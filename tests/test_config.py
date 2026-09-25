"""T5-T6: Config loading and validation tests."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from fasta_lake.config import (
    ComplexityConfig,
    DatabaseSource,
    EntrapmentConfig,
    FastaLakeConfig,
    InferenceConfig,
    SequenceQualityConfig,
    create_example_config,
    load_config,
)

# ── T5: Config defaults ─────────────────────────────────────────────────────


class TestConfigDefaults:
    """T5: Minimal JSON loads with correct defaults."""

    def test_minimal_config_loads(self, tmp_path: Path, sample_config_dict: dict):
        config_file = tmp_path / "minimal.json"
        config_file.write_text(json.dumps(sample_config_dict))

        cfg = load_config(config_file)

        assert len(cfg.database.sources) == 1
        assert cfg.database.sources[0].tag == "REF"
        assert cfg.inference.strategy == "species_budget"
        assert cfg.inference.min_peptides == 1
        assert cfg.peptide_evidence.min_length == 9  # v5 default
        assert cfg.peptide_evidence.max_length == 50
        assert cfg.peptide_evidence.normalize_il is True
        assert cfg.clustering.enabled is True
        assert cfg.clustering.identity == 0.97
        assert cfg.output.write_stats is True

    def test_full_config_loads(self, tmp_path: Path, full_config_dict: dict):
        config_file = tmp_path / "full.json"
        config_file.write_text(json.dumps(full_config_dict))

        cfg = load_config(config_file)
        assert cfg.inference.strategy == "species_budget"
        assert cfg.clustering.identity == 0.97

    def test_config_not_found_raises(self):
        with pytest.raises(FileNotFoundError):
            load_config("/nonexistent/path.json")

    def test_create_example_minimal(self):
        config_str = create_example_config("minimal")
        data = json.loads(config_str)
        assert "database" in data
        assert "sources" in data["database"]

    def test_create_example_full(self):
        config_str = create_example_config("full")
        data = json.loads(config_str)
        assert "inference" in data
        assert "clustering" in data

    def test_validate_passes(self, tmp_path: Path, sample_config_dict: dict):
        config_file = tmp_path / "valid.json"
        config_file.write_text(json.dumps(sample_config_dict))
        cfg = load_config(config_file)
        cfg.validate()  # Should not raise

    def test_validate_no_sources_raises(self, tmp_path: Path):
        data = {
            "database": {"sources": []},
            "peptide_evidence": {"predictions_dir": "./preds"},
            "output": {"directory": "./out"},
        }
        config_file = tmp_path / "nosrc.json"
        config_file.write_text(json.dumps(data))
        cfg = load_config(config_file)
        with pytest.raises(ValueError, match="source"):
            cfg.validate()


# ── T6: Config validation ───────────────────────────────────────────────────


class TestConfigValidation:
    """T6: Invalid strategy name raises ValueError."""

    def test_invalid_strategy_raises(self):
        with pytest.raises(ValueError, match="Invalid strategy"):
            InferenceConfig(strategy="invalid_strategy")

    def test_valid_strategies_accepted(self):
        for strategy in ("species_budget", "razor", "uniform", "uniform_2pep"):
            cfg = InferenceConfig(strategy=strategy)
            assert cfg.strategy == strategy

    def test_invalid_strategy_in_full_config(self, tmp_path: Path, sample_config_dict: dict):
        sample_config_dict["inference"] = {"strategy": "bogus"}
        config_file = tmp_path / "bad.json"
        config_file.write_text(json.dumps(sample_config_dict))
        with pytest.raises(ValueError, match="Invalid strategy"):
            load_config(config_file)


# ── v3 Config features ────────────────────────────────────────────────────────


class TestV3ConfigFeatures:
    """Test v3 config: source types, quality, complexity, entrapment."""

    def test_database_source_types(self):
        for src_type in ("reference", "variant", "assembly", "contaminant", "entrapment"):
            s = DatabaseSource(tag="T", path="/p", type=src_type)
            assert s.type == src_type

    def test_invalid_source_type_raises(self):
        with pytest.raises(ValueError, match="Invalid source type"):
            DatabaseSource(tag="T", path="/p", type="invalid")

    def test_get_sources_by_type(self):
        cfg = FastaLakeConfig(
            database={
                "sources": [
                    {"tag": "SP", "path": "/sp", "type": "reference"},
                    {"tag": "cRAP", "path": "/crap", "type": "contaminant"},
                    {"tag": "Entrap", "path": "/e", "type": "entrapment"},
                ]
            }
        )
        assert len(cfg.get_sources_by_type("reference")) == 1
        assert len(cfg.get_sources_by_type("contaminant")) == 1
        assert len(cfg.get_sources_by_type("entrapment")) == 1
        assert len(cfg.get_sources_by_type("variant")) == 0

    def test_has_entrapment(self):
        cfg = FastaLakeConfig(
            database={
                "sources": [
                    {"tag": "E", "path": "/e", "type": "entrapment"},
                ]
            }
        )
        assert cfg.has_entrapment() is True

    def test_has_entrapment_from_config(self):
        cfg = FastaLakeConfig(entrapment={"enabled": True})
        assert cfg.has_entrapment() is True

    def test_has_variants(self):
        cfg = FastaLakeConfig(
            database={
                "sources": [
                    {"tag": "gnomAD", "path": "/g", "type": "variant"},
                ]
            }
        )
        assert cfg.has_variants() is True

    def test_has_contaminants(self):
        cfg = FastaLakeConfig(
            database={
                "sources": [
                    {"tag": "cRAP", "path": "/c", "type": "contaminant"},
                ]
            }
        )
        assert cfg.has_contaminants() is True

    def test_sequence_quality_defaults(self):
        cfg = SequenceQualityConfig()
        # v5.1: min_length 6→30 (drops fragmented ORFs), poly-runs on by default
        assert cfg.min_length == 30
        assert cfg.reject_poly_runs is True
        assert cfg.reject_low_complexity is False
        assert cfg.reject_fragments is False
        assert cfg.strip_noncanonical is True  # clean ambiguous B/J/Z; preserve U/O
        assert cfg.reject_if_noncanonical_fraction == 0.05

    def test_complexity_defaults(self):
        cfg = ComplexityConfig()
        assert cfg.enabled is False
        assert cfg.auto_recommend is True
        assert "hit_rate" in cfg.report_metrics

    def test_entrapment_defaults(self):
        cfg = EntrapmentConfig()
        assert cfg.enabled is False
        assert cfg.method == "phylogenetic"
        assert cfg.expected_fdr == 0.01

    def test_inference_v3_defaults(self):
        cfg = InferenceConfig()
        assert cfg.strain_resolution is False
        assert cfg.keep_contaminants is True
        assert cfg.keep_entrapment is True

    def test_clinical_template(self):
        config_str = create_example_config("clinical")
        data = json.loads(config_str)
        sources = data["database"]["sources"]
        types = {s["type"] for s in sources}
        assert "reference" in types
        assert "variant" in types
        assert "contaminant" in types
        assert "entrapment" in types
        assert data["inference"]["strain_resolution"] is True

    def test_metaproteomics_template(self):
        config_str = create_example_config("metaproteomics")
        data = json.loads(config_str)
        assert data["inference"]["strategy"] == "species_budget"
        assert data["entrapment"]["enabled"] is True
        assert data["complexity"]["enabled"] is True

    def test_full_config_with_v3_fields(self, tmp_path: Path):
        data = {
            "database": {
                "sources": [
                    {"tag": "SP", "path": "/sp", "type": "reference"},
                    {"tag": "gnomAD", "path": "/g", "type": "variant"},
                    {"tag": "cRAP", "path": "/c", "type": "contaminant"},
                ]
            },
            "sequence_quality": {
                "reject_poly_runs": True,
                "poly_run_length": 8,
                "reject_low_complexity": True,
                "low_complexity_threshold": 0.25,
            },
            "peptide_evidence": {
                "predictions_dir": "/preds",
                "min_score": 0.5,
            },
            "complexity": {
                "enabled": True,
                "auto_recommend": True,
            },
            "inference": {
                "strategy": "razor",
                "strain_resolution": True,
                "keep_contaminants": True,
            },
            "entrapment": {
                "enabled": True,
                "method": "shuffled",
                "n_proteins": 500,
            },
            "output": {
                "directory": "/out",
                "write_audit_log": True,
                "provenance_tracking": True,
            },
        }
        config_file = tmp_path / "v3.json"
        config_file.write_text(json.dumps(data))

        cfg = load_config(config_file)
        cfg.validate()

        assert cfg.sequence_quality.reject_poly_runs is True
        assert cfg.sequence_quality.poly_run_length == 8
        assert cfg.peptide_evidence.min_score == 0.5
        assert cfg.complexity.enabled is True
        assert cfg.inference.strain_resolution is True
        assert cfg.entrapment.enabled is True
        assert cfg.entrapment.method == "shuffled"
        assert cfg.output.write_audit_log is True
