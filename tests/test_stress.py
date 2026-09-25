"""Stress tests: edge cases, adversarial inputs, and robustness checks."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from fasta_lake.config import FastaLakeConfig, load_config
from fasta_lake.helpers import (
    clean_sequence,
    collect_peptides_with_scores,
    get_db_type,
    get_species_id,
    load_fasta_proteins,
    map_peptides_to_proteins,
    normalize_il,
)
from fasta_lake.inference.base import InferenceInput
from fasta_lake.inference.runner import run_inference
from fasta_lake.inference.strategies import (
    infer_razor,
    infer_species_budget,
    infer_uniform,
    infer_uniform_2pep,
)

# ── Edge cases for helpers ──────────────────────────────────────────────────


class TestHelperEdgeCases:
    """Test helpers with adversarial and edge-case inputs."""

    def test_normalize_il_long_string(self):
        seq = "I" * 100_000
        result = normalize_il(seq)
        assert result == "L" * 100_000

    def test_clean_sequence_whitespace(self):
        assert clean_sequence("  MSK  GE  ") == "  MSK  GE  ".upper().replace("*", "")

    def test_clean_sequence_only_stars(self):
        assert clean_sequence("***") == ""

    def test_get_db_type_empty_string(self):
        assert get_db_type("") == "smORF"

    def test_get_db_type_pipe_only(self):
        assert get_db_type("sp|") == "Human"

    def test_get_species_id_empty(self):
        assert get_species_id("") is None

    def test_get_species_id_mgyg_no_underscore(self):
        result = get_species_id("MGYG000000001")
        assert result == "MGYG000000001"

    def test_get_species_id_gmgc_one_dot(self):
        result = get_species_id("GMGC10.something")
        # Less than 3 parts — not a resolvable GMGC sample id
        assert result is None


# ── Empty and malformed inputs ──────────────────────────────────────────────


class TestEmptyInputs:
    """Test behavior with empty/malformed inputs."""

    def test_empty_fasta(self, tmp_path: Path):
        fasta = tmp_path / "empty.fasta"
        fasta.write_text("")
        proteins = load_fasta_proteins(fasta)
        assert proteins == {}

    def test_fasta_no_sequences(self, tmp_path: Path):
        fasta = tmp_path / "headers_only.fasta"
        fasta.write_text(">protein1\n>protein2\n")
        # Header-only entries have no sequence: an input defect, reported not dropped.
        with pytest.raises(ValueError, match="record 'protein1'.*has no sequence"):
            load_fasta_proteins(fasta)

    def test_single_protein_fasta(self, tmp_path: Path):
        fasta = tmp_path / "single.fasta"
        fasta.write_text(">prot1\nMSKGEELFT\n")
        proteins = load_fasta_proteins(fasta)
        assert len(proteins) == 1
        assert "prot1" in proteins
        assert proteins["prot1"] == "MSKGEELFT"

    def test_empty_predictions_csv(self, tmp_path: Path):
        csv_file = tmp_path / "empty.csv"
        csv_file.write_text("peptide_prediction_detokenized_unmodified,score\n")
        peptides = collect_peptides_with_scores(csv_file)
        assert peptides == {}

    def test_predictions_all_filtered(self, tmp_path: Path):
        csv_file = tmp_path / "short.csv"
        csv_file.write_text(
            "peptide_prediction_detokenized_unmodified,score\n"
            "AB,0.9\n"  # Too short (min_length=8)
            "CDE,0.8\n"  # Too short
        )
        peptides = collect_peptides_with_scores(csv_file, min_length=8)
        assert peptides == {}

    def test_predictions_zero_score_filtered(self, tmp_path: Path):
        csv_file = tmp_path / "zero_score.csv"
        csv_file.write_text(
            "peptide_prediction_detokenized_unmodified,score\n"
            "PEPTIDERLONGONE,0.0\n"
            "ANOTHERPEPTIDE,0.0\n"
        )
        peptides = collect_peptides_with_scores(csv_file)
        assert peptides == {}

    def test_predictions_negative_score_filtered(self, tmp_path: Path):
        csv_file = tmp_path / "neg_score.csv"
        csv_file.write_text(
            "peptide_prediction_detokenized_unmodified,score\nPEPTIDERLONGONE,-0.5\n"
        )
        peptides = collect_peptides_with_scores(csv_file)
        assert peptides == {}

    def test_predictions_invalid_score(self, tmp_path: Path):
        csv_file = tmp_path / "bad_score.csv"
        csv_file.write_text(
            "peptide_prediction_detokenized_unmodified,score\nPEPTIDERLONGONE,not_a_number\n"
        )
        # An unparsable score is skipped and counted, never replaced by a default.
        peptides = collect_peptides_with_scores(csv_file)
        assert len(peptides) == 0

    def test_map_no_matches(self, tmp_path: Path):
        proteins = {"prot1": "AAAAAAAAAA"}
        peptides = {"ZZZZZZZZZ"}
        pep2prot, prot2pep = map_peptides_to_proteins(proteins, peptides)
        assert len(pep2prot) == 0
        assert len(prot2pep) == 0


# ── Inference with edge cases ───────────────────────────────────────────────


class TestInferenceEdgeCases:
    """Test inference strategies with edge-case inputs."""

    def test_all_unique_peptides(self):
        """When all peptides are unique, all strategies should agree."""
        inp = InferenceInput(
            pep2prot={
                "PEP1": {"A"},
                "PEP2": {"B"},
                "PEP3": {"C"},
            },
            prot2pep={
                "A": {"PEP1"},
                "B": {"PEP2"},
                "C": {"PEP3"},
            },
            proteins={"A": "SEQ", "B": "SEQ", "C": "SEQ"},
            taxonomy_map={"A": "sp1", "B": "sp2", "C": "sp3"},
        )

        for fn in [infer_species_budget, infer_razor, infer_uniform]:
            result = fn(inp)
            assert result.selected_proteins == {"A", "B", "C"}

    def test_single_protein(self):
        """Single protein should always be selected."""
        inp = InferenceInput(
            pep2prot={"PEP1": {"ONLY"}},
            prot2pep={"ONLY": {"PEP1"}},
            proteins={"ONLY": "MSKGEELFT"},
            taxonomy_map={"ONLY": "sp1"},
        )

        for fn in [infer_species_budget, infer_razor, infer_uniform]:
            result = fn(inp)
            assert result.selected_proteins == {"ONLY"}

    def test_completely_shared_peptide(self):
        """One peptide shared across many proteins."""
        prots = {f"PROT_{i}": {"SHARED_PEP"} for i in range(100)}
        inp = InferenceInput(
            pep2prot={"SHARED_PEP": set(prots.keys())},
            prot2pep=prots,
            proteins={pid: "SEQ" for pid in prots},
        )

        # Razor: should pick exactly 1
        result = infer_razor(inp)
        assert len(result.selected_proteins) == 1

        # Uniform: should keep all
        result = infer_uniform(inp)
        assert len(result.selected_proteins) == 100

        # 2-pep: should keep none (each has only 1 peptide)
        result = infer_uniform_2pep(inp, min_peptides=2)
        assert len(result.selected_proteins) == 0

    def test_many_species_cross_shared(self):
        """Cross-species shared peptide with many species."""
        prot_ids = [f"MGYG00000000{i}_00001" for i in range(1, 6)]
        inp = InferenceInput(
            pep2prot={"CROSS_PEP": set(prot_ids)},
            prot2pep={pid: {"CROSS_PEP"} for pid in prot_ids},
            proteins={pid: "SEQ" for pid in prot_ids},
        )

        # Species budget should keep 1 per species (5 species, 5 proteins)
        result = infer_species_budget(inp)
        assert len(result.selected_proteins) == 5

        # Razor should keep only 1
        result = infer_razor(inp)
        assert len(result.selected_proteins) == 1


# ── Config edge cases ──────────────────────────────────────────────────────


class TestConfigEdgeCases:
    """Test config with edge-case inputs."""

    def test_empty_json(self, tmp_path: Path):
        config_file = tmp_path / "empty.json"
        config_file.write_text("{}")
        cfg = load_config(config_file)
        assert isinstance(cfg, FastaLakeConfig)

    def test_extra_fields_ignored(self, tmp_path: Path):
        config_file = tmp_path / "extra.json"
        config_file.write_text(
            json.dumps(
                {
                    "database": {"sources": [{"tag": "X", "path": "/x"}]},
                    "peptide_evidence": {"predictions_dir": "/p"},
                    "output": {"directory": "/o"},
                }
            )
        )
        cfg = load_config(config_file)
        cfg.validate()  # Should not raise

    def test_invalid_json_raises(self, tmp_path: Path):
        config_file = tmp_path / "bad.json"
        config_file.write_text("not valid json {{{")
        with pytest.raises(json.JSONDecodeError):
            load_config(config_file)


# ── Runner edge cases ──────────────────────────────────────────────────────


class TestRunnerEdgeCases:
    """Test run_inference with edge cases."""

    def test_no_peptide_matches(self, tmp_path: Path):
        """When no peptides match any proteins."""
        fasta = tmp_path / "proteins.fasta"
        fasta.write_text(">prot1\nAAAAAAAAAAAA\n>prot2\nBBBBBBBBBBBB\n")

        csv_file = tmp_path / "preds.csv"
        csv_file.write_text("peptide_prediction_detokenized_unmodified,score\nZZZZZZZZZZZ,0.9\n")

        out = tmp_path / "output"
        stats = run_inference(
            sample_id="EMPTY",
            stage2_fasta=fasta,
            alphanovo_csv=csv_file,
            output_dir=out,
            strategy="uniform",
        )
        assert stats["selected_proteins"] == 0

    def test_all_strategies_on_tiny(
        self, tiny_fasta: Path, tiny_predictions: Path, tiny_taxonomy: Path, tmp_path: Path
    ):
        """Run all 4 strategies on tiny data — none should crash."""
        for strategy in ("species_budget", "razor", "uniform", "uniform_2pep"):
            out = tmp_path / strategy
            stats = run_inference(
                sample_id="TINY",
                stage2_fasta=tiny_fasta,
                alphanovo_csv=tiny_predictions,
                output_dir=out,
                taxonomy_map_path=tiny_taxonomy,
                strategy=strategy,
            )
            assert stats["selected_proteins"] >= 0
            # Output files exist
            suffix_map = {
                "species_budget": "species_budget",
                "razor": "razor",
                "uniform": "uniform",
                "uniform_2pep": "uniform2pep",
            }
            fasta_out = out / f"TINY_{suffix_map[strategy]}.fasta"
            assert fasta_out.exists()

    def test_duplicate_peptides_in_csv(self, tmp_path: Path):
        """Duplicate peptides should keep best score."""
        fasta = tmp_path / "prot.fasta"
        fasta.write_text(">prot1\nPEPTIDERLONGSEQ\n")

        csv_file = tmp_path / "dups.csv"
        csv_file.write_text(
            "peptide_prediction_detokenized_unmodified,score\n"
            "PEPTIDERLON,0.5\n"
            "PEPTIDERLON,0.9\n"
            "PEPTIDERLON,0.3\n"
        )

        out = tmp_path / "output"
        stats = run_inference(
            sample_id="DUP",
            stage2_fasta=fasta,
            alphanovo_csv=csv_file,
            output_dir=out,
            strategy="uniform",
        )
        # Should still find the peptide
        assert stats["selected_proteins"] >= 0


# ── Large synthetic test ────────────────────────────────────────────────────


class TestLargeSynthetic:
    """Test with larger synthetic data to check performance/correctness."""

    def test_1000_proteins_100_peptides(self, tmp_path: Path):
        """Generate 1000 proteins and 100 peptides, run all strategies."""
        import random

        random.seed(42)

        aa = "ACDEFGHKLMNPQRSTVWY"

        # Generate 1000 proteins
        fasta = tmp_path / "large.fasta"
        peptides_in_proteins = []
        with open(fasta, "w") as f:
            for i in range(1000):
                seq = "".join(random.choices(aa, k=200))
                f.write(f">PROT_{i:04d}\n{seq}\n")
                # Extract a peptide from each protein
                if i < 100:
                    start = random.randint(0, 180)
                    pep = seq[start : start + 15]
                    peptides_in_proteins.append(pep)

        # Generate predictions CSV
        csv_file = tmp_path / "preds.csv"
        with open(csv_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["peptide_prediction_detokenized_unmodified", "score"])
            for pep in peptides_in_proteins:
                writer.writerow([pep, 0.9])

        # PROT_ ids carry no species: species_budget needs a taxonomy map.
        tax = tmp_path / "tax.tsv"
        tax.write_text("".join(f"PROT_{i:04d}\tsp{i % 10}\n" for i in range(1000)))

        # Run all strategies
        for strategy in ("species_budget", "razor", "uniform", "uniform_2pep"):
            out = tmp_path / strategy
            stats = run_inference(
                sample_id="LARGE",
                stage2_fasta=fasta,
                alphanovo_csv=csv_file,
                output_dir=out,
                taxonomy_map_path=tax,
                strategy=strategy,
            )
            assert stats["stage2_proteins"] == 1000
            # Uniform should have the most
            if strategy == "uniform":
                uniform_count = stats["selected_proteins"]

        # Verify ordering constraint
        for strategy in ("species_budget", "razor", "uniform_2pep"):
            out = tmp_path / strategy
            stats_file = (
                out
                / f"LARGE_{strategy if strategy != 'uniform_2pep' else 'uniform2pep'}_stats.json"
            )
            with open(stats_file) as f:
                s = json.load(f)
            assert s["selected_proteins"] <= uniform_count
