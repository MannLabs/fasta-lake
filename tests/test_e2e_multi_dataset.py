"""Multi-dataset e2e tests: MicrobPredict + CAPSCAN + Janko.

Tests the pipeline across 3 fundamentally different experiment types:
  - MicrobPredict: gut metaproteomics (233M lake, species_budget)
  - CAPSCAN: human plasma DDA (clinical proteomics)
  - Janko: clinical bacterial isolates (species ID)

Each dataset exercises different aspects of the pipeline and catches
issues specific to its header formats, database sizes, and biology.
"""

import csv
import os
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

FASTA_EXPORT = (
    Path(__file__).parent.parent / "rust" / "fasta_export" / "target" / "release" / "fasta_export"
)


def _count_proteins(path):
    with open(path) as f:
        return sum(1 for line in f if line.startswith(">"))


def _count_headers_with(path, substring):
    with open(path) as f:
        return sum(1 for line in f if line.startswith(">") and substring in line)


# ═══════════════════════════════════════════════════════════════════════════════
# Dataset configs
# ═══════════════════════════════════════════════════════════════════════════════

DATASETS = {
    "microbpredict": {
        "stage2_fasta": Path(
            os.environ.get(
                "FASTALAKE_TEST_E2E_MULTI_DATASET_1", "tests/optional_data/test_e2e_multi_dataset_1"
            )
        ),
        "predictions": Path(
            os.environ.get(
                "FASTALAKE_TEST_E2E_MULTI_DATASET_2", "tests/optional_data/test_e2e_multi_dataset_2"
            )
        ),
        "sample_id": "MicrobPred_11",
        "experiment_type": "metaproteomics",
        "expected_strategy": "species_budget",
        "expected_header_formats": ["GMGC", "MGYG", "sp|"],
        "min_proteins": 1000,
        "max_proteins": 50000,
    },
    "capscan": {
        "stage2_fasta": Path(
            os.environ.get(
                "FASTALAKE_TEST_E2E_MULTI_DATASET_3", "tests/optional_data/test_e2e_multi_dataset_3"
            )
        ),
        "predictions": Path(
            os.environ.get(
                "FASTALAKE_TEST_E2E_MULTI_DATASET_4", "tests/optional_data/test_e2e_multi_dataset_4"
            )
        ),
        "sample_id": "CAPSCAN_1249",
        "experiment_type": "proteomics",
        "expected_strategy": "razor",
        "expected_header_formats": ["sp|", "GMGC"],
        "min_proteins": 1000,
        "max_proteins": 60000,
    },
}


def _dataset_available(name):
    ds = DATASETS[name]
    return ds["stage2_fasta"].exists() and ds["predictions"].exists()


# ═══════════════════════════════════════════════════════════════════════════════
# Per-dataset inference tests
# ═══════════════════════════════════════════════════════════════════════════════


# Manuscript cohorts only. `janko` (pOXA-48 plasmid) and `edwin_mice_gi` (mouse GI)
# were removed 2026-07-31: neither is a manuscript cohort, and edwin_mice_gi's
# predictions are in a 4-column format this pipeline does not read. Healthy stool and
# SIHUMIx are covered by tests/test_e2e_real.py and tests/test_empirical.py.
@pytest.mark.parametrize("dataset_name", ["microbpredict", "capscan"])
class TestMultiDatasetInference:
    """Inference works correctly on 3 fundamentally different datasets."""

    def test_species_budget(self, dataset_name, tmp_path):
        if not _dataset_available(dataset_name):
            pytest.skip(f"{dataset_name} data not available")

        from fasta_lake.inference.runner import run_inference

        ds = DATASETS[dataset_name]

        stats = run_inference(
            sample_id=ds["sample_id"],
            stage2_fasta=ds["stage2_fasta"],
            alphanovo_csv=ds["predictions"],
            output_dir=tmp_path,
            strategy="species_budget",
        )

        assert stats["selected_proteins"] >= ds["min_proteins"], (
            f"{dataset_name}: too few proteins ({stats['selected_proteins']})"
        )
        assert stats["selected_proteins"] <= ds["max_proteins"], (
            f"{dataset_name}: too many proteins ({stats['selected_proteins']})"
        )
        assert stats["unique_peptides"] > 0

    def test_razor(self, dataset_name, tmp_path):
        if not _dataset_available(dataset_name):
            pytest.skip(f"{dataset_name} data not available")

        from fasta_lake.inference.runner import run_inference

        ds = DATASETS[dataset_name]

        stats = run_inference(
            sample_id=ds["sample_id"],
            stage2_fasta=ds["stage2_fasta"],
            alphanovo_csv=ds["predictions"],
            output_dir=tmp_path,
            strategy="razor",
        )

        assert stats["selected_proteins"] > 0
        # Razor should always select fewer than the input
        assert stats["selected_proteins"] < stats["stage2_proteins"]


# ═══════════════════════════════════════════════════════════════════════════════
# Cross-dataset export validation
# ═══════════════════════════════════════════════════════════════════════════════


# Manuscript cohorts only. `janko` (pOXA-48 plasmid) and `edwin_mice_gi` (mouse GI)
# were removed 2026-07-31: neither is a manuscript cohort, and edwin_mice_gi's
# predictions are in a 4-column format this pipeline does not read. Healthy stool and
# SIHUMIx are covered by tests/test_e2e_real.py and tests/test_empirical.py.
@pytest.mark.parametrize("dataset_name", ["microbpredict", "capscan"])
@pytest.mark.parametrize("engine", ["sage", "diann", "msfragger"])
class TestMultiDatasetExport:
    """Export works for all engine × dataset combinations."""

    def test_export(self, dataset_name, engine, tmp_path):
        if not _dataset_available(dataset_name):
            pytest.skip(f"{dataset_name} data not available")
        if not FASTA_EXPORT.exists():
            pytest.skip("fasta_export binary not built")

        ds = DATASETS[dataset_name]
        out = tmp_path / f"{dataset_name}_{engine}.fasta"

        result = subprocess.run(
            [
                str(FASTA_EXPORT),
                "export",
                "-i",
                str(ds["stage2_fasta"]),
                "-o",
                str(out),
                "-e",
                engine,
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"STDERR: {result.stderr}"

        n_input = _count_proteins(ds["stage2_fasta"])
        n_output = _count_proteins(out)

        if engine == "msfragger":
            assert n_output == 2 * n_input, f"{dataset_name}/msfragger: {n_output} != 2×{n_input}"
        else:
            assert n_output == n_input

        # GN= must be present in all headers
        n_gn = _count_headers_with(out, "GN=")
        assert n_gn == n_output, f"{dataset_name}/{engine}: {n_gn}/{n_output} headers have GN="


# ═══════════════════════════════════════════════════════════════════════════════
# Header format coverage
# ═══════════════════════════════════════════════════════════════════════════════


# Manuscript cohorts only. `janko` (pOXA-48 plasmid) and `edwin_mice_gi` (mouse GI)
# were removed 2026-07-31: neither is a manuscript cohort, and edwin_mice_gi's
# predictions are in a 4-column format this pipeline does not read. Healthy stool and
# SIHUMIx are covered by tests/test_e2e_real.py and tests/test_empirical.py.
@pytest.mark.parametrize("dataset_name", ["microbpredict", "capscan"])
class TestHeaderCoverage:
    """Each dataset's FASTA contains the expected header formats."""

    def test_expected_formats_present(self, dataset_name):
        if not _dataset_available(dataset_name):
            pytest.skip(f"{dataset_name} data not available")

        ds = DATASETS[dataset_name]
        for fmt in ds["expected_header_formats"]:
            n = _count_headers_with(ds["stage2_fasta"], fmt)
            assert n > 0, f"{dataset_name}: no headers with '{fmt}' found"


# ═══════════════════════════════════════════════════════════════════════════════
# Prediction CSV validation
# ═══════════════════════════════════════════════════════════════════════════════


# Manuscript cohorts only. `janko` (pOXA-48 plasmid) and `edwin_mice_gi` (mouse GI)
# were removed 2026-07-31: neither is a manuscript cohort, and edwin_mice_gi's
# predictions are in a 4-column format this pipeline does not read. Healthy stool and
# SIHUMIx are covered by tests/test_e2e_real.py and tests/test_empirical.py.
@pytest.mark.parametrize("dataset_name", ["microbpredict", "capscan"])
class TestPredictionCSVs:
    """Prediction CSVs have real peptides (not corrupted)."""

    def test_csv_has_real_peptides(self, dataset_name):
        if not _dataset_available(dataset_name):
            pytest.skip(f"{dataset_name} data not available")

        ds = DATASETS[dataset_name]
        with open(ds["predictions"]) as f:
            reader = csv.DictReader(f)
            peptides = []
            for i, row in enumerate(reader):
                if i >= 500:
                    break
                pep = row.get("peptide_prediction_detokenized_unmodified", "")
                peptides.append(pep)

        lengths = [len(p) for p in peptides]
        avg_len = sum(lengths) / len(lengths)
        long_count = sum(1 for length in lengths if length >= 7)

        assert avg_len > 5.0, f"{dataset_name}: avg peptide length {avg_len:.1f} — likely corrupted"
        assert long_count > 100, f"{dataset_name}: only {long_count}/500 peptides >= 7 AA"
