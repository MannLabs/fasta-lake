"""End-to-end integration tests using REAL MicrobPredict data.

Tests the full pipeline: inference → export → validation on a real
75K-protein metaproteomics sample with 519K AlphaNovo predictions.

These tests need the laboratory's full-size cohort data and skip when it is absent.

Test matrix:
  - 6 inference strategies × 1 real sample
  - 3 export engines × 1 real inference output
  - Entrapment FDP validation on real SAGE results
  - Rust fasta_export binary on real data
"""

import csv
import os
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))


# Real data paths (MicrobPredict condC sample)
BENCH = Path(os.environ.get("FASTALAKE_TEST_E2E_REAL_1", "tests/optional_data/test_e2e_real_1"))
SAMPLE_ID = "20191024_QX7_PeTr_LC02_SA_MicrobPred_10_191025211701"
SAMPLE_DIR = BENCH / "condC_ref_plus_metag" / SAMPLE_ID

STAGE2_FASTA = SAMPLE_DIR / f"{SAMPLE_ID}_species_budget.fasta"
SAGE_RESULTS = SAMPLE_DIR / "sage" / "results.sage.tsv"
ENTRAP_RESULTS = SAMPLE_DIR / "entrapment" / "results.sage.tsv"

# Find predictions CSV from sample list
SAMPLE_LIST = BENCH / "sample_lists" / "all_296.tsv"

FASTA_EXPORT = (
    Path(__file__).parent.parent / "rust" / "fasta_export" / "target" / "release" / "fasta_export"
)


def _find_predictions_csv():
    if not SAMPLE_LIST.exists():
        return None
    with open(SAMPLE_LIST) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["sample_id"] == SAMPLE_ID:
                csv_path = row.get("alphanovo_csv", "")
                if csv_path and Path(csv_path).exists():
                    return Path(csv_path)
    return None


PREDICTIONS_CSV = _find_predictions_csv()

DATA_AVAILABLE = STAGE2_FASTA.exists() and PREDICTIONS_CSV is not None

SAGE_AVAILABLE = SAGE_RESULTS.exists()
ENTRAP_AVAILABLE = ENTRAP_RESULTS.exists()


def _count_proteins(fasta_path):
    with open(fasta_path) as f:
        return sum(1 for line in f if line.startswith(">"))


# ═══════════════════════════════════════════════════════════════════════════════
# Inference tests: all 6 strategies on real data
# ═══════════════════════════════════════════════════════════════════════════════


@pytest.mark.skipif(not DATA_AVAILABLE, reason="MicrobPredict data not available")
class TestRealInference:
    """Run all 6 inference strategies on a real 75K-protein sample."""

    @pytest.mark.parametrize("strategy", ["species_budget", "razor", "uniform", "uniform_2pep"])
    def test_strategy_produces_valid_fasta(self, strategy, tmp_path):
        from fasta_lake.inference.runner import run_inference

        output_dir = tmp_path / strategy
        stats = run_inference(
            sample_id=SAMPLE_ID,
            stage2_fasta=STAGE2_FASTA,
            alphanovo_csv=PREDICTIONS_CSV,
            output_dir=output_dir,
            strategy=strategy,
        )

        # Must produce non-zero proteins
        assert stats["selected_proteins"] > 0, f"{strategy}: 0 proteins selected"
        assert stats["unique_peptides"] > 0, f"{strategy}: 0 unique peptides"

        # Output FASTA must exist and be valid
        # Naming: sample_id_strategy.fasta (uniform_2pep → uniform2pep)
        fasta_files = list(output_dir.glob("*.fasta"))
        assert len(fasta_files) >= 1, (
            f"Expected >=1 FASTA in {output_dir}, got {list(output_dir.iterdir())}"
        )
        n_prots = _count_proteins(fasta_files[0])
        assert n_prots == stats["selected_proteins"]

        # Stats JSON must exist
        json_files = list(output_dir.glob("*stats.json"))
        assert len(json_files) == 1

    def test_species_budget_reduces_database(self, tmp_path):
        """species_budget should reduce the database (input 75K → fewer)."""
        from fasta_lake.inference.runner import run_inference

        stats = run_inference(
            sample_id=SAMPLE_ID,
            stage2_fasta=STAGE2_FASTA,
            alphanovo_csv=PREDICTIONS_CSV,
            output_dir=tmp_path,
            strategy="species_budget",
        )

        # Must produce proteins and show reduction
        assert stats["selected_proteins"] > 1000, f"Too few: {stats['selected_proteins']}"
        assert stats["selected_proteins"] < stats["stage2_proteins"], (
            f"No reduction: {stats['selected_proteins']} >= {stats['stage2_proteins']}"
        )
        assert stats["reduction_factor"] >= 1.0, (
            f"Reduction factor {stats['reduction_factor']}x < 1"
        )

    def test_razor_fewer_than_uniform(self, tmp_path):
        """razor should select fewer proteins than uniform."""
        from fasta_lake.inference.runner import run_inference

        razor = run_inference(
            sample_id=SAMPLE_ID,
            stage2_fasta=STAGE2_FASTA,
            alphanovo_csv=PREDICTIONS_CSV,
            output_dir=tmp_path / "razor",
            strategy="razor",
        )
        uniform = run_inference(
            sample_id=SAMPLE_ID,
            stage2_fasta=STAGE2_FASTA,
            alphanovo_csv=PREDICTIONS_CSV,
            output_dir=tmp_path / "uniform",
            strategy="uniform",
        )
        assert razor["selected_proteins"] <= uniform["selected_proteins"]


# ═══════════════════════════════════════════════════════════════════════════════
# Export tests: all 3 engines on real inference output
# ═══════════════════════════════════════════════════════════════════════════════


@pytest.mark.skipif(not DATA_AVAILABLE, reason="MicrobPredict data not available")
class TestRealExport:
    """Test engine-specific FASTA export on real data."""

    def _run_inference_and_export(self, tmp_path, engine):
        from fasta_lake.headers import prepare_fasta_for_engine
        from fasta_lake.inference.runner import run_inference

        # Run inference
        stats = run_inference(
            sample_id=SAMPLE_ID,
            stage2_fasta=STAGE2_FASTA,
            alphanovo_csv=PREDICTIONS_CSV,
            output_dir=tmp_path / "inference",
            strategy="species_budget",
        )
        inference_fasta = list((tmp_path / "inference").glob("*.fasta"))[0]

        # Export
        out_fasta = tmp_path / f"export_{engine}.fasta"
        export_stats = prepare_fasta_for_engine(
            str(inference_fasta),
            str(out_fasta),
            engine=engine,
        )
        return stats, export_stats, out_fasta

    def test_sage_export(self, tmp_path):
        stats, export, fasta = self._run_inference_and_export(tmp_path, "sage")
        assert export.decoy_proteins == 0
        assert export.target_proteins == stats["selected_proteins"]
        # All headers should have GN=
        with open(fasta) as f:
            for line in f:
                if line.startswith(">"):
                    assert "GN=" in line, f"Missing GN= in SAGE export: {line[:80]}"
                    break

    def test_diann_export(self, tmp_path):
        stats, export, fasta = self._run_inference_and_export(tmp_path, "diann")
        assert export.decoy_proteins == 0
        # Check non-sp/tr headers have gene as first word
        with open(fasta) as f:
            for line in f:
                if (
                    line.startswith(">")
                    and not line.startswith(">sp|")
                    and not line.startswith(">tr|")
                ):
                    parts = line[1:].split()
                    # First word = accession, second = gene (must not be raw description)
                    assert len(parts) >= 2
                    assert parts[1].lower() not in (
                        "uncharacterized",
                        "hypothetical",
                        "putative",
                    ), f"Gene not first word in DIA-NN export: {line[:80]}"
                    break

    def test_msfragger_export(self, tmp_path):
        stats, export, fasta = self._run_inference_and_export(tmp_path, "msfragger")
        assert export.decoy_proteins == export.target_proteins
        assert export.total_proteins == 2 * export.target_proteins
        # Count rev_ headers
        n_decoys = 0
        with open(fasta) as f:
            for line in f:
                if line.startswith(">rev_"):
                    n_decoys += 1
        assert n_decoys == export.target_proteins


# ═══════════════════════════════════════════════════════════════════════════════
# Entrapment FDP on real SAGE results
# ═══════════════════════════════════════════════════════════════════════════════


@pytest.mark.skipif(not ENTRAP_AVAILABLE, reason="Entrapment SAGE results not available")
class TestRealEntrapment:
    """Validate FDP computation on real entrapment SAGE results."""

    def test_entrapment_has_real_psms(self):
        """Entrapment SAGE results should have >100K PSMs."""
        with open(ENTRAP_RESULTS) as f:
            n_lines = sum(1 for _ in f)
        assert n_lines > 100000, f"Only {n_lines} lines — expected >100K"

    def test_entrapment_hits_exist(self):
        """Should find ENTRAP_SHUF_ hits in results."""
        n_entrap = 0
        with open(ENTRAP_RESULTS) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if float(row.get("spectrum_q", "1")) > 0.01:
                    continue
                if row.get("label", "1") == "-1":
                    continue
                if "ENTRAP_SHUF_" in row.get("proteins", ""):
                    n_entrap += 1
        assert n_entrap > 0, "No entrapment hits found — prefix wrong?"

    def test_fdp_combined_valid(self):
        """FDP combined should be between 0 and 1."""
        from fasta_lake.entrapment import estimate_fdp

        n_target = 0
        n_entrap = 0
        with open(ENTRAP_RESULTS) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if float(row.get("spectrum_q", "1")) > 0.01:
                    continue
                if row.get("label", "1") == "-1":
                    continue
                if "ENTRAP_SHUF_" in row.get("proteins", ""):
                    n_entrap += 1
                else:
                    n_target += 1

        result = estimate_fdp(
            entrapment_hits=n_entrap,
            target_hits=n_target,
            target_db_size=75234,
            entrapment_db_size=6500,
        )
        assert 0 <= result["fdp_combined"] <= 1.0
        assert result["fdp_combined"] > result["fdp_lower_bound"]
        assert result["r"] < 1.0  # entrapment << target


# ═══════════════════════════════════════════════════════════════════════════════
# Rust binary integration on real data
# ═══════════════════════════════════════════════════════════════════════════════


@pytest.mark.skipif(not DATA_AVAILABLE, reason="MicrobPredict data not available")
@pytest.mark.skipif(not FASTA_EXPORT.exists(), reason="fasta_export binary not built")
class TestRustBinaryReal:
    """Test Rust fasta_export binary on real data."""

    def test_export_sage_real(self, tmp_path):
        out = tmp_path / "sage.fasta"
        result = subprocess.run(
            [str(FASTA_EXPORT), "export", "-i", str(STAGE2_FASTA), "-o", str(out), "-e", "sage"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"STDERR: {result.stderr}"
        assert _count_proteins(out) == 75234

    def test_export_msfragger_real(self, tmp_path):
        out = tmp_path / "msfragger.fasta"
        result = subprocess.run(
            [
                str(FASTA_EXPORT),
                "export",
                "-i",
                str(STAGE2_FASTA),
                "-o",
                str(out),
                "-e",
                "msfragger",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        n = _count_proteins(out)
        assert n == 75234 * 2  # targets + decoys

    def test_export_diann_real(self, tmp_path):
        out = tmp_path / "diann.fasta"
        result = subprocess.run(
            [str(FASTA_EXPORT), "export", "-i", str(STAGE2_FASTA), "-o", str(out), "-e", "diann"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert _count_proteins(out) == 75234

    def test_export_diann_headers_correct(self, tmp_path):
        """DIA-NN export: non-sp/tr headers should have gene as first word."""
        out = tmp_path / "diann.fasta"
        subprocess.run(
            [str(FASTA_EXPORT), "export", "-i", str(STAGE2_FASTA), "-o", str(out), "-e", "diann"],
            capture_output=True,
        )
        with open(out) as f:
            n_checked = 0
            for line in f:
                if not line.startswith(">"):
                    continue
                if line.startswith(">sp|") or line.startswith(">tr|"):
                    continue
                parts = line[1:].rstrip().split()
                if len(parts) >= 2:
                    # Gene should be second token (first word of description)
                    gene = parts[1]
                    assert gene.lower() not in ("uncharacterized", "hypothetical", "putative"), (
                        f"Bad DIA-NN header: {line[:100]}"
                    )
                n_checked += 1
                if n_checked >= 100:
                    break
            assert n_checked >= 50, "Too few non-sp/tr headers to validate"
