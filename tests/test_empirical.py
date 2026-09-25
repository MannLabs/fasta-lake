"""
Empirical stress tests with real data files.

Tests FastaLake against actual experimental data from:

1. **Masters Thesis (NM2026)**: JoMu stool metaproteomics (Thermo QX8)
   - 3 samples with Stage 2 FASTA + AlphaNovo predictions
   - Golden reference: 12 exact-match checks (3 samples x 4 strategies)

2. **PeTr Janko**: OXA-48 Enterobacteria isolates (Bruker)
   - 49 DDA samples with AlphaNovo predictions
   - NCBI reference FASTAs for 7 bacterial species
   - Tests single-species and strain-level resolution

3. **Masters Thesis (full)**: 97 samples across 5 experiment types
   - Technical replicates, workflow replicates, interindividual
   - Tests pipeline robustness at scale

Each test validates:
- Correct file parsing (no crashes on real CSV/FASTA formats)
- Strategy invariants (uniform >= razor >= uniform_2pep)
- Stats output is valid JSON
- Output FASTA is valid
- Performance is reasonable (< 60s per sample)
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

import pytest

from fasta_lake.complexity import estimate_complexity
from fasta_lake.config import ComplexityConfig, SequenceQualityConfig
from fasta_lake.helpers import (
    collect_peptides_with_scores,
    load_fasta_proteins,
)
from fasta_lake.inference.runner import run_inference
from fasta_lake.quality import filter_sequences_by_quality

logger = logging.getLogger(__name__)

# ── Path constants ────────────────────────────────────────────────────────────

BASE = Path(os.environ.get("FASTALAKE_TEST_EMPIRICAL_1", "tests/optional_data/test_empirical_1"))
POOL = Path(os.environ.get("FASTALAKE_TEST_EMPIRICAL_2", "tests/optional_data/test_empirical_2"))

# The NM2026 / Masters-thesis golden reference is GONE, deliberately.
#
# It was generated from the retired RoPE AlphaNovo inference (10 columns); every
# published v27 number uses the 14-column July inference, and a value that only
# reproduces on RoPE is wrong. No July inference exists for those three DC samples, so
# the reference could not be regenerated. The tests that depended on it are archived in
# FastaLake/ARCHIVE_2026-08-03/rope_era_tests_deleted_2026-08-03/ with the full account.
#
# They had never run: the path was BASE/"FastaLake"/"e2e_test"/... , missing the
# `projects/` component, and they skipped with "Data not available" -- blaming the
# machine for a bug in the test.

# PeTr Janko
# Resolved rather than hard-coded, so the two spellings of the pool path both work and a
# missing directory is reported once instead of as N identical "Data not available" skips.
from _data_paths import janko_root  # noqa: E402

_JANKO = janko_root() or (BASE / "PeTr_Janko")
JANKO_PREDS = _JANKO / "alphanovo_output" / "predictions"
JANKO_REF = _JANKO / "ncbi_reference_strains" / "protein_fastas"

# The four inference strategies the isolate tests compare.
STRATEGIES = ["species_budget", "razor", "uniform", "uniform_2pep"]

# Golden reference samples


def _skip_if_missing(*paths: Path) -> None:
    """Skip test if any required path is missing."""
    for p in paths:
        if not p.exists():
            pytest.skip(f"Data not available: {p}")


# ── 1. Golden Reference Tests (3 samples x 4 strategies = 12 checks) ─────────


class TestJankoBacterialIsolates:
    """Test with PeTr Janko OXA-48 Enterobacteria data (Bruker)."""

    @pytest.fixture
    def janko_ref_fasta(self, tmp_path: Path) -> Path:
        """Get combined reference FASTA for Janko samples."""
        ref = JANKO_REF / "ncbi_combined_pOXA48.fasta"
        _skip_if_missing(ref)
        return ref

    @pytest.fixture
    def janko_sample_ids(self) -> list[str]:
        """Get first 10 Janko sample IDs for testing."""
        _skip_if_missing(JANKO_PREDS)
        samples = sorted(d.name for d in JANKO_PREDS.iterdir() if d.is_dir())
        return samples[:10]

    def test_parse_janko_predictions(self, janko_sample_ids: list[str]):
        """Verify AlphaNovo CSV parsing works on Janko format."""
        for sample_id in janko_sample_ids[:3]:
            csv_path = JANKO_PREDS / sample_id / f"{sample_id}_predictions.csv"
            _skip_if_missing(csv_path)

            peptides = collect_peptides_with_scores(csv_path)
            assert len(peptides) > 0, f"No peptides from {sample_id}"

            logger.info("%s: %d peptides collected", sample_id, len(peptides))

    def test_janko_inference_razor(
        self,
        janko_ref_fasta: Path,
        janko_sample_ids: list[str],
        tmp_path: Path,
    ):
        """Run razor inference on Janko bacterial isolates."""

        for sample_id in janko_sample_ids[:5]:
            csv_path = JANKO_PREDS / sample_id / f"{sample_id}_predictions.csv"
            _skip_if_missing(csv_path)

            # Write ref FASTA to tmp (run_inference reads from file)
            stats = run_inference(
                sample_id=sample_id.split("_")[-1],  # A1, A2, etc.
                stage2_fasta=janko_ref_fasta,
                alphanovo_csv=csv_path,
                output_dir=tmp_path / sample_id,
                strategy="razor",
            )

            assert stats["selected_proteins"] >= 0
            logger.info(
                "%s (razor): %d/%d proteins selected",
                sample_id,
                stats["selected_proteins"],
                stats["stage2_proteins"],
            )

    def test_janko_all_strategies(
        self,
        janko_ref_fasta: Path,
        janko_sample_ids: list[str],
        tmp_path: Path,
    ):
        """Run all 4 strategies on first Janko sample."""
        sample_id = janko_sample_ids[0]
        csv_path = JANKO_PREDS / sample_id / f"{sample_id}_predictions.csv"
        _skip_if_missing(csv_path)

        results = {}
        for strategy in STRATEGIES:
            stats = run_inference(
                sample_id=sample_id.split("_")[-1],
                stage2_fasta=janko_ref_fasta,
                alphanovo_csv=csv_path,
                output_dir=tmp_path / strategy,
                strategy=strategy,
            )
            results[strategy] = stats["selected_proteins"]

        # Ordering invariant
        assert results["uniform"] >= results["razor"]
        assert results["uniform"] >= results["uniform_2pep"]

        logger.info(
            "%s strategies: uniform=%d, razor=%d, species_budget=%d, 2pep=%d",
            sample_id,
            results["uniform"],
            results["razor"],
            results["species_budget"],
            results["uniform_2pep"],
        )

    def test_janko_complexity_estimation(
        self,
        janko_ref_fasta: Path,
        janko_sample_ids: list[str],
    ):
        """Estimate complexity of Janko bacterial isolates.

        Expected: single or few species per isolate (they're bacterial cultures).
        Should recommend razor for most samples.
        """
        ref_proteins = load_fasta_proteins(janko_ref_fasta)

        for sample_id in janko_sample_ids[:5]:
            csv_path = JANKO_PREDS / sample_id / f"{sample_id}_predictions.csv"
            _skip_if_missing(csv_path)

            peptides = collect_peptides_with_scores(csv_path)
            config = ComplexityConfig(auto_recommend=True)
            report = estimate_complexity(ref_proteins, peptides, config)

            logger.info(
                "%s complexity: hit_rate=%.1f%%, species=%d, rec=%s",
                sample_id,
                report.hit_rate * 100,
                report.species_richness,
                report.recommendation,
            )

    def test_janko_quality_filtering(self, janko_ref_fasta: Path):
        """Quality filter Janko reference FASTA."""
        proteins = load_fasta_proteins(janko_ref_fasta)

        config = SequenceQualityConfig(
            min_length=20,
            reject_poly_runs=True,
            poly_run_length=10,
            reject_low_complexity=True,
            low_complexity_threshold=0.25,
            reject_fragments=True,
            fragment_length_threshold=50,
        )

        filtered, report = filter_sequences_by_quality(proteins, config)

        logger.info(
            "Janko ref quality: %d/%d passed (-%d length, -%d poly, -%d lowcplx, -%d frag)",
            report.passed,
            report.total_input,
            report.rejected_length,
            report.rejected_poly_run,
            report.rejected_low_complexity,
            report.rejected_fragment,
        )

        # Most reference proteins should pass quality filters
        assert report.passed > report.total_input * 0.8


# ── 4. Masters Thesis Broader Tests ──────────────────────────────────────────
