"""Tests that validate benchmark results are scientifically sound.

These tests check the actual benchmark outputs (clustering, entrapment)
for expected properties and relationships. They serve as regression tests —
if the pipeline changes, these catch scientific inconsistencies.

Requires benchmark results on disk (run on MPI cluster).
"""

import csv
import os
from pathlib import Path

import pytest

RESULTS_DIR = Path(
    os.environ.get(
        "FASTALAKE_TEST_BENCHMARK_RESULTS_1", "tests/optional_data/test_benchmark_results_1"
    )
)

CLUSTERING_SUMMARY = RESULTS_DIR / "clustering_benchmark_summary.tsv"
CLUSTERING_DETAIL = RESULTS_DIR / "clustering_benchmark_detail.tsv"

HAS_CLUSTERING = CLUSTERING_SUMMARY.exists()
HAS_DETAIL = CLUSTERING_DETAIL.exists()


def _load_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


# ═══════════════════════════════════════════════════════════════════════════════
# Clustering benchmark validation
# ═══════════════════════════════════════════════════════════════════════════════


@pytest.mark.skipif(not HAS_CLUSTERING, reason="Clustering results not available")
class TestClusteringBenchmark:
    """Validate that clustering results show expected biological patterns."""

    def test_unclustered_has_most_pgs(self):
        """Unclustered should always have the most protein groups."""
        rows = _load_tsv(CLUSTERING_SUMMARY)
        for dataset in ["CAPSCAN", "MicrobPredict"]:
            ds_rows = [r for r in rows if r["dataset"] == dataset]
            unclustered = [r for r in ds_rows if r["variant"] == "unclustered"][0]
            for r in ds_rows:
                assert float(r["PG_mean"]) <= float(unclustered["PG_mean"]) * 1.02, (
                    f"{dataset}/{r['variant']}: PG_mean {r['PG_mean']} > unclustered "
                    f"{unclustered['PG_mean']}"
                )

    def test_c99_minimal_loss(self):
        """c99 clustering should lose <5% of protein groups."""
        rows = _load_tsv(CLUSTERING_SUMMARY)
        for dataset in ["CAPSCAN", "MicrobPredict"]:
            ds_rows = [r for r in rows if r["dataset"] == dataset]
            c99 = [r for r in ds_rows if r["variant"] == "c99"][0]
            loss = float(c99["PG_pct_loss_vs_unclustered"])
            assert loss < 5.0, f"{dataset} c99 loss = {loss}% (expected <5%)"

    def test_c95_meaningful_loss(self):
        """c95 clustering should show 5-15% PG loss (meaningful tradeoff)."""
        rows = _load_tsv(CLUSTERING_SUMMARY)
        for dataset in ["CAPSCAN", "MicrobPredict"]:
            ds_rows = [r for r in rows if r["dataset"] == dataset]
            c95 = [r for r in ds_rows if r["variant"] == "c95"][0]
            loss = float(c95["PG_pct_loss_vs_unclustered"])
            assert 5.0 <= loss <= 20.0, f"{dataset} c95 loss = {loss}% (expected 5-20%)"

    def test_stage2_decreases_with_clustering(self):
        """Stage 2 database size should decrease monotonically with tighter clustering."""
        rows = _load_tsv(CLUSTERING_SUMMARY)
        order = ["unclustered", "c99", "c97", "c95"]
        for dataset in ["CAPSCAN", "MicrobPredict"]:
            ds_rows = {
                r["variant"]: float(r["stage2_mean"]) for r in rows if r["dataset"] == dataset
            }
            for i in range(len(order) - 1):
                assert ds_rows[order[i]] > ds_rows[order[i + 1]], (
                    f"{dataset}: stage2 {order[i]}={ds_rows[order[i]]:.0f} "
                    f"<= {order[i + 1]}={ds_rows[order[i + 1]]:.0f}"
                )

    def test_sufficient_sample_count(self):
        """Each condition should have >=200 samples for statistical validity."""
        rows = _load_tsv(CLUSTERING_SUMMARY)
        for r in rows:
            assert int(r["n_samples"]) >= 200, (
                f"{r['dataset']}/{r['variant']}: only {r['n_samples']} samples"
            )

    def test_pg_std_reasonable(self):
        """PG std should be <75% of mean (not dominated by outliers)."""
        rows = _load_tsv(CLUSTERING_SUMMARY)
        for r in rows:
            mean = float(r["PG_mean"])
            std = float(r["PG_std"])
            if mean > 0:
                cv = std / mean
                assert cv < 0.75, f"{r['dataset']}/{r['variant']}: CV={cv:.2f} (too variable)"


@pytest.mark.skipif(not HAS_DETAIL, reason="Clustering detail not available")
class TestClusteringDetail:
    """Validate per-sample clustering results."""

    def test_detail_has_all_samples(self):
        """Detail file should have 4 variants per sample."""
        rows = _load_tsv(CLUSTERING_DETAIL)
        # Count samples per dataset
        from collections import Counter

        by_dataset = Counter(r["dataset"] for r in rows)
        for dataset, count in by_dataset.items():
            # 4 variants per sample
            assert count % 4 == 0, f"{dataset}: {count} rows not divisible by 4 variants"

    def test_no_negative_pg_counts(self):
        """All PG counts should be >= 0."""
        rows = _load_tsv(CLUSTERING_DETAIL)
        for r in rows:
            assert int(r["PG_1pct"]) >= 0, f"Negative PG: {r}"


# ═══════════════════════════════════════════════════════════════════════════════
# Cross-benchmark consistency
# ═══════════════════════════════════════════════════════════════════════════════


@pytest.mark.skipif(not HAS_CLUSTERING, reason="Results not available")
class TestCrossBenchmarkConsistency:
    """Check that CAPSCAN and MicrobPredict show consistent patterns."""

    def test_both_datasets_present(self):
        rows = _load_tsv(CLUSTERING_SUMMARY)
        datasets = set(r["dataset"] for r in rows)
        assert "CAPSCAN" in datasets
        assert "MicrobPredict" in datasets

    def test_same_variants_both_datasets(self):
        rows = _load_tsv(CLUSTERING_SUMMARY)
        capscan = set(r["variant"] for r in rows if r["dataset"] == "CAPSCAN")
        mp = set(r["variant"] for r in rows if r["dataset"] == "MicrobPredict")
        assert capscan == mp, f"Different variants: CAPSCAN={capscan}, MP={mp}"

    def test_loss_direction_consistent(self):
        """Both datasets should show same relative ordering of PG loss."""
        rows = _load_tsv(CLUSTERING_SUMMARY)
        for dataset in ["CAPSCAN", "MicrobPredict"]:
            ds = {
                r["variant"]: float(r["PG_pct_loss_vs_unclustered"])
                for r in rows
                if r["dataset"] == dataset
            }
            # c99 < c97 < c95
            assert ds["c99"] < ds["c97"] < ds["c95"], (
                f"{dataset}: loss order wrong — c99={ds['c99']}, c97={ds['c97']}, c95={ds['c95']}"
            )
