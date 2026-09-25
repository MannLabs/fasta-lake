"""Tests for cross-sample aggregation.

Key principle: intensities in the aggregated matrix must EXACTLY match
the per-sample search results they came from. No silent value changes.
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from fasta_lake.aggregate import (
    SampleResult,
    build_quant_matrix,
    compute_cross_sample_stats,
    load_diann_sample,
)


def _make_sample(sample_id, peptides, genes):
    """Helper to create a SampleResult."""
    return SampleResult(
        sample_id=sample_id,
        peptides=peptides,
        genes=genes,
        n_precursors=sum(1 for _ in peptides),
    )


def test_build_quant_matrix_genes():
    """Gene-level matrix has correct dimensions."""
    s1 = _make_sample("S1", {"A": 100, "B": 200}, {"G1": 300, "G2": 100})
    s2 = _make_sample("S2", {"A": 150, "C": 50}, {"G1": 200, "G3": 50})

    sample_ids, feature_ids, matrix = build_quant_matrix([s1, s2], level="gene")

    assert sample_ids == ["S1", "S2"]
    assert len(feature_ids) == 3  # G1, G2, G3
    assert len(matrix) == 2
    assert len(matrix[0]) == 3


def test_quant_matrix_intensities_match_source():
    """CRITICAL: matrix values must exactly match source sample values."""
    s1 = _make_sample("S1", {}, {"ALB": 1e9, "HP": 5e7, "APOA1": 2e8})
    s2 = _make_sample("S2", {}, {"ALB": 8e8, "HP": 6e7, "GC": 1e6})

    sample_ids, feature_ids, matrix = build_quant_matrix([s1, s2], level="gene")

    # Find ALB column
    alb_idx = feature_ids.index("ALB")
    hp_idx = feature_ids.index("HP")
    gc_idx = feature_ids.index("GC")

    # S1 values must match exactly
    assert matrix[0][alb_idx] == 1e9
    assert matrix[0][hp_idx] == 5e7

    # S2 values must match exactly
    assert matrix[1][alb_idx] == 8e8
    assert matrix[1][gc_idx] == 1e6

    # Missing values must be 0
    apoa1_idx = feature_ids.index("APOA1")
    assert matrix[1][apoa1_idx] == 0.0  # S2 doesn't have APOA1


def test_quant_matrix_no_silent_value_changes():
    """No rounding, normalization, or transformation applied silently."""
    precise_value = 123456789.123456
    s1 = _make_sample("S1", {}, {"GENE": precise_value})

    _, feature_ids, matrix = build_quant_matrix([s1], level="gene")
    gene_idx = feature_ids.index("GENE")

    # Must be BIT-IDENTICAL
    assert matrix[0][gene_idx] == precise_value


def test_quant_matrix_peptide_level():
    """Peptide-level matrix works."""
    s1 = _make_sample("S1", {"AAVLIDK": 1000, "DVFLGMFLYEYAR": 500}, {})
    s2 = _make_sample("S2", {"AAVLIDK": 1200, "NQVALNPQNTVFDAK": 300}, {})

    _, feature_ids, matrix = build_quant_matrix([s1, s2], level="peptide")

    assert "AAVLIDK" in feature_ids
    assert "DVFLGMFLYEYAR" in feature_ids
    assert "NQVALNPQNTVFDAK" in feature_ids

    aav_idx = feature_ids.index("AAVLIDK")
    assert matrix[0][aav_idx] == 1000
    assert matrix[1][aav_idx] == 1200


def test_cross_sample_stats():
    """Cross-sample statistics are computed correctly."""
    s1 = _make_sample("S1", {}, {"G1": 100, "G2": 200, "G3": 300})
    s2 = _make_sample("S2", {}, {"G1": 110, "G2": 190, "G4": 50})
    s3 = _make_sample("S3", {}, {"G1": 105, "G2": 205, "G3": 310, "G4": 55})

    stats = compute_cross_sample_stats([s1, s2, s3], "gene")

    assert stats["n_samples"] == 3
    assert stats["n_features"] == 4  # G1, G2, G3, G4
    assert stats["features_in_all_samples"] == 2  # G1, G2
    assert stats["features_in_one_sample"] == 0  # G3 in 2 samples, G4 in 2 samples
    assert stats["completeness_pct"] > 0


def test_cross_sample_stats_single_sample():
    """Single sample: no CV possible."""
    s1 = _make_sample("S1", {}, {"G1": 100})
    stats = compute_cross_sample_stats([s1], "gene")
    assert stats["n_samples"] == 1
    assert stats["n_features"] == 1
    assert stats["median_cv_pct"] == 0  # can't compute CV with 1 sample


def test_load_diann_roundtrip():
    """Load a DIA-NN report and verify intensities match."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("File.Name\tStripped.Sequence\tGenes\tProtein.Group\tQ.Value\tPrecursor.Quantity\n")
        f.write("s1.mzML\tAAVLIDK\tALB\tP02768\t0.001\t1000000\n")
        f.write("s1.mzML\tDVFLGMFLYEYAR\tALB\tP02768\t0.005\t500000\n")
        f.write("s1.mzML\tVTSIQDWVQK\tHP\tP00738\t0.003\t200000\n")
        # Duplicate peptide with higher intensity — should keep max
        f.write("s1.mzML\tAAVLIDK\tALB\tP02768\t0.002\t1500000\n")
        path = f.name

    try:
        result = load_diann_sample(path, "S1")

        # Peptide intensities: AAVLIDK should be 1500000 (max of two rows)
        assert result.peptides["AAVLIDK"] == 1500000
        assert result.peptides["DVFLGMFLYEYAR"] == 500000
        assert result.peptides["VTSIQDWVQK"] == 200000

        # Gene intensities: ALB = sum of all ALB precursors
        assert result.genes["ALB"] == 1000000 + 500000 + 1500000  # all ALB rows summed
        assert result.genes["HP"] == 200000
    finally:
        os.unlink(path)


def test_load_diann_qvalue_filter():
    """Only rows passing q-value filter are included."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("File.Name\tStripped.Sequence\tGenes\tProtein.Group\tQ.Value\tPrecursor.Quantity\n")
        f.write("s1.mzML\tGOODPEP\tGENE1\tP1\t0.005\t1000\n")
        f.write("s1.mzML\tBADPEP\tGENE2\tP2\t0.05\t2000\n")  # above 1% FDR
        path = f.name

    try:
        result = load_diann_sample(path, "S1", qvalue=0.01)
        assert "GOODPEP" in result.peptides
        assert "BADPEP" not in result.peptides
        assert result.n_precursors == 1
    finally:
        os.unlink(path)


def test_multiple_samples_aggregation():
    """Build matrix from multiple DIA-NN reports."""
    paths = []
    for i, (peps, genes) in enumerate(
        [
            ({"A": 100, "B": 200}, {"G1": 300}),
            ({"A": 150, "C": 50}, {"G1": 200, "G2": 50}),
            ({"A": 120, "B": 180, "C": 60}, {"G1": 300, "G2": 60}),
        ]
    ):
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
        f.write("File.Name\tStripped.Sequence\tGenes\tProtein.Group\tQ.Value\tPrecursor.Quantity\n")
        for pep, intensity in peps.items():
            gene = list(genes.keys())[0]
            f.write(f"s.mzML\t{pep}\t{gene}\tPG\t0.001\t{intensity}\n")
        f.close()
        paths.append(f.name)

    try:
        samples = [load_diann_sample(p, f"S{i}") for i, p in enumerate(paths)]
        _, feature_ids, matrix = build_quant_matrix(samples, level="peptide")

        # Every sample should have peptide A
        a_idx = feature_ids.index("A")
        assert matrix[0][a_idx] == 100
        assert matrix[1][a_idx] == 150
        assert matrix[2][a_idx] == 120

        # B is in S1 and S3 only
        b_idx = feature_ids.index("B")
        assert matrix[0][b_idx] == 200
        assert matrix[1][b_idx] == 0  # missing
        assert matrix[2][b_idx] == 180
    finally:
        for p in paths:
            os.unlink(p)


if __name__ == "__main__":
    test_functions = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    failed = 0
    for test in test_functions:
        try:
            test()
            print(f"  PASS  {test.__name__}")
            passed += 1
        except Exception as e:
            print(f"  FAIL  {test.__name__}: {e}")
            failed += 1
    print(f"\n  {passed} passed, {failed} failed")
