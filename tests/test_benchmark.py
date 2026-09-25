"""Tests for the benchmark comparison module."""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from fasta_lake.benchmark import (
    SearchResult,
    _pearson,
    compare,
    load_diann_report,
    load_sage_report,
)


def test_compare_basic():
    """Basic comparison: FastaLake finds more than baseline."""
    fl = SearchResult(
        peptides={"A", "B", "C", "D", "E"},
        genes={"G1", "G2", "G3"},
        protein_groups={"P1", "P2", "P3"},
        total_precursors=100,
    )
    bl = SearchResult(
        peptides={"A", "B", "C"},
        genes={"G1", "G2"},
        protein_groups={"P1", "P2"},
        total_precursors=80,
    )

    bench = compare(fl, bl)
    assert bench.fl_peptides == 5
    assert bench.bl_peptides == 3
    assert bench.peptide_gain_pct > 0  # +66.7%
    assert bench.shared_peptides == 3
    assert bench.fl_only_peptides == 2
    assert bench.bl_only_peptides == 0
    assert bench.gene_coverage_pct == 100.0  # all baseline genes retained


def test_compare_lost_genes():
    """FastaLake loses some baseline genes."""
    fl = SearchResult(
        peptides={"A", "B", "D", "E"},
        genes={"G1", "G3", "G4"},
        protein_groups=set(),
    )
    bl = SearchResult(
        peptides={"A", "B", "C"},
        genes={"G1", "G2"},
        protein_groups=set(),
    )

    bench = compare(fl, bl)
    assert bench.shared_genes == 1  # G1
    assert bench.fl_only_genes == 2  # G3, G4
    assert bench.bl_only_genes == 1  # G2 lost
    assert bench.gene_coverage_pct == 50.0  # only 1/2 baseline genes retained
    assert any("coverage" in w.lower() for w in bench.warnings)


def test_compare_broken_genes():
    """Broken gene names should be flagged."""
    fl = SearchResult(
        peptides={"A"},
        genes={"sp|P02768|ALBU_HUMAN"},
        protein_groups=set(),
        broken_genes=5,
    )
    bl = SearchResult(
        peptides={"A"},
        genes={"ALB"},
        protein_groups=set(),
    )

    bench = compare(fl, bl)
    assert bench.broken_genes == 5
    assert any("broken" in w.lower() for w in bench.warnings)


def test_compare_entrapment_fdr():
    """Entrapment FDR calculation."""
    fl = SearchResult(
        peptides=set(),
        genes=set(),
        protein_groups=set(),
        total_precursors=1000,
        entrapment_peptides=60,
    )
    bl = SearchResult(
        peptides=set(),
        genes=set(),
        protein_groups=set(),
        total_precursors=800,
        entrapment_peptides=45,
    )

    bench = compare(fl, bl)
    assert abs(bench.fl_entrapment_precursor_fdr - 6.0) < 0.1
    assert abs(bench.bl_entrapment_precursor_fdr - 5.625) < 0.1


def test_intensity_correlation():
    """Intensity correlation for shared peptides."""
    shared_peps = {f"PEP{i}" for i in range(20)}
    fl = SearchResult(
        peptides=shared_peps | {"EXTRA1"},
        genes=set(),
        protein_groups=set(),
        peptide_intensities={f"PEP{i}": 10 ** (i * 0.5) for i in range(20)},
    )
    bl = SearchResult(
        peptides=shared_peps,
        genes=set(),
        protein_groups=set(),
        peptide_intensities={
            f"PEP{i}": 10 ** (i * 0.5) * 1.1 for i in range(20)
        },  # slightly scaled
    )

    bench = compare(fl, bl)
    assert bench.n_shared_for_correlation == 20
    assert bench.intensity_correlation > 0.99  # should be near-perfect


def test_pearson():
    """Pearson correlation."""
    assert abs(_pearson([1, 2, 3, 4, 5], [1, 2, 3, 4, 5]) - 1.0) < 0.001
    assert abs(_pearson([1, 2, 3, 4, 5], [5, 4, 3, 2, 1]) - (-1.0)) < 0.001
    assert abs(_pearson([1, 1, 1], [1, 1, 1]) - 0.0) < 0.001  # zero variance


def test_load_diann_report():
    """Load a DIA-NN report TSV."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("File.Name\tStripped.Sequence\tGenes\tProtein.Group\tQ.Value\tPrecursor.Quantity\n")
        f.write("sample.mzML\tAAVLIDK\tALB\tP02768\t0.001\t1000000\n")
        f.write("sample.mzML\tDVFLGMFLYEYAR\tALB\tP02768\t0.005\t500000\n")
        f.write("sample.mzML\tBADPEPTIDE\tENTRAP_001\tENTRAP_001\t0.009\t1000\n")
        f.write("sample.mzML\tFILTERED\tX\tX\t0.05\t100\n")  # above qvalue
        path = f.name

    try:
        result = load_diann_report(path, qvalue=0.01)
        assert len(result.peptides) == 3  # AAVLIDK, DVFLGMFLYEYAR, BADPEPTIDE
        assert "ALB" in result.genes
        assert result.total_precursors == 3
        assert result.entrapment_proteins == 1
    finally:
        os.unlink(path)


def test_load_sage_report_counts_entrapment(tmp_path):
    """A SAGE results file with entrapment hits must not report 0 % entrapment FDR."""
    path = tmp_path / "results.sage.tsv"
    path.write_text(
        "peptide\tproteins\tlabel\tspectrum_q\tpeptide_q\tms1_intensity\n"
        "AAVLIDK\tP02768\t1\t0.001\t0.001\t1000\n"
        "DVFLGM[+15.9949]FLYEYAR\tP02768;P02769\t1\t0.005\t0.005\t500\n"
        "BADPEPTIDE\tENTRAP_SHUF_000123\t1\t0.009\t0.009\t10\n"
        "SHAREDPEP\tP02768;ENTRAP_PLANT_7\t1\t0.002\t0.002\t20\n"
        "DECOYPEP\trev_P02768\t-1\t0.001\t0.001\t5\n"
        "FILTERED\tENTRAP_SHUF_9\t1\t0.05\t0.05\t1\n"
    )
    result = load_sage_report(path, qvalue=0.01)
    assert result.total_precursors == 4  # decoy and q>0.01 rows excluded
    assert result.peptides == {"AAVLIDK", "DVFLGMFLYEYAR", "BADPEPTIDE", "SHAREDPEP"}
    assert result.entrapment_peptides == 2  # pure entrapment + shared, not filtered/decoy
    bench = compare(result, SearchResult())
    assert abs(bench.fl_entrapment_precursor_fdr - 50.0) < 1e-9


def test_compare_empty():
    """Handle empty results gracefully."""
    fl = SearchResult()
    bl = SearchResult()
    bench = compare(fl, bl)
    assert bench.fl_peptides == 0
    assert bench.gene_coverage_pct == 0.0
    assert bench.intensity_correlation == 0.0


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

    print(f"\n{'=' * 40}")
    print(f"  {passed} passed, {failed} failed")
    print(f"{'=' * 40}")
