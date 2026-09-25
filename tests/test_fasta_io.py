"""Tests for header-preserving FASTA I/O."""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from fasta_lake.headers.fasta_io import (
    load_fasta_full,
    load_fasta_proteins,
    write_fasta,
    write_fasta_simple,
)


def _make_test_fasta():
    """Create a test FASTA with multiple header formats."""
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
    f.write(">sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2\n")
    f.write("MKWVTFISLLFLFSSAYS\n")
    f.write(">MGYG000000001_00001 Ribosomal protein S1\n")
    f.write("MKVLWAALLVTFLAGCQA\n")
    f.write(">var-nfe|P02768|ALB_Glu177Asp|4-73409403 ALB Albumin variant GN=ALB\n")
    f.write("MKWVTFISLLFLFSSAYD\n")
    f.close()
    return f.name


def test_load_fasta_full():
    """Full loader preserves headers."""
    path = _make_test_fasta()
    try:
        proteins = load_fasta_full(path)
        assert len(proteins) == 3

        # Check that full header is preserved
        seq, header = proteins["sp|P02768|ALBU_HUMAN"]
        assert "GN=ALB" in header
        assert "Albumin" in header
        assert seq == "MKWVTFISLLFLFSSAYS"

        # Check MGYG
        seq, header = proteins["MGYG000000001_00001"]
        assert "Ribosomal protein S1" in header

        # Check variant
        seq, header = proteins["var-nfe|P02768|ALB_Glu177Asp|4-73409403"]
        assert "ALB" in header
    finally:
        os.unlink(path)


def test_load_fasta_proteins_backward_compat():
    """Backward-compatible loader returns {id: sequence} only."""
    path = _make_test_fasta()
    try:
        proteins = load_fasta_proteins(path)
        assert len(proteins) == 3
        assert isinstance(proteins["sp|P02768|ALBU_HUMAN"], str)
        assert proteins["sp|P02768|ALBU_HUMAN"] == "MKWVTFISLLFLFSSAYS"
    finally:
        os.unlink(path)


def test_write_fasta_preserves_headers():
    """Write FASTA with full headers."""
    path = _make_test_fasta()
    try:
        proteins = load_fasta_full(path)
        out = path + ".out"
        count = write_fasta(proteins, out)
        assert count == 3

        # Read back and verify headers are preserved
        with open(out) as f:
            headers = [line.rstrip() for line in f if line.startswith(">")]

        assert any("GN=ALB" in h for h in headers)
        assert any("MGYG" in h for h in headers)
        assert any("var-nfe" in h for h in headers)
        os.unlink(out)
    finally:
        os.unlink(path)


def test_write_fasta_simple_subset():
    """Write a subset of proteins."""
    path = _make_test_fasta()
    try:
        proteins = load_fasta_full(path)
        out = path + ".subset"
        count = write_fasta_simple(
            ["sp|P02768|ALBU_HUMAN", "MGYG000000001_00001"],
            proteins,
            out,
        )
        assert count == 2

        # Verify only 2 entries
        with open(out) as f:
            n = sum(1 for line in f if line.startswith(">"))
        assert n == 2
        os.unlink(out)
    finally:
        os.unlink(path)


def test_roundtrip_header_preservation():
    """Load → write → load must preserve all headers exactly."""
    path = _make_test_fasta()
    try:
        proteins1 = load_fasta_full(path)
        out = path + ".roundtrip"
        write_fasta(proteins1, out)
        proteins2 = load_fasta_full(out)

        for pid in proteins1:
            seq1, hdr1 = proteins1[pid]
            seq2, hdr2 = proteins2[pid]
            assert seq1 == seq2, f"Sequence mismatch for {pid}"
            assert hdr1 == hdr2, f"Header mismatch for {pid}: {hdr1!r} != {hdr2!r}"

        os.unlink(out)
    finally:
        os.unlink(path)


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


class TestLoadFastaFullRefusesDefectiveInput:
    def test_duplicate_and_empty_are_refused(self, tmp_path):
        import pytest

        from fasta_lake.headers.fasta_io import load_fasta_full

        p = tmp_path / "dup.fasta"
        p.write_text(">X a\nMKV\n>X b\nMLL\n")
        with pytest.raises(ValueError, match="duplicate accession 'X'"):
            load_fasta_full(p)
        p.write_text(">X a\n>Y b\nMLL\n")
        with pytest.raises(ValueError, match="record 'X'.*has no sequence"):
            load_fasta_full(p)
        p.write_text(">X a\nMKV\n>Y b\nMLL\n")
        assert load_fasta_full(p) == {"X": ("MKV", "X a"), "Y": ("MLL", "Y b")}
