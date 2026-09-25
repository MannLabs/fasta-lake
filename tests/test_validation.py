"""Tests for the pipeline validation module."""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from fasta_lake.headers.validation import (
    validate_denovo_csv,
    validate_evidence_lake,
    validate_inference_output,
    validate_source_fasta,
)

# ═══════════════════════════════════════════════════════════════
# DE NOVO CSV VALIDATION
# ═══════════════════════════════════════════════════════════════


def test_validate_denovo_missing_file():
    results = validate_denovo_csv("/nonexistent/file.csv")
    assert any(not r.passed and r.severity == "error" for r in results)


def test_validate_denovo_alphanovo():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("peptide_prediction_detokenized_unmodified,score\n")
        for i in range(200):
            f.write(f"AAVLIDKTNVKAEFAEVSK,{0.5 + i * 0.001}\n")
        path = f.name
    try:
        results = validate_denovo_csv(path)
        engines = [r for r in results if r.check == "engine_detected"]
        assert len(engines) == 1
        assert engines[0].passed
        assert engines[0].details["engine"] == "alphanovo"
    finally:
        os.unlink(path)


def test_validate_denovo_dianovo():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("graph_idx,pred_path,pred_prob,pred_seq\n")
        for i in range(200):
            f.write(f'{i},"path",0.5 0.3 0.2,A V L I D K\n')
        path = f.name
    try:
        results = validate_denovo_csv(path)
        engines = [r for r in results if r.check == "engine_detected"]
        assert len(engines) == 1
        assert engines[0].passed
        assert engines[0].details["engine"] == "dianovo"
    finally:
        os.unlink(path)


def test_validate_denovo_unknown_format():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("weird_column,other_column\n")
        f.write("data,more_data\n")
        path = f.name
    try:
        results = validate_denovo_csv(path)
        engines = [r for r in results if r.check == "engine_detected"]
        assert len(engines) == 1
        assert not engines[0].passed
    finally:
        os.unlink(path)


def test_validate_denovo_windows_line_endings():
    with tempfile.NamedTemporaryFile(mode="wb", suffix=".csv", delete=False) as f:
        f.write(b"peptide_prediction_detokenized_unmodified,score\r\n")
        f.write(b"AAVLIDKTNVK,0.9\r\n")
        path = f.name
    try:
        results = validate_denovo_csv(path)
        crlf = [r for r in results if r.check == "line_endings"]
        assert len(crlf) == 1
        assert not crlf[0].passed  # should warn
    finally:
        os.unlink(path)


# ═══════════════════════════════════════════════════════════════
# SOURCE FASTA VALIDATION
# ═══════════════════════════════════════════════════════════════


def test_validate_fasta_missing():
    results = validate_source_fasta("/nonexistent.fasta")
    assert any(not r.passed and r.severity == "error" for r in results)


def test_validate_fasta_normal():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        for i in range(50):
            f.write(f">sp|P{i:05d}|PROT_{i} Description GN=GENE{i}\n")
            f.write("MKVLWAALLVTFLAGCQA" * 5 + "\n")
        path = f.name
    try:
        results = validate_source_fasta(path)
        counts = [r for r in results if r.check == "entry_count"]
        assert counts[0].passed
        assert counts[0].details["n_entries"] == 50
    finally:
        os.unlink(path)


def test_validate_fasta_empty():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write("")
        path = f.name
    try:
        results = validate_source_fasta(path)
        counts = [r for r in results if r.check == "entry_count"]
        assert not counts[0].passed
    finally:
        os.unlink(path)


def test_validate_fasta_nucleotide():
    """Should flag nucleotide sequences."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        for i in range(100):
            f.write(f">seq_{i}\n")
            f.write("ATCGATCGATCGATCGATCG\n")
        path = f.name
    try:
        results = validate_source_fasta(path)
        nuc = [r for r in results if r.check == "nucleotide_check"]
        assert len(nuc) == 1
        assert not nuc[0].passed
    finally:
        os.unlink(path)


# ═══════════════════════════════════════════════════════════════
# EVIDENCE LAKE VALIDATION
# ═══════════════════════════════════════════════════════════════


def test_validate_evidence_lake():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        for i in range(100):
            f.write(f">MGYG000000001_{i:05d} protein {i}\n")
            f.write("MKVLWAALLVTFLAGCQA\n")
        path = f.name
    try:
        results = validate_evidence_lake(path)
        counts = [r for r in results if r.check == "entry_count"]
        assert counts[0].passed
        assert counts[0].details["n_entries"] == 100
    finally:
        os.unlink(path)


# ═══════════════════════════════════════════════════════════════
# INFERENCE OUTPUT VALIDATION
# ═══════════════════════════════════════════════════════════════


def test_validate_inference_truncated_headers():
    """Detect Stage 3 header truncation bug."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        # Truncated headers (only protein ID, no description)
        f.write(">sp|P02768|ALBU_HUMAN\nMKWVTF\n")
        f.write(">MGYG000000001_00001\nMKWVTF\n")
        path = f.name
    try:
        results = validate_inference_output(path)
        headers = [r for r in results if r.check == "header_preservation"]
        assert len(headers) == 1
        assert not headers[0].passed  # should warn about truncation
        assert headers[0].details["truncated"] == 2
    finally:
        os.unlink(path)


def test_validate_inference_full_headers():
    """Full headers should pass."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(">sp|P02768|ALBU_HUMAN Albumin GN=ALB\nMKWVTF\n")
        f.write(">MGYG000000001_00001 hypothetical protein\nMKWVTF\n")
        path = f.name
    try:
        results = validate_inference_output(path)
        headers = [r for r in results if r.check == "header_preservation"]
        assert headers[0].passed
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
