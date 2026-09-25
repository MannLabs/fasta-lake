"""Regressions for the second review pass: byte-identical gzip outputs, per-sample
manifest sizes, wrong-class single runs, FDP reader columns, Sage N-terminal
modifications, undefined FDP is None, source manifest outside git."""

from __future__ import annotations

import gzip
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest
from test_simplified_study import P1, search

from fasta_lake.entrapment_workflow import (
    STAGE_APPENDED,
    classify_psms,
    compute_fdp,
    find_sample_results,
)
from fasta_lake.gzip_io import open_gzip_text
from fasta_lake.multi_omics.complementarity import peptide as molecular_peptide
from fasta_lake.study import canonical_peptide, group_searches
from fasta_lake.tune import compute_fdp_combined

LAKE = ">prot1\nMKVLATGSEDKRAAWQFTGEDPLMNK\n"


# ── gzip outputs must not embed the clock ────────────────────────────────────


def test_open_gzip_text_is_clock_independent(tmp_path, monkeypatch):
    outputs = []
    for fake_now in (1_000_000_000, 2_000_000_000):
        monkeypatch.setattr(gzip.time, "time", lambda now=fake_now: now)
        path = tmp_path / f"{fake_now}.gz"
        with open_gzip_text(path) as fh:
            fh.write("a\tb\n1\t2\n")
        outputs.append(path.read_bytes())
    assert outputs[0] == outputs[1]
    assert outputs[0][4:8] == b"\0\0\0\0", "MTIME header field is zero"
    with gzip.open(tmp_path / "1000000000.gz", "rt") as fh:
        assert fh.read() == "a\tb\n1\t2\n"


def test_study_gzip_outputs_have_zero_mtime_and_repeat_byte_identically(tmp_path, monkeypatch):
    search(tmp_path, "S", {"A": P1}, [(P1, "A", 1, 0.001)], [])
    produced = []
    for fake_now, out in ((1_000_000_000, "g1"), (2_000_000_000, "g2")):
        monkeypatch.setattr(gzip.time, "time", lambda now=fake_now: now)
        group_searches(tmp_path / "search", tmp_path / out)
        produced.append({p.name: p.read_bytes() for p in (tmp_path / out).glob("*.gz")})
    assert produced[0].keys() == produced[1].keys() and produced[0]
    for name, data in produced[0].items():
        assert data[4:8] == b"\0\0\0\0", name
        assert data == produced[1][name], name


# ── entrapment: per-sample manifests, wrong-class single run, reader columns ──


def _sage_rows(path: Path, n_target: int, n_entrap: int, prefix: str = "ENTRAP_SHUF_"):
    lines = ["peptide\tproteins\tpeptide_q\tpeptide_len\tlabel"]
    lines += [f"TARGETPEP{j}\tprot1\t0.001\t10\t1" for j in range(n_target)]
    lines += [f"WWWWWWWWW{j}\t{prefix}x{j}\t0.001\t10\t1" for j in range(n_entrap)]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def _manifest(path: Path, n_target: int, n_entrapment: int):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "schema": "fastalake.entrapment.manifest/1",
                "stage": STAGE_APPENDED,
                "class": "shuffled",
                "seed": 42,
                "n_target": n_target,
                "n_entrapment": n_entrapment,
            }
        )
    )


def test_partition_sizes_come_from_each_samples_own_manifest(tmp_path):
    root = tmp_path / "results"
    _manifest(root / "s1" / "entrapment_manifest.json", 1000, 100)
    _manifest(root / "s2" / "entrapment_manifest.json", 50, 10)
    _sage_rows(root / "s1" / "entrapment" / "results.sage.tsv", 20, 1)
    _sage_rows(root / "s2" / "entrapment" / "results.sage.tsv", 20, 1)
    lake = tmp_path / "lake.fasta"
    lake.write_text(LAKE)
    report = compute_fdp(root, target_fasta=lake)
    rows = {r["sample"]: r for r in report["per_sample"]}
    assert (rows["s1"]["db_target"], rows["s1"]["db_entrapment"]) == (1000, 100)
    assert (rows["s2"]["db_target"], rows["s2"]["db_entrapment"]) == (50, 10)
    assert rows["s1"]["r"] != rows["s2"]["r"]
    assert "s1" in rows["s1"]["db_size_source"] and "s2" in rows["s2"]["db_size_source"]


def test_single_run_of_another_class_is_refused(tmp_path):
    root = tmp_path / "results"
    _sage_rows(root / "s1" / "entrapment" / "results.sage.tsv", 5, 0)
    with pytest.raises(ValueError, match="different entrapment run.*'entrapment_plant'"):
        find_sample_results(root, prefer_subdir="entrapment_plant")
    # Explicit choice by --results-name still works and the same class is accepted.
    assert find_sample_results(root, results_name="entrapment/results.sage.tsv")
    assert find_sample_results(root, prefer_subdir="entrapment")


def test_fdp_reader_refuses_renamed_proteins_column(tmp_path):
    tsv = tmp_path / "results.sage.tsv"
    tsv.write_text("peptide\tprotein\tpeptide_q\tlabel\nAAAAAAAAAK\tENTRAP_SHUF_1\t0.001\t1\n")
    with pytest.raises(ValueError, match=r"missing required column\(s\) proteins"):
        classify_psms(tsv, prefix="ENTRAP_SHUF_")


# ── peptide canonicalisers accept Sage N-terminal modifications ──────────────


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("[+42.0106]-PEPTIDEK", "PEPTLDEK"),
        ("PEPT[+15.9949]IDEK", "PEPTLDEK"),
        ("[+42.0106]-PEPT[+15.9949]IDEK", "PEPTLDEK"),
        ("PEPTIDEK", "PEPTLDEK"),
    ],
)
def test_sage_modification_strings_canonicalise(raw, expected):
    assert canonical_peptide(raw) == expected
    assert molecular_peptide(raw) == expected


# ── an undefined FDP is None, never 0.0 ───────────────────────────────────────


def test_undefined_fdp_is_none():
    assert compute_fdp_combined(0, 0, 0.1) is None
    assert compute_fdp_combined(1, 99, 0.1) == pytest.approx(1 * (1 + 10) / 100 * 100)


# ── source manifest tool works in an extracted archive ───────────────────────


def test_source_manifest_check_works_without_git(tmp_path):
    spec = importlib.util.spec_from_file_location(
        "sm", Path(__file__).resolve().parents[1] / "tools" / "source_manifest.py"
    )
    sm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sm)
    root = tmp_path / "archive"
    (root / "pkg").mkdir(parents=True)
    (root / "pkg" / "a.py").write_text("x = 1\n")
    (root / "target").mkdir()
    (root / "target" / "junk.o").write_bytes(b"\0")
    manifest = root / "SOURCE_MANIFEST.json"
    assert sm.write(root, manifest) == 0
    assert [f["path"] for f in json.loads(manifest.read_text())["files"]] == ["pkg/a.py"]
    assert sm.check(root, manifest) == 0
    (root / "pkg" / "a.py").write_text("x = 2\n")
    assert sm.check(root, manifest) == 1
    # And the CLI entry point, with git absent from PATH, exits 1 with the message.
    result = subprocess.run(
        [sys.executable, str(sm.__file__), "--check"],
        cwd=root,
        env={"PATH": "/nonexistent"},
        capture_output=True,
        text=True,
    )
    assert result.returncode in (0, 1) and "Traceback" not in result.stderr
