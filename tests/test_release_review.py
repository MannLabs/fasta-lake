"""Adversarial release regressions: file preservation and scientific accounting."""

from __future__ import annotations

import csv
import importlib.util
import subprocess
from pathlib import Path

import pytest

from fasta_lake.entrapment import estimate_fdp
from fasta_lake.rust.binaries import find_binary, run_fasta_extractor


@pytest.fixture
def extractor():
    try:
        return find_binary("fasta_extractor")
    except FileNotFoundError:
        pytest.skip("Build the Rust binaries to run the release boundary checks")


@pytest.fixture
def inputs(tmp_path):
    db = tmp_path / "db.fasta"
    csv_path = tmp_path / "pred.csv"
    db.write_text(">P1\nMPEPTIDEAKOTHERPEPK\n")
    csv_path.write_text("sequence,score\nMPEPTIDEAK,0.9\n")
    return db, csv_path


def run_extract(extractor, inputs, output, *extra):
    db, predictions = inputs
    return subprocess.run(
        [
            str(extractor),
            "--database",
            str(db),
            "--peptides-dir",
            str(predictions),
            "--output",
            str(output),
            "-t",
            "1",
            *extra,
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )


@pytest.mark.parametrize(
    "extra",
    [
        ("--top-percent", "30"),
        ("--top-percent", "NaN"),
        ("--chunk-size", "0"),
        ("--min-length", "0"),
    ],
)
def test_invalid_options_fail_without_outputs(extractor, inputs, tmp_path, extra):
    output = tmp_path / "result.fasta"
    result = run_extract(extractor, inputs, output, *extra)
    assert result.returncode != 0
    assert not output.exists()
    assert not output.with_suffix(".stats.json").exists()


def test_stats_sidecar_is_returned(extractor, inputs, tmp_path):
    stats = run_fasta_extractor(*inputs, tmp_path / "result.fasta", binary_path=extractor)
    assert stats["proteins_matched"] == 1


def test_existing_output_is_preserved(extractor, inputs, tmp_path):
    output = tmp_path / "result.fasta"
    output.write_bytes(b"original result must survive\n")
    result = run_extract(extractor, inputs, output)
    assert result.returncode != 0
    assert output.read_bytes() == b"original result must survive\n"


def test_missing_include_is_fatal(extractor, inputs, tmp_path):
    output = tmp_path / "result.fasta"
    result = run_extract(extractor, inputs, output, "--include", str(tmp_path / "missing.fa"))
    assert result.returncode != 0
    assert not output.exists()


def test_empty_entrapment_cannot_pass(extractor, inputs, tmp_path):
    entrapment = tmp_path / "empty.fasta"
    entrapment.write_text(">empty\n")
    output = tmp_path / "result.fasta"
    result = run_extract(extractor, inputs, output, "--entrapment", str(entrapment))
    assert result.returncode != 0
    assert not output.with_suffix(".stats.json").exists()


def test_sample_breadth_unions_disjoint_peptide_evidence(extractor, inputs, tmp_path):
    db, _ = inputs
    pool = tmp_path / "pool"
    pool.mkdir()
    (pool / "A.csv").write_text("sequence,score\nMPEPTIDEAK,0.9\nMPEPTIDEAK,0.8\n")
    (pool / "B.csv").write_text("sequence,score\nOTHERPEPK,0.9\n")
    evidence = tmp_path / "evidence.tsv"
    result = run_extract(
        extractor,
        (db, pool),
        tmp_path / "result.fasta",
        "--min-length",
        "8",
        "--output-evidence",
        str(evidence),
    )
    assert result.returncode == 0, result.stderr
    with evidence.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 1
    assert int(rows[0]["n_samples"]) == 2


def test_fdp_rejects_negative_counts_and_marks_empty_target_undefined():
    with pytest.raises(ValueError, match="nonnegative"):
        estimate_fdp(-1, 10, 100, 100)
    result = estimate_fdp(0, 10, 0, 100)
    assert result["evaluable"] is False
    assert result["fdp_combined"] is None


def test_path_discovery_works_outside_checkout(monkeypatch, tmp_path):
    binary = tmp_path / "fasta_extractor"
    binary.write_text("#!/bin/sh\nexit 0\n")
    binary.chmod(0o755)
    monkeypatch.delenv("FASTALAKE_BIN_DIR", raising=False)
    monkeypatch.setenv("PATH", str(tmp_path))
    assert find_binary("fasta_extractor", search_dir=tmp_path / "no_checkout") == binary


def load_runner():
    path = Path(__file__).resolve().parents[1] / "tools/run_manifest.py"
    spec = importlib.util.spec_from_file_location("reviewed_runner", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_manifest_rejects_duplicate_ids_and_resolves_relative_paths(tmp_path):
    runner = load_runner()
    (tmp_path / "S.csv").write_text("sequence,score\nMPEPTIDEAK,0.9\n")
    (tmp_path / "S.mzML").write_text("<mzML/>")
    manifest = tmp_path / "samples.tsv"
    header = "sample\tpredictions\tmzml\n"
    row = "S\tS.csv\tS.mzML\n"
    manifest.write_text(header + row)
    assert runner.load_manifest(manifest)[0]["predictions"] == str(tmp_path / "S.csv")
    manifest.write_text(header + row + row)
    with pytest.raises(ValueError, match="duplicate"):
        runner.load_manifest(manifest)


def test_search_config_exposes_stage_specific_length_difference():
    config = load_runner().sage_config("db.fasta", "S.mzML", "out")
    assert config["database"]["enzyme"]["min_len"] == 5
    assert config["database"]["generate_decoys"] is True
    assert config["quant"]["lfq"] is True


@pytest.mark.parametrize("invalid_score", ["NaN", "inf", "broken"])
def test_parsimony_does_not_promote_invalid_scores(tmp_path, invalid_score):
    try:
        binary = find_binary("parsimony_engine")
    except FileNotFoundError:
        pytest.skip("Build parsimony_engine for the malformed-score regression")
    fasta = tmp_path / "source.fasta"
    fasta.write_text(">P1\nMPEPTIDEAK\n>P2\nOTHERPEPKK\n")
    predictions = tmp_path / "pred.csv"
    predictions.write_text("sequence,score\nMPEPTIDEAK,0.9\nOTHERPEPKK," + invalid_score + "\n")
    output = tmp_path / "inference"
    result = subprocess.run(
        [
            str(binary),
            "--fasta",
            str(fasta),
            "--peptides",
            str(predictions),
            "--sample",
            "S",
            "--output-dir",
            str(output),
            "--strategy",
            "razor",
            "--tiebreak",
            "hash-acc",
            "--min-length",
            "8",
            "-t",
            "1",
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr
    selected = (output / "S_razor.fasta").read_text()
    assert ">P1" in selected
    assert ">P2" not in selected
