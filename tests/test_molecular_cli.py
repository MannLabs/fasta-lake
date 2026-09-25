"""Installed molecular commands preserve the validated evidence contract."""

import gzip
import json

import pytest
from click.testing import CliRunner
from test_molecular_exact import searches, source

from fasta_lake.cli import cli
from fasta_lake.multi_omics.exact import read_fasta


def test_cli_union_and_specimen_rejection(tmp_path):
    manifest = source(tmp_path)
    args = [
        "molecular",
        "add",
        "--baseline",
        str(tmp_path / "baseline.fasta"),
        "--source",
        str(manifest),
        "--specimen",
        "s1",
        "--out",
        str(tmp_path / "out"),
        "--dna-percent",
        "50",
        "--rna-absolute",
        "1",
    ]
    runner = CliRunner()
    wrong = args.copy()
    wrong[wrong.index("s1")] = "foreign"
    failed = runner.invoke(cli, wrong)
    assert failed.exit_code != 0 and "specimen mismatch" in failed.output
    assert not (tmp_path / "out").exists()
    result = runner.invoke(cli, args)
    assert result.exit_code == 0, result.output
    assert json.loads(result.output)["additions"] == 2
    assert {r[1] for r in read_fasta(tmp_path / "out/search.fasta").values()} == {
        "MQQQQK",
        "MACDEK",
        "MLLK",
    }
    before = (tmp_path / "out/search.fasta").read_bytes()
    assert runner.invoke(cli, args).exit_code != 0
    assert (tmp_path / "out/search.fasta").read_bytes() == before


def test_cli_source_integrity(tmp_path):
    manifest = source(tmp_path)
    args = ["molecular", "validate-source", "--source", str(manifest)]
    result = CliRunner().invoke(cli, args)
    assert result.exit_code == 0, result.output
    assert json.loads(result.output)["genes"] == 4
    (tmp_path / "tpm.tsv").write_text("changed input")
    failed = CliRunner().invoke(cli, args)
    assert failed.exit_code != 0 and "checksum mismatch" in failed.output


def test_cli_comparison_labels_and_compressed_evidence(tmp_path):
    left, right = searches(tmp_path)
    result = CliRunner().invoke(
        cli,
        [
            "molecular",
            "compare",
            "--left",
            str(left),
            "--right",
            str(right),
            "--out",
            str(tmp_path / "comparison"),
            "--left-label",
            "baseline",
            "--right-label",
            "union",
            "--compress-tables",
        ],
    )
    assert result.exit_code == 0, result.output
    summary = json.loads(result.output)
    assert summary["protein_cells"] == {
        "both": 2,
        "baseline_only": 0,
        "union_only": 1,
        "neither": 0,
    }
    with gzip.open(tmp_path / "comparison/peptide_overlap.tsv.gz", "rt") as stream:
        assert "QQQQK" in stream.read()


@pytest.mark.parametrize("rejected_label", ["1", "-1"])
def test_foreign_acquisition_rejected_even_when_psm_fails(tmp_path, rejected_label):
    left, right = searches(tmp_path)
    with (left / "results.sage.tsv").open("a") as stream:
        stream.write(f"ACDEK\ta;c\t{rejected_label}\t0.9\tforeign.mzML\n")
    result = CliRunner().invoke(
        cli,
        [
            "molecular",
            "compare",
            "--left",
            str(left),
            "--right",
            str(right),
            "--out",
            str(tmp_path / "comparison"),
        ],
    )
    assert result.exit_code != 0 and "different acquisition" in result.output
    assert not (tmp_path / "comparison").exists()


def test_hidden_extra_lfq_acquisition_rejected(tmp_path):
    left, right = searches(tmp_path)
    path = left / "lfq.tsv"
    rows = path.read_text().splitlines()
    path.write_text(rows[0] + "\tforeign.mzML\n" + "\n".join(r + "\t0" for r in rows[1:]) + "\n")
    result = CliRunner().invoke(
        cli,
        [
            "molecular",
            "compare",
            "--left",
            str(left),
            "--right",
            str(right),
            "--out",
            str(tmp_path / "comparison"),
        ],
    )
    assert result.exit_code != 0 and "LFQ intensity column" in result.output
    assert not (tmp_path / "comparison").exists()
