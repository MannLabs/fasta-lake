"""Adversarial tests for reference-wide protein evidence and real file/CLI accounting."""

import csv
import gzip
import hashlib
import json

import pytest
from click.testing import CliRunner

from fasta_lake.multi_omics.cli import molecular
from fasta_lake.multi_omics.exact import read_fasta
from fasta_lake.multi_omics.increments import evidence_tables, increment_report, reference_map


def digest(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


@pytest.fixture
def searches(tmp_path):
    reference = tmp_path / "reference.fasta"
    reference.write_text(">A\nAAAAKCCCCK\n>B\nAAAAKDDDDK\n>C\nEEEEKFFFFK\n")
    spectrum = tmp_path / "run.mzML"
    spectrum.write_text("immutable test spectrum")
    for name, fasta, peptides in (
        ("left", ">A\nAAAAKCCCCK\n>C\nEEEEKFFFFK\n", [("AAAAK", "A"), ("EEEEK", "C")]),
        ("right", reference.read_text(), [("AAAAK", "A;B"), ("DDDDK", "B")]),
    ):
        folder = tmp_path / name
        folder.mkdir()
        (folder / "db.fasta").write_text(fasta)
        (folder / "config.json").write_text(
            json.dumps(
                {
                    "database": {"fasta": str(folder / "db.fasta"), "decoy_tag": "rev_"},
                    "mzml_paths": [str(spectrum)],
                    "output_directory": str(folder),
                }
            )
        )
        with (folder / "results.sage.tsv").open("w") as f:
            f.write("peptide\tproteins\tlabel\tpeptide_q\tfilename\n")
            for p, proteins in peptides:
                f.write(f"{p}\t{proteins}\t1\t0.005\t{spectrum.name}\n")
        with (folder / "lfq.tsv").open("w") as f:
            f.write(f"peptide\tcharge\tproteins\t{spectrum.name}\n")
            for p, proteins in peptides:
                f.write(f"{p}\t2\t{proteins}\t100\n")
    return tmp_path / "left", tmp_path / "right", reference


def test_shared_support_and_threshold_crossing(searches, tmp_path):
    result = increment_report(*searches, tmp_path / "report")
    with gzip.open(tmp_path / "report/protein_support.tsv.gz", "rt") as f:
        rows = {r["protein_sha256"]: r for r in csv.DictReader(f, delimiter="\t")}
    b = rows[digest("AAAAKDDDDK")]
    assert b["one_peptide_change"] == "retained"
    assert b["two_with_specific_change"] == "gained"
    assert b["left_reference_specific_peptides"] == "0"
    assert b["in_left_candidates"] == "0"  # Shared support outside selected FASTA is explicit.
    assert b["left_reported_peptides"] == "0"
    assert b["gained_peptide_sequences"] == "DDDDK"
    assert result["candidate_sequences_added"] == 1
    assert result["one_peptide_gained"] == 0  # Adding a candidate is not new peptide support.
    assert result["two_with_specific_gained_on_added_candidates"] == 1
    assert result["peptides_gained"] == result["peptides_lost"] == 1
    assert result["lfq_shared"] == result["lfq_gained"] == result["lfq_lost"] == 1
    assert result["right_identity_recovered_by_left"] == 0.5
    assert result["lfq_ratio"] == 1.0  # Equal counts can hide different identities.


def test_reversing_arms_reverses_gain_and_loss(searches, tmp_path):
    a, b, reference = searches
    forward = increment_report(a, b, reference, tmp_path / "f")
    reverse = increment_report(b, a, reference, tmp_path / "r")
    for endpoint in ("peptides", "one_peptide", "two_peptides", "two_with_specific", "lfq"):
        assert forward[endpoint + "_gained"] == reverse[endpoint + "_lost"]
        assert forward[endpoint + "_lost"] == reverse[endpoint + "_gained"]
        retained = "lfq_shared" if endpoint == "lfq" else endpoint + "_retained"
        assert forward[retained] == reverse[retained]


def test_il_equivalence_does_not_collapse_protein_identity(tmp_path):
    p = tmp_path / "il.fasta"
    p.write_text(">one\nMKIIAK\n>alias\nMKIIAK\n>two\nMKLLAK\n")
    sequences, support, matches = reference_map(read_fasta(p), {"MKLLAK"})
    assert len(sequences) == len(matches["MKLLAK"]) == 2
    assert len(support) == 2


def test_reference_omission_is_rejected(searches, tmp_path):
    a, b, ref = searches
    ref.write_text(">A\nAAAAKCCCCK\n")
    with pytest.raises(ValueError, match="absent from declared reference"):
        increment_report(a, b, ref, tmp_path / "bad")
    assert not (tmp_path / "bad").exists()


def test_host_ambiguity_uses_whole_reference(searches, tmp_path):
    host = tmp_path / "host.fasta"
    host.write_text(">B\nAAAAKDDDDK\n")
    result = increment_report(*searches, tmp_path / "host_report", host_reference=host)
    assert result["peptides_retained_host_nonhost_ambiguous"] == 1
    assert result["peptides_gained_host_only"] == 1
    assert result["peptides_lost_nonhost_only"] == 1


def test_settings_mismatch_and_no_overwrite(searches, tmp_path):
    increment_report(*searches, tmp_path / "complete")
    with pytest.raises(FileExistsError):
        increment_report(*searches, tmp_path / "complete")
    config = searches[1] / "config.json"
    settings = json.loads(config.read_text())
    settings["changed_tolerance"] = 1
    config.write_text(json.dumps(settings))
    with pytest.raises(ValueError, match="settings differ"):
        increment_report(*searches, tmp_path / "mismatch")


def test_cli_and_receipt(searches, tmp_path):
    left, right, reference = searches
    output = tmp_path / "cli"
    result = CliRunner().invoke(
        molecular,
        [
            "increments",
            "--left",
            str(left),
            "--right",
            str(right),
            "--reference",
            str(reference),
            "--out",
            str(output),
        ],
    )
    assert result.exit_code == 0, result.output
    receipt = json.loads((output / "COMPLETE.json").read_text())
    assert receipt["status"] == "PASS"
    for name, expected in receipt["output_sha256"].items():
        assert hashlib.sha256((output / name).read_bytes()).hexdigest() == expected


def test_zero_denominators_are_undefined():
    arm = dict(accepted=set(), aliases={}, support={}, lfq={})
    result, proteins, peptides = evidence_tables(arm, arm, {}, {}, {})
    assert result["lfq_ratio"] is None
    assert result["right_identity_recovered_by_left"] is None
    assert proteins == peptides == []


def test_empty_report_keeps_the_complete_export_schema(searches, tmp_path):
    increment_report(*searches, tmp_path / "nonempty")
    for folder in searches[:2]:
        for name in ("results.sage.tsv", "lfq.tsv"):
            path = folder / name
            path.write_text(path.read_text().splitlines()[0] + "\n")
    summary = increment_report(*searches, tmp_path / "empty")
    assert summary["peptides_left"] == summary["peptides_right"] == 0
    assert summary["lfq_ratio"] is None
    for name in ("protein_support.tsv.gz", "peptide_changes.tsv.gz"):
        with gzip.open(tmp_path / "empty" / name, "rt") as empty:
            with gzip.open(tmp_path / "nonempty" / name, "rt") as nonempty:
                assert empty.readline() == nonempty.readline()
            assert empty.read() == ""


def test_reference_change_during_parsing_cannot_receive_a_pass(searches, tmp_path, monkeypatch):
    from fasta_lake.multi_omics import increments

    original = increments.read_fasta
    reference = searches[2]

    def changed(path):
        result = original(path)
        if path == reference:
            reference.write_text(reference.read_text() + ">new_alternative\nAAAAKCCCCKFFFFK\n")
        return result

    monkeypatch.setattr(increments, "read_fasta", changed)
    with pytest.raises(ValueError, match="Source changed"):
        increment_report(*searches, tmp_path / "changed")
    assert not (tmp_path / "changed").exists()


def test_dangling_output_link_is_rejected_before_reading_inputs(tmp_path):
    output = tmp_path / "linked"
    output.symlink_to(tmp_path / "missing")
    with pytest.raises(FileExistsError):
        increment_report(None, None, None, output)
    assert output.is_symlink()


@pytest.mark.parametrize("q", [-1, 2, float("nan")])
def test_invalid_threshold(searches, tmp_path, q):
    with pytest.raises(ValueError, match="Peptide q"):
        increment_report(*searches, tmp_path / "q", q=q)
