"""DIANovo formats must preserve confidence, sequence identity and selected targets."""

import csv
import subprocess

import pytest

from fasta_lake.dianovo import parse_dianovo_prediction
from fasta_lake.headers.validation import validate_denovo_csv
from fasta_lake.helpers import collect_peptides_with_scores
from fasta_lake.rust.binaries import find_binary


@pytest.mark.parametrize("vector", ["[0, .9, .9]", "0,.9,.9", "0 .9 .9", "0 0 .9 .9"])
def test_zero_residues_lowercase_and_path_marker(vector):
    assert parse_dianovo_prediction("A c D", vector) == ("ACD", 0.6)


@pytest.mark.parametrize(
    "sequence,vector",
    [
        ("AC113.084D", "[1, 1, 1]"),
        ("AC[+57]D", "[1, 1, 1]"),
        ("AXB", "[1, 1, 1]"),
        ("AßD", "[1, 1, 1, 1]"),
        ("ACD", "[0, 1, 1, 1]"),
        ("ACD", ".1 1 1 1"),
        ("ACD", "[1, bogus, 1]"),
        ("ACD", "[0_0, 1, 1]"),
        ("ACD", "[٠, 1, 1]"),
        ("ACD", "[1, NaN, 1]"),
        ("ACD", "[1, Infinity, 1]"),
        ("ACD", "[1, -0.1, 1]"),
        ("ACD", "[1, 1.1, 1]"),
        ("ACD", "[1, 1, 1"),
        ("ACD", "[1, 1, 1,]"),
        ("ACD", "[1 1 1]"),
        ("ACD", "[]"),
    ],
)
def test_incomplete_sequences_and_invalid_probabilities_are_rejected(sequence, vector):
    with pytest.raises(ValueError):
        parse_dianovo_prediction(sequence, vector)


def write_predictions(path, rows, headers=("pred_seq", "pred_prob")):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerows(rows)


def test_validation_and_python_helper_share_the_contract(tmp_path, caplog):
    path = tmp_path / "predictions.csv"
    write_predictions(
        path, [("ACDEFGHLK", "[0,1,1,1,1,1,1,1,1]"), ("ACD113.0EFGHLK", "[1,1,1,1,1,1,1,1,1]")]
    )
    scores = collect_peptides_with_scores(path)
    assert scores == {"ACDEFGHLK": 8 / 9}
    assert "excluded 1" in caplog.text
    result = next(r for r in validate_denovo_csv(path) if r.check == "dianovo_complete_predictions")
    assert result.details == {
        "excluded": 1,
        "reasons": {"unresolved_gap_or_noncanonical_sequence": 1},
    }
    assert result.severity == "warning"


@pytest.mark.parametrize("path_format", [False, True])
def test_rust_extraction_and_razor_agree_on_zero_aware_ranking(tmp_path, path_format):
    try:
        extractor = find_binary("fasta_extractor")
        razor = find_binary("parsimony_engine")
    except FileNotFoundError:
        pytest.skip("Build the two Rust selection binaries; CI supplies them")
    # Discarding the real zero would incorrectly rank PEPTIDEAK above ACDEFGHIK.
    rows = [
        ("PEPTIDEAK", [0] + [0.9] * 8),
        ("AcDEFGHIK", [0.85] * 9),
        ("MNPEPTIDK", [0.2] * 9),
        ("ACD113.084MNPQKR", [1] * 9),
    ]
    native = tmp_path / "native"
    generic = tmp_path / "generic"
    native.mkdir()
    generic.mkdir()
    write_predictions(
        native / "sample.csv",
        [
            (
                " ".join(seq) if path_format else seq,
                " ".join(map(str, [0] + prob)) if path_format else str(prob),
            )
            for seq, prob in rows
        ],
    )
    write_predictions(
        generic / "sample.csv",
        [("PEPTIDEAK", 0.8), ("ACDEFGHIK", 0.85), ("MNPEPTIDK", 0.2)],
        headers=("sequence", "score"),
    )
    reference = tmp_path / "reference.fasta"
    reference.write_text(
        ">low_confidence_zero\nMPEPTIDEAK\n>correct_top\nMACDEFGHIK\n"
        ">low\nMMNPEPTIDK\n>must_not_join_gap\nMACDMNPQKR\n"
    )
    selected = []
    reduced = []
    for name, predictions in [("native", native), ("generic", generic)]:
        output = tmp_path / (name + ".fasta")
        result = subprocess.run(
            [
                str(extractor),
                "-p",
                str(predictions),
                "-d",
                str(reference),
                "-o",
                str(output),
                "--source-tag",
                "TEST",
                "--min-length",
                "9",
                "--top-percent",
                "0.3",
                "-t",
                "1",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        selected.append(output.read_bytes())
        result = subprocess.run(
            [
                str(razor),
                "--strategy",
                "razor",
                "--fasta",
                str(reference),
                "--peptides",
                str(predictions / "sample.csv"),
                "--output-dir",
                str(tmp_path / name / "razor"),
                "--sample",
                "sample",
                "--tiebreak",
                "hash-acc",
                "--min-length",
                "9",
                "--top-percent",
                "0.3",
                "-t",
                "1",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        files = list((tmp_path / name / "razor").glob("*.fasta"))
        assert len(files) == 1
        reduced.append(files[0].read_bytes())
    assert selected[0] == selected[1]
    assert selected[0].count(b">") == 1 and b"MACDEFGHIK" in selected[0]
    assert reduced[0] == reduced[1]
    assert reduced[0].count(b">") == 1 and b"MACDEFGHIK" in reduced[0]
