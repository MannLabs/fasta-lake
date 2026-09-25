"""Reference-piece boundaries must not change evidence or prediction ranking."""

import json
import subprocess

import pytest
from click.testing import CliRunner

from fasta_lake.chunked import extract_chunked
from fasta_lake.cli import cli
from fasta_lake.rust.binaries import find_binary


@pytest.fixture
def extractor():
    try:
        return find_binary("fasta_extractor")
    except FileNotFoundError:
        pytest.skip("Build Rust extraction fixtures; CI provides them")


@pytest.fixture
def inputs(tmp_path):
    database = tmp_path / "lake.fasta"
    database.write_bytes(
        b">first retained description\r\nMPEPTI\r\nDEAKOTHERPEPK\r\n"
        b">no_hit\nMSSSSSSSSK\n"
        b">second_same_sequence_different_header\nMPEPTIDEAKOTHERPEPK\n"
        b">I_L_variant\nMPEPTLDEAK\n"
        b">tie_a\nMALTERNPEPK\n"
        b">tie_b\nMHIGHPEPTK\n"
        b">low_score_only\nMLOWPEPTK\n"
        b">nonutf8_\xff_header\nMPEPTIDEAK"
    )
    predictions = tmp_path / "predictions"
    predictions.mkdir()
    (predictions / "A.csv").write_text(
        "sequence,score\nMPEPTIDEAK,1\nALTERNPEPK,0.9\nHIGHPEPTK,0.9\nLOWPEPTK,0.1\n"
    )
    (predictions / "B.csv").write_text("sequence,score\nOTHERPEPK,1\nOTHERPEPK,1\n")
    return database, predictions


def whole(extractor, inputs, output, normalize):
    output.mkdir()
    command = [
        str(extractor),
        "--database",
        str(inputs[0]),
        "--peptides-dir",
        str(inputs[1]),
        "--output",
        str(output / "evidence.fasta"),
        "--output-evidence",
        str(output / "evidence.tsv"),
        "--output-hits",
        str(output / "hits.tsv"),
        "--top-percent",
        "0.5",
        "--min-length",
        "8",
        "--normalize-il",
        str(normalize).lower(),
        "--threads",
        "1",
    ]
    result = subprocess.run(command, capture_output=True)
    assert result.returncode == 0, result.stderr
    return json.loads((output / "evidence.stats.json").read_text())


@pytest.mark.parametrize("byte_limit,record_limit", [(1, 500000), (100, 500000), (10**6, 2)])
@pytest.mark.parametrize("normalize", [True, False])
def test_piece_boundaries_preserve_complete_evidence(
    extractor, inputs, tmp_path, byte_limit, record_limit, normalize
):
    expected = whole(extractor, inputs, tmp_path / "whole", normalize)
    result = extract_chunked(
        *inputs,
        tmp_path / "pieces",
        chunk_bytes=byte_limit,
        chunk_records=record_limit,
        top_fraction=0.5,
        min_length=8,
        normalize_il=normalize,
        threads=2,
        output_hits=True,
        binary_path=extractor,
    )
    for name in ("evidence.fasta", "evidence.tsv", "hits.tsv"):
        assert (tmp_path / "pieces" / name).read_bytes() == (tmp_path / "whole" / name).read_bytes()
    for key, value in expected.items():
        if key == "total_seconds":
            continue
        assert result["stats"][key] == value, key
    assert result["stats"]["pieces"] > 1
    assert sum(p["reference_records"] for p in result["pieces"]) == 8
    assert not (tmp_path / "pieces/reference_piece.fasta").exists()
    if byte_limit == 1:
        assert any(p["stats"]["proteins_matched"] == 0 for p in result["pieces"])
    assert b">low_score_only" not in (tmp_path / "pieces/evidence.fasta").read_bytes()
    assert b">tie_a" in (tmp_path / "pieces/evidence.fasta").read_bytes()
    assert b">tie_b" in (tmp_path / "pieces/evidence.fasta").read_bytes()


def test_installed_cli_and_existing_output_protection(extractor, inputs, tmp_path, monkeypatch):
    monkeypatch.setenv("FASTALAKE_BIN_DIR", str(extractor.parent))
    out = tmp_path / "cli"
    args = [
        "extract-chunked",
        "--database",
        str(inputs[0]),
        "--predictions",
        str(inputs[1]),
        "--out",
        str(out),
        "--chunk-records",
        "1",
        "--min-length",
        "8",
    ]
    result = CliRunner().invoke(cli, args)
    assert result.exit_code == 0, result.output
    before = (out / "COMPLETE.json").read_bytes()
    result = CliRunner().invoke(cli, args)
    assert result.exit_code != 0 and "Output exists" in result.output
    assert (out / "COMPLETE.json").read_bytes() == before


@pytest.mark.parametrize(
    "options",
    [{"chunk_bytes": 0}, {"chunk_records": -1}, {"top_fraction": float("nan")}, {"threads": 0}],
)
def test_invalid_options_do_not_create_output(inputs, tmp_path, options):
    out = tmp_path / "invalid"
    with pytest.raises(ValueError):
        extract_chunked(*inputs, out, **options)
    assert not out.exists()


def test_all_no_hit_pieces_produce_an_explicit_empty_result(extractor, inputs, tmp_path):
    inputs[0].write_bytes(b">none1\nMSSSSSSSSK\n>none2\nMQQQQQQQQK\n")
    result = extract_chunked(*inputs, tmp_path / "empty", chunk_records=1, binary_path=extractor)
    assert result["stats"]["proteins_matched"] == 0
    assert result["stats"]["proteins_searched"] == 2
    assert (tmp_path / "empty/evidence.fasta").read_bytes() == b""
    assert len((tmp_path / "empty/evidence.tsv").read_text().splitlines()) == 1


@pytest.mark.parametrize(
    "raw", [b"", b">empty\n", b">\nPEPTIDEK\n", b"not a fasta\n", b">first\nPEPT>IDEK\n"]
)
def test_malformed_reference_cannot_receive_a_completion_receipt(extractor, inputs, tmp_path, raw):
    inputs[0].write_bytes(raw)
    out = tmp_path / "bad"
    with pytest.raises(ValueError):
        extract_chunked(*inputs, out, binary_path=extractor)
    assert not (out / "COMPLETE.json").exists()


def test_prediction_changes_are_rejected_between_pieces(extractor, inputs, tmp_path, monkeypatch):
    from fasta_lake import chunked

    real_run = chunked.run_recorded

    def changed(*args, **kwargs):
        result = real_run(*args, **kwargs)
        (inputs[1] / "B.csv").write_text("sequence,score\nMPEPTIDEAK,1\n")
        return result

    monkeypatch.setattr(chunked, "run_recorded", changed)
    out = tmp_path / "changed"
    with pytest.raises(RuntimeError, match="Prediction inputs changed"):
        extract_chunked(*inputs, out, chunk_records=1, binary_path=extractor)
    assert not (out / "COMPLETE.json").exists()
    assert (out / "pieces/000001/RUN.json").is_file()


@pytest.mark.parametrize("keep", [False, True])
def test_piece_storage_choice_preserves_merged_products(extractor, inputs, tmp_path, keep):
    out = tmp_path / "kept"
    result = extract_chunked(*inputs, out, chunk_records=1, keep_pieces=keep, binary_path=extractor)
    for piece in result["pieces"]:
        folder = out / "pieces" / f"{piece['piece']:06d}"
        assert (folder / "matched.fasta").exists() == keep
        assert (folder / "evidence.tsv").exists() == keep
        assert (folder / "matched.stats.json").is_file()
        assert piece["resources"]["exit_code"] == 0
    assert (out / "evidence.fasta").stat().st_size > 0


def test_reference_changes_cannot_receive_completion(extractor, inputs, tmp_path, monkeypatch):
    from fasta_lake import chunked

    real_run = chunked.run_recorded
    altered = False

    def changed(*args, **kwargs):
        nonlocal altered
        result = real_run(*args, **kwargs)
        if not altered:
            with inputs[0].open("ab") as stream:
                stream.write(b"\n>new_record\nMPEPTIDEAK\n")
            altered = True
        return result

    monkeypatch.setattr(chunked, "run_recorded", changed)
    out = tmp_path / "changed_reference"
    with pytest.raises(RuntimeError, match="Reference FASTA changed"):
        extract_chunked(*inputs, out, chunk_records=1, binary_path=extractor)
    assert not (out / "COMPLETE.json").exists()
