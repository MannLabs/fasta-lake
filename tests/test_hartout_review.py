"""Reachable review regressions: preservation, acquisition identity and provenance."""

import hashlib
import json
import os
import subprocess
from pathlib import Path

import pytest
from click.testing import CliRunner

from fasta_lake.cli import cli
from fasta_lake.presets import PRESETS
from fasta_lake.run_config import RunConfig
from fasta_lake.study import group_searches


def rust_bin(name):
    root = Path(__file__).resolve().parents[1]
    crate = "fasta_export" if name == "fasta_export" else "fasta_extractor_v2"
    path = Path(os.environ.get("FASTALAKE_BIN_DIR", root / "rust" / crate / "target/release"))
    path /= name
    if not path.is_file():
        pytest.skip(f"Build {name} for the Rust regression")
    return str(path)


@pytest.mark.parametrize("mode", ["same", "symlink", "hardlink", "existing"])
def test_contaminant_tag_preserves_existing_destinations(tmp_path, mode):
    source = tmp_path / "source.fasta"
    source.write_text(">P00761\nPEPTIDEAK\n")
    output = source if mode == "same" else tmp_path / "output.fasta"
    if mode == "symlink":
        output.symlink_to(source)
    elif mode == "hardlink":
        os.link(source, output)
    elif mode == "existing":
        output.write_text("unrelated sentinel")
    before = source.read_bytes(), output.read_bytes()
    result = CliRunner().invoke(cli, ["contaminants", "tag", "-i", str(source), "-o", str(output)])
    assert result.exit_code != 0
    assert (source.read_bytes(), output.read_bytes()) == before


@pytest.mark.parametrize("engine", ["sage", "diann", "msfragger"])
@pytest.mark.parametrize("mode", ["same", "symlink", "hardlink", "existing"])
def test_export_preserves_existing_destinations(tmp_path, engine, mode):
    source = tmp_path / "source.fasta"
    source.write_text(">A\nPEPTIDEAK\n")
    output = source if mode == "same" else tmp_path / "output.fasta"
    if mode == "symlink":
        output.symlink_to(source)
    elif mode == "hardlink":
        os.link(source, output)
    elif mode == "existing":
        output.write_text("unrelated sentinel")
    before = source.read_bytes(), output.read_bytes()
    result = subprocess.run(
        [rust_bin("fasta_export"), "export", "-i", str(source), "-o", str(output), "-e", engine],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0, result.stdout
    assert (source.read_bytes(), output.read_bytes()) == before


def test_msfragger_export_preserves_unrelated_temporary_name(tmp_path):
    source, output = tmp_path / "input.fasta", tmp_path / "export.fasta"
    source.write_text(">A\nPEPTIDEAK\n")
    sentinel = tmp_path / "export.tmp.fasta"
    sentinel.write_bytes(b"unrelated existing file")
    result = subprocess.run(
        [
            rust_bin("fasta_export"),
            "export",
            "-i",
            str(source),
            "-o",
            str(output),
            "-e",
            "msfragger",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert sentinel.read_bytes() == b"unrelated existing file"
    assert output.read_text().count(">") == 2
    assert "KAEDITPEP" in output.read_text()


@pytest.mark.parametrize(
    "option", [[], ["--reassign-header-tags", "ASM"], ["--uniform-tag", "SEQ"]]
)
def test_lake_sidecar_names_the_emitted_header(tmp_path, option):
    source, output, sidecar = [tmp_path / p for p in ["input.fasta", "lake.fasta", "headers.tsv"]]
    source.write_text(">Z_first\nPEPTIDEAK\n>A_second\nPEPTIDEAK\n")
    subprocess.run(
        [
            rust_bin("lake_builder"),
            "-s",
            f"ASM:1:{source}",
            "-o",
            str(output),
            "--output-headers",
            str(sidecar),
            "-t",
            "2",
            *option,
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    emitted = output.read_text().splitlines()[0][1:]
    # The legacy sidecar stores source/header pairs with embedded tabs in its
    # fourth field; its established reader splits at most three times.
    header, *lines = sidecar.read_text().splitlines()
    row = [dict(zip(header.split("\t"), line.split("\t", 3))) for line in lines]
    assert len(row) == 1
    assert row[0]["kept_header"] == emitted
    assert "Z_first" in row[0]["all_source_headers"]
    assert "A_second" in row[0]["all_source_headers"]


def test_preset_commands_cover_every_registered_preset():
    runner = CliRunner()
    result = runner.invoke(cli, ["preset", "list"])
    assert result.exit_code == 0, result.output
    for name in PRESETS:
        assert name in result.output
        for command in ["show", "sage-config"]:
            r = runner.invoke(cli, ["preset", command, name])
            assert r.exit_code == 0, r.output
    assert "0.001" in result.output
    assert "whole-workflow" in result.output


def test_run_config_hash_covers_changes_beyond_first_mib(tmp_path):
    prefix = b">A\n" + b"A" * (1 << 20)
    hashes = []
    for suffix in [b"PEPTIDEAK\n", b"OTHERPEPK\n"]:
        path = tmp_path / (str(len(hashes)) + ".fasta")
        path.write_bytes(prefix + suffix)
        cfg = RunConfig(lake_path=str(path))
        cfg.capture_environment()
        assert cfg.lake_md5 == hashlib.md5(path.read_bytes()).hexdigest()
        hashes.append(cfg.lake_md5)
    assert hashes[0] != hashes[1]


def acquisition(root, name, intensity, *, raw=None, search_name=None, lfq_name=None):
    sample = root / name
    sample.mkdir(parents=True)
    fasta = sample / "input.fasta"
    fasta.write_text(">A\nPEPTIDEAK\n")
    raw = raw or name + ".mzML"
    (sample / "config.json").write_text(
        json.dumps(
            {
                "database": {"fasta": str(fasta)},
                "mzml_paths": ["/recorded/raw/" + raw],
            }
        )
    )
    (sample / "results.sage.tsv").write_text(
        "peptide\tproteins\tlabel\tpeptide_q\tfilename\n"
        f"PEPTIDEAK\tA\t1\t0.001\t{search_name or raw}\n"
    )
    (sample / "lfq.tsv").write_text(
        f"peptide\tproteins\t{lfq_name or raw}\nPEPTIDEAK\tA\t{intensity}\n"
    )
    return sample


def test_grouping_rejects_swapped_lfq_even_when_intensity_conserves(tmp_path):
    root = tmp_path / "search"
    acquisition(root, "case", 10, lfq_name="control.mzML")
    acquisition(root, "control", 100, lfq_name="case.mzML")
    with pytest.raises(ValueError, match="LFQ acquisition identity"):
        group_searches(root, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_grouping_rejects_swapped_search_file(tmp_path):
    root = tmp_path / "search"
    acquisition(root, "case", 10, search_name="control.mzML")
    with pytest.raises(ValueError, match="search acquisition identity"):
        group_searches(root, tmp_path / "out")


def test_explicit_alias_and_sage_full_path_are_valid(tmp_path):
    root = tmp_path / "search"
    acquisition(root, "case", 10, raw="acq001.mzML", lfq_name="C:\\raw\\acq001.mzML")
    result = group_searches(root, tmp_path / "out")
    assert result["intensity_in"] == result["intensity_out"] == 10
    assert result["acquisition_identity_checked"]


def test_empty_lfq_still_validates_acquisition_header(tmp_path):
    root = tmp_path / "search"
    sample = acquisition(root, "case", 10)
    (sample / "lfq.tsv").write_text("peptide\tproteins\tcontrol.mzML\n")
    with pytest.raises(ValueError, match="LFQ acquisition identity"):
        group_searches(root, tmp_path / "out")


def test_evidence_metadata_is_repeatable(tmp_path):
    source, predictions = tmp_path / "lake.fasta", tmp_path / "pred.csv"
    peptides = ["ACDEFGHIK", "LMNPQRSTV", "WYACDEFGH"]
    source.write_text(">A\n" + "".join(peptides) + "\n")
    predictions.write_text(
        "sequence,score\n"
        + "\n".join(f"{p},{score}" for p, score in zip(peptides, [1e16, 1.0, 1.0]))
        + "\n"
    )
    evidence, fastas = set(), set()
    for i in range(32):
        output, table = tmp_path / f"o{i}.fasta", tmp_path / f"e{i}.tsv"
        subprocess.run(
            [
                rust_bin("fasta_extractor"),
                "-d",
                str(source),
                "-p",
                str(predictions),
                "-o",
                str(output),
                "--output-evidence",
                str(table),
                "--min-length",
                "9",
                "--min-score",
                "0",
                "-t",
                "1",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        evidence.add(table.read_bytes())
        fastas.add(output.read_bytes())
    assert len(fastas) == len(evidence) == 1
