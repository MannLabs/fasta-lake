"""Molecular additions must preserve the baseline and the benchmark's evidence boundaries."""

import csv
import gzip
import json
import shutil

import pytest
from click.testing import CliRunner

from fasta_lake.cli import cli
from fasta_lake.multi_omics.exact import build_union, read_fasta, sha256
from fasta_lake.multi_omics.selection import build_selection
from fasta_lake.multi_omics.sources import load_reference, prepare_source


def inputs(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "proteins.faa").write_text(
        ">g1\nMACDEK\n>g2 duplicate\nMACDEK\n>g3\nMILK\n>g4\nMLLK\n>g5\nMWWWWK\n>foreign\nMHHHHK\n"
    )
    (raw / "tpm.tsv").write_text(
        "prodigal_protein\tmetaG_tpm\tmetaT_tpm\n"
        "g1\t1\tNA\ng2\t10\t0\ng3\t5\t3\ng4\t0\t3\ng5\t0\tnan\n"
        "missing_high\t10000\t10000\n"
    )
    (raw / "annotations.tsv").write_text(
        "prodigal_protein\tKEGG_ko\n"
        "g1\t-\ng2\tK00001\ng3\tK00002\ng4\tK00002\ng5\tko:K00003\n"
        "missing_high\tK00003\n"
    )
    (raw / "metabolites.tsv").write_text(
        "specimen\tcompound\tvalue\ns1\tC00001\t2\ns1\tC00002\t0\ns1\tC00003\tNA\n"
    )
    (raw / "links.tsv").write_text(
        "compound\treaction\tKO\nC00001\tR00001\tK00003\n"
        "C00002\tR00002\tK00001\nC00003\tR00003\tK00002\n"
    )
    (raw / "provenance.json").write_text('{"kind":"synthetic test; no measured TPM"}\n')
    (raw / "baseline.fasta").write_bytes(b">de_novo keep header\r\nMQQQQK\r\n")
    return raw


def prepared(tmp_path):
    raw = inputs(tmp_path)
    prepare_source(
        tmp_path / "source",
        specimen="s1",
        reference_id="assembly_1",
        reference_scope="available",
        tpm=raw / "tpm.tsv",
        proteins=raw / "proteins.faa",
        provenance=raw / "provenance.json",
        annotations=raw / "annotations.tsv",
        metabolites=raw / "metabolites.tsv",
        compound_links=raw / "links.tsv",
    )
    return raw, tmp_path / "source/source.json"


def sequences(path):
    return {r[1] for r in read_fasta(path).values()}


def test_available_reference_retains_missingness_and_is_portable(tmp_path):
    _, manifest = prepared(tmp_path)
    reference = load_reference(manifest)
    assert reference.coverage == {
        "reference_scope": "available",
        "raw_tpm_genes": 6,
        "available_genes": 5,
        "unavailable_genes": 1,
        "proteins_without_tpm": 1,
        "unavailable_positive_dna": 1,
        "unavailable_positive_rna": 1,
        "missing_dna_in_available": 0,
        "missing_rna_in_available": 2,
        "available_genes_with_ko": 4,
        "positive_compounds": 1,
    }
    assert sequences(manifest.parent / "reference.fasta") == {
        "MACDEK",
        "MILK",
        "MLLK",
        "MWWWWK",
    }
    with (manifest.parent / "unavailable_genes.tsv").open() as stream:
        missing = list(csv.DictReader(stream, delimiter="\t"))
    assert missing == [{"gene_id": "missing_high", "metaG_tpm": "10000", "metaT_tpm": "10000"}]
    moved = tmp_path / "moved"
    shutil.copytree(manifest.parent, moved)
    assert load_reference(moved / "source.json").coverage == reference.coverage


@pytest.mark.parametrize("layer", ["dna", "rna", "both", "metabolites", "all"])
def test_every_addition_keeps_every_denovo_protein_and_header(tmp_path, layer):
    raw, manifest = prepared(tmp_path)
    kwargs = dict(dna_percent=0)
    expected = {"MQQQQK"}
    if layer in ("dna", "both", "all"):
        kwargs["dna_percent"] = 50
        expected.add("MACDEK")
    if layer in ("rna", "both", "all"):
        kwargs["rna_percent"] = 50
        expected.update(("MILK", "MLLK"))  # Inclusive RNA boundary tie; I and L stay distinct.
    if layer in ("metabolites", "all"):
        kwargs["metabolites"] = True
        expected.add("MWWWWK")
    result = build_union(raw / "baseline.fasta", manifest, "s1", tmp_path / "out", **kwargs)
    assert sequences(tmp_path / "out/search.fasta") == expected
    assert (
        (tmp_path / "out/search.fasta")
        .read_bytes()
        .startswith((raw / "baseline.fasta").read_bytes())
    )
    assert result["final_sequences"] == len(expected)
    with (tmp_path / "out/gene_evidence.tsv").open() as stream:
        ledger = {r["gene_id"]: r for r in csv.DictReader(stream, delimiter="\t")}
    assert "missing_high" not in ledger and "foreign" not in ledger
    if layer in ("dna", "both", "all"):
        assert ledger["g1"]["dna_pass"] == "0" and ledger["g1"]["included"] == "1"
        assert ledger["g1"]["final_accession"] == ledger["g2"]["final_accession"]
    if layer in ("metabolites", "all"):
        with (tmp_path / "out/metabolite_links.tsv").open() as stream:
            assert list(csv.DictReader(stream, delimiter="\t")) == [
                {"compound": "C00001", "reaction": "R00001", "KO": "K00003"}
            ]
    for name, digest in result["output_sha256"].items():
        assert sha256(tmp_path / "out" / name) == digest


def test_standalone_and_addition_use_the_same_gene_rules_and_background(tmp_path):
    raw, manifest = prepared(tmp_path)
    background = raw / "host.fasta"
    background.write_text(">host\nMSSSSK\n")
    options = dict(dna_percent=100, rna_percent=100, metabolites=True, background=background)
    build_selection(manifest, "s1", tmp_path / "only", **options)
    build_union(raw / "baseline.fasta", manifest, "s1", tmp_path / "plus", **options)
    assert sequences(tmp_path / "plus/search.fasta") == (
        sequences(tmp_path / "only/search.fasta") | {"MQQQQK"}
    )
    assert sequences(tmp_path / "only/search.fasta") == {
        "MACDEK",
        "MILK",
        "MLLK",
        "MWWWWK",
        "MSSSSK",
    }


@pytest.mark.parametrize("defect", ["complete", "alternate", "specimen", "annotation_tpm"])
def test_bad_source_fails_before_creating_output(tmp_path, defect):
    raw = inputs(tmp_path)
    kwargs = dict(
        specimen="s1",
        reference_id="r1",
        reference_scope="available",
        tpm=raw / "tpm.tsv",
        proteins=raw / "proteins.faa",
        provenance=raw / "provenance.json",
        annotations=raw / "annotations.tsv",
        metabolites=raw / "metabolites.tsv",
        compound_links=raw / "links.tsv",
    )
    if defect == "complete":
        kwargs["reference_scope"] = "complete"
    elif defect == "alternate":
        (raw / "alternate.faa").write_text(">g1\nMCHANGEDK\n")
        kwargs["alternate_proteins"] = raw / "alternate.faa"
    elif defect == "specimen":
        (raw / "metabolites.tsv").write_text("specimen\tcompound\tvalue\nother\tC00001\t2\n")
    else:
        (raw / "annotations.tsv").write_text(
            "prodigal_protein\tKEGG_ko\tmetaG_tpm\ng1\tK00001\t1000\n"
        )
    with pytest.raises(ValueError):
        prepare_source(tmp_path / "bad", **kwargs)
    assert not (tmp_path / "bad").exists()


def test_missing_rna_never_passes_and_empty_standalone_is_explicit(tmp_path):
    _, manifest = prepared(tmp_path)
    result = build_selection(manifest, "s1", tmp_path / "rna", rna_absolute=0)
    assert result["rna_measured_genes"] == 3
    assert sequences(tmp_path / "rna/search.fasta") == {"MILK", "MLLK"}
    with pytest.raises(ValueError, match="no target"):
        build_selection(manifest, "s1", tmp_path / "empty", rna_absolute=3)
    assert not (tmp_path / "empty").exists()


@pytest.mark.parametrize(
    "options",
    [
        {"rna_percent": 1, "rna_absolute": 1},
        {"rna_percent": -1},
        {"dna_absolute": float("nan")},
        {"rna_percent": float("inf")},
    ],
)
def test_invalid_combination_fails_before_output(tmp_path, options):
    raw, manifest = prepared(tmp_path)
    with pytest.raises(ValueError):
        build_union(
            raw / "baseline.fasta", manifest, "s1", tmp_path / "bad", dna_percent=0, **options
        )
    assert not (tmp_path / "bad").exists()


def test_installed_cli_contract_for_available_rna_and_metabolites(tmp_path):
    raw, manifest = prepared(tmp_path)
    args = [
        "molecular",
        "add",
        "--source",
        str(manifest),
        "--specimen",
        "s1",
        "--baseline",
        str(raw / "baseline.fasta"),
        "--out",
        str(tmp_path / "out"),
        "--dna-percent",
        "0",
        "--rna-percent",
        "50",
        "--metabolites",
    ]
    result = CliRunner().invoke(cli, args)
    assert result.exit_code == 0, result.output
    assert sequences(tmp_path / "out/search.fasta") == {"MQQQQK", "MILK", "MLLK", "MWWWWK"}
    assert CliRunner().invoke(cli, args).exit_code != 0
    assert json.loads(result.output)["reference_coverage"]["unavailable_genes"] == 1


def test_compressed_tpm_input(tmp_path):
    raw = inputs(tmp_path)
    with gzip.open(raw / "tpm.tsv.gz", "wb") as stream:
        stream.write((raw / "tpm.tsv").read_bytes())
    prepare_source(
        tmp_path / "source",
        specimen="s1",
        reference_id="r",
        reference_scope="available",
        tpm=raw / "tpm.tsv.gz",
        proteins=raw / "proteins.faa",
        provenance=raw / "provenance.json",
    )
    assert load_reference(tmp_path / "source/source.json").coverage["raw_tpm_genes"] == 6


def test_manifest_binds_each_acquisition_to_the_declared_specimen(tmp_path):
    import importlib.util
    from pathlib import Path

    raw, manifest = prepared(tmp_path)
    (raw / "predictions.csv").write_text("sequence,score\nACDEK,1\n")
    (raw / "spectra.mzML").write_text("format fixture; not measured spectra")
    path = Path(__file__).resolve().parents[1] / "tools/run_manifest.py"
    spec = importlib.util.spec_from_file_location("molecular_runner", path)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    acquisition = tmp_path / "samples.tsv"
    header = "sample\tpredictions\tmzml\tmolecular_source\tspecimen\n"
    prefix = "run1\traw/predictions.csv\traw/spectra.mzML\tsource/source.json\t"
    acquisition.write_text(header + prefix + "s1\n")
    rows = runner.load_manifest(acquisition)
    assert rows[0]["molecular_source"] == str(manifest.resolve())
    assert rows[0]["molecular_source_sha256"] == sha256(manifest)
    acquisition.write_text(header + prefix + "another_specimen\n")
    with pytest.raises(ValueError, match="specimen mismatch"):
        runner.load_manifest(acquisition)
    acquisition.write_text(header + prefix + "\n")
    with pytest.raises(ValueError, match="supplied together"):
        runner.load_manifest(acquisition)
