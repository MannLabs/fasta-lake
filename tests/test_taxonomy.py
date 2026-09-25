"""Protect reference-dependent taxonomy, rank limits and quantitative denominators."""

import csv
import gzip
import json

import pytest
from click.testing import CliRunner
from test_simplified_study import P1, P2, P3, search

from fasta_lake.cli import cli
from fasta_lake.study import group_searches
from fasta_lake.taxonomy import (
    build_taxonomy_tree,
    canonical_peptide,
    read_unipept_output,
    taxonomy_report,
)


def row(peptide, *, taxon=1678, name="Bifidobacterium", rank="genus", cutoff=False):
    return dict(
        peptide=peptide,
        taxon_id=taxon,
        taxon_name=name,
        taxon_rank=rank,
        cutoff_used=cutoff,
        domain_id=2,
        domain_name="Bacteria",
        phylum_id=201174,
        phylum_name="Actinomycetota",
        genus_id=taxon,
        genus_name=name,
        species_id=None,
        species_name="",
    )


@pytest.fixture
def study(tmp_path):
    for sample, values in [("S1", [0.5, 3, 6]), ("S2", [0, 2, 0]), ("empty", [0, 0, 0])]:
        search(
            tmp_path,
            sample,
            dict(A=P1, B=P2, C=P3),
            [(p, g, 1, 0.001) for p, g in zip([P1, P2, P3], "ABC")],
            [(p, 2, g, v) for p, g, v in zip([P1, P2, P3], "ABC", values)],
        )
    groups = tmp_path / "groups"
    group_searches(tmp_path / "search", groups)
    annotations = tmp_path / "unipept.json"
    host = row(P1, taxon=9606, name="Homo sapiens", rank="species")
    host.update(
        domain_id=2759,
        domain_name="Eukaryota",
        phylum_id=7711,
        phylum_name="Chordata",
        genus_id=9605,
        genus_name="Homo",
        species_id=9606,
        species_name="Homo sapiens",
    )
    annotations.write_text(json.dumps([host, row(P2), row(P3, cutoff=True)]))
    return groups, annotations, tmp_path / "out"


def test_report_conserves_signal_and_keeps_host_distinct_from_contaminants(study, tmp_path):
    groups, source, out = study
    contaminants = tmp_path / "contaminants.txt"
    contaminants.write_text(P1 + "\n")
    original = source.read_bytes()
    result = taxonomy_report(
        groups,
        source,
        out,
        database="Fixture UniProt",
        equate_il=True,
        contaminant_peptides=contaminants,
    )
    data = json.loads((out / "taxonomy.json").read_text())
    first = data["samples"][0]
    assert first["sample"] == "S1"
    assert first["trees"]["signal"]["count"] == pytest.approx(9.5)
    assert first["balance"]["Host"]["signal"] == 0.5
    assert first["balance"]["Bacteria"]["signal"] == 3
    assert first["balance"]["Unresolved"]["signal"] == 6
    assert first["overlap"]["signal"] == 0.5
    assert first["trees"]["peptides"]["count"] == 3
    assert data["samples"][2]["totals"]["signal"] == 0
    assert (out / "unipept_original.json").read_bytes() == original
    assert result["visualizations"] == "Unipept Visualizations 3.0.0 (MIT)"
    document = (out / "index.html").read_text()
    assert "10.1093/bioinformatics/btab590" in document
    assert "10.1093/bioinformatics/btw039" in document
    assert "data:text/javascript;base64," in document
    assert '<script src="https://' not in document
    if result["static_plots"]["status"] == "PASS":
        for name in ["Holobiont_balance", "Contaminant_reference_overlap"]:
            svg = (out / (name + ".svg")).read_text()
            assert "10.1093/bioinformatics/btw039" in svg
            assert "Fixture UniProt" in svg
            assert (out / (name + ".pdf")).stat().st_size > 1000


def test_genus_assignment_is_not_distributed_into_species(tmp_path):
    source = tmp_path / "native.json"
    source.write_text(json.dumps([row(P1)]))
    records = read_unipept_output(source, {canonical_peptide(P1), canonical_peptide(P2)})
    assert records[canonical_peptide(P2)]["status"] == "not_returned"
    tree = build_taxonomy_tree(records, {canonical_peptide(P1): 2, canonical_peptide(P2): 3})
    assert tree["count"] == 5
    assert not any(
        node.get("rank") == "species" for node in records[canonical_peptide(P1)]["lineage"]
    )
    assert any(c["name"] == "Not returned by Unipept" and c["count"] == 3 for c in tree["children"])


def test_native_input_limits_fail_before_processing_large_files(tmp_path):
    source = tmp_path / "input.json"
    source.write_text(json.dumps([row(P1), row(P2)]))
    with pytest.raises(ValueError, match="file-size limit"):
        read_unipept_output(source, {P1, P2}, max_bytes=1)
    with pytest.raises(ValueError, match="row limit"):
        read_unipept_output(source, {P1, P2}, max_rows=1)


@pytest.mark.parametrize(
    "defect,match",
    [
        ("duplicate", "duplicate"),
        ("foreign", "Foreign"),
        ("descendant", "descendants"),
        ("fractional_taxon", "taxonomy ID"),
        ("missing_lineage", "--all"),
        ("invalid_cutoff", "cutoff"),
        ("wrong_lca", "disagrees"),
    ],
)
def test_malformed_or_misleading_native_results_fail(tmp_path, defect, match):
    records = [row(P1)]
    if defect == "duplicate":
        records *= 2
    elif defect == "foreign":
        records = [row("ACDEFGHIK")]
    elif defect == "descendant":
        records[0].update(species_id=123, species_name="Wrong species")
    elif defect == "fractional_taxon":
        records[0]["taxon_id"] = 2.5
    elif defect == "missing_lineage":
        del records[0]["genus_id"]
    elif defect == "invalid_cutoff":
        records[0]["cutoff_used"] = "maybe"
    elif defect == "wrong_lca":
        records[0]["genus_id"] = 999
    source = tmp_path / "native.json"
    source.write_text(json.dumps(records))
    with pytest.raises(ValueError, match=match):
        read_unipept_output(source, {canonical_peptide(P1)})


def test_legacy_csv_blank_false_and_json_false_agree(tmp_path):
    native = row(P1)
    native["cutoff_used"] = ""
    csvpath = tmp_path / "native.csv"
    with csvpath.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=native)
        writer.writeheader()
        writer.writerow(native)
    result = read_unipept_output(csvpath, {canonical_peptide(P1)})
    assert result[canonical_peptide(P1)]["status"] == "assigned"
    del native["cutoff_used"]
    csvpath.with_suffix(".json").write_text(json.dumps([native]))
    result = read_unipept_output(csvpath.with_suffix(".json"), {canonical_peptide(P1)})
    assert result[canonical_peptide(P1)]["status"] == "cutoff_not_reported"


@pytest.mark.parametrize("defect", ["duplicate", "foreign", "negative", "wrong_group"])
def test_feature_ledger_cannot_change_quantitative_identity(study, defect):
    groups, source, out = study
    path = groups / "feature_assignments.tsv.gz"
    with gzip.open(path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames
        rows = list(reader)
    if defect == "duplicate":
        rows.append(rows[0].copy())
    elif defect == "foreign":
        rows[0]["sample"] = "foreign"
    elif defect == "negative":
        rows[0]["intensity"] = "-1"
    elif defect == "wrong_group":
        rows[0]["study_group"] = "wrong"
    with gzip.open(path, "wt") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError):
        taxonomy_report(groups, source, out, database="Fixture", equate_il=True)
    assert not out.exists()


def test_cli_requires_explicit_il_setting_and_preserves_existing_output(study):
    groups, source, out = study
    args = [
        "taxonomy",
        "report",
        "--groups",
        str(groups),
        "--unipept-output",
        str(source),
        "--out",
        str(out),
        "--database",
        "Fixture",
    ]
    runner = CliRunner()
    result = runner.invoke(cli, args)
    assert result.exit_code != 0 and "--equate" in result.output
    result = runner.invoke(cli, args + ["--equate-il"])
    assert result.exit_code == 0, result.output
    before = (out / "COMPLETE.json").read_bytes()
    assert runner.invoke(cli, args + ["--equate-il"]).exit_code != 0
    assert (out / "COMPLETE.json").read_bytes() == before


def test_mass_gap_is_not_silently_removed():
    with pytest.raises(ValueError):
        canonical_peptide("PEP[113]TIDEK")


def test_negative_or_nonfinite_tree_weights_fail(tmp_path):
    path = tmp_path / "native.json"
    path.write_text(json.dumps([row(P1)]))
    records = read_unipept_output(path, {canonical_peptide(P1)})
    for value in [-1, float("nan"), float("inf")]:
        with pytest.raises(ValueError):
            build_taxonomy_tree(records, {canonical_peptide(P1): value})
