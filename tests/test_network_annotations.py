"""Protect exact sequence annotation, unknown coverage and functional graph provenance."""

import csv
import json

import pytest
from test_networks import study, write_table  # noqa: F401

from fasta_lake.network_annotations import (
    annotate_evidence_network,
    function_terms,
    load_eggnog_annotations,
)
from fasta_lake.networks import (
    _sha,
    add_function_annotations,
    export_evidence_network,
    load_evidence_network,
    load_network_annotations,
    select_network_view,
)


@pytest.fixture
def annotated_study(study, tmp_path):  # noqa: F811
    root, _, _ = study
    fasta = root / "representatives.fasta"
    fasta.write_text(">PEPTLDEK\nACDEKPEPTLDEK\n")
    summary = json.loads((root / "summary.json").read_text())
    summary["annotation_targets"] = {"fasta": fasta.name, "sha256": _sha(fasta)}
    (root / "summary.json").write_text(json.dumps(summary))
    graph = tmp_path / "graph"
    export_evidence_network(root, graph)
    query = tmp_path / "queries.fasta"
    query.write_text(">alternative_name\nACDEKPEPTLDEK\n>unused\nACDEFGK\n")
    annotation = tmp_path / "input.emapper.annotations"
    annotation.write_text(
        "#query\tseed_ortholog\tKEGG_ko\tEC\tGOs\tPFAMs\n"
        "alternative_name\tseed\tko:K00001,ko:K00001\t1.1.1.1\tGO:0000001\tDomainA\n"
    )
    return graph, annotation, fasta, query, tmp_path / "functions"


def run_annotation(paths, **options):
    graph, annotation, fasta, query, out = paths
    return annotate_evidence_network(graph, annotation, fasta, out, query_fasta=query, **options)


def test_exact_sequence_aliases_retain_shared_peptide_ambiguity(annotated_study):
    paths = annotated_study
    before = (paths[0] / "network.json").read_bytes()
    receipt = run_annotation(paths)
    data = load_network_annotations(paths[-1], receipt["graph_id"])
    assert len(data["proteins"]) == 2
    matched = next(r for r in data["proteins"] if r["protein_id"] == "PEPTLDEK")
    assert matched["match_rule"] == "exact_sequence_alias"
    assert matched["raw_annotation"]["seed_ortholog"] == "seed"
    assert data["status_counts"] == {"assigned": 1, "sequence_not_supplied": 1}
    assert data["unused_annotation_queries"] == ["unused"]
    assert len(data["edges"]) == 4
    assert all(
        e["source"].startswith("protein:") and e["relation"] == "orthology_annotation"
        for e in data["edges"]
    )
    assert (paths[0] / "network.json").read_bytes() == before
    with (paths[-1] / "group_terms_KEGG_ko.tsv").open() as f:
        assert list(csv.DictReader(f, delimiter="\t")) == [
            {"study_group": "PEPTLDEK", "term": "K00001"}
        ]


def test_no_hits_are_preserved_and_do_not_become_function_nodes(annotated_study):
    paths = annotated_study
    paths[1].write_text(paths[1].read_text().splitlines()[0] + "\n")
    receipt = run_annotation(paths)
    assert receipt["status_counts"] == {"no_hit": 1, "sequence_not_supplied": 1}
    assert receipt["annotation_edges"] == 0 and receipt["function_nodes"] == 0


def test_different_sequence_under_same_accession_fails(annotated_study):
    paths = annotated_study
    paths[3].write_text(">PEPTLDEK\nACDEFGK\n")
    paths[1].write_text(paths[1].read_text().replace("alternative_name", "PEPTLDEK"))
    with pytest.raises(ValueError, match="different sequence"):
        run_annotation(paths)
    assert not paths[-1].exists()


def test_i_and_l_protein_variants_are_not_annotation_aliases(annotated_study):
    paths = annotated_study
    paths[3].write_text(paths[3].read_text().replace("PEPTLDEK", "PEPTIDEK"))
    with pytest.raises(ValueError, match="No annotation query matches"):
        run_annotation(paths)


def test_unrecorded_fasta_fails_even_with_matching_accession(annotated_study):
    paths = annotated_study
    paths[2].write_text(paths[2].read_text() + "\n")
    with pytest.raises(ValueError, match="recorded search FASTA"):
        run_annotation(paths)


def test_foreign_annotation_query_fails(annotated_study):
    paths = annotated_study
    paths[1].write_text(paths[1].read_text().replace("alternative_name", "foreign"))
    with pytest.raises(ValueError, match="Foreign or duplicate"):
        run_annotation(paths)


def test_conflicting_annotations_for_same_sequence_fail(annotated_study):
    paths = annotated_study
    paths[3].write_text(paths[3].read_text() + ">duplicate_sequence\nACDEKPEPTLDEK\n")
    with paths[1].open("a") as f:
        f.write("duplicate_sequence\tseed\tko:K00002\t-\t-\t-\n")
    with pytest.raises(ValueError, match="conflicting annotations"):
        run_annotation(paths)


def test_annotation_receipt_rejects_other_graph_and_changed_file(annotated_study):
    paths = annotated_study
    receipt = run_annotation(paths)
    with pytest.raises(ValueError, match="another evidence graph"):
        load_network_annotations(paths[-1], "foreign")
    with (paths[-1] / "protein_annotations.tsv").open("a") as f:
        f.write("changed\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        load_network_annotations(paths[-1], receipt["graph_id"])


def test_existing_output_is_not_overwritten(annotated_study):
    run_annotation(annotated_study)
    with pytest.raises(FileExistsError):
        run_annotation(annotated_study)


def test_annotation_view_preserves_memberships_and_bounds(annotated_study):
    paths = annotated_study
    receipt = run_annotation(paths)
    base = load_evidence_network(paths[0])
    original = select_network_view(base, max_nodes=5, max_edges=5)
    layer = load_network_annotations(paths[-1], receipt["graph_id"])
    view = add_function_annotations(original, layer, max_nodes=6, max_edges=7)
    assert len(view["nodes"]) <= 6 and len(view["edges"]) <= 7
    assert view["edges"][: len(original["edges"])] == original["edges"]
    assert len(original["nodes"]) == 4
    ids = {n["id"] for n in view["nodes"]}
    assert all(e["source"] in ids and e["target"] in ids for e in view["edges"])
    assert any(n["node_type"] == "annotated_function" for n in view["nodes"])
    assert view["hidden_edges"] == view["total_edges"] - len(view["edges"])


def test_eggnog_cli_matches_python_export(annotated_study):
    from click.testing import CliRunner

    from fasta_lake.cli import cli

    paths = annotated_study
    graph, annotation, fasta, query, out = paths
    result = CliRunner().invoke(
        cli,
        [
            "network",
            "eggnog",
            "--bundle",
            str(graph),
            "--annotation",
            str(annotation),
            "--fasta",
            str(fasta),
            "--query-fasta",
            str(query),
            "--out",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.output)["function_nodes"] == 4


@pytest.mark.parametrize(
    "value,field,expected",
    [
        ("-", "KEGG_ko", []),
        ("ko:K00001,K00001", "KEGG_ko", ["K00001"]),
        ("1.2.3.-", "EC", ["1.2.3.-"]),
        ("GO:0000001,GO:0000001", "GOs", ["GO:0000001"]),
        ("Enolase_N,Enolase", "PFAMs", ["Enolase", "Enolase_N"]),
        ("EG", "COG_category", ["E", "G"]),
    ],
)
def test_namespace_and_missing_terms(value, field, expected):
    assert function_terms(value, field) == expected


@pytest.mark.parametrize("value,field", [("K1", "KEGG_ko"), ("1.2", "EC"), ("GO:12", "GOs")])
def test_malformed_term_identifiers_fail(value, field):
    with pytest.raises(ValueError):
        function_terms(value, field)


def test_native_annotations_require_original_query_fasta(annotated_study):
    with pytest.raises(ValueError, match="query-fasta"):
        load_eggnog_annotations(annotated_study[1])
