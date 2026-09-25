"""Protect peptide membership, study assignment and portable graph identity."""

import copy
import csv
import gzip
import json

import pytest

from fasta_lake._ckg_worker import _json_for_script
from fasta_lake.networks import (
    export_evidence_network,
    load_evidence_network,
    read_evidence_network,
    select_network_view,
)


def write_table(path, fields, rows):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


@pytest.fixture
def study(tmp_path):
    root = tmp_path / "study"
    root.mkdir()
    summary = dict(
        schema_version=1,
        acquisition_names=["one", "empty", "two"],
        acquisitions=3,
        inferred_study_groups=1,
        pooled_canonical_peptides=2,
        proteins_with_observed_evidence=2,
    )
    (root / "summary.json").write_text(json.dumps(summary))
    (root / "input_manifest.json").write_text("{}")
    write_table(
        root / "study_groups.tsv",
        ["study_group", "selection_rank", "observed_support_peptides", "assigned_peptides"],
        [
            dict(
                study_group="PEPTLDEK",
                selection_rank=0,
                observed_support_peptides=2,
                assigned_peptides=2,
            )
        ],
    )
    write_table(
        root / "peptide_to_study_group.tsv",
        ["canonical_peptide", "study_group"],
        [dict(canonical_peptide=p, study_group="PEPTLDEK") for p in ["PEPTLDEK", "ACDEK"]],
    )
    fields = [
        "sample",
        "study_group",
        "canonical_peptide",
        "reported_proteins",
        "n_reported_proteins_sample",
        "n_reported_proteins_study",
    ]
    rows = [
        dict(zip(fields, values))
        for values in [
            ["one", "PEPTLDEK", "PEPTLDEK", "PEPTLDEK;protein<script>", 2, 2],
            ["two", "PEPTLDEK", "PEPTLDEK", "protein<script>", 1, 2],
            ["one", "PEPTLDEK", "ACDEK", "PEPTLDEK", 1, 1],
        ]
    ]
    write_table(root / "peptide_evidence.tsv.gz", fields, rows)
    return root, fields, rows


def test_shared_memberships_are_separate_from_reporting_assignment(study):
    data = read_evidence_network(study[0])
    assert len(data["nodes"]) == 4 and len(data["edges"]) == 3
    same_label = [n for n in data["nodes"] if n["label"] == "PEPTLDEK"]
    assert len({n["id"] for n in same_label}) == 2
    assert sum(e["assigned_representative"] for e in data["edges"]) == 2
    assert sorted(e["observed_acquisitions"] for e in data["edges"]) == [1, 1, 2]
    assert data["acquisitions"][1] == dict(
        sample="empty", canonical_peptides=0, reported_proteins=0
    )


def test_ckg_script_data_escapes_html_breakout_characters():
    serialized = _json_for_script({"label": "</script><img src=x onerror=alert(1)> &"})
    assert "</script" not in serialized.lower()
    assert all(character not in serialized for character in "<>&")


@pytest.mark.parametrize(
    "field,value,message",
    [
        ("sample", "foreign", "foreign acquisition"),
        ("study_group", "unknown", "reporting assignment"),
        ("n_reported_proteins_sample", "3", "membership count"),
        ("n_reported_proteins_study", "3", "membership"),
        ("reported_proteins", "PEPTLDEK;PEPTLDEK", "membership count"),
    ],
)
def test_inconsistent_evidence_fails(study, field, value, message):
    root, fields, rows = study
    rows[0][field] = value
    write_table(root / "peptide_evidence.tsv.gz", fields, rows)
    with pytest.raises(ValueError, match=message):
        read_evidence_network(root)


def test_duplicate_acquisition_peptide_rejected(study):
    root, fields, rows = study
    write_table(root / "peptide_evidence.tsv.gz", fields, rows + [rows[0]])
    with pytest.raises(ValueError, match="Duplicate acquisition"):
        read_evidence_network(root)


def test_missing_peptide_rejected(study):
    root, fields, rows = study
    write_table(root / "peptide_evidence.tsv.gz", fields, rows[:-1])
    with pytest.raises(ValueError, match="complete assignment"):
        read_evidence_network(root)


@pytest.mark.parametrize(
    "option,value", [("max_rows", 1), ("max_edges", 2), ("max_rows", 0), ("max_edges", False)]
)
def test_limits_fail_instead_of_truncating(study, option, value):
    with pytest.raises(ValueError):
        read_evidence_network(study[0], **{option: value})


def test_export_checksums_and_no_overwrite(study, tmp_path):
    out = tmp_path / "export"
    result = export_evidence_network(study[0], out)
    assert result["status"] == "PASS"
    data = load_evidence_network(out)
    assert data["graph_id"] == result["graph_id"]
    with pytest.raises(FileExistsError):
        export_evidence_network(study[0], out)
    (out / "edges.tsv").write_text("edited")
    with pytest.raises(ValueError, match="checksum"):
        load_evidence_network(out)


def test_view_limits_preserve_endpoints_and_source(study):
    data = read_evidence_network(study[0])
    before = copy.deepcopy(data)
    view = select_network_view(data, max_nodes=3, max_edges=2)
    assert len(view["nodes"]) == 3 and len(view["edges"]) == 2
    assert view["hidden_nodes"] == 1 and view["hidden_edges"] == 1
    ids = {n["id"] for n in view["nodes"]}
    assert all(e["source"] in ids and e["target"] in ids for e in view["edges"])
    assert view["edges"][0]["assigned_representative"]
    assert data == before
    assert select_network_view(data, max_nodes=3, max_edges=2) == view


def test_empty_study_preserves_empty_acquisition(tmp_path):
    root = tmp_path / "empty"
    root.mkdir()
    (root / "summary.json").write_text(
        json.dumps(
            dict(
                schema_version=1,
                acquisition_names=["empty"],
                acquisitions=1,
                inferred_study_groups=0,
                pooled_canonical_peptides=0,
                proteins_with_observed_evidence=0,
            )
        )
    )
    (root / "input_manifest.json").write_text("{}")
    write_table(
        root / "study_groups.tsv",
        ["study_group", "selection_rank", "observed_support_peptides", "assigned_peptides"],
        [],
    )
    write_table(root / "peptide_to_study_group.tsv", ["canonical_peptide", "study_group"], [])
    write_table(
        root / "peptide_evidence.tsv.gz",
        [
            "sample",
            "study_group",
            "canonical_peptide",
            "reported_proteins",
            "n_reported_proteins_sample",
            "n_reported_proteins_study",
        ],
        [],
    )
    data = read_evidence_network(root)
    assert data["nodes"] == [] and data["edges"] == []
    assert select_network_view(data)["hidden_nodes"] == 0
