"""Failures that would change biological evidence or quantification silently."""

import gzip
import importlib.metadata
import json

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner
from test_simplified_study import P1, P2, search, table

from fasta_lake.cli import cli
from fasta_lake.quantification import quantify_study
from fasta_lake.study import group_searches


@pytest.fixture
def grouped(tmp_path):
    for sample, factor in [("NA", 1), ("S2", 2), ("empty", 0)]:
        search(
            tmp_path,
            sample,
            {"A": P1 + P2, "B": P1},
            [(P1, "A;B", 1, 0.001), (P2, "A", 1, 0.001)],
            [
                (P1, 2, "A;B", factor * 100),
                (P1, 2, "A;B", factor * 25),
                (P1, 3, "A;B", factor * 30),
                (P2, 2, "A", factor * 10),
            ],
        )
    group_searches(tmp_path / "search", tmp_path / "groups")
    return tmp_path / "groups"


def require_backend():
    pytest.importorskip("directlfq")
    if importlib.metadata.version("directlfq") != "0.3.3":
        pytest.skip("Integration targets the audited directlfq 0.3.3 runtime")


def test_sum_is_byte_identical_and_cli_works_without_loading_backend(
    grouped, tmp_path, monkeypatch
):
    import sys

    monkeypatch.setitem(sys.modules, "directlfq", None)
    result = CliRunner().invoke(
        cli,
        ["quantify-study", "--groups", str(grouped), "--out", str(tmp_path / "sum"), "--skip-qc"],
    )
    assert result.exit_code == 0, result.output
    assert (tmp_path / "sum/study_group_matrix.tsv").read_bytes() == (
        grouped / "study_group_matrix.tsv"
    ).read_bytes()
    assert json.loads(result.output)["method"] == "sum"


def test_identification_evidence_retains_unquantified_and_ambiguity(grouped):
    evidence = table(grouped / "peptide_evidence.tsv.gz")
    assert len(evidence) == 6
    assert len([r for r in evidence if r["sample"] == "empty"]) == 2
    first = next(r for r in evidence if r["canonical_peptide"] == P1)
    assert first["n_reported_proteins_sample"] == first["n_reported_proteins_study"] == "2"
    assert first["study_group"] == "A"
    assert table(grouped / "peptide_ambiguity.tsv") == [
        {"n_reported_proteins_study": "1", "n_canonical_peptides": "1"},
        {"n_reported_proteins_study": "2", "n_canonical_peptides": "1"},
    ]


def test_local_and_study_ambiguity_are_different(tmp_path):
    search(tmp_path, "S1", {"B": P1}, [(P1, "B", 1, 0.001)], [(P1, 2, "B", 10)])
    search(tmp_path, "S2", {"A": P1 + P2}, [(P1, "A", 1, 0.001), (P2, "A", 1, 0.001)], [])
    group_searches(tmp_path / "search", tmp_path / "groups")
    first = table(tmp_path / "groups/peptide_evidence.tsv.gz")[0]
    assert (first["reported_proteins"], first["study_group"]) == ("B", "A")
    assert (first["n_reported_proteins_sample"], first["n_reported_proteins_study"]) == ("1", "2")


def test_real_backend_preserves_duplicates_empty_sample_and_peptide_support(grouped, tmp_path):
    require_backend()
    out = tmp_path / "lfq"
    record = quantify_study(grouped, out, method="directlfq")
    ions = pd.read_csv(out / "ions.aq_reformat.tsv", sep="\t")
    assert list(ions.columns) == ["protein", "ion", "NA", "S2", "empty"]
    assert len(ions) == 3
    assert ions.set_index("ion").loc[P1 + "|charge=2", "NA"] == 125
    assert not ions["empty"].any()
    result = pd.read_csv(out / "study_group_matrix.tsv", sep="\t").set_index("study_group")
    assert result["empty"].isna().all()
    np.testing.assert_allclose(result["NA"], result["S2"])
    support = table(out / "peptide_support.tsv.gz")
    p1 = next(r for r in support if r["sample"] == "NA" and r["canonical_peptide"] == P1)
    assert p1["input_ions"] == p1["retained_ions"] == "2"
    assert float(p1["input_intensity"]) == 155
    assert record["input_positive_ion_cells"] == record["retained_positive_ion_cells"] == 6
    assert record["baseline_sums_reproduced"]


@pytest.mark.parametrize(
    "defect,match",
    [
        ("duplicate", "counted more than once"),
        ("charge", "charge"),
        ("group", "dictionary"),
        ("intensity", "reproduce"),
    ],
)
def test_corrupt_ledger_is_rejected_before_output(grouped, tmp_path, defect, match):
    require_backend()
    path = grouped / "feature_assignments.tsv.gz"
    with gzip.open(path, "rt") as stream:
        lines = stream.readlines()
    if defect == "duplicate":
        lines.append(lines[1])
    else:
        header, row = lines[0].strip().split("\t"), lines[1].strip().split("\t")
        field = {"charge": "charge", "group": "study_group", "intensity": "intensity"}[defect]
        row[header.index(field)] = {"charge": "", "group": "B", "intensity": "12345"}[defect]
        lines[1] = "\t".join(row) + "\n"
    with gzip.open(path, "wt") as stream:
        stream.writelines(lines)
    with pytest.raises(ValueError, match=match):
        quantify_study(grouped, tmp_path / "lfq", method="directlfq")
    assert not (tmp_path / "lfq").exists()


def test_existing_output_and_unknown_method_fail(grouped, tmp_path):
    with pytest.raises(FileExistsError):
        quantify_study(grouped, grouped)
    with pytest.raises(ValueError, match="method"):
        quantify_study(grouped, tmp_path / "out", method="unknown")


def test_all_unquantified_is_explicit(tmp_path):
    require_backend()
    search(tmp_path, "S", {"A": P1}, [(P1, "A", 1, 0.001)], [])
    group_searches(tmp_path / "search", tmp_path / "groups")
    result = quantify_study(tmp_path / "groups", tmp_path / "lfq", method="directlfq")
    assert result["backend_ran"] is False
    assert (tmp_path / "lfq/study_group_matrix.tsv").read_text() == "study_group\tS\n"


def test_sage_combined_charge_is_retained_as_one_trace(tmp_path):
    require_backend()
    for sample in ["S1", "S2"]:
        search(
            tmp_path,
            sample,
            {"A": P1 + P2},
            [(P1, "A", 1, 0.001), (P2, "A", 1, 0.001)],
            [(P1, -1, "A", 100), (P2, -1, "A", 10)],
        )
    group_searches(tmp_path / "search", tmp_path / "groups")
    quantify_study(tmp_path / "groups", tmp_path / "lfq", method="directlfq")
    ions = table(tmp_path / "lfq/ion_mapping.tsv")
    assert len(ions) == 2
    assert all(r["ion"].endswith("|charge=-1") for r in ions)


def test_mixed_combined_and_resolved_charges_fail(tmp_path):
    require_backend()
    for sample, charge in [("S1", -1), ("S2", 2)]:
        search(tmp_path, sample, {"A": P1}, [(P1, "A", 1, 0.001)], [(P1, charge, "A", 100)])
    group_searches(tmp_path / "search", tmp_path / "groups")
    with pytest.raises(ValueError, match="mixes combined-charge"):
        quantify_study(tmp_path / "groups", tmp_path / "lfq", method="directlfq")


def test_il_equivalent_ions_merge_but_modifications_remain_separate(tmp_path):
    require_backend()
    for sample in ["S1", "S2"]:
        search(
            tmp_path,
            sample,
            {"A": "MPEPTIDEK"},
            [("MPEPTIDEK", "A", 1, 0.001)],
            [
                ("M[+15.9949]PEPTIDEK", -1, "A", 100),
                ("M[+15.9949]PEPTLDEK", -1, "A", 25),
                ("MPEPTIDEK", -1, "A", 10),
            ],
        )
    group_searches(tmp_path / "search", tmp_path / "groups")
    quantify_study(tmp_path / "groups", tmp_path / "lfq", method="directlfq")
    ions = table(tmp_path / "lfq/ions.aq_reformat.tsv")
    assert {r["ion"] for r in ions} == {"M[+15.9949]PEPTLDEK|charge=-1", "MPEPTLDEK|charge=-1"}
    assert sorted(float(r["S1"]) for r in ions) == [10, 125]
    support = table(tmp_path / "lfq/peptide_support.tsv.gz")
    assert len(support) == 2
    assert all(r["input_ions"] == r["retained_ions"] == "2" for r in support)


def test_importer_cannot_silently_change_measurements(grouped, tmp_path, monkeypatch):
    require_backend()
    from directlfq import utils

    original = utils.import_data

    def corrupt(*args, **kwargs):
        frame = original(*args, **kwargs)
        frame.loc[0, "S2"] *= 1.01
        return frame

    monkeypatch.setattr(utils, "import_data", corrupt)
    with pytest.raises(ValueError, match="importer changed ion"):
        quantify_study(grouped, tmp_path / "lfq", method="directlfq")
    assert not (tmp_path / "lfq").exists()
