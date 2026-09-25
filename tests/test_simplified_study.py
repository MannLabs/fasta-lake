"""Scientific regression cases for direct peptide assignment and one reporting table."""

import csv
import gzip
import json
from pathlib import Path

import pytest

from fasta_lake.study import group_searches, study_dictionary

P1, P2, P3, P4, P5 = "AAAAAAK", "CCCCCCK", "DDDDDDK", "EEEEEEK", "FFFFFFK"


def search(tmp_path, name, sequences, observations, features, config_name="config.json"):
    folder = tmp_path / "search" / name
    folder.mkdir(parents=True)
    fasta = tmp_path / f"{name}.fasta"
    fasta.write_text("".join(f">{p}\n{s}\n" for p, s in sequences.items()))
    (folder / config_name).write_text(
        json.dumps({"database": {"fasta": str(fasta)}, "mzml_paths": [name + ".mzML"]})
    )
    with (folder / "results.sage.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["peptide", "proteins", "label", "peptide_q", "filename"])
        w.writerows([*row, name + ".mzML"] for row in observations)
    with (folder / "lfq.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            ["peptide", "charge", "proteins", "q_value", "score", "spectral_angle", name + ".mzML"]
        )
        w.writerows(
            (p, c, members, 0.5, 1, 0.9, intensity) for p, c, members, intensity in features
        )
    return folder


def table(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_shared_peptide_follows_explaining_group_not_whole_protein_plurality(tmp_path):
    graph = {"A": {P1, P2, P3}, "B": {P1, P2, P4}, "C": {P4, P5}}
    dictionary = study_dictionary(graph, set(graph))
    assert dictionary["selected"] == [("A", 3), ("C", 2)]
    assert dictionary["peptide_to_group"][P4] == "C"
    search(
        tmp_path,
        "S",
        {p: "".join(sorted(v)) for p, v in graph.items()},
        [
            (P1, "A;B", 1, 0.001),
            (P2, "A;B", 1, 0.001),
            (P3, "A", 1, 0.001),
            (P4, "B;C", 1, 0.001),
            (P5, "C", 1, 0.001),
        ],
        [
            (P1, 2, "A;B", 10),
            (P4, 2, "B;C", 20),
            (P4[:1] + "[+1.0]" + P4[1:], 3, "B;C", 30),
            (P5, 2, "C", 40),
        ],
    )
    out = tmp_path / "grouped"
    result = group_searches(tmp_path / "search", out)
    rows = {r["study_group"]: r for r in table(out / "study_group_long.tsv")}
    assert len(rows) == 2
    assert float(rows["C"]["intensity"]) == 90
    assert (rows["C"]["n_peptides"], rows["C"]["n_features"]) == ("2", "3")
    features = table(out / "feature_assignments.tsv.gz")
    assert [r["study_group"] for r in features if r["canonical_peptide"] == P4] == ["C", "C"]
    assert all(r["reported_proteins"] == "B;C" for r in features if r["canonical_peptide"] == P4)
    assert result["intensity_in"] == result["intensity_out"] == 100
    assert result["lfq_features"] == len(features) == 4
    assert table(out / "study_group_matrix.tsv") == [
        {"study_group": "A", "S": "10.0"},
        {"study_group": "C", "S": "90.0"},
    ]
    assert table(out / "acquisition_summary.tsv")[0]["quantified_study_groups"] == "2"
    assert not (out / "dual_key_long.tsv").exists()
    assert not (out / "protein_to_study_group.tsv").exists()


def test_local_member_sets_no_longer_partition_the_same_group(tmp_path):
    search(
        tmp_path,
        "S",
        {"A": P1, "B": P1 + P2},
        [(P1, "A;B", 1, 0.001), (P2, "B", 1, 0.001)],
        [(P1, 2, "A;B", 10), (P2, 2, "B", 40)],
    )
    out = tmp_path / "out"
    group_searches(tmp_path / "search", out)
    assert table(out / "study_group_long.tsv") == [
        {
            "sample": "S",
            "study_group": "B",
            "n_peptides": "2",
            "n_features": "2",
            "intensity": "50.0",
        },
    ]
    assert {r["member_set"] for r in table(out / "feature_assignments.tsv.gz")} == {"A;B", "B"}


def test_pooled_representative_absent_locally_is_exposed_not_called_local_identification(tmp_path):
    search(tmp_path, "S1", {"B": P1}, [(P1, "B", 1, 0.001)], [(P1, 2, "B", 10)])
    search(
        tmp_path,
        "S2",
        {"A": P1 + P2},
        [(P1, "A", 1, 0.001), (P2, "A", 1, 0.001)],
        [(P2, 2, "A", 20)],
    )
    out = tmp_path / "out"
    result = group_searches(tmp_path / "search", out)
    first = table(out / "feature_assignments.tsv.gz")[0]
    assert (first["study_group"], first["member_set"]) == ("A", "B")
    assert first["representative_in_local_fasta"] == "0"
    assert first["representative_in_local_members"] == "0"
    assert result["features_representative_outside_local_fasta"] == 1


def test_sequence_identity_conflicts_fail_before_output_creation(tmp_path):
    search(tmp_path, "S1", {"A": P1}, [(P1, "A", 1, 0.001)], [(P1, 2, "A", 10)])
    search(tmp_path, "S2", {"A": P2}, [(P2, "A", 1, 0.001)], [(P2, 2, "A", 20)])
    with pytest.raises(ValueError, match="accession-to-sequence conflict"):
        group_searches(tmp_path / "search", tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_representative_sequence_must_explain_the_assigned_peptide(tmp_path):
    search(tmp_path, "S", {"A": P1}, [(P2, "A", 1, 0.001)], [(P2, 2, "A", 10)])
    with pytest.raises(ValueError, match="does not contain assigned peptide"):
        group_searches(tmp_path / "search", tmp_path / "out")
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("where", ["search", "lfq"])
def test_other_samples_fasta_cannot_legitimise_crosswired_local_members(tmp_path, where):
    search(
        tmp_path,
        "S1",
        {"A": P1},
        [(P1, "B" if where == "search" else "A", 1, 0.001)],
        [(P1, 2, "B" if where == "lfq" else "A", 10)],
    )
    search(tmp_path, "S2", {"B": P1}, [(P1, "B", 1, 0.001)], [(P1, 2, "B", 10)])
    with pytest.raises(ValueError, match="S1.*searched FASTA"):
        group_searches(tmp_path / "search", tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_acceptance_and_finite_positive_intensity_match_the_declared_population(tmp_path):
    peptide = "PEPTIDEK"
    search(
        tmp_path,
        "S",
        {"A": peptide + P2 + P3 + P4},
        [(peptide, "A", 1, 0.001), (P2, "A", 1, "nan"), (P3, "A", -1, 0.001), (P4, "A", 1, -0.1)],
        [
            (peptide, 2, "A", 10),
            (peptide, 3, "A", "nan"),
            (peptide, 4, "A", "inf"),
            (peptide, 5, "A", -1),
            (P2, 2, "A", 20),
            (P3, 2, "A", 30),
            (P4, 2, "A", 40),
        ],
        config_name="results.json",
    )
    result = group_searches(tmp_path / "search", tmp_path / "out", config_name="results.json")
    assert result["intensity_out"] == 10
    assert result["lfq_features"] == 1
    assert table(tmp_path / "out/peptide_to_study_group.tsv")[0]["canonical_peptide"] == "PEPTLDEK"


def test_aggregate_overflow_leaves_no_finished_output(tmp_path):
    search(
        tmp_path, "S", {"A": P1}, [(P1, "A", 1, 0.001)], [(P1, 2, "A", 1e308), (P1, 3, "A", 1e308)]
    )
    with pytest.raises(ValueError, match="Non-finite aggregate intensity"):
        group_searches(tmp_path / "search", tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_existing_output_and_nested_output_are_rejected(tmp_path):
    search(tmp_path, "S", {"A": P1}, [(P1, "A", 1, 0.001)], [(P1, 2, "A", 10)])
    out = tmp_path / "out"
    out.mkdir()
    (out / "keep").write_text("original")
    with pytest.raises(FileExistsError):
        group_searches(tmp_path / "search", out)
    assert (out / "keep").read_text() == "original"
    with pytest.raises(ValueError, match="outside the search input"):
        group_searches(tmp_path / "search", tmp_path / "search/new")


def test_excluded_acquisitions_cannot_influence_the_study_dictionary(tmp_path):
    search(tmp_path, "included", {"A": P1}, [(P1, "A", 1, 0.001)], [(P1, 2, "A", 10)])
    search(
        tmp_path,
        "excluded",
        {"B": P1 + P2},
        [(P1, "B", 1, 0.001), (P2, "B", 1, 0.001)],
        [(P2, 2, "B", 20)],
    )
    result = group_searches(tmp_path / "search", tmp_path / "out", acquisition_names=["included"])
    assert result["acquisitions"] == 1
    assert result["acquisition_names"] == ["included"]
    assert result["searched_universe"] == 1
    assert table(tmp_path / "out/peptide_to_study_group.tsv") == [
        {"canonical_peptide": P1, "study_group": "A"}
    ]
    paths = [r["path"] for r in json.loads((tmp_path / "out/input_manifest.json").read_text())]
    assert all(
        "excluded" not in Path(path).parts and Path(path).name != "excluded.fasta" for path in paths
    )


@pytest.mark.parametrize("names", [[], ["S", "S"], ["missing"], ["../S"]])
def test_invalid_explicit_acquisition_population_fails(tmp_path, names):
    search(tmp_path, "S", {"A": P1}, [(P1, "A", 1, 0.001)], [(P1, 2, "A", 10)])
    with pytest.raises(ValueError, match="[Aa]cquisition"):
        group_searches(tmp_path / "search", tmp_path / "out", acquisition_names=names)


def test_output_cannot_enter_search_tree_via_parent_symlink(tmp_path):
    search(tmp_path, "S", {"A": P1}, [(P1, "A", 1, 0.001)], [(P1, 2, "A", 10)])
    (tmp_path / "alias").symlink_to(tmp_path / "search", target_is_directory=True)
    with pytest.raises(ValueError, match="outside the search input"):
        group_searches(tmp_path / "search", tmp_path / "alias/new")
