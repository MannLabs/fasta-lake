"""Worked scientific grouping examples, including ambiguity and conservation."""

import json

import pytest

from fasta_lake.study import group_searches, member_set_key, study_dictionary


def _write_sage(sample, text):
    lines = text.splitlines()
    lines[0] += "\tfilename"
    lines[1:] = [line + "\t" + sample.name + ".mzML" for line in lines[1:]]
    (sample / "results.sage.tsv").write_text("\n".join(lines) + "\n")


def _write_lfq(sample, text):
    text = text.replace("\tintensity\n", "\t" + sample.name + ".mzML\n")
    (sample / "lfq.tsv").write_text(text)


def test_member_sets_are_not_study_groups():
    assert member_set_key("B;A;B") == "A;B"
    result = study_dictionary({"A": {"p1"}, "B": {"p1", "p2"}, "C": {"p2"}}, {"A", "B", "C", "D"})
    assert result["selected"] == [("B", 2)]
    assert result["peptide_to_group"] == {"p1": "B", "p2": "B"}
    assert result["unobserved"] == {"D"}


def test_study_grouping_rejects_unsearched_representatives():
    with pytest.raises(ValueError, match="outside"):
        study_dictionary({"A": {"p1"}}, {"B"})


def test_greedy_is_deterministic_and_covers_every_observed_peptide():
    graph = {"C": {"p2", "p3"}, "A": {"p1", "p2"}, "B": {"p1"}}
    a = study_dictionary(graph, set(graph))
    b = study_dictionary(dict(reversed(list(graph.items()))), set(graph))
    assert a == b
    assert a["selected"] == [("A", 2), ("C", 1)]
    assert a["covered"] == {"p1", "p2", "p3"}


def test_direct_groups_conserve_quantification_and_distinguish_features(tmp_path):
    search = tmp_path / "search"
    sample = search / "S"
    sample.mkdir(parents=True)
    fasta = tmp_path / "database.fasta"
    fasta.write_text(">A\nPEPTIDEAK\n>B\nPEPTIDEAKOTHERPEPK\n")
    (sample / "config.json").write_text(
        json.dumps({"database": {"fasta": str(fasta)}, "mzml_paths": [sample.name + ".mzML"]})
    )
    _write_sage(
        sample,
        "peptide\tproteins\tlabel\tpeptide_q\nPEPTIDEAK\tA;B\t1\t0.001\nOTHERPEPK\tB\t1\t0.001\n",
    )
    _write_lfq(
        sample,
        "peptide\tcharge\tproteins\tq_value\tscore\tspectral_angle\tS.mzML\n"
        "PEPTIDEAK\t2\tA;B\t0.5\t1\t0.9\t10\n"
        "PEPTIDEAK\t3\tA;B\t0.5\t1\t0.9\t20\n"
        "OTHERPEPK\t2\tB\t0.5\t1\t0.9\t40\n",
    )
    result = group_searches(search, tmp_path / "grouped")
    assert result["inferred_study_groups"] == 1
    assert result["intensity_in"] == result["intensity_out"] == 70
    text = (tmp_path / "grouped/study_group_long.tsv").read_text()
    assert "S\tB\t2\t3\t70.0" in text
    assert (tmp_path / "grouped/representatives.fasta").read_text() == ">B\nPEPTIDEAKOTHERPEPK\n"
    assert result["annotation_targets"]["fasta"] == "representatives.fasta"


def test_each_acquisition_must_use_its_own_searched_accessions(tmp_path):
    """A protein searched elsewhere must not legitimize a mismatched sample table."""
    search = tmp_path / "search"
    for sample_name, accession in [("S1", "A"), ("S2", "B")]:
        sample = search / sample_name
        sample.mkdir(parents=True)
        fasta = tmp_path / (sample_name + ".fasta")
        fasta.write_text(f">{accession}\nPEPTIDEAK\n")
        (sample / "config.json").write_text(
            json.dumps({"database": {"fasta": str(fasta)}, "mzml_paths": [sample.name + ".mzML"]})
        )
        # S1 is deliberately cross-wired to S2's protein. B exists globally.
        _write_sage(sample, "peptide\tproteins\tlabel\tpeptide_q\nPEPTIDEAK\tB\t1\t0.001\n")
        _write_lfq(
            sample,
            "peptide\tcharge\tproteins\tq_value\tscore\tspectral_angle\tintensity\n"
            "PEPTIDEAK\t2\tB\t0.5\t1\t0.9\t10\n",
        )
    with pytest.raises(ValueError, match="S1.*searched FASTA"):
        group_searches(search, tmp_path / "grouped")


def test_finite_features_cannot_overflow_into_an_infinite_aggregate(tmp_path):
    sample = tmp_path / "search/S"
    sample.mkdir(parents=True)
    fasta = tmp_path / "database.fasta"
    fasta.write_text(">A\nPEPTIDEAK\n")
    (sample / "config.json").write_text(
        json.dumps({"database": {"fasta": str(fasta)}, "mzml_paths": [sample.name + ".mzML"]})
    )
    _write_sage(sample, "peptide\tproteins\tlabel\tpeptide_q\nPEPTIDEAK\tA\t1\t0.001\n")
    _write_lfq(
        sample,
        "peptide\tcharge\tproteins\tq_value\tscore\tspectral_angle\tintensity\n"
        "PEPTIDEAK\t2\tA\t0.5\t1\t0.9\t1e308\n"
        "PEPTIDEAK\t3\tA\t0.5\t1\t0.9\t1e308\n",
    )
    with pytest.raises(ValueError, match="Non-finite aggregate intensity"):
        group_searches(tmp_path / "search", tmp_path / "grouped")
    assert not (tmp_path / "grouped").exists()


def test_quantified_accessions_must_belong_to_the_acquisition(tmp_path):
    """A valid peptide gate cannot rescue an LFQ row from a different database."""
    search = tmp_path / "search"
    for sample_name, accession in [("S1", "A"), ("S2", "B")]:
        sample = search / sample_name
        sample.mkdir(parents=True)
        fasta = tmp_path / (sample_name + ".fasta")
        fasta.write_text(f">{accession}\nPEPTIDEAK\n")
        (sample / "config.json").write_text(
            json.dumps({"database": {"fasta": str(fasta)}, "mzml_paths": [sample.name + ".mzML"]})
        )
        _write_sage(
            sample, f"peptide\tproteins\tlabel\tpeptide_q\nPEPTIDEAK\t{accession}\t1\t0.001\n"
        )
        _write_lfq(
            sample,
            "peptide\tcharge\tproteins\tq_value\tscore\tspectral_angle\tintensity\n"
            "PEPTIDEAK\t2\tB\t0.5\t1\t0.9\t10\n",
        )
    with pytest.raises(ValueError, match="S1.*searched FASTA"):
        group_searches(search, tmp_path / "grouped")
