"""Readable study groups: gene, protein name and organism from the lake source headers."""

import csv
import json

import pytest

from fasta_lake.study_annotation import annotate_study, load_header_sidecar
from fasta_lake.study_cli import main as group_study

P1, P2 = "AAAAAAK", "CCCCCCK"
SP = "sp|P0A7V0|RS2_ECOLI 30S ribosomal protein S2 OS=Escherichia coli (strain K12) GN=rpsB"
TR = "tr|Q8A1|Q8A1_BACTN Elongation factor Tu OS=Bacteroides thetaiotaomicron GN=tuf"
TR2 = "tr|Q9X2|Q9X2_BACFR Elongation factor Tu OS=Bacteroides fragilis GN=tuf"


def write_search(tmp_path, headers):
    """One acquisition whose two proteins each explain one accepted peptide."""
    folder = tmp_path / "search" / "S"
    folder.mkdir(parents=True)
    fasta = tmp_path / "S.fasta"
    fasta.write_text(f">{headers[0]}\n{P1}\n>{headers[1]}\n{P2}\n")
    accessions = [h.split()[0] for h in headers]
    (folder / "config.json").write_text(
        json.dumps({"database": {"fasta": str(fasta)}, "mzml_paths": ["S.mzML"]})
    )
    with (folder / "results.sage.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["peptide", "proteins", "label", "peptide_q", "filename", "rank"])
        w.writerow([P1, accessions[0], 1, 0.001, "S.mzML", 1])
        w.writerow([P2, accessions[1], 1, 0.001, "S.mzML", 1])
    with (folder / "lfq.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            ["peptide", "charge", "proteins", "q_value", "score", "spectral_angle", "S.mzML"]
        )
        w.writerow([P1, 2, accessions[0], 0.5, 1, 0.9, 10])
        w.writerow([P2, 2, accessions[1], 0.5, 1, 0.9, 20])
    return accessions


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_hash_accessions_are_named_from_the_lake_sidecar(tmp_path):
    accessions = write_search(tmp_path, ["LAKE_aaa", "LAKE_bbb"])
    sidecar = tmp_path / "headers.tsv"
    sidecar.write_text(
        "sha256_hash\tn_sources\tkept_header\tall_source_headers\n"
        f"aaa\t1\tLAKE_aaa\tUNIPROT {SP}\n"
        # One sequence carried by two catalogue entries keeps both organisms.
        f"bbb\t2\tLAKE_bbb\tUNIPROT {TR}|||UNIPROT {TR2}\n"
        "ccc\t1\tLAKE_ccc\tUNIPROT sp|X|Y_HUMAN unrelated OS=Homo sapiens\n"
    )
    assert set(load_header_sidecar(sidecar, set(accessions))) == set(accessions)
    status = group_study(
        [
            "--search",
            str(tmp_path / "search"),
            "--out",
            str(tmp_path / "groups"),
            "--lake-headers",
            str(sidecar),
            "--skip-qc",
        ]
    )
    assert status == 0
    rows = {r["study_group"]: r for r in read(tmp_path / "groups/study_group_annotation.tsv")}
    assert rows["LAKE_aaa"]["gene"] == "rpsB"
    assert rows["LAKE_aaa"]["protein_name"] == "30S ribosomal protein S2"
    assert rows["LAKE_aaa"]["organism"] == "Escherichia coli (strain K12)"
    assert rows["LAKE_bbb"]["gene"] == "tuf"
    assert rows["LAKE_bbb"]["organism"] == ("Bacteroides thetaiotaomicron; Bacteroides fragilis")
    assert rows["LAKE_bbb"]["n_source_headers"] == "2"
    matrix = read(tmp_path / "groups/annotated_study_group_matrix.tsv")
    assert [r["study_group"] for r in matrix] == ["LAKE_aaa", "LAKE_bbb"]
    assert matrix[0]["gene"] == "rpsB" and matrix[0]["S"] == "10.0"
    # The unannotated matrix is unchanged.
    assert read(tmp_path / "groups/study_group_matrix.tsv")[0] == {
        "study_group": "LAKE_aaa",
        "S": "10.0",
    }
    record = json.loads((tmp_path / "groups/WORKFLOW.json").read_text())
    assert record["annotation"]["named_study_groups"] == 2


def test_without_sidecar_the_searched_fasta_headers_name_the_groups(tmp_path):
    write_search(tmp_path, [SP, TR])
    assert (
        group_study(
            ["--search", str(tmp_path / "search"), "--out", str(tmp_path / "groups"), "--skip-qc"]
        )
        == 0
    )
    rows = {r["study_group"]: r for r in read(tmp_path / "groups/study_group_annotation.tsv")}
    assert rows["sp|P0A7V0|RS2_ECOLI"]["organism"] == "Escherichia coli (strain K12)"
    assert rows["tr|Q8A1|Q8A1_BACTN"]["protein_name"] == "Elongation factor Tu"


def test_unnamed_groups_stay_blank_rather_than_invented(tmp_path):
    write_search(tmp_path, ["LAKE_aaa", "LAKE_bbb"])
    assert (
        group_study(
            ["--search", str(tmp_path / "search"), "--out", str(tmp_path / "groups"), "--skip-qc"]
        )
        == 0
    )
    for row in read(tmp_path / "groups/study_group_annotation.tsv"):
        assert row["gene"] == row["protein_name"] == row["organism"] == ""
        assert row["n_source_headers"] == "0"


def test_a_file_that_is_not_a_sidecar_is_rejected(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("accession\theader\n")
    with pytest.raises(ValueError, match="not a lake_builder"):
        load_header_sidecar(bad)
    write_search(tmp_path, [SP, TR])
    groups = tmp_path / "groups"
    from fasta_lake.study import group_searches

    group_searches(tmp_path / "search", groups)
    with pytest.raises(ValueError, match="not a lake_builder"):
        annotate_study(groups, headers_tsv=bad)
