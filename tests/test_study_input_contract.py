"""Prevent external result tables from silently changing accepted evidence."""
import json

import pytest

from fasta_lake.study import group_searches


def inputs(tmp_path, *, rank="1", duplicate=False, duplicate_intensity="10"):
    sample = tmp_path / "search/S"
    sample.mkdir(parents=True)
    fasta = tmp_path / "database.fasta"
    fasta.write_text(">A\nPEPTIDEAKOTHERPEPK\n>B\nPEPTIDEAK\n")
    (sample / "config.json").write_text(json.dumps({
        "database": {"fasta": str(fasta)}, "mzml_paths": ["S.mzML"], "report_psms": 2,
    }))
    (sample / "results.sage.tsv").write_text(
        "peptide\tproteins\tlabel\tpeptide_q\tfilename\trank\n"
        "PEPTIDEAK\tA\t1\t0.001\tS.mzML\t1\n"
        f"OTHERPEPK\tA\t1\t0.001\tS.mzML\t{rank}\n"
        "PEPTIDEAK\tB\t1\t0.001\tS.mzML\t2\n"
    )
    lfq = "peptide\tcharge\tproteins\tS.mzML\nPEPTIDEAK\t2\tA\t10\n"
    lfq += "OTHERPEPK\t2\tA\t1000\n"
    if duplicate:
        lfq += f"PEPTIDEAK\t2\tA\t{duplicate_intensity}\n"
    (sample / "lfq.tsv").write_text(lfq)
    return sample


def test_rank_two_cannot_open_lfq_gate_or_expand_evidence(tmp_path):
    import gzip
    sample = inputs(tmp_path, rank="2")
    result = group_searches(sample.parent, tmp_path / "out")
    assert result["intensity_out"] == 10
    assert result["pooled_canonical_peptides"] == 1
    with gzip.open(tmp_path / "out/peptide_evidence.tsv.gz", "rt") as stream:
        evidence = stream.read()
    assert "\tA\t1\t1\n" in evidence
    assert "OTHERPEPK" not in evidence


@pytest.mark.parametrize("rank", ["", "0", "-1", "1.5", "NaN", "first"])
def test_malformed_rank_fails_before_output(tmp_path, rank):
    sample = inputs(tmp_path, rank=rank)
    with pytest.raises(ValueError, match="rank"):
        group_searches(sample.parent, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_missing_rank_is_not_assumed_to_be_best_hit(tmp_path):
    sample = inputs(tmp_path)
    path = sample / "results.sage.tsv"
    lines = [line.rsplit("\t", 1)[0] for line in path.read_text().splitlines()]
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="rank"):
        group_searches(sample.parent, tmp_path / "out")
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("intensity", ["10", "11"])
def test_duplicate_feature_fails_even_when_intensity_conservation_would_pass(tmp_path, intensity):
    sample = inputs(tmp_path, duplicate=True, duplicate_intensity=intensity)
    with pytest.raises(ValueError, match="Duplicate LFQ feature"):
        group_searches(sample.parent, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_distinct_modifications_and_charges_are_retained(tmp_path):
    sample = inputs(tmp_path)
    with (sample / "lfq.tsv").open("a") as stream:
        stream.write("PEPTIDEAK\t3\tA\t20\nPEPT[+1.0]IDEAK\t2\tA\t30\n")
    result = group_searches(sample.parent, tmp_path / "out")
    assert result["intensity_out"] == 1060


def test_charge_alias_cannot_bypass_duplicate_check(tmp_path):
    sample = inputs(tmp_path)
    with (sample / "lfq.tsv").open("a") as stream:
        stream.write("PEPTIDEAK\t02\tA\t10\n")
    with pytest.raises(ValueError, match="Duplicate LFQ feature"):
        group_searches(sample.parent, tmp_path / "out")


def test_combined_charge_sentinel_remains_supported(tmp_path):
    sample = inputs(tmp_path)
    path = sample / "lfq.tsv"
    path.write_text(path.read_text().replace("\t2\t", "\t-1\t"))
    assert group_searches(sample.parent, tmp_path / "out")["intensity_out"] == 1010
