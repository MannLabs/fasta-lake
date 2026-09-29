"""Invalid q values and lower-ranked hits must not enter control counts."""
import pytest

from fasta_lake.entrapment_workflow import classify_psms


@pytest.mark.parametrize("q", ["NaN", "-inf", "-0.1", "inf", "nonsense"])
def test_invalid_confidence_does_not_enter_either_partition(tmp_path, q):
    path = tmp_path / "results.sage.tsv"
    path.write_text(
        "peptide\tproteins\tlabel\tpeptide_q\trank\n"
        f"PEPTIDEAK\tTARGET\t1\t{q}\t1\n"
        f"OTHERPEPKK\tENTRAP_TEST\t1\t{q}\t1\n"
    )
    result = classify_psms(path)
    assert result.target == result.entrap == result.shared == 0


def test_higher_ranks_do_not_change_control_ratio(tmp_path):
    path = tmp_path / "results.sage.tsv"
    path.write_text(
        "peptide\tproteins\tlabel\tpeptide_q\trank\n"
        "PEPTIDEAK\tTARGET\t1\t0.001\t1\n"
        "OTHERPEPKK\tENTRAP_TEST\t1\t0.001\t2\n"
    )
    result = classify_psms(path)
    assert result.target == 1 and result.entrap == 0
