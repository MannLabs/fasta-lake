"""A renamed input column must be an error naming the column, never an empty result."""

from __future__ import annotations

import pytest

from fasta_lake.aggregate import load_diann_sample
from fasta_lake.benchmark import load_diann_report, load_sage_report
from fasta_lake.canonical_comparison import _load_diann_gene_peptides
from fasta_lake.helpers import collect_peptides_with_scores, require_columns
from fasta_lake.tune import parse_tsv

DIANN = "File.Name\tStripped.Sequence\tGenes\tProtein.Group\tQ.Value\tPrecursor.Quantity\n"
DIANN_ROW = "s.mzML\tAAVLIDK\tALB\tP02768\t0.001\t1000\n"


def test_require_columns_names_missing_and_found():
    with pytest.raises(ValueError, match=r"missing required column\(s\) b, c or d; found: a, x"):
        require_columns(["a", "x"], ("a", "b"), "t.tsv", any_of=("c", "d"))
    require_columns(["a", "b", "c"], ("a", "b"), "t.tsv", any_of=("c", "d"))  # no raise


class TestPredictionCsv:
    def test_alphanovo_and_generic_columns_both_accepted(self, tmp_path):
        a = tmp_path / "a.csv"
        a.write_text("peptide_prediction_detokenized_unmodified,score\nPEPTIDEKAAAR,0.9\n")
        g = tmp_path / "g.csv"
        g.write_text("sequence,score\nPEPTIDEKAAAR,0.9\n")
        assert collect_peptides_with_scores(a) == collect_peptides_with_scores(g)
        assert collect_peptides_with_scores(a) == {"PEPTLDEKAAAR": 0.9}

    def test_renamed_peptide_column_is_refused(self, tmp_path):
        p = tmp_path / "p.csv"
        p.write_text("peptide,score\nPEPTIDEKAAAR,0.9\n")
        with pytest.raises(
            ValueError, match="peptide_prediction_detokenized_unmodified or sequence"
        ):
            collect_peptides_with_scores(p)

    def test_missing_score_column_is_refused(self, tmp_path):
        p = tmp_path / "p.csv"
        p.write_text("sequence,confidence\nPEPTIDEKAAAR,0.9\n")
        with pytest.raises(ValueError, match=r"missing required column\(s\) score"):
            collect_peptides_with_scores(p)

    def test_unparsable_score_is_skipped_not_defaulted(self, tmp_path):
        p = tmp_path / "p.csv"
        p.write_text("sequence,score\nPEPTIDEKAAAR,n/a\nANOTHERPEPK,0.7\n")
        assert collect_peptides_with_scores(p) == {"ANOTHERPEPK": 0.7}


@pytest.mark.parametrize(
    "loader",
    [load_diann_sample, load_diann_report, lambda p: _load_diann_gene_peptides(p, 0.01)],
    ids=["aggregate", "benchmark", "canonical_comparison"],
)
def test_diann_readers_refuse_renamed_q_column(tmp_path, loader):
    ok = tmp_path / "ok.tsv"
    ok.write_text(DIANN + DIANN_ROW)
    loader(ok)  # complete header: fine
    bad = tmp_path / "bad.tsv"
    bad.write_text((DIANN + DIANN_ROW).replace("Q.Value", "QValue"))
    with pytest.raises(ValueError, match=r"missing required column\(s\) Q\.Value"):
        loader(bad)


def test_sage_benchmark_reader_needs_a_q_column(tmp_path):
    bad = tmp_path / "r.tsv"
    bad.write_text("peptide\tproteins\tlabel\tq\nAAVLIDK\tP1\t1\t0.001\n")
    with pytest.raises(ValueError, match="spectrum_q or peptide_q"):
        load_sage_report(bad)


def test_tune_parser_refuses_missing_score_column(tmp_path):
    bad = tmp_path / "r.tsv"
    bad.write_text("spectrum_q\tpeptide_len\tproteins\n0.001\t10\tP1\n")
    with pytest.raises(ValueError, match="sage_discriminant_score"):
        parse_tsv(bad)
