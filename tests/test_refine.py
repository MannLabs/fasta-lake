"""refine.py is proposal-only: it must say so, behave so, and refuse to pretend otherwise."""

from __future__ import annotations

from pathlib import Path

import pytest

from fasta_lake import refine
from fasta_lake.blosum import generate_variants

PROTEIN = "MKTAYIAKQRQISFVKSHFSRQ"
PEPTIDE = "QISFVK"


def _write_fasta(path, records):
    with open(path, "w") as fh:
        for acc, seq in records:
            fh.write(f">{acc}\n{seq}\n")


def _read_fasta(path):
    records = {}
    acc = None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            acc = line[1:]
            records[acc] = ""
        elif acc is not None:
            records[acc] += line
    return records


class TestGenerateVariantFasta:
    def test_modified_record_is_single_substitution_and_renamed(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        _write_fasta(fasta, [("P1", PROTEIN), ("P2", "MMMMMMMMMMMM")])
        out = tmp_path / "out.fasta"

        stats = refine.generate_variant_fasta(
            fasta, [{"peptide": PEPTIDE, "protein_id": "P1"}], output_path=out
        )

        records = _read_fasta(out)
        assert set(records) == {"P1__VAR", "P2"}, "same record count, modified one renamed"
        assert records["P2"] == "MMMMMMMMMMMM", "untouched protein is byte-identical"
        var = records["P1__VAR"]
        assert len(var) == len(PROTEIN)
        diff = [i for i, (a, b) in enumerate(zip(PROTEIN, var)) if a != b]
        assert len(diff) == 1, "exactly one substitution"
        start = PROTEIN.index(PEPTIDE)
        assert start <= diff[0] < start + len(PEPTIDE), "substitution lies inside the peptide"
        assert stats == {
            "candidates": 1,
            "mapped": 1,
            "proteins_modified": 1,
            "total_proteins": 2,
            "matrix": "blosum95",
        }

    def test_best_scoring_substitution_is_chosen(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        _write_fasta(fasta, [("P1", PROTEIN)])
        out = tmp_path / "out.fasta"
        refine.generate_variant_fasta(
            fasta, [{"peptide": PEPTIDE, "protein_id": "P1"}], output_path=out
        )
        start = PROTEIN.index(PEPTIDE)
        best = max(
            (
                v
                for pos in range(start, start + len(PEPTIDE))
                for v in generate_variants(PROTEIN, pos)
            ),
            key=lambda v: v[2],
        )
        # The module keeps the first variant reaching the maximum score, so the
        # chosen sequence must carry that maximal score.
        chosen = _read_fasta(out)["P1__VAR"]
        chosen_score = next(
            v[2]
            for pos in range(start, start + len(PEPTIDE))
            for v in generate_variants(PROTEIN, pos)
            if v[0] == chosen
        )
        assert chosen_score == best[2]

    def test_unmapped_candidates_leave_fasta_unchanged(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        _write_fasta(fasta, [("P1", PROTEIN)])
        out = tmp_path / "out.fasta"
        stats = refine.generate_variant_fasta(
            fasta,
            [
                {"peptide": "WWWWWW", "protein_id": "P1"},  # not in sequence
                {"peptide": PEPTIDE, "protein_id": "MISSING"},  # unknown accession
            ],
            output_path=out,
        )
        assert _read_fasta(out) == {"P1": PROTEIN}
        assert stats["mapped"] == 0 and stats["proteins_modified"] == 0

    def test_il_folded_match_is_found(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        _write_fasta(fasta, [("P1", PROTEIN)])
        stats = refine.generate_variant_fasta(
            fasta, [{"peptide": PEPTIDE.replace("I", "L"), "protein_id": "P1"}]
        )
        assert stats["mapped"] == 1

    def test_output_is_deterministic(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        _write_fasta(fasta, [("P1", PROTEIN), ("P0", PROTEIN)])
        cands = [
            {"peptide": PEPTIDE, "protein_id": "P1"},
            {"peptide": PEPTIDE, "protein_id": "P0"},
        ]
        a, b = tmp_path / "a.fasta", tmp_path / "b.fasta"
        refine.generate_variant_fasta(fasta, cands, output_path=a)
        refine.generate_variant_fasta(fasta, list(reversed(cands)), output_path=b)
        assert a.read_bytes() == b.read_bytes()


class TestProposalOnly:
    def test_run_refinement_refuses(self, tmp_path):
        with pytest.raises(NotImplementedError, match="proposal-only"):
            refine.run_refinement(tmp_path / "r.tsv", tmp_path / "s.fasta", tmp_path / "o.fasta")

    def test_module_docstring_states_status(self):
        assert "proposal-only" in refine.__doc__
        assert "does not validate" in refine.__doc__
