"""T7-T10: Synthetic unit tests for the 4 inference strategies."""

from __future__ import annotations

from fasta_lake.inference.base import InferenceInput
from fasta_lake.inference.strategies import (
    infer_razor,
    infer_species_budget,
    infer_uniform,
    infer_uniform_2pep,
)


def _make_synthetic_input() -> InferenceInput:
    """
    Create a synthetic 5-protein, 2-species scenario with 3 peptide types.

    Proteins:
        MGYG000000001_00001 (species A) - has 2 unique peptides + 1 shared
        MGYG000000001_00002 (species A) - has 1 shared peptide only
        MGYG000000002_00001 (species B) - has 1 unique peptide + 1 cross-species shared
        SMORF_00001           (smORF)    - has 1 unique peptide
        SMORF_00002           (smORF)    - has 1 shared peptide only

    Peptides:
        PEP_UNIQUE_A1 → only MGYG000000001_00001
        PEP_UNIQUE_A2 → only MGYG000000001_00001
        PEP_UNIQUE_B  → only MGYG000000002_00001
        PEP_UNIQUE_S  → only SMORF_00001
        PEP_WITHIN_A  → MGYG000000001_00001 + MGYG000000001_00002 (within species A)
        PEP_CROSS_AB  → MGYG000000001_00002 + MGYG000000002_00001 (cross species)
        PEP_SHARED_SM → SMORF_00001 + SMORF_00002 (within smORF)
    """
    pep2prot = {
        "PEP_UNIQUE_A1": {"MGYG000000001_00001"},
        "PEP_UNIQUE_A2": {"MGYG000000001_00001"},
        "PEP_UNIQUE_B": {"MGYG000000002_00001"},
        "PEP_UNIQUE_S": {"SMORF_00001"},
        "PEP_WITHIN_A": {"MGYG000000001_00001", "MGYG000000001_00002"},
        "PEP_CROSS_AB": {"MGYG000000001_00002", "MGYG000000002_00001"},
        "PEP_SHARED_SM": {"SMORF_00001", "SMORF_00002"},
    }

    prot2pep = {
        "MGYG000000001_00001": {"PEP_UNIQUE_A1", "PEP_UNIQUE_A2", "PEP_WITHIN_A"},
        "MGYG000000001_00002": {"PEP_WITHIN_A", "PEP_CROSS_AB"},
        "MGYG000000002_00001": {"PEP_UNIQUE_B", "PEP_CROSS_AB"},
        "SMORF_00001": {"PEP_UNIQUE_S", "PEP_SHARED_SM"},
        "SMORF_00002": {"PEP_SHARED_SM"},
    }

    proteins = {pid: "MFAKESEQ" * 10 for pid in prot2pep}
    # smORF ids carry no species; the map places both in one genome.
    taxonomy_map = {"SMORF_00001": "SMORF_GENOME", "SMORF_00002": "SMORF_GENOME"}

    return InferenceInput(
        pep2prot=pep2prot, prot2pep=prot2pep, proteins=proteins, taxonomy_map=taxonomy_map
    )


# ── T7: Species Budget ──────────────────────────────────────────────────────


class TestSpeciesBudget:
    """T7: Species budget with synthetic multi-species data."""

    def test_selected_proteins(self):
        inp = _make_synthetic_input()
        result = infer_species_budget(inp)

        # MGYG000000001_00001: has unique peptides → selected
        assert "MGYG000000001_00001" in result.selected_proteins

        # MGYG000000002_00001: has unique peptide → selected
        assert "MGYG000000002_00001" in result.selected_proteins

        # SMORF_00001: has unique peptide → selected
        assert "SMORF_00001" in result.selected_proteins

    def test_within_species_shared_goes_to_best(self):
        inp = _make_synthetic_input()
        result = infer_species_budget(inp)

        # PEP_WITHIN_A shared within species A:
        # MGYG000000001_00001 has 2 unique peptides, _00002 has 0
        # → razor to _00001
        # _00002 only has PEP_WITHIN_A and PEP_CROSS_AB
        # _00002 should still get PEP_CROSS_AB (cross-species)
        assert "MGYG000000001_00001" in result.selected_proteins

    def test_cross_species_keeps_both_species(self):
        inp = _make_synthetic_input()
        result = infer_species_budget(inp)

        # PEP_CROSS_AB is cross-species between A and B
        # Species budget should keep best from EACH species
        # Species A best = _00002 (only one in species A with this peptide)
        # Species B best = _00001 (only one in species B)
        assert "MGYG000000001_00002" in result.selected_proteins
        assert "MGYG000000002_00001" in result.selected_proteins

    def test_stats_structure(self):
        inp = _make_synthetic_input()
        result = infer_species_budget(inp)

        assert result.stats["strategy"] == "species_budget"
        assert result.stats["unique_peptides"] == 4
        assert result.stats["shared_peptides"] == 3
        assert result.stats["species_detected"] >= 2
        assert "within_species_shared" in result.stats
        assert "cross_species_shared" in result.stats


# ── T8: Razor ────────────────────────────────────────────────────────────────


class TestRazor:
    """T8: Shared peptide → protein with most unique peptides."""

    def test_shared_goes_to_most_unique(self):
        inp = _make_synthetic_input()
        result = infer_razor(inp)

        # MGYG000000001_00001 has 2 unique peptides — should win razor
        assert "MGYG000000001_00001" in result.selected_proteins
        assert "SMORF_00001" in result.selected_proteins

    def test_proteins_without_razor_win_excluded(self):
        inp = _make_synthetic_input()
        result = infer_razor(inp)

        # SMORF_00002 has only PEP_SHARED_SM and 0 unique peptides
        # SMORF_00001 has 1 unique → wins the shared peptide
        # So SMORF_00002 gets nothing
        assert "SMORF_00002" not in result.selected_proteins

    def test_razor_fewer_than_uniform(self):
        inp = _make_synthetic_input()
        razor_result = infer_razor(inp)
        uniform_result = infer_uniform(inp)

        assert len(razor_result.selected_proteins) <= len(uniform_result.selected_proteins)

    def test_stats_structure(self):
        inp = _make_synthetic_input()
        result = infer_razor(inp)
        assert result.stats["strategy"] == "razor"
        assert result.stats["proteins_with_unique_peptides"] > 0


# ── T9: Uniform ─────────────────────────────────────────────────────────────


class TestUniform:
    """T9: All evidence proteins kept."""

    def test_selects_all_with_evidence(self):
        inp = _make_synthetic_input()
        result = infer_uniform(inp)

        assert result.selected_proteins == set(inp.prot2pep.keys())
        assert len(result.selected_proteins) == 5

    def test_stats_correct(self):
        inp = _make_synthetic_input()
        result = infer_uniform(inp)
        assert result.stats["selected_proteins"] == 5
        assert result.stats["strategy"] == "uniform"


# ── T10: Uniform 2-peptide ──────────────────────────────────────────────────


class TestUniform2Pep:
    """T10: 1-pep excluded, 2-pep kept."""

    def test_filters_by_peptide_count(self):
        inp = _make_synthetic_input()
        result = infer_uniform_2pep(inp, min_peptides=2)

        # Proteins with >= 2 peptides:
        # MGYG000000001_00001: 3 peptides → kept
        # MGYG000000001_00002: 2 peptides → kept
        # MGYG000000002_00001: 2 peptides → kept
        # SMORF_00001: 2 peptides → kept
        # SMORF_00002: 1 peptide → EXCLUDED
        assert "MGYG000000001_00001" in result.selected_proteins
        assert "MGYG000000001_00002" in result.selected_proteins
        assert "MGYG000000002_00001" in result.selected_proteins
        assert "SMORF_00001" in result.selected_proteins
        assert "SMORF_00002" not in result.selected_proteins

    def test_excluded_count(self):
        inp = _make_synthetic_input()
        result = infer_uniform_2pep(inp, min_peptides=2)
        assert result.stats["excluded_proteins"] == 1

    def test_min_3_peptides(self):
        inp = _make_synthetic_input()
        result = infer_uniform_2pep(inp, min_peptides=3)

        # Only MGYG000000001_00001 has 3 peptides
        assert result.selected_proteins == {"MGYG000000001_00001"}

    def test_uniform_2pep_subset_of_uniform(self):
        inp = _make_synthetic_input()
        uniform = infer_uniform(inp)
        two_pep = infer_uniform_2pep(inp, min_peptides=2)
        assert two_pep.selected_proteins.issubset(uniform.selected_proteins)
