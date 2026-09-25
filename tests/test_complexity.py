"""Tests for sample complexity estimation (v3 feature)."""

from __future__ import annotations

from fasta_lake.complexity import estimate_complexity
from fasta_lake.config import ComplexityConfig


class TestEstimateComplexity:
    """Test complexity estimation with synthetic data."""

    def test_single_species_high_hit_rate(self):
        """Single-species sample with good coverage."""
        proteins = {
            "sp|P12345|GFP_HUMAN": "MSKGEELFTPEPTLDERLONGSEQACDE",
            "sp|P12346|UBQ_HUMAN": "ANOTHERPEPTLDESEQWLTHPEPTLDE",
        }
        peptides = {"PEPTLDER": 0.9, "ANOTHERPEPTLDE": 0.8}

        report = estimate_complexity(proteins, peptides)

        assert report.hit_rate > 0.0
        assert report.species_richness == 1  # All Human
        assert report.proteins_with_evidence > 0
        assert report.matched_peptides > 0

    def test_multi_species(self):
        """Multi-species sample should detect multiple species."""
        # Use a shared peptide that exists in proteins from different species
        # Note: I/L normalization is applied, so use L not I
        shared_pep = "SHAREDPEPTLDE"
        proteins = {
            "sp|P12345|GFP_HUMAN": f"MSK{shared_pep}RLONGSEQ",
            "MGYG000000001_00001": f"ANOTHER{shared_pep}SEQ",
            "MGYG000000002_00001": f"THLRD{shared_pep}SEQRES",
        }
        peptides = {shared_pep: 0.9}

        report = estimate_complexity(proteins, peptides)

        # Should detect 3 species: Human + 2 MGYG genomes
        assert report.species_richness == 3

    def test_empty_peptides(self):
        """Empty peptides -> zero metrics."""
        proteins = {"p1": "MSKGEELFT"}
        report = estimate_complexity(proteins, set())
        assert report.hit_rate == 0.0
        assert report.species_richness == 0
        assert report.total_peptides == 0

    def test_empty_proteins(self):
        """Empty protein DB -> zero metrics."""
        report = estimate_complexity({}, {"PEPTLDER": 0.9})
        assert report.hit_rate == 0.0
        assert report.db_coverage == 0.0

    def test_no_matches(self):
        """Peptides don't match any protein."""
        proteins = {"p1": "MSKGEELFT"}
        peptides = {"ZZZZZZZZZZZ": 0.9}

        report = estimate_complexity(proteins, peptides)

        assert report.hit_rate == 0.0
        assert report.proteins_with_evidence == 0

    def test_auto_recommend_single_species(self):
        """Single species -> recommend razor."""
        proteins = {
            "sp|P12345|GFP_HUMAN": "MSKPEPTLDERLONGSEQ",
            "sp|P12346|UBQ_HUMAN": "ANOTHERPEPTLDESEQ",
        }
        peptides = {"PEPTLDER": 0.9}
        config = ComplexityConfig(auto_recommend=True)

        report = estimate_complexity(proteins, peptides, config)

        assert report.recommendation == "razor"
        assert (
            "single" in report.recommendation_reason.lower()
            or "species" in report.recommendation_reason.lower()
        )

    def test_auto_recommend_disabled(self):
        """No recommendation when auto_recommend=False."""
        proteins = {"sp|P12345|GFP_HUMAN": "MSKPEPTLDERLONGSEQ"}
        peptides = {"PEPTLDER": 0.9}
        config = ComplexityConfig(auto_recommend=False)

        report = estimate_complexity(proteins, peptides, config)

        assert report.recommendation == ""

    def test_peptide_set_input(self):
        """Accept set[str] as peptide input."""
        proteins = {"p1": "MSKPEPTLDERLONGSEQ"}
        peptides = {"PEPTLDER"}

        report = estimate_complexity(proteins, peptides)

        assert report.total_peptides == 1

    def test_db_coverage(self):
        """DB coverage = fraction of proteins with evidence."""
        proteins = {
            "p1": "MSKPEPTLDERLONGSEQ",
            "p2": "ANOTHERPEPTLDESEQ",
            "p3": "NOPEPTLDEMATCH",
            "p4": "ZZZZZZZZZZZZZZZZZ",
        }
        peptides = {"PEPTLDER": 0.9}

        report = estimate_complexity(proteins, peptides)

        # p1, p2, p3 match (PEPTLDE is in all 3), p4 doesn't
        assert 0.0 < report.db_coverage <= 1.0

    def test_report_summary(self):
        """Summary dict has required keys."""
        proteins = {"p1": "MSKPEPTLDERLONGSEQ"}
        peptides = {"PEPTLDER": 0.9}

        report = estimate_complexity(proteins, peptides)
        summary = report.summary()

        assert "hit_rate" in summary
        assert "species_richness" in summary
        assert "peptide_diversity" in summary
        assert "db_coverage" in summary
        assert "total_peptides" in summary
        assert "matched_peptides" in summary

    def test_species_counts_populated(self):
        """Species counts track proteins per species."""
        proteins = {
            "MGYG000000001_00001": "MSKPEPTLDERLONGSEQ",
            "MGYG000000001_00002": "ANOTHERPEPTLDESEQ",
            "MGYG000000002_00001": "THEPEPTLDERLONGSEQ",
        }
        peptides = {"PEPTLDER": 0.9}

        report = estimate_complexity(proteins, peptides)

        assert "MGYG000000001" in report.species_counts


class TestAutoRecommend:
    """Test strategy auto-recommendation logic."""

    def _make_proteins(self, species_ids: list[str]) -> dict[str, str]:
        """Generate proteins for given species IDs."""
        proteins = {}
        for i, sid in enumerate(species_ids):
            pid = f"{sid}_{i:05d}" if not sid.startswith("sp|") else sid
            proteins[pid] = "MSKPEPTLDERLONGSEQ" + f"UNIQUE{i:05d}LONG"
        return proteins

    def test_low_complexity_recommends_razor(self):
        """<10 species -> razor."""
        species = [f"MGYG{i:09d}" for i in range(5)]
        protein_ids = [f"{s}_00001" for s in species]
        proteins = {pid: "MSKPEPTLDERLONGSEQ" for pid in protein_ids}
        peptides = {"PEPTLDER": 0.9}

        report = estimate_complexity(proteins, peptides)

        assert report.recommendation == "razor"

    def test_high_complexity_recommends_species_budget(self):
        """10-100 species -> species_budget."""
        protein_ids = [f"MGYG{i:09d}_00001" for i in range(50)]
        proteins = {pid: "MSKPEPTLDERLONGSEQ" for pid in protein_ids}
        peptides = {"PEPTLDER": 0.9}

        report = estimate_complexity(proteins, peptides)

        assert report.recommendation == "species_budget"

    def test_very_high_complexity_recommends_species_budget(self):
        """100+ species -> species_budget (with clustering note)."""
        protein_ids = [f"MGYG{i:09d}_00001" for i in range(150)]
        proteins = {pid: "MSKPEPTLDERLONGSEQ" for pid in protein_ids}
        peptides = {"PEPTLDER": 0.9}

        report = estimate_complexity(proteins, peptides)

        assert report.recommendation == "species_budget"
        assert "clustering" in report.recommendation_reason.lower()
