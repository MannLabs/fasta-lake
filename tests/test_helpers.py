"""T1-T4: Unit tests for the 7 shared helper functions."""

from __future__ import annotations

from pathlib import Path

import pytest

from fasta_lake.helpers import (
    clean_sequence,
    collect_peptides_with_scores,
    get_db_type,
    get_species_id,
    load_fasta_proteins,
    map_peptides_to_proteins,
    normalize_il,
)

# ── T1: normalize_il ────────────────────────────────────────────────────────


class TestNormalizeIL:
    """T1: I→L replacement and no-change cases."""

    def test_replaces_i_with_l(self):
        assert normalize_il("PEPTIDER") == "PEPTLDER"

    def test_multiple_replacements(self):
        assert normalize_il("IIIII") == "LLLLL"

    def test_no_change_without_i(self):
        assert normalize_il("PEPTLDER") == "PEPTLDER"

    def test_empty_string(self):
        assert normalize_il("") == ""

    def test_preserves_other_chars(self):
        assert normalize_il("ACDEFGHIKLMNPQRSTVWY") == "ACDEFGHLKLMNPQRSTVWY"

    def test_mixed_case_only_upper_i(self):
        # Only uppercase I is replaced (sequences should be uppercased first)
        assert normalize_il("PiPTIDER") == "PiPTLDER"


# ── T2: clean_sequence ──────────────────────────────────────────────────────


class TestCleanSequence:
    """T2: Uppercase, remove *, strip whitespace."""

    def test_uppercase(self):
        assert clean_sequence("mskgeelft") == "MSKGEELFT"

    def test_remove_stop_codon(self):
        assert clean_sequence("MSKGE*ELFT*") == "MSKGEELFT"

    def test_mixed_case_with_star(self):
        assert clean_sequence("Msk*Ge") == "MSKGE"

    def test_empty_string(self):
        assert clean_sequence("") == ""

    def test_already_clean(self):
        assert clean_sequence("MSKGEELFT") == "MSKGEELFT"


# ── T3: get_db_type ─────────────────────────────────────────────────────────


class TestGetDbType:
    """T3: All 5 database types."""

    def test_human_swissprot(self):
        assert get_db_type("sp|P12345|GFP_HUMAN") == "Human"

    def test_human_trembl(self):
        assert get_db_type("tr|Q99999|UBQ_HUMAN") == "Human"

    def test_uhgg(self):
        assert get_db_type("MGYG000000001_00001") == "UHGG"

    def test_gmgc(self):
        assert get_db_type("GMGC10.001_001_001.PROT01") == "GMGC"

    def test_assembly(self):
        assert get_db_type("MuPr_Assembly_SAMP1_00001") == "Assembly"

    def test_smorf(self):
        assert get_db_type("SMORF_00001") == "smORF"

    def test_unknown_prefix(self):
        assert get_db_type("UNKNOWN_PROTEIN") == "smORF"


# ── T4: get_species_id ──────────────────────────────────────────────────────


class TestGetSpeciesId:
    """T4: Species extraction for all protein ID formats."""

    def test_uhgg_genome(self):
        assert get_species_id("MGYG000000001_00001") == "MGYG000000001"

    def test_uhgg_genome_different(self):
        assert get_species_id("MGYG000000002_00001") == "MGYG000000002"

    def test_gmgc_unigene_ids_carry_no_species(self):
        # A unigene id names one gene, not an organism: never a pseudo-species.
        assert get_species_id("GMGC10.000_001_411.RBR") is None
        assert get_species_id("GMGC10.001.002.003.GENE") is None

    def test_uniprot_entry_name_organism(self):
        assert get_species_id("sp|P12345|GFP_HUMAN") == "HUMAN"
        assert get_species_id("tr|Q99999|UBQ_MOUSE") == "MOUSE"
        assert get_species_id("sp|P12345|GFP_HUMAN Green fluorescent") == "HUMAN"

    def test_uniprot_without_entry_name_is_unresolved(self):
        assert get_species_id("sp|P12345") is None

    def test_lake_headers_resolve_through_the_source_tag(self):
        assert get_species_id("uh|MGYG000000001_00001|MGYG000000001_00001") == "MGYG000000001"
        assert get_species_id("gm|GMGC10.001_001_001.PROT|x") is None
        assert get_species_id("uh_v2|MGYG000000001_00001|x") == "MGYG000000001"
        assert get_species_id("lake-3|sp|P02768|ALBU_HUMAN") == "HUMAN"
        assert get_species_id("sp|P02768|ALBU_HUMAN") == "HUMAN"

    def test_unresolvable_ids_are_none_not_a_pseudo_species(self):
        # Assembly proteins and smORFs need a taxonomy map; content-hash ids too.
        assert get_species_id("MuPr_Assembly_SAMP1_00001") is None
        assert get_species_id("SMORF_00001") is None
        assert get_species_id("PREDICT_MG_0123abcd") is None
        assert get_species_id("k141_12345_1") is None


# ── Helpers with file I/O ───────────────────────────────────────────────────


class TestLoadFastaProteins:
    """Test load_fasta_proteins with the tiny fixture."""

    def test_loads_correct_count(self, tiny_fasta: Path):
        proteins = load_fasta_proteins(tiny_fasta)
        assert len(proteins) == 20

    def test_protein_ids_correct(self, tiny_fasta: Path):
        proteins = load_fasta_proteins(tiny_fasta)
        assert "sp|P12345|GFP_HUMAN" in proteins
        assert "MGYG000000001_00001" in proteins
        assert "GMGC10.001_001_001.PROT01" in proteins
        assert "MuPr_Assembly_SAMP1_00001" in proteins
        assert "SMORF_00001" in proteins

    def test_sequences_not_empty(self, tiny_fasta: Path):
        proteins = load_fasta_proteins(tiny_fasta)
        for seq in proteins.values():
            assert len(seq) > 0


class TestCollectPeptides:
    """Test collect_peptides_with_scores with the tiny fixture."""

    def test_loads_correct_count(self, tiny_predictions: Path):
        peptides = collect_peptides_with_scores(tiny_predictions)
        assert len(peptides) == 15

    def test_scores_positive(self, tiny_predictions: Path):
        peptides = collect_peptides_with_scores(tiny_predictions)
        for score in peptides.values():
            assert score > 0.0

    def test_il_normalization_applied(self, tiny_predictions: Path):
        peptides = collect_peptides_with_scores(tiny_predictions)
        for pep in peptides:
            assert "I" not in pep  # All I should be replaced with L

    def test_length_filter(self, tiny_predictions: Path):
        peptides = collect_peptides_with_scores(tiny_predictions, min_length=16)
        # Only peptides >= 16 chars should pass
        assert all(len(p) >= 16 for p in peptides)
        assert len(peptides) < 15  # Some should be filtered


class TestMapPeptidesToProteins:
    """Test map_peptides_to_proteins with tiny fixtures."""

    def test_returns_two_dicts(self, tiny_fasta: Path, tiny_predictions: Path):
        proteins = load_fasta_proteins(tiny_fasta)
        peptides = collect_peptides_with_scores(tiny_predictions)
        pep2prot, prot2pep = map_peptides_to_proteins(proteins, peptides)

        assert isinstance(pep2prot, dict)
        assert isinstance(prot2pep, dict)
        assert len(pep2prot) > 0
        assert len(prot2pep) > 0

    def test_bidirectional_consistency(self, tiny_fasta: Path, tiny_predictions: Path):
        proteins = load_fasta_proteins(tiny_fasta)
        peptides = collect_peptides_with_scores(tiny_predictions)
        pep2prot, prot2pep = map_peptides_to_proteins(proteins, peptides)

        # Every protein in prot2pep should appear in pep2prot values
        for prot, peps in prot2pep.items():
            for pep in peps:
                assert prot in pep2prot[pep]


class TestLoadFastaRefusesDefectiveInput:
    """Duplicate accessions, empty records and bare '>' lines are errors, not silent losses."""

    def test_duplicate_accession_is_refused(self, tmp_path):
        p = tmp_path / "dup.fasta"
        p.write_text(">X first\nMKV\n>Y\nMAA\n>X second copy\nMLL\n")
        with pytest.raises(ValueError, match="duplicate accession 'X'"):
            load_fasta_proteins(p)

    def test_empty_sequence_is_refused(self, tmp_path):
        p = tmp_path / "empty.fasta"
        p.write_text(">X\nMKV\n>EMPTY\n>Y\nMAA\n")
        with pytest.raises(ValueError, match="record 'EMPTY'.*has no sequence"):
            load_fasta_proteins(p)

    def test_trailing_empty_record_is_refused(self, tmp_path):
        p = tmp_path / "trail.fasta"
        p.write_text(">X\nMKV\n>LAST\n")
        with pytest.raises(ValueError, match="record 'LAST'"):
            load_fasta_proteins(p)

    def test_bare_header_is_refused(self, tmp_path):
        p = tmp_path / "bare.fasta"
        p.write_text(">\nMKV\n")
        with pytest.raises(ValueError, match="empty FASTA header"):
            load_fasta_proteins(p)

    def test_clean_file_unchanged(self, tmp_path):
        p = tmp_path / "ok.fasta"
        p.write_text(">X desc\nMKV\nLLK\n>Y\nMAA\n")
        assert load_fasta_proteins(p) == {"X": "MKVLLK", "Y": "MAA"}
