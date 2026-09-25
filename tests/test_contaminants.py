"""Tests for fasta_lake.contaminants (v5.1)."""

from fasta_lake.contaminants import (
    DEFAULT_CONTAMINANTS,
    ContaminantDB,
    _extract_accession,
    exclude_contaminants,
    is_contaminant_id,
    tag_contaminants,
)


class TestExtractAccession:
    def test_swissprot_format(self):
        assert _extract_accession("sp|P00761|TRYP_BOVIN") == "P00761"

    def test_trembl_format(self):
        assert _extract_accession("tr|Q12345|NAME_HUMAN") == "Q12345"

    def test_raw_accession(self):
        assert _extract_accession("P00761") == "P00761"

    def test_con_prefix_stripped(self):
        assert _extract_accession("CON_P00761") == "P00761"

    def test_rev_prefix_stripped(self):
        assert _extract_accession("rev_P00761") == "P00761"

    def test_non_uniprot_returns_none(self):
        assert _extract_accession("MGYG00000001") is None  # bacterial ID
        assert _extract_accession("ENSP00000123456") is None  # Ensembl


class TestContaminantDB:
    def test_default_contaminants_present(self):
        assert len(DEFAULT_CONTAMINANTS) > 10
        assert "P00761" in DEFAULT_CONTAMINANTS  # trypsin
        assert "P02769" in DEFAULT_CONTAMINANTS  # BSA

    def test_db_init_uses_defaults(self):
        db = ContaminantDB()
        assert "P00761" in db.accessions

    def test_db_init_custom(self):
        db = ContaminantDB(accessions={"X12345"})
        assert db.accessions == {"X12345"}
        assert "P00761" not in db.accessions  # defaults replaced

    def test_is_contaminant_direct(self):
        db = ContaminantDB()
        assert db.is_contaminant("P00761") is True
        assert db.is_contaminant("sp|P00761|TRYP") is True

    def test_is_contaminant_already_tagged(self):
        db = ContaminantDB()
        assert db.is_contaminant("CON_P00761") is True
        assert db.is_contaminant("CON_SOMETHING_NEW") is True  # prefix alone suffices

    def test_is_contaminant_false(self):
        db = ContaminantDB()
        assert db.is_contaminant("P04637") is False  # p53, not a contaminant
        assert db.is_contaminant("MGYG00000001") is False  # bacterial

    def test_load_from_text_file(self, tmp_path):
        p = tmp_path / "custom.txt"
        p.write_text("# comment\nX99999\n\nZ88888\n")
        db = ContaminantDB()
        db.load(p)
        assert "X99999" in db.accessions
        assert "Z88888" in db.accessions

    def test_load_from_fasta(self, tmp_path):
        p = tmp_path / "custom.fasta"
        p.write_text(">sp|X99999|TEST protein\nMAKS\n>tr|Z88888|TEST2\nKRMA\n")
        db = ContaminantDB()
        db.load_from_fasta(p)
        assert "X99999" in db.accessions
        assert "Z88888" in db.accessions


class TestTagContaminants:
    def test_tag_matching(self):
        proteins = {
            "sp|P00761|TRYP_BOVIN": "IVGGYTCGANTVP",  # trypsin
            "sp|P04637|P53_HUMAN": "MEEPQSDPSVEPP",  # not a contaminant
        }
        tagged, contam_ids = tag_contaminants(proteins)
        assert "CON_sp|P00761|TRYP_BOVIN" in tagged
        assert "sp|P04637|P53_HUMAN" in tagged
        assert len(contam_ids) == 1

    def test_already_tagged_not_doubled(self):
        proteins = {"CON_sp|P00761|X": "MAKS"}
        tagged, _ = tag_contaminants(proteins)
        assert "CON_sp|P00761|X" in tagged
        assert "CON_CON_sp|P00761|X" not in tagged

    def test_custom_db(self):
        db = ContaminantDB(accessions={"P99999"})
        proteins = {"sp|P99999|CUSTOM": "MAKS", "sp|P00761|TRYP": "MAKS"}
        tagged, _ = tag_contaminants(proteins, db)
        assert "CON_sp|P99999|CUSTOM" in tagged
        # P00761 not in custom db → not tagged
        assert "sp|P00761|TRYP" in tagged


class TestExcludeContaminants:
    def test_exclude_filters_psms(self):
        rows = [
            {"proteins": "sp|P04637|P53_HUMAN", "value": 1},
            {"proteins": "CON_sp|P00761|TRYP_BOVIN", "value": 2},
            {"proteins": "sp|P00761|TRYP_BOVIN", "value": 3},  # untagged but in default list
            {"proteins": "MGYG0001", "value": 4},
        ]
        out = exclude_contaminants(rows)
        values = [r["value"] for r in out]
        assert 1 in values
        assert 2 not in values  # CON_ filtered
        assert 3 not in values  # P00761 matched as default contaminant
        assert 4 in values  # bacterial, kept

    def test_multi_protein_string_uses_first(self):
        rows = [
            {"proteins": "sp|P04637|P53;sp|P00761|TRYP", "value": 1},  # first = p53 → kept
            {"proteins": "sp|P00761|TRYP;sp|P04637|P53", "value": 2},  # first = trypsin → filtered
        ]
        out = exclude_contaminants(rows)
        values = [r["value"] for r in out]
        assert 1 in values
        assert 2 not in values


def test_is_contaminant_id_helper():
    assert is_contaminant_id("P00761") is True
    assert is_contaminant_id("P04637") is False
