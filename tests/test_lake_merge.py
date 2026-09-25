"""Unit tests for fasta_lake.lake_merge (v5.2 annotation-aware dedup)."""

import json
from pathlib import Path

import pytest

from fasta_lake.lake_merge import (
    SourceRecord,
    format_fasta_header,
    merge_records,
    priority_win,
    read_provenance_sidecar,
    sequence_hash,
    write_provenance_sidecar,
)

SEQ = "MEEPQSDPSVEPPLSQETFSDLWKLLPEN"


def _sp() -> SourceRecord:
    return SourceRecord(
        tag="SP",
        accession="P04637",
        gene="TP53",
        description="Cellular tumor antigen p53",
        organism="Homo sapiens",
        tax_id=9606,
        pe=1,
        sv=4,
    )


def _tr() -> SourceRecord:
    return SourceRecord(
        tag="TR",
        accession="H8LKR3",
        gene="TP53",
        description="Cellular tumor antigen p53",
        organism="Homo sapiens",
        tax_id=9606,
        pe=2,
        sv=1,
    )


def _uh() -> SourceRecord:
    return SourceRecord(
        tag="UH",
        accession="UHGP_000123",
        gene=None,
        description="hypothetical protein",
        organism="Unknown bacterium",
        tax_id=2,
        pe=None,
        sv=None,
    )


class TestMergeRecords:
    def test_empty_raises(self):
        with pytest.raises(ValueError):
            merge_records(SEQ, [])

    def test_single_source_passthrough(self):
        r = merge_records(SEQ, [_sp()])
        assert r.primary_tag == "SP"
        assert r.primary_accession == "P04637"
        assert r.gene == "TP53"
        assert r.organism == "Homo sapiens"
        assert r.tax_ids == [9606]
        assert r.pe == 1
        assert r.sv == 4
        assert not r.conflicts

    def test_primary_follows_priority(self):
        """SP beats TR beats UH."""
        r = merge_records(SEQ, [_uh(), _tr(), _sp()])
        assert r.primary_tag == "SP"
        assert r.primary_accession == "P04637"

    def test_custom_priority_flips_primary(self):
        """--priority UH,SP,TR → UH becomes primary."""
        r = merge_records(
            SEQ,
            [_sp(), _tr(), _uh()],
            priority=("UH", "SP", "TR"),
        )
        assert r.primary_tag == "UH"

    def test_gene_merges_from_highest_priority_with_gene(self):
        """UH has null gene, SP has TP53 → TP53 wins."""
        r = merge_records(SEQ, [_uh(), _sp()])
        assert r.gene == "TP53"

    def test_description_prefers_informative(self):
        """'hypothetical protein' loses to anything informative."""
        r = merge_records(SEQ, [_uh(), _sp()])
        assert r.description == "Cellular tumor antigen p53"

    def test_description_falls_back_to_longest_when_all_uninformative(self):
        a = SourceRecord(tag="SP", accession="X", description="hypothetical protein X")
        b = SourceRecord(
            tag="TR", accession="Y", description="uncharacterized protein with longer name"
        )
        r = merge_records(SEQ, [a, b])
        assert r.description == "uncharacterized protein with longer name"

    def test_organism_agreement(self):
        r = merge_records(SEQ, [_sp(), _tr()])  # both Homo sapiens
        assert r.organism == "Homo sapiens"
        assert "organism" not in r.conflicts

    def test_organism_conflict_marks_ambiguous(self):
        r = merge_records(SEQ, [_sp(), _uh()])  # Homo sapiens vs Unknown bacterium
        assert r.organism == "ambiguous"
        assert sorted(r.conflicts["organism"]) == ["Homo sapiens", "Unknown bacterium"]

    def test_tax_ids_union(self):
        r = merge_records(SEQ, [_sp(), _uh()])
        assert r.tax_ids == [2, 9606]

    def test_pe_takes_min_best_evidence(self):
        """PE 1 (experimental) beats PE 2 (transcript-level)."""
        r = merge_records(SEQ, [_sp(), _tr()])  # pe=1 vs pe=2
        assert r.pe == 1

    def test_sv_takes_max(self):
        r = merge_records(SEQ, [_sp(), _tr()])  # sv=4 vs sv=1
        assert r.sv == 4

    def test_contaminant_or_semantics(self):
        a = SourceRecord(tag="SP", accession="X", is_contaminant=True)
        b = SourceRecord(tag="TR", accession="Y", is_contaminant=False)
        r = merge_records(SEQ, [a, b])
        assert r.is_contaminant is True

    def test_gene_conflict_flagged(self):
        a = SourceRecord(tag="SP", accession="X", gene="TP53")
        b = SourceRecord(tag="TR", accession="Y", gene="TP63")
        r = merge_records(SEQ, [a, b])
        assert r.gene == "TP53"  # SP wins
        assert r.conflicts["gene"] == ["TP53", "TP63"]

    def test_length_and_hash(self):
        r = merge_records(SEQ, [_sp()])
        assert r.length == len(SEQ)
        assert r.seq_hash == sequence_hash(SEQ)
        assert len(r.seq_hash) == 64  # SHA256 hex digest (v5.3; was 32 for MD5)

    def test_sequence_hash_case_insensitive(self):
        assert sequence_hash(SEQ) == sequence_hash(SEQ.lower())


class TestPriorityWin:
    def test_picks_highest_priority(self):
        w = priority_win([_uh(), _tr(), _sp()])
        assert w.tag == "SP"
        assert w.accession == "P04637"

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            priority_win([])


class TestFastaHeader:
    def test_basic_header(self):
        r = merge_records(SEQ, [_sp()])
        h = format_fasta_header(r)
        assert h.startswith("sp|P04637")
        assert "GN=TP53" in h
        assert "OS=Homo sapiens" in h
        assert "PE=1" in h
        assert "[FLlake" not in h

    def test_embed_sources_opt_in(self):
        r = merge_records(SEQ, [_sp(), _uh()])
        h = format_fasta_header(r, embed_sources=True)
        assert "[FLlake:sources=" in h
        # sources sorted alphabetically
        assert "SP" in h and "UH" in h

    def test_entry_name_preserved_when_present(self):
        """If the primary SP/TR source has an original entry name
        (e.g. "P53_HUMAN"), the downstream header emits three pipe fields:
        sp|{acc}|{entry_name}. No synthesised accession_gene form."""
        sp = SourceRecord(
            tag="SP",
            accession="P04637",
            gene="TP53",
            description="Cellular tumor antigen p53",
            organism="Homo sapiens",
            tax_id=9606,
            pe=1,
            sv=4,
            extra={"entry_name": "P53_HUMAN"},
        )
        r = merge_records(SEQ, [sp])
        h = format_fasta_header(r)
        assert h.startswith("sp|P04637|P53_HUMAN ")
        # The synthesised form (the old bug) must not appear.
        assert "P04637_TP53" not in h

    def test_no_entry_name_emits_two_fields_not_synthesised(self):
        """Without an entry name (e.g. non-UniProt sources), the header
        must NOT synthesise a third pipe field like accession_gene. Two
        fields is the correct degenerate form — downstream tools key off
        field 2 (accession)."""
        r = merge_records(SEQ, [_uh()])  # UHGP, no entry_name in extra
        h = format_fasta_header(r)
        # Exactly 2 pipes in the leading token (before whitespace).
        head = h.split()[0]
        assert head.count("|") == 1, f"expected 2-field head, got: {head}"
        assert head == "uh|UHGP_000123"

    def test_organism_ambiguous_falls_back_to_primary(self):
        """When sources disagree on organism, the MergedRecord has
        organism='ambiguous' (for the sidecar's benefit). But the
        downstream FASTA header must show the PRIMARY source's organism,
        not the literal word 'ambiguous' — search engines can't parse it."""
        sp_bad_neighbour = SourceRecord(
            tag="SP",
            accession="P04637",
            gene="TP53",
            description="Cellular tumor antigen p53",
            organism="Homo sapiens",
            tax_id=9606,
            pe=1,
            sv=4,
        )
        uh_wrong_org = SourceRecord(
            tag="UH",
            accession="UHGP_BOGUS",
            gene=None,
            description="hypothetical",
            organism="Unknown bacterium",
            tax_id=2,
        )
        r = merge_records(SEQ, [sp_bad_neighbour, uh_wrong_org])
        assert r.organism == "ambiguous"  # sidecar view
        h = format_fasta_header(r)
        assert "OS=Homo sapiens" in h  # header view (primary wins)
        assert "OS=ambiguous" not in h

    def test_tax_id_follows_primary_not_sorted_min(self):
        """When sources contribute different tax_ids, the header's OX=
        must be the PRIMARY source's tax_id, not min() of all of them.
        (Pre-fix bug: OX=2 (Bacteria) overrode OX=9606 (Human)."""
        sp = SourceRecord(
            tag="SP",
            accession="P04637",
            gene="TP53",
            description="Cellular tumor antigen p53",
            organism="Homo sapiens",
            tax_id=9606,
            pe=1,
            sv=4,
        )
        uh = SourceRecord(
            tag="UH",
            accession="UHGP_X",
            gene=None,
            description="hypothetical",
            organism="Unknown bacterium",
            tax_id=2,
        )
        r = merge_records(SEQ, [sp, uh])
        assert r.tax_ids == [2, 9606]  # sidecar: full set
        h = format_fasta_header(r)
        assert "OX=9606" in h  # header: primary only
        assert "OX=2 " not in h and not h.endswith("OX=2")


class TestDerivedSidecarFields:
    """The sidecar JSONL embeds two convenience fields (n_sources,
    alt_accessions) on top of the full sources list, so that the most
    common user query — 'show me every sequence where >1 source
    contributed the same bytes' — is a trivial grep/jq."""

    def test_n_sources_single_source(self, tmp_path):
        rec = merge_records(SEQ, [_sp()])
        p = tmp_path / "p.jsonl"
        write_provenance_sidecar(p, [rec])
        import json as _json

        line = p.read_text().strip()
        d = _json.loads(line)
        assert d["n_sources"] == 1
        assert d["alt_accessions"] == []

    def test_n_sources_multi_source_collapse(self, tmp_path):
        """Three sources contribute byte-identical sequence. The primary
        accession stays in `primary_accession`; the other two appear in
        `alt_accessions` — this is the "also known as" the user will grep."""
        rec = merge_records(SEQ, [_sp(), _tr(), _uh()])
        p = tmp_path / "p.jsonl"
        write_provenance_sidecar(p, [rec])
        import json as _json

        d = _json.loads(p.read_text().strip())
        assert d["n_sources"] == 3
        assert d["primary_accession"] == "P04637"
        assert set(d["alt_accessions"]) == {"H8LKR3", "UHGP_000123"}

    def test_derived_fields_stripped_on_roundtrip(self, tmp_path):
        """The derived fields are purely for downstream users — they must
        not break read_provenance_sidecar (which otherwise would get an
        unknown-kwarg error when constructing MergedRecord)."""
        rec = merge_records(SEQ, [_sp(), _tr()])
        p = tmp_path / "p.jsonl"
        write_provenance_sidecar(p, [rec])
        loaded = read_provenance_sidecar(p)
        assert len(loaded) == 1
        assert loaded[0].primary_accession == rec.primary_accession
        assert len(loaded[0].sources) == 2


class TestSidecarRoundtrip:
    def test_write_and_read(self, tmp_path: Path):
        recs = [
            merge_records(SEQ, [_sp(), _tr()]),
            merge_records(SEQ + "X", [_sp(), _uh()]),
        ]
        p = tmp_path / "lake.provenance.jsonl"
        n = write_provenance_sidecar(p, recs)
        assert n == 2
        loaded = read_provenance_sidecar(p)
        assert len(loaded) == 2
        assert loaded[0].primary_accession == recs[0].primary_accession
        assert loaded[0].tax_ids == recs[0].tax_ids
        assert len(loaded[1].sources) == 2

    def test_roundtrip_preserves_conflicts(self, tmp_path: Path):
        rec = merge_records(SEQ, [_sp(), _uh()])  # organism conflict
        p = tmp_path / "p.jsonl"
        write_provenance_sidecar(p, [rec])
        back = read_provenance_sidecar(p)[0]
        assert back.organism == "ambiguous"
        assert sorted(back.conflicts["organism"]) == ["Homo sapiens", "Unknown bacterium"]

    def test_gzip_roundtrip(self, tmp_path: Path):
        rec = merge_records(SEQ, [_sp()])
        p = tmp_path / "p.jsonl.gz"
        write_provenance_sidecar(p, [rec], compress=True)
        assert p.exists()
        # first two bytes of gzip magic
        assert p.read_bytes()[:2] == b"\x1f\x8b"
        back = read_provenance_sidecar(p)
        assert back[0].primary_accession == "P04637"

    def test_compact_json_one_line_per_record(self, tmp_path: Path):
        recs = [merge_records(SEQ, [_sp()]), merge_records(SEQ + "X", [_uh()])]
        p = tmp_path / "p.jsonl"
        write_provenance_sidecar(p, recs)
        lines = p.read_text().strip().split("\n")
        assert len(lines) == 2
        for ln in lines:
            assert json.loads(ln)  # valid JSON


class TestIsoformHeaderRoundTrip:
    def test_isoform_record_keeps_one_prefix_and_structured_tags(self):
        from fasta_lake.lake_builder import parse_header

        src = parse_header(
            "sp|P04637-2|P53_HUMAN Isoform 2 of Cellular tumor antigen p53 "
            "OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=2"
        )
        rec = merge_records("MEEPQSDPSVEPPLSQETFSDLWK", [src])
        header = format_fasta_header(rec)
        assert header.startswith("sp|P04637-2|P53_HUMAN ")
        assert "sp|sp|" not in header
        assert "OS=Homo sapiens" in header and "OX=9606" in header and "GN=TP53" in header
