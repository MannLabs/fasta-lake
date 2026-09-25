"""Integration tests for fasta_lake.lake_builder (v5.2)."""

from pathlib import Path

from fasta_lake.lake_builder import (
    build_lake,
    iter_fasta,
    parse_header,
)
from fasta_lake.lake_merge import read_provenance_sidecar


def _write(p: Path, lines: list[str]) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


class TestParseHeader:
    def test_swissprot(self):
        r = parse_header(
            "sp|P04637|P53_HUMAN Cellular tumor antigen p53 "
            "OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=4"
        )
        assert r.tag == "SP"
        assert r.accession == "P04637"
        assert r.gene == "TP53"
        assert r.description == "Cellular tumor antigen p53"
        assert r.organism == "Homo sapiens"
        assert r.tax_id == 9606
        assert r.pe == 1
        assert r.sv == 4

    def test_trembl(self):
        r = parse_header(
            "tr|H8LKR3|H8LKR3_HUMAN Putative uncharacterized protein OS=Homo sapiens OX=9606 PE=2 "
            "SV=1"
        )
        assert r.tag == "TR"
        assert r.accession == "H8LKR3"
        assert r.pe == 2

    def test_uhgp(self):
        r = parse_header("MGYP000123456 hypothetical protein")
        assert r.tag == "UH"
        assert r.accession == "MGYP000123456"
        assert r.description == "hypothetical protein"

    def test_gencode(self):
        r = parse_header(
            "ENSP00000269305.4|ENST00000269305.9|ENSG00000141510.19|TP53|cellular tumor antigen p53"
        )
        assert r.tag == "GENCODE"
        assert r.accession == "ENSP00000269305.4"
        assert r.gene == "TP53"

    def test_generic_fallback(self):
        r = parse_header("PROT_xyz123 some random header", tag_hint="NCBI")
        assert r.tag == "NCBI"
        assert r.accession == "PROT_xyz123"


class TestIterFasta:
    def test_reads_multiple_records(self, tmp_path: Path):
        p = tmp_path / "sp.fasta"
        _write(
            p,
            [
                ">sp|P04637|P53_HUMAN p53 OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=4",
                "MEEPQSDPSV",
                "EPPLSQETFSDLWKLLPEN",
                ">sp|P00533|EGFR_HUMAN EGFR OS=Homo sapiens OX=9606 GN=EGFR PE=1 SV=2",
                "MRPSGTAGAAL",
            ],
        )
        recs = iter_fasta(p, tag="SP")
        assert len(recs) == 2
        r0, s0 = recs[0]
        assert r0.accession == "P04637"
        assert s0.startswith("MEEPQSDPSV")
        assert "EPPLSQETFSDLWKLLPEN" in s0  # continues on next line


class TestBuildLake:
    def _setup_mini_lake(self, tmp: Path) -> dict[str, Path]:
        sp = tmp / "sp.fasta"
        _write(
            sp,
            [
                ">sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens OX=9606 GN=TP53 "
                "PE=1 SV=4",
                "MEEPQSDPSVEPPLSQETFSDLWKLLPEN",
                ">sp|P00533|EGFR_HUMAN EGFR OS=Homo sapiens OX=9606 GN=EGFR PE=1 SV=2",
                "MRPSGTAGAALLALLAALCPASRA",
            ],
        )
        # TrEMBL: shares P53 sequence with SP (same AA, different accession)
        tr = tmp / "tr.fasta"
        _write(
            tr,
            [
                ">tr|H8LKR3|H8LKR3_HUMAN Putative uncharacterized protein OS=Homo sapiens OX=9606 "
                "GN=TP53 PE=2 SV=1",
                "MEEPQSDPSVEPPLSQETFSDLWKLLPEN",
                ">tr|H8LKR4|H8LKR4_HUMAN Another OS=Homo sapiens OX=9606 PE=4 SV=1",
                "AAACCCGGGTTT",
            ],
        )
        # UHGP: shares P53 sequence too (a cross-kingdom conflict is contrived but
        # mimics real ribosomal-style identical-across-kingdoms cases)
        uh = tmp / "uh.fasta"
        _write(
            uh,
            [
                ">MGYP000000123 hypothetical protein",
                "MEEPQSDPSVEPPLSQETFSDLWKLLPEN",
                ">MGYP000000124 ribosomal protein",
                "MKKLGTPLVM",
            ],
        )
        return {"SP": sp, "TR": tr, "UH": uh}

    def test_merge_mode(self, tmp_path: Path):
        srcs = self._setup_mini_lake(tmp_path / "src")
        out_fa = tmp_path / "out" / "lake.fasta"
        out_pv = tmp_path / "out" / "lake.provenance.jsonl"
        stats = build_lake(srcs, out_fa, out_pv, dedup_mode="merge")

        # 2 + 2 + 2 = 6 inputs; P53 dedup across SP/TR/UH → 4 unique
        assert stats["n_input_sequences"] == 6
        assert stats["n_unique_sequences"] == 4
        assert stats["n_shared_across_sources"] == 1  # just P53
        # UHGP says "hypothetical protein" (no organism) so only 2 orgs disagree
        # when SP(Homo sapiens) and UHGP(no org) combine → no conflict
        # Actually UHGP has organism=None, so only SP+TR have organism
        # and they agree → 0 ambiguous
        assert stats["n_ambiguous_organism"] == 0

        # Sidecar round-trips
        recs = read_provenance_sidecar(out_pv)
        assert len(recs) == 4
        shared = [r for r in recs if len(r.sources) > 1]
        assert len(shared) == 1
        p53 = shared[0]
        tags = sorted({s.tag for s in p53.sources})
        assert tags == ["SP", "TR", "UH"]
        assert p53.primary_tag == "SP"  # SP beats TR beats UH
        assert p53.gene == "TP53"
        assert p53.pe == 1  # min (best)
        assert p53.sv == 4  # max

    def test_priority_win_mode_discards_other_sources(self, tmp_path: Path):
        srcs = self._setup_mini_lake(tmp_path / "src")
        out_fa = tmp_path / "out" / "lake.fasta"
        out_pv = tmp_path / "out" / "lake.provenance.jsonl"
        stats = build_lake(srcs, out_fa, out_pv, dedup_mode="priority-win")
        assert stats["dedup_mode"] == "priority-win"

        recs = read_provenance_sidecar(out_pv)
        p53 = next(r for r in recs if r.primary_accession == "P04637")
        # priority-win: ONLY SP retained, TR and UH discarded
        assert [s.tag for s in p53.sources] == ["SP"]

    def test_fasta_output_valid(self, tmp_path: Path):
        srcs = self._setup_mini_lake(tmp_path / "src")
        out_fa = tmp_path / "out.fasta"
        out_pv = tmp_path / "out.jsonl"
        build_lake(srcs, out_fa, out_pv)

        text = out_fa.read_text()
        assert text.count(">") == 4
        # Each sequence wrapped at 60 cols — verify no header collisions
        lines = text.splitlines()
        headers = [ln for ln in lines if ln.startswith(">")]
        accs = [h.split("|")[1] if "|" in h else h for h in headers]
        assert len(set(accs)) == 4

    def test_custom_priority_flips_primary(self, tmp_path: Path):
        srcs = self._setup_mini_lake(tmp_path / "src")
        out_fa = tmp_path / "out.fasta"
        out_pv = tmp_path / "out.jsonl"
        build_lake(srcs, out_fa, out_pv, priority=("UH", "SP", "TR"))
        recs = read_provenance_sidecar(out_pv)
        p53 = next(r for r in recs if "MEEPQSD" in str(r.seq_hash) or len(r.sources) > 1)
        assert p53.primary_tag == "UH"

    def test_per_source_contribution_stats(self, tmp_path: Path):
        srcs = self._setup_mini_lake(tmp_path / "src")
        stats = build_lake(srcs, tmp_path / "out.fa", tmp_path / "out.jsonl")
        # All SP records retained in some merged record, etc.
        assert stats["per_source_input"] == {"SP": 2, "TR": 2, "UH": 2}
        # SP contributed to 2 merged records (P53 shared + EGFR unique)
        assert stats["per_source_contributed"]["SP"] == 2
        assert stats["per_source_contributed"]["TR"] == 2
        assert stats["per_source_contributed"]["UH"] == 2

    def test_compress_sidecar(self, tmp_path: Path):
        srcs = self._setup_mini_lake(tmp_path / "src")
        out_pv = tmp_path / "out.jsonl.gz"
        build_lake(srcs, tmp_path / "out.fa", out_pv, compress=True)
        assert out_pv.exists()
        assert out_pv.read_bytes()[:2] == b"\x1f\x8b"  # gzip magic
        recs = read_provenance_sidecar(out_pv)
        assert len(recs) == 4

    def test_embed_sources_header(self, tmp_path: Path):
        srcs = self._setup_mini_lake(tmp_path / "src")
        out_fa = tmp_path / "out.fa"
        build_lake(srcs, out_fa, tmp_path / "out.jsonl", embed_sources=True)
        text = out_fa.read_text()
        assert "[FLlake:sources=" in text


class TestUniprotIsoformAccessions:
    """SwissProt varsplic isoforms (P04637-2) are UniProt entries, not generic headers."""

    def test_isoform_header_parses_as_uniprot(self):
        r = parse_header(
            "sp|P04637-2|P53_HUMAN Isoform 2 of Cellular tumor antigen p53 "
            "OS=Homo sapiens OX=9606 GN=TP53"
        )
        assert r.tag == "SP"
        assert r.accession == "P04637-2"
        assert r.gene == "TP53"
        assert r.organism == "Homo sapiens"
        assert r.tax_id == 9606
        assert r.extra["entry_name"] == "P53_HUMAN"
        assert r.description == "Isoform 2 of Cellular tumor antigen p53"

    def test_trembl_isoform_too(self):
        r = parse_header("tr|A0A024R161-3|A0A024R161_HUMAN Some protein OS=Homo sapiens OX=9606")
        assert r.tag == "TR" and r.accession == "A0A024R161-3"
