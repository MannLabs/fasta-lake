"""Tests for the v4 header parsing and formatting system."""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from fasta_lake.headers import heal_fasta, parse_header

# ═══════════════════════════════════════════════════════════════
# PARSING TESTS
# ═══════════════════════════════════════════════════════════════


def test_parse_swissprot():
    p = parse_header("sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2")
    assert p.source_db == "swissprot"
    assert p.gene == "ALB"
    assert p.organism == "Homo sapiens"
    assert p.accession == "sp|P02768|ALBU_HUMAN"


def test_parse_swissprot_no_gn():
    p = parse_header(
        "sp|P0DOX5|IGG1_HUMAN Immunoglobulin gamma-1 heavy chain OS=Homo sapiens OX=9606 PE=1 SV=2"
    )
    assert p.source_db == "swissprot"
    # No GN= means no gene: the entry-name mnemonic (IGG1_HUMAN) is not a symbol.
    assert p.gene is None


def test_parse_trembl_with_gn():
    p = parse_header(
        (
            "tr|A0A218KGR2|A0A218KGR2_HUMAN Amyloid-beta A4 "
            "protein OS=Homo sapiens OX=9606 GN=APP PE=3 SV=1"
        )
    )
    assert p.source_db == "trembl"
    assert p.gene == "APP"


def test_parse_trembl_no_gn():
    p = parse_header(
        "tr|B2R6Q2|B2R6Q2_HUMAN Alkaline phosphatase OS=Homo sapiens OX=9606 PE=2 SV=1"
    )
    assert p.source_db == "trembl"
    # TrEMBL entry names are ACCESSION_SPECIES; using them made every accession a "gene".
    assert p.gene is None


def test_uniprot_without_gene_is_written_without_gn():
    """No resolved gene: no GN= is invented for sp|/tr| in any engine format."""
    raw = "tr|B2R6Q2|B2R6Q2_HUMAN Alkaline phosphatase OS=Homo sapiens OX=9606 PE=2 SV=1"
    p = parse_header(raw)
    assert p.format_for_diann() == raw
    assert p.format_for_sage() == raw
    assert p.format_for_msfragger() == raw
    p.gene = "ALPL"  # a resolved gene is added, before PE=
    assert p.format_for_diann().endswith("OX=9606 GN=ALPL PE=2 SV=1")
    assert "GN=B2R6Q2" not in p.format_for_diann()


def test_heal_resolves_trembl_gene_from_canonical_description(tmp_path):
    """With no fabricated gene, heal_fasta consults the canonical map (dead before)."""
    canonical = tmp_path / "canonical.fasta"
    canonical.write_text(
        ">sp|P05186|PPBT_HUMAN Alkaline phosphatase, tissue-nonspecific isozyme "
        "OS=Homo sapiens OX=9606 GN=ALPL PE=1 SV=4\nMISPFLVLAIGTCLTNSLVPEK\n"
    )
    inp = tmp_path / "in.fasta"
    inp.write_text(
        ">tr|B2R6Q2|B2R6Q2_HUMAN Alkaline phosphatase, tissue-nonspecific isozyme "
        "OS=Homo sapiens OX=9606 PE=2 SV=1\nMISPFLVLAIGTCLTNSLVPEK\n"
        ">tr|A0A024R161|A0A024R161_HUMAN Uncharacterized protein "
        "OS=Homo sapiens OX=9606 PE=4 SV=1\nMKWVTFISLL\n"
    )
    out = tmp_path / "out.fasta"
    heal_fasta(inp, out, target_engine="diann", canonical_fasta=canonical)
    headers = [ln for ln in out.read_text().splitlines() if ln.startswith(">")]
    assert headers[0].endswith("OX=9606 GN=ALPL PE=2 SV=1")
    assert "GN=" not in headers[1], headers[1]  # unresolved: honest, not fabricated


def test_parse_gnomad():
    p = parse_header(
        "var-nfe|P02768|ALB_Glu177Asp|4-73409403-A-T sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens "
        "OX=9606 GN=ALB PE=1 SV=4 [VARIANT: p.Glu177Asp] [AF_NFE: 0.000897]"
    )
    assert p.source_db == "gnomad"
    assert p.gene == "ALB"
    assert p.population == "nfe"
    assert p.variant_info == "p.Glu177Asp"
    assert "var-nfe|P02768|ALB_Glu177Asp|4-73409403-A-T" in p.accession


def test_parse_gencode():
    p = parse_header(
        "ENSP00000493376.2|ENST00000641515.2|ENSG00000186092.7|OTTHUMG00000001094.4|OTTHUMT00000003223.4|OR4F5-201|OR4F5|326"
    )
    assert p.source_db == "gencode"
    assert p.gene == "OR4F5"


def test_parse_ens_variant():
    p = parse_header(
        "ens|ENST00000393469.8_Arg664Gln|SLC12A8 SLC12A8 p.Arg664Gln (missense_variant)"
    )
    assert p.source_db == "personal_variant"
    assert p.gene == "SLC12A8"
    assert p.variant_info == "p.Arg664Gln"


def test_parse_mgyg():
    p = parse_header("MGYG000071032_00033 Adaptive-response sensory-kinase SasA")
    assert p.source_db == "uhgp"
    assert p.gene == "MGYG000071032_00033"
    assert p.description == "Adaptive-response sensory-kinase SasA"


def test_parse_gmgc():
    p = parse_header("GMGC10.286_822_098.UNKNOWN")
    assert p.source_db == "gmgc"
    assert p.gene == "GMGC10.286_822_098.UNKNOWN"


def test_parse_smorf():
    p = parse_header("SMORF_00001 Small ORF protein 1")
    assert p.source_db == "smorf"
    assert p.gene == "SMORF_00001"


def test_parse_assembly():
    p = parse_header("MuPr_Assembly_SAMP1_00001 Novel assembly protein")
    assert p.source_db == "assembly"


def test_parse_ncbi_with_gene():
    p = parse_header(
        (
            "CIFR|WP_003023861.1|JMV71_RS00005 gene=dnaA chromosomal "
            "replication initiator protein DnaA OS=Citrobacter freundii"
        )
    )
    assert p.source_db == "ncbi_microbe"
    assert p.gene == "dnaA"
    assert p.organism == "Citrobacter freundii"
    assert p.locus_tag == "JMV71_RS00005"


def test_parse_ncbi_no_gene():
    p = parse_header("KLPN|YP_005224301.1|KPHS_00010 flavodoxin OS=Klebsiella pneumoniae")
    assert p.source_db == "ncbi_microbe"
    assert p.gene == "KPHS_00010"  # falls back to locus tag


def test_parse_genscan():
    p = parse_header("GENSCAN00000033973 pep scaffold:GRCh38:1:100:200:1")
    assert p.source_db == "genscan"
    assert p.gene is None


def test_parse_isoform_pipe():
    p = parse_header("Q8IZF2|ADGRF5|ISOFORM_2|liver_hepatocyte")
    assert p.source_db == "isoform_db"
    assert p.gene == "ADGRF5"


def test_parse_ensembl_pep():
    p = parse_header(
        "ENSP00000481738.1 pep chromosome:GRCh38:21:10649400:10649835:-1 gene:ENSG00000277282.1 "
        "transcript:ENST00000622028.1 gene_biotype:IG_V_gene transcript_biotype:IG_V_gene "
        "gene_symbol:IGHV1OR21-1 description:immunoglobulin heavy variable 1/OR21-1"
    )
    assert p.source_db == "ensembl_pep"
    assert p.gene == "IGHV1OR21-1"
    assert p.description is not None


def test_parse_prodigal():
    p = parse_header(
        "1_1 # 1 # 1404 # 1 # "
        "ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.541"
    )
    assert p.source_db == "prodigal"
    assert p.gene == "1_1"


def test_parse_sample_ref():
    p = parse_header(
        "SAMPLE|REF|MGYG000012471_00850 Inosine-5'-monophosphate dehydrogenase [peptides=1]"
    )
    assert p.gene is not None
    assert "MGYG" in p.gene or "MGYG" in p.accession


def test_parse_lake_hybrid():
    p = parse_header(
        "ABR_Project_1 [SWISSPROT] sp|A0A087X1C5|CP2D7_HUMAN Cytochrome P450 2D7 OS=Homo sapiens "
        "OX=9606 GN=CYP2D7 PE=1 SV=1 md5=abc123 n_sources=1"
    )
    assert p.gene == "CYP2D7"
    assert p.source_db == "lake_hybrid"


def test_parse_md5_hash():
    p = parse_header("001fa9937b0a9b73013df2c361fc65bd")
    assert p.source_db == "md5_hash"
    assert p.gene == "001fa9937b0a9b73013df2c361fc65bd"


def test_parse_disease_variant():
    p = parse_header("sp|P02649_v21K|APOE_VAR E21K in ApoE5; associated with hyperlipoproteinemia")
    assert p.source_db == "disease_variant"
    assert p.gene == "APOE"
    assert p.variant_info == "21K"


# ═══════════════════════════════════════════════════════════════
# DIA-NN FORMATTING TESTS
# ═══════════════════════════════════════════════════════════════


def test_diann_swissprot_unchanged():
    """SwissProt with GN= should pass through unchanged."""
    h = "sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2"
    p = parse_header(h)
    assert p.format_for_diann() == h


def test_diann_gnomad_gene_first_word():
    """gnomAD variant must have gene as first word of description."""
    h = (
        "var-nfe|P02768|ALB_Glu177Asp|4-73409403-A-T sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens "
        "OX=9606 GN=ALB PE=1 SV=4 [VARIANT: p.Glu177Asp]"
    )
    p = parse_header(h)
    diann = p.format_for_diann()
    # Accession must be preserved
    assert diann.startswith("var-nfe|P02768|ALB_Glu177Asp|4-73409403-A-T ")
    # First word of description must be ALB
    desc = diann.split(None, 1)[1]
    assert desc.startswith("ALB ")
    # GN= must be present
    assert "GN=ALB" in diann
    # sp|P02768|ALBU_HUMAN must NOT be in description
    assert "sp|P02768|ALBU_HUMAN" not in desc


def test_diann_gencode_gene_first_word():
    """GENCODE must have gene as first word of description."""
    h = (
        "ENSP00000493376.2|ENST00000641515.2|ENSG00000186092.7|"
        "OTTHUMG00000001094.4|OTTHUMT00000003223.4|OR4F5-201|OR4F5|326"
    )
    p = parse_header(h)
    diann = p.format_for_diann()
    desc = diann.split(None, 1)[1]
    assert desc.startswith("OR4F5 ")
    assert "GN=OR4F5" in diann


def test_diann_mgyg_gene_first_word():
    """MGYG must have accession as gene (first word of description)."""
    h = "MGYG000071032_00033 Adaptive-response sensory-kinase SasA"
    p = parse_header(h)
    diann = p.format_for_diann()
    desc = diann.split(None, 1)[1]
    assert desc.startswith("MGYG000071032_00033 ")
    assert "GN=MGYG000071032_00033" in diann


def test_diann_ncbi_gene_first_word():
    """NCBI microbe with gene= must have gene as first word."""
    h = (
        "CIFR|WP_003023861.1|JMV71_RS00005 gene=dnaA chromosomal "
        "replication initiator protein DnaA OS=Citrobacter freundii"
    )
    p = parse_header(h)
    diann = p.format_for_diann()
    desc = diann.split(None, 1)[1]
    assert desc.startswith("dnaA ")
    assert "GN=dnaA" in diann


def test_diann_accession_never_modified():
    """Accession (first token) must NEVER be modified."""
    headers = [
        "sp|P02768|ALBU_HUMAN Albumin GN=ALB",
        "var-nfe|P02768|ALB_Glu177Asp|4-73409403-A-T sp|P02768|ALBU_HUMAN Albumin GN=ALB",
        "ENSP00000493376.2|ENST00000641515.2|ENSG00000186092.7|OTTHUMG00000001094.4|OTTHUMT00000003223.4|OR4F5-201|OR4F5|326",
        "MGYG000071032_00033 Adaptive-response sensory-kinase SasA",
        "CIFR|WP_003023861.1|JMV71_RS00005 gene=dnaA description",
    ]
    for h in headers:
        p = parse_header(h)
        diann = p.format_for_diann()
        assert diann.split()[0] == p.accession, f"Accession modified for: {h[:50]}"


# ═══════════════════════════════════════════════════════════════
# HEAL FASTA TESTS
# ═══════════════════════════════════════════════════════════════


def test_heal_fasta_removes_genscan():
    """GENSCAN entries should be removed."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(">sp|P02768|ALBU_HUMAN Albumin GN=ALB\nMKWVTFISLLFLFSSAYS\n")
        f.write(">GENSCAN00000033973 pep scaffold:GRCh38\nMKWVTF\n")
        f.write(">MGYG000071032_00033 kinase\nMKWVTFISL\n")
        inp = f.name

    out = inp + ".healed"
    try:
        stats = heal_fasta(inp, out, target_engine="diann")
        assert stats.total == 3
        assert stats.kept == 2
        assert stats.skipped == 1
        # Verify GENSCAN not in output
        with open(out) as f:
            content = f.read()
        assert "GENSCAN" not in content
        assert "ALBU_HUMAN" in content
        assert "MGYG" in content
    finally:
        os.unlink(inp)
        os.unlink(out)


def test_heal_fasta_all_have_gn():
    """Every entry in healed output must have GN=."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(">sp|P02768|ALBU_HUMAN Albumin GN=ALB\nMKWVTF\n")
        f.write(">MGYG000071032_00033 kinase\nMKWVTF\n")
        f.write(
            ">ENSP00000493376.2|ENST00000641515.2|ENSG00000186092.7|A|B|OR4F5-201|OR4F5|326\nMKWVTF\n"
        )
        inp = f.name

    out = inp + ".healed"
    try:
        heal_fasta(inp, out, target_engine="diann")
        with open(out) as f:
            for line in f:
                if line.startswith(">"):
                    assert "GN=" in line, f"Missing GN= in: {line.rstrip()}"
    finally:
        os.unlink(inp)
        os.unlink(out)


# ═══════════════════════════════════════════════════════════════
# RUN ALL TESTS
# ═══════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════
# ENGINE-AWARE FASTA PREPARATION TESTS
# ═══════════════════════════════════════════════════════════════


def _write_test_fasta(path):
    """Write a small multi-format FASTA for engine tests."""
    with open(path, "w") as f:
        f.write(">sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens GN=ALB\n")
        f.write("MKWVTFISLLFLFSSAYSRGVFRRDAHK\n")
        f.write(">MGYG000001373_00501 Uncharacterized protein\n")
        f.write("MVKRTILGDTF\n")
        f.write(">GMGC10.000_200_695.UNKNOWN hypothetical protein\n")
        f.write("ACDEFGHIKLM\n")


try:
    from fasta_lake.headers import (
        generate_decoys,
        prepare_fasta_for_engine,
    )

    _HAS_PACKAGE_HEADERS = True
except ImportError:
    _HAS_PACKAGE_HEADERS = False


def test_prepare_sage():
    """SAGE: headers healed (GN= added), no decoys."""
    if not _HAS_PACKAGE_HEADERS:
        return  # skip if package not installed
    with tempfile.TemporaryDirectory() as d:
        inp = os.path.join(d, "input.fasta")
        out = os.path.join(d, "sage.fasta")
        _write_test_fasta(inp)
        stats = prepare_fasta_for_engine(inp, out, engine="sage")
        assert stats.target_proteins == 3
        assert stats.decoy_proteins == 0
        assert stats.total_proteins == 3
        with open(out) as f:
            content = f.read()
        assert "GN=ALB" in content
        assert "GN=MGYG000001373_00501" in content
        assert "rev_" not in content  # no decoys


def test_prepare_msfragger():
    """MSFragger: headers healed + reversed decoys with rev_ prefix."""
    if not _HAS_PACKAGE_HEADERS:
        return
    with tempfile.TemporaryDirectory() as d:
        inp = os.path.join(d, "input.fasta")
        out = os.path.join(d, "msfragger.fasta")
        _write_test_fasta(inp)
        stats = prepare_fasta_for_engine(inp, out, engine="msfragger")
        assert stats.target_proteins == 3
        assert stats.decoy_proteins == 3
        assert stats.total_proteins == 6

        headers = []
        sequences = {}
        with open(out) as f:
            current = None
            for line in f:
                if line.startswith(">"):
                    current = line[1:].rstrip()
                    headers.append(current)
                    sequences[current] = ""
                elif current:
                    sequences[current] += line.rstrip()

        assert len(headers) == 6
        # Check rev_ prefix on decoys
        targets = [h for h in headers if not h.startswith("rev_")]
        decoys = [h for h in headers if h.startswith("rev_")]
        assert len(targets) == 3
        assert len(decoys) == 3
        # Check sequences are reversed
        for t_hdr in targets:
            acc = t_hdr.split()[0]
            d_hdr = [h for h in decoys if h.startswith(f"rev_{acc}")][0]
            assert sequences[d_hdr] == sequences[t_hdr][::-1], f"Decoy of {acc} not reversed"


def test_prepare_diann():
    """DIA-NN: gene as first word of description for non-sp/tr."""
    if not _HAS_PACKAGE_HEADERS:
        return
    with tempfile.TemporaryDirectory() as d:
        inp = os.path.join(d, "input.fasta")
        out = os.path.join(d, "diann.fasta")
        _write_test_fasta(inp)
        stats = prepare_fasta_for_engine(inp, out, engine="diann")
        assert stats.target_proteins == 3
        assert stats.decoy_proteins == 0

        with open(out) as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                hdr = line[1:].rstrip()
                if "MGYG" in hdr:
                    # Gene should be first word of description
                    parts = hdr.split()
                    assert parts[1] == "MGYG000001373_00501", f"Gene not first word: {parts[1]}"
                if "GMGC" in hdr:
                    parts = hdr.split()
                    assert parts[1] == "GMGC10.000_200_695.UNKNOWN", (
                        f"Gene not first word: {parts[1]}"
                    )


def test_generate_decoys_standalone():
    """generate_decoys() reverses sequences with prefix."""
    if not _HAS_PACKAGE_HEADERS:
        return
    with tempfile.TemporaryDirectory() as d:
        inp = os.path.join(d, "input.fasta")
        out = os.path.join(d, "decoy.fasta")
        _write_test_fasta(inp)
        n = generate_decoys(inp, out, prefix="rev_")
        assert n == 3
        with open(out) as f:
            lines = f.readlines()
        # 3 targets + 3 decoys = 6 entries = 12 lines
        assert len(lines) == 12
        # First 6 lines = targets, last 6 = decoys
        assert lines[6].startswith(">rev_sp|P02768|ALBU_HUMAN")


def test_prepare_unknown_engine():
    """Unknown engine raises ValueError."""
    if not _HAS_PACKAGE_HEADERS:
        return
    with tempfile.TemporaryDirectory() as d:
        inp = os.path.join(d, "input.fasta")
        out = os.path.join(d, "out.fasta")
        _write_test_fasta(inp)
        try:
            prepare_fasta_for_engine(inp, out, engine="maxquant")
            assert False, "Should have raised ValueError"
        except ValueError as e:
            assert "maxquant" in str(e)


if __name__ == "__main__":
    test_functions = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    failed = 0
    for test in test_functions:
        try:
            test()
            print(f"  PASS  {test.__name__}")
            passed += 1
        except Exception as e:
            print(f"  FAIL  {test.__name__}: {e}")
            failed += 1

    print(f"\n{'=' * 40}")
    print(f"  {passed} passed, {failed} failed, {passed + failed} total")
    if failed == 0:
        print("  ALL TESTS PASSED")
    print(f"{'=' * 40}")
