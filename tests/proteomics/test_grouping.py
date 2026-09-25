"""Unit tests for proteomics gene-rollup (peptide-voting). Load-bearing — keep strict."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from fasta_lake.proteomics import grouping as g

# Two canonical genes with deliberately distinct sequences (yield discriminative peptides).
G1 = "MEEISAMPLEPEPTIDEKWITNESSFRAGMENTKQUIETENDINGSEQ"
G2 = "MQRDIFFERENTGENEKUNIQUESTRETCHAANRNOVELDOMAINKLAST"
CANON = [
    ("sp|P00001|G1_HUMAN Sample protein one OS=Homo sapiens OX=9606 GN=G1 PE=1 SV=1", G1),
    ("sp|P00002|G2_HUMAN Different protein two OS=Homo sapiens OX=9606 GN=G2 PE=1 SV=1", G2),
]


def _idx():
    import os
    import tempfile

    fd, p = tempfile.mkstemp(suffix=".fasta")
    os.close(fd)
    with open(p, "w") as f:
        for h, s in CANON:
            f.write(f">{h}\n{s}\n")
    return g.build_gene_index(p, min_len=7, max_missed=1)


def test_tryptic_cuts_after_KR_and_folds_IL():
    peps = g.tryptic_peptides("AAAAAAAKCCCCCCCRDDDDDDD", min_len=3, max_missed=0)
    assert "AAAAAAAK" in peps and "CCCCCCCR" in peps and "DDDDDDD" in peps
    # I/L folded: a peptide with L is stored with I
    peps2 = g.tryptic_peptides("LLLLLLLK", min_len=3, max_missed=0)
    assert "IIIIIIIK" in peps2


def test_accession_detection():
    assert g.is_accession_gene("A8K5A4")  # old 6-char
    assert g.is_accession_gene("A0A385KNS5")  # new 10-char A0A...
    assert g.is_accession_gene("P00450")
    assert not g.is_accession_gene("CP")
    assert not g.is_accession_gene("CD44")
    assert not g.is_accession_gene("HLA-A")


def test_protein_name_parsing():
    h = "tr|A8K5A4|A8K5A4_HUMAN Ceruloplasmin OS=Homo sapiens OX=9606 GN=A8K5A4 PE=2 SV=1"
    assert g.parse_protein_name(h) == "Ceruloplasmin"
    assert g.parse_gene(h) == "A8K5A4"


def test_isoform_votes_correct_gene():
    idx = _idx()
    # near-copy of G1 with GN=<accession> -> must roll up to G1 by peptide vote
    iso_header = "tr|A9XXX9|A9XXX9_HUMAN Sample protein one isoform GN=A9XXX9"
    gene, method = g.assign_gene(G1, iso_header, idx)
    assert gene == "G1", (gene, method)
    assert method == "peptide_vote"


def test_name_fallback_when_no_peptide_votes():
    idx = _idx()
    # a sequence with no shared peptides, but name matches a canonical name exactly
    gene, method = g.assign_gene("ZZZZZZZWZZZZ", "tr|Q9|Q9_HUMAN Sample protein one GN=Q9", idx)
    assert gene == "G1" and method == "name_match", (gene, method)


def test_unmapped_when_nothing_matches():
    idx = _idx()
    # no peptide votes, name doesn't match, and existing GN is itself an accession (fake)
    h = "tr|A0A385KNS5|A0A385KNS5_HUMAN Totally novel thing GN=A0A385KNS5"
    gene, method = g.assign_gene("ZZZZZZZWZZZZ", h, idx)
    assert gene is None and method == "unmapped", (gene, method)


def test_real_gene_is_kept():
    idx = _idx()
    gene, method = g.assign_gene("ZZZZZZZWZZZZ", "tr|Q9|Q9_HUMAN Novel GN=REALSYM", idx)
    assert gene == "REALSYM" and method == "kept_existing"


def test_ig_repertoire_class():
    idx = _idx()
    # no canonical votes + Ig variable-region description -> ig_repertoire (gene None, not forced)
    h = "tr|A0A068LKQ8|A0A068LKQ8_HUMAN Ig heavy chain variable region (Fragment) GN=A0A068LKQ8"
    gene, method = g.assign_gene("ZZWZZWZZWZZ", h, idx)
    assert gene is None and method == "ig_repertoire", (gene, method)
    assert g.is_ig(h) and not g.is_ig("tr|X|X_HUMAN Ceruloplasmin GN=CP")


def test_paralog_ambiguity_flagged():
    idx = _idx()
    # chimeric query = G1 + G2 -> votes split across both -> flagged ambiguous
    gene, method = g.assign_gene(G1 + G2, "tr|CHIM|CHIM_HUMAN Chimera GN=CHIM", idx)
    assert method == "peptide_vote_ambiguous", (gene, method)
    assert gene in ("G1", "G2")


if __name__ == "__main__":
    import traceback

    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    ok = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS {fn.__name__}")
            ok += 1
        except Exception:
            print(f"FAIL {fn.__name__}")
            traceback.print_exc()
    print(f"\n{ok}/{len(fns)} passed")
    sys.exit(0 if ok == len(fns) else 1)
