"""Gene rollup for proteomics mode (Point 6) — accession -> real gene symbol.

THE LOAD-BEARING STEP. TrEMBL isoform entries frequently carry GN=<accession>
(no real gene symbol), so a search engine counts each accession as a distinct
"gene" — inflating gene counts and breaking protein grouping (e.g. three
Ceruloplasmin isoforms A8K5A4/B7Z5Q2/A5PL27 should ALL roll up to CP).

Method = peptide-voting (sequence-grounded, not name-guessing):
  1. Digest every canonical SwissProt protein (GN= reliable) into tryptic peptides.
     Build {peptide -> gene} keeping ONLY peptides discriminative for a single gene.
  2. For a query protein, digest it, tally votes from discriminative peptides,
     assign the plurality gene. Deterministic tie-break: more votes, then more
     total index peptides for that gene, then lexicographic.
  3. Fallbacks when no peptide votes: exact protein-name match against canonical
     name->gene; else leave the existing GN= (flagged unmapped).

Why peptide-voting: a TrEMBL isoform is near-identical to its canonical protein,
so the overwhelming majority of its tryptic peptides are gene-unique and vote the
correct gene. Paralog/shared regions are excluded automatically (non-discriminative
peptides never enter the index). No external idmapping file needed.

Deterministic: same inputs -> same output. No randomness.
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

_GN_RE = re.compile(r"\bGN=(\S+)")
# UniProt accession (old 6-char + new 10-char A0A...); used to detect GN=<accession>
_ACC_RE = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")
# Immunoglobulin variable-region / Ig-domain entries (the somatic repertoire) — a CLASS,
# not a fake gene. These legitimately don't map to a canonical SwissProt gene symbol.
_IG_RE = re.compile(
    r"(immunoglobulin|\bIg[- ](?:like|heavy|kappa|lambda|gamma|mu|alpha)|"
    r"\bIg heavy chain|variable region|\bIG[HKL][VDJ]?\b)",
    re.I,
)


def is_ig(header: str) -> bool:
    """Test whether the header matches the configured immunoglobulin pattern."""
    return bool(_IG_RE.search(header))


def iter_fasta(path: str | Path):
    """Yield (header_without_gt, sequence)."""
    h, buf = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(buf)
                h, buf = line[1:].rstrip("\n"), []
            else:
                buf.append(line.strip())
    if h is not None:
        yield h, "".join(buf)


def parse_gene(header: str) -> str | None:
    """Return the first GN= gene token, or None when the header has none."""
    m = _GN_RE.search(header)
    return m.group(1) if m else None


def parse_protein_name(header: str) -> str | None:
    """Protein description between the id token and the first KEY= tag.

    e.g. 'sp|P00450|CERU_HUMAN Ceruloplasmin OS=Homo sapiens ...' -> 'Ceruloplasmin'
    """
    parts = header.split(None, 1)
    if len(parts) < 2:
        return None
    rest = parts[1]
    rest = re.split(r"\s+[A-Z][A-Za-z]*=", rest, maxsplit=1)[0]
    return rest.strip() or None


def is_accession_gene(gene: str | None) -> bool:
    """True if GN= is actually a UniProt accession (i.e. a fake gene)."""
    return bool(gene) and bool(_ACC_RE.match(gene))


def tryptic_peptides(
    seq: str, min_len: int = 7, max_len: int = 40, max_missed: int = 1
) -> set[str]:
    """Cut after K/R (no proline rule — matches the run's --cut 'K*,R*'). I/L folded.

    I/L folding (both -> 'I') makes voting robust to isoleucine/leucine ambiguity,
    which MS cannot distinguish anyway.
    """
    seq = seq.replace("L", "I")
    # cut sites: after every K or R
    cuts = [0]
    for i, aa in enumerate(seq):
        if aa in "KR":
            cuts.append(i + 1)
    if cuts[-1] != len(seq):
        cuts.append(len(seq))
    peps: set[str] = set()
    n = len(cuts) - 1
    for i in range(n):
        for m in range(max_missed + 1):
            j = i + 1 + m
            if j > n:
                break
            pep = seq[cuts[i] : cuts[j]]
            if min_len <= len(pep) <= max_len:
                peps.add(pep)
    return peps


@dataclass
class GeneIndex:
    pep2gene: dict[str, str]  # discriminative peptide -> gene
    name2gene: dict[str, str]  # normalized protein name -> gene (unambiguous only)
    gene_pep_count: dict[str, int]  # gene -> #discriminative peptides (tie-break weight)


def build_gene_index(
    canonical_fasta: str | Path, min_len: int = 7, max_missed: int = 1
) -> GeneIndex:
    """Build the discriminative peptide->gene index from canonical SwissProt."""
    pep_genes: dict[str, set[str]] = defaultdict(set)
    name_genes: dict[str, set[str]] = defaultdict(set)
    for header, seq in iter_fasta(canonical_fasta):
        gene = parse_gene(header)
        if not gene or is_accession_gene(gene):
            continue  # only trust real symbols as ground truth
        for pep in tryptic_peptides(seq, min_len=min_len, max_missed=max_missed):
            pep_genes[pep].add(gene)
        name = parse_protein_name(header)
        if name:
            name_genes[_norm_name(name)].add(gene)
    pep2gene = {p: next(iter(g)) for p, g in pep_genes.items() if len(g) == 1}
    name2gene = {n: next(iter(g)) for n, g in name_genes.items() if len(g) == 1}
    gene_pep_count: Counter = Counter(pep2gene.values())
    return GeneIndex(pep2gene, name2gene, dict(gene_pep_count))


def _norm_name(name: str) -> str:
    # strip isoform/fragment qualifiers so 'Ceruloplasmin (Fragment)' == 'Ceruloplasmin'
    """Lowercase a protein name and remove fragment and isoform qualifiers."""
    name = re.sub(r"\s*\((Fragment|Isoform[^)]*)\)\s*", "", name, flags=re.I)
    return name.strip().lower()


def assign_gene(
    seq: str,
    header: str,
    idx: GeneIndex,
    min_len: int = 7,
    max_missed: int = 1,
    ambig_ratio: float = 0.5,
) -> tuple[str | None, str]:
    """Return (gene, method).

    method in {peptide_vote, peptide_vote_ambiguous, ig_repertoire, name_match,
    kept_existing, unmapped}. For ig_repertoire / unmapped, gene is None (NOT
    counted as a real gene — prevents accession-as-gene inflation; Ig is tracked
    as its own repertoire class).
    """
    votes: Counter = Counter()
    for pep in tryptic_peptides(seq, min_len=min_len, max_missed=max_missed):
        g = idx.pep2gene.get(pep)
        if g:
            votes[g] += 1
    if votes:
        ranked = sorted(
            votes.items(), key=lambda kv: (-kv[1], -idx.gene_pep_count.get(kv[0], 0), kv[0])
        )
        winner, top = ranked[0]
        # paralog ambiguity: runner-up wins a comparable share of votes (e.g. C2 vs CFB)
        if len(ranked) > 1 and ranked[1][1] >= ambig_ratio * top:
            return winner, "peptide_vote_ambiguous"
        return winner, "peptide_vote"
    # No discriminative votes. Immunoglobulin variable regions land here legitimately —
    # the somatic repertoire isn't in canonical SwissProt with a standard gene symbol.
    if is_ig(header):
        return None, "ig_repertoire"
    name = parse_protein_name(header)
    if name:
        g = idx.name2gene.get(_norm_name(name))
        if g:
            return g, "name_match"
    existing = parse_gene(header)
    if existing and not is_accession_gene(existing):
        return existing, "kept_existing"
    return None, "unmapped"


def heal_fasta_genes(
    input_fasta: str | Path,
    canonical_fasta: str | Path,
    output_fasta: str | Path,
    mapping_tsv: str | Path | None = None,
    min_len: int = 7,
    max_missed: int = 1,
) -> dict:
    """Rewrite GN= in input_fasta with peptide-voted real gene symbols. Returns stats."""
    idx = build_gene_index(canonical_fasta, min_len=min_len, max_missed=max_missed)
    stats: Counter = Counter()
    rows = []
    with open(output_fasta, "w") as out:
        for header, seq in iter_fasta(input_fasta):
            old = parse_gene(header)
            gene, method = assign_gene(seq, header, idx, min_len=min_len, max_missed=max_missed)
            stats[method] += 1
            if gene:
                new_header = _set_gene(header, gene)
                changed = old != gene
                stats["changed"] += int(changed)
            else:
                new_header = header
            out.write(f">{new_header}\n{seq}\n")
            acc = header.split("|")[1] if header.startswith(("sp|", "tr|")) else header.split()[0]
            rows.append((acc, old or "", gene or "", method))
    if mapping_tsv:
        with open(mapping_tsv, "w") as m:
            m.write("accession\told_gene\tnew_gene\tmethod\n")
            for r in rows:
                m.write("\t".join(r) + "\n")
    stats["total"] = sum(
        stats[k]
        for k in (
            "peptide_vote",
            "peptide_vote_ambiguous",
            "ig_repertoire",
            "name_match",
            "kept_existing",
            "unmapped",
        )
    )
    stats["index_peptides"] = len(idx.pep2gene)
    stats["index_genes"] = len(idx.gene_pep_count)
    return dict(stats)


def _set_gene(header: str, gene: str) -> str:
    """Replace the first GN= token, or append it when absent."""
    if _GN_RE.search(header):
        return _GN_RE.sub(f"GN={gene}", header, count=1)
    return f"{header} GN={gene}"
