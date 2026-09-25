"""Layered tier builder for proteomics mode (Point 2) — canonical-anchored by construction.

Builds 'canonical + smart-isoforms + Ig + variants', healed to real gene symbols.
Guarantees:
  * FULL canonical SwissProt backbone is always present (fixes the lost-backbone bug).
  * Isoform layer = only TrEMBL entries that (a) roll up to a real gene AND (b) add >=1
    tryptic peptide NOT already in the canonical proteome (genuine new sequence). This
    trims the low-confidence redundant-TrEMBL bulk that drives FDR-boundary clogging.
  * Ig repertoire kept as its own class (not force-mapped, not counted as fake genes).
  * Personal WGS variants (ens|) kept as the specificity class.
All entries emitted with a corrected GN= (real symbol where determinable), so gene
counts and protein grouping are honest. Deterministic; dedup by sequence.
"""

from __future__ import annotations

from collections import Counter

from fasta_lake.proteomics import grouping as gr


def _canonical_peptide_set(canonical_fasta, min_len=7, max_missed=1) -> set[str]:
    """Pool tryptic peptides from the canonical FASTA under the supplied digest limits."""
    peps: set[str] = set()
    for _, seq in gr.iter_fasta(canonical_fasta):
        peps |= gr.tryptic_peptides(seq, min_len=min_len, max_missed=max_missed)
    return peps


def build_smart_tier(
    canonical_fasta, candidate_fasta, output_fasta, mapping_tsv=None, min_len=7, max_missed=1
) -> dict:
    """canonical backbone + smart isoforms + Ig + variants, gene-healed.

    candidate_fasta = the de-novo-evidence-supported pool (TrEMBL isoforms, ens| variants,
    Ig entries) to layer on top of the full canonical backbone.
    """
    idx = gr.build_gene_index(canonical_fasta, min_len=min_len, max_missed=max_missed)
    canon_peps = _canonical_peptide_set(canonical_fasta, min_len=min_len, max_missed=max_missed)

    seen: set[str] = set()
    stats: Counter = Counter()
    rows = []
    with open(output_fasta, "w") as out:
        # 1) full canonical backbone (always)
        for h, s in gr.iter_fasta(canonical_fasta):
            if s and s not in seen:
                seen.add(s)
                out.write(f">{h}\n{s}\n")
                stats["canonical_backbone"] += 1

        # 2) candidate layers
        for h, s in gr.iter_fasta(candidate_fasta):
            if not s or s in seen:
                stats["dup_or_empty"] += 1
                continue
            acc = h.split("|")[1] if h.startswith(("sp|", "tr|")) else h.split()[0]
            is_variant = (
                acc.lower().startswith("ens") or h.lower().startswith(">ens") or "ens|" in h.lower()
            )
            gene, method = gr.assign_gene(s, h, idx, min_len=min_len, max_missed=max_missed)

            keep, klass = False, None
            if is_variant:
                keep, klass = True, "variant"  # always keep personal WGS variants
            elif method == "ig_repertoire":
                keep, klass = True, "ig_repertoire"  # keep Ig V-region content as a class
            elif gene is not None:
                # smart isoform: must add >=1 genuinely new tryptic peptide
                novel = any(
                    p not in canon_peps
                    for p in gr.tryptic_peptides(s, min_len=min_len, max_missed=max_missed)
                )
                if novel:
                    keep, klass = True, "smart_isoform"
                else:
                    stats["dropped_redundant_isoform"] += 1
            else:
                stats["dropped_unmapped"] += 1

            if keep:
                seen.add(s)
                new_h = gr._set_gene(h, gene) if gene else h
                out.write(f">{new_h}\n{s}\n")
                stats[f"kept_{klass}"] += 1
                rows.append((acc, klass, gene or "", method))

    if mapping_tsv:
        with open(mapping_tsv, "w") as m:
            m.write("accession\tclass\tgene\tmethod\n")
            for r in rows:
                m.write("\t".join(r) + "\n")
    stats["total_written"] = (
        stats["canonical_backbone"]
        + stats["kept_smart_isoform"]
        + stats["kept_ig_repertoire"]
        + stats["kept_variant"]
    )
    return dict(stats)
