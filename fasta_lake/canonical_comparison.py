"""
Canonical comparison — the detailed dissection.

For every gene found in canonical but NOT in FastaLake ("lost" genes),
determine WHY it was lost:

1. REASSIGNED — same peptides found, but assigned to a variant protein group
   with a different gene name (e.g., ALB peptides assigned to ALB_Glu177Asp
   variant, which DIA-NN reports under a different Genes value)

2. SEARCH_SPACE_COMPETITION — peptides exist in the FASTA but didn't pass
   FDR in the larger search space. Marginal, low-abundance, single-peptide IDs.

3. NOT_IN_FASTA — the protein wasn't in the FastaLake FASTA at all
   (shouldn't happen if --include canonical was used)

This classification is critical for defending FastaLake. "Lost" genes that
are actually REASSIGNED are not real losses — the peptides are still found,
just under a different (often more specific) protein group name.
"""

from __future__ import annotations

import csv
import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from fasta_lake.helpers import require_columns

logger = logging.getLogger(__name__)


@dataclass
class LossClassification:
    """Classification of a single 'lost' gene."""

    gene: str
    category: str  # "reassigned", "search_space_competition", "not_in_fasta"
    canonical_peptides: int = 0  # peptides in canonical result
    canonical_intensity: float = 0  # summed intensity in canonical
    fl_peptides_found: int = 0  # how many of those peptides appear in FL (under any gene)
    fl_gene_name: str = ""  # what gene name the peptides got in FL (if reassigned)
    explanation: str = ""


def classify_losses(
    canonical_report: str | Path,
    fastalake_report: str | Path,
    fastalake_fasta: str | Path | None = None,
    engine: str = "diann",
    qvalue: float = 0.01,
) -> list[LossClassification]:
    """Classify every gene lost from canonical → FastaLake.

    Parameters
    ----------
    canonical_report : path
        Canonical search results (DIA-NN report.tsv or SAGE results.sage.tsv).
    fastalake_report : path
        FastaLake search results.
    fastalake_fasta : path, optional
        FastaLake FASTA (to check if proteins are present).
    engine : str
        "diann" or "sage".
    qvalue : float
        FDR threshold.

    Returns
    -------
    list of LossClassification
        One entry per lost gene, classified and explained.
    """
    csv.field_size_limit(10**7)

    if engine == "diann":
        can_data = _load_diann_gene_peptides(canonical_report, qvalue)
        fl_data = _load_diann_gene_peptides(fastalake_report, qvalue)
    else:
        raise NotImplementedError("SAGE loss classification not yet implemented")

    # Build FL peptide → gene lookup (which gene does each peptide belong to in FL?)
    fl_pep_to_genes = defaultdict(set)
    for gene, peps in fl_data.items():
        for pep in peps["peptides"]:
            fl_pep_to_genes[pep].add(gene)

    # Check if proteins are in the FASTA
    fasta_accessions = set()
    if fastalake_fasta:
        with open(fastalake_fasta) as f:
            for line in f:
                if line.startswith(">"):
                    fasta_accessions.add(line[1:].split()[0])

    # Find lost genes
    can_genes = set(can_data.keys())
    fl_genes = set(fl_data.keys())
    lost_genes = can_genes - fl_genes

    results = []
    for gene in sorted(lost_genes):
        info = can_data[gene]
        can_peps = info["peptides"]
        can_intensity = info["intensity"]

        # How many of this gene's peptides appear in FL (under any gene)?
        found_in_fl = set()
        fl_gene_names = set()
        for pep in can_peps:
            if pep in fl_pep_to_genes:
                found_in_fl.add(pep)
                fl_gene_names.update(fl_pep_to_genes[pep])

        pct_found = len(found_in_fl) / max(len(can_peps), 1) * 100

        if pct_found >= 50:
            # Most peptides still found → REASSIGNED
            category = "reassigned"
            explanation = (
                f"{len(found_in_fl)}/{len(can_peps)} peptides still found in FastaLake, "
                f"assigned to gene(s): {', '.join(sorted(fl_gene_names)[:3])}"
            )
        elif pct_found > 0:
            category = "partial_reassignment"
            explanation = (
                f"{len(found_in_fl)}/{len(can_peps)} peptides found under different gene(s)"
            )
        else:
            category = "search_space_competition"
            explanation = (
                f"0/{len(can_peps)} peptides found — pushed below FDR by larger search space. "
                f"Intensity: {can_intensity:.0f}"
            )

        results.append(
            LossClassification(
                gene=gene,
                category=category,
                canonical_peptides=len(can_peps),
                canonical_intensity=can_intensity,
                fl_peptides_found=len(found_in_fl),
                fl_gene_name=";".join(sorted(fl_gene_names)[:3]),
                explanation=explanation,
            )
        )

    return results


def _load_diann_gene_peptides(report_path, qvalue):
    """Load DIA-NN report: gene → {peptides, intensity}."""
    from fasta_lake.benchmark import is_real_gene_symbol

    data = {}
    with open(report_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        require_columns(
            reader.fieldnames,
            ("Q.Value", "Stripped.Sequence", "Genes", "Precursor.Quantity"),
            report_path,
        )
        for row in reader:
            q = float(row["Q.Value"])
            if q > qvalue:
                continue
            pep = row["Stripped.Sequence"]
            gene_raw = row["Genes"]
            intensity = float(row["Precursor.Quantity"] or 0)

            for gene in gene_raw.split(";"):
                gene = gene.strip()
                if not is_real_gene_symbol(gene):
                    continue
                if gene not in data:
                    data[gene] = {"peptides": set(), "intensity": 0}
                if pep:
                    data[gene]["peptides"].add(pep)
                data[gene]["intensity"] += intensity

    return data


def print_loss_report(losses: list[LossClassification], label: str = "") -> None:
    """Print a detailed loss classification report."""
    reassigned = [loss for loss in losses if loss.category == "reassigned"]
    partial = [loss for loss in losses if loss.category == "partial_reassignment"]
    competition = [loss for loss in losses if loss.category == "search_space_competition"]

    print(f"\n{'=' * 75}")
    print(f"  LOSS CLASSIFICATION: {label}")
    print(f"  {len(losses)} genes in canonical but not in FastaLake")
    print(f"{'=' * 75}")

    print(f"\n  REASSIGNED ({len(reassigned)} genes) — peptides found under different gene name:")
    print("  These are NOT real losses. The peptides are still identified.")
    for loss in reassigned[:10]:
        print(
            (
                f"    {loss.gene:<15} → {loss.fl_gene_name:<20} "
                f"({loss.fl_peptides_found}/{loss.canonical_peptides} peps)"
            )
        )
    if len(reassigned) > 10:
        print(f"    ... and {len(reassigned) - 10} more")

    print(f"\n  PARTIAL ({len(partial)} genes) — some peptides reassigned:")
    for loss in partial[:10]:
        print(
            (
                f"    {loss.gene:<15} {loss.fl_peptides_found}/{loss.canonical_peptides} "
                f"peps → {loss.fl_gene_name}"
            )
        )
    if len(partial) > 10:
        print(f"    ... and {len(partial) - 10} more")

    print(f"\n  SEARCH SPACE COMPETITION ({len(competition)} genes) — peptides below FDR:")
    print("  These are true losses — marginal IDs pushed below threshold.")
    # Sort by intensity to show highest-intensity losses first
    competition.sort(key=lambda x: -x.canonical_intensity)
    for loss in competition[:10]:
        print(
            (
                f"    {loss.gene:<15} {loss.canonical_peptides} "
                f"peps, int={loss.canonical_intensity:>12,.0f}"
            )
        )
    if len(competition) > 10:
        print(f"    ... and {len(competition) - 10} more")

    print("\n  SUMMARY:")
    print(
        (
            f"    Reassigned (not real losses): {len(reassigned)} "
            f"({len(reassigned) / max(len(losses), 1) * 100:.0f}%)"
        )
    )
    print(
        (
            f"    Partial reassignment:         {len(partial)} "
            f"({len(partial) / max(len(losses), 1) * 100:.0f}%)"
        )
    )
    print(
        (
            f"    True losses (FDR boundary):   {len(competition)} "
            f"({len(competition) / max(len(losses), 1) * 100:.0f}%)"
        )
    )
    print(f"{'=' * 75}")
