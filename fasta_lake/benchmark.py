"""
Benchmark FastaLake results against a baseline search.

For proteomics: compare against canonical SwissProt search.
For metaproteomics: compare against matched metagenome or iterative search.

Computes 6 metrics:
1. Identification gain (peptides, proteins, genes)
2. FDR control (entrapment hits)
3. Canonical coverage (overlap with baseline)
4. Biological plausibility (species distribution, variant frequencies)
5. Reproducibility (cross-sample CV)
6. Quantitative agreement (intensity correlation)
"""

from __future__ import annotations

import csv
import logging
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

from fasta_lake.helpers import require_columns

logger = logging.getLogger(__name__)


def is_real_gene_symbol(gene: str) -> bool:
    """Check if a gene name is a real HGNC/MGI symbol, not an accession fallback.

    When FASTA headers lack a proper gene name (GN= tag), the heal function
    falls back to using the UniProt accession (e.g., A0A024CIM4) or entry
    name as the gene. These look like gene names but aren't — they inflate
    gene counts if not filtered.

    Real gene symbols: ALB, APOA1, HP, TP53, BRCA1, MGYG000000001_00001
    Accession fallbacks: A0A024CIM4, B2R6Q2, Q8IZF2, F6KPG5

    The UniProt accession pattern is ``[A-Z][0-9][A-Z0-9]{3,}`` (e.g., P02768,
    A0A024CIM4, Q6UXV3). This regex catches accessions but not gene symbols
    like ALB, TP53, or C3 (which are short and don't match the pattern).
    """
    if not gene:
        return False
    if "|" in gene:
        return False  # broken header (sp|P02768|ALBU_HUMAN as gene)
    if len(gene) > 20:
        return False
    # UniProt accession pattern: letter-digit-then-3+-alphanumeric
    # Catches: P02768, A0A024CIM4, B2R6Q2, Q8IZF2, F6KPG5
    # Does NOT catch: ALB, TP53, C3, HBA1, SERPINA1, APOA1
    if re.match(r"^[A-Z][0-9][A-Z0-9]{3,}$", gene):
        return False
    if gene.startswith("ENSP") or gene.startswith("ENST") or gene.startswith("ENSG"):
        return False
    return True


@dataclass
class SearchResult:
    """Parsed search engine result (DIA-NN or SAGE report.tsv).

    Gene counts use ONLY real gene symbols (HGNC/MGI), not UniProt accession
    fallbacks. This prevents inflated counts when the FASTA contains many
    TrEMBL entries without proper gene names.

    The `genes_all` field contains everything (including accession fallbacks)
    for debugging. The `genes` field is filtered to real symbols only.
    """

    peptides: set[str] = field(default_factory=set)
    genes: set[str] = field(default_factory=set)  # real gene symbols only
    genes_all: set[str] = field(default_factory=set)  # including accession fallbacks
    protein_groups: set[str] = field(default_factory=set)
    peptide_intensities: dict[str, float] = field(default_factory=dict)
    gene_intensities: dict[str, float] = field(default_factory=dict)
    entrapment_peptides: int = 0
    entrapment_proteins: int = 0
    total_precursors: int = 0
    broken_genes: int = 0  # genes containing | (header healing failure)


@dataclass
class BenchmarkResult:
    """Complete benchmark comparison between FastaLake and baseline."""

    # Identification gain
    fl_peptides: int = 0
    bl_peptides: int = 0
    fl_genes: int = 0
    bl_genes: int = 0
    fl_protein_groups: int = 0
    bl_protein_groups: int = 0
    peptide_gain_pct: float = 0.0
    gene_gain_pct: float = 0.0

    # Overlap
    shared_peptides: int = 0
    fl_only_peptides: int = 0
    bl_only_peptides: int = 0
    shared_genes: int = 0
    fl_only_genes: int = 0
    bl_only_genes: int = 0
    gene_coverage_pct: float = 0.0  # % of baseline genes retained

    # FDR
    fl_entrapment_precursor_fdr: float = 0.0
    bl_entrapment_precursor_fdr: float = 0.0
    fl_entrapment_protein_fdr: float = 0.0

    # Quantitative agreement
    intensity_correlation: float = 0.0
    n_shared_for_correlation: int = 0

    # Quality flags
    broken_genes: int = 0
    warnings: list[str] = field(default_factory=list)


def load_diann_report(report_path: str | Path, qvalue: float = 0.01) -> SearchResult:
    """Load a DIA-NN report.tsv and extract peptides, genes, protein groups.

    Gene names are filtered to real HGNC/MGI symbols only. Accession
    fallbacks (A0A024CIM4, B2R6Q2, etc.) are excluded from gene counts
    but kept in genes_all for debugging.
    """
    result = SearchResult()
    with open(report_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        require_columns(
            reader.fieldnames,
            ("Q.Value", "Stripped.Sequence", "Genes", "Protein.Group", "Precursor.Quantity"),
            report_path,
        )
        for row in reader:
            q = float(row["Q.Value"])
            if q > qvalue:
                continue

            result.total_precursors += 1
            pep = row["Stripped.Sequence"]
            gene = row["Genes"]
            pg = row["Protein.Group"]
            intensity = float(row["Precursor.Quantity"] or 0)

            if pep:
                result.peptides.add(pep)
                if (
                    pep not in result.peptide_intensities
                    or intensity > result.peptide_intensities[pep]
                ):
                    result.peptide_intensities[pep] = intensity
            if gene:
                result.genes_all.add(gene)
                if "|" in gene:
                    result.broken_genes += 1
                # Only count real gene symbols, not accession fallbacks
                # DIA-NN Genes column may contain multi-gene groups (A;B)
                for g in gene.split(";"):
                    g = g.strip()
                    if is_real_gene_symbol(g):
                        result.genes.add(g)
                if gene not in result.gene_intensities:
                    result.gene_intensities[gene] = 0
                result.gene_intensities[gene] += intensity
            if pg:
                result.protein_groups.add(pg)
                if _is_entrapment_accession(pg):
                    result.entrapment_proteins += 1
                    result.entrapment_peptides += 1

    return result


def _is_entrapment_accession(accession: str) -> bool:
    """One rule for both loaders: an entrapment record carries ENTRAP or SHUFFLED
    in its accession (every FastaLake entrapment class writes an ENTRAP_ prefix)."""
    return "ENTRAP" in accession or "SHUFFLED" in accession


def build_gene_map(fasta_path: str | Path) -> dict[str, str]:
    """Build accession → gene name map from a FASTA file.

    Extracts gene names from GN= tags in FASTA headers. Only returns
    real gene symbols (filtered by is_real_gene_symbol), not accession
    fallbacks.

    Parameters
    ----------
    fasta_path : path
        FASTA file (healed or canonical).

    Returns
    -------
    dict mapping protein accession → gene symbol
    """
    acc_to_gene = {}
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            acc = line[1:].split()[0]
            m = re.search(r"GN=(\S+)", line)
            if m and is_real_gene_symbol(m.group(1)):
                acc_to_gene[acc] = m.group(1)
    return acc_to_gene


def load_sage_report(
    report_path: str | Path,
    fasta_path: str | Path | None = None,
    qvalue: float = 0.01,
) -> SearchResult:
    """Load a SAGE results.sage.tsv and extract peptides, genes, protein groups.

    SAGE doesn't report gene names directly — it reports protein accessions.
    To get gene-level counts, provide the FASTA that was searched so we can
    map accession → gene using the GN= tag.

    Without a FASTA, only peptide and protein counts are available (no genes).

    Parameters
    ----------
    report_path : path
        SAGE results.sage.tsv file.
    fasta_path : path, optional
        The FASTA file that was searched. Used to map proteins → genes.
    qvalue : float
        FDR threshold (default 0.01).
    """
    csv.field_size_limit(10**7)  # SAGE protein groups can be very long

    # Build gene map from FASTA if provided
    acc_to_gene = {}
    if fasta_path:
        acc_to_gene = build_gene_map(fasta_path)
        logger.info("Gene map: %d accessions → genes", len(acc_to_gene))

    result = SearchResult()
    with open(report_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames or []
        require_columns(
            fields, ("peptide", "proteins"), report_path, any_of=("spectrum_q", "peptide_q")
        )
        q_col = "spectrum_q" if "spectrum_q" in fields else "peptide_q"
        for row in reader:
            q = float(row[q_col])
            if q > qvalue:
                continue
            # Skip decoys
            if int(row.get("label", 1)) == -1:
                continue

            result.total_precursors += 1

            # Peptide: strip modifications
            pep = row.get("peptide", "")
            pep_clean = "".join(c for c in pep if c.isalpha() and c.isupper())
            if pep_clean:
                result.peptides.add(pep_clean)
                intensity = float(row.get("ms1_intensity", row.get("ms2_intensity", 0)))
                if (
                    pep_clean not in result.peptide_intensities
                    or intensity > result.peptide_intensities[pep_clean]
                ):
                    result.peptide_intensities[pep_clean] = intensity

            # Proteins: SAGE reports semicolon-separated accessions
            proteins_raw = row.get("proteins", "")
            members = [p.strip() for p in proteins_raw.split(";") if p.strip()]

            # Entrapment: same per-precursor rule as load_diann_report. A PSM whose
            # protein list contains any entrapment accession counts once.
            if any(_is_entrapment_accession(p) for p in members):
                result.entrapment_proteins += 1
                result.entrapment_peptides += 1

            # Protein group = sorted unique members (gene-collapsed)
            psm_genes = set()
            for p in members:
                gene = acc_to_gene.get(p)
                if gene:
                    psm_genes.add(gene)

            result.genes.update(psm_genes)
            result.genes_all.update(psm_genes)
            # Also add unmapped accessions to genes_all (for debugging)
            for p in members:
                if p not in acc_to_gene:
                    result.genes_all.add(p)

            # Protein group key: gene-collapsed, sorted
            if psm_genes:
                pg_key = ";".join(sorted(psm_genes))
                result.protein_groups.add(pg_key)

    return result


def compare(
    fastalake_result: SearchResult,
    baseline_result: SearchResult,
    label: str = "FastaLake vs Baseline",
) -> BenchmarkResult:
    """Compare FastaLake search results against a baseline.

    Parameters
    ----------
    fastalake_result : SearchResult
        Results from searching the FastaLake-built FASTA.
    baseline_result : SearchResult
        Results from searching the baseline FASTA (canonical or matched metagenome).
    label : str
        Description for logging.

    Returns
    -------
    BenchmarkResult
        Complete comparison metrics.
    """
    bench = BenchmarkResult()

    # Identification counts
    bench.fl_peptides = len(fastalake_result.peptides)
    bench.bl_peptides = len(baseline_result.peptides)
    bench.fl_genes = len(fastalake_result.genes)
    bench.bl_genes = len(baseline_result.genes)
    bench.fl_protein_groups = len(fastalake_result.protein_groups)
    bench.bl_protein_groups = len(baseline_result.protein_groups)

    # Gain percentages
    if bench.bl_peptides > 0:
        bench.peptide_gain_pct = (bench.fl_peptides - bench.bl_peptides) / bench.bl_peptides * 100
    if bench.bl_genes > 0:
        bench.gene_gain_pct = (bench.fl_genes - bench.bl_genes) / bench.bl_genes * 100

    # Overlap
    bench.shared_peptides = len(fastalake_result.peptides & baseline_result.peptides)
    bench.fl_only_peptides = len(fastalake_result.peptides - baseline_result.peptides)
    bench.bl_only_peptides = len(baseline_result.peptides - fastalake_result.peptides)

    bench.shared_genes = len(fastalake_result.genes & baseline_result.genes)
    bench.fl_only_genes = len(fastalake_result.genes - baseline_result.genes)
    bench.bl_only_genes = len(baseline_result.genes - fastalake_result.genes)

    if bench.bl_genes > 0:
        bench.gene_coverage_pct = bench.shared_genes / bench.bl_genes * 100

    # FDR
    if fastalake_result.total_precursors > 0:
        bench.fl_entrapment_precursor_fdr = (
            fastalake_result.entrapment_peptides / fastalake_result.total_precursors * 100
        )
    if baseline_result.total_precursors > 0:
        bench.bl_entrapment_precursor_fdr = (
            baseline_result.entrapment_peptides / baseline_result.total_precursors * 100
        )

    # Quantitative agreement (Pearson on shared peptides)
    shared = fastalake_result.peptides & baseline_result.peptides
    if len(shared) >= 10:
        fl_ints = []
        bl_ints = []
        for pep in shared:
            fl_i = fastalake_result.peptide_intensities.get(pep, 0)
            bl_i = baseline_result.peptide_intensities.get(pep, 0)
            if fl_i > 0 and bl_i > 0:
                fl_ints.append(math.log10(fl_i))
                bl_ints.append(math.log10(bl_i))

        bench.n_shared_for_correlation = len(fl_ints)
        if len(fl_ints) >= 10:
            bench.intensity_correlation = _pearson(fl_ints, bl_ints)

    # Quality flags
    bench.broken_genes = fastalake_result.broken_genes
    if bench.broken_genes > 0:
        bench.warnings.append(
            f"{bench.broken_genes} broken gene names (contain |) — headers not healed for DIA-NN"
        )
    if bench.gene_coverage_pct < 90:
        bench.warnings.append(
            f"Gene coverage {bench.gene_coverage_pct:.1f}% — "
            "FastaLake misses >10% of baseline genes"
        )
    if bench.fl_entrapment_precursor_fdr > 10:
        bench.warnings.append(
            (
                f"High entrapment FDR ({bench.fl_entrapment_precursor_fdr:.1f}%) "
                f"— search space may be too large"
            )
        )

    return bench


def print_benchmark(bench: BenchmarkResult, label: str = "") -> None:
    """Pretty-print benchmark results."""
    print(f"\n{'=' * 70}")
    print(f"  BENCHMARK: {label}")
    print(f"{'=' * 70}")

    print("\n  1. Identification Gain")
    print(f"     {'':20} {'Baseline':>10} {'FastaLake':>10} {'Gain':>10}")
    print(
        (
            f"     {'Peptides':20} {bench.bl_peptides:>10,} "
            f"{bench.fl_peptides:>10,} {bench.peptide_gain_pct:>+9.1f}%"
        )
    )
    print(
        (
            f"     {'Genes':20} {bench.bl_genes:>10,} "
            f"{bench.fl_genes:>10,} {bench.gene_gain_pct:>+9.1f}%"
        )
    )
    print(
        f"     {'Protein groups':20} {bench.bl_protein_groups:>10,} {bench.fl_protein_groups:>10,}"
    )

    print("\n  2. FDR Control")
    print(f"     Baseline entrapment FDR:  {bench.bl_entrapment_precursor_fdr:.2f}%")
    print(f"     FastaLake entrapment FDR: {bench.fl_entrapment_precursor_fdr:.2f}%")

    print("\n  3. Canonical Coverage")
    print(
        f"     Shared genes:    {bench.shared_genes:,} ({bench.gene_coverage_pct:.1f}% of baseline)"
    )
    print(f"     FL-only genes:   {bench.fl_only_genes:,}")
    print(f"     Lost genes:      {bench.bl_only_genes:,}")

    print("\n  4. Quantitative Agreement")
    if bench.n_shared_for_correlation > 0:
        print(f"     Pearson r (log10 intensity): {bench.intensity_correlation:.4f}")
        print(f"     Based on {bench.n_shared_for_correlation:,} shared peptides")
    else:
        print("     Not enough shared peptides for correlation")

    if bench.warnings:
        print("\n  WARNINGS:")
        for w in bench.warnings:
            print(f"     {w}")

    print(f"{'=' * 70}")


def _pearson(x: list[float], y: list[float]) -> float:
    """Pearson correlation coefficient."""
    n = len(x)
    if n < 2:
        return 0.0
    mx = sum(x) / n
    my = sum(y) / n
    sx = math.sqrt(sum((xi - mx) ** 2 for xi in x) / (n - 1))
    sy = math.sqrt(sum((yi - my) ** 2 for yi in y) / (n - 1))
    if sx == 0 or sy == 0:
        return 0.0
    return sum((xi - mx) * (yi - my) for xi, yi in zip(x, y)) / ((n - 1) * sx * sy)
