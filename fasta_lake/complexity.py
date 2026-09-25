"""
Sample complexity estimation from de novo peptide results.

Maps de novo peptides to a reference database to estimate sample complexity
BEFORE building the per-sample database. Reports metrics that help choose
optimal inference parameters.

Metrics
-------
- **hit_rate**: Fraction of de novo peptides matching the reference DB.
  High hit rate (>0.3) = well-characterized sample, low (<0.1) = novel organisms.
- **species_richness**: Number of distinct species/organisms detected.
  1 = single-species, >100 = complex metaproteomics.
- **peptide_diversity**: Mean unique peptides per protein with evidence.
  High diversity = good coverage, low = sparse evidence.
- **db_coverage**: Fraction of DB proteins with any peptide evidence.
  Indicates how much of the DB is relevant to this sample.

Recommendations
---------------
Based on detected complexity, the module can auto-recommend:

- Single-species (richness=1): ``razor`` strategy
- Low complexity (<10 species): ``razor`` or ``uniform_2pep``
- High complexity (>10 species): ``species_budget``
- Very high complexity (>100 species): ``species_budget`` with clustering
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

from fasta_lake.config import ComplexityConfig
from fasta_lake.helpers import (
    get_species_id,
    map_peptides_to_proteins,
)

logger = logging.getLogger(__name__)


@dataclass
class ComplexityReport:
    """
    Sample complexity estimation report.

    Parameters
    ----------
    hit_rate : float
        Fraction of input peptides that matched the reference DB (0-1).
    species_richness : int
        Number of distinct species/organisms detected.
    peptide_diversity : float
        Mean unique peptides per protein with evidence.
    db_coverage : float
        Fraction of DB proteins with any peptide evidence (0-1).
    total_peptides : int
        Total input peptides analyzed.
    matched_peptides : int
        Peptides that matched at least one protein.
    proteins_with_evidence : int
        Proteins with at least one peptide match.
    species_counts : dict[str, int]
        Proteins per species (species_id -> count).
    recommendation : str
        Auto-recommended inference strategy.
    recommendation_reason : str
        Explanation for the recommendation.
    """

    hit_rate: float = 0.0
    species_richness: int = 0
    peptide_diversity: float = 0.0
    db_coverage: float = 0.0
    total_peptides: int = 0
    matched_peptides: int = 0
    proteins_with_evidence: int = 0
    species_counts: dict[str, int] = field(default_factory=dict)
    unresolved_proteins: int = 0  # proteins with evidence whose species is unknown
    recommendation: str = ""
    recommendation_reason: str = ""

    def summary(self) -> dict:
        """Return summary dict for JSON stats output."""
        result: dict = {
            "hit_rate": round(self.hit_rate, 4),
            "species_richness": self.species_richness,
            "unresolved_proteins": self.unresolved_proteins,
            "peptide_diversity": round(self.peptide_diversity, 2),
            "db_coverage": round(self.db_coverage, 6),
            "total_peptides": self.total_peptides,
            "matched_peptides": self.matched_peptides,
            "proteins_with_evidence": self.proteins_with_evidence,
        }
        if self.recommendation:
            result["recommendation"] = self.recommendation
            result["recommendation_reason"] = self.recommendation_reason
        if self.species_counts:
            # Top 10 species by protein count
            top = sorted(self.species_counts.items(), key=lambda x: -x[1])[:10]
            result["top_species"] = dict(top)
        return result


def estimate_complexity(
    proteins: dict[str, str],
    peptides: dict[str, float] | set[str],
    config: ComplexityConfig | None = None,
) -> ComplexityReport:
    """
    Estimate sample complexity by mapping peptides to a reference database.

    Parameters
    ----------
    proteins : dict[str, str]
        Reference database proteins (protein_id -> sequence).
    peptides : dict[str, float] or set[str]
        De novo peptides (I/L-normalized). If dict, values are scores.
    config : ComplexityConfig or None
        Configuration. If None, uses defaults.

    Returns
    -------
    ComplexityReport
        Complexity metrics and optional strategy recommendation.
    """
    if config is None:
        config = ComplexityConfig()

    peptide_keys = peptides if isinstance(peptides, set) else set(peptides.keys())
    total_peptides = len(peptide_keys)

    if total_peptides == 0 or len(proteins) == 0:
        return ComplexityReport(total_peptides=total_peptides)

    # Map peptides to proteins
    pep2prot, prot2pep = map_peptides_to_proteins(proteins, peptide_keys)

    matched_peptides = len(pep2prot)
    proteins_with_evidence = len(prot2pep)

    # Hit rate
    hit_rate = matched_peptides / total_peptides if total_peptides > 0 else 0.0

    # DB coverage
    db_coverage = proteins_with_evidence / len(proteins) if len(proteins) > 0 else 0.0

    # Peptide diversity (mean peptides per protein)
    peptide_diversity = (
        sum(len(peps) for peps in prot2pep.values()) / proteins_with_evidence
        if proteins_with_evidence > 0
        else 0.0
    )

    # Species richness, over proteins whose species can be resolved from the
    # identifier. Unresolved proteins are counted separately and never pooled
    # into a pseudo-species.
    species_counts: dict[str, int] = {}
    unresolved_proteins = 0
    for protein_id in prot2pep:
        species = get_species_id(protein_id)
        if species is None:
            unresolved_proteins += 1
            continue
        species_counts[species] = species_counts.get(species, 0) + 1

    species_richness = len(species_counts)

    report = ComplexityReport(
        hit_rate=hit_rate,
        species_richness=species_richness,
        peptide_diversity=peptide_diversity,
        db_coverage=db_coverage,
        total_peptides=total_peptides,
        matched_peptides=matched_peptides,
        proteins_with_evidence=proteins_with_evidence,
        species_counts=species_counts,
        unresolved_proteins=unresolved_proteins,
    )

    # Auto-recommend strategy
    if config.auto_recommend:
        _recommend_strategy(report)

    logger.info(
        "Complexity: %.1f%% hit rate, %d species, %.1f pep/prot, %.4f%% DB coverage",
        hit_rate * 100,
        species_richness,
        peptide_diversity,
        db_coverage * 100,
    )

    return report


def _recommend_strategy(report: ComplexityReport) -> None:
    """
    Auto-recommend inference strategy based on complexity metrics.

    Rules (in priority order):

    1. Single species (richness=1) -> razor
    2. Low complexity (<10 species) -> razor
    3. High complexity (10-100 species) -> species_budget
    4. Very high complexity (>100 species) -> species_budget + clustering
    """
    if report.unresolved_proteins > 0:
        report.recommendation = "razor"
        report.recommendation_reason = (
            f"{report.unresolved_proteins} proteins with evidence have no resolvable "
            "species from their identifiers; species_budget would refuse to run. "
            "Use razor, or provide --taxonomy-map for species_budget."
        )
    elif report.species_richness <= 1:
        report.recommendation = "razor"
        report.recommendation_reason = (
            f"Single-species sample ({report.species_richness} species detected). "
            "Razor parsimony is optimal for single-species proteomics."
        )
    elif report.species_richness < 10:
        report.recommendation = "razor"
        report.recommendation_reason = (
            f"Low-complexity sample ({report.species_richness} species). "
            "Razor parsimony works well with limited cross-species ambiguity."
        )
    elif report.species_richness <= 100:
        report.recommendation = "species_budget"
        report.recommendation_reason = (
            f"Moderate-complexity metaproteomics ({report.species_richness} species). "
            "Species-budget parsimony preserves cross-species diversity."
        )
    else:
        report.recommendation = "species_budget"
        report.recommendation_reason = (
            f"High-complexity metaproteomics ({report.species_richness} species). "
            "Species-budget parsimony with clustering recommended for "
            "databases of this size."
        )

    if report.hit_rate < 0.05:
        report.recommendation_reason += (
            f" Warning: very low hit rate ({report.hit_rate:.1%}). "
            "Consider adding more reference databases or assemblies."
        )
