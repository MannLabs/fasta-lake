"""Entrapment control metadata with explicit interpretation limits.

Control membership does not establish absence or exchangeability. Construction
matches are diagnostics unless a defensible estimator and counting-level ratio
are supplied. Source FASTAs must be provided explicitly for copied controls.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

#: Difficulty tiers, ordered easiest to hardest. A class's tier says how hard a
#: test it poses, which bounds how much confidence its passing FDP earns.
TIER_TRIVIAL = "trivial"
TIER_REMOTE = "phylogenetically-remote"
TIER_HOMOLOGOUS = "cross-kingdom-homologous"
TIER_WITHIN_COHORT = "within-cohort"
TIER_USER = "user-defined"

#: What quantity an FDP computed from a class actually estimates.
EST_CONSTRUCTION = "construction-control match diagnostic; no automatic FDR calibration"
EST_SEARCH = "search/scoring false-discovery proportion"
EST_INCOMPLETENESS = "cross-sample sequence-sharing rate (database incompleteness, NOT FDP)"
EST_UNKNOWN = "depends on the user-supplied sequences"


@dataclass(frozen=True)
class EntrapmentClass:
    """One entrapment class and everything needed to interpret its FDP."""

    name: str
    prefix: str
    generation: str  # 'shuffle' (derive from target) or 'copy' (tag a source FASTA)
    difficulty: str
    estimates: str
    description: str
    caveat: str
    default_source: str | None = None
    requires_source: bool = False
    is_definitional_null: bool = False
    recommended_headline: bool = False
    poolable: bool = True

    def as_metadata(self) -> dict:
        """Serialisable metadata for embedding in an FDP report."""
        return asdict(self)


ENTRAPMENT_CLASSES: dict[str, EntrapmentClass] = {
    "shuffled": EntrapmentClass(
        name="shuffled",
        prefix="ENTRAP_SHUF_",
        generation="shuffle",
        difficulty=TIER_TRIVIAL,
        estimates=EST_CONSTRUCTION,
        description=(
            "Target-derived shuffled sequences with K/R residues held at their "
            "original positions. This is a construction diagnostic."
        ),
        caveat=(
            "Holding K/R positions does not preserve all tryptic boundaries when "
            "proline moves. Shuffling does not guarantee absence of shared target "
            "peptides or implement unique paired-peptide entrapment. Check overlap, "
            "effective partition sizes and uncertainty before interpreting a rate."
        ),
        is_definitional_null=False,
        recommended_headline=False,
    ),
    "thermophilic": EntrapmentClass(
        name="thermophilic",
        prefix="ENTRAP_PFUR_",
        generation="copy",
        difficulty=TIER_REMOTE,
        estimates=EST_CONSTRUCTION,
        description="A user-supplied Pyrococcus furiosus foreign-organism control.",
        caveat=(
            "Absence must be justified for the sample. Small or phylogenetically "
            "remote controls test limited failure modes. An empty effective partition"
            " is undefined, not zero FDP."
        ),
        requires_source=True,
    ),
    "archaea": EntrapmentClass(
        name="archaea",
        prefix="ENTRAP_ARCH_",
        generation="copy",
        difficulty=TIER_REMOTE,
        estimates=EST_CONSTRUCTION,
        description="A user-supplied thermophilic-archaeal sequence panel.",
        caveat=(
            "Control size and phylogenetic distance affect interpretation. Document "
            "the actual panel and establish absence at the peptide or protein "
            "counting level."
        ),
        requires_source=True,
    ),
    "plant": EntrapmentClass(
        name="plant",
        prefix="ENTRAP_PLANT_",
        generation="copy",
        difficulty=TIER_HOMOLOGOUS,
        estimates=EST_CONSTRUCTION,
        description="A user-supplied Arabidopsis thaliana sequence panel.",
        caveat=(
            "Dietary material and conserved peptides can violate control-absence "
            "assumptions. Inspect target overlap and biological context before "
            "interpreting apparent control matches as false discoveries."
        ),
        requires_source=True,
    ),
    "non_self": EntrapmentClass(
        name="non_self",
        prefix="ENTRAP_NONSELF_",
        generation="copy",
        difficulty=TIER_WITHIN_COHORT,
        estimates=EST_INCOMPLETENESS,
        description="Sequences from another acquisition's database within the declared cohort.",
        caveat=(
            "Shared communities can genuinely contain these sequences. A match cannot"
            " automatically be labelled a false discovery or a recall failure. Report"
            " cross-sample sharing as a separate diagnostic."
        ),
        requires_source=True,
        poolable=False,
    ),
    "user": EntrapmentClass(
        name="user",
        prefix="ENTRAP_USER_",
        generation="copy",
        difficulty=TIER_USER,
        estimates=EST_UNKNOWN,
        description="A user-supplied entrapment FASTA, tagged and used verbatim.",
        caveat=(
            "Interpretation is the user's responsibility. The entrapment axiom "
            "requires that no sequence in this file is genuinely present in the "
            "sample; if that fails, hits are true identifications and the FDP is "
            "meaningless. Run the conserved-peptide report to check."
        ),
        requires_source=True,
    ),
}


def get_class(name: str) -> EntrapmentClass:
    """Look up an entrapment class by name."""
    try:
        return ENTRAPMENT_CLASSES[name]
    except KeyError:
        raise ValueError(
            f"Unknown entrapment class {name!r}. Available: {', '.join(sorted(ENTRAPMENT_CLASSES))}"
        ) from None


def class_names() -> list[str]:
    """Registered class names, preserving the historical display order."""
    names = sorted(ENTRAPMENT_CLASSES)
    names.remove("shuffled")
    return ["shuffled", *names]
