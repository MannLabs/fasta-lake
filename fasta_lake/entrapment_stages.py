"""
Entrapment injection-stage registry.

Where entrapment enters the pipeline decides what its FDP measures. This module
is the single source of truth for the stage vocabulary, and it deliberately uses
the names the run scripts and figure loaders already use so nothing downstream
breaks.

The canonical chain
-------------------
::

    stage 0  lake build
    stage 1  cohort evidence filter   (de novo peptides admit sequences to the lake)
    stage 1.5 cluster97               (MMseqs2 --min-seq-id 0.97 -c 0.8 --cov-mode 0)
    stage 2  per-sample extraction    (de novo string matching = "selection")
    stage 3  parsimony / razor
    stage 4  SAGE search

Four injection points exist in this repository, each verified against the run
script that implements it. They are ordered here from earliest to latest, which is
also from hardest to easiest for an entrapment sequence to survive.

``prelake``
    Entrapment is concatenated into the **input lake, before the cohort evidence
    filter**. It must then survive evidence filtering, clustering, per-sample
    selection, parsimony and the search. This is the fullest construction test --
    entrapment faces every curation step a real sequence faces.
    Implemented by ``projects/nearneighbour_entrapment_2026-07-20/scripts/
    00b_build_lakes.sh`` (``cat TARGET_FASTA TAGGED > lake_NN.fasta``), consumed by
    ``01_evidence_cluster.sh``.

``prebuild``
    Entrapment is concatenated onto the **already-clustered lake and the result is
    re-clustered at 97%**, so entrapment must also survive clustering -- it can be
    absorbed into a target cluster, or absorb targets into itself. It does **not**
    face the cohort evidence filter.
    Implemented by ``projects/masters_circularity_F1/scripts/inject_prep.sh:54-56``
    and ``projects/entrapment_rerun_2026-07-09/scripts/01_prebuild_prep.sh``
    (``mmseqs easy-cluster`` over the augmented lake).

``preselect``
    Entrapment is concatenated onto the **already-clustered lake with no
    re-clustering**. It faces per-sample selection, parsimony and the search.
    Cluster97 is *upstream* of this injection point, not downstream.
    Implemented by ``masters_circularity_F1/scripts/f1b_prep.sh:25-26`` and the
    ``preselect`` arm of ``inject_prep.sh``.

``appended``
    Entrapment is concatenated onto the **finished per-sample database after
    parsimony**, so it enters the search untested by any curation step. Its FDP
    measures ordinary search/scoring error, **not** construction error.
    Implemented by ``inject_persample.sh:72-74`` and ``03_appended_sage.sh``.

Why the distinction matters
---------------------------
``prebuild`` and ``preselect`` are **not** the same stage under two names: they
differ by whether entrapment had to survive clustering. Collapsing them loses the
strongest construction claim available. ``prelake`` is strictly earlier than both.

Numbers are deliberately absent from this module. Stage metadata describes the
design; measured FDPs belong to a run's output record, not to the schema.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

#: Curation steps an entrapment sequence can be required to survive, in order.
STEP_EVIDENCE_FILTER = "cohort evidence filter"
STEP_CLUSTER97 = "cluster97"
STEP_SELECTION = "per-sample de novo selection"
STEP_PARSIMONY = "parsimony/razor"
STEP_SEARCH = "search (SAGE)"

MEASURES_CONSTRUCTION = "database-construction error"
MEASURES_SEARCH = "search error (NOT construction error)"


@dataclass(frozen=True)
class EntrapmentStage:
    """One injection point and what its FDP therefore estimates."""

    name: str
    order: int
    injection_point: str
    faces: tuple[str, ...]
    measures: str
    is_construction_test: bool
    description: str
    reference_script: str

    def as_metadata(self) -> dict:
        """Serialisable metadata for embedding in a manifest or FDP record."""
        d = asdict(self)
        d["faces"] = list(self.faces)
        return d


ENTRAPMENT_STAGES: dict[str, EntrapmentStage] = {
    "prelake": EntrapmentStage(
        name="prelake",
        order=0,
        injection_point="input lake, before the cohort evidence filter",
        faces=(STEP_EVIDENCE_FILTER, STEP_CLUSTER97, STEP_SELECTION, STEP_PARSIMONY, STEP_SEARCH),
        measures=MEASURES_CONSTRUCTION,
        is_construction_test=True,
        description=(
            "Fullest construction test: entrapment faces every curation step a "
            "real sequence faces, starting with the cohort evidence filter."
        ),
        reference_script=(
            "projects/nearneighbour_entrapment_2026-07-20/scripts/00b_build_lakes.sh"
        ),
    ),
    "prebuild": EntrapmentStage(
        name="prebuild",
        order=1,
        injection_point="clustered lake + entrapment, then re-clustered at 97%",
        faces=(STEP_CLUSTER97, STEP_SELECTION, STEP_PARSIMONY, STEP_SEARCH),
        measures=MEASURES_CONSTRUCTION,
        is_construction_test=True,
        description=(
            "Entrapment must also survive clustering: it can be absorbed into a "
            "target cluster or absorb targets. Does NOT face the cohort evidence "
            "filter, which ran before injection."
        ),
        reference_script="projects/masters_circularity_F1/scripts/inject_prep.sh",
    ),
    "preselect": EntrapmentStage(
        name="preselect",
        order=2,
        injection_point="clustered lake + entrapment, no re-clustering",
        faces=(STEP_SELECTION, STEP_PARSIMONY, STEP_SEARCH),
        measures=MEASURES_CONSTRUCTION,
        is_construction_test=True,
        description=(
            "Entrapment must be nominated by the sample's own de novo peptides. "
            "cluster97 is UPSTREAM of this injection point, so entrapment skips "
            "both the evidence filter and clustering."
        ),
        reference_script="projects/masters_circularity_F1/scripts/f1b_prep.sh",
    ),
    "appended": EntrapmentStage(
        name="appended",
        order=3,
        injection_point="finished per-sample database, after parsimony",
        faces=(STEP_SEARCH,),
        measures=MEASURES_SEARCH,
        is_construction_test=False,
        description=(
            "Entrapment enters the search untested by any curation step, so a hit "
            "means the search engine wrongly scored it. Legitimate, but it is a "
            "different quantity from the construction stages and must never be "
            "pooled with or compared against them."
        ),
        reference_script=("projects/masters_circularity_F1/scripts/inject_persample.sh"),
    ),
}

#: Historical / alternative spellings accepted on input and normalised away.
#: ``pre-selection`` was used by an earlier version of the CLI and is the label
#: ``manuscript/figures/fig5_versions/fig5_build.py`` prints for ``preselect``.
STAGE_ALIASES: dict[str, str] = {
    "pre-selection": "preselect",
    "pre_selection": "preselect",
    "pre-build": "prebuild",
    "pre_build": "prebuild",
    "pre-lake": "prelake",
    "pre_lake": "prelake",
}

#: Canonical stage names, earliest injection point first.
VALID_ENTRAPMENT_STAGES: tuple[str, ...] = tuple(
    s.name for s in sorted(ENTRAPMENT_STAGES.values(), key=lambda s: s.order)
)

#: Stages whose FDP is a database-construction measurement.
CONSTRUCTION_STAGES: tuple[str, ...] = tuple(
    s.name
    for s in sorted(ENTRAPMENT_STAGES.values(), key=lambda s: s.order)
    if s.is_construction_test
)


def normalize_stage(stage: str) -> str:
    """Map an accepted spelling onto a canonical stage name.

    Raises ``ValueError`` on anything unrecognised -- the stage determines what a
    number means and is never guessed.
    """
    if stage in ENTRAPMENT_STAGES:
        return stage
    if stage in STAGE_ALIASES:
        return STAGE_ALIASES[stage]
    raise ValueError(
        f"Unknown entrapment stage {stage!r}. Valid stages, earliest injection "
        f"point first: {', '.join(VALID_ENTRAPMENT_STAGES)}"
    )


def get_stage(stage: str) -> EntrapmentStage:
    """Look up a stage, accepting aliases."""
    return ENTRAPMENT_STAGES[normalize_stage(stage)]


def stages_comparable(a: str, b: str) -> bool:
    """True only when two stages are the same injection point.

    Two construction stages are still not interchangeable: ``prebuild`` entrapment
    had to survive clustering and ``preselect`` entrapment did not, so their FDPs
    are not the same quantity.
    """
    return normalize_stage(a) == normalize_stage(b)
