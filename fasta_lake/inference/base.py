"""
Data structures for the inference engine.

Defines the input/output contracts for all inference strategies.
Each strategy function takes an ``InferenceInput`` and returns
an ``InferenceResult``.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class InferenceInput:
    """
    Pre-computed data passed to every inference strategy.

    Parameters
    ----------
    pep2prot : dict[str, set[str]]
        Peptide → set of protein IDs.
    prot2pep : dict[str, set[str]]
        Protein ID → set of peptides.
    proteins : dict[str, str]
        Protein ID → sequence (from stage 2 FASTA).
    peptide_scores : dict[str, float]
        Peptide → best de novo confidence score. Used for score-weighted
        tie-breaking in inference strategies (higher score = stronger evidence).
        If empty, strategies fall back to count-based tie-breaking.
    """

    pep2prot: dict[str, set[str]]
    prot2pep: dict[str, set[str]]
    proteins: dict[str, str]
    peptide_scores: dict[str, float] = field(default_factory=dict)
    taxonomy_map: dict[str, str] = field(default_factory=dict)
    taxonomy_pruning: bool = False


@dataclass
class InferenceResult:
    """
    Output from an inference strategy.

    Parameters
    ----------
    selected_proteins : set[str]
        Set of selected protein IDs.
    stats : dict
        Strategy-specific statistics.
    """

    selected_proteins: set[str] = field(default_factory=set)
    stats: dict = field(default_factory=dict)
