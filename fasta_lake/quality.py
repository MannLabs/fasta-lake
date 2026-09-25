"""
Protein sequence quality filtering.

Optional filters that remove low-quality sequences BEFORE peptide matching.
Each filter is independently toggleable via ``SequenceQualityConfig``.

These filters target common database artifacts:

- **Poly-runs**: Assembly artifacts producing ``AAAAAAAAAA...`` sequences
- **Low-complexity**: Sequences dominated by 1-2 amino acids (low entropy)
- **Fragments**: Truncated ORFs lacking start methionine

Design Principles
-----------------
- No hidden transformations: every rejection is logged with reason
- Ambiguous residues may be replaced with X when explicitly configured
- Defaults and individual filter switches are defined by SequenceQualityConfig
"""

from __future__ import annotations

import logging
import math
from collections import Counter
from dataclasses import dataclass, field

from fasta_lake.config import SequenceQualityConfig

logger = logging.getLogger(__name__)

# Genetically encoded residues, including Sage-supported Sec (U) and Pyl (O).
AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYUO")

# Non-canonical residues commonly found in databases:
#   B = Asx (D or N)      — disambiguated during translation
#   J = Xle (I or L)      — ambiguity
#   Z = Glx (E or Q)      — disambiguated during translation
#   X = unknown residue   — kept (represents any residue)
#   * = stop codon symbol — stripped (terminator)
NONCANONICAL_REPLACE = {"B": "X", "J": "X", "Z": "X"}
STOP_SYMBOL = "*"


def normalise_residues(sequence: str) -> tuple[str, int]:
    """
    Clean a protein sequence: strip trailing stop symbols, replace non-canonical
    ambiguous residues (B/J/Z) with X. Preserve selenocysteine and pyrrolysine.

    Returns
    -------
    cleaned : str
        Normalised uppercase sequence.
    n_replaced : int
        Count of ambiguous residues and internal stop symbols replaced.
        Trailing stop symbols are stripped and are not counted.
    """
    seq = sequence.upper().rstrip(STOP_SYMBOL)  # drop trailing stop
    n_replaced = sum(seq.count(c) for c in NONCANONICAL_REPLACE)
    # Also count internal stops if any
    n_replaced += seq.count(STOP_SYMBOL)
    seq = seq.replace(STOP_SYMBOL, "X")
    for orig, repl in NONCANONICAL_REPLACE.items():
        if orig in seq:
            seq = seq.replace(orig, repl)
    return seq, n_replaced


def noncanonical_fraction(sequence: str) -> float:
    """Fraction of residues outside the 22 genetically encoded amino acids."""
    if not sequence:
        return 0.0
    non = sum(1 for c in sequence.upper() if c not in AA_ALPHABET)
    return non / len(sequence)


def has_poly_run(sequence: str, run_length: int = 10) -> bool:
    """
    Check if a sequence contains a homopolymeric run.

    A homopolymeric run is ``run_length`` or more consecutive identical
    amino acids (e.g., ``AAAAAAAAAA``). These are typically assembly
    artifacts or low-quality ORF predictions.

    Parameters
    ----------
    sequence : str
        Protein sequence (uppercase).
    run_length : int
        Minimum consecutive identical residues to trigger. Default: 10.

    Returns
    -------
    bool
        True if the sequence contains a poly-run.

    Examples
    --------
    >>> has_poly_run("MSKAAAAAAAAAAAGELFT", run_length=10)
    True
    >>> has_poly_run("MSKGEELFT", run_length=10)
    False
    """
    if run_length < 2 or len(sequence) < run_length:
        return False

    count = 1
    for i in range(1, len(sequence)):
        if sequence[i] == sequence[i - 1]:
            count += 1
            if count >= run_length:
                return True
        else:
            count = 1

    return False


def compute_sequence_entropy(sequence: str) -> float:
    """
    Compute normalized Shannon entropy of a protein sequence.

    Measures amino acid diversity on a 0-1 scale:

    - 0.0: completely uniform (e.g., ``AAAAAAA``)
    - 1.0: maximum diversity (all 20 amino acids equally represented)

    Uses log base 20 (number of canonical amino acids) for normalization,
    so the twenty-residue alphabet has maximum 1. Values can exceed 1
    when additional residue symbols occur; the filtering scale is unchanged.

    Parameters
    ----------
    sequence : str
        Protein sequence (uppercase).

    Returns
    -------
    float
        Shannon entropy divided by log2(20).

    Examples
    --------
    >>> compute_sequence_entropy("AAAAAAAAAA")
    0.0
    >>> compute_sequence_entropy("ACDEFGHIKLMNPQRSTVWY")  # all 20 AAs
    1.0
    """
    if len(sequence) < 2:
        return 0.0

    counts = Counter(sequence)
    n = len(sequence)
    entropy = 0.0

    for count in counts.values():
        if count > 0:
            p = count / n
            entropy -= p * math.log2(p)

    # Normalize by maximum possible entropy (log2(20) for 20 amino acids)
    max_entropy = math.log2(20)
    if max_entropy == 0:
        return 0.0

    return entropy / max_entropy


def is_fragment(
    sequence: str,
    length_threshold: int = 100,
) -> bool:
    """
    Check if a protein is likely a fragment.

    A fragment is defined as a sequence that:
    1. Does NOT start with methionine (M)
    2. Is shorter than ``length_threshold`` amino acids

    Long sequences without M-start are less likely to be fragments
    (could be legitimate mature proteins with signal peptide removed).

    Parameters
    ----------
    sequence : str
        Protein sequence (uppercase).
    length_threshold : int
        Maximum length for fragment classification. Sequences shorter
        than this AND lacking M-start are flagged. Default: 100.

    Returns
    -------
    bool
        True if the sequence is likely a fragment.

    Examples
    --------
    >>> is_fragment("SKGEELFT", length_threshold=100)
    True
    >>> is_fragment("MSKGEELFT", length_threshold=100)
    False
    >>> is_fragment("SKGEELFT" * 20, length_threshold=100)  # 160 aa
    False
    """
    if not sequence:
        return True
    return not sequence.startswith("M") and len(sequence) < length_threshold


@dataclass
class QualityReport:
    """
    Report from sequence quality filtering.

    Tracks how many sequences were removed by each filter and why.

    Parameters
    ----------
    total_input : int
        Total sequences before filtering.
    passed : int
        Sequences that passed all filters.
    rejected_length : int
        Rejected for length (too short or too long).
    rejected_poly_run : int
        Rejected for homopolymeric runs.
    rejected_low_complexity : int
        Rejected for low sequence entropy.
    rejected_fragment : int
        Rejected as likely fragments.
    rejection_details : dict[str, str]
        Mapping of rejected protein_id -> rejection reason.
    """

    total_input: int = 0
    passed: int = 0
    rejected_length: int = 0
    rejected_poly_run: int = 0
    rejected_low_complexity: int = 0
    rejected_fragment: int = 0
    rejected_noncanonical: int = 0  # too many ambiguous/unknown residues
    residues_normalised: int = 0  # B/J/Z → X replacements made
    rejection_details: dict[str, str] = field(default_factory=dict)

    def summary(self) -> dict:
        """Return summary dict for JSON stats output."""
        total_rejected = (
            self.rejected_length
            + self.rejected_poly_run
            + self.rejected_low_complexity
            + self.rejected_fragment
            + self.rejected_noncanonical
        )
        return {
            "total_input": self.total_input,
            "passed": self.passed,
            "total_rejected": total_rejected,
            "rejected_length": self.rejected_length,
            "rejected_poly_run": self.rejected_poly_run,
            "rejected_low_complexity": self.rejected_low_complexity,
            "rejected_fragment": self.rejected_fragment,
            "rejected_noncanonical": self.rejected_noncanonical,
            "residues_normalised": self.residues_normalised,
        }


def filter_sequences_by_quality(
    proteins: dict[str, str],
    config: SequenceQualityConfig | None = None,
) -> tuple[dict[str, str], QualityReport]:
    """
    Apply sequence quality filters to a protein dictionary.

    Filters are applied in order: length -> poly-run -> low-complexity ->
    fragment. A sequence rejected by an earlier filter is not checked by
    later filters (first rejection wins).

    Parameters
    ----------
    proteins : dict[str, str]
        Mapping of protein ID -> sequence.
    config : SequenceQualityConfig or None
        Quality filter configuration. If None, uses defaults (no
        optional filters enabled, only length filtering).

    Returns
    -------
    filtered : dict[str, str]
        Proteins that passed all filters.
    report : QualityReport
        Detailed rejection report.

    Examples
    --------
    >>> from fasta_lake.config import SequenceQualityConfig
    >>> cfg = SequenceQualityConfig(reject_poly_runs=True, poly_run_length=5)
    >>> proteins = {"p1": "MSKGEELFT", "p2": "MAAAAAAAAAGEELFT"}
    >>> filtered, report = filter_sequences_by_quality(proteins, cfg)
    >>> "p1" in filtered
    True
    >>> "p2" in filtered
    False
    """
    if config is None:
        config = SequenceQualityConfig()

    report = QualityReport(total_input=len(proteins))
    filtered: dict[str, str] = {}

    for protein_id, sequence in proteins.items():
        seq_upper = sequence.upper()

        # Non-canonical fraction check on ORIGINAL (before stripping), so that
        # proteins with many ambiguous residues get rejected regardless of
        # whether strip_noncanonical is on.
        if config.reject_if_noncanonical_fraction > 0:
            frac = noncanonical_fraction(seq_upper)
            if frac > config.reject_if_noncanonical_fraction:
                report.rejected_noncanonical += 1
                report.rejection_details[protein_id] = (
                    f"noncanonical:frac={frac:.2%}>{config.reject_if_noncanonical_fraction:.2%}"
                )
                continue

        # Ambiguous residue normalisation (B/J/Z → X); retain encoded U/O.
        if config.strip_noncanonical:
            seq_upper, n_replaced = normalise_residues(seq_upper)
            report.residues_normalised += n_replaced

        # Length filter (always applied)
        if len(seq_upper) < config.min_length:
            report.rejected_length += 1
            report.rejection_details[protein_id] = f"too_short:{len(seq_upper)}<{config.min_length}"
            continue
        if len(seq_upper) > config.max_length:
            report.rejected_length += 1
            report.rejection_details[protein_id] = f"too_long:{len(seq_upper)}>{config.max_length}"
            continue

        # Poly-run filter (optional)
        if config.reject_poly_runs and has_poly_run(seq_upper, config.poly_run_length):
            report.rejected_poly_run += 1
            report.rejection_details[protein_id] = f"poly_run:>={config.poly_run_length}"
            continue

        # Low-complexity filter (optional)
        if config.reject_low_complexity:
            entropy = compute_sequence_entropy(seq_upper)
            if entropy < config.low_complexity_threshold:
                report.rejected_low_complexity += 1
                report.rejection_details[protein_id] = (
                    f"low_complexity:entropy={entropy:.3f}<{config.low_complexity_threshold}"
                )
                continue

        # Fragment filter (optional)
        if config.reject_fragments and is_fragment(seq_upper, config.fragment_length_threshold):
            report.rejected_fragment += 1
            report.rejection_details[protein_id] = (
                f"fragment:no_M_start,len={len(seq_upper)}<{config.fragment_length_threshold}"
            )
            continue

        # Use the cleaned sequence, not the original, so downstream sees normalised form
        filtered[protein_id] = seq_upper

    report.passed = len(filtered)

    logger.info(
        "Quality filtering: %d/%d passed (-%d length, -%d poly, -%d complexity, -%d fragment)",
        report.passed,
        report.total_input,
        report.rejected_length,
        report.rejected_poly_run,
        report.rejected_low_complexity,
        report.rejected_fragment,
    )

    return filtered, report
