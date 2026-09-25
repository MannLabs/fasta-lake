"""
BLOSUM-guided amino acid substitution scoring.

BLOSUM95 is optimized for closely related sequences (>95% identity),
appropriate for single amino acid variant (SAV) discovery within a species.

For de novo error correction, we also support:
- BLOSUM62: general-purpose, good for cross-species comparisons
- De novo error matrix: empirical substitution rates from de novo prediction
  errors (mass-equivalent substitutions, terminal errors, I/L confusion)

The de novo error matrix captures what de novo sequencing engines actually
get wrong — which is different from evolutionary substitution patterns.
Common de novo errors:
  - I/L confusion (isobaric, 100% confusion rate)
  - N/D confusion (deamidation, +0.98 Da)
  - Q/E confusion (deamidation, +0.98 Da)
  - K/Q confusion (near-isobaric, 0.036 Da difference)
  - Adjacent AA swaps (mass-preserving)
  - F/oxidized-M (near-isobaric)

Reference:
  - Henikoff & Henikoff (1992) BLOSUM matrices
  - Miura et al. (2025) ICS-guided refinement
  - Lebedev, Treit, Mann (2026) de novo error taxonomy
"""

from __future__ import annotations

from typing import Dict, List, Tuple

# Standard amino acids
AMINO_ACIDS = "ARNDCQEGHILKMFPSTWYV"

# ============================================================================
# BLOSUM95
# ============================================================================
_BLOSUM95_RAW = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
A  5 -2 -1 -2 -1 -1 -1  0 -2 -1 -2 -1 -1 -3 -1  1  0 -3 -2  0
R -2  6 -1 -2 -4  1  0 -3  0 -3 -3  2 -2 -4 -2 -1 -1 -3 -2 -3
N -1 -1  6  1 -3  0 -1 -1  0 -4 -4  0 -3 -4 -3  0  0 -4 -3 -4
D -2 -2  1  6 -4 -1  1 -2 -2 -4 -5 -1 -4 -4 -2 -1 -1 -5 -4 -4
C -1 -4 -3 -4  9 -4 -5 -4 -4 -2 -2 -4 -2 -3 -4 -2 -1 -3 -3 -1
Q -1  1  0 -1 -4  6  2 -2  1 -3 -3  1  0 -4 -2  0 -1 -3 -2 -3
E -1  0 -1  1 -5  2  6 -3  0 -4 -4  1 -2 -4 -2  0 -1 -4 -3 -3
G  0 -3 -1 -2 -4 -2 -3  6 -3 -5 -4 -2 -4 -4 -3 -1 -2 -4 -4 -4
H -2  0  0 -2 -4  1  0 -3  8 -4 -3 -1 -2 -2 -3 -1 -2 -3  1 -4
I -1 -3 -4 -4 -2 -3 -4 -5 -4  5  1 -3  1 -1 -4 -3 -1 -3 -2  3
L -2 -3 -4 -5 -2 -3 -4 -4 -3  1  4 -3  2  0 -4 -3 -2 -2 -2  1
K -1  2  0 -1 -4  1  1 -2 -1 -3 -3  5 -2 -4 -1 -1 -1 -4 -3 -3
M -1 -2 -3 -4 -2  0 -2 -4 -2  1  2 -2  6  0 -3 -2 -1 -2 -2  1
F -3 -4 -4 -4 -3 -4 -4 -4 -2 -1  0 -4  0  6 -4 -3 -2  0  3 -1
P -1 -2 -3 -2 -4 -2 -2 -3 -3 -4 -4 -1 -3 -4  8 -1 -2 -4 -4 -3
S  1 -1  0 -1 -2  0  0 -1 -1 -3 -3 -1 -2 -3 -1  5  1 -4 -2 -2
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -2 -1 -1 -2 -2  1  5 -4 -2  0
W -3 -3 -4 -5 -3 -3 -4 -4 -3 -3 -2 -4 -2  0 -4 -4 -4 11  2 -3
Y -2 -2 -3 -4 -3 -2 -3 -4  1 -2 -2 -3 -2  3 -4 -2 -2  2  7 -2
V  0 -3 -4 -4 -1 -3 -3 -4 -4  3  1 -3  1 -1 -3 -2  0 -3 -2  5
"""

# ============================================================================
# De novo error matrix — empirical confusion rates from mass spec de novo
# Scores represent log-likelihood of confusion (higher = more likely to confuse)
# Based on error taxonomy from Lebedev, Treit, Mann (2026)
# ============================================================================
_DENOVO_CONFUSIONS: Dict[Tuple[str, str], int] = {
    # Isobaric (always confused)
    ("I", "L"): 10,
    ("L", "I"): 10,
    # Deamidation-like (N->D, Q->E: +0.98 Da)
    ("N", "D"): 6,
    ("D", "N"): 6,
    ("Q", "E"): 6,
    ("E", "Q"): 6,
    # Near-isobaric (K/Q: 0.036 Da, mass-equivalent at low resolution)
    ("K", "Q"): 4,
    ("Q", "K"): 4,
    # Oxidized-M ~ F (147.035 vs 147.068, 0.033 Da)
    ("M", "F"): 3,
    ("F", "M"): 3,
    # Common adjacent swaps (mass-preserving, frequent de novo error)
    # These don't change the total mass but reorder residues
    # Represented as high self-scores for flanking AAs
    # GG/AG/GA type swaps
    ("A", "G"): 2,
    ("G", "A"): 2,
    ("A", "S"): 1,
    ("S", "A"): 1,
    ("V", "A"): 1,
    ("A", "V"): 1,
    # W ~ DA, W ~ AD (tryptophan = Asp + Ala mass)
    ("W", "D"): 2,
    ("D", "W"): 2,
    ("W", "A"): 2,
    ("A", "W"): 2,
}


def _parse_matrix(raw: str) -> Dict[Tuple[str, str], int]:
    """Parse a whitespace-delimited substitution matrix into residue-pair scores."""
    lines = [ln.strip() for ln in raw.strip().split("\n") if ln.strip()]
    header = lines[0].split()
    m: Dict[Tuple[str, str], int] = {}
    for line in lines[1:]:
        parts = line.split()
        aa1 = parts[0]
        for i, aa2 in enumerate(header):
            m[(aa1, aa2)] = int(parts[i + 1])
    return m


BLOSUM95 = _parse_matrix(_BLOSUM95_RAW)


def score(aa1: str, aa2: str, matrix: str = "blosum95") -> int:
    """
    Substitution score for aa1 -> aa2.

    Parameters
    ----------
    matrix : str
        ``'blosum95'`` or ``'denovo'``.
    """
    if matrix == "denovo":
        return _DENOVO_CONFUSIONS.get((aa1.upper(), aa2.upper()), -5)
    return BLOSUM95.get((aa1.upper(), aa2.upper()), -4)


def substitutions(
    aa: str,
    min_score: int = 0,
    matrix: str = "blosum95",
) -> List[Tuple[str, int]]:
    """
    Possible substitutions for an amino acid, sorted by score descending.

    Parameters
    ----------
    aa : str
        Source amino acid.
    min_score : int
        Minimum score to include.
    matrix : str
        ``'blosum95'`` or ``'denovo'``.

    Returns
    -------
    list of (target_aa, score)
    """
    aa = aa.upper()
    if aa not in AMINO_ACIDS:
        return []

    m = _DENOVO_CONFUSIONS if matrix == "denovo" else BLOSUM95
    subs = []
    for target in AMINO_ACIDS:
        if target != aa:
            s = m.get((aa, target), -5 if matrix == "denovo" else -4)
            if s >= min_score:
                subs.append((target, s))
    subs.sort(key=lambda x: -x[1])
    return subs


def generate_variants(
    sequence: str,
    position: int,
    min_score: int = -1,
    matrix: str = "blosum95",
    max_per_position: int = 3,
) -> List[Tuple[str, str, int]]:
    """
    Generate single-AA variants at a position.

    Returns
    -------
    list of (variant_sequence, mutation_notation, score)
    """
    if position < 0 or position >= len(sequence):
        return []
    aa = sequence[position].upper()
    subs = substitutions(aa, min_score=min_score, matrix=matrix)
    variants = []
    for target, s in subs[:max_per_position]:
        var_seq = sequence[:position] + target + sequence[position + 1 :]
        notation = f"{aa}{position + 1}{target}"
        variants.append((var_seq, notation, s))
    return variants


def generate_all_variants(
    sequence: str,
    min_score: int = -1,
    matrix: str = "blosum95",
    max_per_position: int = 3,
) -> List[Tuple[str, str, int]]:
    """
    Generate all single-substitution variants for a sequence.

    Parameters
    ----------
    matrix : str
        ``'blosum95'`` for evolutionary variants, ``'denovo'`` for
        de novo error correction.
    """
    all_vars = []
    for pos in range(len(sequence)):
        all_vars.extend(generate_variants(sequence, pos, min_score, matrix, max_per_position))
    return all_vars
