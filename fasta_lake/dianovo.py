"""Read complete DIANovo predictions without guessing across unresolved mass gaps.

DIANovo: https://github.com/hearthewind/dianovo (cc413f25140e47a5f6a4234bae15da0ef525ac3a),
``denovo_from_feature_detection/rnova/task/sequence_generation_inference_utils.py``
and ``optimal_path_inference_classic.py``. These are FastaLake input checks, not
DIANovo inference. Sequence generation writes one probability per residue; the
space-delimited path output prepends a zero marker. See docs/how-to/de-novo-inputs.md.
The matching Rust implementation is rust/common/dianovo.rs.
"""

import math

_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def parse_dianovo_prediction(sequence: str, probabilities: str) -> tuple[str, float]:
    """Return a complete uppercase peptide and mean confidence, including zeros.

    Accept lowercase residues and whitespace between residues. Reject mass gaps,
    noncanonical residues, malformed lists, nonfinite/out-of-range probabilities
    and length mismatches with ValueError. Only a space-delimited vector with
    exactly one extra leading zero is interpreted as a path-marker vector.
    All-zero residue vectors are valid here; positive-score selection excludes them.
    """
    raw_sequence = "".join(c for c in sequence if c not in " \t\r\n\v\f")
    if not raw_sequence.isascii():
        raise ValueError("unresolved_gap_or_noncanonical_sequence")
    peptide = raw_sequence.upper()
    if not peptide or any(aa not in _AMINO_ACIDS for aa in peptide):
        raise ValueError("unresolved_gap_or_noncanonical_sequence")
    raw = probabilities.strip()
    bracketed = raw.startswith("[") and raw.endswith("]")
    if bracketed:
        raw = raw[1:-1].strip()
    comma_format = bracketed or "," in raw
    tokens = raw.split(",") if comma_format else raw.split()
    if any(not token.isascii() or "_" in token for token in tokens):
        raise ValueError("invalid_probability_vector")
    try:
        values = [float(token.strip()) for token in tokens]
    except ValueError as exc:
        raise ValueError("invalid_probability_vector") from exc
    if not values or any(not math.isfinite(p) or not 0 <= p <= 1 for p in values):
        raise ValueError("invalid_probability_vector")
    if not comma_format and len(values) == len(peptide) + 1 and values[0] == 0:
        values = values[1:]
    if len(values) != len(peptide):
        raise ValueError("sequence_probability_length_mismatch")
    # Explicit left-to-right f64 addition matches Rust's ranking at cutoff ties.
    total = 0.0
    for value in values:
        total += value
    return peptide, total / len(values)
