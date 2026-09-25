"""
Deterministic sequence hashing for FastaLake (v5.3).

One frozen contract used by every stage that identifies a protein by its
amino-acid sequence. Python (`fasta_lake.lake_merge`, inference, lineage
aggregation) and Rust (`rust/fasta_extractor_v2/src/bin/lake_builder.rs`) both
implement this contract; `tests/test_hashing_crosslang.py` asserts byte-exact
agreement.

Normalization contract
----------------------
Given an input string, `normalize_sequence` yields the bytes actually fed to
SHA256:

1. Keep only ASCII letters [A-Z a-z]. Everything else is dropped — whitespace,
   digits, punctuation, stop codons (`*`), newlines, non-ASCII characters.
2. Lowercase letters are folded to uppercase.
3. Ambiguity codes (B, Z, X, U, O) are preserved verbatim because they are
   letters. Two sequences differing only in ambiguity-code placement are
   considered distinct.

Rationale
---------
- Matches the Rust lake_builder's `extract_clean_sequence` behavior. Any
  divergence between stages is a silent dedup bug, so the contract is
  deliberately permissive about input noise (whitespace, stop codons) and
  strict about output bytes (A-Z only, uppercase).
- Stop codons (`*`) represent translation artifacts, not biology. Dropping
  them means `ARG*` and `ARG` dedup to the same sequence. Internal `*` (a
  mid-sequence stop, typically a data-quality issue in a bad translation) is
  also stripped. Lake-build stats log the count of sequences containing `*`.
- Ambiguity codes are kept literal to avoid guessing biology. Expanding them
  (X → any, B → D/N, Z → E/Q) would invent identity relationships that do not
  exist in the source database.
- SHA256 is chosen for reviewer-defensibility over MD5 (cryptographically
  broken) and non-cryptographic hashes (xxhash, FNV — adversarial collisions
  trivial). Hashing is I/O-bound at our scale, so the ~2x cost vs MD5 is not
  observable in end-to-end timings.

Versioning
----------
A lake built with the v5.3 hash function is not compatible with pre-v5.3
lakes (which used MD5). The `hashes.version` field in lake stats records
this; rebuilding is the only migration path.
"""

from __future__ import annotations

import hashlib

__all__ = [
    "HASH_VERSION",
    "HASH_ALGORITHM",
    "normalize_sequence",
    "sequence_hash",
]

# Bumped whenever the normalization contract or hash algorithm changes.
# Consumers that persist hashes should record this alongside them.
HASH_VERSION = "5.3"
HASH_ALGORITHM = "sha256"


def normalize_sequence(seq: str) -> bytes:
    """Normalize an amino-acid sequence per the frozen hashing contract.

    Returns the exact bytes that `sequence_hash` feeds to SHA256, so callers
    that need to log or round-trip the canonical form can do so.

    See module docstring for the contract. This implementation must stay
    byte-identical to Rust's `normalize_sequence` in lake_builder.rs —
    `tests/test_hashing_crosslang.py` enforces that.
    """
    out = bytearray()
    for ch in seq:
        c = ord(ch)
        if 0x41 <= c <= 0x5A:  # A-Z
            out.append(c)
        elif 0x61 <= c <= 0x7A:  # a-z → uppercase
            out.append(c - 0x20)
    return bytes(out)


def sequence_hash(seq: str) -> str:
    """Return the v5.3 SHA256 hex digest of a normalized amino-acid sequence.

    64 lowercase hex characters. Deterministic across Python/Rust.
    """
    return hashlib.sha256(normalize_sequence(seq)).hexdigest()
