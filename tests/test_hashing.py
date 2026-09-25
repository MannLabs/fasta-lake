"""
Unit tests for the frozen sequence-hashing contract.

These tests pin the v5.3 behavior. If any of them change meaning, the
normalization contract has changed, and that requires:
  - Bumping HASH_VERSION in fasta_lake/hashing.py
  - Updating the Rust counterpart in lake_builder.rs
  - Rebuilding any persisted lake (hashes are not backward-compatible)
  - Re-running tests/test_hashing_crosslang.py

The cross-language equality test lives in test_hashing_crosslang.py.
"""

from __future__ import annotations

import hashlib

import pytest

from fasta_lake.hashing import (
    HASH_ALGORITHM,
    HASH_VERSION,
    normalize_sequence,
    sequence_hash,
)

# ---------------------------------------------------------------------------
# Contract metadata — these surface in lake stats; changing them is a breaking
# change for persisted lakes.
# ---------------------------------------------------------------------------


def test_hash_algorithm_is_sha256():
    assert HASH_ALGORITHM == "sha256"


def test_hash_version_v53():
    assert HASH_VERSION == "5.3"


def test_sha256_digest_length():
    # SHA256 → 32 bytes → 64 hex chars. Distinguishes from MD5 (32 hex).
    assert len(sequence_hash("MAAVR")) == 64


# ---------------------------------------------------------------------------
# Normalization rules — the six edge cases the contract commits to.
# ---------------------------------------------------------------------------


class TestNormalizationRules:
    def test_trailing_stop_codon_stripped(self):
        # ARG* and ARG must hash identically (stop codon is a translation artifact).
        assert sequence_hash("ARG*") == sequence_hash("ARG")

    def test_internal_stop_codon_stripped(self):
        # Internal * (mid-sequence stop from bad translation / frame-shift) is
        # stripped — matches Rust's behavior. A sequence of interleaved valid
        # residues and *s hashes to the concatenation of the valid residues.
        assert sequence_hash("AR*G") == sequence_hash("ARG")
        assert sequence_hash("A*R*G*") == sequence_hash("ARG")

    def test_whitespace_stripped(self):
        # Spaces, tabs, and newlines (common in multi-line FASTA lines) must
        # not change the hash.
        assert sequence_hash("AR G") == sequence_hash("ARG")
        assert sequence_hash("AR\tG") == sequence_hash("ARG")
        assert sequence_hash("AR\nG") == sequence_hash("ARG")
        assert sequence_hash("  ARG  ") == sequence_hash("ARG")

    def test_case_folding(self):
        # Lowercase folds to uppercase; mixed case hashes like uppercase.
        assert sequence_hash("arg") == sequence_hash("ARG")
        assert sequence_hash("ArG") == sequence_hash("ARG")

    def test_ambiguity_codes_preserved(self):
        # B/Z/X/U/O are amino-acid letters: kept literal, not expanded or dropped.
        # Two sequences differing only in an ambiguity code are distinct.
        assert sequence_hash("ARGX") != sequence_hash("ARG")
        assert sequence_hash("ARGX") != sequence_hash("ARGA")
        assert sequence_hash("ARGB") != sequence_hash("ARGZ")
        # But case-folding still applies to them.
        assert sequence_hash("argx") == sequence_hash("ARGX")

    def test_digits_and_punctuation_stripped(self):
        # FASTA rarely contains these, but defensive: they drop out cleanly.
        assert sequence_hash("AR1G") == sequence_hash("ARG")
        assert sequence_hash("AR-G") == sequence_hash("ARG")
        assert sequence_hash("AR.G") == sequence_hash("ARG")

    def test_empty_input(self):
        # Empty string is a valid (degenerate) input. Hashes to SHA256("").
        assert normalize_sequence("") == b""
        assert sequence_hash("") == hashlib.sha256(b"").hexdigest()

    def test_noise_only_input_equals_empty(self):
        # A "sequence" that is pure noise normalizes to empty, hashes as empty.
        # This is a footgun — callers should check length before hashing. The
        # hash function itself is not responsible for rejecting junk.
        assert sequence_hash("***") == sequence_hash("")
        assert sequence_hash("   ") == sequence_hash("")
        assert sequence_hash("123") == sequence_hash("")


# ---------------------------------------------------------------------------
# normalize_sequence returns bytes — callers that want to log/round-trip the
# canonical form need this.
# ---------------------------------------------------------------------------


class TestNormalizeSequence:
    def test_returns_bytes(self):
        assert isinstance(normalize_sequence("ARG"), bytes)

    def test_canonical_form_is_uppercase_only_AZ(self):
        result = normalize_sequence("  ar*g-1  ")
        assert result == b"ARG"

    def test_round_trip_through_ascii(self):
        # Decoding the normalized bytes as ASCII must succeed — we do not emit
        # any byte outside 0x41..=0x5A.
        result = normalize_sequence("The Quick Brown Fox 123 *")
        assert result.decode("ascii") == "THEQUICKBROWNFOX"


# ---------------------------------------------------------------------------
# Known-answer tests — pin a handful of hex digests. If these ever change,
# someone edited the contract.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "seq, expected_hex",
    [
        ("", "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"),
        ("M", "08f271887ce94707da822d5263bae19d5519cb3614e0daedc4c7ce5dab7473f1"),
        ("ARG", "204ca7e7cf7d3e9cd8905f36d40edbe3910102d8ae4dcdb9475b6d5fe959be3d"),
    ],
)
def test_known_answer_vectors(seq, expected_hex):
    # These are what Python produces today. The cross-language test will
    # verify Rust produces the same hex for the same inputs. If any of these
    # ever need updating, the contract has changed — bump HASH_VERSION.
    assert sequence_hash(seq) == expected_hex
