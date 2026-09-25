"""Tests for sequence quality filtering (v3 feature)."""

from __future__ import annotations

import math

from fasta_lake.config import SequenceQualityConfig
from fasta_lake.quality import (
    compute_sequence_entropy,
    filter_sequences_by_quality,
    has_poly_run,
    is_fragment,
)

# ── has_poly_run ──────────────────────────────────────────────────────────────


class TestHasPolyRun:
    """Test homopolymeric run detection."""

    def test_poly_a_detected(self):
        assert has_poly_run("MSK" + "A" * 10 + "GEELFT", run_length=10) is True

    def test_poly_a_below_threshold(self):
        assert has_poly_run("MSK" + "A" * 9 + "GEELFT", run_length=10) is False

    def test_exact_threshold(self):
        assert has_poly_run("A" * 10, run_length=10) is True

    def test_no_poly_run(self):
        assert has_poly_run("MSKGEELFTACDEFGH", run_length=10) is False

    def test_empty_string(self):
        assert has_poly_run("", run_length=10) is False

    def test_short_sequence(self):
        assert has_poly_run("AAA", run_length=10) is False

    def test_multiple_short_runs(self):
        # Multiple short runs that individually don't trigger
        assert has_poly_run("AAAAA" + "G" + "BBBBB", run_length=6) is False

    def test_poly_at_start(self):
        assert has_poly_run("AAAAAAAAAAMSKGE", run_length=10) is True

    def test_poly_at_end(self):
        assert has_poly_run("MSKGEAAAAAAAAAA", run_length=10) is True

    def test_different_amino_acids(self):
        for aa in "ACDEFGHKLMNPQRSTVWY":
            assert has_poly_run(aa * 10, run_length=10) is True

    def test_run_length_2(self):
        assert has_poly_run("AA", run_length=2) is True
        assert has_poly_run("AB", run_length=2) is False


# ── compute_sequence_entropy ──────────────────────────────────────────────────


class TestComputeSequenceEntropy:
    """Test Shannon entropy calculation."""

    def test_uniform_sequence(self):
        """All same amino acid = 0 entropy."""
        assert compute_sequence_entropy("AAAAAAAAAA") == 0.0

    def test_max_diversity(self):
        """All 20 amino acids equally = entropy of 1.0."""
        seq = "ACDEFGHIKLMNPQRSTVWY"
        entropy = compute_sequence_entropy(seq)
        assert abs(entropy - 1.0) < 0.01

    def test_two_amino_acids(self):
        """Two equally frequent AAs -> entropy = log2(2)/log2(20)."""
        seq = "AAAAAAAAAALLLLLLLLL"  # ~half A, ~half L
        entropy = compute_sequence_entropy(seq)
        expected = math.log2(2) / math.log2(20)
        assert abs(entropy - expected) < 0.05

    def test_empty_string(self):
        assert compute_sequence_entropy("") == 0.0

    def test_single_char(self):
        assert compute_sequence_entropy("A") == 0.0

    def test_realistic_protein(self):
        """A realistic protein sequence should have moderate-high entropy."""
        seq = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGK"
        entropy = compute_sequence_entropy(seq)
        assert 0.5 < entropy < 1.0

    def test_low_complexity(self):
        """90% one amino acid = very low entropy."""
        seq = "A" * 90 + "BCDEFGHIKL"
        entropy = compute_sequence_entropy(seq)
        assert entropy < 0.3

    def test_entropy_increases_with_diversity(self):
        """More diverse sequences have higher entropy."""
        seq_low = "A" * 100
        seq_mid = "A" * 50 + "L" * 50
        seq_high = "ACDEFGHIKLMNPQRSTVWY" * 5
        e_low = compute_sequence_entropy(seq_low)
        e_mid = compute_sequence_entropy(seq_mid)
        e_high = compute_sequence_entropy(seq_high)
        assert e_low < e_mid < e_high


# ── is_fragment ───────────────────────────────────────────────────────────────


class TestIsFragment:
    """Test fragment detection."""

    def test_short_no_m_start(self):
        assert is_fragment("SKGEELFT", length_threshold=100) is True

    def test_short_with_m_start(self):
        assert is_fragment("MSKGEELFT", length_threshold=100) is False

    def test_long_no_m_start(self):
        """Long sequences without M-start are NOT fragments."""
        seq = "SKGEELFT" * 20  # 160 aa
        assert is_fragment(seq, length_threshold=100) is False

    def test_empty_string(self):
        assert is_fragment("", length_threshold=100) is True

    def test_exact_threshold(self):
        # 99 chars, no M-start -> fragment
        assert is_fragment("S" * 99, length_threshold=100) is True
        # 100 chars, no M-start -> NOT fragment (>= threshold)
        assert is_fragment("S" * 100, length_threshold=100) is False

    def test_m_start_always_passes(self):
        assert is_fragment("M", length_threshold=100) is False
        assert is_fragment("MSKGE", length_threshold=100) is False


# ── filter_sequences_by_quality ───────────────────────────────────────────────


class TestFilterSequencesByQuality:
    """Test batch quality filtering."""

    def test_defaults_no_optional_filters(self):
        """v5.1 default config: min_length=30, poly_runs on. Use custom config for these
        short test peptides."""
        # v5.1 bumped defaults — use explicit short-friendly config here
        config = SequenceQualityConfig(min_length=1, reject_poly_runs=False)
        proteins = {
            "p1": "MSKGEELFT",
            "p2": "MSKGEELFTACDEFGHKLMNPQRSTVWY",
        }
        filtered, report = filter_sequences_by_quality(proteins, config)
        assert len(filtered) == 2
        assert report.passed == 2
        assert report.rejected_poly_run == 0

    def test_length_filter_min(self):
        config = SequenceQualityConfig(min_length=10)
        proteins = {"short": "MSKGE", "long": "MSKGEELFTACDE"}
        filtered, report = filter_sequences_by_quality(proteins, config)
        assert "short" not in filtered
        assert "long" in filtered
        assert report.rejected_length == 1

    def test_length_filter_max(self):
        config = SequenceQualityConfig(min_length=1, max_length=10)
        proteins = {"short": "MSKGE", "long": "MSKGEELFTACDE"}
        filtered, report = filter_sequences_by_quality(proteins, config)
        assert "short" in filtered
        assert "long" not in filtered

    def test_poly_run_filter(self):
        config = SequenceQualityConfig(min_length=1, reject_poly_runs=True, poly_run_length=5)
        proteins = {
            "good": "MSKGEELFT",
            "poly": "MSKAAAAAAGEELFT",
        }
        filtered, report = filter_sequences_by_quality(proteins, config)
        assert "good" in filtered
        assert "poly" not in filtered
        assert report.rejected_poly_run == 1

    def test_low_complexity_filter(self):
        config = SequenceQualityConfig(
            min_length=1,
            reject_poly_runs=False,
            reject_low_complexity=True,
            low_complexity_threshold=0.3,
        )
        proteins = {
            "good": "MSKGEELFTACDEFGHKLMNPQRSTVWY",
            "low": "A" * 90 + "BCDEFGHIKL",
        }
        filtered, report = filter_sequences_by_quality(proteins, config)
        assert "good" in filtered
        assert "low" not in filtered
        assert report.rejected_low_complexity == 1

    def test_fragment_filter(self):
        config = SequenceQualityConfig(
            min_length=1,
            reject_poly_runs=False,
            reject_fragments=True,
            fragment_length_threshold=50,
        )
        proteins = {
            "good": "MSKGEELFTACDE",
            "frag": "SKGEELFT",  # no M, < 50 aa
        }
        filtered, report = filter_sequences_by_quality(proteins, config)
        assert "good" in filtered
        assert "frag" not in filtered
        assert report.rejected_fragment == 1

    def test_all_filters_combined(self):
        config = SequenceQualityConfig(
            min_length=5,
            max_length=200,
            reject_poly_runs=True,
            poly_run_length=5,
            reject_low_complexity=True,
            low_complexity_threshold=0.3,
            reject_fragments=True,
            fragment_length_threshold=50,
        )
        proteins = {
            "good": "MSKGEELFTACDEFGHKLMNPQRSTVWY",
            "too_short": "MSK",
            "poly": "MSKAAAAAAGEELFT",
            "low_entropy": "M" + "A" * 90 + "B" * 9,
            "fragment": "SKGEELFT",
        }
        filtered, report = filter_sequences_by_quality(proteins, config)
        assert filtered.keys() == {"good"}
        assert report.passed == 1
        assert report.total_input == 5

    def test_first_rejection_wins(self):
        """A sequence rejected by length is NOT checked for poly-runs."""
        config = SequenceQualityConfig(
            min_length=100,
            reject_poly_runs=True,
            poly_run_length=5,
        )
        # This sequence has poly-run AND is too short
        proteins = {"both": "MAAAAAAGEELFT"}
        filtered, report = filter_sequences_by_quality(proteins, config)
        # Should be counted as length rejection, not poly
        assert report.rejected_length == 1
        assert report.rejected_poly_run == 0

    def test_empty_input(self):
        filtered, report = filter_sequences_by_quality({})
        assert len(filtered) == 0
        assert report.total_input == 0
        assert report.passed == 0

    def test_report_summary(self):
        config = SequenceQualityConfig(min_length=1, reject_poly_runs=True, poly_run_length=5)
        proteins = {"good": "MSKGEELFT", "poly": "MSKAAAAAAGEELFT"}
        _, report = filter_sequences_by_quality(proteins, config)
        summary = report.summary()
        assert summary["total_input"] == 2
        assert summary["passed"] == 1
        assert summary["total_rejected"] == 1

    def test_rejection_details_populated(self):
        config = SequenceQualityConfig(min_length=1, reject_poly_runs=True, poly_run_length=5)
        proteins = {"poly": "MSKAAAAAAGEELFT"}
        _, report = filter_sequences_by_quality(proteins, config)
        assert "poly" in report.rejection_details
        assert "poly_run" in report.rejection_details["poly"]

    def test_none_config_uses_defaults(self):
        # v5.1 default min_length=30 — use a long enough sequence
        proteins = {"p1": "MSKGEELFTACDEFGHIKLMNPQRSTVWYACDEFG"}  # 35 AA
        filtered, _ = filter_sequences_by_quality(proteins, None)
        assert len(filtered) == 1

    # v5.1 additions
    def test_noncanonical_residues_stripped(self):
        config = SequenceQualityConfig(
            min_length=1,
            reject_poly_runs=False,
            strip_noncanonical=True,
            reject_if_noncanonical_fraction=0.5,  # permissive for test
        )
        proteins = {"p1": "MSKBUZGEELFT"}  # B and Z are ambiguous; U is selenocysteine
        filtered, report = filter_sequences_by_quality(proteins, config)
        assert "p1" in filtered
        # Replace B and Z while retaining U.
        assert filtered["p1"] == "MSKXUXGEELFT"
        assert report.residues_normalised == 2

    def test_noncanonical_reject_high_fraction(self):
        config = SequenceQualityConfig(
            min_length=1,
            reject_poly_runs=False,
            strip_noncanonical=False,  # keep original so fraction is high
            reject_if_noncanonical_fraction=0.1,
        )
        # 3 X in 10 = 30% non-canonical
        proteins = {"bad": "MSKXXXGEEL", "good": "MSKGEELFTA"}
        filtered, report = filter_sequences_by_quality(proteins, config)
        assert "good" in filtered
        assert "bad" not in filtered
        assert report.rejected_noncanonical == 1


def test_quality_filters_preserve_encoded_rare_residues():
    # Both rare residues exceed the configured ambiguity-fraction cutoff here.
    # They must remain searchable instead of causing rejection or becoming X.
    sequence = "MACDEFGUOKLMNPQRSTVWYACDEFGUOKLMNPQRSTVWY"
    filtered, report = filter_sequences_by_quality({"selenoprotein": sequence})
    assert filtered == {"selenoprotein": sequence}
    assert report.rejected_noncanonical == 0
    assert report.residues_normalised == 0
