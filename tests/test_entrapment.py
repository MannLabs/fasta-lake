"""Tests for entrapment database generation (v3 feature)."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pytest

from fasta_lake.config import EntrapmentConfig
from fasta_lake.entrapment import (
    estimate_fdp_paired,
    generate_noble_paired,
    generate_shuffled_entrapment,
    load_entrapment_fasta,
    prepare_entrapment,
    validate_entrapment_fdr,
)

# ── generate_shuffled_entrapment ──────────────────────────────────────────────


class TestGenerateShuffledEntrapment:
    """Test shuffled entrapment generation."""

    def test_correct_count(self):
        source = {f"p{i}": "MSKGEELFTACDE" for i in range(100)}
        entrapment = generate_shuffled_entrapment(source, n_proteins=50)
        assert len(entrapment) == 50

    def test_fewer_than_requested(self):
        """If source has fewer proteins than requested, use all."""
        source = {"p1": "MSKGEELFT", "p2": "ACDEFGHIKL"}
        entrapment = generate_shuffled_entrapment(source, n_proteins=100)
        assert len(entrapment) == 2

    def test_aa_composition_preserved(self):
        """Shuffled sequences have the same AA composition as source."""
        source = {"p1": "MSKGEELFTACDEFGHIKLMNPQRSTVWY"}
        entrapment = generate_shuffled_entrapment(source, n_proteins=1)

        source_counts = Counter(source["p1"])
        entrap_key = list(entrapment.keys())[0]
        entrap_counts = Counter(entrapment[entrap_key])

        assert source_counts == entrap_counts

    def test_sequence_is_shuffled(self):
        """Shuffled sequence should (almost certainly) differ from original."""
        source = {"p1": "MSKGEELFTACDEFGHIKLMNPQRSTVWY" * 5}
        entrapment = generate_shuffled_entrapment(source, n_proteins=1)
        entrap_seq = list(entrapment.values())[0]
        # With 140 chars, probability of getting exact same order is ~1/140!
        assert entrap_seq != source["p1"]

    def test_deterministic_with_seed(self):
        """Same seed -> same result."""
        source = {f"p{i}": f"MSKGEELFT{i}" * 10 for i in range(20)}
        e1 = generate_shuffled_entrapment(source, n_proteins=10, seed=42)
        e2 = generate_shuffled_entrapment(source, n_proteins=10, seed=42)
        assert e1 == e2

    def test_different_seeds(self):
        """Different seeds -> different results."""
        source = {f"p{i}": "MSKGEELFTACDE" * 10 for i in range(20)}
        e1 = generate_shuffled_entrapment(source, n_proteins=10, seed=42)
        e2 = generate_shuffled_entrapment(source, n_proteins=10, seed=99)
        # At least some sequences should differ
        assert e1 != e2

    def test_prefix_in_ids(self):
        source = {"p1": "MSKGEELFT"}
        entrapment = generate_shuffled_entrapment(source, n_proteins=1, prefix="MY_PREFIX")
        assert all(k.startswith("MY_PREFIX_") for k in entrapment)

    def test_empty_source(self):
        entrapment = generate_shuffled_entrapment({}, n_proteins=10)
        assert entrapment == {}


# ── load_entrapment_fasta ─────────────────────────────────────────────────────


class TestLoadEntrapmentFasta:
    """Test loading user-supplied entrapment FASTA."""

    def test_loads_with_prefix(self, tmp_path: Path):
        fasta = tmp_path / "entrap.fasta"
        fasta.write_text(">prot1\nMSKGEELFT\n>prot2\nACDEFGHIKL\n")

        result = load_entrapment_fasta(fasta, prefix="ENTRAP")

        assert len(result) == 2
        assert "ENTRAP_prot1" in result
        assert "ENTRAP_prot2" in result

    def test_loads_without_prefix(self, tmp_path: Path):
        fasta = tmp_path / "entrap.fasta"
        fasta.write_text(">prot1\nMSKGEELFT\n")

        result = load_entrapment_fasta(fasta, prefix="")

        assert "prot1" in result

    def test_sequences_preserved(self, tmp_path: Path):
        fasta = tmp_path / "entrap.fasta"
        fasta.write_text(">prot1\nMSKGEELFT\n")

        result = load_entrapment_fasta(fasta, prefix="E")

        assert result["E_prot1"] == "MSKGEELFT"


# ── generate_noble_paired ─────────────────────────────────────────────────────


class TestGenerateNoblePaired:
    """Test Noble paired entrapment partitioning."""

    def test_partition_sizes(self):
        proteins = {f"p{i}": f"SEQ{i}" for i in range(100)}
        target, entrapment = generate_noble_paired(proteins)

        assert len(target) == 50
        assert len(entrapment) == 50

    def test_odd_number(self):
        proteins = {f"p{i}": f"SEQ{i}" for i in range(101)}
        target, entrapment = generate_noble_paired(proteins)

        assert len(target) + len(entrapment) == 101

    def test_no_overlap(self):
        """Target and entrapment should have no overlapping sequences."""
        proteins = {f"p{i}": f"UNIQUE_SEQ_{i}" for i in range(50)}
        target, entrapment = generate_noble_paired(proteins)

        target_seqs = set(target.values())
        entrap_seqs = set(entrapment.values())
        assert len(target_seqs & entrap_seqs) == 0

    def test_entrapment_prefixed(self):
        proteins = {"p1": "SEQ1", "p2": "SEQ2"}
        _, entrapment = generate_noble_paired(proteins)

        for key in entrapment:
            assert key.startswith("ENTRAP_NOBLE_")

    def test_deterministic(self):
        proteins = {f"p{i}": f"SEQ{i}" for i in range(20)}
        t1, e1 = generate_noble_paired(proteins, seed=42)
        t2, e2 = generate_noble_paired(proteins, seed=42)
        assert t1 == t2
        assert e1 == e2

    def test_all_proteins_accounted(self):
        """All original sequences appear in either target or entrapment."""
        proteins = {f"p{i}": f"SEQ{i}" for i in range(20)}
        target, entrapment = generate_noble_paired(proteins)

        all_seqs = set(target.values()) | set(entrapment.values())
        assert all_seqs == set(proteins.values())


# ── validate_entrapment_fdr ───────────────────────────────────────────────────


class TestValidateEntrapmentFDR:
    """Test FDR validation using entrapment hits.

    Uses the combined method (Eq 1, Wen/Keich/Noble 2025 Nat Methods):
        FDP = N_E × (1 + 1/r) / (N_T + N_E)
    where r = entrapment_db_size / target_db_size.
    """

    def test_passes_below_threshold(self):
        # r=1 (equal DBs): combined = 2 × lower = 2 × 5/1000 = 0.01
        result = validate_entrapment_fdr(
            entrapment_hits=5,
            target_hits=995,
            target_db_size=10000,
            entrapment_db_size=10000,
            expected_fdr=0.01,
        )
        assert result["passed"] is True
        assert result["fdp_combined"] <= 0.01

    def test_fails_above_threshold(self):
        result = validate_entrapment_fdr(
            entrapment_hits=50,
            target_hits=950,
            target_db_size=10000,
            entrapment_db_size=10000,
            expected_fdr=0.01,
        )
        assert result["passed"] is False
        assert "warning" in result

    def test_zero_hits(self):
        result = validate_entrapment_fdr(
            entrapment_hits=0,
            target_hits=1000,
            target_db_size=10000,
            entrapment_db_size=5000,
            expected_fdr=0.01,
        )
        assert result["passed"] is True
        assert result["fdp_combined"] == 0.0

    def test_zero_total(self):
        """No discoveries at all means the test never ran -- UNDEFINED, not a pass.

        This previously asserted ``passed is True`` / ``fdp_combined == 0.0``,
        which pinned a vacuous pass: 0.0 is the most favourable value the
        statistic can take, standing in for a measurement that was never
        possible. The assertion now requires the honest verdict.
        """
        result = validate_entrapment_fdr(
            entrapment_hits=0,
            target_hits=0,
            target_db_size=10000,
            entrapment_db_size=5000,
            expected_fdr=0.01,
        )
        assert result["evaluable"] is False
        assert result["passed"] is None
        assert result["verdict"] == "UNDEFINED"
        assert result["fdp_combined"] is None
        assert "warning" in result

    def test_empty_entrapment_partition_is_undefined(self):
        """Nothing to hit means nothing was tested, however many targets there are."""
        result = validate_entrapment_fdr(
            entrapment_hits=0,
            target_hits=5000,
            target_db_size=10000,
            entrapment_db_size=0,
            expected_fdr=0.01,
        )
        assert result["evaluable"] is False
        assert result["passed"] is None
        assert result["fdp_combined"] is None

    def test_genuine_zero_still_passes(self):
        """Zero hits against a populated entrapment partition is a real 0% pass."""
        result = validate_entrapment_fdr(
            entrapment_hits=0,
            target_hits=1000,
            target_db_size=10000,
            entrapment_db_size=5000,
            expected_fdr=0.01,
        )
        assert result["evaluable"] is True
        assert result["passed"] is True
        assert result["verdict"] == "PASS"
        assert result["fdp_combined"] == 0.0

    def test_exact_threshold(self):
        # r=1: combined = 2 × 10/1000 = 0.02 > 0.01, so fails
        result = validate_entrapment_fdr(
            entrapment_hits=10,
            target_hits=990,
            target_db_size=10000,
            entrapment_db_size=10000,
            expected_fdr=0.01,
        )
        # At r=1, 10 hits gives FDP=0.02, which exceeds 0.01
        assert result["passed"] is False

    def test_r_scaling(self):
        """When entrapment DB is small relative to target (r < 1),
        combined FDP is larger than lower bound."""
        result = validate_entrapment_fdr(
            entrapment_hits=5,
            target_hits=995,
            target_db_size=100000,  # large target
            entrapment_db_size=5000,  # small entrapment
            expected_fdr=0.10,
        )
        # r = 5000/100000 = 0.05, so (1+1/r) = 21
        # FDP = 5 × 21 / 1000 = 0.105
        assert result["r"] < 1.0
        assert result["fdp_combined"] > result["fdp_lower_bound"]


# ── estimate_fdp_paired (Wen 2025 Eq 4) ───────────────────────────────────────


class TestEstimateFdpPaired:
    def test_formula_matches_wen_eq4(self):
        # FDP* = (N_E + N_{E>=s>T} + 2*N_{E>T>=s}) / (N_T + N_E)
        r = estimate_fdp_paired(
            n_target=980,
            n_entrapment=20,
            n_entrap_target_undiscovered=6,
            n_entrap_gt_target=4,
        )
        # numerator = 20 + 6 + 2*4 = 34 ; denom = 980 + 20 = 1000
        assert r["numerator"] == 34
        assert r["fdp_paired"] == pytest.approx(0.034)
        assert r["evaluable"] is True

    def test_no_target_score_advantage_reduces_to_lower_bound(self):
        # When every discovered entrapment peptide's paired target was also
        # discovered and scored at least as high (both sub-counts zero), the
        # numerator collapses to N_E and the paired estimate equals the
        # lower bound N_E/(N_T+N_E). This is the tightest, best-case reading.
        r = estimate_fdp_paired(
            n_target=990,
            n_entrapment=10,
            n_entrap_target_undiscovered=0,
            n_entrap_gt_target=0,
        )
        assert r["numerator"] == 10
        assert r["fdp_paired"] == pytest.approx(10 / 1000)

    def test_upper_bound_never_below_lower_bound(self):
        # Paired numerator (N_E + a + 2b) >= N_E, so paired >= lower bound.
        n_t, n_e = 900, 100
        lower = n_e / (n_t + n_e)
        r = estimate_fdp_paired(
            n_target=n_t,
            n_entrapment=n_e,
            n_entrap_target_undiscovered=30,
            n_entrap_gt_target=20,
        )
        assert r["fdp_paired"] >= lower

    def test_zero_discoveries_is_undefined_not_zero(self):
        r = estimate_fdp_paired(
            n_target=0,
            n_entrapment=0,
            n_entrap_target_undiscovered=0,
            n_entrap_gt_target=0,
        )
        assert r["evaluable"] is False
        assert r["fdp_paired"] is None

    def test_subcount_exceeds_entrapment_raises(self):
        with pytest.raises(ValueError):
            estimate_fdp_paired(
                n_target=100,
                n_entrapment=5,
                n_entrap_target_undiscovered=6,  # > N_E
                n_entrap_gt_target=0,
            )

    def test_disjoint_subcount_sum_exceeds_entrapment_raises(self):
        with pytest.raises(ValueError):
            estimate_fdp_paired(
                n_target=100,
                n_entrapment=10,
                n_entrap_target_undiscovered=7,
                n_entrap_gt_target=5,  # 7 + 5 > 10
            )


# ── prepare_entrapment ────────────────────────────────────────────────────────


class TestPrepareEntrapment:
    """Test entrapment dispatcher."""

    def test_disabled_returns_empty(self):
        config = EntrapmentConfig(enabled=False)
        result = prepare_entrapment(config)
        assert result == {}

    def test_shuffled_method(self):
        config = EntrapmentConfig(
            enabled=True,
            method="shuffled",
            n_proteins=5,
        )
        source = {f"p{i}": "MSKGEELFTACDE" for i in range(10)}
        result = prepare_entrapment(config, source_proteins=source)
        assert len(result) == 5

    def test_shuffled_no_source_raises(self):
        config = EntrapmentConfig(enabled=True, method="shuffled")
        with pytest.raises(ValueError, match="source_proteins required"):
            prepare_entrapment(config)

    def test_shuffled_uses_configured_seed(self):
        """config.seed must reach the generator: a different seed gives a different
        shuffle, the same seed reproduces it, and seed 7 is not silently 42."""
        source = {f"p{i}": "MSKGEELFTACDEKRAAWQFTGEDPLMNK" for i in range(4)}
        cfg7 = EntrapmentConfig(enabled=True, method="shuffled", n_proteins=4, seed=7)
        cfg42 = EntrapmentConfig(enabled=True, method="shuffled", n_proteins=4, seed=42)
        seed7_a = prepare_entrapment(cfg7, source_proteins=source)
        seed7_b = prepare_entrapment(cfg7, source_proteins=source)
        seed42 = prepare_entrapment(cfg42, source_proteins=source)
        assert seed7_a == seed7_b
        assert seed7_a != seed42
        assert seed42 == generate_shuffled_entrapment(source, n_proteins=4, seed=42)
        assert seed7_a == generate_shuffled_entrapment(source, n_proteins=4, seed=7)

    def test_user_method(self, tmp_path: Path):
        fasta = tmp_path / "entrap.fasta"
        fasta.write_text(">prot1\nMSKGEELFT\n>prot2\nACDEFGHIKL\n")

        config = EntrapmentConfig(
            enabled=True,
            method="user",
            source=str(fasta),
        )
        result = prepare_entrapment(config)
        assert len(result) == 2

    def test_phylogenetic_method(self, tmp_path: Path):
        fasta = tmp_path / "arabidopsis.fasta"
        fasta.write_text(">AT1G01010\nMSKGEELFT\n>AT1G01020\nACDEFGHIKL\n>AT1G01030\nMNPQRSTVWY\n")

        config = EntrapmentConfig(
            enabled=True,
            method="phylogenetic",
            source=str(fasta),
            n_proteins=2,
        )
        result = prepare_entrapment(config)
        assert len(result) == 2

    def test_missing_source_raises(self):
        config = EntrapmentConfig(enabled=True, method="phylogenetic", source="")
        with pytest.raises(ValueError, match="source.*required"):
            prepare_entrapment(config)

    def test_unknown_method_raises(self):
        config = EntrapmentConfig(enabled=True, method="invalid")
        with pytest.raises(ValueError, match="Unknown entrapment method"):
            prepare_entrapment(config)
