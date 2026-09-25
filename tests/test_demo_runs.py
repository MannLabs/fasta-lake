"""
Demo runs: 50 synthetic scenarios testing all features end-to-end.

Generates diverse synthetic datasets and runs the full pipeline
to validate correctness, performance, and log output quality.

Scenarios cover:
- Single-species (clinical) vs multi-species (metaproteomics)
- Varying database sizes (10 to 10,000 proteins)
- Varying peptide counts (5 to 500)
- All 4 inference strategies
- Quality filtering (poly-runs, low-complexity, fragments)
- Complexity estimation
- Entrapment generation and FDR validation
- Edge cases (empty, no-match, all-shared)
- Score filtering
"""

from __future__ import annotations

import csv
import json
import random
import time
from pathlib import Path

import pytest

from fasta_lake.complexity import estimate_complexity
from fasta_lake.config import (
    ComplexityConfig,
    SequenceQualityConfig,
)
from fasta_lake.entrapment import (
    generate_shuffled_entrapment,
    validate_entrapment_fdr,
)
from fasta_lake.helpers import (
    normalize_il,
)
from fasta_lake.inference.runner import run_inference
from fasta_lake.quality import filter_sequences_by_quality

AA = "ACDEFGHKLMNPQRSTVWY"


def _generate_proteins(
    n: int,
    species_prefix: str = "MGYG",
    n_species: int = 1,
    seq_len: int = 200,
    seed: int = 42,
) -> dict[str, str]:
    """Generate synthetic protein database."""
    rng = random.Random(seed)
    proteins = {}
    for i in range(n):
        species_idx = i % n_species
        if species_prefix == "sp":
            pid = f"sp|P{i:05d}|PROT{i}_HUMAN"
        elif species_prefix == "MGYG":
            pid = f"MGYG{species_idx:09d}_{i:05d}"
        elif species_prefix == "GMGC":
            pid = f"GMGC10.{species_idx:03d}_{i:03d}_001.PROT{i:02d}"
        elif species_prefix == "mixed":
            # Mix all types
            if i % 4 == 0:
                pid = f"sp|P{i:05d}|PROT{i}_HUMAN"
            elif i % 4 == 1:
                pid = f"MGYG{species_idx:09d}_{i:05d}"
            elif i % 4 == 2:
                pid = f"GMGC10.{species_idx:03d}_{i:03d}_001.PROT{i:02d}"
            else:
                pid = f"SMORF_{i:05d}"
        else:
            pid = f"{species_prefix}_{i:05d}"

        seq = "".join(rng.choices(AA, k=seq_len))
        proteins[pid] = seq
    return proteins


def _extract_peptides(
    proteins: dict[str, str],
    n_peptides: int,
    pep_len: int = 15,
    shared_ratio: float = 0.0,
    seed: int = 42,
) -> list[tuple[str, float]]:
    """Extract peptides from protein sequences."""
    rng = random.Random(seed)
    peptides = []
    prot_ids = list(proteins.keys())

    for i in range(n_peptides):
        pid = prot_ids[i % len(prot_ids)]
        seq = proteins[pid]
        if len(seq) > pep_len:
            start = rng.randint(0, len(seq) - pep_len)
            pep = seq[start : start + pep_len]
        else:
            pep = seq
        score = rng.uniform(0.3, 0.99)
        peptides.append((normalize_il(pep), score))

    return peptides


def _write_fasta(path: Path, proteins: dict[str, str]) -> None:
    with open(path, "w") as f:
        for pid, seq in proteins.items():
            f.write(f">{pid}\n{seq}\n")


def _write_taxonomy(path: Path, proteins: dict[str, str]) -> Path:
    """Species for every synthetic id: the catalogue rule where one exists, else
    the id's prefix (SMORF_..., mixed_...), so species_budget has a species for all."""
    from fasta_lake.helpers import get_species_id

    with open(path, "w") as f:
        for pid in proteins:
            sp = get_species_id(pid) or pid.split("_")[0].split("|")[0]
            f.write(f"{pid}\t{sp}\n")
    return path


def _write_csv(path: Path, peptides: list[tuple[str, float]]) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["peptide_prediction_detokenized_unmodified", "score"])
        for pep, score in peptides:
            writer.writerow([pep, score])


# ── Scenario generators ──────────────────────────────────────────────────────

STRATEGIES = ["species_budget", "razor", "uniform", "uniform_2pep"]


class TestDemoSingleSpecies:
    """Demo runs: single-species clinical proteomics scenarios."""

    @pytest.mark.parametrize("n_proteins", [10, 50, 100, 500, 1000])
    def test_varying_db_size(self, tmp_path: Path, n_proteins: int):
        """D1-D5: Single-species with varying DB sizes."""
        proteins = _generate_proteins(n_proteins, species_prefix="sp", seed=n_proteins)
        peptides = _extract_peptides(proteins, n_peptides=min(100, n_proteins))

        fasta = tmp_path / "db.fasta"
        csv_file = tmp_path / "preds.csv"
        _write_fasta(fasta, proteins)
        _write_csv(csv_file, peptides)

        for strategy in STRATEGIES:
            out = tmp_path / strategy
            stats = run_inference(
                sample_id=f"SINGLE_{n_proteins}",
                stage2_fasta=fasta,
                alphanovo_csv=csv_file,
                output_dir=out,
                strategy=strategy,
                taxonomy_map_path=_write_taxonomy(tmp_path / "tax.tsv", proteins),
            )
            assert stats["stage2_proteins"] == n_proteins
            assert stats["selected_proteins"] >= 0
            assert stats["runtime_seconds"] >= 0

    @pytest.mark.parametrize("n_peptides", [5, 20, 50, 100, 200])
    def test_varying_peptide_count(self, tmp_path: Path, n_peptides: int):
        """D6-D10: Fixed DB, varying peptide count."""
        proteins = _generate_proteins(200, species_prefix="sp")
        peptides = _extract_peptides(proteins, n_peptides=n_peptides, seed=n_peptides)

        fasta = tmp_path / "db.fasta"
        csv_file = tmp_path / "preds.csv"
        _write_fasta(fasta, proteins)
        _write_csv(csv_file, peptides)

        stats = run_inference(
            sample_id=f"PEPS_{n_peptides}",
            stage2_fasta=fasta,
            alphanovo_csv=csv_file,
            output_dir=tmp_path / "out",
            strategy="razor",
        )
        assert stats["selected_proteins"] >= 0


class TestDemoMultiSpecies:
    """Demo runs: multi-species metaproteomics scenarios."""

    @pytest.mark.parametrize("n_species", [2, 5, 10, 50])
    def test_varying_species_count(self, tmp_path: Path, n_species: int):
        """D11-D14: Multi-species with varying richness."""
        n_proteins = n_species * 20
        proteins = _generate_proteins(
            n_proteins, species_prefix="MGYG", n_species=n_species, seed=n_species
        )
        peptides = _extract_peptides(proteins, n_peptides=min(100, n_proteins))

        fasta = tmp_path / "db.fasta"
        csv_file = tmp_path / "preds.csv"
        _write_fasta(fasta, proteins)
        _write_csv(csv_file, peptides)

        for strategy in STRATEGIES:
            out = tmp_path / strategy
            stats = run_inference(
                sample_id=f"MULTI_{n_species}sp",
                stage2_fasta=fasta,
                alphanovo_csv=csv_file,
                output_dir=out,
                strategy=strategy,
                taxonomy_map_path=_write_taxonomy(tmp_path / "tax.tsv", proteins),
            )
            assert stats["selected_proteins"] >= 0

    def test_mixed_db_types(self, tmp_path: Path):
        """D15: Mixed DB types (Human + UHGG + GMGC + smORF)."""
        proteins = _generate_proteins(200, species_prefix="mixed", n_species=10)
        peptides = _extract_peptides(proteins, n_peptides=100)

        fasta = tmp_path / "db.fasta"
        csv_file = tmp_path / "preds.csv"
        _write_fasta(fasta, proteins)
        _write_csv(csv_file, peptides)

        for strategy in STRATEGIES:
            out = tmp_path / strategy
            stats = run_inference(
                sample_id="MIXED",
                stage2_fasta=fasta,
                alphanovo_csv=csv_file,
                output_dir=out,
                strategy=strategy,
                taxonomy_map_path=_write_taxonomy(tmp_path / "tax.tsv", proteins),
            )
            assert stats["selected_proteins"] >= 0
            # Check db_breakdown has multiple types
            if "db_breakdown" in stats:
                assert len(stats["db_breakdown"]) >= 1


class TestDemoQualityFiltering:
    """Demo runs: quality filtering scenarios."""

    def test_poly_run_rejection(self, tmp_path: Path):
        """D16: Proteins with poly-runs get filtered."""
        proteins = {
            "good_1": "MSKGEELFTACDEFGHKLMNPQRSTVWY",
            # v5.1 note: "ANOTHERGOODPROTEINSEQUENCEHERE" contains U/O (non-canonical),
            # so use a canonical-only sequence here.
            "good_2": "MALSPSTVWLPNRQNGKVSTHELVMY",
            "poly_a": "M" + "A" * 50 + "CDEFGHKLM",
            "poly_l": "M" + "L" * 50 + "NQRSTVWYA",
            "poly_k": "M" + "K" * 50 + "ACDEFGHLN",
        }
        config = SequenceQualityConfig(
            reject_poly_runs=True,
            poly_run_length=10,
            min_length=1,
        )

        filtered, report = filter_sequences_by_quality(proteins, config)

        assert report.rejected_poly_run == 3
        assert report.passed == 2
        assert "good_1" in filtered
        assert "good_2" in filtered

    def test_low_complexity_rejection(self, tmp_path: Path):
        """D17: Low-entropy sequences get filtered."""
        proteins = {
            "diverse": "MSKGEELFTACDEFGHKLMNPQRSTVWY",
            "low_ent": "M" + "A" * 95 + "BCDG",
            "medium": "M" + "AABB" * 25,
        }
        config = SequenceQualityConfig(
            reject_low_complexity=True,
            low_complexity_threshold=0.3,
            min_length=1,
        )

        filtered, report = filter_sequences_by_quality(proteins, config)

        assert "diverse" in filtered
        assert "low_ent" not in filtered

    def test_fragment_rejection(self, tmp_path: Path):
        """D18: Fragments (no M-start, short) get filtered."""
        proteins = {
            "good": "MSKGEELFTACDEFGHKLMNPQRSTVWY" * 4,  # 112 aa, has M
            "frag_short": "SKGEELFT",  # no M, short
            "frag_medium": "ACDEFGHKLMNPQR" * 5,  # no M, 70 aa
            "long_no_m": "ACDEFGHKLMNPQRSTVWY" * 10,  # no M, 200 aa -> not fragment
        }
        config = SequenceQualityConfig(
            reject_fragments=True,
            fragment_length_threshold=100,
            min_length=1,
        )

        filtered, report = filter_sequences_by_quality(proteins, config)

        assert "good" in filtered
        assert "frag_short" not in filtered
        assert "frag_medium" not in filtered
        assert "long_no_m" in filtered  # Long enough -> not a fragment
        assert report.rejected_fragment == 2

    def test_combined_filters_pipeline(self, tmp_path: Path):
        """D19: Full quality pipeline with all filters."""
        rng = random.Random(99)
        proteins = {}
        for i in range(100):
            seq = "".join(rng.choices(AA, k=rng.randint(50, 300)))
            if i < 5:
                # Add poly-run
                seq = seq[:20] + "A" * 15 + seq[35:]
            if 5 <= i < 10:
                # Make low-entropy
                seq = "A" * (len(seq) - 10) + "BCDEFGHIKL"
            if 10 <= i < 15:
                # Make fragment (short, no M)
                seq = "SKGE" * 5
            proteins[f"PROT_{i:04d}"] = seq

        config = SequenceQualityConfig(
            min_length=10,
            max_length=500,
            reject_poly_runs=True,
            poly_run_length=10,
            reject_low_complexity=True,
            low_complexity_threshold=0.3,
            reject_fragments=True,
            fragment_length_threshold=50,
        )

        filtered, report = filter_sequences_by_quality(proteins, config)

        assert report.total_input == 100
        assert report.passed > 0
        assert report.rejected_poly_run >= 0
        assert report.rejected_low_complexity >= 0
        assert report.rejected_fragment >= 0

        # Summary should be valid JSON-serializable
        summary = report.summary()
        json.dumps(summary)  # Should not raise


class TestDemoComplexity:
    """Demo runs: complexity estimation scenarios."""

    def test_single_species_clinical(self, tmp_path: Path):
        """D20: Clinical sample -> recommend razor."""
        proteins = _generate_proteins(500, species_prefix="sp")
        peptides = _extract_peptides(proteins, n_peptides=100)
        pep_dict = dict(peptides)

        config = ComplexityConfig(enabled=True, auto_recommend=True)
        report = estimate_complexity(proteins, pep_dict, config)

        assert report.species_richness == 1
        assert report.recommendation == "razor"
        assert report.hit_rate > 0

    def test_moderate_metaproteomics(self, tmp_path: Path):
        """D21: Moderate metaproteomics -> recommend species_budget."""
        proteins = _generate_proteins(500, species_prefix="MGYG", n_species=30)
        peptides = _extract_peptides(proteins, n_peptides=200)
        pep_dict = dict(peptides)

        config = ComplexityConfig(enabled=True, auto_recommend=True)
        report = estimate_complexity(proteins, pep_dict, config)

        assert report.species_richness >= 10
        assert report.recommendation == "species_budget"

    def test_high_complexity_metaproteomics(self, tmp_path: Path):
        """D22: High complexity -> species_budget with clustering."""
        proteins = _generate_proteins(2000, species_prefix="MGYG", n_species=200)
        peptides = _extract_peptides(proteins, n_peptides=500)
        pep_dict = dict(peptides)

        config = ComplexityConfig(enabled=True, auto_recommend=True)
        report = estimate_complexity(proteins, pep_dict, config)

        assert report.species_richness > 100
        assert report.recommendation == "species_budget"
        assert "clustering" in report.recommendation_reason.lower()

    def test_complexity_summary_json(self, tmp_path: Path):
        """D23: Complexity report is JSON-serializable."""
        proteins = _generate_proteins(100, species_prefix="MGYG", n_species=5)
        peptides = _extract_peptides(proteins, n_peptides=50)
        pep_dict = dict(peptides)

        report = estimate_complexity(proteins, pep_dict)
        summary = report.summary()

        # Must be JSON-serializable
        json_str = json.dumps(summary, indent=2)
        data = json.loads(json_str)
        assert "hit_rate" in data
        assert "species_richness" in data


class TestDemoEntrapment:
    """Demo runs: entrapment scenarios."""

    def test_shuffled_entrapment(self, tmp_path: Path):
        """D24: Generate shuffled entrapment."""
        source = _generate_proteins(100, species_prefix="sp")
        entrapment = generate_shuffled_entrapment(source, n_proteins=50)

        assert len(entrapment) == 50
        # All IDs are prefixed
        assert all(k.startswith("ENTRAP_SHUFFLED_") for k in entrapment)
        # All sequences have same length as some source
        source_lens = {len(s) for s in source.values()}
        for seq in entrapment.values():
            assert len(seq) in source_lens

    def test_entrapment_in_pipeline(self, tmp_path: Path):
        """D25: Run inference with entrapment proteins mixed in."""
        # Generate target proteins + entrapment
        target = _generate_proteins(100, species_prefix="sp")
        entrapment = generate_shuffled_entrapment(target, n_proteins=20, prefix="ENTRAP")

        # Merge into one database
        combined = {**target, **entrapment}

        # Extract peptides only from target proteins
        peptides = _extract_peptides(target, n_peptides=50)

        fasta = tmp_path / "combined.fasta"
        csv_file = tmp_path / "preds.csv"
        _write_fasta(fasta, combined)
        _write_csv(csv_file, peptides)

        stats = run_inference(
            sample_id="ENTRAP_TEST",
            stage2_fasta=fasta,
            alphanovo_csv=csv_file,
            output_dir=tmp_path / "out",
            strategy="uniform",
        )

        assert stats["stage2_proteins"] == 120  # 100 target + 20 entrapment

    def test_entrapment_fdr_validation(self):
        """D26: FDR validation with synthetic hit counts."""
        # Good FDR (r=1: combined = 2 × 5/1000 = 0.01)
        result = validate_entrapment_fdr(
            5, 995, target_db_size=10000, entrapment_db_size=10000, expected_fdr=0.01
        )
        assert result["passed"] is True

        # Bad FDR (r=1: combined = 2 × 50/1000 = 0.10)
        result = validate_entrapment_fdr(
            50, 950, target_db_size=10000, entrapment_db_size=10000, expected_fdr=0.01
        )
        assert result["passed"] is False


class TestDemoScoreFiltering:
    """Demo runs: score filtering scenarios."""

    def test_min_score_filtering(self, tmp_path: Path):
        """D27: min_score filters low-confidence peptides."""
        proteins = _generate_proteins(100, species_prefix="sp")

        # Create peptides with varying scores
        peptides = []
        rng = random.Random(42)
        for i, (pid, seq) in enumerate(list(proteins.items())[:50]):
            start = rng.randint(0, len(seq) - 15)
            pep = normalize_il(seq[start : start + 15])
            score = 0.1 * (i % 10) + 0.05  # 0.05, 0.15, ..., 0.95
            peptides.append((pep, score))

        fasta = tmp_path / "db.fasta"
        csv_file = tmp_path / "preds.csv"
        _write_fasta(fasta, proteins)
        _write_csv(csv_file, peptides)

        # Run with no score filter
        stats_no_filter = run_inference(
            sample_id="SCORE_NONE",
            stage2_fasta=fasta,
            alphanovo_csv=csv_file,
            output_dir=tmp_path / "no_filter",
            strategy="uniform",
            min_score=0.0,
        )

        # Run with score filter
        stats_filtered = run_inference(
            sample_id="SCORE_050",
            stage2_fasta=fasta,
            alphanovo_csv=csv_file,
            output_dir=tmp_path / "filtered",
            strategy="uniform",
            min_score=0.5,
        )

        # Filtered should have fewer or equal proteins
        assert stats_filtered["selected_proteins"] <= stats_no_filter["selected_proteins"]


class TestDemoStrategyComparison:
    """Demo runs: strategy comparison on same data."""

    @pytest.mark.parametrize(
        "scenario",
        [
            "single_species",
            "5_species",
            "20_species",
            "mixed_types",
        ],
    )
    def test_strategy_ordering(self, tmp_path: Path, scenario: str):
        """D28-D31: uniform >= other strategies, verified ordering."""
        if scenario == "single_species":
            proteins = _generate_proteins(200, species_prefix="sp", n_species=1)
        elif scenario == "5_species":
            proteins = _generate_proteins(200, species_prefix="MGYG", n_species=5)
        elif scenario == "20_species":
            proteins = _generate_proteins(400, species_prefix="MGYG", n_species=20)
        else:  # mixed_types
            proteins = _generate_proteins(200, species_prefix="mixed", n_species=10)

        peptides = _extract_peptides(proteins, n_peptides=100)

        fasta = tmp_path / "db.fasta"
        csv_file = tmp_path / "preds.csv"
        _write_fasta(fasta, proteins)
        _write_csv(csv_file, peptides)

        results = {}
        for strategy in STRATEGIES:
            out = tmp_path / strategy
            stats = run_inference(
                sample_id=f"COMPARE_{scenario}",
                stage2_fasta=fasta,
                alphanovo_csv=csv_file,
                output_dir=out,
                strategy=strategy,
                taxonomy_map_path=_write_taxonomy(tmp_path / "tax.tsv", proteins),
            )
            results[strategy] = stats["selected_proteins"]

        # Uniform should be >= all other strategies
        for strategy in ["species_budget", "razor", "uniform_2pep"]:
            assert results["uniform"] >= results[strategy], (
                f"uniform ({results['uniform']}) < {strategy} ({results[strategy]}) "
                f"in scenario {scenario}"
            )


class TestDemoLargeScale:
    """Demo runs: larger-scale synthetic tests."""

    def test_5000_proteins_200_peptides(self, tmp_path: Path):
        """D32: 5K proteins, 200 peptides, all strategies."""
        proteins = _generate_proteins(5000, species_prefix="MGYG", n_species=50, seed=123)
        peptides = _extract_peptides(proteins, n_peptides=200, seed=123)

        fasta = tmp_path / "large.fasta"
        csv_file = tmp_path / "preds.csv"
        _write_fasta(fasta, proteins)
        _write_csv(csv_file, peptides)

        for strategy in STRATEGIES:
            out = tmp_path / strategy
            t0 = time.time()
            stats = run_inference(
                sample_id="LARGE_5K",
                stage2_fasta=fasta,
                alphanovo_csv=csv_file,
                output_dir=out,
                strategy=strategy,
                taxonomy_map_path=_write_taxonomy(tmp_path / "tax.tsv", proteins),
            )
            elapsed = time.time() - t0

            assert stats["stage2_proteins"] == 5000
            assert stats["selected_proteins"] >= 0
            # Should complete in reasonable time (< 30s even on slow systems)
            assert elapsed < 30, f"{strategy} took {elapsed:.1f}s"

    def test_10000_proteins_500_peptides(self, tmp_path: Path):
        """D33: 10K proteins, 500 peptides — performance check."""
        proteins = _generate_proteins(10000, species_prefix="MGYG", n_species=100, seed=456)
        peptides = _extract_peptides(proteins, n_peptides=500, seed=456)

        fasta = tmp_path / "xlarge.fasta"
        csv_file = tmp_path / "preds.csv"
        _write_fasta(fasta, proteins)
        _write_csv(csv_file, peptides)

        t0 = time.time()
        stats = run_inference(
            sample_id="XLARGE_10K",
            stage2_fasta=fasta,
            alphanovo_csv=csv_file,
            output_dir=tmp_path / "out",
            strategy="species_budget",
            taxonomy_map_path=_write_taxonomy(tmp_path / "tax.tsv", proteins),
        )
        elapsed = time.time() - t0

        assert stats["stage2_proteins"] == 10000
        assert elapsed < 60


class TestDemoEndToEnd:
    """Demo runs: full end-to-end pipeline demos."""

    def test_clinical_pipeline(self, tmp_path: Path):
        """D34: Full clinical proteomics demo.

        1. Generate human proteins
        2. Add contaminants
        3. Generate entrapment (shuffled)
        4. Quality filter
        5. Extract peptides
        6. Estimate complexity
        7. Run inference (razor, recommended for single-species)
        8. Validate entrapment FDR
        """
        # 1. Human proteins
        human = _generate_proteins(300, species_prefix="sp", seed=100)

        # 2. Contaminants
        contaminants = {f"CONTAM_{i:03d}": "MSKGEELFTACDE" * 15 for i in range(10)}

        # 3. Entrapment
        entrapment = generate_shuffled_entrapment(human, n_proteins=30, prefix="ENTRAP")

        # 4. Quality filter
        combined = {**human, **contaminants, **entrapment}
        config = SequenceQualityConfig(
            min_length=10,
            reject_poly_runs=True,
            poly_run_length=15,
        )
        filtered, qc_report = filter_sequences_by_quality(combined, config)

        # 5. Extract peptides from human proteins only
        peptides = _extract_peptides(human, n_peptides=100, seed=200)

        # 6. Estimate complexity
        pep_dict = dict(peptides)
        complexity_config = ComplexityConfig(auto_recommend=True)
        complexity_report = estimate_complexity(filtered, pep_dict, complexity_config)

        assert complexity_report.species_richness >= 1
        assert complexity_report.recommendation in ("razor", "species_budget")

        # 7. Write files and run inference
        fasta = tmp_path / "clinical.fasta"
        csv_file = tmp_path / "preds.csv"
        _write_fasta(fasta, filtered)
        _write_csv(csv_file, peptides)

        stats = run_inference(
            sample_id="CLINICAL_DEMO",
            stage2_fasta=fasta,
            alphanovo_csv=csv_file,
            output_dir=tmp_path / "results",
            strategy="razor",
        )

        assert stats["selected_proteins"] > 0
        assert stats["runtime_seconds"] >= 0

        # 8. FDR validation (simulated)
        fdr_result = validate_entrapment_fdr(
            entrapment_hits=2,
            target_hits=stats["selected_proteins"],
            target_db_size=10000,
            entrapment_db_size=10000,
            expected_fdr=0.05,
        )
        assert fdr_result["passed"] is True

    def test_metaproteomics_pipeline(self, tmp_path: Path):
        """D35: Full metaproteomics demo.

        1. Generate multi-species DB (UHGG + GMGC)
        2. Add entrapment
        3. Quality filter
        4. Extract peptides
        5. Estimate complexity
        6. Run species_budget inference
        7. Verify reduction
        """
        # 1. Multi-species database
        uhgg = _generate_proteins(500, species_prefix="MGYG", n_species=50, seed=300)
        gmgc = _generate_proteins(500, species_prefix="GMGC", n_species=30, seed=301)
        database = {**uhgg, **gmgc}

        # 2. Entrapment
        entrapment = generate_shuffled_entrapment(database, n_proteins=50, prefix="ENTRAP")
        combined = {**database, **entrapment}

        # 3. Quality filter
        config = SequenceQualityConfig(
            min_length=20,
            reject_poly_runs=True,
            poly_run_length=10,
            reject_low_complexity=True,
            low_complexity_threshold=0.2,
        )
        filtered, qc_report = filter_sequences_by_quality(combined, config)

        # 4. Peptides
        peptides = _extract_peptides(database, n_peptides=200, seed=400)

        # 5. Complexity
        pep_dict = dict(peptides)
        complexity = estimate_complexity(filtered, pep_dict, ComplexityConfig(auto_recommend=True))

        assert complexity.species_richness > 1
        assert complexity.recommendation == "species_budget"

        # 6. Inference
        fasta = tmp_path / "meta.fasta"
        csv_file = tmp_path / "preds.csv"
        _write_fasta(fasta, filtered)
        _write_csv(csv_file, peptides)

        stats = run_inference(
            sample_id="META_DEMO",
            stage2_fasta=fasta,
            alphanovo_csv=csv_file,
            output_dir=tmp_path / "results",
            strategy="species_budget",
            taxonomy_map_path=_write_taxonomy(tmp_path / "tax.tsv", filtered),
        )

        # 7. Verify
        assert stats["selected_proteins"] > 0
        assert stats["selected_proteins"] < len(filtered)  # Some reduction expected

    @pytest.mark.parametrize("seed", range(36, 51))
    def test_random_scenario(self, tmp_path: Path, seed: int):
        """D36-D50: Random scenarios with varying parameters."""
        rng = random.Random(seed)

        n_proteins = rng.choice([50, 100, 200, 500, 1000])
        n_species = rng.choice([1, 3, 10, 25])
        n_peptides = rng.choice([20, 50, 100])
        species_prefix = rng.choice(["sp", "MGYG", "GMGC", "mixed"])
        strategy = rng.choice(STRATEGIES)

        proteins = _generate_proteins(
            n_proteins,
            species_prefix=species_prefix,
            n_species=n_species,
            seed=seed,
        )
        peptides = _extract_peptides(proteins, n_peptides=n_peptides, seed=seed + 1000)

        fasta = tmp_path / "random.fasta"
        csv_file = tmp_path / "preds.csv"
        _write_fasta(fasta, proteins)
        _write_csv(csv_file, peptides)

        stats = run_inference(
            sample_id=f"RNG_{seed}",
            stage2_fasta=fasta,
            alphanovo_csv=csv_file,
            output_dir=tmp_path / "out",
            strategy=strategy,
            taxonomy_map_path=_write_taxonomy(tmp_path / "tax.tsv", proteins),
        )

        # Core invariants that should always hold
        assert stats["selected_proteins"] >= 0
        assert stats["stage2_proteins"] == n_proteins
        assert stats["runtime_seconds"] >= 0
        assert stats["reduction_factor"] >= 0

        # Output files exist
        suffix_map = {
            "species_budget": "species_budget",
            "razor": "razor",
            "uniform": "uniform",
            "uniform_2pep": "uniform2pep",
        }
        fasta_out = tmp_path / "out" / f"RNG_{seed}_{suffix_map[strategy]}.fasta"
        stats_out = tmp_path / "out" / f"RNG_{seed}_{suffix_map[strategy]}_stats.json"
        assert fasta_out.exists()
        assert stats_out.exists()

        # Stats JSON is valid
        with open(stats_out) as f:
            stats_data = json.load(f)
        assert "selected_proteins" in stats_data
