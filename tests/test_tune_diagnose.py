"""Unit tests for fasta_lake.tune and fasta_lake.diagnose (v5)."""

from pathlib import Path

from fasta_lake.diagnose import (
    CALIBRATION_Q,
    calibration_curve,
    outlier_samples,
    per_sample_fdp_at_q,
)
from fasta_lake.tune import (
    FilterSpec,
    apply_filter,
    compute_fdp_combined,
    default_filter_grid,
    tune,
    write_reports,
)


def _mk_fake_tsv(path: Path, psms: list[tuple[float, int, float, str]]) -> None:
    """Write a minimal SAGE TSV for testing.
    psms: list of (q, len, score, proteins)
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("spectrum_q\tpeptide_len\tsage_discriminant_score\tproteins\n")
        for q, L, s, p in psms:
            f.write(f"{q}\t{L}\t{s}\t{p}\n")


def _mk_dataset(base: Path, n_samples: int = 3) -> Path:
    """Create a fake results dir with n_samples × 3 entrap types."""
    base.mkdir(parents=True, exist_ok=True)
    for i in range(n_samples):
        sample_dir = base / f"sample_{i:03d}"
        # shuffled: 100 targets @ q<=0.01, 2 entrap @ q<=0.01
        _mk_fake_tsv(
            sample_dir / "entrapment" / "results.sage.tsv",
            [
                *[(0.005, 14, 0.9, "PROT_x")] * 100,
                *[(0.008, 8, 0.6, "SHUF_y")] * 2,
                *[(0.03, 12, 0.4, "PROT_z")] * 50,  # won't pass q<=0.01
            ],
        )
        _mk_fake_tsv(
            sample_dir / "entrapment_plant" / "results.sage.tsv",
            [
                *[(0.005, 14, 0.9, "PROT_x")] * 100,
                *[(0.008, 9, 0.6, "PLANT_y")] * 5,
            ],
        )
        _mk_fake_tsv(
            sample_dir / "entrapment_archaea" / "results.sage.tsv",
            [
                *[(0.005, 14, 0.9, "PROT_x")] * 100,
            ],
        )
    return base


class TestTune:
    def test_filter_grid_default(self):
        grid = default_filter_grid()
        assert len(grid) == 3 * 4 * 3
        assert all(isinstance(s, FilterSpec) for s in grid)

    def test_apply_filter(self):
        rows = [
            {"q": 0.001, "len": 10, "score": 0.9, "proteins": "PROT_a", "is_decoy": False},
            {"q": 0.005, "len": 8, "score": 0.5, "proteins": "SHUF_b", "is_decoy": False},
            {"q": 0.02, "len": 12, "score": 0.8, "proteins": "PROT_c", "is_decoy": False},
            {"q": 0.005, "len": 10, "score": 0.9, "proteins": "rev_d", "is_decoy": True},
        ]
        spec = FilterSpec("t", q_threshold=0.01, min_length=9, min_score=0.0)
        c = apply_filter(rows, spec, entrap_prefix="SHUF_")
        # len>=9 drops the SHUF one (len=8). q<=0.01 drops the len=12 one.
        assert c["target"] == 1
        assert c["decoy"] == 1
        assert c["entrap"] == 0

    def test_compute_fdp_combined(self):
        # 2 entrap, 100 target+decoy, r=0.1 → 2 × (1+1/0.1) / 102 = 2 × 11 / 102 ≈ 21.57%
        fdp = compute_fdp_combined(2, 100, 0.1)
        assert 21.0 < fdp < 22.0

    def test_tune_full(self, tmp_path):
        results_dir = _mk_dataset(tmp_path / "results", n_samples=5)
        report = tune(
            results_dir,
            target_fdp_pct=5.0,
            target_db_size_fn=lambda _sid: 75000,
        )
        assert report["n_samples"] == 5
        assert len(report["all_results"]) > 0
        # Winner must exist and satisfy threshold
        if report["winner"]:
            assert report["winner"]["median_fdp_pct"] <= 5.0

    def test_write_reports(self, tmp_path):
        results_dir = _mk_dataset(tmp_path / "r", n_samples=3)
        report = tune(results_dir, target_fdp_pct=5.0)
        out = tmp_path / "out"
        write_reports(report, out)
        assert (out / "tuning_report.csv").exists()
        assert (out / "recommended.json").exists()


class TestDiagnose:
    def test_per_sample_fdp_at_q(self):
        rows = [
            {"q": 0.001, "len": 10, "score": 0.9, "proteins": "PROT_a", "is_decoy": False},
            {"q": 0.005, "len": 9, "score": 0.8, "proteins": "SHUF_b", "is_decoy": False},
        ]
        c = per_sample_fdp_at_q(rows, "SHUF_", q_thr=0.01, min_length=9)
        assert c["target"] == 1
        assert c["entrap"] == 1
        # v5.1: also returns protein-level counts
        assert c["prot_target"] == 1
        assert c["prot_entrap"] == 1
        assert c["contaminant"] == 0

    def test_per_sample_fdp_protein_level(self):
        # 3 PSMs all mapping to PROT_a → 1 protein (picked)
        rows = [
            {"q": 0.001, "len": 10, "score": 0.9, "proteins": "PROT_a", "is_decoy": False},
            {"q": 0.002, "len": 10, "score": 0.8, "proteins": "PROT_a", "is_decoy": False},
            {"q": 0.003, "len": 10, "score": 0.7, "proteins": "PROT_a", "is_decoy": False},
        ]
        c = per_sample_fdp_at_q(rows, "SHUF_", q_thr=0.01, min_length=9)
        assert c["target"] == 3
        assert c["prot_target"] == 1  # picked down to one protein

    def test_per_sample_excludes_contaminants(self):
        # Trypsin P00761 is in the default contaminant list
        rows = [
            {"q": 0.001, "len": 10, "score": 0.9, "proteins": "sp|P04637|P53", "is_decoy": False},
            {"q": 0.005, "len": 10, "score": 0.8, "proteins": "sp|P00761|TRYP", "is_decoy": False},
        ]
        c = per_sample_fdp_at_q(rows, "SHUF_", q_thr=0.01, min_length=9, exclude_contaminants=True)
        assert c["target"] == 1  # p53 counted
        assert c["contaminant"] == 1  # trypsin excluded from target

    def test_calibration_curve(self, tmp_path):
        results_dir = _mk_dataset(tmp_path / "results", n_samples=5)
        cal = calibration_curve(results_dir, min_length=9)
        assert cal["n_samples"] == 5
        entrap_types = set(p["entrap_type"] for p in cal["curve"])
        assert entrap_types == {"shuffled", "plant", "archaea"}
        # Curve has 9 q-values × 3 entrap types = 27 points max
        assert len(cal["curve"]) <= len(CALIBRATION_Q) * 3

    def test_outliers_empty_when_small_dataset(self, tmp_path):
        results_dir = _mk_dataset(tmp_path / "r", n_samples=3)
        out = outlier_samples(results_dir, fdp_z_threshold=2.0, min_length=9)
        # With 3 identical samples, no outliers
        assert out == []
