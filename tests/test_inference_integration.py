"""T11-T12: Integration tests for inference with file I/O."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from fasta_lake.inference.runner import run_inference

# ── T11: Full inference with file I/O ────────────────────────────────────────


class TestInferenceWithFileIO:
    """T11: Full run_inference() with tiny.fasta → output FASTA + stats JSON."""

    def test_produces_fasta_and_stats(
        self, tiny_fasta: Path, tiny_predictions: Path, tiny_taxonomy: Path, tmp_output: Path
    ):
        run_inference(
            sample_id="TEST_SAMPLE",
            stage2_fasta=tiny_fasta,
            alphanovo_csv=tiny_predictions,
            output_dir=tmp_output,
            strategy="species_budget",
            taxonomy_map_path=tiny_taxonomy,
        )

        # Check output FASTA exists
        fasta_out = tmp_output / "TEST_SAMPLE_species_budget.fasta"
        assert fasta_out.exists()
        assert fasta_out.stat().st_size > 0

        # Check stats JSON exists
        stats_out = tmp_output / "TEST_SAMPLE_species_budget_stats.json"
        assert stats_out.exists()

        with open(stats_out) as f:
            saved_stats = json.load(f)

        assert saved_stats["sample_id"] == "TEST_SAMPLE"
        assert saved_stats["strategy"] == "species_budget"
        assert saved_stats["selected_proteins"] > 0
        assert saved_stats["stage2_proteins"] == 20

    def test_output_fasta_valid(
        self, tiny_fasta: Path, tiny_predictions: Path, tiny_taxonomy: Path, tmp_output: Path
    ):
        run_inference(
            sample_id="TEST_SAMPLE",
            stage2_fasta=tiny_fasta,
            alphanovo_csv=tiny_predictions,
            output_dir=tmp_output,
            strategy="species_budget",
            taxonomy_map_path=tiny_taxonomy,
        )

        fasta_out = tmp_output / "TEST_SAMPLE_species_budget.fasta"
        # Verify it's valid FASTA
        with open(fasta_out) as f:
            lines = f.readlines()

        headers = [line for line in lines if line.startswith(">")]
        assert len(headers) > 0
        # All non-header lines should be sequence
        for line in lines:
            if not line.startswith(">"):
                assert line.strip().isalpha()

    def test_species_budget_refuses_unresolved_species(
        self, tiny_fasta: Path, tiny_predictions: Path, tmp_output: Path
    ):
        """Without a taxonomy map the smORF/assembly ids have no species: refuse, name them."""
        with pytest.raises(
            ValueError,
            match=r"9 of 14 proteins .* no resolvable species \(e\.g\. GMGC10\.001_001_001",
        ):
            run_inference(
                sample_id="TEST_SAMPLE",
                stage2_fasta=tiny_fasta,
                alphanovo_csv=tiny_predictions,
                output_dir=tmp_output,
                strategy="species_budget",
            )

    def test_missing_fasta_raises(self, tiny_predictions: Path, tmp_output: Path):
        with pytest.raises(FileNotFoundError):
            run_inference(
                sample_id="X",
                stage2_fasta="/nonexistent.fasta",
                alphanovo_csv=tiny_predictions,
                output_dir=tmp_output,
            )

    def test_invalid_strategy_raises(
        self, tiny_fasta: Path, tiny_predictions: Path, tmp_output: Path
    ):
        with pytest.raises(ValueError, match="Unknown strategy"):
            run_inference(
                sample_id="X",
                stage2_fasta=tiny_fasta,
                alphanovo_csv=tiny_predictions,
                output_dir=tmp_output,
                strategy="nonexistent",
            )


# ── T12: All strategies produce output ──────────────────────────────────────


class TestAllStrategiesProduceOutput:
    """T12: All 4 strategies → valid outputs, uniform >= species_budget."""

    @pytest.fixture
    def all_results(
        self, tiny_fasta: Path, tiny_predictions: Path, tiny_taxonomy: Path, tmp_output: Path
    ) -> dict[str, dict]:
        results = {}
        for strategy in ("species_budget", "razor", "uniform", "uniform_2pep"):
            out_dir = tmp_output / strategy
            out_dir.mkdir()
            stats = run_inference(
                sample_id="TEST",
                stage2_fasta=tiny_fasta,
                alphanovo_csv=tiny_predictions,
                output_dir=out_dir,
                strategy=strategy,
                taxonomy_map_path=tiny_taxonomy,
            )
            results[strategy] = stats
        return results

    def test_all_produce_output(self, all_results: dict):
        for strategy, stats in all_results.items():
            assert stats["selected_proteins"] > 0, f"{strategy} selected 0 proteins"

    def test_uniform_is_superset(self, all_results: dict):
        # Uniform keeps everything; others should be <= uniform
        uniform_count = all_results["uniform"]["selected_proteins"]
        for strategy in ("species_budget", "razor", "uniform_2pep"):
            assert all_results[strategy]["selected_proteins"] <= uniform_count, (
                f"{strategy} ({all_results[strategy]['selected_proteins']}) "
                f"> uniform ({uniform_count})"
            )

    def test_stats_have_required_keys(self, all_results: dict):
        required = {"sample_id", "strategy", "selected_proteins", "stage2_proteins"}
        for strategy, stats in all_results.items():
            assert required.issubset(stats.keys()), (
                f"{strategy} missing keys: {required - stats.keys()}"
            )


# ── k-mer rescue: failures are hard errors, zero is recorded ─────────────────


class TestKmerRescueFailureIsLoud:
    """A requested rescue that cannot run must stop the run, not shrink the database."""

    @staticmethod
    def _predictions_with_failed_peptide(tiny_predictions: Path, tmp_path: Path) -> Path:
        p = tmp_path / "preds.csv"
        p.write_text(tiny_predictions.read_text() + "WWWWWWWWWWWWK,0.99\n")
        return p

    @staticmethod
    def _rescue_db(tmp_path: Path) -> Path:
        db = tmp_path / "rescue.fasta"
        db.write_text(">RESCUED_1 from lake\nMWWWWWWWWWWWWKLLL\n>OTHER_1\nMAAAAK\n")
        return db

    @staticmethod
    def _fake_binary(tmp_path: Path, body: str) -> Path:
        script = tmp_path / "fake_mass_kmer_anchor.sh"
        script.write_text("#!/bin/bash\n" + body)
        script.chmod(0o755)
        return script

    def _run(self, tiny_fasta, preds, tmp_output, db):
        return run_inference(
            sample_id="S",
            stage2_fasta=tiny_fasta,
            alphanovo_csv=preds,
            output_dir=tmp_output,
            strategy="razor",
            kmer_rescue_db=db,
        )

    def test_missing_binary_raises(
        self, tiny_fasta, tiny_predictions, tmp_output, tmp_path, monkeypatch
    ):
        preds = self._predictions_with_failed_peptide(tiny_predictions, tmp_path)
        monkeypatch.setenv("KMER_BINARY", str(tmp_path / "no_such_binary"))
        with pytest.raises(FileNotFoundError, match="mass_kmer_anchor not found"):
            self._run(tiny_fasta, preds, tmp_output, self._rescue_db(tmp_path))

    def test_missing_rescue_db_raises(self, tiny_fasta, tiny_predictions, tmp_output, tmp_path):
        with pytest.raises(FileNotFoundError, match="rescue DB not found"):
            self._run(tiny_fasta, tiny_predictions, tmp_output, tmp_path / "absent.fasta")

    def test_nonzero_exit_raises(
        self, tiny_fasta, tiny_predictions, tmp_output, tmp_path, monkeypatch
    ):
        preds = self._predictions_with_failed_peptide(tiny_predictions, tmp_path)
        fake = self._fake_binary(tmp_path, "echo 'index exploded' >&2\nexit 3\n")
        monkeypatch.setenv("KMER_BINARY", str(fake))
        with pytest.raises(RuntimeError, match=r"failed \(exit 3\): index exploded"):
            self._run(tiny_fasta, preds, tmp_output, self._rescue_db(tmp_path))
        assert not list(tmp_output.glob("*_stats.json")), "no output on a failed rescue"

    def test_successful_rescue_is_injected_and_recorded(
        self, tiny_fasta, tiny_predictions, tmp_output, tmp_path, monkeypatch
    ):
        preds = self._predictions_with_failed_peptide(tiny_predictions, tmp_path)
        # The fake binary writes one hit to whatever --output it is given.
        fake = self._fake_binary(
            tmp_path,
            'while [ $# -gt 0 ]; do [ "$1" = "--output" ] && OUT="$2"; shift; done\n'
            'printf "query_sequence\\tprotein_id\\nWWWWWWWWWWWWK\\tRESCUED_1\\n" > "$OUT"\n',
        )
        monkeypatch.setenv("KMER_BINARY", str(fake))
        stats = self._run(tiny_fasta, preds, tmp_output, self._rescue_db(tmp_path))
        assert stats["kmer_proteins_rescued"] == 1
        assert stats["kmer_proteins_new"] == 1
        assert "RESCUED_1" in (tmp_output / "S_razor.fasta").read_text()

    def test_zero_rescue_is_recorded_not_silent(
        self, tiny_fasta, tiny_predictions, tmp_output, tmp_path, monkeypatch
    ):
        preds = self._predictions_with_failed_peptide(tiny_predictions, tmp_path)
        fake = self._fake_binary(tmp_path, "exit 0\n")  # exits 0, writes nothing
        monkeypatch.setenv("KMER_BINARY", str(fake))
        stats = self._run(tiny_fasta, preds, tmp_output, self._rescue_db(tmp_path))
        assert stats["kmer_proteins_rescued"] == 0
        saved = json.loads((tmp_output / "S_razor_stats.json").read_text())
        assert saved["kmer_rescue_db"].endswith("rescue.fasta")
        assert saved["kmer_proteins_rescued"] == 0
