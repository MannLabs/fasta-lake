"""T13: CLI smoke tests."""

from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from fasta_lake.cli import cli


class TestCLISmoke:
    """T13: Click CliRunner invokes fasta-lake commands → exit 0."""

    def test_help(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "fasta-lake" in result.output.lower() or "metaproteomics" in result.output.lower()

    def test_version(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["--version"])
        assert result.exit_code == 0
        assert "fasta-lake, version" in result.output

    def test_infer_help(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["infer", "--help"])
        assert result.exit_code == 0
        assert "strategy" in result.output.lower()

    def test_infer_runs(self, tiny_fasta: Path, tiny_predictions: Path, tmp_output: Path):
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "infer",
                "--sample",
                "CLI_TEST",
                "--stage2-fasta",
                str(tiny_fasta),
                "--peptides",
                str(tiny_predictions),
                "-o",
                str(tmp_output),
                "--strategy",
                "uniform",
            ],
        )
        assert result.exit_code == 0, f"CLI failed: {result.output}"
        assert "Selected proteins" in result.output

    def test_validate_command(self, tmp_path: Path):
        config = {
            "database": {"sources": [{"tag": "X", "path": "/fake"}]},
            "peptide_evidence": {"predictions_dir": "/fake"},
            "output": {"directory": "/tmp/out"},
        }
        config_file = tmp_path / "test.json"
        config_file.write_text(json.dumps(config))

        runner = CliRunner()
        result = runner.invoke(cli, ["validate", str(config_file)])
        assert result.exit_code == 0
        assert "valid" in result.output.lower()

    def test_init_command(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["init"])
        assert result.exit_code == 0
        # Should output valid JSON
        data = json.loads(result.output)
        assert "database" in data

    def test_init_minimal(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["init", "--template", "minimal"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "database" in data
