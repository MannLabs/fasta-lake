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


def test_run_and_group_study_are_installed_commands():
    """The main workflow is reachable as ``fasta-lake run``, not only as a repository script."""
    import subprocess
    import sys

    for command, flag in (("run", "--manifest"), ("group-study", "--search")):
        result = subprocess.run(
            [sys.executable, "-m", "fasta_lake.cli", command, "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        assert f"fasta-lake {command}" in result.stdout
        assert flag in result.stdout
    listing = CliRunner().invoke(cli, ["--help"])
    assert "run" in listing.output and "group-study" in listing.output


def test_repository_scripts_remain_compatibility_wrappers():
    import subprocess
    import sys

    root = Path(__file__).resolve().parents[1]
    for script in ("tools/run_manifest.py", "tools/group_study.py"):
        result = subprocess.run(
            [sys.executable, str(root / script), "--help"], capture_output=True, text=True
        )
        assert result.returncode == 0, result.stderr


def test_member_set_aggregate_is_opt_in_and_needs_searches(tmp_path, capsys):
    """The default run reports one protein-group matrix; member sets are a diagnostic."""
    import pytest

    from fasta_lake import workflow

    base = ["--manifest", str(tmp_path / "m.tsv"), "--out", str(tmp_path / "out")]
    with pytest.raises(SystemExit) as error:
        workflow.main(base + ["--databases-only", "--member-set-aggregate"])
    assert error.value.code == 2
    assert "member-set-aggregate requires searches" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        workflow.main(["--help"])
    assert "Diagnostic only" in capsys.readouterr().out
    assert not (tmp_path / "out").exists()


LEGACY = ("infer", "convert", "tune", "diagnose", "preset", "lake-build")
INTERNAL = (
    "DESIGN_V5",
    "PREFIX_CENSUS",
    "razor_tiebreak_bias",
    "COMPARABILITY_RESOLUTION",
    "MicrobPredict",
    "PREDICT_<",
    "pre-v26",
    "For v26",
)


def test_legacy_commands_are_hidden_but_still_run_with_a_notice():
    runner = CliRunner()
    listing = runner.invoke(cli, ["--help"]).output
    commands = {line.split()[0] for line in listing.split("Commands:")[1].splitlines() if line}
    assert not commands & set(LEGACY)
    assert {"run", "group-study"} <= commands
    result = runner.invoke(cli, ["preset", "list"])
    assert result.exit_code == 0
    assert "legacy route" in result.stderr and "fasta-lake run" in result.stderr
    assert "legacy route" not in result.stdout


def test_visible_help_has_no_internal_references():
    """Help text must not point users at private files, cohorts or internal versions."""
    runner = CliRunner()

    def walk(group, prefix):
        for name, command in group.commands.items():
            if command.hidden:
                continue
            text = runner.invoke(cli, [*prefix, name, "--help"]).output
            yield " ".join([*prefix, name]), text
            if isinstance(command, click.Group):
                yield from walk(command, [*prefix, name])

    import click

    texts = dict(walk(cli, []))
    assert "run" in texts and "entrapment fdp" in texts
    for name, text in texts.items():
        for marker in INTERNAL:
            assert marker not in text, f"{name} --help mentions {marker}"


def test_rust_help_has_no_internal_references():
    import os
    import subprocess

    import pytest

    folder = os.environ.get("FASTALAKE_BIN_DIR")
    if not folder:
        pytest.skip("Set FASTALAKE_BIN_DIR to the built Rust executables")
    for name in ("fasta_extractor", "parsimony_engine", "aggregation_engine", "lake_builder"):
        text = subprocess.run(
            [str(Path(folder) / name), "--help"], capture_output=True, text=True
        ).stdout
        assert text, name
        for marker in INTERNAL + ("Measured on PREDICT", "CAPSCAN"):
            assert marker not in text, f"{name} --help mentions {marker}"
