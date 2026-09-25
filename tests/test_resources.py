"""Stage measurements retain output and failures from real subprocesses."""

import sys

import pytest

from fasta_lake.resources import run_recorded


@pytest.mark.parametrize("exit_code", [0, 7])
def test_command_status_output_and_resource_units(tmp_path, exit_code):
    path = tmp_path / "stdout"
    with path.open("w") as stream:
        result = run_recorded(
            [
                sys.executable,
                "-c",
                f"import sys; data=bytearray(8*1024*1024); print(len(data)); sys.exit({exit_code})",
            ],
            stdout=stream,
        )
    assert result["exit_code"] == exit_code
    assert path.read_text().strip() == str(8 * 1024 * 1024)
    assert result["seconds"] > 0
    if sys.platform in {"linux", "darwin"}:
        assert result["max_process_rss_bytes"] >= 8 * 1024 * 1024
        assert result["user_cpu_seconds"] >= 0
