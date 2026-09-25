"""T14: Rust binary integration tests."""

from __future__ import annotations

from pathlib import Path

import pytest

from fasta_lake.rust import binaries
from fasta_lake.rust.binaries import find_binary, run_fasta_extractor

REPO_ROOT = Path(__file__).parent.parent


def _capture_extractor_cmd(monkeypatch, tmp_path, **kwargs):
    """Run run_fasta_extractor with find_binary/subprocess stubbed, return argv."""
    captured = {}

    def fake_find_binary(name, search_dir=None, explicit_path=None):
        return Path("/fake/fasta_extractor")

    class FakeResult:
        stdout = ""

    def fake_run(cmd, **_):
        captured["cmd"] = cmd
        return FakeResult()

    monkeypatch.setattr(binaries, "find_binary", fake_find_binary)
    monkeypatch.setattr(binaries.subprocess, "run", fake_run)

    (tmp_path / "db.fasta").write_text(">P1\nPEPTIDEAK\n")
    (tmp_path / "pred.csv").write_text("sequence,score\nPEPTIDEAK,0.9\n")
    run_fasta_extractor(
        database=tmp_path / "db.fasta",
        predictions=tmp_path / "pred.csv",
        output=tmp_path / "out.fasta",
        **kwargs,
    )
    return captured["cmd"]


class TestNormalizeIlFlag:
    """Regression: rebuilt v2 binary uses clap ArgAction::Set, so --normalize-il
    REQUIRES a value. A bare --normalize-il exits 2 (broke jC1_stage2, job 5340153)."""

    def test_normalize_il_true_passes_value(self, monkeypatch, tmp_path):
        cmd = _capture_extractor_cmd(monkeypatch, tmp_path, normalize_il=True)
        assert "--normalize-il" in cmd
        idx = cmd.index("--normalize-il")
        # Must be followed by an explicit value, never left bare.
        assert idx + 1 < len(cmd), "--normalize-il passed bare (would exit 2)"
        assert cmd[idx + 1] == "true"

    def test_normalize_il_false_passes_value(self, monkeypatch, tmp_path):
        cmd = _capture_extractor_cmd(monkeypatch, tmp_path, normalize_il=False)
        idx = cmd.index("--normalize-il")
        assert cmd[idx + 1] == "false"

    def test_normalize_il_never_bare(self, monkeypatch, tmp_path):
        """The value must not be another flag (i.e. the flag is never bare)."""
        for val in (True, False):
            cmd = _capture_extractor_cmd(monkeypatch, tmp_path, normalize_il=val)
            idx = cmd.index("--normalize-il")
            assert not cmd[idx + 1].startswith("--")


class TestRustBinary:
    """T14: find_binary() locates fasta_extractor."""

    def test_find_fasta_extractor(self):
        """Test that fasta_extractor binary can be found."""
        try:
            binary = find_binary("fasta_extractor", search_dir=REPO_ROOT)
            assert binary.exists()
            assert binary.is_file()
        except FileNotFoundError:
            pytest.skip(
                "Rust binary not compiled. Run: cd rust/fasta_extractor_v2 && cargo build --release"
            )

    def test_find_lake_builder(self):
        """Test that lake_builder binary can be found."""
        try:
            binary = find_binary("lake_builder", search_dir=REPO_ROOT)
            assert binary.exists()
            assert binary.is_file()
        except FileNotFoundError:
            pytest.skip(
                "Rust binary not compiled. Run: cd rust/fasta_extractor_v2 && cargo build --release"
            )

    def test_explicit_path_not_found(self):
        """Test that explicit nonexistent path raises."""
        with pytest.raises(FileNotFoundError):
            find_binary("fasta_extractor", explicit_path="/nonexistent/binary")

    def test_binary_not_found_message(self, monkeypatch):
        """Test helpful error message when binary not found."""
        monkeypatch.delenv("FASTALAKE_BIN_DIR", raising=False)
        monkeypatch.setattr(binaries.shutil, "which", lambda _: None)
        with pytest.raises(FileNotFoundError, match="cargo build"):
            find_binary("fasta_extractor", search_dir="/tmp/empty_dir_12345")
