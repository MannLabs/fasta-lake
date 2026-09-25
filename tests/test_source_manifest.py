"""tools/source_manifest.py must describe exactly the tracked tree and notice drift."""

from __future__ import annotations

import importlib.util
import json
import subprocess
from pathlib import Path

import pytest

TOOL = Path(__file__).resolve().parents[1] / "tools" / "source_manifest.py"


def _load():
    spec = importlib.util.spec_from_file_location("source_manifest", TOOL)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _git(root: Path, *args: str) -> None:
    subprocess.run(
        ["git", "-C", str(root), *args],
        check=True,
        capture_output=True,
        env={
            "GIT_AUTHOR_NAME": "t",
            "GIT_AUTHOR_EMAIL": "t@t",
            "GIT_COMMITTER_NAME": "t",
            "GIT_COMMITTER_EMAIL": "t@t",
            "HOME": str(root),
            "PATH": "/usr/bin:/bin",
        },
    )


def _repo(tmp_path: Path) -> Path:
    root = tmp_path / "software"
    (root / "pkg").mkdir(parents=True)
    (root / "pkg" / "a.py").write_text("print('a')\n")
    (root / "README.md").write_text("hello\n")
    (root / "untracked.tmp").write_text("not in git\n")
    _git(root, "init", "-q")
    _git(root, "add", "pkg/a.py", "README.md")
    _git(root, "commit", "-q", "-m", "init")
    return root


def test_write_lists_tracked_files_only_and_check_passes(tmp_path):
    sm = _load()
    root = _repo(tmp_path)
    manifest = root / "SOURCE_MANIFEST.json"
    assert sm.write(root, manifest) == 0
    listed = [f["path"] for f in json.loads(manifest.read_text())["files"]]
    assert listed == ["README.md", "pkg/a.py"], listed  # sorted, no untracked, not itself
    assert sm.check(root, manifest) == 0


def test_check_fails_on_content_change_and_on_new_tracked_file(tmp_path, capsys):
    sm = _load()
    root = _repo(tmp_path)
    manifest = root / "SOURCE_MANIFEST.json"
    sm.write(root, manifest)
    (root / "pkg" / "a.py").write_text("print('b')\n")
    assert sm.check(root, manifest) == 1
    assert "content differs: pkg/a.py" in capsys.readouterr().out
    (root / "pkg" / "a.py").write_text("print('a')\n")
    (root / "pkg" / "new.py").write_text("\n")
    _git(root, "add", "pkg/new.py")
    assert sm.check(root, manifest) == 1
    assert "tracked but not listed: pkg/new.py" in capsys.readouterr().out


def test_write_is_deterministic(tmp_path):
    sm = _load()
    root = _repo(tmp_path)
    m1, m2 = root / "m1.json", root / "m2.json"
    sm.write(root, m1)
    sm.write(root, m2)
    assert m1.read_bytes() == m2.read_bytes()


@pytest.mark.parametrize("with_git", [True, False], ids=["checkout", "source-archive"])
def test_nested_manifest_excludes_itself_and_keeps_other_same_named_files(tmp_path, with_git):
    sm = _load()
    root = _repo(tmp_path) if with_git else tmp_path / "archive"
    root.mkdir(exist_ok=True)
    # A different file sharing the manifest's basename is still release input.
    other = root / "SOURCE_MANIFEST.json"
    other.write_text('{"historical": true}\n')
    if with_git:
        _git(root, "add", other.name)
    manifest = root / "release" / "SOURCE_MANIFEST.json"
    assert sm.write(root, manifest) == 0  # creates release/ in a fresh archive
    if with_git:
        _git(root, "add", "release/SOURCE_MANIFEST.json")
    assert sm.check(root, manifest) == 0
    before = manifest.read_bytes()
    assert sm.write(root, manifest) == 0
    assert manifest.read_bytes() == before
    listed = {f["path"] for f in json.loads(manifest.read_text())["files"]}
    assert "SOURCE_MANIFEST.json" in listed
    assert "release/SOURCE_MANIFEST.json" not in listed
