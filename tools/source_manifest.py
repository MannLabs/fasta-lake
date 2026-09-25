#!/usr/bin/env python3
"""Write or check release/SOURCE_MANIFEST.json for the tracked repository tree.

    python tools/source_manifest.py --write   # regenerate from the tracked tree
    python tools/source_manifest.py --check   # exit 1 if the manifest and the tree differ

The manifest lists every git-tracked file in the repository (paths relative to
the root, sorted) except itself. It is regenerated with --write on the release
branch and verified with --check in CI on release branches and tags, so a
release cannot ship a manifest that describes another tree. Between releases
the tracked manifest may lag the tree; it is not enforced on ordinary pull
requests, because automated dependency updates cannot regenerate it.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "release" / "SOURCE_MANIFEST.json"


EXCLUDED = {
    ".git",
    "target",
    "__pycache__",
    ".venv",
    "build",
    "dist",
    ".pytest_cache",
    ".ruff_cache",
}


def tracked_files(
    root: Path = ROOT, manifest_name: str = MANIFEST.relative_to(ROOT).as_posix()
) -> list[Path]:
    """Files below root, relative to root, sorted, without the manifest itself.

    Inside a git checkout these are the tracked files, so run outputs written
    next to the source do not count. In an extracted source archive (no git
    metadata, or no git on PATH) every file is taken, minus build directories.
    """
    try:
        out = subprocess.run(
            ["git", "-C", str(root), "ls-files", "-z", "--", "."],
            check=True,
            capture_output=True,
        ).stdout
        rel = sorted(p for p in out.decode().split("\0") if p)
    except (OSError, subprocess.CalledProcessError):
        rel = sorted(
            p.relative_to(root).as_posix()
            for p in root.rglob("*")
            if p.is_file() and not EXCLUDED.intersection(p.relative_to(root).parts)
        )
    return [Path(p) for p in rel if (root / p).is_file() and p != manifest_name]


def build_manifest(root: Path = ROOT, manifest: Path = MANIFEST) -> dict:
    files = []
    for rel in tracked_files(root, manifest.relative_to(root).as_posix()):
        data = (root / rel).read_bytes()
        files.append(
            {"path": rel.as_posix(), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
        )
    return {"generator": "tools/source_manifest.py", "files": files}


def write(root: Path = ROOT, manifest: Path = MANIFEST) -> int:
    built = build_manifest(root, manifest)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(built, indent=2) + "\n")
    print(f"wrote {manifest.relative_to(root)}: {len(built['files'])} files")
    return 0


def check(root: Path = ROOT, manifest: Path = MANIFEST) -> int:
    if not manifest.exists():
        print(f"missing {manifest}")
        return 1
    recorded = {f["path"]: f for f in json.loads(manifest.read_text())["files"]}
    actual = {f["path"]: f for f in build_manifest(root, manifest)["files"]}
    problems = []
    for path in sorted(set(recorded) | set(actual)):
        if path not in actual:
            problems.append(f"listed but not tracked: {path}")
        elif path not in recorded:
            problems.append(f"tracked but not listed: {path}")
        elif recorded[path]["sha256"] != actual[path]["sha256"]:
            problems.append(f"content differs: {path}")
    if problems:
        print("\n".join(problems))
        print(f"{manifest.relative_to(root)} is stale ({len(problems)} problems); run --write")
        return 1
    print(f"OK: {manifest.relative_to(root)} matches {len(actual)} tracked files")
    return 0


def main(argv: list[str]) -> int:
    if argv == ["--write"]:
        return write()
    if argv == ["--check"]:
        return check()
    print(__doc__)
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
