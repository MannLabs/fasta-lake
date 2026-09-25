#!/usr/bin/env python3
"""Check a source release for hardcoded laboratory paths, including source archives.

Scans the git-tracked files when run inside a checkout, so run outputs written
next to the source (example results, run records with absolute paths) do not
count; in an extracted archive without git it scans every text file instead.
"""

import gzip
import io
import re
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PATTERNS = [
    re.compile(
        r"(?<![A-Za-z0-9_./])/(?:fs|home|Users|scratch|ptmp|lustre|gpfs|nfs|mnt|cluster)/[^\s\x00<>\"']+"
    ),
    re.compile(r"(?i)\b[A-Z]:\\Users\\[^\s\x00<>\"']+"),
    re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----"),
    re.compile(
        r"\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}|sk-(?:proj-|ant-)?[A-Za-z0-9_-]{30,}|xox[baprs]-[A-Za-z0-9-]{20,})"
    ),
]
# notebook stored outputs, data tables and browser assets leak paths as readily as code
TEXT = {
    ".py",
    ".rs",
    ".sh",
    ".sbatch",
    ".json",
    ".toml",
    ".md",
    ".yml",
    ".yaml",
    ".txt",
    ".cff",
    ".ipynb",
    ".tsv",
    ".csv",
    ".html",
    ".mjs",
}
EXCLUDED = {".git", "target", "__pycache__", ".venv", "build", "dist", "share"}


def candidate_files() -> list[Path]:
    """Tracked files in a git checkout; every file otherwise."""
    try:
        out = subprocess.run(
            ["git", "-C", str(REPO), "ls-files", "-z", "--", "."],
            check=True,
            capture_output=True,
        ).stdout
        files = [REPO / p for p in out.decode().split("\0") if p]
    except (OSError, subprocess.CalledProcessError):
        files = list(REPO.rglob("*"))
    return [
        p for p in files if p.is_file() and not EXCLUDED.intersection(p.relative_to(REPO).parts)
    ]


def inspect_payload(data: bytes, label: str, depth: int = 0) -> list[str]:
    """Inspect metadata and archives; never print matched values."""
    if depth > 4:
        return [f"{label}: archive nesting exceeds inspection limit"]
    content = data
    if label.lower().endswith(".mzml"):
        content = re.sub(rb"(<binary>)[^<]*(</binary>)", rb"\1\2", content)
    if label.lower().endswith((".fasta", ".faa", ".fa")):
        content = b"\n".join(
            line if line.startswith(b">") else b"" for line in content.splitlines()
        )
    text = content.decode("utf-8", errors="replace")
    hits = []
    for pattern in PATTERNS:
        match = pattern.search(text)
        if match:
            line = text.count("\n", 0, match.start()) + 1
            hits.append(f"{label}:{line}: private path or credential pattern")
    limit = 512 * 1024 * 1024
    try:
        if data.startswith(b"PK\x03\x04"):
            with zipfile.ZipFile(io.BytesIO(data)) as archive:
                if sum(x.file_size for x in archive.infolist()) > limit:
                    return hits + [f"{label}: archive exceeds inspection limit"]
                for member in archive.infolist():
                    if not member.is_dir():
                        hits += inspect_payload(
                            archive.read(member), label + "!" + member.filename, depth + 1
                        )
        elif data.startswith(b"\x1f\x8b"):
            with gzip.GzipFile(fileobj=io.BytesIO(data)) as stream:
                expanded = stream.read(limit + 1)
            if len(expanded) > limit:
                return hits + [f"{label}: gzip exceeds inspection limit"]
            if label.endswith((".tar.gz", ".tgz")):
                with tarfile.open(fileobj=io.BytesIO(expanded)) as archive:
                    for member in archive:
                        if member.isfile():
                            hits += inspect_payload(
                                archive.extractfile(member).read(),
                                label + "!" + member.name,
                                depth + 1,
                            )
            else:
                hits += inspect_payload(expanded, label.removesuffix(".gz"), depth + 1)
    except (OSError, EOFError, ValueError, zipfile.BadZipFile, tarfile.TarError):
        hits.append(f"{label}: cannot inspect archive")
    return hits


def main():
    if "--selftest" in sys.argv:
        dirty = "/" + "fs/pool/fixture"
        assert any(p.search(dirty) for p in PATTERNS)
        assert not any(p.search("data/reference.fasta") for p in PATTERNS)
        print("selftest: PASS")
        return 0
    paths = candidate_files()
    hits = [
        hit
        for p in paths
        for hit in inspect_payload(p.read_bytes(), p.relative_to(REPO).as_posix())
    ]
    print("\n".join(hits) if hits else f"OK: {len(paths)} tracked files, including archive members")
    return int(bool(hits))


if __name__ == "__main__":
    raise SystemExit(main())
