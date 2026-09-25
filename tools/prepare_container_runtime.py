#!/usr/bin/env python3
"""Prepare the existing example runtime contract from installed container binaries."""

from __future__ import annotations

import json
import platform
import shlex
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "examples/campi"))
from run_example import (  # noqa: E402
    BINARIES,
    SAGE_ARCHIVE,
    SAGE_SHA256,
    sha,
    source_identity,
)


def prepare(runtime: Path) -> None:
    runtime = runtime.resolve()
    archive_path = ROOT / "examples/campi/runtime" / SAGE_ARCHIVE
    if sha(archive_path) != SAGE_SHA256:
        raise ValueError("Bundled Sage archive checksum mismatch")
    runtime.mkdir(parents=True, exist_ok=False)
    (runtime / "bin").mkdir()
    (runtime / "venv/bin").mkdir(parents=True)
    (runtime / "venv/bin/python").write_text(
        "#!/bin/sh\nexec " + shlex.quote(sys.executable) + ' "$@"\n'
    )
    (runtime / "venv/bin/python").chmod(0o755)
    for name in BINARIES:
        binary = shutil.which(name)
        if binary is None:
            raise ValueError("Missing prebuilt binary: " + name)
        (runtime / "bin" / name).symlink_to(Path(binary).resolve())
    with tarfile.open(archive_path) as archive:
        member = archive.getmember("sage-v0.14.6-x86_64-unknown-linux-gnu/sage")
        if not member.isfile():
            raise ValueError("Sage archive entry is not a regular file")
        with archive.extractfile(member) as source, (runtime / "bin/sage").open("xb") as target:
            shutil.copyfileobj(source, target)
    (runtime / "bin/sage").chmod(0o755)
    version = subprocess.check_output([runtime / "bin/sage", "--version"], text=True).strip()
    if version != "sage 0.14.6":
        raise ValueError("Unexpected Sage version: " + version)
    (runtime / "RUNTIME.json").write_text(
        json.dumps(
            {
                "source_sha256": source_identity(),
                "platform": platform.platform(),
                "binaries_sha256": {
                    name: sha(runtime / "bin" / name) for name in (*BINARIES, "sage")
                },
            },
            indent=2,
        ) + "\n"
    )


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: prepare_container_runtime.py NEW_RUNTIME_DIRECTORY")
    prepare(Path(sys.argv[1]))
