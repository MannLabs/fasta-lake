#!/usr/bin/env python3
"""Render the generated parts of the documentation reference from the code.

Writes docs/reference/generated/run_manifest_help.txt (the study runner's --help)
and one example configuration per `fasta-lake init` template. `--check` fails when
a file differs from what the current code produces, so CI keeps the site honest.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "reference" / "generated"
TEMPLATES = ("minimal", "clinical", "metaproteomics", "full")


def render() -> dict[str, str]:
    env = dict(os.environ, COLUMNS="88", PYTHONUTF8="1")
    files: dict[str, str] = {}
    help_text = subprocess.run(
        [sys.executable, str(ROOT / "tools" / "run_manifest.py"), "--help"],
        check=True,
        capture_output=True,
        text=True,
        env=env,
        cwd=ROOT,
    ).stdout
    files["run_manifest_help.txt"] = help_text
    for template in TEMPLATES:
        config = subprocess.run(
            [sys.executable, "-m", "fasta_lake.cli", "init", "--template", template],
            check=True,
            capture_output=True,
            text=True,
            env=env,
            cwd=ROOT,
        ).stdout
        files[f"fastalake_{template}.json"] = config
    return files


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true", help="write the generated files")
    mode.add_argument("--check", action="store_true", help="fail if any generated file is stale")
    args = parser.parse_args()
    files = render()
    OUT.mkdir(parents=True, exist_ok=True)
    stale = []
    for name, text in files.items():
        path = OUT / name
        if args.write:
            path.write_text(text)
        elif not path.exists() or path.read_text() != text:
            stale.append(name)
    if args.write:
        print(f"wrote {len(files)} files to {OUT.relative_to(ROOT)}")
        return 0
    if stale:
        print("stale generated reference files: " + ", ".join(stale))
        print("regenerate with: python tools/render_reference.py --write")
        return 1
    print(f"OK: {len(files)} generated reference files are current")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
