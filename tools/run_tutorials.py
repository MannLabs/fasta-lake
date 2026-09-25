#!/usr/bin/env python3
"""Execute the tutorial notebooks against real example outputs.

Each notebook under docs/tutorials/ reads the outputs of the bundled examples
from the directory named by FASTALAKE_TUTORIAL_DATA (default: tutorial_data/,
see --layout). `--write` executes every notebook in place, so the stored
outputs on the documentation site come from a real run; `--check` executes
into a temporary directory and fails on the first error, which is what CI runs
after the examples have produced their outputs.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import tempfile
from pathlib import Path

import nbformat
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = ROOT / "docs" / "tutorials"
LAYOUT = {
    "campi": "examples/campi/run_example.py --out <data>/campi",
    "multiomics": (
        "examples/multiomics/run_example.py --campi <data>/campi --out <data>/multiomics"
    ),
    "community": (
        "examples/community/run_example.py --groups <data>/campi/study --out <data>/community"
    ),
    "analysis": "examples/analysis/run_example.py --out <data>/analysis",
    "functional": "examples/functional/run_example.py --out <data>/functional",
}


def execute(path: Path, data: Path, out: Path) -> None:
    nb = nbformat.read(path, as_version=4)
    client = NotebookClient(
        nb,
        timeout=600,
        kernel_name="python3",
        resources={"metadata": {"path": str(path.parent)}},
        allow_errors=False,
    )
    os.environ["FASTALAKE_TUTORIAL_DATA"] = str(data)
    os.environ.setdefault("MPLBACKEND", "Agg")
    client.execute()
    nbformat.write(nb, out)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--write", action="store_true", help="execute in place, updating stored outputs"
    )
    mode.add_argument("--check", action="store_true", help="execute into a temporary directory")
    mode.add_argument(
        "--layout", action="store_true", help="print the commands that produce the data"
    )
    parser.add_argument(
        "--data", default=os.environ.get("FASTALAKE_TUTORIAL_DATA", "tutorial_data")
    )
    parser.add_argument("notebooks", nargs="*", help="subset of notebooks (default: all)")
    args = parser.parse_args()
    if args.layout:
        for name, command in LAYOUT.items():
            print(f"{name:12s} python {command}")
        return 0
    data = Path(args.data).resolve()
    missing = [name for name in LAYOUT if not (data / name).is_dir()]
    if missing:
        print(f"missing example outputs under {data}: {', '.join(missing)}; see --layout")
        return 2
    paths = [Path(n) for n in args.notebooks] or sorted(NOTEBOOKS.glob("*.ipynb"))
    failed = 0
    tmp = None if args.write else Path(tempfile.mkdtemp(prefix="fl_tutorials_"))
    for path in paths:
        out = path if tmp is None else tmp / path.name
        try:
            execute(path, data, out)
            print(f"ok      {path.relative_to(ROOT)}")
        except Exception as exc:  # noqa: BLE001 - report every notebook, then fail
            failed += 1
            print(f"FAILED  {path.relative_to(ROOT)}: {type(exc).__name__}: {str(exc)[:400]}")
    if tmp is not None:
        shutil.rmtree(tmp, ignore_errors=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
