#!/usr/bin/env python3
"""Run the two bundled HSTOOL technical-acquisition windows and verify results."""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "campi"))
from run_example import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main(example_dir=HERE, cohort="HSTOOL"))
