#!/usr/bin/env python3
"""Compatibility wrapper; the installed command is ``fasta-lake group-study``."""

from fasta_lake.study_cli import *  # noqa: F401,F403
from fasta_lake.study_cli import main

if __name__ == "__main__":
    raise SystemExit(main())
