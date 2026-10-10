#!/usr/bin/env python3
"""Compatibility wrapper; the installed command is ``fasta-lake run``."""

from fasta_lake.workflow import *  # noqa: F401,F403
from fasta_lake.workflow import entry

if __name__ == "__main__":
    raise SystemExit(entry())
