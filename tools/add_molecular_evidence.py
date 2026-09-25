#!/usr/bin/env python3
"""Append verified, specimen-matched TPM-supported exact proteins after razor."""

import argparse
import json

from fasta_lake.multi_omics.exact import build_union


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--baseline", required=True)
    p.add_argument("--source", required=True)
    p.add_argument("--specimen", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--dna-percent", type=float, default=1.0)
    p.add_argument("--dna-absolute", type=float)
    p.add_argument("--rna-absolute", type=float)
    p.add_argument("--rank-rounding", choices=["floor", "ceil"], default="floor")
    a = p.parse_args()
    print(
        json.dumps(
            build_union(
                a.baseline,
                a.source,
                a.specimen,
                a.out,
                dna_percent=a.dna_percent,
                dna_absolute=a.dna_absolute,
                rna_absolute=a.rna_absolute,
                rounding=a.rank_rounding,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
