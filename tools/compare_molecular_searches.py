#!/usr/bin/env python3
"""Exact-sequence, peptide and quantitative overlap for two matched Sage searches."""

import argparse
import json

from fasta_lake.multi_omics.complementarity import compare


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--left", required=True)
    p.add_argument("--right", required=True)
    p.add_argument("--left-label", default="metaG")
    p.add_argument("--right-label", default="de_novo")
    p.add_argument("--out", required=True)
    p.add_argument("--peptide-q", type=float, default=0.01)
    p.add_argument("--protein-q", type=float)
    p.add_argument("--compress-tables", action="store_true")
    a = p.parse_args()
    print(
        json.dumps(
            compare(
                a.left,
                a.right,
                a.out,
                q=a.peptide_q,
                protein_q=a.protein_q,
                labels=(a.left_label, a.right_label),
                compress_tables=a.compress_tables,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
