"""
Reproduce the published FL_MO headline numbers from committed artifacts.

The end-to-end FL_MO arm (per-sample SAGE searches over 306 samples) is too
heavy to rerun here, so its final per-sample depth was frozen into committed
per-sample TSVs by `manuscript/verification/scripts/fdr_gate_benchmark.py`:

    manuscript/verification/fdr_gate_benchmark_FL_MO_razor_persample.tsv
    manuscript/verification/fdr_gate_benchmark_FL_razor_persample.tsv

The published razor depth numbers (Methods §M5.5 / Results, defect-ledger D3)
are the per-sample MEDIAN of the corrected-gate (protein_q <= 0.01) column:

    FL_MO razor : 13,776 protein groups   (this module reproduces this)
    FL    razor : 12,275 protein groups   (de-novo-only counterpart)

This module recomputes that median from the committed TSV so the headline
number is regenerable from committed code, and asserts it matches 13,776.

Run:
    python -m fasta_lake.multi_omics reproduce
or:
    python -m fasta_lake.multi_omics.reproduce
"""

from __future__ import annotations

import csv
import statistics
from pathlib import Path

# Committed verification artifacts (frozen final per-sample depth).
#
# These SHIP WITH THE PACKAGE. They used to be read from `manuscript/verification/`
# outside the package, where they were never committed to any repository -- so
# `reproduce()` and its test silently skipped on every clone, including this one.
# A verification that cannot run where it is published verifies nothing.
_VERIF = Path(__file__).resolve().parent / "verification"
FL_MO_RAZOR_TSV = _VERIF / "fdr_gate_benchmark_FL_MO_razor_persample.tsv"
FL_RAZOR_TSV = _VERIF / "fdr_gate_benchmark_FL_razor_persample.tsv"

# The corrected gate the manuscript reports (protein-level search FDR <= 1%).
DEPTH_COLUMN = "grp_protein_q_0.01"

PUBLISHED_FL_MO_RAZOR_MEDIAN = 13776


def per_sample_median(tsv_path, column: str = DEPTH_COLUMN) -> tuple[float, int]:
    """Return (median, n_samples) of an integer depth column in a per-sample TSV."""
    with open(tsv_path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    vals = [int(r[column]) for r in rows]
    return statistics.median(vals), len(vals)


def reproduce(verbose: bool = True) -> dict:
    """Recompute and check the published FL_MO razor median."""
    fl_mo_median, n_mo = per_sample_median(FL_MO_RAZOR_TSV)
    fl_median, n_fl = per_sample_median(FL_RAZOR_TSV)
    lift_pct = 100.0 * (fl_mo_median - fl_median) / fl_median if fl_median else 0.0
    match = round(fl_mo_median) == PUBLISHED_FL_MO_RAZOR_MEDIAN

    result = {
        "fl_mo_razor_median": fl_mo_median,
        "fl_razor_median": fl_median,
        "n_samples_fl_mo": n_mo,
        "n_samples_fl": n_fl,
        "lift_pct": lift_pct,
        "published_fl_mo_razor_median": PUBLISHED_FL_MO_RAZOR_MEDIAN,
        "matches_published": match,
    }
    if verbose:
        print("FL_MO reproduction (gate: protein_q <= 0.01, razor):")
        print(f"  FL    razor per-sample median = {fl_median:.0f}  (n={n_fl})")
        print(f"  FL_MO razor per-sample median = {fl_mo_median:.0f}  (n={n_mo})")
        print(f"  per-sample lift               = +{lift_pct:.1f}%")
        print(f"  published FL_MO razor median  = {PUBLISHED_FL_MO_RAZOR_MEDIAN}")
        print(f"  REPRODUCED: {'MATCH' if match else 'MISMATCH'}")
    return result


if __name__ == "__main__":
    r = reproduce()
    raise SystemExit(0 if r["matches_published"] else 1)
