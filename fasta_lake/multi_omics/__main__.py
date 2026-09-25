"""
CLI for the FL_MO multi-omics evidence gate.

    python -m fasta_lake.multi_omics gate    --sample S --denovo-fasta ... [--gate or|and]
    python -m fasta_lake.multi_omics reproduce

This is also wired into the main `fasta-lake` CLI as `fasta-lake fl-mo`.
"""

from __future__ import annotations

import argparse
import sys

from fasta_lake.multi_omics.evidence_gate import Gate
from fasta_lake.multi_omics.reproduce import reproduce
from fasta_lake.multi_omics.runner import run_fl_mo_sample


def _add_gate_args(p: argparse.ArgumentParser) -> None:
    """Add options for the historical clustered FL_MO evidence-gate command."""
    p.add_argument("--sample", required=True)
    p.add_argument(
        "--denovo-fasta", required=True, help="Stage-3 de-novo extracted FASTA (evidence (a))"
    )
    p.add_argument(
        "--tpm-tsv", required=True, help="Per-sample TPM TSV with metaG_tpm/metaT_tpm columns"
    )
    p.add_argument(
        "--assembly-fasta", required=True, help="Per-sample MEGAHIT predicted-genes FASTA"
    )
    p.add_argument(
        "--hash-to-rep", required=True, help="FL_MO hash_to_rep.tsv (sha256 -> cluster rep id)"
    )
    p.add_argument(
        "--cluster-rep-fasta", required=True, help="FL_MO cluster97 representative FASTA"
    )
    p.add_argument("--output-fasta", required=True)
    p.add_argument("--output-lookup", required=True)
    p.add_argument("--output-stats", default=None)
    p.add_argument(
        "--gate",
        choices=[Gate.OR, Gate.AND],
        default=Gate.OR,
        help="or = FL_MO default; and = low-FDR option",
    )
    p.add_argument("--metag-top-pct", type=float, default=1.0)
    p.add_argument("--metat-abs", type=float, default=1.0)


def main(argv=None) -> int:
    """Dispatch the historical FL_MO gate or frozen-median reproduction command."""
    ap = argparse.ArgumentParser(
        prog="fasta_lake.multi_omics", description="FL_MO multi-omics evidence gate (M5.5)"
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("gate", help="Run the gate for one sample at Stage 3")
    _add_gate_args(g)

    sub.add_parser("reproduce", help="Reproduce published FL_MO razor median (13,776)")

    args = ap.parse_args(argv)

    if args.cmd == "reproduce":
        r = reproduce()
        return 0 if r["matches_published"] else 1

    if args.cmd == "gate":
        stats = run_fl_mo_sample(
            sample=args.sample,
            denovo_fasta=args.denovo_fasta,
            tpm_tsv=args.tpm_tsv,
            assembly_fasta=args.assembly_fasta,
            hash_to_rep=args.hash_to_rep,
            cluster_rep_fasta=args.cluster_rep_fasta,
            output_fasta=args.output_fasta,
            output_lookup=args.output_lookup,
            output_stats=args.output_stats,
            gate=args.gate,
            metag_top_pct=args.metag_top_pct,
            metat_abs=args.metat_abs,
        )
        print(
            f"[FL_MO {args.gate}-gate] sample={stats['sample']} "
            f"denovo={stats['fasta_denovo_written']:,} "
            f"rescued={stats['fasta_rescued_written']:,} "
            f"total={stats['fasta_total']:,} "
            f"(OR={stats['or_gate_members']:,} AND={stats['and_gate_members']:,})"
        )
        return 0

    return 2


if __name__ == "__main__":
    sys.exit(main())
