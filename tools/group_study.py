#!/usr/bin/env python3
"""Write one study-group matrix by assigning LFQ features directly through peptides."""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fasta_lake.study import group_searches  # noqa: E402


def main():
    """Parse grouping options, validate optional backends, then group and report the study."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search", required=True, help="Completed acquisition search directories")
    parser.add_argument("--out", required=True, help="New output directory outside the inputs")
    parser.add_argument("--peptide-q", type=float, default=0.01)
    parser.add_argument("--samples-file", type=Path, help="Exact acquisition names, one per line")
    parser.add_argument(
        "--config-name", choices=["config.json", "results.json"], default="config.json"
    )
    parser.add_argument(
        "--quantification",
        choices=["sum", "directlfq"],
        help="Also write selected quantities under out/quantification",
    )
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument(
        "--skip-qc",
        action="store_true",
        help="Explicitly omit AlphaPeptTools QC/PCA for a core-only installation",
    )
    args = parser.parse_args()
    try:
        if args.threads < 1:
            raise ValueError("threads must be positive")
        if args.quantification == "directlfq":
            from fasta_lake.quantification import require_directlfq

            require_directlfq()
        if not args.skip_qc:
            from fasta_lake.downstream import require_analysis

            require_analysis()
        names = None
        if args.samples_file:
            names = [
                line.strip() for line in args.samples_file.read_text().splitlines() if line.strip()
            ]
        result = group_searches(
            args.search,
            args.out,
            args.peptide_q,
            config_name=args.config_name,
            acquisition_names=names,
        )
        if args.quantification:
            from fasta_lake.quantification import quantify_study

            result["quantification"] = quantify_study(
                args.out,
                Path(args.out) / "quantification",
                method=args.quantification,
                threads=args.threads,
            )
        if not args.skip_qc:
            from fasta_lake.downstream import analyze_study

            result["qc"] = analyze_study(
                args.out,
                Path(args.out) / "qc",
                quantification=Path(args.out) / "quantification" if args.quantification else None,
            )
        else:
            result["qc"] = {"status": "SKIPPED: requested with --skip-qc"}
        (Path(args.out) / "WORKFLOW.json").write_text(json.dumps(result, indent=2) + "\n")
    except (OSError, ValueError, RuntimeError, KeyError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
