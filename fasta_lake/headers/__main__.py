"""Allow running as: python -m fasta_lake.headers heal/audit/prepare"""

import argparse

from . import audit_fasta, heal_fasta, prepare_fasta_for_engine, print_audit


def main():
    """Parse and dispatch the FASTA header healing, audit and preparation commands."""
    parser = argparse.ArgumentParser(prog="fastalake-headers")
    sub = parser.add_subparsers(dest="cmd")

    p = sub.add_parser("heal")
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--engine", default="diann", choices=["diann", "sage", "msfragger"])
    p.add_argument("--canonical", default=None)

    p = sub.add_parser("audit")
    p.add_argument("input")

    p = sub.add_parser("prepare", help="Prepare FASTA for a specific search engine")
    p.add_argument("input", help="Input inference FASTA from FastaLake")
    p.add_argument("output", help="Engine-ready output FASTA")
    p.add_argument(
        "--engine",
        required=True,
        choices=["sage", "diann", "msfragger"],
        help="Target search engine",
    )
    p.add_argument("--canonical", default=None, help="SwissProt FASTA for gene resolution (DIA-NN)")

    args = parser.parse_args()
    if args.cmd == "heal":
        stats = heal_fasta(
            args.input, args.output, target_engine=args.engine, canonical_fasta=args.canonical
        )
        print(f"Kept: {stats.kept:,}, Healed: {stats.healed:,}, Skipped: {stats.skipped:,}")
    elif args.cmd == "audit":
        stats = audit_fasta(args.input)
        print_audit(stats)
    elif args.cmd == "prepare":
        stats = prepare_fasta_for_engine(
            args.input,
            args.output,
            engine=args.engine,
            canonical_fasta=args.canonical,
        )
        print(f"Engine: {stats.engine}")
        print(f"Target proteins: {stats.target_proteins:,}")
        if stats.decoy_proteins:
            print(f"Decoy proteins: {stats.decoy_proteins:,}")
        print(f"Total proteins: {stats.total_proteins:,}")
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
