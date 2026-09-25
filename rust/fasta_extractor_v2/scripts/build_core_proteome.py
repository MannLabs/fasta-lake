#!/usr/bin/env python3
"""
Build core proteome FASTA from evidence manifest.

Filters the evidence lake to proteins detected in >= min_samples samples
(or >= min_fraction of total samples). These "core" proteins are injected
into every per-sample database via fasta_extractor --include, ensuring
cross-sample comparability without fabricating identifications.

Usage:
  python build_core_proteome.py \
    --evidence evidence_manifest.tsv \
    --fasta evidence_lake.fasta \
    --output core_proteome.fasta \
    --min-fraction 0.10  # proteins in >=10% of samples

The search engine (SAGE/DIA-NN) still controls FDR — this only expands
the search space with high-confidence proteins.
"""

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Build core proteome from evidence manifest")
    parser.add_argument("--evidence", required=True, help="Evidence manifest TSV")
    parser.add_argument("--fasta", required=True, help="Evidence lake FASTA")
    parser.add_argument("--output", required=True, help="Output core proteome FASTA")
    parser.add_argument(
        "--min-samples", type=int, default=0, help="Min samples a protein must appear in (absolute)"
    )
    parser.add_argument(
        "--min-fraction",
        type=float,
        default=0.10,
        help="Min fraction of samples (default: 0.10 = 10%%)",
    )
    parser.add_argument(
        "--min-peptides", type=int, default=2, help="Min unique peptides (default: 2)"
    )
    parser.add_argument(
        "--min-support", type=float, default=0.0, help="Min support_score (default: 0.0)"
    )
    args = parser.parse_args()

    evidence_path = Path(args.evidence)
    fasta_path = Path(args.fasta)
    output_path = Path(args.output)

    # First pass: read all entries and find max n_samples (= total samples)
    entries = []
    max_samples = 0

    with open(evidence_path) as f:
        header = f.readline().strip().split("\t")
        idx = {col: i for i, col in enumerate(header)}

        for line in f:
            parts = line.strip().split("\t")
            n_peptides = int(parts[idx["n_unique_peptides"]])
            n_samples = int(parts[idx["n_samples"]])
            support = float(parts[idx["support_score"]])
            protein_id = parts[idx["protein_id"]]
            entries.append((protein_id, n_peptides, n_samples, support))
            if n_samples > max_samples:
                max_samples = n_samples

    total_proteins = len(entries)

    # Determine effective min_samples threshold
    # Priority: explicit --min-samples > --min-fraction > default (no filter)
    effective_min_samples = args.min_samples
    if effective_min_samples == 0 and args.min_fraction > 0 and max_samples > 1:
        effective_min_samples = max(1, int(max_samples * args.min_fraction))
        print(f"  Inferred total samples: {max_samples}")
        print(f"  min_fraction {args.min_fraction:.0%} -> min_samples = {effective_min_samples}")

    # Second pass: filter
    core_ids = set()
    for protein_id, n_peptides, n_samples, support in entries:
        if n_peptides < args.min_peptides:
            continue
        if effective_min_samples > 0 and n_samples < effective_min_samples:
            continue
        if args.min_support > 0 and support < args.min_support:
            continue
        core_ids.add(protein_id)

    print(f"Evidence manifest: {total_proteins} proteins ({max_samples} samples detected)")
    print(
        f"Filters: min_peptides={args.min_peptides}, min_samples={effective_min_samples}, "
        f"min_support={args.min_support}"
    )
    print(f"Core proteome candidates: {len(core_ids)} proteins")

    # Extract core proteins from evidence FASTA
    written = 0
    writing = False

    with open(fasta_path) as f_in, open(output_path, "w") as f_out:
        for line in f_in:
            if line.startswith(">"):
                protein_id = line[1:].split()[0]
                writing = protein_id in core_ids
                if writing:
                    written += 1
            if writing:
                f_out.write(line)

    print(f"Written: {written} proteins to {output_path}")
    print(f"Coverage: {written}/{total_proteins} ({written / total_proteins * 100:.1f}%)")


if __name__ == "__main__":
    main()
