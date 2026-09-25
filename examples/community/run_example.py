#!/usr/bin/env python3
"""Join recorded public CAMPI taxonomy and eggNOG output without network access."""

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path


def main():
    """Check source identities, run installed commands and reconcile every sample."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    inputs = Path(__file__).resolve().parent / "inputs"
    provenance = json.loads((inputs / "PROVENANCE.json").read_text())
    for name, digest in provenance["files"].items():
        assert hashlib.sha256((inputs / name).read_bytes()).hexdigest() == digest, name
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)

    def run(*words):
        subprocess.run([sys.executable, "-m", "fasta_lake.cli", *map(str, words)], check=True)

    run("taxonomy", "prepare", "--groups", args.groups, "--out", out / "queries")
    assert (out / "queries/peptides.txt").read_bytes() == (inputs / "peptides.txt").read_bytes()
    run(
        "taxonomy",
        "report",
        "--groups",
        args.groups,
        "--unipept-output",
        inputs / "campi_unipept.json",
        "--out",
        out / "taxonomy",
        "--database",
        "UniProt 2026.02; Unipept CLI 4.2.2",
        "--equate-il",
    )
    run("network", "export", "--groups", args.groups, "--out", out / "network")
    run(
        "network",
        "eggnog",
        "--bundle",
        out / "network",
        "--annotation",
        inputs / "eggnog.annotations.tsv",
        "--query-fasta",
        inputs / "eggnog_query.fasta",
        "--fasta",
        args.groups / "representatives.fasta",
        "--out",
        out / "functions",
    )
    data = json.loads((out / "taxonomy/taxonomy.json").read_text())
    for sample in data["samples"]:
        for measure in ["peptides", "signal"]:
            total = sample["totals"][measure]
            assert abs(sum(c[measure] for c in sample["balance"].values()) - total) <= (
                1e-10 * max(1, total)
            )
            assert abs(sample["trees"][measure]["count"] - total) <= 1e-10 * max(1, total)
    functions = json.loads((out / "functions/COMPLETE.json").read_text())
    assert functions["status_counts"]["assigned"] == 1
    assert "K01689" in (out / "functions/group_terms_KEGG_ko.tsv").read_text()
    for name in ["taxonomy", "network", "functions"]:
        receipt = json.loads((out / name / "COMPLETE.json").read_text())
        assert receipt["status"] == "PASS"
        for filename, digest in receipt["files"].items():
            assert hashlib.sha256((out / name / filename).read_bytes()).hexdigest() == digest
    receipt = json.loads((out / "taxonomy/COMPLETE.json").read_text())
    assert receipt["static_plots"]["status"] == "PASS"
    for suffix in ["pdf", "png", "svg"]:
        assert (out / "taxonomy" / ("Holobiont_balance." + suffix)).stat().st_size > 1000
    (out / "VALIDATION.json").write_text(
        json.dumps(
            {
                "status": "PASS",
                "samples": len(data["samples"]),
                "network_access": False,
                "taxonomy": "Recorded native Unipept output; reference-dependent LCA",
                "function": "One real sequence-matched enolase annotation; partial coverage",
                "ckg_renderer": "Separate optional environment; not executed by this example",
                "counts_and_signal_conserved": True,
            },
            indent=2,
        )
        + "\n"
    )
    print("PASS: offline Unipept report, holobiont balance, evidence graph and exact eggNOG link")


if __name__ == "__main__":
    main()
