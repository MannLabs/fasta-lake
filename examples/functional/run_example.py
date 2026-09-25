#!/usr/bin/env python3
"""Run experimental functional exploration on an explicitly synthetic paired study."""

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    names = [f"donor{d}_{condition}" for d in range(3) for condition in ("A", "B")]
    values = [
        [4, 7, 5, 8, 6, 9],
        [9, 5, 11, 7, 10, 8],
        [3, 6, 4, 5, 5, 7],
        [100, 20, 80, 40, 90, 30],
    ]
    with (out / "matrix.tsv").open("x") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["study_group", *names])
        writer.writerows(
            [[group, *row] for group, row in zip(["g1", "g2", "g3", "unknown"], values)]
        )
    (out / "terms.tsv").write_text("study_group\tterm\ng1\tK1\ng2\tK2\ng3\tK3\n")
    with (out / "metadata.tsv").open("x") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["sample", "condition", "donor"])
        writer.writerows([[name, name.rsplit("_", 1)[1], name.rsplit("_", 1)[0]] for name in names])
    subprocess.run(
        [
            sys.executable,
            "-m",
            "fasta_lake.cli",
            "explore-study",
            "--matrix",
            str(out / "matrix.tsv"),
            "--annotations",
            str(out / "terms.tsv"),
            "--metadata",
            str(out / "metadata.tsv"),
            "--group-column",
            "condition",
            "--block-column",
            "donor",
            "--permutations",
            "99",
            "--seed",
            "23",
            "--level",
            "synthetic_terms",
            "--out",
            str(out / "results"),
        ],
        check=True,
    )
    result = json.loads((out / "results/COMPLETE.json").read_text())
    assert result["experimental"] and result["outputs"]["pca"] == "PASS"
    for name, digest in result["files"].items():
        assert hashlib.sha256((out / "results" / name).read_bytes()).hexdigest() == digest
    with (out / "results/annotation_coverage.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    for i, row in enumerate(rows):
        assert abs(float(row["observed_signal"]) - sum(x[i] for x in values)) < 1e-10
        expected = sum(x[i] for x in values[:3]) / sum(x[i] for x in values)
        assert abs(float(row["annotation_fraction"]) - expected) < 1e-12
    with (out / "results/diversity.tsv").open() as stream:
        assert all(
            float(row["observed_richness"]) == 3 for row in csv.DictReader(stream, delimiter="\t")
        )
    test = json.loads((out / "results/PERMANOVA.json").read_text())
    assert test["blocked"] and test["permutations"] == 99
    (out / "VALIDATION.json").write_text(
        json.dumps(
            {
                "status": "PASS",
                "scope": "Synthetic paired functional software example",
                "biological_result": False,
                "signal_conserved": True,
                "unknown_excluded_from_diversity": True,
                "pca": True,
                "blocked_permutations": 99,
            },
            indent=2,
        )
        + "\n"
    )
    print("PASS: synthetic functional coverage, diversity, CLR PCA and blocked permutations")


if __name__ == "__main__":
    main()
