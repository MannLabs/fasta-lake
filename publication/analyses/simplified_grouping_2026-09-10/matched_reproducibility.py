#!/usr/bin/env python3
"""Compare correlations using exactly the same positive rows in both methods."""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd


def common_pair(original, simplified, a, b):
    shared = original[[a, b]].join(simplified[[a, b]], how="inner", rsuffix="_direct")
    shared = shared.loc[(shared > 0).all(axis=1) & np.isfinite(shared).all(axis=1)]
    values = np.log2(shared.to_numpy())
    result = {"common_positive_rows": len(shared)}
    for method, columns in [("original", (0, 1)), ("simplified", (2, 3))]:
        x, y = (values[:, i] for i in columns)
        result[method] = (
            float(np.corrcoef(x, y)[0, 1])
            if (len(x) >= 3 and np.std(x) > 0 and np.std(y) > 0)
            else None
        )
    result["change"] = (
        result["simplified"] - result["original"]
        if None not in (result["original"], result["simplified"])
        else None
    )
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--comparison", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    design_path = args.comparison / "fixed_pipeline_metrics.tsv"
    design = pd.read_csv(design_path, sep="\t")
    if len(design) != 58 or design["sample"].nunique() != 58:
        raise ValueError("Exactly 58 HSTOOL acquisitions required")
    paths = [
        design_path,
        args.comparison / "functional_profiles.tsv",
        *(args.comparison / f"{m}_study_group_matrix.tsv" for m in ["original", "simplified"]),
    ]
    functions = pd.read_csv(paths[1], sep="\t")
    matrices = {
        "groups": {
            m: pd.read_csv(args.comparison / f"{m}_study_group_matrix.tsv", sep="\t", index_col=0)
            for m in ["original", "simplified"]
        },
        "KO": {
            m: functions.query("method == @m").pivot(
                index="ko", columns="sample", values="normalised_function_intensity"
            )
            for m in ["original", "simplified"]
        },
    }
    rows = []
    for tier, population in design.groupby("tier", sort=False):
        for a, b in itertools.combinations(population["sample"], 2):
            for space, pair in matrices.items():
                rows.append(
                    dict(
                        tier=tier,
                        space=space,
                        sample_a=a,
                        sample_b=b,
                        **common_pair(pair["original"], pair["simplified"], a, b),
                    )
                )
    pairs = pd.DataFrame(rows)
    medians = (
        pairs.groupby(["tier", "space"], sort=False)
        .agg(
            pairs=("sample_a", "size"),
            common_rows=("common_positive_rows", "median"),
            original=("original", "median"),
            simplified=("simplified", "median"),
            median_paired_change=("change", "median"),
        )
        .reset_index()
    )
    args.out.mkdir(parents=True, exist_ok=False)
    pairs.to_csv(args.out / "pairs.tsv", sep="\t", index=False)
    medians.to_csv(args.out / "medians.tsv", sep="\t", index=False)
    paths.append(Path(__file__))
    provenance = [
        {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths
    ]
    (args.out / "PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(medians.to_string(index=False))


if __name__ == "__main__":
    main()
