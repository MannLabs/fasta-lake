"""Compare chunked and unchunked construction on the measured teaching inputs."""

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campi", required=True, type=Path)
    parser.add_argument("--hstool", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument(
        "--laptop", action="store_true", help="Also exercise automatic resource planning"
    )
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    checked = []
    for cohort in ("campi", "hstool"):
        source = getattr(args, cohort)
        if source is None:
            continue
        source = source.resolve()
        destination = out / cohort
        argv = [
            sys.executable,
            str(ROOT / "tools/run_manifest.py"),
            "--manifest",
            str(ROOT / "examples" / cohort / "inputs/samples.tsv"),
            "--lake",
            str(source / "lake.fasta"),
            "--out",
            str(destination),
            "--reference-chunk-mib",
            "1" if cohort == "campi" else "16",
        ]
        argv += ["--laptop"] if args.laptop else ["--threads", str(args.threads)]
        with (out / (cohort + ".log")).open("x") as log:
            completed = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT)
        if completed.returncode:
            raise RuntimeError(f"Chunked {cohort} workflow failed; inspect retained log")
        records = json.loads((destination / "evidence_pieces/COMPLETE.json").read_text())
        assert records["stats"]["pieces"] > 1
        exact = ["evidence.fasta"]
        with (ROOT / "examples" / cohort / "inputs/samples.tsv").open() as stream:
            samples = list(csv.DictReader(stream, delimiter="\t"))
        for row in samples:
            sample = row["sample"]
            exact.extend(
                [
                    f"samples/{sample}/sample_lake.fasta",
                    f"samples/{sample}/inference/{sample}_razor.fasta",
                ]
            )
        for name in exact:
            assert (destination / name).read_bytes() == (source / "workflow" / name).read_bytes()
        before = pd.read_csv(source / "study/study_group_matrix.tsv", sep="\t", index_col=0)
        after = pd.read_csv(
            destination / "study_groups/study_group_matrix.tsv", sep="\t", index_col=0
        )
        assert before.index.equals(after.index) and before.columns.equals(after.columns)
        np.testing.assert_allclose(before.to_numpy(), after.to_numpy(), rtol=1e-9, equal_nan=True)
        qc = json.loads((destination / "study_groups/qc/summary.json").read_text())
        assert qc["status"] == "PASS" and qc["method"] == "sum"
        assert (destination / "study_groups/qc/QC_overview.png").stat().st_size > 1000
        plan = None
        if args.laptop:
            plan = json.loads((destination / "RESOURCE_PLAN.json").read_text())
            assert plan["threads_per_process"] <= plan["capacity"]["cpu_slots"]
            assert plan["memory_budget_bytes"] <= plan["capacity"]["memory_available_bytes"]
            assert plan["parallel_processes"] == 1
        checked.append(
            {
                "cohort": cohort,
                "pieces": records["stats"]["pieces"],
                "exact_fasta_comparisons": exact,
                "study_matrix_axes_and_values": "PASS",
                "automatic_qc": "PASS",
                "resource_plan": plan,
            }
        )
    receipt = {
        "status": "PASS",
        "cohorts": checked,
        "scope": "Measured cropped acquisitions; execution parity, not full-lake capacity",
    }
    (out / "VALIDATION.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
