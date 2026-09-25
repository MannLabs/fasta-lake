#!/usr/bin/env python3
"""Run the experimental AlphaNovo bag/mass comparison from a frozen study manifest."""

import argparse
import csv
import gzip
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fasta_lake.bags import ARMS, compile_queries, dump, metadata, select_databases, sha
from fasta_lake.quantification import quantify_study
from fasta_lake.study import canonical_peptide, group_searches


def verify(record):
    path = Path(record["path"])
    stat = path.stat()
    for key, actual in (("bytes", stat.st_size), ("mtime_ns", stat.st_mtime_ns)):
        if key in record and record[key] != actual:
            raise ValueError("Changed input " + str(path))
    if "sha256" in record and sha(path) != record["sha256"]:
        raise ValueError("Input hash mismatch " + str(path))


def fasta(path):
    header, seq = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:], []
            else:
                seq.append(line)
    if header is not None:
        yield header, "".join(seq)


def case(manifest, index):
    task = manifest["tasks"][index]
    return Path(manifest["output"]) / task["cohort"] / task["sample"]


def candidate(manifest, index):
    started = time.monotonic()
    task = manifest["tasks"][index]
    for field in ("predictions", "baseline", "candidates"):
        verify(task[field])
    output = case(manifest, index)
    output.mkdir(parents=True, exist_ok=False)
    dump(output / "task.json", task)
    compiled = compile_queries(task["predictions"]["path"], manifest["alphanovo_source"], output)
    scan_start = time.monotonic()
    command = [
        manifest["bag_engine"]["path"],
        str(output / "queries.jsonl"),
        task["candidates"]["path"],
        str(output / "scan.json"),
    ]
    dump(output / "scan_command.json", command)
    with (
        open(output / "scan.stderr", "w") as log,
        gzip.open(output / "candidates.jsonl.gz", "wt", compresslevel=1) as result,
    ):
        process = subprocess.Popen(
            ["/usr/bin/time", "-v", *command], stdout=subprocess.PIPE, stderr=log, text=True
        )
        try:
            for line in process.stdout:
                result.write(line)
        except BaseException:
            process.terminate()
            process.wait()
            raise
        if process.wait() != 0:
            raise RuntimeError("Full-reservoir scan failed")
    scan_seconds = time.monotonic() - scan_start
    scan = json.loads((output / "scan.json").read_text())
    if scan["reservoir_records"] != task["expected_reservoir_records"]:
        raise ValueError("Incomplete reservoir scan")
    verify(task["candidates"])
    selection = select_databases(output / "candidates.jsonl.gz", output, compiled["query_units"])
    expected = {header.split()[0]: sequence for header, sequence in fasta(task["baseline"]["path"])}
    actual = {header.split()[0]: sequence for header, sequence in fasta(output / "exact.fasta")}
    baseline_check = {
        "identical_accessions_and_sequences": expected == actual,
        "expected": len(expected),
        "actual": len(actual),
        "missing_accessions": sorted(set(expected) - set(actual)),
        "added_accessions": sorted(set(actual) - set(expected)),
        "sequence_disagreements": [
            key for key in expected.keys() & actual.keys() if expected[key] != actual[key]
        ],
    }
    dump(output / "baseline_check.json", baseline_check)
    if expected != actual:
        raise ValueError(
            "Exact pipeline differs from the frozen baseline; diagnose before cohort execution"
        )
    receipt = {
        "status": "CANDIDATES_COMPLETE",
        "job_id": os.environ.get("SLURM_JOB_ID"),
        "compiled": compiled,
        "scan": scan,
        "scan_seconds": scan_seconds,
        "selection": selection,
        "elapsed_seconds": time.monotonic() - started,
        "baseline_reproduced": True,
    }
    dump(output / "CANDIDATES_COMPLETE.json", receipt)


def accepted(path):
    peptides, spectra = set(), set()
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            q = float(row["peptide_q"])
            if int(row["label"]) == 1 and math.isfinite(q) and q <= 0.01:
                if "rank" in row and int(row["rank"]) != 1:
                    raise ValueError("Unexpected accepted non-rank-one PSM")
                peptides.add(canonical_peptide(row["peptide"]))
                spectra.add((row["filename"], row["scannr"]))
    return peptides, spectra


def search(manifest, index, arms=None, baseline_arm="exact"):
    arms = ARMS if arms is None else arms
    task = manifest["tasks"][index]
    output = case(manifest, index)
    candidates = json.loads((output / "CANDIDATES_COMPLETE.json").read_text())
    assert candidates["status"] == "CANDIDATES_COMPLETE"
    for record in (task["search_config"], *task["mzml"], manifest["sage"]):
        verify(record)
    version = subprocess.check_output([manifest["sage"]["path"], "--version"], text=True).strip()
    if version.split()[-1] != "0.14.6":
        raise ValueError("Sage version mismatch")
    original = json.loads(Path(task["search_config"]["path"]).read_text())
    root = output / "searches"
    root.mkdir(exist_ok=False)
    results, executions, seen = {}, {}, {}
    env = dict(
        os.environ,
        RAYON_NUM_THREADS="4",
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
    )
    for arm in arms:
        fa = candidates["selection"]["arms"][arm]["fasta"]
        verify(fa)
        if fa["bytes"] == 0:
            results[arm] = (set(), set())
            executions[arm] = {"status": "EMPTY_SELECTED_DATABASE", "seconds": 0}
            continue
        if fa["sha256"] in seen:
            other = seen[fa["sha256"]]
            (root / arm).symlink_to(other, target_is_directory=True)
            results[arm] = results[other]
            executions[arm] = {"reused_identical_fasta": other, "seconds": 0}
            continue
        directory = root / arm
        directory.mkdir()
        config = json.loads(json.dumps(original))
        config.pop("version", None)
        config.pop("output_paths", None)
        config["database"]["fasta"] = fa["path"]
        config["mzml_paths"] = [record["path"] for record in task["mzml"]]
        config["output_directory"] = str(directory)
        dump(directory / "config.json", config)
        command = [
            manifest["sage"]["path"],
            str(directory / "config.json"),
            "--write-pin",
            "--disable-telemetry-i-dont-want-to-improve-sage",
        ]
        dump(directory / "command.json", command)
        started = time.monotonic()
        with (
            open(directory / "stdout.log", "w") as stdout,
            open(directory / "stderr.log", "w") as stderr,
        ):
            subprocess.run(
                ["/usr/bin/time", "-v", *command], env=env, stdout=stdout, stderr=stderr, check=True
            )
        results[arm] = accepted(directory / "results.sage.tsv")
        executions[arm] = {
            "seconds": time.monotonic() - started,
            "results": metadata(directory / "results.sage.tsv"),
        }
        seen[fa["sha256"]] = arm
    baseline = results[baseline_arm][0]
    rows = []
    with gzip.open(output / "peptide_changes.tsv.gz", "wt") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["arm", "change", "peptide_IL"])
        for arm, (peptides, spectra) in results.items():
            gained, lost = peptides - baseline, baseline - peptides
            rows.append(
                {
                    "cohort": task["cohort"],
                    "sample": task["sample"],
                    "arm": arm,
                    "peptides": len(peptides),
                    "spectra": len(spectra),
                    "gained": len(gained),
                    "lost": len(lost),
                    "retained": len(peptides & baseline),
                    "net": len(gained) - len(lost),
                }
            )
            for label, values in (("gained", gained), ("lost", lost)):
                writer.writerows((arm, label, value) for value in sorted(values))
    dump(
        output / "SEARCHES_COMPLETE.json",
        {
            "status": "MATCHED_SEARCHES_COMPLETE",
            "rows": rows,
            "comparison_reference_arm": baseline_arm,
            "executions": executions,
            "novel_subset_error_rate": "NOT_ESTABLISHED",
            "workflow_fdr": "UNVERIFIED",
            "acceptance_rule": "Sage target peptide_q <= 0.01; nominal search threshold",
            "same_fdr_sensitivity_claim_supported": False,
            "job_id": os.environ.get("SLURM_JOB_ID"),
        },
    )


def gate(manifest):
    records = []
    for index in manifest["paired_pilot_indices"]:
        output = case(manifest, index)
        candidate_record = json.loads((output / "CANDIDATES_COMPLETE.json").read_text())
        search_record = json.loads((output / "SEARCHES_COMPLETE.json").read_text())
        assert candidate_record["baseline_reproduced"]
        assert search_record["status"] == "MATCHED_SEARCHES_COMPLETE"
        assert candidate_record["elapsed_seconds"] <= manifest["pilot_candidate_limit_seconds"]
        base = search_record["rows"][0]["peptides"]
        for row in search_record["rows"]:
            assert row["retained"] + row["lost"] == base
            assert row["retained"] + row["gained"] == row["peptides"]
        grouped = {}
        task = manifest["tasks"][index]
        for arm in list(ARMS)[:4]:
            stage = output / "pilot_group_inputs" / arm
            stage.mkdir(parents=True, exist_ok=False)
            (stage / task["sample"]).symlink_to(output / "searches" / arm, target_is_directory=True)
            target = output / "pilot_groups" / arm
            summary = group_searches(stage, target, 0.01, acquisition_names=[task["sample"]])
            summary["quantification"] = quantify_study(
                target, target / "quantification", method="sum", threads=1
            )
            dump(target / "WORKFLOW.json", summary)
            grouped[arm] = metadata(target / "WORKFLOW.json")
        records.append(
            {
                "index": index,
                "candidate_seconds": candidate_record["elapsed_seconds"],
                "searches": search_record["executions"],
                "pilot_grouping_and_sum": grouped,
            }
        )
    dump(
        Path(manifest["output"]).parent / "PILOT_GATE.json",
        {
            "status": "PASS",
            "gate_scope": "EXECUTION_AND_ACCOUNTING_ONLY",
            "workflow_fdr": "UNVERIFIED",
            "same_fdr_sensitivity_claim_supported": False,
            "pilots": records,
            "positive_gain_required": False,
        },
    )


def group(manifest, index):
    cohort = list(manifest["expected_cohort_counts"])[index // 4]
    arm = list(ARMS)[index % 4]
    root = Path(manifest["output"]).parent
    stage = root / "group_inputs" / cohort / arm
    stage.mkdir(parents=True, exist_ok=False)
    names = []
    for i, task in enumerate(manifest["tasks"]):
        if task["cohort"] != cohort:
            continue
        receipt = json.loads((case(manifest, i) / "SEARCHES_COMPLETE.json").read_text())
        assert receipt["status"] == "MATCHED_SEARCHES_COMPLETE"
        (stage / task["sample"]).symlink_to(
            case(manifest, i) / "searches" / arm, target_is_directory=True
        )
        names.append(task["sample"])
    assert len(names) == manifest["expected_cohort_counts"][cohort]
    output = root / "groups" / cohort / arm
    summary = group_searches(stage, output, 0.01, acquisition_names=names)
    summary["quantification"] = quantify_study(
        output, output / "quantification", method="sum", threads=4
    )
    summary["qc"] = (
        "Core coverage/missingness recorded; optional AlphaPeptTools PCA is a separate endpoint"
    )
    dump(output / "WORKFLOW.json", summary)


def collect(manifest):
    rows, errors, timing = [], [], []
    for index, task in enumerate(manifest["tasks"]):
        try:
            output = case(manifest, index)
            cand = json.loads((output / "CANDIDATES_COMPLETE.json").read_text())
            result = json.loads((output / "SEARCHES_COMPLETE.json").read_text())
            rows.extend(result["rows"])
            timing.append(
                {
                    "cohort": task["cohort"],
                    "sample": task["sample"],
                    "compile_seconds": cand["compiled"]["elapsed_seconds"],
                    "scan_seconds": cand["scan_seconds"],
                    "selection_seconds": cand["selection"]["elapsed_seconds"],
                    "total_candidate_seconds": cand["elapsed_seconds"],
                    "searches": result["executions"],
                }
            )
        except (OSError, KeyError, ValueError) as error:
            errors.append({"index": index, "sample": task["sample"], "error": str(error)})
    root = Path(manifest["output"]).parent
    if rows:
        with (root / "RESULTS.tsv").open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    dump(root / "TIMINGS.json", timing)
    dump(
        root / "COLLECTION.json",
        {
            "status": "INCOMPLETE" if errors else "SEARCH_ACCOUNTING_COMPLETE",
            "errors": errors,
            "completed": len(timing),
            "expected": len(manifest["tasks"]),
            "participant_statistics": "PENDING",
            "scientific_review": "PENDING",
            "workflow_fdr": "UNVERIFIED",
            "same_fdr_sensitivity_claim_supported": False,
        },
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("stage", choices=["candidate", "search", "gate", "group", "collect"])
    parser.add_argument("index", type=int, nargs="?", default=0)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    for record in manifest["code"] + [manifest["bag_engine"]]:
        verify(record)
    if args.stage in {"candidate", "search"} and args.index not in manifest["paired_pilot_indices"]:
        gate_record = json.loads((args.manifest.parent / "PILOT_GATE.json").read_text())
        assert gate_record["status"] == "PASS"
    functions = {"candidate": candidate, "search": search, "group": group}
    if args.stage in functions:
        functions[args.stage](manifest, args.index)
    elif args.stage == "gate":
        gate(manifest)
    else:
        collect(manifest)


if __name__ == "__main__":
    main()
