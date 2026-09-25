#!/usr/bin/env python3
"""Run a prespecified absent-species diagnostic through bag selection and Sage.

The control FASTA is included before retrieval, and shares selection and search
with ordinary targets. Observed control fractions are diagnostics, never a
calibrated estimate of target-only, protein-group or whole-workflow FDR.
"""

import argparse
import csv
import gzip
import json
import math
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fasta_lake.bags import (
    ARMS,
    compile_queries,
    dump,
    evidence_for,
    metadata,
    norm,
    select_databases,
)
from fasta_lake.study import canonical_peptide
from tools.run_bags import case, fasta, search, verify

THRESHOLDS = (0.001, 0.005, 0.01, 0.02, 0.05)


def validate_absence(evidence, sample):
    """Require an attributable, acquisition-specific absence assertion, not a guess."""
    for key in ("source", "reviewed_by", "species", "contaminant_policy"):
        if not evidence.get(key):
            raise ValueError(f"Absence evidence requires {key}")
    if evidence.get("absence_confirmed") is not True or sample not in evidence.get("samples", []):
        raise ValueError("No confirmed absence record for this acquisition")


def augment(target, entrapment, output):
    """Include complete controls before selection; keep IDs and tie-breaking unchanged."""
    controls = {}
    for header, sequence in fasta(entrapment):
        accession = header.split()[0]
        if accession in controls or not sequence:
            raise ValueError("Duplicate control accession or empty control protein")
        controls[accession] = sequence
    if not controls:
        raise ValueError("Empty entrapment FASTA")
    count = 0
    with output.open("x") as handle:
        for source in (target, entrapment):
            for header, sequence in fasta(source):
                accession = header.split()[0]
                if source == target and accession in controls:
                    raise ValueError("Target/control accession collision: " + accession)
                if not sequence or set(sequence.upper()) - set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*"):
                    raise ValueError("Invalid protein sequence: " + accession)
                # Class labels are kept outside the FASTA and never enter selection.
                handle.write(">" + header + "\n" + sequence + "\n")
                count += 1
    return controls, count


def scan(engine, queries, reservoir, output):
    """Retain every scan hit with bounded output buffering and preserve failures."""
    command = [str(engine), str(queries), str(reservoir), str(output / "scan.json")]
    dump(output / "scan_command.json", command)
    with (
        open(output / "scan.stderr", "w") as log,
        gzip.open(output / "candidates.jsonl.gz", "wt", compresslevel=1) as dest,
    ):
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log, text=True)
        try:
            for line in process.stdout:
                dest.write(line)
            if process.wait() != 0:
                raise RuntimeError("Entrapment reservoir scan failed; see scan.stderr")
        finally:
            process.stdout.close()
            if process.poll() is None:
                process.terminate()
                process.wait()


def read_psms(path):
    """Keep rank-one target PSMs and both nominal thresholds for separate endpoints."""
    rows = []
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if int(row["label"]) != 1:
                continue
            if int(row["rank"]) != 1:
                raise ValueError("Expected rank-one search output")
            pq, sq = float(row["peptide_q"]), float(row["spectrum_q"])
            if not all(math.isfinite(q) and 0 <= q <= 1 for q in (pq, sq)):
                raise ValueError("Invalid target q-value")
            rows.append((canonical_peptide(row["peptide"]), row["filename"], row["scannr"], pq, sq))
    keys = [(r[1], r[2]) for r in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Multiple rank-one targets for one spectrum")
    return rows


def membership(peptides, target, entrapment, possible_present):
    """Classify exact I/L peptide containment against complete preselection sources."""
    import ahocorasick

    if not peptides:
        return {}
    matcher = ahocorasick.Automaton()
    for peptide in peptides:
        matcher.add_word(norm(peptide), peptide)
    matcher.make_automaton()
    present, control = set(), set()
    for path, found in ((target, present), (possible_present, present), (entrapment, control)):
        for _, sequence in fasta(path):
            found.update(peptide for _, peptide in matcher.iter(norm(sequence)))
    classes = {}
    for peptide in peptides:
        if peptide in present and peptide in control:
            classes[peptide] = "shared"
        elif peptide in present:
            classes[peptide] = "target_compatible"
        elif peptide in control:
            classes[peptide] = "entrapment_only"
        else:
            raise ValueError("Accepted peptide missing from input sources: " + peptide)
    return classes


def summarize(rows, classes, baseline):
    """Report canonical peptides and PSMs separately, with exact-arm gains and losses."""
    result = []
    for q in THRESHOLDS:
        for level, column in (("canonical_peptide", 3), ("PSM", 4)):

            def identities(values):
                """Use peptide_q for peptide sets and spectrum_q for PSM assignments."""
                return {
                    (row[0] if level == "canonical_peptide" else row[:3]): row[0]
                    for row in values
                    if row[column] <= q
                }

            current, exact = identities(rows), identities(baseline)
            for population, keys in (
                ("all", current.keys()),
                ("gained_vs_exact", current.keys() - exact.keys()),
                ("lost_vs_exact", exact.keys() - current.keys()),
            ):
                source = exact if population == "lost_vs_exact" else current
                counts = Counter(classes[source[key]] for key in keys)
                n = sum(counts.values())
                e = counts["entrapment_only"]
                result.append(
                    {
                        "q": q,
                        "level": level,
                        "population": population,
                        "accepted": n,
                        "entrapment_only": e,
                        "shared": counts["shared"],
                        "target_compatible": counts["target_compatible"],
                        "observed_control_fraction": e / n if n else None,
                    }
                )
    return result


def run(manifest, index, control_spec, output):
    """Execute a fresh augmented-reservoir experiment; never alter the frozen study."""
    if index < 0 or index >= len(manifest["tasks"]):
        raise ValueError("Acquisition index is outside the frozen roster")
    if set(control_spec) != {"entrapment_fasta", "possible_present_fasta", "absence_evidence"}:
        raise ValueError("Expected exactly three control input records")
    root = Path(__file__).resolve().parents[1]
    required_code = {
        str((root / name).resolve())
        for name in ("fasta_lake/bags.py", "tools/run_bags.py", "tools/run_bag_entrapment.py")
    }
    if not required_code <= {str(Path(r["path"]).resolve()) for r in manifest["code"]}:
        raise ValueError("Manifest must freeze the executing matching module and both runners")
    task = manifest["tasks"][index]
    for record in (
        manifest["bag_engine"],
        manifest["sage"],
        *manifest["code"],
        task["candidates"],
        task["predictions"],
        task["search_config"],
        *task["mzml"],
        *control_spec.values(),
    ):
        if not record.get("sha256"):
            raise ValueError("Every input and code record must have a frozen SHA-256")
        verify(record)
    absence = json.loads(Path(control_spec["absence_evidence"]["path"]).read_text())
    validate_absence(absence, task["sample"])
    config = json.loads(Path(task["search_config"]["path"]).read_text())
    if config["database"].get("generate_decoys", True) is not True:
        raise ValueError("Entrapments must be targets alongside Sage-generated decoys")
    output.mkdir(parents=True, exist_ok=False)
    dump(output / "INPUTS.json", {"manifest": manifest, "controls": control_spec, "index": index})
    target = Path(task["candidates"]["path"])
    entrapment = Path(control_spec["entrapment_fasta"]["path"])
    possible = Path(control_spec["possible_present_fasta"]["path"])
    augmented = output / "augmented.fasta"
    controls, total = augment(target, entrapment, augmented)
    if total - len(controls) != task["expected_reservoir_records"] or total == len(controls):
        raise ValueError("Original reservoir record count disagrees with frozen manifest")
    decoy_tag = config["database"].get("decoy_tag", "rev_")
    if not decoy_tag or any(decoy_tag in accession for accession in controls):
        raise ValueError("Control accessions conflict with Sage decoy labeling")
    augmented_record = metadata(augmented)
    dump(
        output / "AUGMENTATION.json",
        {
            "reservoir": augmented_record,
            "records": total,
            "control_records": len(controls),
            "entry_stage": "BEFORE_RETRIEVAL_AND_SELECTION",
            "absence_evidence": control_spec["absence_evidence"],
        },
    )
    local = dict(manifest, output=str(output / "acquisitions"))
    dest = case(local, index)
    dest.mkdir(parents=True)
    compiled = compile_queries(task["predictions"]["path"], manifest["alphanovo_source"], dest)
    scan(manifest["bag_engine"]["path"], dest / "queries.jsonl", augmented, dest)
    scanned = json.loads((dest / "scan.json").read_text())
    if (
        scanned["reservoir_records"] != total
        or scanned["reservoir_sha256"] != augmented_record["sha256"]
    ):
        raise ValueError("Incomplete augmented reservoir scan")
    selection = select_databases(dest / "candidates.jsonl.gz", dest, compiled["query_units"])
    stages = {arm: Counter() for arm in ARMS}
    with gzip.open(dest / "candidates.jsonl.gz", "rt") as handle:
        for line in handle:
            record = json.loads(line)
            for arm, mask in ARMS.items():
                if evidence_for(record, mask):
                    label = "entrapment" if record["header"].split()[0] in controls else "target"
                    stages[arm]["retrieved_" + label] += 1
    for arm in ARMS:
        for header, _ in fasta(dest / (arm + ".fasta")):
            label = "entrapment" if header.split()[0] in controls else "target"
            stages[arm]["selected_" + label] += 1
    dump(
        dest / "CANDIDATES_COMPLETE.json",
        {
            "status": "CANDIDATES_COMPLETE",
            "selection": selection,
            "comparison": "fresh exact versus relaxed selection from the same augmented reservoir",
        },
    )
    # Uses the unchanged settings, all six arms and original spectra, just as the study does.
    search(local, index)
    rows = {
        arm: read_psms(dest / "searches" / arm / "results.sage.tsv")
        if selection["arms"][arm]["fasta"]["bytes"]
        else []
        for arm in ARMS
    }
    classes = membership(
        {row[0] for values in rows.values() for row in values}, target, entrapment, possible
    )
    with gzip.open(output / "identifications.tsv.gz", "wt") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["arm", "peptide_IL", "filename", "scan", "peptide_q", "spectrum_q", "class"]
        )
        for arm, values in rows.items():
            writer.writerows((arm, *row, classes[row[0]]) for row in values)
    dump(
        output / "ENTRAPMENT_DIAGNOSTIC.json",
        {
            "status": "DIAGNOSTIC_COMPLETE_REVIEW_REQUIRED",
            "sample": task["sample"],
            "workflow_fdr": "UNVERIFIED",
            "same_fdr_sensitivity_claim_supported": False,
            "interpretation": "Control fractions are not target-only FDR estimates or upper bounds",
            "protein_group_fdr": "NOT_EVALUATED",
            "absence_record": absence,
            "thresholds": THRESHOLDS,
            "stage_counts": stages,
            "arms": {
                arm: summarize(values, classes, rows["exact"]) for arm, values in rows.items()
            },
            "job_id": os.environ.get("SLURM_JOB_ID"),
        },
    )


def main():
    """Require a frozen study, control specification, acquisition index and new output."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("index", type=int)
    parser.add_argument("controls", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    run(
        json.loads(args.manifest.read_text()),
        args.index,
        json.loads(args.controls.read_text()),
        args.output.resolve(),
    )


if __name__ == "__main__":
    main()
