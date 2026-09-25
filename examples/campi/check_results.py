#!/usr/bin/env python3
"""Check database identities and scientific accounting in the CAMPI example."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fasta_count(path):
    with path.open() as stream:
        return sum(line.startswith(">") for line in stream)


def canonical(sequence):
    return re.sub(r"\[[^\]]*\]", "", sequence).replace("I", "L")


def inspect(run):
    workflow = run / "workflow"
    if not (workflow / "COMPLETE.json").is_file():
        raise ValueError("Workflow did not finish")
    summary = {
        "reference_sequences": fasta_count(run / "lake.fasta"),
        "reference_sha256": sha(run / "lake.fasta"),
        "evidence_sequences": fasta_count(workflow / "evidence.fasta"),
        "evidence_sha256": sha(workflow / "evidence.fasta"),
        "samples": {},
    }
    peptides = {}
    provenance = json.loads((workflow / "PROVENANCE.json").read_text())
    samples = [row["sample"] for row in provenance["samples"]]
    if len(samples) != 2 or len(set(samples)) != 2:
        raise ValueError("Expected exactly two distinct acquisition IDs")
    for sample in samples:
        drawn = workflow / "samples" / sample / "sample_lake.fasta"
        selected = workflow / "samples" / sample / "inference" / (sample + "_razor.fasta")
        accepted = set()
        psms = 0
        with (workflow / "search" / sample / "results.sage.tsv").open() as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                if row["label"] != "1":
                    continue
                sq, pq = float(row["spectrum_q"]), float(row["peptide_q"])
                if math.isfinite(sq) and 0 <= sq <= 0.01:
                    psms += 1
                if math.isfinite(pq) and 0 <= pq <= 0.01:
                    accepted.add(canonical(row["peptide"]))
        if not accepted or not psms:
            raise ValueError(sample + ": no accepted target identifications")
        peptides[sample] = accepted
        summary["samples"][sample] = {
            "drawn_sequences": fasta_count(drawn),
            "drawn_sha256": sha(drawn),
            "selected_sequences": fasta_count(selected),
            "selected_sha256": sha(selected),
            "target_psms_spectrum_q_001": psms,
            "canonical_target_peptides_peptide_q_001": len(accepted),
        }
    study = json.loads((run / "study" / "summary.json").read_text())
    if study["acquisitions"] != 2:
        raise ValueError("Expected exactly two completed study acquisitions")
    if not math.isclose(study["intensity_in"], study["intensity_out"], rel_tol=1e-12):
        raise ValueError("Study assignment did not conserve included intensity")
    if study["intensity_in"] <= 0 or study["quantified_study_groups"] <= 0:
        raise ValueError("No positive study quantification")
    with (run / "study" / "study_group_long.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows:
        raise ValueError("Study reporting table is empty")
    if {row["sample"] for row in rows} != set(samples):
        raise ValueError("Study reporting must contain both declared acquisitions")
    if any(not row["study_group"] for row in rows):
        raise ValueError("Study reporting contains a missing counting key")
    if any(not math.isfinite(float(row["intensity"])) for row in rows):
        raise ValueError("Study reporting contains a non-finite intensity")
    if not math.isclose(
        math.fsum(float(row["intensity"]) for row in rows), study["intensity_out"], rel_tol=1e-12
    ):
        raise ValueError("Study table intensity differs from its summary")
    if len({(row["sample"], row["study_group"]) for row in rows}) != len(rows):
        raise ValueError("Repeated acquisition/study-group row")
    summary["study"] = study
    summary["study_group_rows"] = len(rows)
    return summary, peptides


def validate(run, expected_dir):
    actual, peptides = inspect(run)
    expected = json.loads((expected_dir / "results.json").read_text())
    baseline = expected["measurements"]
    if set(actual["samples"]) != set(baseline["samples"]):
        raise ValueError("Observed acquisition IDs differ from the expected example")
    errors = []
    for key in ("reference_sequences", "reference_sha256", "evidence_sequences", "evidence_sha256"):
        if actual[key] != baseline[key]:
            errors.append(key + " differs from the expected database")
    similarities = {}
    for sample, result in actual["samples"].items():
        reference = baseline["samples"][sample]
        for key in ("drawn_sequences", "drawn_sha256", "selected_sequences", "selected_sha256"):
            if result[key] != reference[key]:
                errors.append(sample + ": " + key + " differs")
        for key in ("target_psms_spectrum_q_001", "canonical_target_peptides_peptide_q_001"):
            tolerance = max(2, reference[key] * expected["search_count_relative_tolerance"])
            if abs(result[key] - reference[key]) > tolerance:
                errors.append(sample + ": " + key + " is outside tolerance")
        reference_peptides = set(
            (expected_dir / (sample + "_canonical_peptides.txt")).read_text().splitlines()
        )
        union = peptides[sample] | reference_peptides
        similarity = len(peptides[sample] & reference_peptides) / len(union)
        similarities[sample] = similarity
        if similarity < expected["minimum_peptide_jaccard"]:
            errors.append(sample + ": accepted peptide identities differ substantially")
    study_tolerance = max(
        2, baseline["study"]["study_groups"] * expected["search_count_relative_tolerance"]
    )
    inferred_count = actual["study"]["inferred_study_groups"]
    if abs(inferred_count - baseline["study"]["study_groups"]) > study_tolerance:
        errors.append("Study group count is outside tolerance")
    result = {
        "status": "passed" if not errors else "failed",
        "errors": errors,
        "peptide_jaccard_to_reference": similarities,
        "measurements": actual,
        "interpretation": (
            "Software example acceptance, not independent FDR calibration "
            "or a full-cohort benchmark."
        ),
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument(
        "--expected", type=Path, default=Path(__file__).resolve().parent / "expected"
    )
    args = parser.parse_args()
    result = validate(args.run, args.expected)
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
