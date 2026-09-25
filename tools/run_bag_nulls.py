#!/usr/bin/env python3
"""Compare matched exact/bag/mass nulls from an existing full, route-tagged scan.

Every shuffled query was already scanned by the native experiment under exact,
bag and mass routes. This adds the missing separate null selection/search arms
without recompiling queries, reshuffling them or scanning the reservoir again.
"""

import argparse
import csv
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fasta_lake.bags import ARMS, dump, evidence_for, metadata, select_databases
from tools.run_bag_entrapment import THRESHOLDS, read_psms
from tools.run_bags import case, search, verify

NULL_ARMS = {
    f"shuffle_{seed}_{method}": mask << (seed * 3)
    for seed in (1, 2)
    for method, mask in (("exact", 1), ("bags", 3), ("mass", 5))
}
ALL_ARMS = dict(ARMS, **NULL_ARMS)


def retrieval(candidates, evidence, output):
    """Count query hits and protein multiplicity with each full queried denominator."""
    units = {}
    with gzip.open(evidence, "rt") as handle:
        for line in handle:
            row = json.loads(line)
            units[row["query"]] = row
    counts = {arm: Counter() for arm in ALL_ARMS}
    proteins = Counter()
    with gzip.open(candidates, "rt") as handle:
        for line in handle:
            record = json.loads(line)
            for arm, mask in ALL_ARMS.items():
                ids = evidence_for(record, mask)
                if ids:
                    proteins[arm] += 1
                for uid in ids:
                    if units[uid]["lookup"]:
                        counts[arm][uid] += 1
    rows = []
    with gzip.open(output / "NULL_QUERY_HITS.tsv.gz", "wt") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["arm", "query", "length", "candidate_proteins"])
        for arm, counter in counts.items():
            seed = int(arm.split("_")[1]) if arm.startswith("shuffle_") else 0
            by_length = {}
            for uid, unit in units.items():
                if unit["control"] != seed or not unit["lookup"]:
                    continue
                length = len(unit["sequence"])
                n = counter[uid]
                writer.writerow([arm, uid, length, n])
                row = by_length.setdefault(length, Counter())
                row["queried"] += 1
                row["matched"] += n > 0
                row["query_protein_edges"] += n
            for length, value in sorted(by_length.items()):
                rows.append(
                    dict(
                        arm=arm,
                        length=length,
                        **value,
                        hit_fraction=value["matched"] / value["queried"],
                    )
                )
    comparisons = []
    for arm, counter in counts.items():
        seed = int(arm.split("_")[1]) if arm.startswith("shuffle_") else 0
        base_arm = f"shuffle_{seed}_exact" if seed else "exact"
        queried = {uid for uid, unit in units.items() if unit["control"] == seed and unit["lookup"]}
        matched = {uid for uid in queried if counter[uid]}
        exact = {uid for uid in queried if counts[base_arm][uid]}
        comparisons.append(
            dict(
                arm=arm,
                reference=base_arm,
                queried=len(queried),
                matched=len(matched),
                candidate_proteins=proteins[arm],
                hit_fraction=len(matched) / len(queried) if queried else None,
                exact_matched=len(exact),
                gained=len(matched - exact),
                lost=len(exact - matched),
                ratio_to_exact=len(matched) / len(exact) if exact else None,
            )
        )
    return dict(
        by_length=rows,
        paired_comparisons=comparisons,
        candidate_proteins=dict(proteins),
        unit="original unique IL-normalized prediction, paired across matching rules",
    )


def candidate(manifest, index, source_manifest):
    """Reuse only a completed baseline-verified scan; select six missing null arms."""
    source = case(source_manifest, index)
    receipt = json.loads((source / "CANDIDATES_COMPLETE.json").read_text())
    if receipt["status"] != "CANDIDATES_COMPLETE" or receipt["baseline_reproduced"] is not True:
        raise ValueError("Source candidate run is not complete and baseline-verified")
    dest = case(manifest, index)
    dest.mkdir(parents=True, exist_ok=False)
    files = {
        name: metadata(source / name)
        for name in (
            "CANDIDATES_COMPLETE.json",
            "candidates.jsonl.gz",
            "evidence.jsonl.gz",
            "scan.json",
        )
    }
    counts = retrieval(source / "candidates.jsonl.gz", source / "evidence.jsonl.gz", dest)
    selected = select_databases(
        source / "candidates.jsonl.gz", dest, receipt["compiled"]["query_units"], arms=NULL_ARMS
    )
    dump(dest / "NULL_RETRIEVAL.json", counts)
    dump(
        dest / "CANDIDATES_COMPLETE.json",
        {
            "status": "CANDIDATES_COMPLETE",
            "selection": selected,
            "source_outputs": files,
            "source_scan": receipt["scan"],
            "workflow_fdr": "UNVERIFIED",
        },
    )


def report(manifest, index, source_manifest):
    """Compare each seed's exact, bag, mass and combined retrieval and final searches."""
    dest, source = case(manifest, index), case(source_manifest, index)
    base_receipt = json.loads((source / "SEARCHES_COMPLETE.json").read_text())
    added_receipt = json.loads((dest / "SEARCHES_COMPLETE.json").read_text())
    if any(r["status"] != "MATCHED_SEARCHES_COMPLETE" for r in (base_receipt, added_receipt)):
        raise ValueError("Both original and added null searches must be complete")
    base_selection = json.loads((source / "selection.json").read_text())["arms"]
    new_selection = json.loads((dest / "selection.json").read_text())["arms"]
    selection = dict(base_selection, **new_selection)
    values, provenance = {}, {}
    for arm in ALL_ARMS:
        root = dest if arm in NULL_ARMS else source
        if selection[arm]["fasta"]["bytes"]:
            path = root / "searches" / arm / "results.sage.tsv"
            # Empty/database-reuse receipts have no result hash of their own.
            receipt = added_receipt if arm in NULL_ARMS else base_receipt
            execution = receipt["executions"][arm]
            if "reused_identical_fasta" in execution:
                execution = receipt["executions"][execution["reused_identical_fasta"]]
            verify(execution["results"])
            values[arm] = read_psms(path)
            provenance[arm] = metadata(path)
        else:
            values[arm] = []
    rows = []
    for seed in (1, 2):
        for q in THRESHOLDS:
            for level, col in (("canonical_peptide", 3), ("PSM", 4)):

                def identities(arm):
                    """Keep endpoint-specific q-values and spectrum-peptide assignments."""
                    return {r[0] if col == 3 else r[:3] for r in values[arm] if r[col] <= q}

                exact = identities(f"shuffle_{seed}_exact")
                for method in ("exact", "bags", "mass", "combined"):
                    arm = f"shuffle_{seed}" if method == "combined" else f"shuffle_{seed}_{method}"
                    current = identities(arm)
                    rows.append(
                        dict(
                            seed=seed,
                            method=method,
                            q=q,
                            level=level,
                            accepted=len(current),
                            exact_accepted=len(exact),
                            gained=len(current - exact),
                            lost=len(exact - current),
                            net=len(current) - len(exact),
                            ratio_to_exact=len(current) / len(exact) if exact else None,
                            selected_proteins=selection[arm]["selected_proteins"],
                        )
                    )
    dump(
        dest / "NULL_COMPARISON.json",
        {
            "status": "NULL_DIAGNOSTIC_COMPLETE_REVIEW_REQUIRED",
            "rows": rows,
            "results": provenance,
            "workflow_fdr": "UNVERIFIED",
            "interpretation": (
                "Shuffled-query selection can retrieve real proteins; "
                "accepted hits are not known false IDs"
            ),
            "same_fdr_sensitivity_claim_supported": False,
        },
    )


def main():
    """Run one scheduled stage from a frozen manifest, preserving existing outputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("stage", choices=["candidate", "search", "report"])
    parser.add_argument("index", type=int)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    for record in [manifest["source_manifest"], *manifest["code"]]:
        verify(record)
    source = json.loads(Path(manifest["source_manifest"]["path"]).read_text())
    if args.stage == "candidate":
        candidate(manifest, args.index, source)
    elif args.stage == "search":
        search(manifest, args.index, arms=NULL_ARMS, baseline_arm="shuffle_1_exact")
    else:
        report(manifest, args.index, source)


if __name__ == "__main__":
    main()
