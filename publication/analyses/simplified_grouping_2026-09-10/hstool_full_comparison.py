#!/usr/bin/env python3
"""Paired HSTOOL comparison on exactly 58 declared acquisitions.

Rebuild the original study rule on the same evidence as direct peptide assignment.
Inventory the retained upstream search/selection outputs, then compare assigned
peptide support, coverage, quantification, reproducibility and representative
annotation. No new database selection or Sage search is performed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import itertools
import json
import math
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from fasta_lake.study import _read_fasta, canonical_peptide  # noqa: E402

TIERS = ["TechnicalReplicates", "WorkflowReplicates", "OnePiece", "Intra", "Interindividual"]
EXPECTED = dict(zip(TIERS, [12, 15, 10, 11, 10]))
METHODS = ["original", "simplified"]


def read_rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        yield from csv.DictReader(f, delimiter="\t")


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1048576), b""):
            h.update(b)
    return h.hexdigest()


def write_rows(path, rows):
    if not rows:
        raise ValueError(f"Empty result table: {path}")
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def median(values):
    values = list(values)
    return float(statistics.median(values)) if values else None


def passing(q):
    q = float(q)
    return math.isfinite(q) and 0 <= q <= 0.01


def coverage(sequence, peptides):
    """Union covered residues over all exact I/L-equivalent occurrences.

    Unsupported peptides add no coverage and are counted separately.
    """
    sequence = sequence.replace("I", "L")
    if not sequence:
        raise ValueError("Cannot calculate coverage of an empty sequence")
    intervals, unsupported = [], 0
    for peptide in peptides:
        if not peptide:
            raise ValueError("Empty peptide cannot supply sequence coverage")
        peptide = peptide.replace("I", "L")
        start = sequence.find(peptide)
        if start < 0:
            unsupported += 1
            continue
        while start >= 0:
            intervals.append((start, start + len(peptide)))
            start = sequence.find(peptide, start + 1)
    end, covered = 0, 0
    for start, stop in sorted(intervals):
        covered += max(0, stop - max(start, end))
        end = max(end, stop)
    return 100 * covered / len(sequence), unsupported


def correlation(a, b):
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def group_class(classes):
    return next(iter(classes)) if len(classes) == 1 else "Mixed"


def metadata_key(key):
    parts = []
    for p in key.split(";"):
        tagged = p.split("|", 2)
        if len(tagged) == 3 and all(re.fullmatch(r"[A-Z0-9_]+", x) for x in tagged[:2]):
            p = tagged[2]
        parts.append(p)
    return ";".join(sorted(set(parts)))


def annotation_maps(eggnog, genomes):
    kos = {}
    for r in read_rows(eggnog):
        kos[r["accession"]] = {
            k.removeprefix("ko:") for k in r["KEGG_ko"].split(",") if k.startswith("ko:K")
        }
    records = list(read_rows(genomes))
    fallback = {r["Species_rep"]: r["Lineage"] for r in records if r["Lineage"]}
    species = {}
    for r in records:
        lineage = r["Lineage"] or fallback.get(r["Species_rep"], "")
        names = [x[3:] for x in lineage.split(";") if x.startswith("s__") and len(x) > 3]
        if names:
            species[r["Genome"]] = names[-1]
    return kos, species


def parse_searches(base, names, design, registry):
    sequences, graph, local_graphs, fixed = {}, defaultdict(set), {}, []
    for number, s in enumerate(names, 1):
        directory = base / "cohorts/HSTOOL/5_sage/razor_hash-acc" / s
        config = directory / "results.json"
        raw = config.read_bytes()
        registry.append({"path": str(config), "sha256": hashlib.sha256(raw).hexdigest()})
        fasta = Path(json.loads(raw)["database"]["fasta"])
        members = _read_fasta(fasta, sequences, registry)
        local, identified = defaultdict(set), defaultdict(set)
        peptides, raw_peptides = set(), set()
        target_q, decoy_q, accepted_psms = 0, 0, 0
        matched, fractions = [], []
        path = directory / "results.sage.tsv"
        digest = hashlib.sha256()
        with path.open("rb") as f:
            header = f.readline()
            digest.update(header)
            fields = header.decode().rstrip().split("\t")
            idx = {
                k: fields.index(k)
                for k in [
                    "peptide",
                    "proteins",
                    "label",
                    "peptide_q",
                    "protein_q",
                    "spectrum_q",
                    "matched_peaks",
                    "matched_intensity_pct",
                ]
            }
            for raw in f:
                digest.update(raw)
                cells = raw.decode().rstrip("\r\n").split("\t")
                label = cells[idx["label"]]
                if passing(cells[idx["spectrum_q"]]):
                    target_q += label == "1"
                    decoy_q += label == "-1"
                if label != "1" or not passing(cells[idx["peptide_q"]]):
                    continue
                peptide = sys.intern(canonical_peptide(cells[idx["peptide"]]))
                ps = {
                    sys.intern(p)
                    for p in cells[idx["proteins"]].split(";")
                    if p and not p.startswith("rev_")
                }
                if not ps <= members:
                    raise ValueError(f"Out-of-database protein in {s}")
                peptides.add(peptide)
                raw_peptides.add(cells[idx["peptide"]])
                accepted_psms += 1
                matched.append(float(cells[idx["matched_peaks"]]))
                fractions.append(float(cells[idx["matched_intensity_pct"]]))
                for p in ps:
                    graph[p].add(peptide)
                    local[p].add(peptide)
                    if passing(cells[idx["protein_q"]]):
                        identified[p].add(peptide)
        registry.append({"path": str(path), "sha256": digest.hexdigest()})
        local_graphs[s] = local
        core = base / "bench/HSTOOL/core" / (s + ".json")
        c = json.loads(core.read_text())
        metrics = {r["metric"]: r["value"] for r in c["metrics"]}
        registry.append({"path": str(core), "sha256": sha256(core)})
        selection_path = base / "cohorts/HSTOOL/3_persample" / s / "sample_lake.stats.json"
        selection = json.loads(selection_path.read_text())
        registry.append({"path": str(selection_path), "sha256": sha256(selection_path)})
        if metrics["db_entries"] != len(members):
            raise ValueError(f"Archived database count does not match actual FASTA in {s}")
        fixed.append(
            dict(
                sample=s,
                tier=design[s]["group"],
                database_sequences=len(members),
                database_residues=sum(len(sequences[p]) for p in members),
                selection_evidence_peptides=metrics["evidence_peptides"],
                candidate_lake_proteins_searched=metrics["lake_proteins_searched"],
                proteins_matching_selection_evidence=selection["proteins_matched"],
                pre_razor_selected_sequences=selection["total_output_proteins"],
                accepted_canonical_peptides=len(peptides),
                accepted_modified_peptide_strings=len(raw_peptides),
                accepted_target_psms_at_peptide_q=accepted_psms,
                target_psms_at_spectrum_q=target_q,
                decoy_psms_at_spectrum_q=decoy_q,
                median_matched_peaks_accepted_psms=median(matched),
                median_matched_intensity_pct_accepted_psms=median(fractions),
                identified_accessions_at_peptide_and_protein_q=len(identified),
                identified_fraction_of_search_database_pct=100 * len(identified) / len(members),
                median_observed_peptides_per_identified_accession=median(
                    map(len, identified.values())
                ),
                mean_observed_peptides_per_identified_accession=statistics.mean(
                    map(len, identified.values())
                ),
                median_observed_coverage_pct_per_identified_accession=median(
                    coverage(sequences[p], peps)[0] for p, peps in identified.items()
                ),
                archived_local_parsimony_groups=metrics["pg_parsimony"],
                archived_local_parsimony_groups_at_least_two_assigned_peptides=metrics[
                    "pg_parsimony_2pep"
                ],
            )
        )
        if number % 10 == 0 or number == len(names):
            print(
                f"Read accepted evidence for {number}/{len(names)} HSTOOL acquisitions", flush=True
            )
    return sequences, graph, local_graphs, fixed


def quantify(args, names, design, sequences, local_graphs, dictionary, class_of, kos, species):
    group_matrix = {m: defaultdict(dict) for m in METHODS}
    classes = {m: defaultdict(set) for m in METHODS}
    ko_matrix = {m: defaultdict(dict) for m in METHODS}
    ko_by_class = {m: defaultdict(dict) for m in METHODS}
    summaries, function_rows, moves, composition = [], [], [], []
    annotation_qc, group_details = [], args.out / "group_evidence.tsv.gz"
    stream = read_rows(args.direct / "feature_assignments.tsv.gz")
    observed_samples = []
    with gzip.open(group_details, "wt", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "method",
                "sample",
                "study_group",
                "assigned_canonical_peptides",
                "lfq_features",
                "intensity",
                "compatible_assigned_coverage_pct",
                "unsupported_assigned_peptides",
                "locally_observed_representative_peptides",
                "feature_membership_class",
            ]
        )
        for s, feature_rows in itertools.groupby(stream, key=lambda r: r["sample"]):
            if s not in names or s in observed_samples:
                raise ValueError("Unexpected or duplicated acquisition in feature ledger")
            observed_samples.append(s)
            stats = {m: defaultdict(lambda: [set(), 0, 0.0, set(), 0, 0.0]) for m in METHODS}
            total, features, changed, changed_intensity = 0.0, 0, 0, 0.0
            outside_fasta, outside_members = 0, 0
            feature_class_mass = defaultdict(float)
            for r in feature_rows:
                pep = sys.intern(r["canonical_peptide"])
                v = float(r["intensity"])
                outside_fasta += r["representative_in_local_fasta"] == "0"
                outside_members += r["representative_in_local_members"] == "0"
                ms = r["member_set"].split(";")
                candidates = {
                    dictionary["protein_to_group"][p]
                    for p in ms
                    if p in dictionary["protein_to_group"]
                }
                original = min(candidates, key=lambda p: (dictionary["rank"][p], p))
                direct = dictionary["peptide_to_group"][pep]
                if direct != r["study_group"]:
                    raise ValueError("Direct output differs from rebuilt exact-cohort dictionary")
                classification = class_of.get(metadata_key(r["member_set"]), "Unknown")
                feature_class_mass[classification] += v
                total += v
                features += 1
                changed += original != direct
                changed_intensity += v * (original != direct)
                for method, g in [("original", original), ("simplified", direct)]:
                    row = stats[method][g]
                    row[0].add(pep)
                    row[1] += 1
                    row[2] += v
                    row[3].add(classification)
                    incompatible = pep not in sequences[g].replace("I", "L")
                    row[4] += incompatible
                    row[5] += v * incompatible
            for c, value in sorted(feature_class_mass.items()):
                composition.append(
                    dict(
                        sample=s,
                        tier=design[s]["group"],
                        classification=c,
                        intensity=value,
                        intensity_pct=100 * value / total,
                    )
                )
            moves.append(
                dict(
                    sample=s,
                    tier=design[s]["group"],
                    features=features,
                    changed_features=changed,
                    changed_features_pct=100 * changed / features,
                    changed_intensity_pct=100 * changed_intensity / total,
                    direct_representative_outside_local_fasta_features=outside_fasta,
                    direct_representative_outside_local_fasta_features_pct=(
                        100 * outside_fasta / features
                    ),
                    direct_representative_outside_local_members_features=outside_members,
                    direct_representative_outside_local_members_features_pct=(
                        100 * outside_members / features
                    ),
                )
            )
            for method in METHODS:
                supports, coverages, compatible_supports, local_supports = [], [], [], []
                missing_seq_features, missing_seq_intensity = 0, 0.0
                functions = defaultdict(lambda: defaultdict(float))
                resolved_mass, ko_mass, joint_mass = 0.0, 0.0, 0.0
                for g, (peps, n, v, cs, incompatible_n, incompatible_v) in stats[method].items():
                    group_matrix[method][s][g] = v
                    classes[method][g].update(cs)
                    cov, unsupported = coverage(sequences[g], peps)
                    supports.append(len(peps))
                    coverages.append(cov)
                    compatible_supports.append(len(peps) - unsupported)
                    local_supports.append(len(local_graphs[s].get(g, set())))
                    missing_seq_features += incompatible_n
                    missing_seq_intensity += incompatible_v
                    c = group_class(cs)
                    w.writerow(
                        [method, s, g, len(peps), n, v, cov, unsupported, local_supports[-1], c]
                    )
                    matches = re.findall(r"MGYG\d+", g)
                    sp = species.get(matches[0], "Unassigned") if matches else "Unassigned"
                    ks = kos.get(g, set())
                    resolved_mass += v * (sp != "Unassigned")
                    ko_mass += v * bool(ks)
                    joint_mass += v * (sp != "Unassigned" and bool(ks))
                    for k in ks:
                        functions[k][sp] += v / total
                        ko_matrix[method][s][k] = ko_matrix[method][s].get(k, 0.0) + v
                        ck = c + "|" + k
                        ko_by_class[method][s][ck] = ko_by_class[method][s].get(ck, 0.0) + v
                if not math.isclose(sum(group_matrix[method][s].values()), total, rel_tol=1e-10):
                    raise ValueError("Intensity conservation failed")
                if method == "simplified" and missing_seq_features:
                    raise ValueError("Direct assignment contains sequence-incompatible features")
                summaries.append(
                    dict(
                        method=method,
                        sample=s,
                        tier=design[s]["group"],
                        quantified_groups=len(supports),
                        assigned_canonical_peptides_total=sum(supports),
                        lfq_features=features,
                        intensity=total,
                        median_assigned_peptides_per_group=median(supports),
                        mean_assigned_peptides_per_group=statistics.mean(supports),
                        groups_with_at_least_two_assigned_peptides=sum(n >= 2 for n in supports),
                        groups_with_at_least_two_assigned_peptides_pct=100
                        * sum(n >= 2 for n in supports)
                        / len(supports),
                        median_compatible_assigned_peptides_per_group=median(compatible_supports),
                        median_assigned_sequence_coverage_pct=median(coverages),
                        median_locally_observed_representative_peptides=median(local_supports),
                        sequence_incompatible_features=missing_seq_features,
                        sequence_incompatible_feature_pct=100 * missing_seq_features / features,
                        sequence_incompatible_intensity_pct=100 * missing_seq_intensity / total,
                    )
                )
                annotation_qc.append(
                    dict(
                        method=method,
                        sample=s,
                        tier=design[s]["group"],
                        ko_annotated_intensity_pct=100 * ko_mass / total,
                        species_annotated_intensity_pct=100 * resolved_mass / total,
                        jointly_annotated_intensity_pct=100 * joint_mass / total,
                    )
                )
                for k, values in functions.items():
                    ranked = sorted(
                        ((sp, v) for sp, v in values.items() if sp != "Unassigned"),
                        key=lambda t: (-t[1], t[0]),
                    )
                    resolved = sum(v for _, v in ranked)
                    amount = sum(values.values())
                    leader, lead = ranked[0] if ranked else ("", 0.0)
                    second = ranked[1][1] if len(ranked) > 1 else 0.0
                    function_rows.append(
                        dict(
                            method=method,
                            sample=s,
                            tier=design[s]["group"],
                            ko=k,
                            normalised_function_intensity=amount,
                            resolved_fraction=resolved / amount,
                            leading_species=leader,
                            leader_share=lead / resolved if resolved else None,
                            leader_margin=(lead - second) / resolved if resolved else None,
                        )
                    )
    if set(observed_samples) != set(names):
        raise ValueError("Feature ledger does not contain exactly the declared HSTOOL acquisitions")
    emitted = {
        (r["sample"], r["study_group"]): float(r["intensity"])
        for r in read_rows(args.direct / "study_group_long.tsv")
    }
    reconstructed = {
        (s, g): v for s, gs in group_matrix["simplified"].items() for g, v in gs.items()
    }
    if emitted.keys() != reconstructed.keys() or any(
        not math.isclose(v, reconstructed[k], rel_tol=1e-12) for k, v in emitted.items()
    ):
        raise ValueError("Independent feature recount disagrees with emitted direct group table")
    return dict(
        group_matrix=group_matrix,
        classes=classes,
        ko_matrix=ko_matrix,
        ko_by_class=ko_by_class,
        summaries=summaries,
        functions=function_rows,
        moves=moves,
        composition=composition,
        annotation_qc=annotation_qc,
    )


def reproducibility(data, names, design):
    result = []
    for method in METHODS:
        cls = {g: group_class(v) for g, v in data["classes"][method].items()}
        for tier in TIERS:
            for a, b in itertools.combinations([n for n in names if design[n]["group"] == tier], 2):
                for space in ["groups", "KO"]:
                    for c in ["All", "Human", "Microbial"]:
                        if space == "groups":
                            aa, bb = (
                                data["group_matrix"][method][a],
                                data["group_matrix"][method][b],
                            )
                            ka = {k for k in aa if c == "All" or cls[k] == c}
                            kb = {k for k in bb if c == "All" or cls[k] == c}
                        else:
                            source = data["ko_matrix"] if c == "All" else data["ko_by_class"]
                            aa, bb = source[method][a], source[method][b]
                            ka = {k for k in aa if c == "All" or k.startswith(c + "|")}
                            kb = {k for k in bb if c == "All" or k.startswith(c + "|")}
                        shared = sorted(ka & kb)
                        x, y = np.array([aa[k] for k in shared]), np.array([bb[k] for k in shared])
                        result.append(
                            dict(
                                method=method,
                                tier=tier,
                                space=space,
                                classification=c,
                                sample_a=a,
                                sample_b=b,
                                shared_features=len(shared),
                                feature_jaccard=len(shared) / len(ka | kb) if ka | kb else None,
                                pearson_log2=correlation(np.log2(x), np.log2(y)),
                                median_raw_intensity_pair_cv_pct=median(
                                    100 * np.abs(x - y) / (x + y)
                                ),
                            )
                        )
    return result


def export_matrices_and_detection(data, names, design, representatives, out):
    """Keep zero/missing distinct; completeness is descriptive within each tier."""
    result = []
    for method in METHODS:
        matrix = data["group_matrix"][method]
        with (out / f"{method}_study_group_matrix.tsv").open("w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["study_group", *names])
            for group in sorted(representatives):
                writer.writerow([group, *(matrix[s].get(group, "") for s in names)])
        for tier in ["All HSTOOL", *TIERS]:
            samples = [s for s in names if tier == "All HSTOOL" or design[s]["group"] == tier]
            counts = Counter(g for s in samples for g in matrix[s])
            result.append(
                dict(
                    method=method,
                    tier=tier,
                    acquisitions=len(samples),
                    groups_detected_anywhere=len(counts),
                    groups_detected_at_least_half=sum(
                        n >= math.ceil(0.5 * len(samples)) for n in counts.values()
                    ),
                    groups_detected_at_least_80pct=sum(
                        n >= math.ceil(0.8 * len(samples)) for n in counts.values()
                    ),
                    groups_detected_in_all=sum(n == len(samples) for n in counts.values()),
                    mean_completeness_pct_among_detected=100
                    * sum(counts.values())
                    / (len(samples) * len(counts)),
                )
            )
    write_rows(out / "detection_completeness.tsv", result)


def functional_comparison(functions, names, design, top12):
    lookup = {(r["method"], r["sample"], r["ko"]): r for r in functions}
    kos = sorted({r["ko"] for r in functions})
    comparisons, profiles = [], []
    for s in names:
        shared = [k for k in kos if ("original", s, k) in lookup and ("simplified", s, k) in lookup]
        x = np.array([lookup[("original", s, k)]["normalised_function_intensity"] for k in shared])
        y = np.array(
            [lookup[("simplified", s, k)]["normalised_function_intensity"] for k in shared]
        )
        profiles.append(
            dict(
                sample=s,
                tier=design[s]["group"],
                shared_KOs=len(shared),
                pearson_log2_KO_profile=correlation(np.log2(x), np.log2(y)),
            )
        )
        for k in kos:
            a, b = lookup.get(("original", s, k), {}), lookup.get(("simplified", s, k), {})
            la, lb = a.get("leading_species", ""), b.get("leading_species", "")
            comparisons.append(
                dict(
                    sample=s,
                    tier=design[s]["group"],
                    ko=k,
                    displayed_top12=k in top12,
                    original_leader=la,
                    simplified_leader=lb,
                    both_resolved=bool(la and lb),
                    changed_when_both_resolved=(la != lb) if la and lb else None,
                    resolution_status_changed=bool(la) != bool(lb),
                    original_intensity=a.get("normalised_function_intensity", 0.0),
                    simplified_intensity=b.get("normalised_function_intensity", 0.0),
                )
            )
    return comparisons, profiles


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in [
        "base",
        "direct",
        "samples",
        "design",
        "baseline-module",
        "baseline-map",
        "metadata",
        "eggnog",
        "genomes",
        "top12",
        "out",
    ]:
        ap.add_argument("--" + name, type=Path, required=True)
    args = ap.parse_args()
    names = args.samples.read_text().splitlines()
    design = {r["sample_id"]: r for r in read_rows(args.design)}
    if len(names) != len(set(names)) or len(names) != 58:
        raise ValueError("HSTOOL must contain exactly 58 unique acquisitions")
    if Counter(design[n]["group"] for n in names) != EXPECTED:
        raise ValueError("HSTOOL tier membership differs from the declared cohort")
    direct_summary = json.loads((args.direct / "summary.json").read_text())
    if set(direct_summary["acquisition_names"]) != set(names):
        raise ValueError("Direct results were not inferred from the exact HSTOOL cohort")
    if args.out.exists():
        raise FileExistsError(args.out)
    args.out.mkdir(parents=True)
    registry = [
        {"path": str(p), "sha256": sha256(p)}
        for p in [
            args.samples,
            args.design,
            args.baseline_module,
            args.baseline_map,
            args.metadata,
            args.eggnog,
            args.genomes,
            args.top12,
            Path(__file__),
        ]
    ]
    spec = importlib.util.spec_from_file_location("original_study_rule", args.baseline_module)
    baseline = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(baseline)
    sequences, graph, local_graphs, fixed = parse_searches(args.base, names, design, registry)
    dictionary = baseline.study_dictionary(graph, set(sequences))
    direct_groups = {r["study_group"] for r in read_rows(args.direct / "study_groups.tsv")}
    if direct_groups != set(dictionary["rank"]):
        raise ValueError("The two approaches do not share exactly the same inferred groups")
    audited_map = {r["protein"]: r["study_group"] for r in read_rows(args.baseline_map)}
    if audited_map != dictionary["protein_to_group"]:
        raise ValueError("Original protein map differs from independent 58-acquisition audit")
    direct_peptides = {
        r["canonical_peptide"]: r["study_group"]
        for r in read_rows(args.direct / "peptide_to_study_group.tsv")
    }
    if direct_peptides != dictionary["peptide_to_group"]:
        raise ValueError("Approaches differ in their peptide dictionary")
    write_rows(args.out / "fixed_pipeline_metrics.tsv", fixed)
    write_rows(
        args.out / "original_protein_to_group.tsv",
        [dict(protein=p, study_group=g) for p, g in sorted(dictionary["protein_to_group"].items())],
    )
    class_of = {r["protein_group"]: r["classification"] for r in read_rows(args.metadata)}
    kos, species = annotation_maps(args.eggnog, args.genomes)
    data = quantify(
        args, names, design, sequences, local_graphs, dictionary, class_of, kos, species
    )
    del graph, local_graphs, class_of
    export_matrices_and_detection(data, names, design, direct_groups, args.out)
    for filename, key in [
        ("per_acquisition.tsv", "summaries"),
        ("assignment_changes.tsv", "moves"),
        ("feature_composition.tsv", "composition"),
        ("annotation_resolution.tsv", "annotation_qc"),
        ("functional_profiles.tsv", "functions"),
    ]:
        write_rows(args.out / filename, data[key])
    pairs = reproducibility(data, names, design)
    write_rows(args.out / "reproducibility_pairs.tsv", pairs)
    top12 = {r["ko"] for r in read_rows(args.top12)}
    comparisons, profiles = functional_comparison(data["functions"], names, design, top12)
    write_rows(args.out / "functional_leader_comparison.tsv", comparisons)
    write_rows(args.out / "functional_profile_agreement.tsv", profiles)
    primary = [r for r in comparisons if r["displayed_top12"]]
    resolved = [r for r in primary if r["both_resolved"]]
    metrics = [k for k in data["summaries"][0] if k not in {"method", "sample", "tier"}]
    cohort = []
    for tier in ["All HSTOOL", *TIERS]:
        for method in METHODS:
            population = [
                r
                for r in data["summaries"]
                if r["method"] == method and (tier == "All HSTOOL" or r["tier"] == tier)
            ]
            cohort.append(
                dict(
                    method=method,
                    tier=tier,
                    acquisitions=len(population),
                    **{k: median(r[k] for r in population) for k in metrics},
                )
            )
    write_rows(args.out / "cohort_medians.tsv", cohort)
    result = dict(
        cohort="HSTOOL",
        acquisitions=58,
        tier_counts=EXPECTED,
        scope=(
            "Exact cohort applied before both groupings; "
            "identical retained searches and feature population"
        ),
        inferred_groups_in_both_methods=len(dictionary["selected"]),
        checks=dict(
            same_cohort=True,
            same_greedy_dictionary=True,
            original_protein_map_matches_independent_58_acquisition_audit=True,
            same_features=True,
            intensity_conserved=True,
            direct_table_independently_recounted=True,
            direct_representative_explains_every_assigned_peptide=True,
        ),
        cohort_medians=[r for r in cohort if r["tier"] == "All HSTOOL"],
        fixed_pipeline_medians={
            k: median(r[k] for r in fixed) for k in fixed[0] if k not in {"sample", "tier"}
        },
        assignment_changed_features_pct=100
        * sum(r["changed_features"] for r in data["moves"])
        / sum(r["features"] for r in data["moves"]),
        assignment_changed_intensity_pct_median=median(
            r["changed_intensity_pct"] for r in data["moves"]
        ),
        functional_top12=dict(
            cells=len(primary),
            both_resolved_cells=len(resolved),
            changed_leader_cells=sum(r["changed_when_both_resolved"] for r in resolved),
            changed_leader_pct_among_both_resolved=100
            * sum(r["changed_when_both_resolved"] for r in resolved)
            / len(resolved),
            changed_resolution_status_cells=sum(r["resolution_status_changed"] for r in primary),
        ),
        median_KO_profile_correlation=median(
            r["pearson_log2_KO_profile"]
            for r in profiles
            if r["pearson_log2_KO_profile"] is not None
        ),
        annotation_rule=(
            "Same representative KO/UHGG annotation in both approaches; "
            "full group intensity per distinct KO"
        ),
        coverage_rule=(
            "Union of exact I/L-equivalent matched residues; "
            "incompatible peptides counted separately"
        ),
        limitations=[
            "Post-search rerun, no new de novo prediction or Sage searching",
            "Group counts describe quantified study groups, not new peptide identifications",
            "Original local pg_parsimony and accession evidence remain separate fixed endpoints",
            "Representative annotation is still an approximation; no group FDR calibration",
            "Acquisition pairs and technical replicates are dependent; comparisons are descriptive",
        ],
    )
    registry += [
        {"path": str(args.direct / f), "sha256": sha256(args.direct / f)}
        for f in [
            "summary.json",
            "feature_assignments.tsv.gz",
            "study_groups.tsv",
            "study_group_long.tsv",
            "peptide_to_study_group.tsv",
            "input_manifest.json",
        ]
    ]
    (args.out / "PROVENANCE.json").write_text(json.dumps(registry, indent=2) + "\n")
    (args.out / "SUMMARY.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
