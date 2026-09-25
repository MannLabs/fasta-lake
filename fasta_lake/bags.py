"""Experimental protein retrieval using pinned AlphaNovo reference-explanation rules.

The matching rules come directly from MannLabs/alphanovo, commit
93d647cf1d996864c3367d5f04e880b5df08d977. Unlike its FASTA-extension caller,
we retain every compatible protein and preserve the original prediction as
one selection-evidence unit. This is not a novel-peptide discovery workflow.
"""

import csv
import gzip
import hashlib
import importlib
import json
import math
import random
import re
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

ALPHANOVO_REVISION = "93d647cf1d996864c3367d5f04e880b5df08d977"
ARMS = {"exact": 1, "bags": 3, "mass": 5, "combined": 7, "shuffle_1": 56, "shuffle_2": 448}
RESIDUES = set("ACDEFGHIKLMNPQRSTVWY")


def norm(sequence):
    """Canonicalize residue case and collapse indistinguishable I/L spellings."""
    return sequence.upper().replace("I", "L")


def sha(path):
    """Hash a complete file while bounding read-buffer memory."""
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def metadata(path):
    """Record the absolute path, byte size and SHA-256 of an input or output."""
    path = Path(path).resolve()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha(path)}


def dump(path, value):
    """Atomically replace a JSON receipt, rejecting non-finite numbers."""
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def upstream(source):
    """Load the actual pinned implementation; never select another installed copy."""
    source = Path(source).resolve()
    revision = subprocess.check_output(
        ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
    ).strip()
    if revision != ALPHANOVO_REVISION:
        raise ValueError(f"Expected AlphaNovo {ALPHANOVO_REVISION}, got {revision}")
    dirty = subprocess.check_output(
        ["git", "-C", str(source), "diff", "HEAD", "--", "alphanovo"], text=True
    )
    if dirty:
        raise ValueError("Pinned AlphaNovo source has tracked changes")
    sys.path.insert(0, str(source))
    mass = importlib.import_module("alphanovo.utils.mass_equivalent")
    bag = importlib.import_module("alphanovo.utils.anchor_bag_matcher")
    for module in (mass, bag):
        if source not in Path(module.__file__).resolve().parents:
            raise ValueError("Another AlphaNovo installation was imported")
    return mass, bag


def low_complexity(sequence):
    """Apply the archived composition and alternating-repeat exclusion rules."""
    n = len(sequence)
    counts = Counter(sequence)
    return (
        max(counts.values()) / n > 0.6
        or (n >= 12 and len(counts) <= 2)
        or (n >= 10 and all(sequence[i] == sequence[i % 2] for i in range(n)))
    )


def load_queries(path):
    """Keep original row ranking/ties; collapse selection evidence by original sequence."""
    rows = []
    rejected = Counter()
    with open(path) as handle:
        for number, row in enumerate(csv.DictReader(handle), 1):
            sequence = row["peptide_prediction_detokenized_unmodified"].strip()
            try:
                score = float(row["score"])
            except ValueError:
                rejected["score"] += 1
                continue
            if not math.isfinite(score) or score <= 0 or not 9 <= len(sequence) <= 50:
                rejected["length_or_score"] += 1
                continue
            if set(sequence) - RESIDUES:
                rejected["unsupported_sequence"] += 1
                continue
            if low_complexity(sequence):
                rejected["complexity"] += 1
                continue
            if int(row["beam_rank"]) != 0 or not row["spec_idx"]:
                raise ValueError(
                    "This study requires archived rank-zero predictions and spectrum IDs"
                )
            rows.append((number, norm(sequence), score, row))
    if not rows:
        raise ValueError("No eligible predictions")
    scores = sorted((row[2] for row in rows), reverse=True)
    cutoff = scores[math.ceil(len(scores) * 0.3) - 1]
    queries = {}
    for number, sequence, score, row in rows:
        q = queries.setdefault(
            sequence, {"sequence": sequence, "lookup": False, "rows": [], "forms": set()}
        )
        q["lookup"] |= score >= cutoff
        q["rows"].append({"row": number, "spectrum": row["spec_idx"], "score": score})
        if score >= cutoff:
            q["forms"].add(
                (
                    norm(row["sequence_with_bags"].strip()),
                    row["peptide_prediction_detokenized"].strip(),
                )
            )
    return list(queries.values()), {
        "eligible_rows": len(rows),
        "retained_rows": sum(row[2] >= cutoff for row in rows),
        "cutoff": cutoff,
        "rejected": dict(rejected),
    }


def bag_plan(pattern, sequence, native):
    """Validate recorded spans and choose an anchor for all-hit streaming lookup."""
    parts = native.parse_pattern(pattern)
    elements, anchors, pos = [], [], 0
    for part in parts:
        unordered = not isinstance(part, str)
        letters = norm("".join(sorted(part.elements())) if unordered else part)
        observed = sequence[pos : pos + len(letters)]
        if (Counter(letters) != Counter(observed)) if unordered else (letters != observed):
            raise ValueError("Recorded bag does not describe its original prediction")
        elements.append({"letters": letters, "unordered": unordered})
        if not unordered:
            anchors.append((len(letters), letters, pos))
        pos += len(letters)
    if pos != len(sequence):
        raise ValueError("Recorded bag length mismatch")
    if not anchors:
        return None
    longest = max(anchors, key=lambda item: item[0])
    if longest[0] < 4 and (longest[0] < 3 or sum(item[0] for item in anchors) < 4):
        return None
    return {"seed": longest[1], "offset": longest[2], "elements": elements}


def tokens(modified, sequence):
    """Separate residue and terminal tokens and check their unmodified spelling."""
    units = re.findall(r"\[[^\]]+\]|[A-Z]", modified)
    if "".join(units) != modified:
        raise ValueError("Unsupported modified prediction syntax")
    bare, residues, prefix, suffix = [], [], [], []
    for unit in units:
        if unit.startswith("["):
            site = unit.rsplit("@", 1)[1].rstrip("]").split("^")[0]
            if site in {"Any_N-term", "Protein_N-term"}:
                if residues:
                    raise ValueError("N-terminal modification follows a residue")
                prefix.append(unit)
                continue
            if site in {"Any_C-term", "Protein_C-term"}:
                suffix.append(unit)
                continue
            if len(site) != 1 or site not in RESIDUES:
                raise ValueError("Unrecognized modification location")
            bare.append(site)
        else:
            bare.append(unit)
        if suffix:
            raise ValueError("Residue follows a C-terminal modification")
        residues.append(unit)
    if norm("".join(bare)) != sequence:
        raise ValueError("Modified and unmodified predictions disagree")
    return residues, "".join(prefix), "".join(suffix)


def control_form(sequence, pattern, modified, seed, native):
    """Shuffle residues reproducibly while preserving bag spans and modifications."""
    permutation = list(range(len(sequence)))
    random.Random(f"{seed}:{sequence}").shuffle(permutation)
    shuffled = "".join(sequence[i] for i in permutation)
    units, prefix, suffix = tokens(modified, sequence)
    shuffled_modified = prefix + "".join(units[i] for i in permutation) + suffix
    pieces, pos = [], 0
    for part in native.parse_pattern(pattern or sequence):
        unordered = not isinstance(part, str)
        length = sum(part.values()) if unordered else len(part)
        piece = shuffled[pos : pos + length]
        pieces.append("[" + "".join(sorted(piece)) + "]" if unordered else piece)
        pos += length
    return shuffled, "".join(pieces), shuffled_modified


def compile_queries(predictions, source, output, max_patterns=40000000):
    """Compile native AlphaNovo alternatives without turning variants into observations.

    The paper's metaproteome setting omits D->N and E->Q modification reversions.
    Expanded lookup uses the original top-30% population. Below-cutoff eligible
    predictions retain their existing exact selection support only.
    """
    started = time.monotonic()
    mass, native = upstream(source)
    table = {
        key: tuple(value for value in values if (key, value) not in {("D", "N"), ("E", "Q")})
        for key, values in mass.build_swap_table().items()
    }
    queries, receipt = load_queries(predictions)
    output = Path(output)
    counts = Counter()
    with (
        open(output / "queries.jsonl", "w") as handle,
        gzip.open(output / "evidence.jsonl.gz", "wt") as evidence,
    ):
        for index, original in enumerate(queries):
            for control, seed in enumerate((None, 20260923, 20260924)):
                uid = index * 3 + control
                forms = original["forms"] or {(original["sequence"], original["sequence"])}
                transformed = []
                for pattern, modified in sorted(forms):
                    tokens(modified, original["sequence"])
                    transformed.append(
                        (original["sequence"], pattern or original["sequence"], modified)
                        if seed is None
                        else control_form(original["sequence"], pattern, modified, seed, native)
                    )
                sequence = transformed[0][0]
                assert all(form[0] == sequence for form in transformed)
                evidence.write(
                    json.dumps(
                        {
                            "query": uid,
                            "original": original["sequence"],
                            "control": control,
                            "sequence": sequence,
                            "lookup": original["lookup"],
                            "rows": original["rows"],
                        }
                    )
                    + "\n"
                )
                shift = 3 * control

                def write(plan, route):
                    """Emit one route-tagged pattern, failing before truncating a query."""
                    counts["patterns"] += 1
                    if counts["patterns"] > max_patterns:
                        raise RuntimeError("Pattern guard exceeded; incomplete run, never truncate")
                    handle.write(
                        json.dumps(
                            dict(plan, query=uid, route=route << shift, lookup=original["lookup"]),
                            separators=(",", ":"),
                        )
                        + "\n"
                    )

                write({"seed": sequence}, 1)
                if original["lookup"]:
                    bags_seen = set()
                    for _, pattern, _ in transformed:
                        if pattern in bags_seen:
                            continue
                        bags_seen.add(pattern)
                        plan = bag_plan(pattern, sequence, native)
                        if plan is not None:
                            write(plan, 2)
                        else:
                            counts["unsearchable_bags"] += 1
                    variants = mass.generate_variants(
                        sequence.replace("L", "I"), table, max_substitutions=2
                    )
                    for _, _, modified in transformed:
                        if "[Carbamidomethyl@C]" in modified:
                            variants.update(
                                mass.generate_modified_mass_equivalent_variants(modified)
                            )
                    for variant in sorted({norm(v) for v in variants} - {sequence}):
                        write({"seed": variant}, 4)
                if index % 1000 == 0 and control == 0:
                    dump(
                        output / "compile_progress.json",
                        {
                            "queries_processed": index,
                            "queries_total": len(queries),
                            "patterns": counts["patterns"],
                        },
                    )
    receipt.update(counts)
    receipt.update(
        query_units=len(queries) * 3,
        original_sequences=len(queries),
        elapsed_seconds=time.monotonic() - started,
        alphanovo_revision=ALPHANOVO_REVISION,
        mass_tolerance=mass.MASS_TOL,
        max_substitutions=2,
        deamidation_reversions=False,
        alternative_spelling_is_additional_evidence=False,
        expanded_support_below_lookup_cutoff=False,
    )
    dump(output / "compile.json", receipt)
    return receipt


def fnv1a(value):
    """Compute the production unsigned 64-bit FNV-1a accession tie-breaker."""
    h = 14695981039346656037
    for byte in value.encode():
        h = ((h ^ byte) * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return h


def evidence_for(record, mask):
    """Return original-query support only when the arm admits this protein."""
    if not record["lookup_routes"] & mask:
        return []
    return [query for query, routes in record["evidence"] if routes & mask]


def select_databases(candidates, output, query_count, arms=None):
    """Apply FastaLake's fixed-count razor rule to original-query/protein edges.

    Stream candidate sequences rather than holding the reservoir in memory.
    Multiple variant hits in one protein are one original-query edge.
    """
    arms = ARMS if arms is None else arms
    started = time.monotonic()
    counts = {arm: [0] * query_count for arm in arms}
    candidates_per_arm = Counter()
    with gzip.open(candidates, "rt") as handle:
        for line in handle:
            record = json.loads(line)
            for arm, mask in arms.items():
                ids = evidence_for(record, mask)
                if ids:
                    candidates_per_arm[arm] += 1
                for query in ids:
                    counts[arm][query] += 1
    best = {arm: {} for arm in arms}
    with gzip.open(candidates, "rt") as handle:
        for line in handle:
            record = json.loads(line)
            accession = record["header"].split()[0]
            for arm, mask in arms.items():
                ids = evidence_for(record, mask)
                unique = sum(counts[arm][query] == 1 for query in ids)
                key = (-unique, -len(ids), fnv1a(accession), accession)
                for query in ids:
                    if query not in best[arm] or key < best[arm][query]:
                        best[arm][query] = key
    chosen = {arm: {key[-1] for key in values.values()} for arm, values in best.items()}
    handles = {arm: open(Path(output) / (arm + ".fasta"), "w") for arm in arms}
    seen = {arm: set() for arm in arms}
    try:
        with gzip.open(candidates, "rt") as handle:
            for line in handle:
                record = json.loads(line)
                accession = record["header"].split()[0]
                for arm in arms:
                    if accession in chosen[arm]:
                        if accession in seen[arm]:
                            raise ValueError("Duplicate selected accession in reservoir")
                        seen[arm].add(accession)
                        handles[arm].write(
                            ">" + record["header"] + "\n" + record["sequence"] + "\n"
                        )
    finally:
        for handle in handles.values():
            handle.close()
    assert seen == chosen
    dump(
        Path(output) / "assignments.json",
        {
            arm: {str(query): key[-1] for query, key in values.items()}
            for arm, values in best.items()
        },
    )
    receipt = {
        "elapsed_seconds": time.monotonic() - started,
        "arms": {
            arm: {
                "candidate_proteins": candidates_per_arm[arm],
                "selected_proteins": len(chosen[arm]),
                "assigned_original_queries": len(best[arm]),
                "fasta": metadata(Path(output) / (arm + ".fasta")),
            }
            for arm in arms
        },
    }
    dump(Path(output) / "selection.json", receipt)
    return receipt
