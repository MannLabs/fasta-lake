"""Compare two Sage searches by exact proteins, canonical peptides and LFQ.

Reported protein evidence means a sequence is compatible with an accepted PSM;
shared-peptide support is retained and is never presented as unique identity.
Quantitative comparisons align peptidoform/charge features, not group labels.
"""

from __future__ import annotations

import csv
import json
import math
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from .exact import read_fasta, sha256, write_table

MOD = re.compile(r"\[[^\[\]]*\]")


def peptide(sequence):
    """Validate a Sage peptide and remove modification syntax before I/L collapse."""
    value = MOD.sub("", sequence).replace("-", "").upper()  # "-" joins an N-terminal mod
    if not value or re.fullmatch("[A-Z]+", value) is None:
        raise ValueError(f"Invalid Sage peptide: {sequence!r}")
    return value.replace("I", "L")


def finite(value, label):
    """Convert a value to float and reject NaN or infinity with a field label."""
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Nonfinite {label}")
    return number


def search_inputs(folder):
    """Resolve a single-acquisition search and separate settings from file paths.

    Return the search folder, existing FASTA, existing mzML, non-file settings
    and configuration path. Relative input paths resolve beside config.json.
    """
    folder = Path(folder).resolve()
    config_path = folder / "config.json"
    config = json.loads(config_path.read_text())
    fasta = Path(config["database"]["fasta"])
    if not fasta.is_absolute():
        fasta = folder / fasta
    mzml = config.get("mzml_paths", [])
    if not isinstance(mzml, list) or len(mzml) != 1 or not isinstance(mzml[0], str):
        raise ValueError("Compare one acquisition per search directory")
    spectrum = Path(mzml[0])
    if not spectrum.is_absolute():
        spectrum = folder / spectrum
    settings = json.loads(json.dumps(config))
    del settings["database"]["fasta"]
    settings.pop("output_directory", None)
    settings.pop("mzml_paths", None)
    return folder, fasta.resolve(strict=True), spectrum.resolve(strict=True), settings, config_path


def load_search(folder, fasta, spectrum, q, protein_q, decoy_tag="rev_"):
    """Read accepted target evidence and positive LFQ features for one acquisition.

    Apply peptide q and optional protein q thresholds to PSMs. Verify reported
    members contain each canonical peptide. Keep exact peptidoform/charge LFQ
    keys and record missing/nonpositive exclusions; do not apply an LFQ-q gate.
    Return sequence aliases, compatible support, accepted peptides and LFQ maps.
    """
    records = read_fasta(fasta)
    hash_records = defaultdict(list)
    for accession, record in records.items():
        if accession.startswith(decoy_tag):
            continue
        hash_records[record[2]].append(accession)
    support = defaultdict(set)
    unique_local = defaultdict(set)
    accepted = set()
    psms = 0
    path = folder / "results.sage.tsv"
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"peptide", "proteins", "label", "peptide_q", "filename"}
        if protein_q is not None:
            required.add("protein_q")
        if len(reader.fieldnames or []) != len(set(reader.fieldnames or [])):
            raise ValueError("Duplicate Sage columns")
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Missing required Sage columns")
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError("Malformed Sage row")
            if Path(row["filename"].replace("\\", "/")).name != spectrum.name:
                raise ValueError("Sage result belongs to a different acquisition")
            if int(row["label"]) != 1:
                continue
            value = finite(row["peptide_q"], "peptide q")
            if not 0 <= value <= 1:
                raise ValueError("Peptide q outside [0, 1]")
            if value > q:
                continue
            if protein_q is not None:
                pq = finite(row["protein_q"], "protein q")
                if not 0 <= pq <= 1:
                    raise ValueError("Protein q outside [0, 1]")
                if pq > protein_q:
                    continue
            seq = peptide(row["peptide"])
            accessions = row["proteins"].split(";")
            if not accessions or any(
                a not in records or a.startswith(decoy_tag) for a in accessions
            ):
                raise ValueError("Accepted target has absent/decoy protein accession")
            hashes = set()
            for accession in accessions:
                _, protein, digest = records[accession]
                if seq not in protein.replace("I", "L"):
                    raise ValueError(
                        f"Accepted peptide not contained in reported protein: {accession}"
                    )
                support[digest].add(seq)
                hashes.add(digest)
            if len(hashes) == 1:
                unique_local[next(iter(hashes))].add(seq)
            accepted.add(seq)
            psms += 1
    lfq = {}
    excluded_lfq = []
    lfq_families = defaultdict(list)
    lfq_path = folder / "lfq.tsv"
    with lfq_path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError("Duplicate LFQ columns")
        fixed = {"peptide", "charge", "proteins", "q_value", "score", "spectral_angle"}
        columns = [f for f in fields if f not in fixed]
        if (
            not {"peptide", "charge", "proteins"} <= set(fields)
            or len(columns) != 1
            or Path(columns[0].replace("\\", "/")).name != spectrum.name
        ):
            raise ValueError("Ambiguous/missing LFQ intensity column")
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError("Malformed LFQ row")
            seq = peptide(row["peptide"])
            if seq not in accepted:
                continue
            text = row[columns[0]]
            if text.strip().lower() in {"", "na", "nan"}:
                excluded_lfq.append(
                    dict(
                        peptidoform=row["peptide"],
                        charge=row["charge"],
                        raw_intensity=text,
                        reason="missing",
                    )
                )
                continue
            value = finite(text, "LFQ intensity")
            if value <= 0:
                excluded_lfq.append(
                    dict(
                        peptidoform=row["peptide"],
                        charge=row["charge"],
                        raw_intensity=text,
                        reason="nonpositive",
                    )
                )
                continue
            # Sage modifications are numeric brackets, kept in this feature key.
            key = (row["peptide"].upper(), int(row["charge"]))
            # Sage 0.14.6 serializes charge-combined LFQ features with -1.
            if key[1] == 0 or key[1] < -1:
                raise ValueError("LFQ charge must be -1 (combined) or positive")
            if key in lfq:
                raise ValueError(
                    "Duplicate accepted LFQ peptidoform/charge; resolve before comparison"
                )
            for accession in row["proteins"].split(";"):
                if accession not in records or seq not in records[accession][1].replace("I", "L"):
                    raise ValueError("LFQ member is absent or lacks its peptide")
            lfq[key] = value
            lfq_families[(key[0].replace("I", "L"), key[1])].append(key)
    charge_modes = {"combined" if k[1] == -1 else "separate" for k in lfq}
    if len(charge_modes) > 1:
        raise ValueError("LFQ mixes combined and separate charge modes")
    return dict(
        charge_modes=sorted(charge_modes),
        excluded_lfq=excluded_lfq,
        records=records,
        aliases=hash_records,
        support=support,
        unique_local=unique_local,
        accepted=accepted,
        psms=psms,
        lfq=lfq,
        lfq_families=lfq_families,
    )


def rank(values):
    """Return one-based ranks in input order, averaging the ranks of tied values."""
    order = sorted(range(len(values)), key=values.__getitem__)
    ranked = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        for index in order[i:j]:
            ranked[index] = (i + j - 1) / 2 + 1
        i = j
    return ranked


def correlation(x, y):
    """Return Pearson correlation, or None for fewer than two or constant values.

    The caller supplies equal-length, finite paired observations.
    """
    if len(x) < 2:
        return None
    mx, my = statistics.mean(x), statistics.mean(y)
    denominator = math.sqrt(sum((v - mx) ** 2 for v in x) * sum((v - my) ** 2 for v in y))
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / denominator if denominator else None


def compare(
    left,
    right,
    output,
    *,
    q=0.01,
    protein_q=None,
    labels=("metaG", "de_novo"),
    compress_tables=False,
):
    """Compare two searches of the same acquisition with identical non-file settings.

    Parameters
    ----------
    left, right : path-like
        Completed Sage search directories containing config.json, results.sage.tsv
        and lfq.tsv. The recorded mzML inputs must resolve identically or hash equally.
    output : path-like
        New directory for overlap tables and COMPLETE.json; existing paths are rejected.
    q : float
        Inclusive peptide-q threshold in [0, 1].
    protein_q : float or None
        Optional additional protein-q threshold; None leaves it unused.
    labels : tuple of str
        Two distinct simple arm names used in the exported columns.
    compress_tables : bool
        Write deterministic gzip TSVs when true.

    Returns
    -------
    dict
        Counts, shared-feature correlations, exclusions, settings and file hashes.
        Protein support is sequence compatibility, including shared peptides;
        it does not identify a unique biological protein.
    """
    if not math.isfinite(q) or not 0 <= q <= 1:
        raise ValueError("Peptide q threshold outside [0, 1]")
    if protein_q is not None and (not math.isfinite(protein_q) or not 0 <= protein_q <= 1):
        raise ValueError("Protein q threshold outside [0, 1]")
    if (
        len(labels) != 2
        or labels[0] == labels[1]
        or any(not re.fullmatch("[A-Za-z0-9_-]+", x) for x in labels)
    ):
        raise ValueError("Two distinct simple arm labels are required")
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    a, b = search_inputs(left), search_inputs(right)
    if a[3] != b[3]:
        raise ValueError("Non-file Sage settings differ between arms")
    if a[2] != b[2] and sha256(a[2]) != sha256(b[2]):
        raise ValueError("Spectra differ between arms")
    decoy_tag = a[3]["database"].get("decoy_tag", "rev_")
    if not isinstance(decoy_tag, str) or not decoy_tag:
        raise ValueError("Sage decoy tag must be a nonempty string")
    arms = [load_search(x[0], x[1], x[2], q, protein_q, decoy_tag) for x in (a, b)]
    left_data, right_data = arms
    universe = set(left_data["aliases"]) | set(right_data["aliases"])
    if not universe:
        raise ValueError("Compared FASTAs contain no target protein sequences")
    cells = {"both": 0, labels[0] + "_only": 0, labels[1] + "_only": 0, "neither": 0}
    rows = []
    for digest in sorted(universe):
        left_value = left_data["support"].get(digest, set())
        r = right_data["support"].get(digest, set())
        category = (
            "both"
            if left_value and r
            else labels[0] + "_only"
            if left_value
            else labels[1] + "_only"
            if r
            else "neither"
        )
        cells[category] += 1
        row = dict(protein_sha256=digest, category=category)
        available = next(arm for arm in arms if digest in arm["aliases"])
        sequence = available["records"][available["aliases"][digest][0]][1]
        row.update(protein_sequence=sequence, protein_length=len(sequence))
        for name, arm, other in zip(labels, arms, reversed(arms)):
            aliases = sorted(arm["aliases"].get(digest, []))
            peptides = arm["support"].get(digest, set())
            row.update(
                {
                    name + "_in_database": int(bool(aliases)),
                    name + "_accessions": ";".join(aliases),
                    name + "_headers": json.dumps([arm["records"][x][0] for x in aliases]),
                    name + "_accepted_peptides": ";".join(sorted(peptides)),
                    name + "_peptide_count": len(peptides),
                    name + "_single_sequence_candidate_peptides": len(
                        arm["unique_local"].get(digest, set())
                    ),
                    name + "_support_peptides_accepted_by_other_arm": len(
                        peptides & other["accepted"]
                    ),
                    name + "_support_peptides_not_accepted_by_other_arm": ";".join(
                        sorted(peptides - other["accepted"])
                    ),
                }
            )
        rows.append(row)
    peptide_rows = []
    for seq in sorted(left_data["accepted"] | right_data["accepted"]):
        peptide_rows.append(
            dict(
                peptide=seq,
                **{
                    labels[0]: int(seq in left_data["accepted"]),
                    labels[1]: int(seq in right_data["accepted"]),
                },
            )
        )
    features = []
    for key in sorted(set(left_data["lfq"]) | set(right_data["lfq"])):
        left_value, r = left_data["lfq"].get(key), right_data["lfq"].get(key)
        family = (key[0].replace("I", "L"), key[1])
        ambiguous = any(len(arm["lfq_families"].get(family, [])) > 1 for arm in arms)
        features.append(
            dict(
                peptidoform=key[0],
                peptidoform_il=family[0],
                charge=key[1],
                il_family_multiple_rows=int(ambiguous),
                **{labels[0] + "_lfq": left_value, labels[1] + "_lfq": r},
                log2_right_over_left=math.log2(r) - math.log2(left_value)
                if left_value is not None and r is not None and not ambiguous
                else None,
            )
        )
    all_shared = set(left_data["lfq"]) & set(right_data["lfq"])
    shared = sorted(
        k
        for k in all_shared
        if all(
            len(arm["lfq_families"].get((k[0].replace("I", "L"), k[1]), [])) == 1 for arm in arms
        )
    )
    x = [math.log2(left_data["lfq"][k]) for k in shared]
    y = [math.log2(right_data["lfq"][k]) for k in shared]
    summary = dict(
        comparison_code_sha256=sha256(__file__),
        spectra_identity_check="same_resolved_path" if a[2] == b[2] else "equal_sha256",
        labels=list(labels),
        peptide_q=q,
        protein_q=protein_q,
        protein_scope=(
            "exact-sequence-compatible accepted PSM support; shared peptides remain ambiguous"
        ),
        universe="union of target protein sequences actually searched in these two databases",
        protein_cells=cells,
        protein_universe=len(universe),
        peptide_counts={name: len(arm["accepted"]) for name, arm in zip(labels, arms)},
        peptide_shared=len(left_data["accepted"] & right_data["accepted"]),
        peptide_left_only=len(left_data["accepted"] - right_data["accepted"]),
        peptide_right_only=len(right_data["accepted"] - left_data["accepted"]),
        accepted_psms={name: arm["psms"] for name, arm in zip(labels, arms)},
        lfq_features={name: len(arm["lfq"]) for name, arm in zip(labels, arms)},
        lfq_charge_modes={name: arm["charge_modes"] for name, arm in zip(labels, arms)},
        lfq_excluded_counts={
            name: dict(Counter(r["reason"] for r in arm["excluded_lfq"]))
            for name, arm in zip(labels, arms)
        },
        lfq_shared_features=len(shared),
        lfq_shared_exact_rows=len(all_shared),
        lfq_shared_rows_excluded_for_il_ambiguity=len(all_shared) - len(shared),
        lfq_canonical_family_counts={
            name: len(arm["lfq_families"]) for name, arm in zip(labels, arms)
        },
        lfq_canonical_family_shared=len(
            set(left_data["lfq_families"]) & set(right_data["lfq_families"])
        ),
        lfq_shared_log2_pearson=correlation(x, y),
        lfq_shared_spearman=correlation(rank(x), rank(y)),
        lfq_median_log2_right_over_left=statistics.median([b - a for a, b in zip(x, y)])
        if x
        else None,
        lfq_rule=(
            "positive finite exact reported peptidoform/charge rows whitelisted by "
            "accepted canonical peptides; no LFQ-q gate; agreement excludes any "
            "I/L family with multiple rows in either arm, without summing or "
            "selecting a representative"
        ),
        hashes={
            name: {
                "fasta": sha256(paths[1]),
                "config": sha256(paths[4]),
                "psms": sha256(paths[0] / "results.sage.tsv"),
                "lfq": sha256(paths[0] / "lfq.tsv"),
            }
            for name, paths in zip(labels, (a, b))
        },
    )
    if sum(cells.values()) != len(universe):
        raise RuntimeError("Protein overlap does not exhaust the declared universe")
    output.mkdir(parents=True, exist_ok=False)
    suffix = ".tsv.gz" if compress_tables else ".tsv"
    table_names = {key: key + "_overlap" + suffix for key in ("protein", "peptide", "lfq")}
    table_names["lfq_excluded"] = "lfq_excluded" + suffix
    summary["tables"] = table_names
    write_table(output / table_names["protein"], list(rows[0]), rows)
    write_table(output / table_names["peptide"], ["peptide", *labels], peptide_rows)
    write_table(
        output / table_names["lfq"],
        [
            "peptidoform",
            "peptidoform_il",
            "charge",
            "il_family_multiple_rows",
            labels[0] + "_lfq",
            labels[1] + "_lfq",
            "log2_right_over_left",
        ],
        features,
    )
    write_table(
        output / table_names["lfq_excluded"],
        ["arm", "peptidoform", "charge", "raw_intensity", "reason"],
        [dict(arm=name, **r) for name, arm in zip(labels, arms) for r in arm["excluded_lfq"]],
    )
    summary["output_sha256"] = {n: sha256(output / n) for n in table_names.values()}
    (output / "COMPLETE.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    return summary
