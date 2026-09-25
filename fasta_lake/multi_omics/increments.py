"""Peptide and exact protein increments against one declared reference.

Reference-wide compatibility prevents a restricted FASTA from making a shared
peptide appear sequence-specific. Evidence counts are not protein inference,
protein-level FDR, species identity or independent LFQ validation.
"""

from __future__ import annotations

import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

from .complementarity import correlation, load_search, rank, search_inputs
from .exact import read_fasta, sha256, write_table

PROTEIN_FIELDS = (
    "protein_sha256",
    "sequence_length",
    "in_left_candidates",
    "in_right_candidates",
    "host_sequence",
    "left_peptides",
    "right_peptides",
    "left_reference_specific_peptides",
    "right_reference_specific_peptides",
    "gained_peptides",
    "lost_peptides",
    "left_peptide_sequences",
    "right_peptide_sequences",
    "gained_peptide_sequences",
    "lost_peptide_sequences",
    "reference_specific_peptide_sequences",
    "left_reported_peptides",
    "right_reported_peptides",
    "one_peptide_change",
    "two_peptides_change",
    "two_with_specific_change",
    "reference_accessions",
)
PEPTIDE_FIELDS = (
    "peptide",
    "change",
    "reference_proteins",
    "reference_specific",
    "reference_source",
    "protein_sha256",
)


def reference_map(records, peptides):
    """Map I/L-canonical peptides to distinct exact protein sequences."""
    import ahocorasick

    sequences = {r[2]: r[1] for r in records.values()}
    peptides = set(peptides)
    support = defaultdict(set)
    matches = {p: set() for p in peptides}
    if not peptides:
        return sequences, dict(support), matches
    automaton = ahocorasick.Automaton()
    for p in sorted(peptides):
        automaton.add_word(p, p)
    automaton.make_automaton()
    for digest, sequence in sequences.items():
        for _, p in automaton.iter(sequence.replace("I", "L")):
            support[digest].add(p)
            matches[p].add(digest)
    missing = [p for p, values in matches.items() if not values]
    if missing:
        raise ValueError(f"Accepted peptides absent from declared reference: {missing[:5]}")
    return sequences, dict(support), matches


def evidence_tables(left, right, sequences, support, matches, host=None):
    """Return exhaustive support rows and bidirectional counts.

    Compatible sequences outside both selected FASTAs are retained and labelled.
    Candidate membership and reported PSM membership remain separate.
    """
    lp, rp = left["accepted"], right["accepted"]
    lc, rc = set(left["aliases"]), set(right["aliases"])
    if not (lc | rc) <= set(sequences):
        raise ValueError("Search FASTA contains sequences outside declared reference")
    if not (lp | rp) <= set(matches):
        raise ValueError("Reference mapping omits accepted peptides")
    if host is not None and not set(host) <= set(sequences):
        raise ValueError("Host reference contains sequences outside declared reference")
    retained = lp & rp
    gains, losses = rp - lp, lp - rp
    unique = {p for p, hashes in matches.items() if len(hashes) == 1}
    summary = dict(
        peptides_left=len(lp),
        peptides_right=len(rp),
        peptides_retained=len(lp & rp),
        peptides_gained=len(gains),
        peptides_lost=len(losses),
        right_identity_recovered_by_left=len(lp & rp) / len(rp) if rp else None,
        left_identity_retained_by_right=len(lp & rp) / len(lp) if lp else None,
        left_candidates=len(lc),
        right_candidates=len(rc),
        candidate_sequences_added=len(rc - lc),
        candidate_sequences_removed=len(lc - rc),
        reference_sequences=len(sequences),
    )
    endpoints = {k: [set(), set()] for k in ("one_peptide", "two_peptides", "two_with_specific")}
    proteins = []
    for digest in sorted(support):
        a, b = support[digest] & lp, support[digest] & rp
        if not (a or b):
            continue
        au, bu = a & unique, b & unique
        row = dict(
            protein_sha256=digest,
            sequence_length=len(sequences[digest]),
            in_left_candidates=int(digest in lc),
            in_right_candidates=int(digest in rc),
            host_sequence="unclassified" if host is None else int(digest in host),
            left_peptides=len(a),
            right_peptides=len(b),
            left_reference_specific_peptides=len(au),
            right_reference_specific_peptides=len(bu),
            gained_peptides=len(b - a),
            lost_peptides=len(a - b),
            left_peptide_sequences=";".join(sorted(a)),
            right_peptide_sequences=";".join(sorted(b)),
            gained_peptide_sequences=";".join(sorted(b - a)),
            lost_peptide_sequences=";".join(sorted(a - b)),
            reference_specific_peptide_sequences=";".join(sorted(au | bu)),
            left_reported_peptides=len(left["support"].get(digest, set())),
            right_reported_peptides=len(right["support"].get(digest, set())),
        )
        for name, condition in (
            ("one_peptide", (bool(a), bool(b))),
            ("two_peptides", (len(a) >= 2, len(b) >= 2)),
            ("two_with_specific", (len(a) >= 2 and bool(au), len(b) >= 2 and bool(bu))),
        ):
            for i, passed in enumerate(condition):
                if passed:
                    endpoints[name][i].add(digest)
            row[name + "_change"] = (
                "retained"
                if all(condition)
                else "gained"
                if condition[1]
                else "lost"
                if condition[0]
                else "below_threshold"
            )
        proteins.append(row)
    for name, (a, b) in endpoints.items():
        for suffix, value in (
            ("left", len(a)),
            ("right", len(b)),
            ("retained", len(a & b)),
            ("gained", len(b - a)),
            ("lost", len(a - b)),
        ):
            summary[name + "_" + suffix] = value
        summary[name + "_gained_on_added_candidates"] = len((b - a) & (rc - lc))
        summary[name + "_gained_on_retained_candidates"] = len((b - a) & lc & rc)
        summary[name + "_gained_outside_right_candidates"] = len((b - a) - rc)
    peptides = []
    for p in sorted(lp | rp):
        hashes = matches[p]
        category = "retained" if p in retained else "gained" if p in gains else "lost"
        if host is None:
            source = "unclassified"
        else:
            source = (
                "host_only"
                if hashes <= host
                else "nonhost_only"
                if not hashes & host
                else "host_nonhost_ambiguous"
            )
        peptides.append(
            dict(
                peptide=p,
                change=category,
                reference_proteins=len(hashes),
                reference_specific=int(len(hashes) == 1),
                reference_source=source,
                protein_sha256=";".join(sorted(hashes)),
            )
        )
    if host is not None:
        for category in ("retained", "gained", "lost"):
            for source in ("host_only", "nonhost_only", "host_nonhost_ambiguous"):
                summary[f"peptides_{category}_{source}"] = sum(
                    r["change"] == category and r["reference_source"] == source for r in peptides
                )
    a, b = left["lfq"], right["lfq"]
    shared = sorted(set(a) & set(b))
    ratios = [math.log2(b[k] / a[k]) for k in shared]
    summary.update(
        lfq_left=len(a),
        lfq_right=len(b),
        lfq_shared=len(shared),
        lfq_gained=len(set(b) - set(a)),
        lfq_lost=len(set(a) - set(b)),
        lfq_ratio=len(b) / len(a) if a else None,
        lfq_spearman=correlation(rank([a[k] for k in shared]), rank([b[k] for k in shared])),
        lfq_median_log2_ratio=statistics.median(ratios) if ratios else None,
        lfq_median_absolute_log2_ratio=statistics.median(map(abs, ratios)) if ratios else None,
    )
    return summary, proteins, peptides


def increment_report(
    left,
    right,
    reference,
    output,
    *,
    host_reference=None,
    q=0.01,
    protein_q=None,
    labels=("de_novo", "with_molecular"),
):
    """Validate matched searches and export complete reference-wide evidence."""
    if not math.isfinite(q) or not 0 <= q <= 1:
        raise ValueError("Peptide q must be in [0, 1]")
    if protein_q is not None and (not math.isfinite(protein_q) or not 0 <= protein_q <= 1):
        raise ValueError("Protein q must be in [0, 1]")
    if len(labels) != 2 or labels[0] == labels[1] or any(not v.strip() for v in labels):
        raise ValueError("Two distinct nonempty labels are required")
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    a, b = search_inputs(left), search_inputs(right)
    if a[3] != b[3]:
        raise ValueError("Non-file Sage settings differ")
    if a[2] != b[2] and sha256(a[2]) != sha256(b[2]):
        raise ValueError("Spectra differ")
    decoy_tag = a[3]["database"].get("decoy_tag", "rev_")
    if not isinstance(decoy_tag, str) or not decoy_tag:
        raise ValueError("Invalid decoy tag")
    paths = {Path(reference)} | ({Path(host_reference)} if host_reference else set())
    for info in (a, b):
        paths.update([info[1], info[2], info[4], info[0] / "results.sage.tsv", info[0] / "lfq.tsv"])
    source_hashes = {str(p.resolve()): sha256(p) for p in paths}
    records = read_fasta(reference)
    if any(k.startswith(decoy_tag) for k in records):
        raise ValueError("Declared reference must contain targets only")
    host = None
    if host_reference is not None:
        host = {r[2] for r in read_fasta(host_reference).values()}
    arms = [load_search(x[0], x[1], x[2], q, protein_q, decoy_tag) for x in (a, b)]
    sequences, support, matches = reference_map(records, arms[0]["accepted"] | arms[1]["accepted"])
    summary, proteins, peptides = evidence_tables(*arms, sequences, support, matches, host)
    if any(sha256(p) != h for p, h in source_hashes.items()):
        raise ValueError("Source changed during increment report")
    aliases = defaultdict(list)
    for accession, record in records.items():
        aliases[record[2]].append(accession)
    for row in proteins:
        row["reference_accessions"] = ";".join(sorted(aliases[row["protein_sha256"]]))
    output.mkdir(parents=True)
    write_table(output / "protein_support.tsv.gz", PROTEIN_FIELDS, proteins)
    write_table(output / "peptide_changes.tsv.gz", PEPTIDE_FIELDS, peptides)
    summary.update(
        labels=list(labels),
        peptide_q=q,
        protein_q=protein_q,
        protein_scope="Compatibility with declared reference; no protein-level FDR inference",
        specificity_scope="One exact sequence in declared reference; not universal uniqueness",
        lfq_scope="Original same-acquisition positive accepted-peptide features; no joint LFQ/MBR",
    )
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    receipt = dict(
        status="PASS",
        input_sha256=source_hashes,
        code_sha256={
            p.name: sha256(p)
            for p in (
                Path(__file__),
                Path(__file__).with_name("complementarity.py"),
                Path(__file__).with_name("exact.py"),
                Path(__file__).parents[1] / "hashing.py",
            )
        },
        output_sha256={p.name: sha256(p) for p in output.iterdir() if p.is_file()},
    )
    (output / "COMPLETE.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return summary
