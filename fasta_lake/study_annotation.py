"""Name each study group with the gene, protein and organism of its source headers.

Lake accessions are often sequence hashes (``TAG_<sha256>``), which make a matrix
unreadable. The original headers survive in the lake builder's ``--output-headers``
sidecar and, for lakes that keep catalogue headers, in the searched FASTA files.
This module reads either source and writes a readable annotation table plus an
annotated copy of the study-group matrix. It does not change any grouping: when
one sequence came from several catalogue entries, every distinct value is listed.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from fasta_lake.headers import parse_header

ANNOTATION_COLUMNS = ["gene", "protein_name", "organism", "source", "n_source_headers"]


def load_header_sidecar(path, wanted=None) -> dict[str, list[str]]:
    """Map each emitted lake accession to its original source headers.

    Parameters
    ----------
    path : path-like
        ``sha256_hash, n_sources, kept_header, all_source_headers`` TSV written by
        ``lake_builder --output-headers``. ``all_source_headers`` holds
        ``"TAG header"`` entries joined by ``|||``.
    wanted : set[str] or None
        Accessions to keep; None keeps all.

    Returns
    -------
    dict[str, list[str]]
        Original headers (source tag removed), in sidecar order.
    """
    headers: dict[str, list[str]] = {}
    with Path(path).open(newline="") as stream:
        columns = stream.readline().rstrip("\n").split("\t")
        if columns[:4] != ["sha256_hash", "n_sources", "kept_header", "all_source_headers"]:
            raise ValueError(f"{path}: not a lake_builder --output-headers sidecar")
        for line in stream:
            fields = line.rstrip("\n").split("\t", 3)
            if len(fields) < 4 or not fields[2].split():
                continue
            accession = fields[2].split()[0]
            if wanted is not None and accession not in wanted:
                continue
            entries = []
            for entry in fields[3].split("|||"):
                # "TAG header"; older sidecars separated the tag with a tab.
                parts = entry.split(None, 1)
                if len(parts) == 2:
                    entries.append(parts[1].strip())
            headers[accession] = entries
    return headers


def fasta_headers(paths, wanted=None) -> dict[str, list[str]]:
    """Read header lines for the wanted accessions from FASTA files."""
    headers: dict[str, list[str]] = {}
    for path in paths:
        with Path(path).open() as stream:
            for line in stream:
                if not line.startswith(">"):
                    continue
                header = line[1:].strip()
                accession = header.split()[0] if header else ""
                if accession and (wanted is None or accession in wanted):
                    known = headers.setdefault(accession, [])
                    if header not in known:
                        known.append(header)
    return headers


def searched_fastas(search_root, config_name="config.json") -> list[Path]:
    """Return the distinct FASTA files named by each acquisition's Sage record."""
    paths = []
    for config in sorted(Path(search_root).glob(f"*/{config_name}")):
        fasta = Path(json.loads(config.read_text())["database"]["fasta"])
        if fasta not in paths:
            paths.append(fasta)
    return paths


def describe(accession: str, headers: list[str]) -> dict[str, str]:
    """Summarise parsed headers; distinct values are joined with ``; ``."""
    values = {"gene": [], "protein_name": [], "organism": [], "source": []}
    for header in headers:
        parsed = parse_header(header)
        # A gene equal to the accession (MGYG, GMGC, smORF) is an identifier, not a name.
        gene = parsed.gene if parsed.gene and parsed.gene != parsed.accession else None
        for key, value in (
            ("gene", gene),
            ("protein_name", parsed.description),
            ("organism", parsed.organism),
            ("source", parsed.source_db if parsed.source_db != "unknown" else None),
        ):
            if value and value not in values[key]:
                values[key].append(value)
    row = {key: "; ".join(found) for key, found in values.items()}
    row["n_source_headers"] = str(len(headers))
    return row


def annotate_study(groups_dir, *, headers_tsv=None, fastas=()) -> dict:
    """Write ``study_group_annotation.tsv`` and ``annotated_study_group_matrix.tsv``.

    Parameters
    ----------
    groups_dir : path-like
        Completed ``group_searches`` output.
    headers_tsv : path-like or None
        Lake header sidecar; preferred because it keeps every source header.
    fastas : iterable of path-like
        FASTA files whose header lines are read for groups the sidecar lacks.

    Returns
    -------
    dict
        Counts recorded in the study summary.
    """
    groups_dir = Path(groups_dir)
    with (groups_dir / "study_groups.tsv").open(newline="") as stream:
        groups = [row["study_group"] for row in csv.DictReader(stream, delimiter="\t")]
    wanted = set(groups)
    headers = load_header_sidecar(headers_tsv, wanted) if headers_tsv else {}
    missing = wanted - headers.keys()
    if missing and fastas:
        for accession, found in fasta_headers(fastas, missing).items():
            # A bare accession line carries no description; keep looking elsewhere.
            if any(header != accession for header in found):
                headers[accession] = found
    rows = {group: describe(group, headers.get(group, [])) for group in groups}
    columns = ["study_group"] + ANNOTATION_COLUMNS
    with (groups_dir / "study_group_annotation.tsv").open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for group in groups:
            writer.writerow([group] + [rows[group][key] for key in ANNOTATION_COLUMNS])
    matrix = groups_dir / "study_group_matrix.tsv"
    annotated = groups_dir / "annotated_study_group_matrix.tsv"
    with matrix.open(newline="") as source, annotated.open("x", newline="") as target:
        reader = csv.reader(source, delimiter="\t")
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        header = next(reader)
        writer.writerow(header[:1] + ["gene", "protein_name", "organism"] + header[1:])
        for row in reader:
            info = rows.get(row[0]) or describe(row[0], [])
            writer.writerow(
                row[:1] + [info["gene"], info["protein_name"], info["organism"]] + row[1:]
            )
    named = sum(1 for row in rows.values() if row["protein_name"] or row["gene"])
    return {
        "annotation_table": "study_group_annotation.tsv",
        "annotated_matrix": "annotated_study_group_matrix.tsv",
        "header_source": str(headers_tsv) if headers_tsv else "searched FASTA headers",
        "study_groups": len(groups),
        "named_study_groups": named,
        "scope": "Source-header labels of the representative; shared peptides can "
        "leave biological origin ambiguous",
    }
