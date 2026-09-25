"""Strict sequence and TPM inputs for optional post-selection molecular unions.

Protein identity follows FastaLake hash version 5.3; peptide I/L equivalence is
deliberately absent from protein hashing. No cluster representative substitution.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

from fasta_lake.gzip_io import open_gzip_text
from fasta_lake.hashing import normalize_sequence


def sha256(path):
    """Hash file bytes in 8 MiB blocks and return the SHA-256 hex digest."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path, *, max_records=None):
    """Return accession -> (full header, normalized sequence, SHA-256)."""
    if max_records is not None and (type(max_records) is not int or max_records < 1):
        raise ValueError("max_records must be a positive integer")
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    result, header, parts = {}, None, []

    def save():
        """Store a normalized record with its exact-sequence hash, enforcing ID limits."""
        if header is None:
            return
        accession = header.split()[0] if header.split() else ""
        sequence = normalize_sequence("".join(parts)).decode("ascii")
        if not accession or not sequence or accession in result:
            raise ValueError(f"Empty/duplicate FASTA accession or empty sequence: {accession!r}")
        result[accession] = (header, sequence, hashlib.sha256(sequence.encode("ascii")).hexdigest())
        if max_records is not None and len(result) > max_records:
            raise ValueError(
                f"FASTA exceeds the {max_records} record limit; select a smaller subset"
            )

    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith(">"):
                save()
                header, parts = line[1:].rstrip("\r\n"), []
            elif line.strip():
                if header is None:
                    raise ValueError("Sequence before first FASTA header")
                parts.append(line.strip())
        save()
    if not result:
        raise ValueError("Empty FASTA")
    return result


def numeric(value, *, nan_is_missing=False):
    """Parse a nonnegative finite TPM while preserving explicit missing tokens.

    Blank, NA, N/A, null and none return None. NaN is missing only when requested.
    A missing TSV field, negative value or nonfinite number raises ValueError.
    """
    if value is None:
        raise ValueError("Malformed TSV row: missing field")
    if value.strip().lower() in {"", "na", "n/a", "null", "none"} or (
        nan_is_missing and value.strip().lower() == "nan"
    ):
        return None
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise ValueError(f"Invalid nonnegative finite TPM: {value!r}")
    return number


def read_tpm(path, *, nan_is_missing=False):
    """Read plain or gzipped gene TPM rows without converting missing values to zero.

    Require unique prodigal_protein IDs and metaG_tpm/metaT_tpm columns.
    Return a gene-keyed dictionary with parsed dna/rna and original raw strings;
    malformed rows, duplicate genes and invalid numeric values raise ValueError.
    """
    required = {"prodigal_protein", "metaG_tpm", "metaT_tpm"}
    rows = {}
    opener = gzip.open if Path(path).suffix == ".gz" else open
    with opener(path, "rt", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or not required <= set(fields):
            raise ValueError("TPM requires unique prodigal_protein, metaG_tpm, metaT_tpm columns")
        for row in reader:
            gene = row["prodigal_protein"]
            if None in row or not gene or gene != gene.strip() or gene in rows:
                raise ValueError(f"Malformed or duplicate TPM gene: {gene!r}")
            rows[gene] = dict(
                dna=numeric(row["metaG_tpm"], nan_is_missing=nan_is_missing),
                rna=numeric(row["metaT_tpm"], nan_is_missing=nan_is_missing),
                raw_dna=row["metaG_tpm"],
                raw_rna=row["metaT_tpm"],
            )
    if not rows:
        raise ValueError("TPM table has no genes")
    return rows


def top_cutoff(values, percent, rounding="floor"):
    """Top percentage among positive GENES; zero disables, ties all included."""
    if not math.isfinite(percent) or not 0 <= percent <= 100:
        raise ValueError("DNA top percent must be finite and in [0, 100]")
    if rounding not in {"floor", "ceil"}:
        raise ValueError("Rank rounding must be floor or ceil")
    values = list(values)
    if any(v is not None and (not math.isfinite(v) or v < 0) for v in values):
        raise ValueError("Invalid TPM values")
    positive = sorted(v for v in values if v is not None and v > 0)
    if not positive or percent == 0:
        return None
    rank = max(1, (math.floor if rounding == "floor" else math.ceil)(len(positive) * percent / 100))
    return positive[-rank]


def write_table(path, fields, rows):
    """Write dictionary rows as a TSV, using deterministic gzip for a .gz path.

    The supplied field order defines the columns. Existing destination files
    are overwritten; callers enforce new output directories where required.
    """
    path = Path(path)
    opener = open_gzip_text if path.suffix == ".gz" else (lambda p: open(p, "w", newline=""))
    with opener(path) as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_source(manifest_path):
    """Load one explicitly specimen-qualified quantification reference.

    Checksums and all-gene agreement detect changed/mismatched inputs. The
    provenance record must establish how the original quantification reference
    was linked to these proteins; the caller's specimen assertion is not inferred
    from filenames. A hash alone cannot establish biological specimen identity.
    """
    manifest_path = Path(manifest_path).resolve()
    obj = json.loads(manifest_path.read_text())
    if obj.get("schema") != "fastalake.molecular-source.v1":
        raise ValueError("Unsupported molecular source schema")
    for field in ("specimen", "reference_id"):
        if not isinstance(obj.get(field), str) or not obj[field].strip():
            raise ValueError(f"Missing {field}")
    inputs = {}
    for field in ("tpm", "proteins", "provenance"):
        record = obj[field]
        path = Path(record["path"])
        if not path.is_absolute():
            path = manifest_path.parent / path
        path = path.resolve(strict=True)
        if sha256(path) != record["sha256"]:
            raise ValueError(f"{field} checksum mismatch")
        inputs[field] = path
    if obj.get("reference_status") != "verified_original_quantification_reference":
        raise ValueError("Original quantification reference is not verified")
    genes = read_tpm(inputs["tpm"])
    proteins = read_fasta(inputs["proteins"])
    if set(genes) != set(proteins):
        missing = sorted(set(genes) - set(proteins))
        extra = sorted(set(proteins) - set(genes))
        raise ValueError(
            f"All-gene reference mismatch: missing={len(missing)}, extra={len(extra)}; "
            f"examples={missing[:5]}"
        )
    return obj, inputs, genes, proteins


def build_union(
    baseline,
    source_manifest,
    specimen,
    output,
    *,
    dna_percent=1.0,
    dna_absolute=None,
    rna_absolute=None,
    rounding="floor",
    rna_percent=0,
    background=None,
    metabolites=False,
):
    """Keep every baseline sequence and append independently selected molecular proteins."""
    from .selection import build_selection

    return build_selection(
        source_manifest,
        specimen,
        output,
        baseline=baseline,
        background=background,
        dna_percent=dna_percent,
        rna_percent=rna_percent,
        dna_absolute=dna_absolute,
        rna_absolute=rna_absolute,
        rounding=rounding,
        metabolites=metabolites,
    )
