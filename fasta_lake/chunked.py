"""Extract a prepared reference lake in bounded pieces, then merge in source order.

Every piece sees the same full prediction roster. This changes storage/execution,
not peptide ranking, sequence identity or downstream inference.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import shutil
import time
from pathlib import Path

from fasta_lake.resources import run_recorded
from fasta_lake.rust.binaries import find_binary


def _sha(path):
    """Stream a file into a lowercase SHA-256 hex digest."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _stamp(path):
    """Record a file's byte count and SHA-256 without loading it into memory."""
    path = Path(path)
    return {"bytes": path.stat().st_size, "sha256": _sha(path)}


def _predictions(path):
    """Mirror Rust's recursive, case-sensitive CSV roster, including symlink files."""
    if path.is_file():
        if path.suffix != ".csv":
            raise ValueError("Chunked extraction requires a predictions CSV or CSV directory")
        files = [path]
    else:
        files = []
        for folder, _, names in os.walk(path, followlinks=True):
            files.extend(Path(folder) / n for n in names if n.endswith(".csv"))
        files.sort()
    if not files:
        raise ValueError("No prediction CSV files found")
    return [{"path": str(p), **_stamp(p)} for p in files]


def _pieces(source, scratch, byte_limit, record_limit, digest):
    """Yield complete-record pieces without buffering a whole FASTA record.

    A piece may exceed the byte target by its last complete record; record_limit
    also bounds Rust's per-record offset table. Only one scratch piece is live.
    """
    current = None
    count = size = number = total = 0
    sequence_bytes = 0
    try:
        with source.open("rb") as stream:
            for line in stream:
                digest.update(line)
                if line.startswith(b">"):
                    if current is not None and sequence_bytes == 0:
                        raise ValueError(f"Empty FASTA record before record {total + 1}")
                    if len(line[1:].strip()) == 0:
                        raise ValueError("Empty FASTA header")
                    if current is not None and (size >= byte_limit or count >= record_limit):
                        current.close()
                        current = None
                        yield number, count, size, total - count
                        count = size = 0
                    if current is None:
                        number += 1
                        current = scratch.open("xb")
                    count += 1
                    total += 1
                    sequence_bytes = 0
                else:
                    if current is None:
                        if line.strip():
                            raise ValueError("Expected plain FASTA beginning with a header")
                        continue
                    if b">" in line:
                        raise ValueError("A FASTA header marker must begin a new line")
                    sequence_bytes += len(line.rstrip(b"\r\n"))
                current.write(line)
                size += len(line)
        if current is None or sequence_bytes == 0:
            raise ValueError("Empty FASTA input or final record")
        current.close()
        current = None
        yield number, count, size, total - count
    finally:
        if current is not None:
            current.close()


def extract_chunked(
    database,
    predictions,
    output,
    *,
    chunk_bytes=256 * 1024 * 1024,
    chunk_records=500_000,
    top_fraction=0.30,
    min_score=0.0,
    min_length=9,
    max_length=50,
    normalize_il=True,
    threads=1,
    output_hits=False,
    keep_pieces=False,
    binary_path=None,
):
    """Extract and merge without altering the prepared input lake's identities.

    Output is a new directory containing evidence.fasta, evidence.tsv, optional
    hits.tsv, per-piece receipts/logs and a checksummed COMPLETE.json. No razor,
    clustering, additional deduplication or independently ranked peptide batches
    are run here. The existing extractor's complete per-sample prediction ranking
    is repeated identically for every reference piece.
    """
    if any(
        not isinstance(n, int) or isinstance(n, bool) or n < 1
        for n in (chunk_bytes, chunk_records, threads, min_length, max_length)
    ):
        raise ValueError("Chunk sizes, thread count and peptide lengths must be positive integers")
    if min_length > max_length or not math.isfinite(min_score):
        raise ValueError("Invalid peptide length range or score")
    if top_fraction is not None and (not math.isfinite(top_fraction) or not 0 < top_fraction <= 1):
        raise ValueError("top_fraction must be in (0, 1], or None for a flat score threshold")
    database = Path(database).resolve(strict=True)
    predictions = Path(predictions).resolve(strict=True)
    out = Path(output).absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"Output exists: {out}; choose a new directory")
    if not database.is_file() or database.stat().st_size == 0:
        raise ValueError("Reference FASTA must be a nonempty file")
    if predictions.is_dir() and predictions in out.resolve().parents:
        raise ValueError("Output must be outside the prediction directory")
    binary = find_binary("fasta_extractor", explicit_path=binary_path)
    inputs = {
        "database": {"path": str(database), **_stamp(database)},
        "predictions": _predictions(predictions),
        "binary": {"path": str(binary), **_stamp(binary)},
    }
    out.mkdir(parents=True, exist_ok=False)
    scratch = out / "reference_piece.fasta"
    started = time.monotonic()
    parameters = {
        "chunk_bytes": chunk_bytes,
        "chunk_records": chunk_records,
        "top_fraction": top_fraction,
        "min_score": min_score,
        "min_length": min_length,
        "max_length": max_length,
        "normalize_il": normalize_il,
        "threads": threads,
        "output_hits": output_hits,
        "keep_pieces": keep_pieces,
    }
    (out / "INPUTS.json").write_text(
        json.dumps({"inputs": inputs, "parameters": parameters}, indent=2) + "\n"
    )
    records = []
    read_digest = hashlib.sha256()
    conserved_keys = ("total_peptides", "score_filter", "normalize_il", "peptide_length_range")
    conserved = None
    total_hits = 0
    evidence_header = hits_header = None
    with (
        (out / "evidence.fasta").open("xb") as merged,
        (out / "evidence.tsv").open("xb") as evidence,
    ):
        # Keep the optional hit stream separate so it cannot allocate a full table.
        hits = (out / "hits.tsv").open("xb") if output_hits else None
        try:
            for number, count, size, offset in _pieces(
                database, scratch, chunk_bytes, chunk_records, read_digest
            ):
                piece = out / "pieces" / f"{number:06d}"
                piece.mkdir(parents=True)
                command = [
                    str(binary),
                    "--database",
                    str(scratch),
                    "--peptides-dir",
                    str(predictions),
                    "--output",
                    str(piece / "matched.fasta"),
                    "--output-evidence",
                    str(piece / "evidence.tsv"),
                    "--min-length",
                    str(min_length),
                    "--max-length",
                    str(max_length),
                    "--normalize-il",
                    str(normalize_il).lower(),
                    "--threads",
                    str(threads),
                ]
                if top_fraction is not None:
                    command += ["--top-percent", str(top_fraction)]
                else:
                    command += ["--min-score", str(min_score)]
                if output_hits:
                    command += ["--output-hits", str(piece / "hits.tsv")]
                reference = _stamp(scratch)
                if reference["bytes"] != size:
                    raise RuntimeError("Scratch reference changed before extraction")
                begin = time.monotonic()
                with (
                    (piece / "stdout.log").open("xb") as stdout,
                    (piece / "stderr.log").open("xb") as stderr,
                ):
                    resources = run_recorded(command, stdout=stdout, stderr=stderr)
                record = {
                    "piece": number,
                    "first_record_zero_based": offset,
                    "reference_records": count,
                    "reference": reference,
                    "command": command,
                    "exit_code": resources["exit_code"],
                    "resources": resources,
                    "seconds": time.monotonic() - begin,
                }
                (piece / "RUN.json").write_text(json.dumps(record, indent=2) + "\n")
                if resources["exit_code"]:
                    raise RuntimeError(f"Piece {number} failed; inspect {piece / 'stderr.log'}")
                if _predictions(predictions) != inputs["predictions"]:
                    raise RuntimeError("Prediction inputs changed during chunked extraction")
                stats = json.loads((piece / "matched.stats.json").read_text())
                if stats["proteins_searched"] != count:
                    raise RuntimeError("Rust and streaming FASTA record counts disagree")
                current = {key: stats[key] for key in conserved_keys}
                if conserved is not None and current != conserved:
                    raise RuntimeError("Peptide selection changed between reference pieces")
                conserved = current
                with (piece / "matched.fasta").open("rb") as source:
                    shutil.copyfileobj(source, merged)
                with (piece / "evidence.tsv").open("rb") as source:
                    header = source.readline()
                    if evidence_header is None:
                        evidence_header = header
                        evidence.write(header)
                    if not header or header != evidence_header:
                        raise RuntimeError("Piece evidence headers differ")
                    shutil.copyfileobj(source, evidence)
                with (piece / "evidence.tsv").open() as source:
                    rows = 0
                    for row in csv.DictReader(source, delimiter="\t"):
                        total_hits += int(row["n_peptides"])
                        rows += 1
                    if rows != stats["proteins_matched"]:
                        raise RuntimeError("Piece evidence count differs from matched proteins")
                if hits is not None:
                    with (piece / "hits.tsv").open("rb") as source:
                        header = source.readline()
                        if hits_header is None:
                            hits_header = header
                            hits.write(header)
                        if not header or header != hits_header:
                            raise RuntimeError("Piece hit headers differ")
                        shutil.copyfileobj(source, hits)
                record.update(
                    stats=stats,
                    outputs={p.name: _stamp(p) for p in sorted(piece.iterdir()) if p.is_file()},
                )
                record["merged_outputs"] = {
                    "matched.fasta": "evidence.fasta",
                    "evidence.tsv": "evidence.tsv",
                }
                if output_hits:
                    record["merged_outputs"]["hits.tsv"] = "hits.tsv"
                if not keep_pieces:
                    # Flush the merged products before reclaiming our generated copies.
                    # Failed runs lack COMPLETE.json and retain the original inputs.
                    merged.flush()
                    evidence.flush()
                    if hits is not None:
                        hits.flush()
                    for name in record["merged_outputs"]:
                        (piece / name).unlink()
                record["retained_outputs"] = sorted(p.name for p in piece.iterdir())
                records.append(record)
                (out / "PROGRESS.json").write_text(json.dumps(records, indent=2) + "\n")
                # Only this invocation's successfully consumed scratch copy is removed.
                scratch.unlink()
        finally:
            if hits is not None:
                hits.close()
    if read_digest.hexdigest() != inputs["database"]["sha256"]:
        raise RuntimeError("Reference FASTA changed during extraction")
    if _stamp(database) != {k: v for k, v in inputs["database"].items() if k != "path"}:
        raise RuntimeError("Reference FASTA changed during extraction")
    if _predictions(predictions) != inputs["predictions"]:
        raise RuntimeError("Prediction inputs changed during extraction")
    if _stamp(binary) != {k: v for k, v in inputs["binary"].items() if k != "path"}:
        raise RuntimeError("Extractor binary changed during extraction")
    searched = sum(r["stats"]["proteins_searched"] for r in records)
    matched = sum(r["stats"]["proteins_matched"] for r in records)
    stats = {
        **conserved,
        "database": str(database),
        "source_tag": "DB",
        "proteins_searched": searched,
        "proteins_matched": matched,
        "core_proteome_added": 0,
        "total_output_proteins": matched,
        "match_rate": matched / searched,
        "avg_peptides_per_protein": (total_hits / matched if matched else 0.0),
        "total_seconds": time.monotonic() - started,
        "extraction_mode": "reference_pieces",
        "pieces": len(records),
    }
    (out / "evidence.stats.json").write_text(json.dumps(stats, indent=2) + "\n")
    complete = {
        "status": "PASS",
        "schema": "fastalake.chunked-extraction.v1",
        "inputs": inputs,
        "parameters": parameters,
        "stats": stats,
        "implementation": _stamp(Path(__file__)),
        "pieces": records,
        "outputs": {p.name: _stamp(p) for p in sorted(out.iterdir()) if p.is_file()},
        "scope": "Exact extraction from a prepared lake; no inference or re-deduplication",
    }
    (out / "COMPLETE.json").write_text(json.dumps(complete, indent=2) + "\n")
    return complete
