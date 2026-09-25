"""
Lake builder: merge multiple source FASTAs into one with provenance (v5.2).

Reads N input FASTAs (each tagged with a source label such as "SP", "TR",
"UH"), parses headers into ``SourceRecord`` objects, groups sequences by
SHA-256 hash (normalized sequence; see ``hashing.py``), and for each group
emits a single merged representative in the output FASTA plus one JSONL line
in the provenance sidecar.

Priority-win mode (back-compat) picks one source and discards the others.
Merge mode (default in v5.2) combines annotations per the rules in
``fasta_lake.lake_merge``.

Header parser is deliberately conservative: it knows UniProt sp|/tr|, UHGP
MGYP, GENCODE ENSP, and falls back to "UNKNOWN" for other formats.
"""

from __future__ import annotations

import logging
import re
import time
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

from fasta_lake.lake_merge import (
    DEFAULT_PRIORITY,
    MergedRecord,
    SourceRecord,
    format_fasta_header,
    merge_records,
    priority_win,
    sequence_hash,
    write_provenance_sidecar,
)

logger = logging.getLogger(__name__)


# Accession may carry a SwissProt isoform suffix (P04637-2); without it the
# varsplic entries fell through to the generic parser and the merged header
# came out as `sp|sp|P04637-2|P53_HUMAN` with OS/OX/GN demoted to free text.
UNIPROT_RE = re.compile(
    r"^(?P<db>sp|tr)\|(?P<acc>[A-Z0-9]+(?:-\d+)?)\|(?P<name>\S+)\s*(?P<desc>.*)$",
    re.IGNORECASE,
)
OS_RE = re.compile(r"OS=([^=]+?)(?=\s+[A-Z]{2}=|$)")
OX_RE = re.compile(r"OX=(\d+)")
GN_RE = re.compile(r"GN=(\S+)")
PE_RE = re.compile(r"PE=(\d+)")
SV_RE = re.compile(r"SV=(\d+)")


def _parse_uniprot_header(header: str, tag: str) -> SourceRecord | None:
    """Parse UniProt metadata into a source record, or return None for other headers."""
    m = UNIPROT_RE.match(header)
    if not m:
        return None
    desc = m.group("desc")
    os_ = (OS_RE.search(desc) or [None, None])[1]
    ox_ = (OX_RE.search(desc) or [None, None])[1]
    gn_ = (GN_RE.search(desc) or [None, None])[1]
    pe_ = (PE_RE.search(desc) or [None, None])[1]
    sv_ = (SV_RE.search(desc) or [None, None])[1]

    # Description = everything before the first OS=/OX=/GN=/...
    core_desc = re.split(r"\s+[A-Z]{2}=", desc, maxsplit=1)[0].strip()

    return SourceRecord(
        tag=tag,
        accession=m.group("acc"),
        gene=gn_,
        description=core_desc,
        organism=os_.strip() if os_ else None,
        tax_id=int(ox_) if ox_ else None,
        pe=int(pe_) if pe_ else None,
        sv=int(sv_) if sv_ else None,
        # Preserve the original SP/TR entry name (e.g. "P53_HUMAN") so the
        # downstream merged header can emit `sp|P04637|P53_HUMAN` rather than
        # a synthesised `sp|P04637|P04637_TP53`. Downstream search engines
        # (DIA-NN, SAGE) key off the accession, but reviewers recognise the
        # canonical entry-name format.
        extra={"entry_name": m.group("name")},
    )


def _parse_uhgp_header(header: str) -> SourceRecord | None:
    """UHGP: >MGYP000123 description"""
    m = re.match(r"^(MGYP\w+)\s*(.*)$", header)
    if not m:
        return None
    return SourceRecord(
        tag="UH",
        accession=m.group(1),
        description=m.group(2).strip(),
    )


def _parse_gencode_header(header: str) -> SourceRecord | None:
    """GENCODE: >ENSP00000269305.4|ENST...|ENSG...|gene_name|desc"""
    if not header.startswith("ENSP"):
        return None
    parts = header.split("|")
    return SourceRecord(
        tag="GENCODE",
        accession=parts[0],
        gene=parts[3] if len(parts) > 3 else None,
        description=parts[-1] if len(parts) > 1 else "",
    )


def _parse_generic_header(header: str, tag: str) -> SourceRecord:
    """Preserve the first header token as accession and the remaining description."""
    first = header.split()[0]
    desc = header[len(first) :].strip()
    return SourceRecord(tag=tag, accession=first, description=desc)


def parse_header(header: str, tag_hint: str | None = None) -> SourceRecord:
    """
    Parse a FASTA header line (without leading ``>``).

    ``tag_hint`` forces the source tag (e.g. "SP" or "UH"); otherwise the tag
    is inferred from the header format.
    """
    header = header.strip()
    if header.startswith("sp|"):
        r = _parse_uniprot_header(header, "SP")
        if r:
            return r
    if header.startswith("tr|"):
        r = _parse_uniprot_header(header, "TR")
        if r:
            return r
    if header.startswith("MGYP"):
        r = _parse_uhgp_header(header)
        if r:
            return r
    if header.startswith("ENSP"):
        r = _parse_gencode_header(header)
        if r:
            return r
    return _parse_generic_header(header, tag_hint or "UNKNOWN")


def iter_fasta(path: Path, tag: str) -> list[tuple[SourceRecord, str]]:
    """
    Yield (SourceRecord, sequence) pairs from a FASTA, tagging every record
    with the given source label.
    """
    out: list[tuple[SourceRecord, str]] = []
    header: str | None = None
    seq_parts: list[str] = []

    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    rec = parse_header(header, tag_hint=tag)
                    rec = replace(rec, tag=tag)  # enforce the caller's tag
                    out.append((rec, "".join(seq_parts).upper()))
                header = line[1:]
                seq_parts = []
            elif line:
                seq_parts.append(line)
        if header is not None:
            rec = replace(parse_header(header, tag_hint=tag), tag=tag)
            out.append((rec, "".join(seq_parts).upper()))
    return out


def build_lake(
    sources: dict[str, Path],
    out_fasta: Path,
    out_provenance: Path,
    dedup_mode: str = "merge",
    priority: tuple[str, ...] = DEFAULT_PRIORITY,
    embed_sources: bool = False,
    compress: bool = False,
) -> dict:
    """
    Build a merged lake FASTA + provenance sidecar from multiple sources.

    Parameters
    ----------
    sources
        ``{tag: fasta_path}``. Tags typically "SP", "TR", "UH", etc.
    out_fasta
        Output merged FASTA.
    out_provenance
        Output provenance JSONL (optionally .gz).
    dedup_mode
        "merge" (default, v5.2) or "priority-win" (legacy).
    priority
        Tag priority tuple. Default: SP, TR, UH, GM, GENCODE, GNOMAD, NCBI.
    embed_sources
        If True, append ``[FLlake:sources=...]`` tag to FASTA headers.
    compress
        If True, write sidecar as JSONL.gz.

    Returns
    -------
    dict
        Statistics of the merge: input, unique and per-source counts and timings.
    """
    if dedup_mode not in ("merge", "priority-win"):
        raise ValueError(f"unknown dedup_mode: {dedup_mode!r}")

    t0 = time.time()
    # seq_hash -> list[SourceRecord] AND the (first) sequence string for writing
    groups: dict[str, list[SourceRecord]] = defaultdict(list)
    seqs: dict[str, str] = {}
    per_source_counts: dict[str, int] = defaultdict(int)
    total_input = 0

    for tag, path in sources.items():
        recs = iter_fasta(Path(path), tag)
        per_source_counts[tag] = len(recs)
        total_input += len(recs)
        for src, seq in recs:
            h = sequence_hash(seq)
            groups[h].append(src)
            seqs.setdefault(h, seq)

    # Merge / pick-one per group
    merged: list[MergedRecord] = []
    for h, srcs in groups.items():
        if dedup_mode == "merge":
            merged.append(merge_records(seqs[h], srcs, priority=priority))
        else:  # priority-win — wrap in a MergedRecord with single source
            w = priority_win(srcs, priority=priority)
            merged.append(merge_records(seqs[h], [w], priority=priority))

    # Write FASTA
    out_fasta = Path(out_fasta)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(out_fasta, "w") as f:
        for rec in merged:
            header = format_fasta_header(rec, embed_sources=embed_sources)
            f.write(f">{header}\n")
            s = seqs[rec.seq_hash]
            for i in range(0, len(s), 60):
                f.write(s[i : i + 60] + "\n")

    # Write provenance sidecar
    out_provenance = Path(out_provenance)
    out_provenance.parent.mkdir(parents=True, exist_ok=True)
    write_provenance_sidecar(out_provenance, merged, compress=compress)

    # Summary stats
    n_shared = sum(1 for rec in merged if len(rec.sources) > 1)
    n_conflicts = sum(1 for rec in merged if rec.conflicts)
    n_ambiguous = sum(1 for rec in merged if rec.organism == "ambiguous")
    # Per-source contribution (how many output records had this source)
    contributed: dict[str, int] = defaultdict(int)
    for rec in merged:
        for s in rec.sources:
            contributed[s.tag] += 1

    stats = {
        "dedup_mode": dedup_mode,
        "n_input_sequences": total_input,
        "n_unique_sequences": len(merged),
        "dedup_ratio": round(total_input / max(len(merged), 1), 3),
        "n_shared_across_sources": n_shared,
        "n_with_conflicts": n_conflicts,
        "n_ambiguous_organism": n_ambiguous,
        "per_source_input": dict(per_source_counts),
        "per_source_contributed": dict(contributed),
        "runtime_seconds": round(time.time() - t0, 2),
    }
    logger.info(
        "Lake built: %d inputs → %d unique "
        "(%s mode, %.1fx dedup, %d collapsed across ≥2 sources) in %.1fs",
        total_input,
        len(merged),
        dedup_mode,
        stats["dedup_ratio"],
        n_shared,
        stats["runtime_seconds"],
    )
    return stats
