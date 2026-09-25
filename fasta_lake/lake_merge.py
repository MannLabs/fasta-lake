"""
Annotation-aware dedup for lake building (v5.2, hash upgraded to SHA256 in v5.3).

Replaces priority-based winner-take-all with field-level merge across sources.
When two sources hold the same sequence (by SHA256 hash of the normalized
amino-acid string), this module merges their annotations under configurable
rules and emits a provenance sidecar so the choice is auditable.

The merged representative is written to the FASTA; the full per-source record
is serialised to JSONL (one line per sequence-hash) alongside it.

The hashing contract (normalization + algorithm) lives in
``fasta_lake.hashing``. It is shared byte-for-byte with the Rust lake_builder
so that all stages of the pipeline agree on identity.

Design doc: ``DESIGN_dedup.md``.
"""

from __future__ import annotations

import gzip
import json
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path

from fasta_lake.gzip_io import open_gzip_text
from fasta_lake.hashing import sequence_hash  # re-exported below

# Source tags in default priority order. Lower index = higher priority.
DEFAULT_PRIORITY: tuple[str, ...] = (
    "SP",  # SwissProt
    "TR",  # TrEMBL
    "UH",  # UHGP (bacterial gut microbiome)
    "GM",  # GMGC (gut microbial genes)
    "GENCODE",
    "GNOMAD",
    "NCBI",
    "UNKNOWN",
)

# Descriptions that are non-informative and should be dropped in favour of
# anything longer. Case-insensitive substring match.
UNINFORMATIVE_DESC_PATTERNS = re.compile(
    r"\b(hypothetical|uncharacteri[sz]ed|putative\s+uncharacteri[sz]ed|predicted\s+protein|"
    r"unknown\s+function|conserved\s+hypothetical|fragment)\b",
    re.IGNORECASE,
)


@dataclass
class SourceRecord:
    """One protein record from one source database."""

    tag: str  # "SP", "TR", "UH", ...
    accession: str
    gene: str | None = None
    description: str = ""
    organism: str | None = None
    tax_id: int | None = None
    pe: int | None = None  # UniProt Protein Existence (1=experimental, 5=uncertain)
    sv: int | None = None  # Sequence Version
    is_contaminant: bool = False
    is_canonical: bool = True
    extra: dict = field(default_factory=dict)


@dataclass
class MergedRecord:
    """Merged representative for one sequence-hash, with full source trail."""

    seq_hash: str
    length: int
    primary_tag: str
    primary_accession: str
    gene: str | None
    description: str
    organism: str | None  # "ambiguous" if sources disagree
    tax_ids: list[int]  # union, sorted
    pe: int | None  # min (= best evidence) across UniProt sources
    sv: int | None  # max (= latest) across sources
    is_contaminant: bool
    is_canonical: bool
    sources: list[SourceRecord]
    conflicts: dict[str, list]  # field -> divergent values, only when disagree


# ``sequence_hash`` is defined in ``fasta_lake.hashing`` (v5.3) and re-exported
# here so existing callers (``fasta_lake.lake_builder`` and tests that import
# from ``fasta_lake.lake_merge``) continue to work unchanged.
__all_hashing_reexports__ = ("sequence_hash",)


def _priority_index(tag: str, priority: tuple[str, ...]) -> int:
    """Return a source's priority position, sorting unknown tags after known ones."""
    try:
        return priority.index(tag)
    except ValueError:
        return len(priority)  # unknown tags sort last


def _merge_gene(sources: list[SourceRecord], priority: tuple[str, ...]) -> str | None:
    """First non-empty gene in priority order."""
    for s in sorted(sources, key=lambda r: _priority_index(r.tag, priority)):
        if s.gene:
            return s.gene
    return None


def _merge_description(sources: list[SourceRecord], priority: tuple[str, ...]) -> str:
    """Prefer the longest *informative* description; fall back to longest."""
    informative = [
        s
        for s in sources
        if s.description and not UNINFORMATIVE_DESC_PATTERNS.search(s.description)
    ]
    pool = informative if informative else [s for s in sources if s.description]
    if not pool:
        return ""
    pool.sort(key=lambda s: (-len(s.description), _priority_index(s.tag, priority)))
    return pool[0].description


def _merge_organism(sources: list[SourceRecord]) -> tuple[str | None, list[str] | None]:
    """If all sources agree → that organism. Else → 'ambiguous', list conflicts."""
    names = sorted({s.organism for s in sources if s.organism})
    if not names:
        return None, None
    if len(names) == 1:
        return names[0], None
    return "ambiguous", names


def _merge_tax_ids(sources: list[SourceRecord]) -> list[int]:
    """Collect distinct recorded taxonomy IDs in numeric order."""
    return sorted({s.tax_id for s in sources if s.tax_id is not None})


def _merge_pe(sources: list[SourceRecord]) -> int | None:
    """UniProt PE: LOWER is stronger evidence. Return min across sources."""
    pes = [s.pe for s in sources if s.pe is not None]
    return min(pes) if pes else None


def _merge_sv(sources: list[SourceRecord]) -> int | None:
    """Return the highest recorded sequence version, or None when all are absent."""
    svs = [s.sv for s in sources if s.sv is not None]
    return max(svs) if svs else None


def merge_records(
    seq: str,
    sources: list[SourceRecord],
    priority: tuple[str, ...] = DEFAULT_PRIORITY,
) -> MergedRecord:
    """
    Merge N source records that share an identical sequence into one.

    The primary representative is the highest-priority source's accession.
    Other annotation fields are merged per the rules in ``DESIGN_dedup.md``.

    Raises
    ------
    ValueError
        If the source list is empty.
    """
    if not sources:
        raise ValueError("merge_records requires at least one source")

    # Primary = highest-priority tag; tie-break by accession
    primary = min(
        sources,
        key=lambda s: (_priority_index(s.tag, priority), s.accession),
    )

    organism, organism_conflicts = _merge_organism(sources)

    conflicts: dict[str, list] = {}
    if organism_conflicts:
        conflicts["organism"] = organism_conflicts

    # Conflict detection for gene: only note if sources gave DIFFERENT non-null genes
    gene_values = sorted({s.gene for s in sources if s.gene})
    if len(gene_values) > 1:
        conflicts["gene"] = gene_values

    return MergedRecord(
        seq_hash=sequence_hash(seq),
        length=len(seq),
        primary_tag=primary.tag,
        primary_accession=primary.accession,
        gene=_merge_gene(sources, priority),
        description=_merge_description(sources, priority),
        organism=organism,
        tax_ids=_merge_tax_ids(sources),
        pe=_merge_pe(sources),
        sv=_merge_sv(sources),
        is_contaminant=any(s.is_contaminant for s in sources),
        is_canonical=any(s.is_canonical for s in sources),
        sources=list(sources),
        conflicts=conflicts,
    )


# ---------------------------------------------------------------------------
# Provenance sidecar I/O
# ---------------------------------------------------------------------------


def _primary_source(rec: MergedRecord) -> SourceRecord | None:
    """Return the SourceRecord that won primary status, or None if not found."""
    for s in rec.sources:
        if s.tag == rec.primary_tag and s.accession == rec.primary_accession:
            return s
    return None


def format_fasta_header(
    rec: MergedRecord,
    embed_sources: bool = False,
) -> str:
    """
    Format the merged FASTA header line (without leading ``>``).

    Rules for the downstream-visible fields (SAGE / DIA-NN see these):
      - Third pipe field (entry name like ``P53_HUMAN``) is preserved verbatim
        from the primary SP/TR source when available; otherwise omitted. We
        do NOT synthesise a pseudo entry name, because non-standard forms
        like ``P04637_TP53`` are neither valid SP entry names nor the
        accession, and confuse header-aware tooling.
      - ``OS=`` and ``OX=`` follow the primary source, even when other
        sources disagreed during the annotation merge. The full multi-source
        view (including the disagreement) stays in the provenance sidecar,
        not in the FASTA header where search engines can't use it.

    The ``[FLlake:...]`` suffix is opt-in — downstream search engines mostly
    tolerate it, but some header healers may strip it.
    """
    primary = _primary_source(rec)

    # Third pipe field: preserve original entry name if we have it.
    entry_name = primary.extra.get("entry_name") if primary else None
    if entry_name:
        head = f"{rec.primary_tag.lower()}|{rec.primary_accession}|{entry_name}"
    else:
        head = f"{rec.primary_tag.lower()}|{rec.primary_accession}"

    # OS / OX: prefer primary source's values. This is what downstream tools
    # expect and what reviewers grep for. Ambiguity lives in the sidecar.
    display_os = primary.organism if primary and primary.organism else None
    if display_os is None and rec.organism and rec.organism != "ambiguous":
        display_os = rec.organism
    display_ox = primary.tax_id if primary and primary.tax_id is not None else None
    if display_ox is None and rec.tax_ids:
        display_ox = rec.tax_ids[0]

    tail_bits: list[str] = []
    if rec.description:
        tail_bits.append(rec.description)
    if display_os:
        tail_bits.append(f"OS={display_os}")
    if display_ox is not None:
        tail_bits.append(f"OX={display_ox}")
    if rec.gene:
        tail_bits.append(f"GN={rec.gene}")
    if rec.pe is not None:
        tail_bits.append(f"PE={rec.pe}")
    if rec.sv is not None:
        tail_bits.append(f"SV={rec.sv}")

    if embed_sources:
        srctags = ",".join(sorted({s.tag for s in rec.sources}))
        txs = ",".join(str(t) for t in rec.tax_ids) or "-"
        tail_bits.append(f"[FLlake:sources={srctags};tax={txs}]")

    tail = " ".join(tail_bits)
    return f"{head} {tail}".strip() if tail else head


def _record_to_jsonable(rec: MergedRecord) -> dict:
    """Convert MergedRecord to a JSON-serialisable dict (flattens dataclasses).

    Adds two derived convenience fields that make the sidecar directly
    grep/jq-friendly for the most common question users ask — "which
    sequences had more than one source contribute the exact same bytes?":

      n_sources      — len(sources); makes `jq 'select(.n_sources > 1)'` work
      alt_accessions — accessions from non-primary sources; the "also known as"
                       list for byte-identical collapses.

    The primary accession stays in ``primary_accession``; these fields are
    additive, not replacements, and are stripped by ``read_provenance_sidecar``
    on round-trip so they don't need to be dataclass fields.
    """
    d = asdict(rec)
    d["n_sources"] = len(rec.sources)
    d["alt_accessions"] = [s.accession for s in rec.sources if s.accession != rec.primary_accession]
    return d


def write_provenance_sidecar(
    out_path: Path,
    records: list[MergedRecord],
    compress: bool = False,
) -> int:
    """
    Write merged records as JSONL (or JSONL.gz if compress=True).

    Returns the number of records written.
    """
    out_path = Path(out_path)
    opener = open_gzip_text if compress else (lambda p: open(p, "w"))
    with opener(out_path) as f:
        for rec in records:
            f.write(json.dumps(_record_to_jsonable(rec), separators=(",", ":")) + "\n")
    return len(records)


def read_provenance_sidecar(in_path: Path) -> list[MergedRecord]:
    """Read a provenance sidecar JSONL (or JSONL.gz) back into MergedRecords."""
    in_path = Path(in_path)
    opener = (
        (lambda p: gzip.open(p, "rt", encoding="utf-8"))
        if str(in_path).endswith(".gz")
        else (lambda p: open(p))
    )
    records: list[MergedRecord] = []
    with opener(in_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            d = json.loads(line)
            src = [SourceRecord(**s) for s in d.pop("sources", [])]
            # Strip derived convenience fields that were added for grep/jq
            # queries but are not dataclass fields.
            d.pop("n_sources", None)
            d.pop("alt_accessions", None)
            rec = MergedRecord(sources=src, **d)
            records.append(rec)
    return records


# ---------------------------------------------------------------------------
# Priority-win back-compat mode
# ---------------------------------------------------------------------------


def priority_win(
    sources: list[SourceRecord],
    priority: tuple[str, ...] = DEFAULT_PRIORITY,
) -> SourceRecord:
    """
    Legacy behaviour: pick one winner by priority, discard other sources.

    Kept for ``--dedup-mode priority-win`` back-compat and regression testing.
    """
    if not sources:
        raise ValueError("priority_win requires at least one source")
    return min(
        sources,
        key=lambda s: (_priority_index(s.tag, priority), s.accession),
    )
