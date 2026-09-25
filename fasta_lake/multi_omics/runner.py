"""
FL_MO Stage-3 runner: wire the evidence gate to on-disk inputs/outputs.

This is the file-level orchestration that takes the verified gate logic in
`evidence_gate` and applies it to a single sample's real inputs:

  - the sample's de-novo-extracted FASTA (Stage-3 output, evidence (a)),
  - the sample's per-assembly TPM TSV (metaG_tpm / metaT_tpm columns),
  - the sample's MEGAHIT predicted-genes FASTA,
  - the FL_MO hash_to_rep table,
  - the FL_MO cluster97 representative FASTA (for verbatim sequence copy).

It writes:
  - <output_fasta>           de-novo proteins UNION gate-admitted rescued reps,
  - <output_lookup>          evidence_lookup.tsv (full attribution),
  - <output_stats>           JSON stats sidecar (optional).

Default gate is OR (the FL_MO default); AND is the low-FDR option. Default runs
of the de-novo-only pipeline are untouched: this code only runs when explicitly
invoked via the CLI flag or this function.

No new dependencies; headers are copied byte-verbatim from the rep FASTA.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from fasta_lake.multi_omics.evidence_gate import (
    EVIDENCE_LOOKUP_COLUMNS,
    Gate,
    build_evidence_lookup,
    gate_members,
)

__all__ = ["run_fl_mo_sample", "load_tpm_tsv", "index_fasta", "load_hash_to_rep"]

# TPM TSV column names (matches the project's *_TPM.tsv / *_TPM_annotated.tsv).
_TPM_ID_COLS = ("prodigal_protein", "protein", "gene_id", "k141")
_TPM_METAG_COL = "metaG_tpm"
_TPM_METAT_COL = "metaT_tpm"


def _f(x: str) -> float:
    """Parse a legacy TPM value, treating blank or malformed text as zero.

    This fallback belongs to the historical gate; the current exact-input
    parser preserves missing measurements instead.
    """
    try:
        return float(x) if x and x.strip() else 0.0
    except ValueError:
        return 0.0


def load_tpm_tsv(path):
    """Yield (assembly_protein_id, metaG_tpm, metaT_tpm) from a TPM TSV.

    The id column is auto-detected from `_TPM_ID_COLS`.
    """
    path = Path(path)
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        id_col = next((c for c in _TPM_ID_COLS if c in idx), None)
        if id_col is None:
            raise ValueError(
                f"{path}: no recognized id column (looked for {_TPM_ID_COLS}); header={header[:6]}"
            )
        if _TPM_METAG_COL not in idx or _TPM_METAT_COL not in idx:
            raise ValueError(f"{path}: missing {_TPM_METAG_COL}/{_TPM_METAT_COL} columns")
        ip, img, imt = idx[id_col], idx[_TPM_METAG_COL], idx[_TPM_METAT_COL]
        n = max(ip, img, imt)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= n:
                continue
            yield p[ip], _f(p[img]), _f(p[imt])


def index_fasta(path):
    """Return {first_whitespace_token_of_header: sequence_str} for a FASTA."""
    out: dict[str, str] = {}
    cur = None
    chunks: list[str] = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur is not None:
                    out[cur] = "".join(chunks)
                cur = line[1:].split(None, 1)[0].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if cur is not None:
            out[cur] = "".join(chunks)
    return out


def load_hash_to_rep(path):
    """Load the FL_MO hash_to_rep table: {sha256_hex: cluster_rep_id}.

    Tab-separated, first column = hash, second = rep id; a header line is
    skipped.
    """
    htr: dict[str, str] = {}
    with open(path) as f:
        first = f.readline()
        # if the first line looks like data (64-hex first token), keep it
        tok = first.split("\t", 1)[0].strip()
        if len(tok) == 64 and all(c in "0123456789abcdef" for c in tok.lower()):
            p = first.rstrip("\n").split("\t")
            if len(p) >= 2:
                htr[p[0]] = p[1]
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                htr[p[0]] = p[1]
    return htr


def _denovo_rep_ids(denovo_fasta):
    """First-token header ids from the de-novo extracted FASTA."""
    ids = set()
    with open(denovo_fasta) as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line[1:].split(None, 1)[0].strip())
    return ids


def _pull_rep_seqs(cluster_rep_fasta, needed):
    """Stream the rep FASTA, returning {rep_id: raw_sequence_bytes} for needed."""
    needed = set(needed)
    seqs: dict[str, bytes] = {}
    cur = None
    chunks: list[bytes] = []
    with open(cluster_rep_fasta, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                if cur is not None and cur in needed:
                    seqs[cur] = b"".join(chunks)
                cur = line[1:].split(None, 1)[0].decode("ascii", "replace")
                chunks = []
            elif cur in needed:
                chunks.append(line.rstrip(b"\n"))
        if cur is not None and cur in needed:
            seqs[cur] = b"".join(chunks)
    return seqs


def run_fl_mo_sample(
    *,
    sample: str,
    denovo_fasta: str,
    tpm_tsv: str,
    assembly_fasta: str,
    hash_to_rep: str,
    cluster_rep_fasta: str,
    output_fasta: str,
    output_lookup: str,
    output_stats: str | None = None,
    gate: str = Gate.OR,
    metag_top_pct: float = 1.0,
    metat_abs: float = 1.0,
) -> dict:
    """Run the FL_MO gate for one sample and write the augmented FASTA + lookup.

    Returns the stats dict (also written to `output_stats` if given).
    """
    t0 = time.time()

    denovo_ids = _denovo_rep_ids(denovo_fasta)
    tpm_rows = list(load_tpm_tsv(tpm_tsv))
    asm = index_fasta(assembly_fasta)
    htr = load_hash_to_rep(hash_to_rep)

    rows, stats = build_evidence_lookup(
        denovo_rep_ids=denovo_ids,
        tpm_rows=tpm_rows,
        assembly_seqs=asm,
        hash_to_rep=htr,
        metag_top_pct=metag_top_pct,
        metat_abs=metat_abs,
    )

    # Write evidence_lookup.tsv (full attribution, both sources)
    with open(output_lookup, "w") as f:
        f.write("\t".join(EVIDENCE_LOOKUP_COLUMNS) + "\n")
        for r in rows:
            f.write(r.to_tsv() + "\n")

    # Members admitted under the chosen gate
    keep = gate_members(rows, gate)

    # The augmented FASTA = de-novo proteins (verbatim) + gate-admitted rescued
    # reps not already in the de-novo set. `keep` is a set, so the rescued block
    # is written in sorted id order to keep the FASTA byte-identical across runs.
    rescued_keep = sorted(rep for rep in keep if rep not in denovo_ids)
    rep_seqs = _pull_rep_seqs(cluster_rep_fasta, rescued_keep) if rescued_keep else {}

    # Two-pass write: de-novo records that pass the gate (verbatim), then the
    # gate-admitted rescued reps. Under OR-gate every de-novo rep is kept
    # (has_denovo => in_gate(OR)); under AND-gate only de-novo reps that also
    # have TPM evidence are kept (`keep` already encodes this).
    n_denovo = 0
    n_rescued = 0
    with open(output_fasta, "wb") as fout:
        cur_keep = False
        with open(denovo_fasta, "rb") as fin:
            for line in fin:
                if line.startswith(b">"):
                    rep = line[1:].split(None, 1)[0].decode("ascii", "replace")
                    cur_keep = rep in keep
                    if cur_keep:
                        n_denovo += 1
                        fout.write(line)
                elif cur_keep:
                    fout.write(line)
        for rep in rescued_keep:
            seq = rep_seqs.get(rep)
            if seq is None:
                continue
            fout.write(f">{rep} [evidence:tpm_rescue]\n".encode("ascii"))
            for i in range(0, len(seq), 60):
                fout.write(seq[i : i + 60])
                fout.write(b"\n")
            n_rescued += 1

    stats.update(
        sample=sample,
        gate=gate,
        fasta_denovo_written=n_denovo,
        fasta_rescued_written=n_rescued,
        fasta_total=n_denovo + n_rescued,
        gate_members=len(keep),
        elapsed_s=round(time.time() - t0, 2),
    )
    if output_stats:
        with open(output_stats, "w") as f:
            json.dump(stats, f, indent=2)
    return stats
