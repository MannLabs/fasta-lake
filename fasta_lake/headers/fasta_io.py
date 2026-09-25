"""
FastaLake v4 — FASTA I/O with full header preservation.

This module provides header-preserving alternatives to the functions in helpers.py.
The original helpers.py functions are kept unchanged for backward compatibility.

Key difference from v3 helpers.py:
- load_fasta_proteins() in helpers.py truncates headers to first token
- load_fasta_full() here preserves the complete header line
- write_fasta_proteins() here writes full headers, not just protein IDs
"""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def load_fasta_full(fasta_file: str | Path) -> dict[str, tuple[str, str]]:
    """Load FASTA file preserving full headers.

    Parameters
    ----------
    fasta_file : path
        Input FASTA file.

    Returns
    -------
    dict mapping protein_id → (sequence, full_header)
        protein_id is the first whitespace-delimited token (accession).
        full_header is the complete header line without '>'.
        sequence is the concatenated sequence.
    """
    proteins: dict[str, tuple[str, str]] = {}
    current_id: str | None = None
    current_header: str | None = None
    current_seq: list[str] = []

    def _flush(line_no: int) -> None:
        # Same contract as helpers.load_fasta_proteins: duplicates and empty
        # records are input defects, reported rather than silently dropped.
        """Store the pending FASTA record, rejecting empty sequences and duplicate IDs."""
        if current_id is None:
            return
        seq = "".join(current_seq)
        if not seq:
            raise ValueError(
                f"{fasta_file}: record {current_id!r} (before line {line_no}) has no sequence"
            )
        if current_id in proteins:
            raise ValueError(
                f"{fasta_file}: duplicate accession {current_id!r} (second occurrence before line "
                f"{line_no}); FASTA accessions must be unique"
            )
        proteins[current_id] = (seq, current_header)

    with open(fasta_file) as f:
        line_no = 0
        for line_no, line in enumerate(f, 1):
            if line.startswith(">"):
                _flush(line_no)
                full = line[1:].strip()
                if not full:
                    raise ValueError(
                        f"{fasta_file}:{line_no}: empty FASTA header ('>' with no accession)"
                    )
                current_id = full.split()[0]
                current_header = full
                current_seq = []
            else:
                current_seq.append(line.strip())
        _flush(line_no)

    logger.info("  Loaded %s proteins with full headers", f"{len(proteins):,}")
    return proteins


def load_fasta_proteins(fasta_file: str | Path) -> dict[str, str]:
    """Load FASTA file returning protein_id → sequence (backward compatible).

    This is equivalent to helpers.load_fasta_proteins() but uses the full
    loader internally. Provided for drop-in compatibility.
    """
    full = load_fasta_full(fasta_file)
    return {pid: seq for pid, (seq, _header) in full.items()}


def write_fasta(
    proteins: dict[str, tuple[str, str]],
    output_path: str | Path,
    line_width: int = 60,
) -> int:
    """Write FASTA file with full headers.

    Parameters
    ----------
    proteins : dict
        Mapping protein_id → (sequence, full_header).
    output_path : path
        Output file path.
    line_width : int
        Sequence line width (default 60).

    Returns
    -------
    int
        Number of proteins written.
    """
    count = 0
    with open(output_path, "w") as f:
        for protein_id in sorted(proteins):
            seq, header = proteins[protein_id]
            f.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i : i + line_width] + "\n")
            count += 1
    return count


def write_fasta_simple(
    protein_ids: list[str],
    proteins: dict[str, tuple[str, str]],
    output_path: str | Path,
    line_width: int = 60,
) -> int:
    """Write a subset of proteins to FASTA with full headers.

    Parameters
    ----------
    protein_ids : list
        List of protein IDs to write (order preserved).
    proteins : dict
        Full protein dict from load_fasta_full().
    output_path : path
        Output file path.

    Returns
    -------
    int
        Number of proteins written.
    """
    count = 0
    with open(output_path, "w") as f:
        for pid in protein_ids:
            if pid not in proteins:
                logger.warning("Protein %s not found in source", pid)
                continue
            seq, header = proteins[pid]
            f.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i : i + line_width] + "\n")
            count += 1
    return count
