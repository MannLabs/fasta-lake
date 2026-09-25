"""
Shared helper functions for FASTA Lake.

These helpers originated in the tested standalone inference scripts
and now include explicit input checks. They form the foundation for all inference strategies.

Notes
-----
The Aho-Corasick automaton performs multi-pattern peptide matching. This
experimental Python interface has different parsing and inference rules from
the supported Rust workflow; sharing the matching algorithm does not establish
end-to-end equivalence. See docs/concepts/algorithms.md, section 4.4.
"""

from __future__ import annotations

import csv
import logging
import re
from collections import defaultdict
from pathlib import Path

import ahocorasick

from .dianovo import parse_dianovo_prediction

logger = logging.getLogger(__name__)


def normalize_il(seq: str) -> str:
    """
    Normalize isoleucine (I) to leucine (L) for mass spec equivalence.

    Parameters
    ----------
    seq : str
        Protein or peptide sequence.

    Returns
    -------
    str
        Sequence with all I replaced by L.
    """
    return seq.replace("I", "L")


def clean_sequence(seq: str) -> str:
    """
    Clean a protein sequence by uppercasing and removing stop codons.

    Parameters
    ----------
    seq : str
        Raw protein sequence.

    Returns
    -------
    str
        Cleaned sequence (uppercase, no ``*``).
    """
    return seq.upper().replace("*", "")


def get_db_type(protein_id: str) -> str:
    """
    Classify a protein by its database source based on ID prefix.

    Parameters
    ----------
    protein_id : str
        Protein identifier (first token of FASTA header).

    Returns
    -------
    str
        One of ``'Human'``, ``'UHGG'``, ``'GMGC'``, ``'Assembly'``, ``'smORF'``.
    """
    if protein_id.startswith("sp|") or protein_id.startswith("tr|"):
        return "Human"
    elif protein_id.startswith("MGYG"):
        return "UHGG"
    elif protein_id.startswith("GMGC"):
        return "GMGC"
    elif protein_id.startswith("MuPr_Assembly_"):
        return "Assembly"
    else:
        return "smORF"


_LAKE_TAG_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_-]{0,15}$")


def get_species_id(protein_id: str) -> str | None:
    """
    Extract a species/genome identifier from a protein ID, or ``None``.

    Only identifiers that carry a species are resolved. Everything else is
    ``None`` so that callers can see, count and refuse unresolved proteins
    instead of pooling them into one pseudo-species (which turned
    ``species_budget`` into razor on every content-hash lake).

    Parameters
    ----------
    protein_id : str
        Protein identifier, either a raw catalogue accession or a FastaLake
        lake header ``<tag>|<accession>[|<entry>]``.

    Returns
    -------
    str or None
        - UHGG: ``MGYG000000001_00001`` → ``MGYG000000001`` (genome)
        - GMGC: ``None``. A GMGC unigene id (``GMGC10.001_001_001.GENE``) names
          the unigene, not an organism; taking its middle field made every
          protein its own "species". Provide a taxonomy map.
        - UniProt: ``sp|P12345|GFP_HUMAN`` → ``HUMAN`` (entry-name organism)
        - Lake header: ``uh|MGYG000000001_00001|…`` → ``MGYG000000001``
        - Anything else (assembly proteins, smORFs, content-hash ids): ``None``;
          provide a taxonomy map for these.
    """
    fields = protein_id.split("|")
    # FastaLake lake headers prefix a short source tag; it is not a species.
    if len(fields) >= 2 and fields[0].lower() not in ("sp", "tr") and _LAKE_TAG_RE.match(fields[0]):
        fields = fields[1:]
    if fields[0].lower() in ("sp", "tr"):
        if len(fields) >= 3:
            entry = fields[2].split()[0] if fields[2].strip() else ""
            if "_" in entry and entry.rsplit("_", 1)[1]:
                return entry.rsplit("_", 1)[1]
        return None
    accession = fields[0]
    if accession.startswith("MGYG"):
        return accession.split("_")[0]
    return None


def count_fasta_headers(path: str | Path) -> int:
    """Number of '>' records in a FASTA file, read with the file closed afterwards."""
    with open(path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def load_fasta_proteins(fasta_file: str | Path) -> dict[str, str]:
    """
    Load proteins from a FASTA file into a dict.

    Parameters
    ----------
    fasta_file : str or Path
        Path to FASTA file.

    Returns
    -------
    dict[str, str]
        Mapping of protein ID (first header token) → sequence.
    """
    fasta_file = Path(fasta_file)
    logger.info("Loading FASTA: %s...", fasta_file.name)

    proteins: dict[str, str] = {}
    current_header: str | None = None
    current_seq: list[str] = []

    def _flush(line_no: int) -> None:
        """Pass the pending FASTA record to the shared record validator."""
        if current_header is None:
            return
        _store_record(proteins, current_header, "".join(current_seq), fasta_file, line_no)

    with open(fasta_file) as f:
        line_no = 0
        for line_no, line in enumerate(f, 1):
            if line.startswith(">"):
                _flush(line_no)
                current_header = _accession_from_header(line, fasta_file, line_no)
                current_seq = []
            else:
                current_seq.append(line.strip())
        _flush(line_no)

    logger.info("  Loaded %s proteins", f"{len(proteins):,}")
    return proteins


def _accession_from_header(line: str, path, line_no: int) -> str:
    """Read the first header token, reporting the source line for an empty header."""
    tokens = line[1:].split()
    if not tokens:
        raise ValueError(f"{path}:{line_no}: empty FASTA header ('>' with no accession)")
    return tokens[0]


def _store_record(proteins: dict, accession: str, sequence: str, path, line_no: int) -> None:
    """Refuse duplicate accessions and empty sequences instead of losing records.

    A dict keyed by accession silently kept the LAST of two records with the
    same accession and dropped empty ones, so a concatenated database lost
    proteins with no count to show for it. Both are input defects that must
    be reported, not repaired.
    """
    if not sequence:
        raise ValueError(f"{path}: record {accession!r} (before line {line_no}) has no sequence")
    if accession in proteins:
        raise ValueError(
            f"{path}: duplicate accession {accession!r} (second occurrence before line "
            f"{line_no}); FASTA accessions must be unique"
        )
    proteins[accession] = sequence


def require_columns(
    fieldnames,
    required: tuple[str, ...] | list[str],
    path,
    *,
    any_of: tuple[str, ...] | list[str] = (),
) -> None:
    """Raise ValueError when a table lacks the columns a reader depends on.

    Readers used to fill a missing column with a default (``1`` for a q value,
    ``0`` for an intensity, ``""`` for a sequence), so a renamed column produced
    an empty result and exit 0. Naming the missing and the present columns is
    the only honest answer.
    """
    fields = list(fieldnames or [])
    missing = [c for c in required if c not in fields]
    if any_of and not any(c in fields for c in any_of):
        missing.append(" or ".join(any_of))
    if missing:
        raise ValueError(
            f"{path}: missing required column(s) {', '.join(missing)}; "
            f"found: {', '.join(fields) if fields else '(no header)'}"
        )


PEPTIDE_COLUMNS = ("peptide_prediction_detokenized_unmodified", "sequence")


def collect_peptides_with_scores(
    alphanovo_csv: str | Path,
    min_length: int = 9,
    max_length: int = 50,
) -> dict[str, float]:
    """
    Collect peptides with their best AlphaNovo scores from a CSV file.

    Parameters
    ----------
    alphanovo_csv : str or Path
        Path to a predictions CSV: AlphaNovo columns
        ``peptide_prediction_detokenized_unmodified`` and ``score``, the generic
        ``sequence`` and ``score`` pair, or DIANovo ``pred_seq`` and ``pred_prob``.
        Any other layout is refused with the missing column named; a row whose
        score cannot be parsed is skipped and counted, never given a default.
    min_length : int
        Minimum peptide length (inclusive).
    max_length : int
        Maximum peptide length (inclusive).

    Returns
    -------
    dict[str, float]
        Mapping of I/L-normalized peptide sequence → best score.
    """
    alphanovo_csv = Path(alphanovo_csv)
    logger.info("Collecting peptides from %s...", alphanovo_csv.name)

    peptide_scores: dict[str, float] = {}
    invalid_dianovo = 0
    invalid_score = 0

    with open(alphanovo_csv) as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames or []
        dianovo = "pred_seq" in fields
        if dianovo:
            require_columns(fields, ("pred_seq", "pred_prob"), alphanovo_csv)
            peptide_col = "pred_seq"
        else:
            require_columns(fields, ("score",), alphanovo_csv, any_of=PEPTIDE_COLUMNS)
            peptide_col = next(c for c in PEPTIDE_COLUMNS if c in fields)
        for row in reader:
            if dianovo:
                try:
                    peptide, score = parse_dianovo_prediction(
                        row.get("pred_seq") or "", row.get("pred_prob") or ""
                    )
                except ValueError:
                    invalid_dianovo += 1
                    continue
            else:
                peptide = (row.get(peptide_col) or "").strip().upper()
                try:
                    score = float(row["score"])
                except (ValueError, TypeError):
                    invalid_score += 1
                    continue

            if not peptide or len(peptide) < min_length or len(peptide) > max_length:
                continue

            if score <= 0.0:
                continue

            peptide_norm = normalize_il(peptide)
            if peptide_norm not in peptide_scores or score > peptide_scores[peptide_norm]:
                peptide_scores[peptide_norm] = score

    if invalid_dianovo:
        logger.warning("DIANovo: excluded %d unresolved/invalid prediction rows", invalid_dianovo)
    if invalid_score:
        logger.warning("Excluded %d rows whose score could not be parsed", invalid_score)
    logger.info("  Collected %s peptides", f"{len(peptide_scores):,}")
    return peptide_scores


def map_peptides_to_proteins(
    proteins: dict[str, str],
    peptides: set[str] | dict[str, float],
) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """
    Map peptides to proteins using Aho-Corasick multi-pattern matching.

    Parameters
    ----------
    proteins : dict[str, str]
        Mapping of protein ID → sequence.
    peptides : set[str] or dict[str, float]
        Peptide sequences to search for (I/L-normalized).

    Returns
    -------
    pep2prot : dict[str, set[str]]
        Mapping of peptide → set of protein IDs containing it.
    prot2pep : dict[str, set[str]]
        Mapping of protein ID → set of peptides found in it.
    """
    peptide_keys = peptides if isinstance(peptides, set) else set(peptides.keys())

    logger.info("Building Aho-Corasick automaton...")
    automaton = ahocorasick.Automaton()
    for peptide in peptide_keys:
        automaton.add_word(peptide, peptide)
    if not peptide_keys:
        logger.warning("No peptides to map — returning empty mappings")
        return {}, {}
    automaton.make_automaton()

    logger.info("Mapping peptides to proteins...")
    pep2prot: dict[str, set[str]] = defaultdict(set)
    prot2pep: dict[str, set[str]] = defaultdict(set)

    for protein_id, sequence in proteins.items():
        seq_norm = normalize_il(clean_sequence(sequence))

        for _end_pos, peptide in automaton.iter(seq_norm):
            pep2prot[peptide].add(protein_id)
            prot2pep[protein_id].add(peptide)

    logger.info("  Mapped %s peptides", f"{len(pep2prot):,}")
    logger.info("  %s proteins with peptide evidence", f"{len(prot2pep):,}")

    return dict(pep2prot), dict(prot2pep)
