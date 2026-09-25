"""
Contaminant handling (v5.1).

Common proteomics contaminants — keratins, trypsin, BSA, caseins, etc. — must
be tagged explicitly so they're identifiable in search results and can be
excluded from target/entrapment FDR accounting.

This module provides:

- A bundled list of common contaminant UniProt accessions (``DEFAULT_CONTAMINANTS``).
- ``tag_contaminants()``: rewrite FASTA headers with ``CON_`` prefix where needed.
- ``is_contaminant()``: O(1) check from protein accession.
- ``strip_contaminants_from_fdr()``: exclusion helper for FDR tables.

The bundled list is NOT a full cRAP — it's a minimal conservative set
(~50 accessions). Users needing the full cRAP should download it from
thegpm.org and pass via ``--contaminant-db``.

All contaminant lookups are by accession, not by sequence. This means a
contaminant from one source (e.g. SwissProt P00761 trypsin) is correctly
tagged even if the same sequence appears in another source.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

# Common proteomics contaminants (UniProt accessions).
# Minimal conservative set — users may extend via --contaminant-db.
# Sourced from the GPM's cRAP minimal set; bundled as ID list (not sequences).
DEFAULT_CONTAMINANTS: frozenset[str] = frozenset(
    {
        # Trypsin — bovine
        "P00761",
        # Bovine serum albumin
        "P02769",
        # Streptavidin (bacterial, high-affinity tag)
        "P22629",
        # Carbonic anhydrase (bovine, standard)
        "P00921",
        # Glyceraldehyde-3-phosphate dehydrogenase (rabbit, standard)
        "P46406",
        # Keratin family — human (skin/hair contamination)
        "P04264",  # K1
        "P35908",  # K2e
        "P13645",  # K10
        "P35527",  # K9
        "P04259",  # K6A
        "P13647",  # K5
        "P05787",  # K8
        "P08727",  # K19
        "P08779",  # K16
        "P19013",  # K4
        "Q04695",  # K17
        # Keratin family — sheep/wool
        "P02535",  # sheep K1
        # Caseins (milk/cheese contamination, handled samples)
        "P02662",  # alpha-S1-casein bovine
        "P02663",  # alpha-S2-casein bovine
        "P02666",  # beta-casein bovine
        "P02668",  # kappa-casein bovine
        # Alpha-lactalbumin, beta-lactoglobulin (milk)
        "P00711",
        "P02754",
        # Ovalbumin (chicken egg)
        "P01012",
        # Lysozyme C (chicken egg white, common buffer)
        "P00698",
        # Myoglobin (standard, bovine)
        "P02192",
        # Immunoglobulin common chains (human IgG)
        "P01857",  # IGHG1
        "P01859",  # IGHG2
        "P01860",  # IGHG3
        "P01861",  # IGHG4
        "P01876",  # IGHA1
        "P01877",  # IGHA2
        "P01834",  # IGKC
        # Serum components (common contaminants in cell culture via FBS)
        "P12763",  # bovine alpha-2-HS-glycoprotein
        "P02070",  # bovine haemoglobin beta
        "P01966",  # bovine haemoglobin alpha
        # ProteaseMAX / cleavable surfactants — often as degraded protein
        # (none)
        # Tag proteins often seen (his-tag fused)
        # (none universal)
    }
)


@dataclass
class ContaminantDB:
    """In-memory contaminant accession database."""

    accessions: set[str] = field(default_factory=set)
    prefix: str = "CON_"

    def __post_init__(self) -> None:
        """Populate the default contaminant accessions when none were supplied."""
        if not self.accessions:
            self.accessions = set(DEFAULT_CONTAMINANTS)

    def is_contaminant(self, protein_id: str) -> bool:
        """True if this protein ID corresponds to a contaminant.

        Matches:
          - exact accession (e.g., ``P00761``)
          - SwissProt-style ``sp|ACC|NAME`` (extracts ACC)
          - already-tagged ``CON_ACC``
        """
        if protein_id.startswith(self.prefix):
            return True
        acc = _extract_accession(protein_id)
        return acc in self.accessions

    def load(self, path: str | Path) -> None:
        """Load additional accessions from a text file (one per line)."""
        p = Path(path)
        if not p.exists():
            return
        for line in p.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            self.accessions.add(line)

    def load_from_fasta(self, fasta: str | Path) -> None:
        """Extract accessions from a FASTA file."""
        p = Path(fasta)
        if not p.exists():
            return
        with open(p) as f:
            for line in f:
                if line.startswith(">"):
                    acc = _extract_accession(line[1:].split()[0])
                    if acc:
                        self.accessions.add(acc)


def _extract_accession(protein_id: str) -> str | None:
    """Extract UniProt-style accession from various header formats.

    >>> _extract_accession('sp|P00761|TRYP_BOVIN')
    'P00761'
    >>> _extract_accession('tr|Q12345|NAME')
    'Q12345'
    >>> _extract_accession('CON_P00761')
    'P00761'
    >>> _extract_accession('P00761')
    'P00761'
    """
    s = protein_id.strip()
    # CON_ / rev_ prefix
    for pref in ("CON_", "rev_", "DECOY_"):
        if s.startswith(pref):
            s = s[len(pref) :]
    # sp| / tr| prefix
    if s.startswith(("sp|", "tr|")):
        parts = s.split("|")
        if len(parts) >= 2:
            return parts[1]
    # Raw accession
    # UniProt accessions: 6-10 alphanumeric, starting with letter
    if 6 <= len(s) <= 10 and s[0].isalpha() and s.isalnum():
        return s
    return None


def tag_contaminants(
    proteins: dict[str, str],
    db: ContaminantDB | None = None,
) -> tuple[dict[str, str], list[str]]:
    """Rewrite protein IDs with ``CON_`` prefix where they match contaminants.

    Returns (tagged_proteins, tagged_accessions).
    """
    if db is None:
        db = ContaminantDB()
    tagged: dict[str, str] = {}
    contam_ids: list[str] = []
    for pid, seq in proteins.items():
        if db.is_contaminant(pid):
            new_id = pid if pid.startswith(db.prefix) else f"{db.prefix}{pid}"
            tagged[new_id] = seq
            contam_ids.append(new_id)
        else:
            tagged[pid] = seq
    return tagged, contam_ids


def is_contaminant_id(protein_id: str, db: ContaminantDB | None = None) -> bool:
    """Check if a protein_id is a contaminant. Cheap wrapper for FDR filters."""
    if db is None:
        db = ContaminantDB()
    return db.is_contaminant(protein_id)


def exclude_contaminants(
    rows: Iterable[dict],
    protein_field: str = "proteins",
    db: ContaminantDB | None = None,
) -> list[dict]:
    """Filter out PSM rows whose protein is a contaminant.

    For use in FDR calculation — contaminants should not count as targets
    or entrapment.
    """
    if db is None:
        db = ContaminantDB()
    out = []
    for r in rows:
        prot = r.get(protein_field, "")
        # Header may contain multiple proteins separated by ; or ,
        first = prot.split(";")[0].split(",")[0].strip()
        if db.is_contaminant(first):
            continue
        out.append(r)
    return out
