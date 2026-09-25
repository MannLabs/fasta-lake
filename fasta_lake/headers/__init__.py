"""
FastaLake v4 — FASTA header parsing, formatting, and healing.

This module handles the complete lifecycle of FASTA headers across the pipeline:
1. Auto-detection of header format from 13+ known patterns
2. Structured parsing into a universal ParsedHeader
3. Reformatting for specific search engines (DIA-NN, SAGE)
4. Healing entire FASTA files with format validation and QC

DIA-NN header rules (empirically verified 2026-03-26):
- For sp| and tr| prefixed accessions: uses GN= tag from description for gene name
- For ALL other accession formats: uses the FIRST WORD of description as gene name
- Ignores GN= for non-sp/tr entries
- Protein.Group column preserves full accession (first whitespace-delimited token)

SAGE header rules:
- Uses the full header line for protein name
- Gene extraction varies by version
- Generally more tolerant of non-standard formats

Design principle: the accession (first field before space) is NEVER modified.
It is the unique identifier and traceability key. Only the description is reformatted.
"""

from __future__ import annotations

import logging
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════════
# DATA MODEL
# ═══════════════════════════════════════════════════════════════════════════════


@dataclass
class ParsedHeader:
    """Universal representation of a parsed FASTA header.

    The `raw` field always contains the original header exactly as it appeared
    in the source FASTA (without the leading '>'). This ensures traceability.

    The `accession` is the first whitespace-delimited token — this is what
    search engines use as the protein identifier. It is NEVER modified.
    """

    raw: str  # original header verbatim (without >)
    accession: str  # first token (protein ID for search engines)
    gene: str | None = None  # gene name (e.g., ALB, dnaA, MGYG000005705_02541)
    description: str | None = None  # protein function/name
    organism: str | None = None  # species name
    source_db: str = "unknown"  # database of origin
    locus_tag: str | None = None  # locus tag (microbes)
    variant_info: str | None = None  # mutation description (variants)
    population: str | None = None  # population code (gnomAD)
    extra: dict = field(default_factory=dict)

    def format_for_diann(self) -> str:
        """Format header for DIA-NN compatibility.

        Rules:
        - sp|/tr| entries: keep as-is (DIA-NN uses GN= from description). A
          resolved gene is added as GN=; an unresolved one is NOT invented, the
          header is written without GN= and DIA-NN keys the entry by accession.
        - All others: gene name must be FIRST WORD of description; the accession
          stands in when the source has no gene concept (contig-local ids).
        - Accession is NEVER modified (traceability)
        """
        # sp| and tr| entries: DIA-NN's UniProt parser handles these correctly
        if self.source_db in ("swissprot", "trembl"):
            if "GN=" in self.raw or self.gene is None:
                return self.raw
            desc = self.raw[len(self.accession) :].strip()
            pe_pos = desc.find(" PE=")
            if pe_pos > 0:
                return f"{self.accession} {desc[:pe_pos]} GN={self.gene}{desc[pe_pos:]}"
            return f"{self.accession} {desc} GN={self.gene}"

        gene = self.gene or self.accession  # non-UniProt: accession is the id

        # All non-sp/tr entries: gene must be FIRST WORD of description
        # Build: ACCESSION GENE description_text GN=GENE
        parts = []
        parts.append(self.accession)

        # Description with gene as first word
        desc_parts = [gene]  # gene as first word
        if self.description:
            desc_parts.append(self.description)
        if self.variant_info:
            desc_parts.append(f"[variant {self.variant_info}]")
        if self.organism:
            desc_parts.append(f"OS={self.organism}")
        desc_parts.append(f"GN={gene}")

        parts.append(" ".join(desc_parts))
        return " ".join(parts)

    def format_for_sage(self) -> str:
        """Format header for SAGE compatibility.

        SAGE is more tolerant — it uses the full accession as protein ID
        and can extract gene names from standard UniProt format.
        For non-UniProt entries, adding GN= helps but is not strictly required.
        """
        # For SAGE, we can use a simpler format — just ensure GN= is present
        if "GN=" in self.raw:
            return self.raw
        if self.source_db in ("swissprot", "trembl") and self.gene is None:
            return self.raw  # never invent a UniProt gene symbol

        gene = self.gene or self.accession
        return f"{self.raw} GN={gene}"

    def format_for_msfragger(self) -> str:
        """Format header for MSFragger/Philosopher compatibility.

        MSFragger uses the accession (first token) as protein ID.
        Philosopher's PeptideProphet needs standard header formatting.
        For non-sp/tr entries, ensure gene is discoverable.
        Headers are kept as-is — the critical MSFragger requirement
        is reversed decoy sequences (rev_ prefix), which is handled
        at the FASTA level by ``prepare_fasta_for_engine()``.
        """
        # MSFragger is tolerant of header format like SAGE.
        # Main thing: ensure GN= present for gene-level reporting.
        if "GN=" in self.raw:
            return self.raw
        if self.source_db in ("swissprot", "trembl") and self.gene is None:
            return self.raw  # never invent a UniProt gene symbol

        gene = self.gene or self.accession
        return f"{self.raw} GN={gene}"


# ═══════════════════════════════════════════════════════════════════════════════
# PARSER: Auto-detect and parse any FASTA header format
# ═══════════════════════════════════════════════════════════════════════════════


def parse_header(raw: str) -> ParsedHeader:
    """Parse a FASTA header line (without leading '>') into a ParsedHeader.

    Auto-detects the format from 13+ known patterns and extracts all available
    metadata. Unknown formats are handled gracefully with the accession as gene.

    Parameters
    ----------
    raw : str
        The header line without the leading '>' character, stripped of newline.

    Returns
    -------
    ParsedHeader
        Structured representation with all extracted fields.
    """
    raw = raw.strip()
    accession = raw.split()[0] if raw else raw

    # ── 0a. Disease variants: sp|ACC_vNNN|GENE_VAR (must check before generic sp|) ──
    if raw.startswith("sp|") and "_VAR" in accession:
        return _parse_disease_variant(raw, accession)

    # ── 1. UniProt SwissProt: sp|ACC|NAME_HUMAN description GN=GENE ──
    if raw.startswith("sp|"):
        return _parse_uniprot(raw, accession, "swissprot")

    # ── 2. UniProt TrEMBL: tr|ACC|NAME_HUMAN description [GN=GENE] ──
    if raw.startswith("tr|"):
        return _parse_uniprot(raw, accession, "trembl")

    # ── 3. gnomAD variants: var-POP|ACC|GENE_Mut|pos description GN=GENE ──
    if raw.startswith("var-"):
        return _parse_gnomad(raw, accession)

    # ── 4a. Ensembl pep: ENSP... pep chromosome:... gene_symbol:XXX ──
    #    (must check before GENCODE — Ensembl pep has "pep" keyword, no pipes in accession)
    if raw.startswith("ENSP") and " pep " in raw:
        return _parse_ensembl_pep(raw, accession)

    # ── 4b. Ensembl mouse pep: ENSMUSP... pep chromosome:... ──
    if raw.startswith("ENSMUSP") and " pep " in raw:
        return _parse_ensembl_pep(raw, accession)

    # ── 4. GENCODE pipe-delimited: ENSP|ENST|ENSG|...|transcript-NNN|GENE|length ──
    if raw.startswith("ENSP") and "|" in accession:
        return _parse_gencode(raw, accession)

    # ── 5. Personal variants: ens|ENST_Mut|GENE description ──
    if raw.startswith("ens|"):
        return _parse_ens_variant(raw, accession)

    # ── 6. NCBI microbe: SPECIES|WP/YP_ACC|LOCUS [gene=X] description OS=... ──
    # Detect by pipe-delimited with WP_, YP_, NP_ in second field
    if "|" in accession:
        fields = accession.split("|")
        if len(fields) >= 2 and any(fields[1].startswith(p) for p in ("WP_", "YP_", "NP_")):
            return _parse_ncbi_microbe(raw, accession, fields)

    # ── 7. MGYG (UHGP): MGYG000000000_00000 description ──
    if raw.startswith("MGYG"):
        return _parse_mgyg(raw, accession)

    # ── 8. GMGC: GMGC10.xxx_xxx_xxx.UNKNOWN ──
    if raw.startswith("GMGC"):
        return _parse_simple_microbe(raw, accession, "gmgc")

    # ── 9. SMORF: SMORF_00000 description ──
    if raw.startswith("SMORF") or raw.startswith("smORF"):
        return _parse_simple_microbe(raw, accession, "smorf")

    # ── 10. Assembly: MuPr_Assembly_SAMP_00000 description ──
    if raw.startswith("MuPr_Assembly"):
        return _parse_simple_microbe(raw, accession, "assembly")

    # ── 12. GENSCAN: GENSCAN00000 pep scaffold:... ──
    if raw.startswith("GENSCAN"):
        return ParsedHeader(
            raw=raw,
            accession=accession,
            gene=None,
            description="computational prediction",
            source_db="genscan",
        )

    # ── 13. Prodigal: 1_1 # 1 # 1404 # 1 # ID=1_1;partial=10;... ──
    if " # " in raw and "ID=" in raw:
        return _parse_prodigal(raw, accession)

    # ── 14. SAMPLE|REF| prefix (lake pipeline artifact) ──
    if accession.startswith("SAMPLE|REF|") or "|REF|" in accession:
        return _parse_sample_ref(raw, accession)

    # ── 15. Lake hybrid: PROJECT [SWISSPROT] sp|ACC|NAME desc ──
    if "[SWISSPROT]" in raw or "[UHGG]" in raw or "[GMGC]" in raw:
        return _parse_lake_hybrid(raw, accession)

    # ── 16. MD5 hash (from clustering/inference output) ──
    if re.match(r"^[0-9a-f]{32}$", accession):
        return ParsedHeader(
            raw=raw,
            accession=accession,
            gene=accession,
            description=None,
            source_db="md5_hash",
        )

    # ── 17. Disease variants: sp|P02649_v21K|APOE_VAR desc ──
    if raw.startswith("sp|") and "_v" in accession and "_VAR" in accession:
        return _parse_disease_variant(raw, accession)

    # ── 18. Pipe-delimited isoforms: ACC|GENE|ISOFORM_N|tissue ──
    if "|" in accession:
        fields = accession.split("|")
        if len(fields) >= 3:
            candidate = fields[1].strip()
            if re.match(r"^[A-Z][A-Z0-9]{1,14}$", candidate):
                desc = raw[len(accession) :].strip() or None
                return ParsedHeader(
                    raw=raw,
                    accession=accession,
                    gene=candidate,
                    description=desc,
                    source_db="isoform_db",
                )

    # ── 13. Population-specific variants: tr|ACC|NAME|POP|MUT ──
    if raw.startswith("tr|") or raw.startswith("sp|"):
        # Already handled above, but catch edge cases
        return _parse_uniprot(raw, accession, "swissprot" if raw.startswith("sp|") else "trembl")

    # ── Fallback: unknown format ──
    desc = raw[len(accession) :].strip() or None
    # Try to extract gene from GN= if present anywhere
    m_gn = re.search(r"GN=(\S+)", raw)
    gene = m_gn.group(1) if m_gn else accession

    return ParsedHeader(
        raw=raw,
        accession=accession,
        gene=gene,
        description=desc,
        source_db="unknown",
    )


# ═══════════════════════════════════════════════════════════════════════════════
# FORMAT-SPECIFIC PARSERS
# ═══════════════════════════════════════════════════════════════════════════════


def _parse_uniprot(raw: str, accession: str, source_db: str) -> ParsedHeader:
    """Parse UniProt sp| or tr| format."""
    # Extract gene from GN= tag
    m_gn = re.search(r"GN=(\S+)", raw)
    gene = m_gn.group(1) if m_gn else None

    # Extract organism
    m_os = re.search(r"OS=(.+?)(?:\s+OX=|\s+GN=|\s+PE=|\s*$)", raw)
    organism = m_os.group(1).strip() if m_os else None

    # Extract description (between entry name and OS=)
    desc = None
    m_desc = re.search(r"\|\S+_\S+\s+(.*?)(?:\s+OS=|\s*$)", raw)
    if m_desc:
        desc = m_desc.group(1).strip()

    # No GN= means no gene. The entry-name mnemonic is NOT a gene symbol
    # (ALBU_HUMAN is gene ALB; a TrEMBL entry name is ACCESSION_SPECIES), so
    # deriving one from it fabricated symbols and, for TrEMBL, turned every
    # accession into its own "gene". Leave it None so heal_fasta can consult the
    # canonical map and downstream code can tell "unknown" from "known".

    return ParsedHeader(
        raw=raw,
        accession=accession,
        gene=gene,
        description=desc,
        organism=organism,
        source_db=source_db,
    )


def _parse_gnomad(raw: str, accession: str) -> ParsedHeader:
    """Parse gnomAD variant format: var-POP|ACC|GENE_Mut|pos [sp|...|... desc GN=GENE]"""
    parts = accession.split("|")
    population = parts[0].replace("var-", "") if parts else None
    parent_acc = parts[1] if len(parts) >= 2 else None
    variant_name = parts[2] if len(parts) >= 3 else None

    # Gene from variant name: ALB_Asp574Gly → ALB
    gene = None
    if variant_name:
        gene_candidate = variant_name.split("_")[0]
        if re.match(r"^[A-Z][A-Z0-9]+$", gene_candidate):
            gene = gene_candidate

    # Also try GN= from description
    m_gn = re.search(r"GN=(\S+)", raw)
    if m_gn:
        gene = m_gn.group(1)  # GN= takes precedence

    # Extract variant info
    variant_info = None
    m_var = re.search(r"\[VARIANT:\s*(.*?)\]", raw)
    if m_var:
        variant_info = m_var.group(1)
    elif variant_name and "_" in variant_name:
        variant_info = "p." + variant_name.split("_", 1)[1]

    # Extract description (remove sp|...|... prefix from description)
    desc_part = raw[len(accession) :].strip()
    desc_clean = re.sub(r"^(?:sp|tr)\|\S+\|\S+_\S+\s+", "", desc_part)
    # Remove GN=, VARIANT tags to get clean description
    desc_clean = re.sub(r"\s*GN=\S+", "", desc_clean)
    desc_clean = re.sub(r"\s*\[VARIANT:.*?\]", "", desc_clean)
    desc_clean = re.sub(r"\s*\[AF_\w+:.*?\]", "", desc_clean)
    desc_clean = re.sub(r"\s*PE=\d+\s*SV=\d+", "", desc_clean)
    desc_clean = re.sub(r"\s*OS=.*?(?=\s+PE=|\s+GN=|\s*$)", "", desc_clean).strip()
    if not desc_clean:
        desc_clean = None

    # Extract organism
    m_os = re.search(r"OS=(.+?)(?:\s+OX=|\s+PE=|\s+GN=|\s*$)", raw)
    organism = m_os.group(1).strip() if m_os else None

    return ParsedHeader(
        raw=raw,
        accession=accession,
        gene=gene,
        description=desc_clean,
        organism=organism,
        source_db="gnomad",
        variant_info=variant_info,
        population=population,
        extra={"parent_accession": parent_acc},
    )


def _parse_gencode(raw: str, accession: str) -> ParsedHeader:
    """Parse GENCODE format: ENSP|ENST|ENSG|OTTHUMG|OTTHUMT|transcript-NNN|GENE|length"""
    fields = accession.split("|")
    gene = fields[6] if len(fields) >= 7 else None
    transcript = fields[5] if len(fields) >= 6 else None

    if gene in ("-", "", None):
        gene = None

    return ParsedHeader(
        raw=raw,
        accession=accession,
        gene=gene,
        description=f"isoform {transcript}" if transcript else None,
        source_db="gencode",
        extra={
            "ensp": fields[0] if fields else None,
            "enst": fields[1] if len(fields) >= 2 else None,
            "ensg": fields[2] if len(fields) >= 3 else None,
            "transcript_name": transcript,
        },
    )


def _parse_ens_variant(raw: str, accession: str) -> ParsedHeader:
    """Parse personal variant format: ens|ENST_Mut|GENE description"""
    fields = accession.split("|")
    gene = None
    variant_info = None

    if len(fields) >= 3:
        gene = fields[2].split()[0] if fields[2] else None
        # Variant info from field 1: ENST00000393469.8_Arg664Gln
        if "_" in fields[1]:
            variant_info = "p." + fields[1].split("_", 1)[1]

    # Also try GN= from description
    m_gn = re.search(r"GN=(\S+)", raw)
    if m_gn:
        gene = m_gn.group(1)

    desc = raw[len(accession) :].strip() or None

    return ParsedHeader(
        raw=raw,
        accession=accession,
        gene=gene,
        description=desc,
        source_db="personal_variant",
        variant_info=variant_info,
    )


def _parse_ncbi_microbe(raw: str, accession: str, fields: list[str]) -> ParsedHeader:
    """Parse NCBI microbe format: SPECIES|WP_ACC|LOCUS [gene=X] description OS=..."""
    species_code = fields[0] if fields else None
    refseq_acc = fields[1] if len(fields) >= 2 else None
    locus_tag = fields[2] if len(fields) >= 3 else None

    # Extract gene from gene= tag in description
    desc = raw[len(accession) :].strip()
    gene = None
    m_gene = re.search(r"gene=(\S+)", desc)
    if m_gene:
        gene = m_gene.group(1)
    elif locus_tag:
        gene = locus_tag  # fall back to locus tag

    # Extract organism
    m_os = re.search(r"OS=(.+?)(?:\s*$)", desc)
    organism = m_os.group(1).strip() if m_os else None

    # Clean description (remove gene= and OS= parts)
    desc_clean = re.sub(r"gene=\S+\s*", "", desc)
    desc_clean = re.sub(r"OS=.*$", "", desc_clean).strip()

    return ParsedHeader(
        raw=raw,
        accession=accession,
        gene=gene,
        description=desc_clean or None,
        organism=organism,
        source_db="ncbi_microbe",
        locus_tag=locus_tag,
        extra={"species_code": species_code, "refseq_accession": refseq_acc},
    )


def _parse_mgyg(raw: str, accession: str) -> ParsedHeader:
    """Parse UHGP/MGYG format: MGYG000000000_00000 description"""
    desc = raw[len(accession) :].strip() or None

    # For MGYG, the accession IS the gene identifier (no separate gene name)
    # The genome ID is the part before the underscore
    genome_id = accession.split("_")[0] if "_" in accession else accession

    return ParsedHeader(
        raw=raw,
        accession=accession,
        gene=accession,
        description=desc,
        source_db="uhgp",
        extra={"genome_id": genome_id},
    )


def _parse_simple_microbe(raw: str, accession: str, source_db: str) -> ParsedHeader:
    """Parse simple microbe formats: GMGC, SMORF, Assembly."""
    desc = raw[len(accession) :].strip() or None
    return ParsedHeader(
        raw=raw,
        accession=accession,
        gene=accession,
        description=desc,
        source_db=source_db,
    )


def _parse_ensembl_pep(raw: str, accession: str) -> ParsedHeader:
    """Parse Ensembl pep format: ENSP... pep chromosome:... gene:ENSG... gene_symbol:XXX
    description:...
    """
    gene = None

    # Extract gene_symbol: tag
    m_gs = re.search(r"gene_symbol:(\S+)", raw)
    if m_gs:
        gene = m_gs.group(1)

    # Fall back to gene: tag (ENSG ID)
    if not gene:
        m_gene = re.search(r"gene:(\S+)", raw)
        if m_gene:
            gene = m_gene.group(1)

    # Extract description: tag
    desc = None
    m_desc = re.search(r"description:(.+?)(?:\s+\[Source:|\s*$)", raw)
    if m_desc:
        desc = m_desc.group(1).strip()

    return ParsedHeader(
        raw=raw,
        accession=accession,
        gene=gene,
        description=desc,
        source_db="ensembl_pep",
    )


def _parse_prodigal(raw: str, accession: str) -> ParsedHeader:
    """Parse Prodigal format: 1_1 # start # end # strand # ID=1_1;partial=..."""
    # Prodigal accessions are just numbers (1_1, 1_2, etc.)
    # No gene name — use the accession as gene
    # Extract any useful metadata
    desc = raw[len(accession) :].strip() or None

    return ParsedHeader(
        raw=raw,
        accession=accession,
        gene=accession,
        description=desc,
        source_db="prodigal",
    )


def _parse_sample_ref(raw: str, accession: str) -> ParsedHeader:
    """Parse SAMPLE|REF| prefix format from Stage 2 output.

    Format: SAMPLE|REF|REAL_ACCESSION description
    The real accession is after the second pipe.
    """
    # Strip the SAMPLE|REF| prefix to get the real accession
    parts = accession.split("|")
    if len(parts) >= 3:
        real_accession = "|".join(parts[2:])
        # Re-parse the real accession part
        real_header = real_accession + raw[len(accession) :]
        parsed = parse_header(real_header)
        # Override accession to keep the full SAMPLE|REF| prefix for traceability
        parsed.accession = accession
        parsed.raw = raw
        parsed.extra["sample_ref_prefix"] = "|".join(parts[:2])
        return parsed

    desc = raw[len(accession) :].strip() or None
    return ParsedHeader(
        raw=raw,
        accession=accession,
        gene=accession,
        description=desc,
        source_db="sample_ref",
    )


def _parse_lake_hybrid(raw: str, accession: str) -> ParsedHeader:
    """Parse lake_builder output: PROJECT [SWISSPROT] sp|ACC|NAME desc GN=GENE md5=...
    n_sources=...
    """
    # Find the embedded sp| or tr| header
    m_sp = re.search(r"((?:sp|tr)\|\S+\|\S+)\s+(.*?)(?:\s+md5=|\s*$)", raw)
    if m_sp:
        # Re-parse the embedded UniProt header
        embedded = m_sp.group(1) + " " + m_sp.group(2)
        parsed = parse_header(embedded)
        parsed.accession = accession  # keep full lake accession
        parsed.raw = raw
        parsed.source_db = "lake_hybrid"
        return parsed

    # If no embedded UniProt found, try to extract gene from GN=
    m_gn = re.search(r"GN=(\S+)", raw)
    gene = m_gn.group(1) if m_gn else accession
    desc = raw[len(accession) :].strip() or None

    return ParsedHeader(
        raw=raw,
        accession=accession,
        gene=gene,
        description=desc,
        source_db="lake_hybrid",
    )


def _parse_disease_variant(raw: str, accession: str) -> ParsedHeader:
    """Parse disease variant format: sp|P02649_v21K|APOE_VAR description."""
    parts = accession.split("|")
    gene = None
    variant_info = None

    if len(parts) >= 3:
        # Entry name like APOE_VAR → gene = APOE
        entry_name = parts[2]
        if "_VAR" in entry_name:
            gene = entry_name.replace("_VAR", "")
        elif "_" in entry_name:
            gene = entry_name.split("_")[0]

        # Variant info from accession: P02649_v21K → variant = v21K
        if "_v" in parts[1]:
            variant_info = parts[1].split("_v", 1)[1]

    # Also try GN= from description
    m_gn = re.search(r"GN=(\S+)", raw)
    if m_gn:
        gene = m_gn.group(1)

    desc = raw[len(accession) :].strip() or None

    return ParsedHeader(
        raw=raw,
        accession=accession,
        gene=gene,
        description=desc,
        source_db="disease_variant",
        variant_info=variant_info,
    )


# ═══════════════════════════════════════════════════════════════════════════════
# FASTA HEALING: Process entire files
# ═══════════════════════════════════════════════════════════════════════════════


@dataclass
class HealStats:
    """Statistics from healing a FASTA file."""

    total: int = 0
    kept: int = 0
    skipped: int = 0
    healed: int = 0
    already_ok: int = 0
    by_source_db: dict = field(default_factory=lambda: defaultdict(int))
    by_skip_reason: dict = field(default_factory=lambda: defaultdict(int))
    genes_found: int = 0
    genes_missing: int = 0


def heal_fasta(
    input_fasta: str | Path,
    output_fasta: str | Path,
    target_engine: str = "diann",
    skip_genscan: bool = True,
    skip_no_gene: bool = False,
    canonical_fasta: str | Path | None = None,
) -> HealStats:
    """Heal a FASTA file for a specific search engine.

    Reads every header, parses it, reformats for the target engine,
    and writes the result. Optionally skips GENSCAN predictions, and, only
    when asked (``skip_no_gene=True``), entries whose gene could not be
    resolved. The default keeps them: a header-healing step must not change
    the search space, and an unresolved UniProt entry is written without
    GN= and counted in ``HealStats.genes_missing`` rather than dropped.

    Parameters
    ----------
    input_fasta : path
        Input FASTA file (evidence FASTA from fasta_extractor).
    output_fasta : path
        Output healed FASTA file.
    target_engine : str
        Target search engine: "diann", "sage", or "msfragger".
    skip_genscan : bool
        If True, remove GENSCAN predictions.
    skip_no_gene : bool
        If True, remove entries where no gene could be determined.
    canonical_fasta : path, optional
        SwissProt canonical FASTA for gene name mapping of tr| entries.

    Returns
    -------
    HealStats
        Statistics about the healing process.
    """
    stats = HealStats()

    # Build gene mapping from canonical if provided
    gene_map = _build_gene_map(canonical_fasta) if canonical_fasta else {}

    with open(input_fasta) as fin, open(output_fasta, "w") as fout:
        current_skip = False

        for line in fin:
            if line.startswith(">"):
                stats.total += 1
                raw = line[1:].rstrip()
                parsed = parse_header(raw)
                stats.by_source_db[parsed.source_db] += 1

                # Try to resolve missing gene from canonical map
                if parsed.gene is None and gene_map:
                    parsed.gene = _resolve_gene(parsed, gene_map)

                # Skip decisions
                should_skip = False
                if skip_genscan and parsed.source_db == "genscan":
                    should_skip = True
                    stats.by_skip_reason["genscan"] += 1
                elif skip_no_gene and parsed.gene is None:
                    should_skip = True
                    stats.by_skip_reason[f"{parsed.source_db}_no_gene"] += 1

                if should_skip:
                    stats.skipped += 1
                    current_skip = True
                    continue

                # Format for target engine
                if target_engine == "diann":
                    formatted = parsed.format_for_diann()
                elif target_engine == "sage":
                    formatted = parsed.format_for_sage()
                elif target_engine == "msfragger":
                    formatted = parsed.format_for_msfragger()
                else:
                    formatted = raw

                if formatted != raw:
                    stats.healed += 1
                else:
                    stats.already_ok += 1

                if parsed.gene:
                    stats.genes_found += 1
                else:
                    stats.genes_missing += 1

                stats.kept += 1
                fout.write(f">{formatted}\n")
                current_skip = False
            else:
                if not current_skip:
                    fout.write(line)

    return stats


# ═══════════════════════════════════════════════════════════════════════════════
# ENGINE-AWARE FASTA PREPARATION
# ═══════════════════════════════════════════════════════════════════════════════


@dataclass
class PrepareStats:
    """Statistics from preparing a FASTA for a search engine."""

    target_proteins: int = 0
    decoy_proteins: int = 0
    total_proteins: int = 0
    engine: str = ""
    healed: bool = False
    heal_stats: HealStats | None = None


def generate_decoys(
    input_fasta: str | Path,
    output_fasta: str | Path,
    prefix: str = "rev_",
) -> int:
    """Generate reversed decoy sequences and append to FASTA.

    For each protein in the input, creates a reversed version with the
    given prefix prepended to the accession. The sequence is reversed
    but the header metadata is preserved.

    This is the standard target-decoy approach used by MSFragger,
    Philosopher, Comet, and other tools that don't generate internal
    decoys (unlike SAGE which does this automatically).

    Parameters
    ----------
    input_fasta : path
        Input FASTA file (target sequences).
    output_fasta : path
        Output FASTA with targets + reversed decoys.
    prefix : str
        Decoy prefix. Default: "rev_" (MSFragger/Philosopher standard).

    Returns
    -------
    int
        Number of decoy sequences generated.
    """
    entries: list[tuple[str, str]] = []  # (header, sequence)
    current_header = ""
    current_seq: list[str] = []

    with open(input_fasta) as f:
        for line in f:
            if line.startswith(">"):
                if current_header:
                    entries.append((current_header, "".join(current_seq)))
                current_header = line.rstrip()
                current_seq = []
            else:
                current_seq.append(line.rstrip())
        if current_header:
            entries.append((current_header, "".join(current_seq)))

    n_decoys = 0
    with open(output_fasta, "w") as fout:
        # Write targets
        for header, seq in entries:
            fout.write(f"{header}\n{seq}\n")

        # Write reversed decoys
        for header, seq in entries:
            # >rev_ACCESSION description
            accession = header[1:].split()[0]
            rest = header[1 + len(accession) :]
            fout.write(f">{prefix}{accession}{rest}\n")
            fout.write(f"{seq[::-1]}\n")
            n_decoys += 1

    return n_decoys


def prepare_fasta_for_engine(
    input_fasta: str | Path,
    output_fasta: str | Path,
    engine: str,
    canonical_fasta: str | Path | None = None,
) -> PrepareStats:
    """Prepare a FASTA file for a specific search engine.

    Each search engine has specific requirements:

    - **SAGE**: Works with raw headers. Generates internal decoys.
      Only needs GN= tag for gene-level reporting.
    - **DIA-NN**: Ignores GN= for non-sp/tr headers. Gene name
      must be FIRST WORD of description. No decoys needed (DIA-NN
      uses its own FDR model).
    - **MSFragger**: Needs reversed decoy sequences with ``rev_``
      prefix in the FASTA. Philosopher uses these for target-decoy
      FDR estimation.

    Parameters
    ----------
    input_fasta : path
        Input inference FASTA from FastaLake.
    output_fasta : path
        Engine-ready FASTA.
    engine : str
        Target engine: "sage", "diann", or "msfragger".
    canonical_fasta : path, optional
        SwissProt FASTA for gene name resolution (DIA-NN only).

    Returns
    -------
    PrepareStats
        Preparation statistics.
    """
    stats = PrepareStats(engine=engine)

    if engine == "sage":
        # SAGE: heal headers (add GN=), no decoys needed
        heal_stats = heal_fasta(
            input_fasta,
            output_fasta,
            target_engine="sage",
            skip_genscan=False,
            skip_no_gene=False,
        )
        stats.target_proteins = heal_stats.kept
        stats.total_proteins = heal_stats.kept
        stats.healed = True
        stats.heal_stats = heal_stats

    elif engine == "diann":
        # DIA-NN: heal headers (gene as first word), no decoys
        heal_stats = heal_fasta(
            input_fasta,
            output_fasta,
            target_engine="diann",
            canonical_fasta=canonical_fasta,
        )
        stats.target_proteins = heal_stats.kept
        stats.total_proteins = heal_stats.kept
        stats.healed = True
        stats.heal_stats = heal_stats

    elif engine == "msfragger":
        # MSFragger: heal headers + generate reversed decoys
        # Step 1: heal headers into a temp file
        import tempfile

        with tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".fasta",
            delete=False,
            dir=Path(output_fasta).parent,
        ) as tmp:
            tmp_path = tmp.name

        heal_stats = heal_fasta(
            input_fasta,
            tmp_path,
            target_engine="msfragger",
            skip_genscan=False,
            skip_no_gene=False,
        )
        stats.target_proteins = heal_stats.kept
        stats.healed = True
        stats.heal_stats = heal_stats

        # Step 2: generate reversed decoys
        n_decoys = generate_decoys(tmp_path, output_fasta, prefix="rev_")
        stats.decoy_proteins = n_decoys
        stats.total_proteins = stats.target_proteins + n_decoys

        # Clean up temp
        Path(tmp_path).unlink(missing_ok=True)

    else:
        raise ValueError(f"Unknown engine: {engine!r}. Supported: 'sage', 'diann', 'msfragger'")

    return stats


def _build_gene_map(canonical_fasta: str | Path) -> dict:
    """Build accession → gene and description → gene maps from canonical FASTA."""
    acc_to_gene = {}
    desc_to_gene = {}

    with open(canonical_fasta) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            m_gn = re.search(r"GN=(\S+)", line)
            if not m_gn:
                continue
            gene = m_gn.group(1)

            # Map accession
            parts = line[1:].split("|")
            if len(parts) >= 2:
                acc_to_gene[parts[1].strip()] = gene

            # Map description
            m_desc = re.search(r"\|\S+_\S+\s+(.*?)\s+OS=", line)
            if m_desc:
                desc_to_gene[m_desc.group(1).strip().lower()] = gene

    return {"acc": acc_to_gene, "desc": desc_to_gene}


def _resolve_gene(parsed: ParsedHeader, gene_map: dict) -> str | None:
    """Try to resolve a missing gene name using canonical mapping."""
    acc_map = gene_map.get("acc", {})
    desc_map = gene_map.get("desc", {})

    # Try accession mapping (for tr| entries that share accession with canonical)
    if parsed.source_db == "trembl":
        parts = parsed.accession.split("|")
        if len(parts) >= 2 and parts[1] in acc_map:
            return acc_map[parts[1]]

    # Try description matching
    if parsed.description:
        desc_lower = parsed.description.lower().strip()
        if desc_lower in desc_map:
            return desc_map[desc_lower]
        # Try prefix matching for truncated descriptions
        desc_clean = re.sub(r"\s*\(fragment\)", "", desc_lower).strip()
        for can_desc, gene in desc_map.items():
            if can_desc.startswith(desc_clean) and len(desc_clean) >= 8:
                return gene

    return None


# ═══════════════════════════════════════════════════════════════════════════════
# QC: Validate and report on FASTA headers
# ═══════════════════════════════════════════════════════════════════════════════


def audit_fasta(fasta_path: str | Path) -> dict:
    """Audit a FASTA file's headers without modifying it.

    Returns a dictionary with format distribution, GN= coverage,
    gene counts, and any detected issues.
    """
    stats = {
        "total": 0,
        "by_source_db": defaultdict(int),
        "has_gn": 0,
        "missing_gn": 0,
        "unique_genes": set(),
        "unique_accessions": set(),
        "issues": [],
    }

    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            stats["total"] += 1
            raw = line[1:].rstrip()
            parsed = parse_header(raw)
            stats["by_source_db"][parsed.source_db] += 1
            stats["unique_accessions"].add(parsed.accession)

            if parsed.gene:
                stats["unique_genes"].add(parsed.gene)
                stats["has_gn"] += 1
            else:
                stats["missing_gn"] += 1

    # Convert sets to counts for serialization
    stats["unique_genes"] = len(stats["unique_genes"])
    stats["unique_accessions"] = len(stats["unique_accessions"])
    stats["by_source_db"] = dict(stats["by_source_db"])
    return stats


def print_audit(stats: dict, label: str = "FASTA Audit") -> None:
    """Pretty-print audit results."""
    print(f"\n{'=' * 60}")
    print(f"  {label}")
    print(f"{'=' * 60}")
    print(f"  Total entries:      {stats['total']:,}")
    print(f"  Unique accessions:  {stats['unique_accessions']:,}")
    print(f"  Unique genes:       {stats['unique_genes']:,}")
    print(
        (
            f"  With gene name:     {stats['has_gn']:,} "
            f"({stats['has_gn'] / max(stats['total'], 1) * 100:.1f}%)"
        )
    )
    print(f"  Missing gene name:  {stats['missing_gn']:,}")
    print("\n  Format breakdown:")
    for db, count in sorted(stats["by_source_db"].items(), key=lambda x: -x[1]):
        pct = count / max(stats["total"], 1) * 100
        print(f"    {db:<20} {count:>8,}  ({pct:.1f}%)")
