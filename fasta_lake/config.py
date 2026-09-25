"""
JSON configuration system for FASTA Lake v3 (SAGE-inspired).

Supports 5 use cases from a single config format:

1. **Metaproteomics** (UHGG + GMGC + assemblies, 233M proteins)
2. **Clinical proteomics** (SwissProt + TrEMBL + gnomAD variants)
3. **Strain-level resolution** (variant-specific peptides from WGS)
4. **Environmental** (GMSC + custom assemblies, unknown organisms)
5. **Single-species DDA** (UniProt organism, standard proteomics)

Only 3 fields required: ``database.sources``, ``peptide_evidence.predictions_dir``,
and ``output.directory``. Everything else has sensible defaults.

Design Principles
-----------------
- **No hidden transformations**: every sequence decision is logged
- **I/L normalization is the ONLY sequence modification** (documented, configurable)
- **Provenance tracking**: every protein in output traceable to source DB + evidence
- **Entrapment is orthogonal**: real proteins from known-absent organisms, not decoys
  (SAGE handles target-decoy FDR internally)
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path

from fasta_lake.entrapment_stages import VALID_ENTRAPMENT_STAGES, normalize_stage

logger = logging.getLogger(__name__)

VALID_STRATEGIES = ("species_budget", "razor", "uniform", "uniform_2pep", "info_score", "bayesian")

# Re-exported for callers that historically imported it from this module. Where
# entrapment is injected, earliest first: prelake, prebuild, preselect, appended.
# The first three are database-construction measurements and differ by how much
# curation the entrapment had to survive; appended measures search error.
# Different quantities -- never pool them. Defined in entrapment_stages.py, which
# also documents the run script implementing each.
__all__ = ["VALID_ENTRAPMENT_STAGES", "normalize_stage"]


# ---------------------------------------------------------------------------
# Database configuration
# ---------------------------------------------------------------------------


@dataclass
class DatabaseSource:
    """
    A single protein database source.

    Parameters
    ----------
    tag : str
        Short identifier for this source (e.g., ``'SwissProt'``, ``'UHGG'``).
        Used in output headers and provenance tracking.
    path : str
        Path to FASTA file (plain or gzipped).
    priority : int
        Priority for deduplication (lower = higher priority). When the same
        sequence appears in multiple sources, the highest-priority source
        header is kept.
    type : str
        Source type for inference logic. One of:
        ``'reference'``, ``'variant'``, ``'assembly'``, ``'contaminant'``,
        ``'entrapment'``.
    """

    tag: str
    path: str
    priority: int = 1
    type: str = "reference"

    def __post_init__(self) -> None:
        """Normalize the source path and reject unsupported source categories."""
        self.path = str(Path(self.path))
        valid_types = ("reference", "variant", "assembly", "contaminant", "entrapment")
        if self.type not in valid_types:
            raise ValueError(
                f"Invalid source type '{self.type}' for '{self.tag}'. "
                f"Must be one of: {', '.join(valid_types)}"
            )


@dataclass
class DatabaseConfig:
    """
    Database configuration section.

    Supports multiple source types for different use cases:

    - ``reference``: canonical protein databases (SwissProt, UHGG, GMGC)
    - ``variant``: protein variants (gnomAD SAVs, patient WGS, TrEMBL isoforms)
    - ``assembly``: de novo assemblies (metagenomic, transcriptomic)
    - ``contaminant``: known contaminants (cRAP, keratin, trypsin)
    - ``entrapment``: known-absent proteins for FDR validation
    """

    sources: list[DatabaseSource] = field(default_factory=list)
    output: str = "reference_lake.fasta"

    def __post_init__(self) -> None:
        """Convert source dictionaries into validated DatabaseSource objects."""
        if isinstance(self.sources, list):
            self.sources = [DatabaseSource(**s) if isinstance(s, dict) else s for s in self.sources]


# ---------------------------------------------------------------------------
# Sequence quality filtering
# ---------------------------------------------------------------------------


@dataclass
class SequenceQualityConfig:
    """
    Protein sequence quality filters (all optional).

    These filters remove low-quality sequences BEFORE peptide matching.
    Each filter can be independently enabled/disabled.

    Parameters
    ----------
    min_length : int
        Minimum protein length in amino acids. Proteins shorter than this
        are excluded. Default: 30.
    max_length : int
        Maximum protein length. Proteins longer than this are excluded.
        Default: 50000 (titin is ~34,000 aa).
    reject_poly_runs : bool
        If True, reject proteins with homopolymeric runs (e.g., AAAAAAAAAA).
        These are typically assembly artifacts. Default: True.
    poly_run_length : int
        Minimum consecutive identical amino acids to trigger rejection.
        Default: 10 (e.g., 10 consecutive A's).
    reject_low_complexity : bool
        If True, reject sequences dominated by a small number of amino acids.
        Uses a simple entropy-based filter. Default: False.
    low_complexity_threshold : float
        Shannon entropy threshold (0-1 normalized). Sequences below this
        are rejected. 0.0 = completely uniform, 1.0 = maximum diversity.
        Default: 0.3 (~sequences where >70% is one amino acid).
    reject_fragments : bool
        If True, flag likely protein fragments: sequences that don't start
        with M and are shorter than ``fragment_length_threshold``.
        Default: False.
    fragment_length_threshold : int
        Proteins without M-start shorter than this are flagged as fragments.
        Default: 100.
    """

    min_length: int = 30  # v5.1: raised 6 → 30 to drop fragmented ORFs / pseudogenes
    max_length: int = 50000
    reject_poly_runs: bool = True  # v5.1: on by default (assembly artifacts)
    poly_run_length: int = 10
    reject_low_complexity: bool = False  # keep OFF by default (takes real signal peptides)
    low_complexity_threshold: float = 0.3
    reject_fragments: bool = False  # keep OFF (some truncated ORFs are real)
    fragment_length_threshold: int = 100
    # v5.1: non-canonical residue handling
    strip_noncanonical: bool = True  # replace ambiguous B/J/Z with X; preserve U/O
    reject_if_noncanonical_fraction: float = 0.05  # drop if >5% of residues are non-canonical


# ---------------------------------------------------------------------------
# Peptide evidence filtering
# ---------------------------------------------------------------------------


@dataclass
class PeptideEvidenceConfig:
    """
    Peptide evidence configuration.

    Controls how AlphaNovo predictions are filtered before mapping
    to the protein database.

    Parameters
    ----------
    predictions_dir : str
        Path to directory containing AlphaNovo prediction CSVs.
    min_length : int
        Minimum peptide length (amino acids). Default: 8.
    max_length : int
        Maximum peptide length (amino acids). Default: 50.
    min_score : float
        Minimum AlphaNovo confidence score. Peptides below this are
        discarded. Default: 0.0 (keep all non-zero scores).
    top_percent : float
        Keep only the top N% of peptides per sample by score.
        Set to 1.0 to disable. Default: 0.30 (top 30%).
    normalize_il : bool
        Replace I with L for mass spec equivalence. This is the ONLY
        sequence transformation applied. Default: True.
    """

    predictions_dir: str = ""
    min_length: int = 9  # v5: raised 8 → 9 to match SAGE search default.
    max_length: int = 50
    min_score: float = 0.0
    top_percent: float = 0.30
    normalize_il: bool = True


# ---------------------------------------------------------------------------
# Sample complexity estimation
# ---------------------------------------------------------------------------


@dataclass
class ComplexityConfig:
    """
    Sample complexity estimation from de novo results.

    Maps de novo peptides to the reference database to estimate sample
    complexity BEFORE building the per-sample database. This helps
    choose optimal parameters.

    Parameters
    ----------
    enabled : bool
        Enable complexity estimation. Default: False.
    reference_db : str
        Path to a reference database for complexity estimation.
        If empty, uses the first database source.
    auto_recommend : bool
        If True, recommend optimal strategy and parameters based on
        detected complexity. Default: True.
    report_metrics : list[str]
        Metrics to compute: ``'hit_rate'`` (fraction of peptides matching),
        ``'species_richness'`` (distinct species/organisms),
        ``'peptide_diversity'`` (unique peptides per protein),
        ``'db_coverage'`` (fraction of DB with evidence).
    """

    enabled: bool = False
    reference_db: str = ""
    auto_recommend: bool = True
    report_metrics: list[str] = field(
        default_factory=lambda: ["hit_rate", "species_richness", "peptide_diversity"]
    )


# ---------------------------------------------------------------------------
# Protein inference
# ---------------------------------------------------------------------------


@dataclass
class InferenceConfig:
    """
    Protein inference strategy configuration.

    Parameters
    ----------
    strategy : str
        Inference strategy. One of:

        - ``'species_budget'``: Cross-species aware (metaproteomics)
        - ``'razor'``: MaxQuant-style (classical proteomics)
        - ``'uniform'``: Keep all with evidence (baseline)
        - ``'uniform_2pep'``: Require >=2 peptides (conservative)

    min_peptides : int
        Minimum peptides per protein for uniform_2pep. Default: 1.
    strain_resolution : bool
        If True, treat variant proteins as distinct strains. Variant-specific
        peptides are used to distinguish between canonical and variant forms.
        Relevant for clinical proteomics with gnomAD/WGS data. Default: False.
    keep_contaminants : bool
        If True, always include contaminant-tagged proteins in output
        regardless of peptide evidence. Standard practice for proteomics
        search. Default: True.
    keep_entrapment : bool
        If True, include entrapment-tagged proteins in output for FDR
        validation. Entrapment hits in SAGE results indicate false
        discoveries. Default: True.
    """

    strategy: str = "species_budget"
    min_peptides: int = 1
    strain_resolution: bool = False
    keep_contaminants: bool = True
    keep_entrapment: bool = True

    def __post_init__(self) -> None:
        """Reject an inference strategy absent from the supported strategy registry."""
        if self.strategy not in VALID_STRATEGIES:
            raise ValueError(
                f"Invalid strategy '{self.strategy}'. Must be one of: {', '.join(VALID_STRATEGIES)}"
            )


# ---------------------------------------------------------------------------
# Entrapment
# ---------------------------------------------------------------------------


@dataclass
class EntrapmentConfig:
    """
    Entrapment database configuration for FDR validation.

    Entrapment proteins are real proteins from organisms guaranteed to
    be absent from the sample. Any SAGE hits to entrapment proteins
    indicate false discoveries, providing an orthogonal FDR estimate
    independent of the target-decoy approach.

    Parameters
    ----------
    enabled : bool
        Enable entrapment. Default: False.
    method : str
        How to generate/source entrapment proteins:

        - ``'phylogenetic'``: Use proteins from phylogenetically distant
          organisms (e.g., Arabidopsis in human sample). Most realistic.
        - ``'shuffled'``: Shuffle real protein sequences preserving amino
          acid composition. Computational approach.
        - ``'user'``: User provides their own entrapment FASTA.

    source : str
        Path to entrapment FASTA (for method ``'user'`` or ``'phylogenetic'``).
        For ``'phylogenetic'``, can be a species name (e.g., ``'arabidopsis'``)
        and FASTA Lake will use a built-in set.
    n_proteins : int
        Number of entrapment proteins to include. Default: 1000.
        For ``'user'`` method, includes all proteins from the FASTA.
    expected_fdr : float
        Expected false discovery rate. Used for reporting -- if entrapment
        hit rate exceeds this, a warning is raised. Default: 0.01 (1%).
    entrap_class : str
        Entrapment class from ``fasta_lake.entrapment_classes``: one of
        ``'shuffled'``, ``'thermophilic'``, ``'archaea'``, ``'plant'``,
        ``'non_self'``, ``'user'``. Each class declares its difficulty tier,
        what quantity its FDP estimates, and its caveat, and that metadata is
        carried into the FDP report. Empty string falls back to ``method``.
    stage : str
        Where entrapment is injected. One of, earliest first:

        - ``'prelake'``: into the input lake before the cohort evidence filter.
          Entrapment faces every curation step; the fullest construction test.
        - ``'prebuild'``: onto the clustered lake, then re-clustered at 97%, so
          entrapment must also survive clustering.
        - ``'preselect'`` (default): onto the clustered lake with no
          re-clustering. Entrapment must be nominated by the sample's own de
          novo peptides.
        - ``'appended'``: onto the finished per-sample database after
          parsimony, entering the search untested.

        The first three measure database-construction error and differ by how
        much curation entrapment had to survive; ``appended`` measures ordinary
        search error. These are DIFFERENT QUANTITIES and must never be pooled.
        ``'pre-selection'`` is accepted as an alias for ``'preselect'``.
    seed : int
        Random seed for entrapment generation. Default: 42. Zero is rejected
        (it is a fixed point of the xorshift64 generator in the Rust extractor).
    """

    enabled: bool = False
    method: str = "phylogenetic"
    source: str = ""
    n_proteins: int = 1000
    expected_fdr: float = 0.01
    entrap_class: str = ""
    stage: str = "preselect"
    seed: int = 42

    def __post_init__(self) -> None:
        # Accepts aliases ('pre-selection') and normalises to the canonical name
        # the run scripts and figure loaders use, so a config and a run cannot
        # disagree about which injection point produced a number.
        """Normalize the injection stage and validate seed, count and FDR settings."""
        try:
            self.stage = normalize_stage(self.stage)
        except ValueError as exc:
            raise ValueError(f"entrapment.stage: {exc}") from None
        if self.seed == 0:
            raise ValueError(
                "entrapment.seed must be nonzero (0 is a fixed point of the "
                "xorshift64 generator used by the Rust extractor)"
            )
        if self.n_proteins < 0:
            raise ValueError(f"entrapment.n_proteins must be >= 0; got {self.n_proteins}")
        if not (0.0 < self.expected_fdr <= 1.0):
            raise ValueError(f"entrapment.expected_fdr must be in (0, 1]; got {self.expected_fdr}")


# ---------------------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------------------


@dataclass
class ClusteringConfig:
    """Clustering configuration section.

    The configured ``identity`` is the *single source of truth* for the MMseqs2
    ``--min-seq-id`` value the pipeline runs. Supported sweep values are
    ``0.90, 0.95, 0.97, 0.99, 1.0``. A value of ``1.0`` means exact dedup /
    skip clustering, matching the historical pipeline semantics (clustering at
    100% identity is degenerate, so it is treated as a no-op pass-through of
    the input FASTA).
    """

    enabled: bool = True
    identity: float = 0.97
    coverage: float = 0.8

    def __post_init__(self) -> None:
        """Require identity and coverage fractions in the interval (0, 1]."""
        if not (0.0 < self.identity <= 1.0):
            raise ValueError(f"clustering.identity must be in (0, 1]; got {self.identity}")
        if not (0.0 < self.coverage <= 1.0):
            raise ValueError(f"clustering.coverage must be in (0, 1]; got {self.coverage}")

    @property
    def skip_clustering(self) -> bool:
        """True when clustering is a no-op for this config.

        Clustering is skipped when explicitly disabled, or when ``identity``
        is ``1.0`` (exact-dedup / pass-through), matching current semantics.
        """
        return (not self.enabled) or self.identity >= 1.0

    @property
    def min_seq_id(self) -> float:
        """The value to pass as MMseqs2 ``--min-seq-id``."""
        return float(self.identity)

    def mmseqs_args(self, *, threads: int = 16) -> list[str]:
        """The MMseqs2 ``easy-cluster`` argument list driven by this config.

        Returns the flags appended after ``easy-cluster <in> <out> <tmp>``.
        Default (identity=0.97, coverage=0.8) reproduces the historical
        hardcoded command exactly:
        ``--min-seq-id 0.97 -c 0.8 --cov-mode 0 --threads 16 --remove-tmp-files 1``.
        """
        return [
            "--min-seq-id",
            f"{self.min_seq_id:g}",
            "-c",
            f"{self.coverage:g}",
            "--cov-mode",
            "0",
            "--threads",
            str(threads),
            "--remove-tmp-files",
            "1",
        ]


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


@dataclass
class OutputConfig:
    """
    Output configuration.

    Parameters
    ----------
    directory : str
        Output directory path. Created if it doesn't exist.
    write_stats : bool
        Write JSON stats file alongside output FASTA. Default: True.
    write_hits_tsv : bool
        Write peptide-protein hits TSV. Default: True.
    write_audit_log : bool
        Write a detailed audit log tracking every protein decision
        (included/excluded, reason, source DB, evidence). Default: False.
    provenance_tracking : bool
        Tag output FASTA headers with source database and evidence count.
        Enables full traceability from output back to input. Default: True.
    """

    directory: str = "./output"
    write_stats: bool = True
    write_hits_tsv: bool = True
    write_audit_log: bool = False
    provenance_tracking: bool = True

    def __post_init__(self) -> None:
        """Normalize the configured output path to a string."""
        self.directory = str(Path(self.directory))


# ---------------------------------------------------------------------------
# Rust binary config
# ---------------------------------------------------------------------------


@dataclass
class RustConfig:
    """Rust binary configuration section."""

    fasta_extractor: str | None = None
    lake_builder: str | None = None
    threads: int | None = None


# ---------------------------------------------------------------------------
# Top-level config
# ---------------------------------------------------------------------------

VALID_ACQUISITIONS = ("dda", "dia")
VALID_ENGINES = ("sage", "diann", "msfragger")


@dataclass
class FastaLakeConfig:
    """
    Complete FASTA Lake v3 configuration.

    Supports all proteomics use cases from a single config format.
    Only 3 fields required: ``database.sources``,
    ``peptide_evidence.predictions_dir``, and ``output.directory``.

    Parameters
    ----------
    acquisition : str
        Data acquisition mode. Determines de novo engine expectations and
        search engine compatibility:

        - ``'dda'``: Data-Dependent Acquisition → AlphaNovo → SAGE or MSFragger
        - ``'dia'``: Data-Independent Acquisition → DIANovo → DIA-NN or SAGE (wide-window)

    search_engine : str
        Target search engine. Determines FASTA export format:

        - ``'sage'``: Add GN= tags. SAGE generates internal decoys.
        - ``'diann'``: Gene name as first word of description (critical
          for non-sp/tr entries). No decoys needed.
        - ``'msfragger'``: Add GN= tags + reversed decoy sequences (rev_ prefix).

    database : DatabaseConfig
        Reference database sources (target + variant + contaminant + entrapment).
    sequence_quality : SequenceQualityConfig
        Protein sequence quality filters (optional, all off by default).
    peptide_evidence : PeptideEvidenceConfig
        Peptide evidence filtering parameters.
    complexity : ComplexityConfig
        Sample complexity estimation (optional).
    inference : InferenceConfig
        Protein inference strategy configuration.
    entrapment : EntrapmentConfig
        Entrapment database for FDR validation (optional).
    clustering : ClusteringConfig
        MMseqs2 clustering parameters.
    output : OutputConfig
        Output directory and format options.
    rust : RustConfig
        Rust binary paths and thread configuration.
    """

    acquisition: str = "dda"
    search_engine: str = "sage"
    database: DatabaseConfig = field(default_factory=DatabaseConfig)
    sequence_quality: SequenceQualityConfig = field(default_factory=SequenceQualityConfig)
    peptide_evidence: PeptideEvidenceConfig = field(default_factory=PeptideEvidenceConfig)
    complexity: ComplexityConfig = field(default_factory=ComplexityConfig)
    inference: InferenceConfig = field(default_factory=InferenceConfig)
    entrapment: EntrapmentConfig = field(default_factory=EntrapmentConfig)
    clustering: ClusteringConfig = field(default_factory=ClusteringConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    rust: RustConfig = field(default_factory=RustConfig)

    def __post_init__(self) -> None:
        # Validate acquisition and search_engine
        """Normalize and validate acquisition and search-engine names."""
        self.acquisition = self.acquisition.lower()
        self.search_engine = self.search_engine.lower()
        if self.acquisition not in VALID_ACQUISITIONS:
            raise ValueError(
                f"Invalid acquisition '{self.acquisition}'. "
                f"Must be one of: {', '.join(VALID_ACQUISITIONS)}"
            )
        if self.search_engine not in VALID_ENGINES:
            raise ValueError(
                f"Invalid search_engine '{self.search_engine}'. "
                f"Must be one of: {', '.join(VALID_ENGINES)}"
            )
        if self.acquisition == "dia" and self.search_engine == "msfragger":
            logger.warning(
                "DIA + MSFragger is an unusual combination. "
                "Consider DIA-NN for DIA data (set search_engine='diann')."
            )

        # Convert nested dicts to dataclass instances
        _sections = {
            "database": DatabaseConfig,
            "sequence_quality": SequenceQualityConfig,
            "peptide_evidence": PeptideEvidenceConfig,
            "complexity": ComplexityConfig,
            "inference": InferenceConfig,
            "entrapment": EntrapmentConfig,
            "clustering": ClusteringConfig,
            "output": OutputConfig,
            "rust": RustConfig,
        }
        for attr, cls in _sections.items():
            val = getattr(self, attr)
            if isinstance(val, dict):
                setattr(self, attr, cls(**val))

    def validate(self) -> None:
        """Validate that required fields are set and paths are sensible."""
        if not self.database.sources:
            raise ValueError("At least one database source is required")
        if not self.peptide_evidence.predictions_dir:
            raise ValueError("peptide_evidence.predictions_dir is required")
        if not self.output.directory:
            raise ValueError("output.directory is required")

    def get_sources_by_type(self, source_type: str) -> list[DatabaseSource]:
        """Return database sources filtered by type."""
        return [s for s in self.database.sources if s.type == source_type]

    def has_entrapment(self) -> bool:
        """Check if entrapment is configured (via config or source tags)."""
        return self.entrapment.enabled or any(s.type == "entrapment" for s in self.database.sources)

    def has_variants(self) -> bool:
        """Check if variant databases are configured."""
        return any(s.type == "variant" for s in self.database.sources)

    def has_contaminants(self) -> bool:
        """Check if contaminant databases are configured."""
        return any(s.type == "contaminant" for s in self.database.sources)


# ---------------------------------------------------------------------------
# Config I/O
# ---------------------------------------------------------------------------


def load_config(config_path: str | Path) -> FastaLakeConfig:
    """
    Load a FASTA Lake configuration from a JSON file.

    Parameters
    ----------
    config_path : str or Path
        Path to JSON configuration file.

    Returns
    -------
    FastaLakeConfig
        Parsed and validated configuration.

    Raises
    ------
    FileNotFoundError
        If config file does not exist.
    json.JSONDecodeError
        If config file is not valid JSON.
    """
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    logger.info("Loading config from %s", config_path)
    with open(config_path) as f:
        data = json.load(f)

    return FastaLakeConfig(**data)


def create_example_config(template: str = "full") -> str:
    """
    Generate an example JSON configuration string.

    Parameters
    ----------
    template : str
        Template type: ``'minimal'``, ``'clinical'``, ``'metaproteomics'``,
        or ``'full'``.

    Returns
    -------
    str
        JSON configuration string.
    """
    if template == "minimal":
        config = {
            "acquisition": "dda",
            "search_engine": "sage",
            "database": {"sources": [{"tag": "REF", "path": "proteins.fasta"}]},
            "peptide_evidence": {"predictions_dir": "./predictions"},
            "output": {"directory": "./output"},
        }
    elif template == "clinical":
        config = {
            "acquisition": "dda",
            "search_engine": "sage",
            "database": {
                "sources": [
                    {
                        "tag": "SwissProt",
                        "priority": 1,
                        "path": "/data/uniprot_sprot_human.fasta",
                        "type": "reference",
                    },
                    {
                        "tag": "TrEMBL",
                        "priority": 2,
                        "path": "/data/uniprot_trembl_human.fasta",
                        "type": "reference",
                    },
                    {
                        "tag": "gnomAD",
                        "priority": 3,
                        "path": "/data/gnomad_variants.fasta",
                        "type": "variant",
                    },
                    {
                        "tag": "cRAP",
                        "priority": 99,
                        "path": "/data/cRAP_contaminants.fasta",
                        "type": "contaminant",
                    },
                    {
                        "tag": "Entrap_Arabidopsis",
                        "priority": 100,
                        "path": "/data/arabidopsis_thaliana.fasta",
                        "type": "entrapment",
                    },
                ],
                "output": "clinical_lake.fasta",
            },
            "sequence_quality": {
                "min_length": 6,
                "max_length": 50000,
                "reject_poly_runs": True,
                "poly_run_length": 10,
                "reject_low_complexity": True,
                "low_complexity_threshold": 0.3,
                "reject_fragments": True,
                "fragment_length_threshold": 100,
            },
            "peptide_evidence": {
                "predictions_dir": "/data/alphanovo_results",
                "min_length": 8,
                "max_length": 50,
                "min_score": 0.5,
                "top_percent": 0.50,
                "normalize_il": True,
            },
            "complexity": {
                "enabled": True,
                "auto_recommend": True,
                "report_metrics": [
                    "hit_rate",
                    "species_richness",
                    "peptide_diversity",
                    "db_coverage",
                ],
            },
            "inference": {
                "strategy": "razor",
                "min_peptides": 1,
                "strain_resolution": True,
                "keep_contaminants": True,
                "keep_entrapment": True,
            },
            "entrapment": {
                "enabled": True,
                "method": "phylogenetic",
                "source": "/data/arabidopsis_thaliana.fasta",
                "n_proteins": 1000,
                "expected_fdr": 0.01,
            },
            "output": {
                "directory": "./output",
                "write_stats": True,
                "write_hits_tsv": True,
                "write_audit_log": True,
                "provenance_tracking": True,
            },
        }
    elif template == "metaproteomics":
        config = {
            "acquisition": "dda",
            "search_engine": "sage",
            "database": {
                "sources": [
                    {
                        "tag": "UHGG",
                        "priority": 1,
                        "path": "/data/uhgp-100.faa",
                        "type": "reference",
                    },
                    {
                        "tag": "GMGC",
                        "priority": 2,
                        "path": "/data/GMGC10.100AA.faa",
                        "type": "reference",
                    },
                    {
                        "tag": "GMSC",
                        "priority": 3,
                        "path": "/data/GMSC10.100AA.faa",
                        "type": "reference",
                    },
                    {
                        "tag": "Entrap_Arabidopsis",
                        "priority": 100,
                        "path": "/data/arabidopsis_thaliana.fasta",
                        "type": "entrapment",
                    },
                ],
                "output": "reference_lake.fasta",
            },
            "sequence_quality": {
                "reject_poly_runs": True,
                "reject_low_complexity": True,
            },
            "peptide_evidence": {
                "predictions_dir": "/data/alphanovo_results",
                "min_length": 8,
                "max_length": 50,
                "min_score": 0.0,
                "top_percent": 0.30,
                "normalize_il": True,
            },
            "complexity": {
                "enabled": True,
                "auto_recommend": True,
            },
            "inference": {
                "strategy": "species_budget",
                "min_peptides": 1,
                "strain_resolution": False,
                "keep_entrapment": True,
            },
            "entrapment": {
                "enabled": True,
                "method": "phylogenetic",
                "source": "/data/arabidopsis_thaliana.fasta",
                "n_proteins": 1000,
                "expected_fdr": 0.01,
            },
            "clustering": {
                "enabled": True,
                "identity": 0.97,
                "coverage": 0.8,
            },
            "output": {
                "directory": "./output",
                "write_stats": True,
                "write_hits_tsv": True,
                "write_audit_log": True,
                "provenance_tracking": True,
            },
        }
    else:  # "full"
        config = {
            "database": {
                "sources": [
                    {"tag": "UHGG", "priority": 1, "path": "/data/uhgp-100.faa"},
                    {"tag": "GMGC", "priority": 2, "path": "/data/GMGC10.100AA.faa"},
                ],
                "output": "reference_lake.fasta",
            },
            "sequence_quality": {
                "min_length": 6,
                "max_length": 50000,
                "reject_poly_runs": False,
                "poly_run_length": 10,
                "reject_low_complexity": False,
                "low_complexity_threshold": 0.3,
                "reject_fragments": False,
                "fragment_length_threshold": 100,
            },
            "peptide_evidence": {
                "predictions_dir": "/data/alphanovo_results",
                "min_length": 8,
                "max_length": 50,
                "min_score": 0.0,
                "top_percent": 0.30,
                "normalize_il": True,
            },
            "complexity": {
                "enabled": False,
                "auto_recommend": True,
                "report_metrics": ["hit_rate", "species_richness", "peptide_diversity"],
            },
            "inference": {
                "strategy": "species_budget",
                "min_peptides": 1,
                "strain_resolution": False,
                "keep_contaminants": True,
                "keep_entrapment": True,
            },
            "entrapment": {
                "enabled": False,
                "method": "phylogenetic",
                "source": "",
                "n_proteins": 1000,
                "expected_fdr": 0.01,
            },
            "clustering": {
                "enabled": True,
                "identity": 0.97,
                "coverage": 0.8,
            },
            "output": {
                "directory": "./output",
                "write_stats": True,
                "write_hits_tsv": True,
                "write_audit_log": False,
                "provenance_tracking": True,
            },
            "rust": {
                "fasta_extractor": None,
                "lake_builder": None,
                "threads": None,
            },
        }
    return json.dumps(config, indent=2)
