"""
FastaLake run configuration — the single file that captures everything.

Every pipeline run produces a run_config.yaml that records:
1. All input paths (predictions, databases, mzML)
2. All parameters (score thresholds, clustering identity, inference strategy)
3. All output paths and metrics (entry counts, hit rates, gene counts)
4. Benchmark results (vs canonical or matched metagenome)
5. Timestamps and software versions

This enables:
- Exact reproduction of any run
- Parameter comparison between runs
- Provenance tracking from result → config → raw data
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

from fasta_lake.helpers import count_fasta_headers

logger = logging.getLogger(__name__)


@dataclass
class RunConfig:
    """Complete configuration + results for a FastaLake pipeline run."""

    # ── Project metadata ──
    project_name: str = ""
    description: str = ""
    date: str = field(default_factory=lambda: time.strftime("%Y-%m-%d"))
    fastalake_version: str = "4.0.0-dev"
    fastalake_git_commit: str = ""  # git rev-parse HEAD
    rust_binary_path: str = ""  # path to fasta_extractor binary
    rust_binary_date: str = ""  # modification date of binary
    python_version: str = ""
    hostname: str = ""

    # ── Lake provenance ──
    lake_build_date: str = ""  # when the lake was built
    lake_build_config: str = ""  # path to the lake_builder config/command used
    lake_md5: str = ""  # md5sum of the lake FASTA (for exact version tracking)

    # ── Input ──
    predictions_dir: str = ""
    denovo_engine: str = "alphanovo"  # alphanovo, dianovo, casanovo
    n_samples: int = 0
    acquisition: str = "dda"  # dda, dia
    host_organism: str = "human"
    experiment_type: str = "proteomics"  # proteomics, metaproteomics

    # ── Database ──
    lake_path: str = ""
    lake_entries: int = 0
    canonical_path: str = ""
    canonical_entries: int = 0
    sources: list[dict] = field(default_factory=list)
    # Each source: {"tag": "UHGG", "path": "...", "entries": 95303635, "priority": 1}

    # ── Pipeline parameters ──
    # DENOVO_EXTRACT — de novo peptide filter (stages 1-3)
    top_percent: float = 0.30
    denovo_min_length: int = 9  # v5: admits proteins into evidence lake. Stages 1-3.
    denovo_max_length: int = 50
    normalize_il: bool = True

    # CLUSTER
    cluster_enabled: bool = True
    cluster_identity: float = 0.97
    cluster_coverage: float = 0.8

    # SAMPLE_EXTRACT
    include_canonical: bool = True

    # INFER
    infer_strategy: str = "species_budget"
    infer_min_peptides: int = 1

    # KMER_RESCUE
    kmer_rescue: bool = False
    kmer_k: int = 5

    # FORMAT_HEADERS
    header_engine: str = "sage"  # sage, diann

    # SEARCH — tryptic digest during search (stage 4)
    search_engine: str = "sage"  # sage, diann
    search_min_length: int = 9  # v5: SAGE enzyme.min_len (theoretical digest floor).
    search_max_length: int = 50
    search_missed_cleavages: int = 2
    search_max_variable_mods: int = 2
    search_precursor_ppm: float = 10.0
    search_fdr: float = 0.01  # q-value threshold for acceptance

    # ── Back-compat aliases ──
    # min_length / max_length used to be single fields (stages 1-4). Keep as
    # properties that mirror denovo_min_length for code that still reads them.
    @property
    def min_length(self) -> int:
        """Read the legacy length alias, or set both de novo and search minima."""
        return self.denovo_min_length

    @min_length.setter
    def min_length(self, v: int) -> None:
        # When callers set .min_length, propagate to both stages (legacy behaviour).
        """Read the legacy length alias, or set both de novo and search minima."""
        self.denovo_min_length = v
        self.search_min_length = v

    @property
    def max_length(self) -> int:
        """Read the legacy length alias, or set both de novo and search maxima."""
        return self.denovo_max_length

    @max_length.setter
    def max_length(self, v: int) -> None:
        """Read the legacy length alias, or set both de novo and search maxima."""
        self.denovo_max_length = v
        self.search_max_length = v

    # ── Results (filled after run) ──
    evidence_entries: int | None = None
    evidence_hit_rate: float | None = None
    clustered_entries: int | None = None
    per_sample_entries: dict[str, int] = field(default_factory=dict)
    inferred_entries: dict[str, int] = field(default_factory=dict)

    # Search results per sample
    search_peptides: dict[str, int] = field(default_factory=dict)
    search_genes: dict[str, int] = field(default_factory=dict)
    search_protein_groups: dict[str, int] = field(default_factory=dict)

    # Benchmark
    benchmark_baseline: str = ""  # path to baseline results
    benchmark_peptide_gain_pct: float | None = None
    benchmark_gene_coverage_pct: float | None = None
    benchmark_entrapment_fdr: float | None = None
    benchmark_intensity_r: float | None = None

    # ── Timestamps ──
    started: str = ""
    completed: str = ""

    def capture_environment(self) -> None:
        """Auto-populate environment fields (git commit, binary paths, etc.)."""
        import platform
        import subprocess

        self.hostname = platform.node()
        self.python_version = platform.python_version()
        self.date = time.strftime("%Y-%m-%d %H:%M:%S")

        # Git commit
        try:
            pkg_dir = Path(__file__).parent.parent
            result = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                cwd=str(pkg_dir),
                capture_output=True,
                text=True,
                timeout=5,
            )
            if result.returncode == 0:
                self.fastalake_git_commit = result.stdout.strip()[:12]
        except Exception:
            pass

        # Rust binary
        candidates = [
            Path(__file__).parent.parent
            / "rust"
            / "fasta_extractor_v2"
            / "target"
            / "release"
            / "fasta_extractor",
        ]
        for c in candidates:
            if c.exists():
                self.rust_binary_path = str(c)
                self.rust_binary_date = time.strftime("%Y-%m-%d", time.localtime(c.stat().st_mtime))
                break

        # Lake MD5: full-file identity, including records beyond the first MiB.
        if self.lake_path and Path(self.lake_path).exists():
            try:
                import hashlib

                h = hashlib.md5(usedforsecurity=False)  # file identity only, not a security hash
                with open(self.lake_path, "rb") as f:
                    for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
                        h.update(block)
                self.lake_md5 = h.hexdigest()
            except Exception:
                pass

    def save(self, path: str | Path) -> None:
        """Save config to YAML."""
        path = Path(path)
        # Use JSON since we don't want to add pyyaml as dependency
        # YAML is a superset of JSON so this is fine
        with open(path, "w") as f:
            json.dump(asdict(self), f, indent=2, default=str)
        logger.info("Config saved: %s", path)

    @classmethod
    def load(cls, path: str | Path) -> RunConfig:
        """Load config from JSON/YAML."""
        path = Path(path)
        with open(path) as f:
            data = json.load(f)
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})

    def update_evidence(self, evidence_fasta: str | Path) -> None:
        """Update config with evidence lake results."""
        n = count_fasta_headers(evidence_fasta)
        self.evidence_entries = n
        if self.lake_entries > 0:
            self.evidence_hit_rate = n / self.lake_entries * 100

    def update_cluster(self, clustered_fasta: str | Path) -> None:
        """Update config with clustering results."""
        self.clustered_entries = count_fasta_headers(clustered_fasta)

    def update_sample(self, sample_id: str, fasta_path: str | Path) -> None:
        """Update config with per-sample extraction results."""
        n = count_fasta_headers(fasta_path)
        self.per_sample_entries[sample_id] = n

    def update_search(self, sample_id: str, peptides: int, genes: int, pgs: int) -> None:
        """Update config with search results for one sample."""
        self.search_peptides[sample_id] = peptides
        self.search_genes[sample_id] = genes
        self.search_protein_groups[sample_id] = pgs

    def summary(self) -> str:
        """Human-readable summary."""
        lines = [
            f"{'=' * 60}",
            f"  FastaLake Run: {self.project_name}",
            f"  {self.description}",
            f"  Date: {self.date}",
            f"{'=' * 60}",
            "",
            "  Input:",
            f"    Predictions: {self.predictions_dir}",
            f"    Engine: {self.denovo_engine}, Acquisition: {self.acquisition}",
            f"    Samples: {self.n_samples}",
            "",
            "  Database:",
            f"    Lake: {self.lake_entries:,} proteins",
            f"    Canonical: {self.canonical_entries:,} proteins",
            f"    Sources: {len(self.sources)}",
            "",
            "  Parameters:",
            f"    Score filter: top {self.top_percent * 100:.0f}%",
            f"    Clustering: {'on' if self.cluster_enabled else 'off'}"
            f" ({self.cluster_identity * 100:.0f}%)"
            if self.cluster_enabled
            else "",
            f"    Inference: {self.infer_strategy}",
            f"    K-mer rescue: {'on' if self.kmer_rescue else 'off'}",
            f"    Search: {self.search_engine} ({self.header_engine} headers)",
        ]

        if self.evidence_entries is not None:
            lines.extend(
                [
                    "",
                    "  Results:",
                    (
                        f"    Evidence lake: {self.evidence_entries:,} "
                        f"({self.evidence_hit_rate:.2f}% hit rate)"
                    ),
                ]
            )
            if self.clustered_entries is not None:
                lines.append(f"    After clustering: {self.clustered_entries:,}")
            if self.per_sample_entries:
                mean_ps = sum(self.per_sample_entries.values()) // len(self.per_sample_entries)
                lines.append(f"    Per-sample mean: {mean_ps:,}")
            if self.search_peptides:
                mean_pep = sum(self.search_peptides.values()) // len(self.search_peptides)
                mean_gene = sum(self.search_genes.values()) // len(self.search_genes)
                lines.append(f"    Search mean: {mean_pep:,} peptides, {mean_gene:,} genes")

        if self.benchmark_peptide_gain_pct is not None:
            lines.extend(
                [
                    "",
                    f"  Benchmark vs {Path(self.benchmark_baseline).name}:",
                    f"    Peptide gain: {self.benchmark_peptide_gain_pct:+.1f}%",
                    f"    Gene coverage: {self.benchmark_gene_coverage_pct:.1f}%",
                    f"    Entrapment FDR: {self.benchmark_entrapment_fdr:.2f}%",
                    f"    Intensity r: {self.benchmark_intensity_r:.4f}",
                ]
            )

        lines.append(f"{'=' * 60}")
        return "\n".join(lines)
