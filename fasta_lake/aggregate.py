"""
Cross-sample aggregation for FastaLake search results.

Wraps the Rust aggregation_engine for dual-index (peptide + protein group)
aggregation, and provides Python utilities for cross-sample comparison,
quantification matrices, and cohort-level statistics.

The Rust aggregation_engine handles:
- Peptide-level aggregation with protein group normalization
- Source tag stripping (REF|, SAMPLE|REF|)
- Protein group classification (human, microbial, mixed)
- Dual-index output (both peptide and protein group level)

This Python module adds:
- Multi-sample loading and merging
- Quantification matrix construction
- Cross-sample statistics (CV, missing values, overlap)
- Comparison with canonical baseline (via benchmark.py)
"""

from __future__ import annotations

import csv
import logging
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from fasta_lake.helpers import require_columns

logger = logging.getLogger(__name__)

DIANN_COLUMNS = ("Q.Value", "Stripped.Sequence", "Genes", "Protein.Group", "Precursor.Quantity")


@dataclass
class SampleResult:
    """Parsed search result for one sample."""

    sample_id: str
    peptides: dict[str, float] = field(default_factory=dict)  # peptide → intensity
    genes: dict[str, float] = field(default_factory=dict)  # gene → summed intensity
    protein_groups: dict[str, float] = field(default_factory=dict)  # pg → intensity
    n_precursors: int = 0


def load_diann_sample(
    report_path: str | Path,
    sample_id: str | None = None,
    qvalue: float = 0.01,
) -> SampleResult:
    """Load a single-sample DIA-NN report."""
    report_path = Path(report_path)
    if sample_id is None:
        sample_id = report_path.parent.name

    result = SampleResult(sample_id=sample_id)

    with open(report_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        require_columns(reader.fieldnames, DIANN_COLUMNS, report_path)
        for row in reader:
            q = float(row["Q.Value"])
            if q > qvalue:
                continue

            result.n_precursors += 1
            pep = row["Stripped.Sequence"]
            gene = row["Genes"]
            pg = row["Protein.Group"]
            intensity = float(row["Precursor.Quantity"] or 0)

            if pep and intensity > 0:
                if pep not in result.peptides or intensity > result.peptides[pep]:
                    result.peptides[pep] = intensity

            if gene and intensity > 0:
                result.genes[gene] = result.genes.get(gene, 0) + intensity

            if pg and intensity > 0:
                if pg not in result.protein_groups or intensity > result.protein_groups[pg]:
                    result.protein_groups[pg] = intensity

    logger.info(
        "  %s: %d peptides, %d genes, %d precursors",
        sample_id,
        len(result.peptides),
        len(result.genes),
        result.n_precursors,
    )
    return result


def build_quant_matrix(
    samples: list[SampleResult],
    level: str = "gene",
) -> tuple[list[str], list[str], list[list[float]]]:
    """Build a quantification matrix across samples.

    Parameters
    ----------
    samples : list of SampleResult
    level : 'gene', 'peptide', or 'protein_group'

    Returns
    -------
    (sample_ids, feature_ids, matrix)
        ``matrix[i][j]`` is the intensity of feature j in sample i; 0 means missing.
    """
    # Collect all features
    all_features: set[str] = set()
    for s in samples:
        if level == "gene":
            all_features.update(s.genes.keys())
        elif level == "peptide":
            all_features.update(s.peptides.keys())
        elif level == "protein_group":
            all_features.update(s.protein_groups.keys())

    feature_ids = sorted(all_features)
    sample_ids = [s.sample_id for s in samples]

    matrix = []
    for s in samples:
        if level == "gene":
            data = s.genes
        elif level == "peptide":
            data = s.peptides
        else:
            data = s.protein_groups

        row = [data.get(f, 0.0) for f in feature_ids]
        matrix.append(row)

    return sample_ids, feature_ids, matrix


def compute_cross_sample_stats(
    samples: list[SampleResult],
    level: str = "gene",
) -> dict:
    """Compute cross-sample statistics.

    Returns dict with:
    - n_samples, n_features
    - completeness: % of matrix cells that are non-zero
    - per_feature_cv: coefficient of variation for each feature (across samples where present)
    - median_cv: median CV across all features
    - features_in_all: features found in every sample
    - features_in_one: features found in only one sample
    """
    sample_ids, feature_ids, matrix = build_quant_matrix(samples, level)

    n_samples = len(samples)
    n_features = len(feature_ids)
    total_cells = n_samples * n_features

    non_zero = sum(1 for row in matrix for v in row if v > 0)
    completeness = non_zero / max(total_cells, 1) * 100

    # Per-feature stats
    feature_cvs = []
    features_in_all = 0
    features_in_one = 0

    for j in range(n_features):
        values = [matrix[i][j] for i in range(n_samples) if matrix[i][j] > 0]
        n_present = len(values)

        if n_present == n_samples:
            features_in_all += 1
        if n_present == 1:
            features_in_one += 1

        if n_present >= 2:
            import math

            mean = sum(values) / n_present
            if mean > 0:
                std = math.sqrt(sum((v - mean) ** 2 for v in values) / (n_present - 1))
                cv = std / mean * 100
                feature_cvs.append(cv)

    median_cv = sorted(feature_cvs)[len(feature_cvs) // 2] if feature_cvs else 0

    return {
        "n_samples": n_samples,
        "n_features": n_features,
        "completeness_pct": completeness,
        "median_cv_pct": median_cv,
        "features_in_all_samples": features_in_all,
        "features_in_one_sample": features_in_one,
    }


def run_rust_aggregation(
    input_tsv: str | Path,
    output_dir: str | Path,
    dataset: str,
    data_type: str = "lfq",
    threads: int = 1,
    binary: str | None = None,
) -> Path:
    """Run the Rust aggregation_engine binary.

    Parameters
    ----------
    input_tsv : path
        Input search result TSV (DIA-NN or SAGE format).
    output_dir : path
        Output directory for aggregated results.
    dataset : str
        Dataset name for labeling.
    data_type : str
        'lfq' (DIA-NN) or 'psm' (SAGE).
    threads : int
        Number of threads.
    binary : str, optional
        Path to aggregation_engine binary. Auto-detected if None.

    Returns
    -------
    Path
        Output directory.
    """
    if binary is None:
        # Try to find the binary relative to this package
        pkg_dir = Path(__file__).parent.parent
        candidates = [
            pkg_dir / "rust" / "aggregation_engine" / "target" / "release" / "aggregation_engine",
        ]
        for c in candidates:
            if c.exists():
                binary = str(c)
                break

    if binary is None:
        raise FileNotFoundError(
            "aggregation_engine binary not found. "
            "Build with: cd rust/aggregation_engine && cargo build --release"
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        binary,
        "--input",
        str(input_tsv),
        "--output",
        str(output_dir),
        "--dataset",
        dataset,
        "--data-type",
        data_type,
        "--threads",
        str(threads),
    ]

    logger.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error("aggregation_engine failed: %s", result.stderr)
        raise RuntimeError(f"aggregation_engine failed: {result.stderr}")

    logger.info("Aggregation complete: %s", output_dir)
    return output_dir


def print_cohort_summary(samples: list[SampleResult]) -> None:
    """Print a summary table of all samples."""
    print(f"\n{'=' * 70}")
    print(f"  Cohort Summary ({len(samples)} samples)")
    print(f"{'=' * 70}")
    print(f"  {'Sample':<25} {'Peptides':>10} {'Genes':>8} {'Precursors':>12}")
    print(f"  {'─' * 55}")

    total_pep = 0
    total_gene = 0
    for s in sorted(samples, key=lambda x: x.sample_id):
        print(
            f"  {s.sample_id:<25} {len(s.peptides):>10,} {len(s.genes):>8,} {s.n_precursors:>12,}"
        )
        total_pep += len(s.peptides)
        total_gene += len(s.genes)

    print(f"  {'─' * 55}")
    print(f"  {'Mean':<25} {total_pep // len(samples):>10,} {total_gene // len(samples):>8,}")

    stats = compute_cross_sample_stats(samples, "gene")
    print("\n  Cross-sample gene stats:")
    print(f"    Total unique genes:  {stats['n_features']:,}")
    print(f"    In all samples:      {stats['features_in_all_samples']:,}")
    print(f"    In only 1 sample:    {stats['features_in_one_sample']:,}")
    print(f"    Completeness:        {stats['completeness_pct']:.1f}%")
    print(f"    Median CV:           {stats['median_cv_pct']:.1f}%")
    print(f"{'=' * 70}")
