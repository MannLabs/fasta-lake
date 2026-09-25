"""
FastaLake v4 — Comprehensive failure mode analysis.

For EVERY step of the pipeline, enumerate:
1. What can go wrong (failure modes)
2. How would we detect it (validation checks)
3. What does it look like when it's right (expected metrics)

This module implements the validation checks as executable functions.
Each returns a ValidationResult with pass/fail + details.
"""

from __future__ import annotations

import csv
import logging
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from fasta_lake.helpers import count_fasta_headers

from ..dianovo import parse_dianovo_prediction

logger = logging.getLogger(__name__)


@dataclass
class ValidationResult:
    """Result of a validation check."""

    stage: str
    check: str
    passed: bool
    message: str
    details: dict = field(default_factory=dict)
    severity: str = "error"  # "error", "warning", "info"


# ═══════════════════════════════════════════════════════════════════════════════
# STAGE -1: INPUT RAW FILES / mzML
# ═══════════════════════════════════════════════════════════════════════════════
# Failure modes:
# - mzML not found / corrupted / truncated
# - Wrong acquisition mode (user says DDA but file is DIA, or vice versa)
# - Wrong instrument / different MS parameters across samples
# - Mixed acquisition modes in the same batch
# - Centroid vs profile mode (some tools need centroid)
# - File size suspiciously small (empty run) or large (wrong settings)
# - Retention time range unexpected (very short run = problem)


def validate_mzml(mzml_path: str | Path) -> list[ValidationResult]:
    """Validate an mzML file before processing."""
    results = []
    p = Path(mzml_path)

    # Check file exists and is readable
    if not p.exists():
        results.append(
            ValidationResult(
                "input",
                "mzml_exists",
                False,
                f"mzML file not found: {p}",
                severity="error",
            )
        )
        return results

    # Check file size
    size_mb = p.stat().st_size / (1024**2)
    if size_mb < 1:
        results.append(
            ValidationResult(
                "input",
                "mzml_size",
                False,
                f"mzML suspiciously small: {size_mb:.1f} MB (empty run?)",
                details={"size_mb": size_mb},
                severity="warning",
            )
        )
    elif size_mb > 50000:
        results.append(
            ValidationResult(
                "input",
                "mzml_size",
                False,
                f"mzML very large: {size_mb:.0f} MB (check settings)",
                details={"size_mb": size_mb},
                severity="warning",
            )
        )
    else:
        results.append(
            ValidationResult(
                "input",
                "mzml_size",
                True,
                f"mzML size OK: {size_mb:.1f} MB",
                details={"size_mb": size_mb},
                severity="info",
            )
        )

    # Check file extension
    if p.suffix.lower() not in (".mzml", ".mzxml"):
        results.append(
            ValidationResult(
                "input",
                "mzml_format",
                False,
                f"Unexpected file extension: {p.suffix} (expected .mzML)",
                severity="warning",
            )
        )

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 0: DE NOVO PREDICTIONS (AlphaNovo / DIANovo / Casanovo)
# ═══════════════════════════════════════════════════════════════════════════════
# Failure modes:
# - CSV file empty or has wrong columns
# - Column names vary by engine version (pred_seq vs peptide_prediction_detokenized_unmodified)
# - Scores in unexpected range (0-1 vs 0-100 vs log-scale)
# - Peptides contain modifications that aren't stripped (brackets, masses)
# - Peptides contain non-standard amino acids (U, X, B, Z, J, O)
# - Windows line endings (\r\n) corrupt header detection
# - Very few predictions (de novo failed silently)
# - All scores identical (model degenerate)
# - Peptide length distribution wrong (all very short = chimeric spectra)
# - File encoding issues (UTF-8 BOM, Latin-1, etc.)


def validate_denovo_csv(csv_path: str | Path) -> list[ValidationResult]:
    """Validate a de novo prediction CSV before processing."""
    results = []
    p = Path(csv_path)

    if not p.exists():
        results.append(
            ValidationResult(
                "denovo",
                "csv_exists",
                False,
                f"CSV not found: {p}",
                severity="error",
            )
        )
        return results

    # Read first lines to detect format
    with open(p, "rb") as f:
        first_bytes = f.read(1000)

    # Check for Windows line endings
    if b"\r\n" in first_bytes:
        results.append(
            ValidationResult(
                "denovo",
                "line_endings",
                False,
                "Windows line endings (\\r\\n) detected — may cause parsing issues",
                severity="warning",
            )
        )

    # Check for BOM
    if first_bytes.startswith(b"\xef\xbb\xbf"):
        results.append(
            ValidationResult(
                "denovo",
                "encoding",
                False,
                "UTF-8 BOM detected — may cause header mismatch",
                severity="warning",
            )
        )

    # Parse CSV and check structure
    try:
        with open(p, newline="") as f:
            reader = csv.reader(f)
            raw_headers = next(reader)
            # Strip \r from headers
            headers = [h.strip().strip("\r") for h in raw_headers]
    except Exception as e:
        results.append(
            ValidationResult(
                "denovo",
                "csv_readable",
                False,
                f"Cannot parse CSV: {e}",
                severity="error",
            )
        )
        return results

    # Detect engine from columns
    alphanovo_cols = {"peptide_prediction_detokenized_unmodified", "score"}
    dianovo_cols = {"pred_seq", "pred_prob"}
    casanovo_cols = {"sequence", "score"}

    headers_set = set(headers)
    engine = "unknown"
    if alphanovo_cols.issubset(headers_set):
        engine = "alphanovo"
    elif dianovo_cols.issubset(headers_set):
        engine = "dianovo"
    elif casanovo_cols.issubset(headers_set):
        engine = "casanovo"

    results.append(
        ValidationResult(
            "denovo",
            "engine_detected",
            engine != "unknown",
            f"Detected engine: {engine}"
            if engine != "unknown"
            else f"Unknown CSV format. Columns: {headers[:5]}",
            details={"engine": engine, "columns": headers},
            severity="error" if engine == "unknown" else "info",
        )
    )

    # Count rows and check data quality
    n_rows = 0
    scores = []
    lengths = []
    non_standard = 0
    empty_peptides = 0
    invalid_dianovo = Counter()

    with open(p) as f:
        reader = csv.DictReader(f)
        # Normalize header names
        if reader.fieldnames:
            reader.fieldnames = [h.strip().strip("\r") for h in reader.fieldnames]

        for row in reader:
            n_rows += 1
            if n_rows > 10000:
                break  # sample first 10K rows

            # Get peptide
            pep = ""
            if engine == "alphanovo":
                pep = row.get("peptide_prediction_detokenized_unmodified", "")
            elif engine == "dianovo":
                try:
                    pep, dianovo_score = parse_dianovo_prediction(
                        row.get("pred_seq") or "", row.get("pred_prob") or ""
                    )
                except ValueError as exc:
                    invalid_dianovo[str(exc)] += 1
                    continue
            elif engine == "casanovo":
                pep = row.get("sequence", "")

            if not pep:
                empty_peptides += 1
                continue

            lengths.append(len(pep))

            # Check for non-standard amino acids
            if re.search(r"[^ACDEFGHIKLMNPQRSTVWY]", pep.upper()):
                non_standard += 1

            # Get score
            try:
                if engine == "dianovo":
                    score = dianovo_score
                else:
                    score = float(row.get("score", 0))
                scores.append(score)
            except (ValueError, TypeError):
                pass

    results.append(
        ValidationResult(
            "denovo",
            "row_count",
            n_rows > 100,
            f"Rows: {n_rows:,}" + (" (sampled first 10K)" if n_rows >= 10000 else ""),
            details={"n_rows": n_rows},
            severity="warning" if n_rows < 100 else "info",
        )
    )

    if engine == "dianovo":
        invalid = sum(invalid_dianovo.values())
        results.append(
            ValidationResult(
                "denovo",
                "dianovo_complete_predictions",
                invalid == 0,
                f"DIANovo: {invalid}/{n_rows} sampled rows excluded as unresolved/invalid; "
                "zero residue probabilities are retained",
                details={"excluded": invalid, "reasons": dict(invalid_dianovo)},
                severity="error" if invalid == n_rows else "warning" if invalid else "info",
            )
        )

    if empty_peptides > 0:
        results.append(
            ValidationResult(
                "denovo",
                "empty_peptides",
                empty_peptides < n_rows * 0.5,
                (
                    f"Empty peptides: {empty_peptides}/{n_rows} "
                    f"({empty_peptides / max(n_rows, 1) * 100:.0f}%)"
                ),
                details={"empty": empty_peptides},
                severity="warning" if empty_peptides < n_rows * 0.5 else "error",
            )
        )

    if non_standard > 0:
        results.append(
            ValidationResult(
                "denovo",
                "non_standard_aa",
                non_standard < n_rows * 0.1,
                f"Non-standard amino acids in {non_standard} peptides",
                details={"count": non_standard},
                severity="warning",
            )
        )

    if lengths:
        import statistics

        med_len = statistics.median(lengths)
        results.append(
            ValidationResult(
                "denovo",
                "peptide_length",
                8 <= med_len <= 25,
                f"Median peptide length: {med_len:.0f} aa (range {min(lengths)}-{max(lengths)})",
                details={"median": med_len, "min": min(lengths), "max": max(lengths)},
                severity="warning" if med_len < 8 or med_len > 25 else "info",
            )
        )

    if scores:
        import statistics

        med_score = statistics.median(scores)
        score_range = (min(scores), max(scores))
        # Check if scores are in expected range
        if score_range[1] > 1.1:
            results.append(
                ValidationResult(
                    "denovo",
                    "score_range",
                    False,
                    f"Scores exceed 1.0 (max={score_range[1]:.2f}) — wrong scale?",
                    severity="warning",
                )
            )
        elif med_score < 0.01:
            results.append(
                ValidationResult(
                    "denovo",
                    "score_range",
                    False,
                    f"Median score very low ({med_score:.4f}) — de novo model may have failed",
                    severity="warning",
                )
            )
        else:
            results.append(
                ValidationResult(
                    "denovo",
                    "score_range",
                    True,
                    f"Score range: {score_range[0]:.3f}-{score_range[1]:.3f}, "
                    f"median={med_score:.3f}",
                    details={"median": med_score, "range": score_range},
                    severity="info",
                )
            )

        # Check if all scores identical (degenerate model)
        if len(set(round(s, 4) for s in scores)) <= 3:
            results.append(
                ValidationResult(
                    "denovo",
                    "score_diversity",
                    False,
                    f"Only {len(set(round(s, 4) for s in scores))} unique scores — degenerate "
                    f"model?",
                    severity="warning",
                )
            )

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 0: LAKE BUILDING (Rust lake_builder)
# ═══════════════════════════════════════════════════════════════════════════════
# Failure modes:
# - Source FASTA file not found / corrupted / truncated
# - FASTA has no > characters (wrong file, maybe nucleotide)
# - FASTA has * (stop codons) in sequences — need clean_sequence
# - Duplicate accessions across sources with different sequences
# - Wildly different sizes between sources (millions vs thousands)
# - Sequences contain non-amino-acid characters (numbers, spaces)
# - Empty sequences (header with no sequence lines)
# - Very long headers that break tools (>10K characters)
# - Mixed nucleotide and protein sequences in same file
# - Gzipped file but not recognized


def validate_source_fasta(fasta_path: str | Path) -> list[ValidationResult]:
    """Validate a source FASTA file before lake building."""
    results = []
    p = Path(fasta_path)

    if not p.exists():
        results.append(
            ValidationResult(
                "lake_build",
                "fasta_exists",
                False,
                f"FASTA not found: {p}",
                severity="error",
            )
        )
        return results

    # Quick scan: count entries, check for issues
    n_entries = 0
    n_empty_seq = 0
    n_stop_codons = 0
    n_non_aa = 0
    n_nucleotide_like = 0
    max_header_len = 0
    min_seq_len = float("inf")
    max_seq_len = 0
    seq_lens = []
    header_prefixes = Counter()

    current_seq_len = 0
    in_seq = False

    # Only scan first 50K entries for speed
    with open(p) as f:
        for line in f:
            if line.startswith(">"):
                if in_seq:
                    if current_seq_len == 0:
                        n_empty_seq += 1
                    seq_lens.append(current_seq_len)
                    min_seq_len = min(min_seq_len, current_seq_len)
                    max_seq_len = max(max_seq_len, current_seq_len)

                n_entries += 1
                if n_entries > 50000:
                    break

                header = line[1:].strip()
                max_header_len = max(max_header_len, len(header))

                # Track prefix
                prefix = header.split("|")[0].split()[0][:10]
                header_prefixes[prefix] += 1

                current_seq_len = 0
                in_seq = True
            else:
                seq_line = line.strip().upper()
                current_seq_len += len(seq_line)

                if "*" in seq_line:
                    n_stop_codons += 1

                # Check for nucleotide-only sequences
                if seq_line and all(c in "ACGTURYSWKMBDHVN" for c in seq_line):
                    n_nucleotide_like += 1

                # Check for non-amino-acid characters
                if re.search(r"[^ACDEFGHIKLMNPQRSTVWY*XU\n\r ]", seq_line):
                    n_non_aa += 1

    results.append(
        ValidationResult(
            "lake_build",
            "entry_count",
            n_entries > 0,
            f"Entries: {n_entries:,}" + (" (sampled first 50K)" if n_entries >= 50000 else ""),
            details={"n_entries": n_entries},
            severity="error" if n_entries == 0 else "info",
        )
    )

    if n_empty_seq > 0:
        results.append(
            ValidationResult(
                "lake_build",
                "empty_sequences",
                False,
                f"Empty sequences: {n_empty_seq} (headers with no sequence)",
                severity="warning",
            )
        )

    if n_stop_codons > n_entries * 0.1:
        results.append(
            ValidationResult(
                "lake_build",
                "stop_codons",
                False,
                f"Stop codons (*) in {n_stop_codons} sequences — need clean_sequence()",
                severity="warning",
            )
        )

    if n_nucleotide_like > n_entries * 0.5:
        results.append(
            ValidationResult(
                "lake_build",
                "nucleotide_check",
                False,
                f"{n_nucleotide_like} sequences look like nucleotide (not protein!) — wrong file?",
                severity="error",
            )
        )

    if max_header_len > 10000:
        results.append(
            ValidationResult(
                "lake_build",
                "header_length",
                False,
                f"Very long header: {max_header_len} chars (may break tools)",
                severity="warning",
            )
        )

    if seq_lens:
        import statistics

        med = statistics.median(seq_lens)
        results.append(
            ValidationResult(
                "lake_build",
                "sequence_length",
                10 < med < 10000,
                f"Sequence length: median={med:.0f}, range={min_seq_len}-{max_seq_len}",
                details={"median": med, "min": min_seq_len, "max": max_seq_len},
                severity="info",
            )
        )

    # Header format distribution
    results.append(
        ValidationResult(
            "lake_build",
            "header_formats",
            True,
            f"Header prefixes: {dict(header_prefixes.most_common(5))}",
            details={"prefixes": dict(header_prefixes)},
            severity="info",
        )
    )

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 1: EVIDENCE LAKE (Rust fasta_extractor)
# ═══════════════════════════════════════════════════════════════════════════════
# Failure modes:
# - 0% hit rate (de novo predictions don't match database at all)
#   → wrong organism? wrong database? wrong de novo settings?
# - Very high hit rate (>50%) on a metaproteomics database
#   → may be matching too many short/ambiguous peptides
# - Evidence lake much larger than expected
#   → top_percent not set? min_length too short?
# - Evidence lake too small (< canonical proteome size)
#   → too aggressive filtering? wrong database?
# - All peptides match one protein (contamination?)
# - Score distribution after filtering is flat (no discrimination)
# - I/L normalization off when it should be on (or vice versa)


def validate_evidence_lake(
    evidence_fasta: str | Path,
    denovo_stats: dict | None = None,
    lake_stats: dict | None = None,
) -> list[ValidationResult]:
    """Validate evidence lake output from fasta_extractor."""
    results = []
    p = Path(evidence_fasta)

    if not p.exists():
        results.append(
            ValidationResult(
                "evidence_lake",
                "exists",
                False,
                f"Evidence lake not found: {p}",
                severity="error",
            )
        )
        return results

    # Count entries
    n_entries = count_fasta_headers(p)

    results.append(
        ValidationResult(
            "evidence_lake",
            "entry_count",
            n_entries > 0,
            f"Evidence lake: {n_entries:,} proteins",
            details={"n_entries": n_entries},
            severity="error" if n_entries == 0 else "info",
        )
    )

    # Check hit rate if we have de novo stats
    if lake_stats and "proteins_searched" in lake_stats:
        hit_rate = n_entries / lake_stats["proteins_searched"] * 100
        results.append(
            ValidationResult(
                "evidence_lake",
                "hit_rate",
                0.1 < hit_rate < 50,
                f"Hit rate: {hit_rate:.2f}% ({n_entries:,}/{lake_stats['proteins_searched']:,})",
                details={"hit_rate": hit_rate},
                severity="warning" if hit_rate < 0.1 or hit_rate > 50 else "info",
            )
        )

    # Check file size relative to entry count
    size_mb = p.stat().st_size / (1024**2)
    avg_entry_size = size_mb * 1024 / max(n_entries, 1)  # KB per entry
    results.append(
        ValidationResult(
            "evidence_lake",
            "file_size",
            True,
            f"File size: {size_mb:.1f} MB ({avg_entry_size:.1f} KB/entry)",
            details={"size_mb": size_mb, "avg_kb_per_entry": avg_entry_size},
            severity="info",
        )
    )

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 2: PER-SAMPLE EXTRACTION
# ═══════════════════════════════════════════════════════════════════════════════
# Failure modes:
# - Sample CSV has 0 matching peptides (wrong sample? wrong format?)
# - Per-sample FASTA much larger or smaller than expected
# - Per-sample FASTA has no proteins from canonical (--include failed?)
# - All samples produce identical FASTAs (no sample specificity)
# - Some samples have 10x fewer proteins than others (outlier detection)


def validate_persample_fasta(
    persample_fasta: str | Path,
    canonical_fasta: str | Path | None = None,
) -> list[ValidationResult]:
    """Validate per-sample FASTA from Stage 2."""
    results = []
    p = Path(persample_fasta)

    if not p.exists():
        results.append(
            ValidationResult(
                "persample",
                "exists",
                False,
                f"Per-sample FASTA not found: {p}",
                severity="error",
            )
        )
        return results

    n_entries = count_fasta_headers(p)
    results.append(
        ValidationResult(
            "persample",
            "entry_count",
            n_entries > 0,
            f"Per-sample: {n_entries:,} proteins",
            details={"n_entries": n_entries},
            severity="error" if n_entries == 0 else "info",
        )
    )

    # Check canonical inclusion
    if canonical_fasta:
        can_path = Path(canonical_fasta)
        if can_path.exists():
            n_canonical = count_fasta_headers(can_path)
            if n_entries < n_canonical:
                results.append(
                    ValidationResult(
                        "persample",
                        "canonical_inclusion",
                        False,
                        f"Per-sample ({n_entries:,}) smaller than canonical ({n_canonical:,}) — "
                        f"--include may have failed",
                        severity="warning",
                    )
                )
            else:
                results.append(
                    ValidationResult(
                        "persample",
                        "canonical_inclusion",
                        True,
                        f"Per-sample ({n_entries:,}) >= canonical ({n_canonical:,})",
                        severity="info",
                    )
                )

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 3: INFERENCE
# ═══════════════════════════════════════════════════════════════════════════════
# Failure modes:
# - Strategy produces 0 proteins (all filtered out)
# - Reduction too aggressive (>100x) or too weak (<2x)
# - Headers truncated to first token (the Stage 3 bug we found)
# - Species budget assigns everything to one species
# - Razor assigns shared peptides illogically


def validate_inference_output(
    output_fasta: str | Path,
    stats_json: str | Path | None = None,
) -> list[ValidationResult]:
    """Validate Stage 3 inference output."""
    results = []
    p = Path(output_fasta)

    if not p.exists():
        results.append(
            ValidationResult(
                "inference",
                "exists",
                False,
                f"Inference output not found: {p}",
                severity="error",
            )
        )
        return results

    # Check if headers are truncated (the Stage 3 bug)
    n_full_headers = 0
    n_truncated = 0
    with open(p) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                if " " in header:
                    n_full_headers += 1
                else:
                    n_truncated += 1

    results.append(
        ValidationResult(
            "inference",
            "header_preservation",
            n_full_headers > n_truncated,
            f"Headers: {n_full_headers} full, {n_truncated} truncated (first-token only)",
            details={"full": n_full_headers, "truncated": n_truncated},
            severity="warning" if n_truncated > n_full_headers else "info",
        )
    )

    # Check stats if available
    if stats_json:
        import json

        sp = Path(stats_json)
        if sp.exists():
            with open(sp) as f:
                stats = json.load(f)
            reduction = stats.get("reduction_factor", 0)
            selected = stats.get("selected_proteins", 0)
            results.append(
                ValidationResult(
                    "inference",
                    "reduction",
                    1 < reduction < 1000,
                    f"Selected {selected:,} proteins ({reduction:.1f}x reduction)",
                    details=stats,
                    severity="warning" if reduction < 1 or reduction > 1000 else "info",
                )
            )

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# FULL PIPELINE VALIDATION
# ═══════════════════════════════════════════════════════════════════════════════


def validate_pipeline(
    mzml_files: list[str | Path] | None = None,
    denovo_csvs: list[str | Path] | None = None,
    source_fastas: list[str | Path] | None = None,
    evidence_fasta: str | Path | None = None,
    persample_fastas: list[str | Path] | None = None,
    inference_fastas: list[str | Path] | None = None,
) -> list[ValidationResult]:
    """Run all validation checks on available pipeline outputs."""
    all_results = []

    if mzml_files:
        for f in mzml_files[:5]:  # sample first 5
            all_results.extend(validate_mzml(f))

    if denovo_csvs:
        for f in denovo_csvs[:5]:
            all_results.extend(validate_denovo_csv(f))

    if source_fastas:
        for f in source_fastas:
            all_results.extend(validate_source_fasta(f))

    if evidence_fasta:
        all_results.extend(validate_evidence_lake(evidence_fasta))

    if persample_fastas:
        for f in persample_fastas[:5]:
            all_results.extend(validate_persample_fasta(f))

    if inference_fastas:
        for f in inference_fastas[:5]:
            all_results.extend(validate_inference_output(f))

    return all_results


def print_validation(results: list[ValidationResult]) -> None:
    """Pretty-print validation results."""
    errors = [r for r in results if not r.passed and r.severity == "error"]
    warnings = [r for r in results if not r.passed and r.severity == "warning"]
    passed = [r for r in results if r.passed]

    print(f"\n{'=' * 60}")
    print(f"  Pipeline Validation: {len(results)} checks")
    print(f"{'=' * 60}")
    print(f"  PASSED:   {len(passed)}")
    print(f"  WARNINGS: {len(warnings)}")
    print(f"  ERRORS:   {len(errors)}")

    if errors:
        print(f"\n  {'─' * 56}")
        print("  ERRORS:")
        for r in errors:
            print(f"    [{r.stage}] {r.check}: {r.message}")

    if warnings:
        print(f"\n  {'─' * 56}")
        print("  WARNINGS:")
        for r in warnings:
            print(f"    [{r.stage}] {r.check}: {r.message}")

    if passed:
        print(f"\n  {'─' * 56}")
        print("  PASSED:")
        for r in passed:
            print(f"    [{r.stage}] {r.check}: {r.message}")


# ═══════════════════════════════════════════════════════════════════════════════
# AUTO-DETECTION: Infer full pipeline config from inputs
# ═══════════════════════════════════════════════════════════════════════════════


@dataclass
class ExperimentDetection:
    """Auto-detected experiment configuration."""

    denovo_engine: str = "unknown"  # alphanovo, dianovo, casanovo
    acquisition: str = "unknown"  # dda, dia
    search_engine: str = "sage"  # sage, diann, msfragger
    inference_strategy: str = "species_budget"
    needs_header_healing: bool = False
    needs_decoys: bool = False
    confidence: str = "low"  # low, medium, high
    reasons: list = field(default_factory=list)


def detect_experiment(
    predictions_dir: str | Path | None = None,
    mzml_dir: str | Path | None = None,
) -> ExperimentDetection:
    """Auto-detect experiment type from de novo CSVs and mzML filenames.

    The detection chain:
    1. CSV columns → de novo engine (AlphaNovo / DIANovo / Casanovo)
    2. De novo engine → acquisition mode (AlphaNovo → DDA, DIANovo → DIA)
    3. Acquisition → search engine (DDA → SAGE, DIA → DIA-NN)
    4. Search engine → FASTA preparation needs

    Optionally uses mzML filename patterns as confirmation:
    - "100SPD", "SWATH", "diaPASEF" → DIA
    - Absence of DIA markers → DDA (default)

    Parameters
    ----------
    predictions_dir : path, optional
        Directory containing de novo prediction CSVs.
    mzml_dir : path, optional
        Directory containing mzML files.

    Returns
    -------
    ExperimentDetection
        Auto-detected configuration with confidence level and reasoning.
    """
    det = ExperimentDetection()

    # Step 1: Detect de novo engine from CSV columns
    if predictions_dir is not None:
        pred_path = Path(predictions_dir)
        csvs = list(pred_path.glob("*.csv")) + list(pred_path.glob("*.CSV"))
        if csvs:
            # Check first CSV
            results = validate_denovo_csv(csvs[0])
            for r in results:
                if r.check == "engine_detected" and r.details.get("engine") != "unknown":
                    det.denovo_engine = r.details["engine"]
                    det.reasons.append(f"De novo engine: {det.denovo_engine} (from CSV columns)")

    # Step 2: De novo engine → acquisition mode
    engine_to_acquisition = {
        "alphanovo": "dda",
        "dianovo": "dia",
        "casanovo": "dda",  # Casanovo is typically DDA
    }
    if det.denovo_engine in engine_to_acquisition:
        det.acquisition = engine_to_acquisition[det.denovo_engine]
        det.confidence = "high"
        det.reasons.append(f"Acquisition: {det.acquisition} (from {det.denovo_engine})")

    # Step 3: Cross-check with mzML filenames
    if mzml_dir is not None:
        mzml_path = Path(mzml_dir)
        mzmls = list(mzml_path.glob("*.mzML")) + list(mzml_path.glob("*.mzml"))
        dia_patterns = ["100SPD", "SWATH", "diaPASEF", "dia_PASEF", "_DIA_", "_dia_"]
        dia_count = sum(1 for m in mzmls if any(p.lower() in m.name.lower() for p in dia_patterns))
        if mzmls:
            fraction_dia = dia_count / len(mzmls)
            if fraction_dia > 0.5:
                if det.acquisition == "dda":
                    det.reasons.append(
                        f"WARNING: mzML names suggest DIA ({dia_count}/{len(mzmls)}) "
                        f"but de novo engine ({det.denovo_engine}) suggests DDA"
                    )
                    det.confidence = "low"
                else:
                    det.reasons.append(f"mzML names confirm DIA ({dia_count}/{len(mzmls)})")

    # Step 4: Acquisition → recommended search engine
    if det.acquisition == "dia":
        det.search_engine = "diann"
        det.needs_header_healing = True
        det.needs_decoys = False
        det.reasons.append("Search engine: DIA-NN (recommended for DIA)")
    else:
        det.search_engine = "sage"
        det.needs_header_healing = False
        det.needs_decoys = False  # SAGE generates internal decoys
        det.reasons.append("Search engine: SAGE (recommended for DDA)")

    return det
