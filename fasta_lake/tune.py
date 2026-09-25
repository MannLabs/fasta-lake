"""
FastaLake tune — parameter tuning from existing SAGE results.

Takes a directory of SAGE entrapment searches, applies a grid of filters
post-hoc, and recommends parameters that achieve a user-specified target FDP.

No re-search required — operates on existing results.sage.tsv files.

CLI:
    fasta-lake tune --results-dir PATH --target-fdp 0.02

Output:
    tuning_report.csv   — per-filter FDP, target count, pareto verdict
    recommended.json    — the winning filter as a SAGE config diff
"""

from __future__ import annotations

import csv
import json
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from fasta_lake.helpers import require_columns

# Entrapment configs — shared with sweep scripts
ENTRAP_PREFIXES = {
    "shuffled": "SHUF_",
    "plant": "PLANT_",
    "archaea": "ARCH_",
}


@dataclass
class FilterSpec:
    """One parameter combination to test."""

    name: str
    q_threshold: float
    min_length: int
    min_score: float


def default_filter_grid() -> list[FilterSpec]:
    """Standard grid. len × q × score cartesian product."""
    grid = []
    for q in [0.01, 0.005, 0.001]:
        for L in [7, 8, 9, 10]:
            for s in [0.0, 0.5, 0.7]:
                name = f"q{q:.3f}_len{L}_score{s:.1f}"
                grid.append(FilterSpec(name, q, L, s))
    return grid


def apply_filter(rows: list[dict], spec: FilterSpec, entrap_prefix: str) -> dict:
    """Count target / decoy / entrapment PSMs under a filter."""
    target = decoy = entrap = 0
    for r in rows:
        if r["q"] > spec.q_threshold:
            continue
        if r["len"] < spec.min_length:
            continue
        if r["score"] < spec.min_score:
            continue
        if r["is_decoy"]:
            decoy += 1
        elif entrap_prefix in r["proteins"]:
            entrap += 1
        else:
            target += 1
    return {"target": target, "decoy": decoy, "entrap": entrap}


def parse_tsv(path: Path) -> list[dict]:
    """Parse required Sage fields for the historical spectrum-q tuning diagnostic.

    Return q, length, discriminant score, protein text and the legacy decoy flag
    for each row. Missing columns and invalid numeric fields raise errors.
    """
    rows = []
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        require_columns(
            r.fieldnames, ("spectrum_q", "peptide_len", "sage_discriminant_score", "proteins"), path
        )
        for row in r:
            p = row["proteins"]
            rows.append(
                {
                    "q": float(row["spectrum_q"]),
                    "len": int(row["peptide_len"]),
                    "score": float(row["sage_discriminant_score"]),
                    "proteins": p,
                    "is_decoy": "rev_" in p,
                }
            )
    return rows


def compute_fdp_combined(entrap: int, target_plus_decoy: int, r: float) -> float | None:
    """Noble 2025 combined-method FDP as percentage; ``None`` when undefined.

    A sample with no accepted PSMs has no FDP, not an FDP of zero; returning
    0.0 pulled the median that picks the winning filter towards zero.
    """
    denom = target_plus_decoy + entrap
    if denom == 0:
        return None
    return entrap * (1 + 1 / r) / denom * 100


def find_sample_results(results_dir: Path) -> Iterator[tuple[str, str, Path]]:
    """Yield (sample_id, entrap_type, results.sage.tsv path) for each sample."""
    # Expected layout: {results_dir}/{sample_id}/{entrap_type}/results.sage.tsv
    # Where entrap_type in {entrapment, entrapment_plant, entrapment_archaea}
    subdir_map = {
        "entrapment": "shuffled",
        "entrapment_plant": "plant",
        "entrapment_archaea": "archaea",
    }
    for sample_dir in sorted(results_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        for subdir_name, entrap_type in subdir_map.items():
            tsv = sample_dir / subdir_name / "results.sage.tsv"
            if tsv.exists():
                yield sample_dir.name, entrap_type, tsv


def tune(
    results_dir: Path,
    target_fdp_pct: float = 2.0,
    entrap_db_sizes: dict[str, int] | None = None,
    target_db_size_fn=None,
    grid: list[FilterSpec] | None = None,
) -> dict:
    """
    Tune parameters by post-hoc filter sweep.

    Parameters
    ----------
    results_dir
        Directory with per-sample subdirs (each containing entrapment/, entrapment_plant/,
        entrapment_archaea/ with results.sage.tsv).
    target_fdp_pct
        Maximum acceptable median shuffled FDP (percent). Winner is the filter
        with highest target PSM count that satisfies this threshold.
    entrap_db_sizes
        Map entrap_type -> protein count. Defaults to FastaLake standard sizes.
    target_db_size_fn
        Callable (sample_id) -> int returning the target FASTA protein count.
        If None, assumes 75000 (typical MP metaproteomics).
    grid
        List of FilterSpec to evaluate. Defaults to `default_filter_grid()`.

    Returns
    -------
    dict with keys: all_results (per-filter per-entrap summary), winner (best filter).
    """
    if grid is None:
        grid = default_filter_grid()
    if entrap_db_sizes is None:
        entrap_db_sizes = {"shuffled": 6500, "plant": 16418, "archaea": 505}
    if target_db_size_fn is None:

        def target_db_size_fn(_sid):
            """Supply the historical 75,000-protein fallback when no size callback exists."""
            return 75000

    # Collect per-sample rows
    per_sample: dict[tuple[str, str], list[dict]] = {}
    for sample_id, entrap_type, tsv in find_sample_results(results_dir):
        per_sample[(sample_id, entrap_type)] = parse_tsv(tsv)

    if not per_sample:
        raise RuntimeError(f"No entrapment result TSVs found under {results_dir}")

    # Evaluate each filter spec × each entrapment type across all samples
    summary = []
    for spec in grid:
        for entrap_type, entrap_db in entrap_db_sizes.items():
            prefix = ENTRAP_PREFIXES[entrap_type]
            target_counts = []
            fdps = []
            for (sid, etype), rows in per_sample.items():
                if etype != entrap_type:
                    continue
                target_db = target_db_size_fn(sid)
                r = entrap_db / target_db if target_db else 0.0
                if r <= 0:
                    continue
                c = apply_filter(rows, spec, prefix)
                fdp = compute_fdp_combined(c["entrap"], c["target"] + c["decoy"], r)
                if fdp is None:
                    continue  # undefined for this sample: not a zero
                target_counts.append(c["target"])
                fdps.append(fdp)

            if not fdps:
                continue
            summary.append(
                {
                    "filter": spec.name,
                    "q_threshold": spec.q_threshold,
                    "min_length": spec.min_length,
                    "min_score": spec.min_score,
                    "entrap_type": entrap_type,
                    "n_samples": len(fdps),
                    "median_target": statistics.median(target_counts),
                    "median_fdp_pct": statistics.median(fdps),
                    "mean_fdp_pct": statistics.mean(fdps),
                    "p90_fdp_pct": sorted(fdps)[int(0.9 * len(fdps))],
                }
            )

    # Pick winner: highest target count where shuffled median FDP ≤ target_fdp_pct
    shuffled_rows = [r for r in summary if r["entrap_type"] == "shuffled"]
    passing = [r for r in shuffled_rows if r["median_fdp_pct"] <= target_fdp_pct]

    if not passing:
        # No filter passes → recommend the tightest
        winner = min(shuffled_rows, key=lambda r: r["median_fdp_pct"]) if shuffled_rows else None
        winner_verdict = "NO_FILTER_PASSES — showing tightest available"
    else:
        winner = max(passing, key=lambda r: r["median_target"])
        winner_verdict = f"PASSES @ target_fdp_pct={target_fdp_pct}"

    return {
        "all_results": summary,
        "winner": winner,
        "winner_verdict": winner_verdict,
        "target_fdp_pct": target_fdp_pct,
        "n_samples": len(set(sid for sid, _ in per_sample.keys())),
        "n_entrap_types": len(set(et for _, et in per_sample.keys())),
    }


def write_reports(report: dict, out_dir: Path) -> None:
    """Write tuning CSV, recommended settings and supporting reports.

    Create the directory as needed and overwrite report files. An empty result
    list produces no new tuning_report.csv.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # CSV report
    csv_path = out_dir / "tuning_report.csv"
    if report["all_results"]:
        cols = list(report["all_results"][0].keys())
        with open(csv_path, "w") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(report["all_results"])

    # Recommended preset / filter
    rec_path = out_dir / "recommended.json"
    recommended = {
        "winner_verdict": report["winner_verdict"],
        "target_fdp_pct": report["target_fdp_pct"],
        "winner": report["winner"],
    }
    rec_path.write_text(json.dumps(recommended, indent=2, default=float))
