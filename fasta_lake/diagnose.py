"""
FastaLake diagnose — FDR calibration curves and per-sample QC.

Given a directory of SAGE entrapment searches, produces:
- Calibration curve: empirical FDP vs nominal q-value, per entrapment type
- Score distribution comparison: target vs decoy vs entrap
- Per-sample QC: flags samples with outlier FDP

CLI:
    fasta-lake diagnose --results-dir PATH -o OUTPUT_DIR [--plot]
"""

from __future__ import annotations

import csv
import json
import statistics
from pathlib import Path

from fasta_lake.tune import (
    ENTRAP_PREFIXES,
    compute_fdp_combined,
    find_sample_results,
    parse_tsv,
)

# Q-values to sweep for calibration curve
CALIBRATION_Q = [0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]


def per_sample_fdp_at_q(
    rows: list[dict],
    entrap_prefix: str,
    q_thr: float,
    min_length: int = 9,
    exclude_contaminants: bool = True,
) -> dict:
    """Compute per-sample PSM- and protein-level FDP at a given q threshold.

    Protein-level FDR uses the picked-protein method (Savitski 2015): for each
    protein, keep only the best PSM; then count target/decoy/entrap at the
    protein level. This is typically 2-5× tighter than PSM-level FDR.

    Contaminants (``CON_`` prefix, or accessions in the default contaminant DB)
    are excluded from both target and entrapment counts — they are neither
    false positives nor true targets.
    """
    from fasta_lake.contaminants import is_contaminant_id

    target = decoy = entrap = contaminant = 0
    # Protein-level: map protein -> (best_score, is_decoy, is_entrap)
    protein_best: dict[str, tuple[float, bool, bool]] = {}

    for r in rows:
        if r["q"] > q_thr:
            continue
        if r["len"] < min_length:
            continue

        first_prot = r["proteins"].split(";")[0].split(",")[0].strip()

        # Contaminants are counted separately and excluded from FDP numerator/denom
        if exclude_contaminants and is_contaminant_id(first_prot):
            contaminant += 1
            continue

        if r["is_decoy"]:
            decoy += 1
        elif entrap_prefix in r["proteins"]:
            entrap += 1
        else:
            target += 1

        # Protein-level picked scoring
        score = r.get("score", 0.0)
        existing = protein_best.get(first_prot)
        if (existing is None) or (score > existing[0]):
            is_entrap = entrap_prefix in r["proteins"]
            protein_best[first_prot] = (score, r["is_decoy"], is_entrap)

    # Protein-level counts
    prot_target = sum(1 for _, d, e in protein_best.values() if not d and not e)
    prot_decoy = sum(1 for _, d, e in protein_best.values() if d)
    prot_entrap = sum(1 for _, d, e in protein_best.values() if e)

    return {
        "target": target,
        "decoy": decoy,
        "entrap": entrap,
        "contaminant": contaminant,
        "prot_target": prot_target,
        "prot_decoy": prot_decoy,
        "prot_entrap": prot_entrap,
    }


def calibration_curve(
    results_dir: Path,
    entrap_db_sizes: dict[str, int] | None = None,
    target_db_size_fn=None,
    min_length: int = 9,
) -> dict:
    """Compute empirical FDP vs nominal q-value curve across samples."""
    if entrap_db_sizes is None:
        entrap_db_sizes = {"shuffled": 6500, "plant": 16418, "archaea": 505}
    if target_db_size_fn is None:

        def target_db_size_fn(_sid):
            """Supply the historical 75,000-protein fallback when no size callback exists."""
            return 75000

    # Collect rows
    per_sample: dict[tuple[str, str], list[dict]] = {}
    for sid, et, tsv in find_sample_results(results_dir):
        per_sample[(sid, et)] = parse_tsv(tsv)

    curve = []
    for entrap_type, db_size in entrap_db_sizes.items():
        prefix = ENTRAP_PREFIXES[entrap_type]
        for q in CALIBRATION_Q:
            psm_fdps = []
            prot_fdps = []
            psm_targets = []
            prot_targets = []
            per_sample_ratios_psm = []
            contaminants = 0
            for (sid, et), rows in per_sample.items():
                if et != entrap_type:
                    continue
                td = target_db_size_fn(sid)
                if td <= 0:
                    continue
                r = db_size / td
                c = per_sample_fdp_at_q(rows, prefix, q, min_length)
                psm_fdp = compute_fdp_combined(c["entrap"], c["target"] + c["decoy"], r)
                prot_fdp = compute_fdp_combined(
                    c["prot_entrap"], c["prot_target"] + c["prot_decoy"], r
                )
                if psm_fdp is None or prot_fdp is None:
                    continue  # undefined for this sample: not a zero
                psm_fdps.append(psm_fdp)
                prot_fdps.append(prot_fdp)
                psm_targets.append(c["target"])
                prot_targets.append(c["prot_target"])
                contaminants += c["contaminant"]
                if q > 0:
                    per_sample_ratios_psm.append(psm_fdp / (q * 100))
            if not psm_fdps:
                continue
            curve.append(
                {
                    "entrap_type": entrap_type,
                    "nominal_q": q,
                    "nominal_fdp_pct": q * 100,
                    "n_samples": len(psm_fdps),
                    "median_psm_fdp_pct": statistics.median(psm_fdps),
                    "mean_psm_fdp_pct": statistics.mean(psm_fdps),
                    "p90_psm_fdp_pct": sorted(psm_fdps)[int(0.9 * len(psm_fdps))],
                    "median_prot_fdp_pct": statistics.median(prot_fdps),
                    "mean_prot_fdp_pct": statistics.mean(prot_fdps),
                    "p90_prot_fdp_pct": sorted(prot_fdps)[int(0.9 * len(prot_fdps))],
                    "median_psm_target": statistics.median(psm_targets),
                    "median_prot_target": statistics.median(prot_targets),
                    "calibration_ratio_psm": statistics.median(per_sample_ratios_psm)
                    if per_sample_ratios_psm
                    else float("nan"),
                    "total_contaminants": contaminants,
                }
            )
    return {
        "curve": curve,
        "n_samples": len(set(s for s, _ in per_sample.keys())),
        "min_length": min_length,
    }


def outlier_samples(
    results_dir: Path, q_thr: float = 0.01, fdp_z_threshold: float = 2.0, min_length: int = 9
) -> list[dict]:
    """Identify samples with FDP > mean + z*std (outliers)."""
    per_sample_rows: dict[tuple[str, str], list[dict]] = {}
    for sid, et, tsv in find_sample_results(results_dir):
        per_sample_rows[(sid, et)] = parse_tsv(tsv)

    # Just shuffled
    entrap_type = "shuffled"
    prefix = ENTRAP_PREFIXES[entrap_type]
    all_fdps = []
    sample_fdps = {}
    for (sid, et), rows in per_sample_rows.items():
        if et != entrap_type:
            continue
        c = per_sample_fdp_at_q(rows, prefix, q_thr, min_length)
        r = 6500 / 75000
        fdp = compute_fdp_combined(c["entrap"], c["target"] + c["decoy"], r)
        if fdp is None:
            continue  # undefined for this sample: not a zero
        all_fdps.append(fdp)
        sample_fdps[sid] = (fdp, c["target"], c["entrap"], c["prot_target"], c["prot_entrap"])

    if len(all_fdps) < 2:
        return []
    mu = statistics.mean(all_fdps)
    sd = statistics.stdev(all_fdps)
    threshold = mu + fdp_z_threshold * sd
    outliers = []
    for sid, (fdp, t, e, pt, pe) in sorted(sample_fdps.items(), key=lambda x: -x[1][0]):
        if fdp > threshold:
            outliers.append(
                {
                    "sample": sid,
                    "shuffled_psm_fdp_pct": round(fdp, 3),
                    "target_psms": t,
                    "entrap_psms": e,
                    "target_proteins": pt,
                    "entrap_proteins": pe,
                    "z_score": round((fdp - mu) / sd if sd else 0, 2),
                }
            )
    return outliers


def write_report(report: dict, out_dir: Path) -> None:
    """Write calibration JSON/CSV and outlier JSON into the requested directory.

    Existing report files are overwritten. An empty calibration curve produces
    no new CSV file; the caller supplies the computed report dictionary.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "calibration_curve.json").write_text(json.dumps(report["calibration"], indent=2))
    if report["calibration"]["curve"]:
        cols = list(report["calibration"]["curve"][0].keys())
        with open(out_dir / "calibration_curve.csv", "w") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(report["calibration"]["curve"])
    (out_dir / "outlier_samples.json").write_text(json.dumps(report["outliers"], indent=2))
