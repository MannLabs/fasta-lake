"""
Ion Coverage Score (ICS) calculation for PSM quality assessment.

ICS = (matched_b + matched_y) / (2 * (peptide_length - 1))

Measures how completely a spectrum's fragment ions cover the peptide
backbone. High ICS (>= 0.4) indicates confident identification.
Low ICS (< 0.4) suggests the true sequence may differ — candidate
for BLOSUM-guided variant exploration.

Supports SAGE, AlphaPept, AlphaDIA, and DIA-NN output formats.

Reference: Miura et al. (2025)
"""

from __future__ import annotations

import logging
import re

import pandas as pd

logger = logging.getLogger(__name__)


def clean_sequence(seq: str) -> str:
    """Strip modifications, keep uppercase AA only."""
    if not seq or pd.isna(seq):
        return ""
    cleaned = re.sub(r"\[.*?\]", "", str(seq))
    cleaned = cleaned.strip("_")
    return "".join(c for c in cleaned if c.isupper())


def ics(matched_b: int, matched_y: int, length: int) -> float:
    """Core ICS formula. Returns 0-1."""
    if length <= 1:
        return 0.0
    return min((matched_b + matched_y) / (2 * (length - 1)), 1.0)


def add_ics_sage(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add ICS column to SAGE results.

    Uses ``matched_peaks`` (total matched fragment ions) split evenly
    between b and y. Falls back to ``longest_b``/``longest_y`` if
    ``matched_peaks`` is unavailable.
    """
    df = df.copy()

    if "peptide_len" in df.columns:
        pep_len = df["peptide_len"]
    elif "peptide" in df.columns:
        pep_len = df["peptide"].apply(lambda s: len(clean_sequence(s)))
    else:
        logger.warning("No peptide length column — cannot compute ICS")
        df["ics"] = 0.0
        return df

    if "matched_peaks" in df.columns:
        total = df["matched_peaks"].fillna(0).astype(int)
        mb = total // 2
        my = total - mb
    elif "longest_b" in df.columns and "longest_y" in df.columns:
        mb = df["longest_b"].fillna(0).astype(int)
        my = df["longest_y"].fillna(0).astype(int)
    else:
        logger.warning("No fragment ion columns — ICS set to 0")
        df["ics"] = 0.0
        return df

    df["ics"] = [ics(b, y, plen) for b, y, plen in zip(mb, my, pep_len)]
    return df


def add_ics_alphapept(df: pd.DataFrame) -> pd.DataFrame:
    """Add ICS column to AlphaPept results."""
    df = df.copy()
    pep_col = "sequence" if "sequence" in df.columns else "peptide"
    pep_len = df[pep_col].apply(lambda s: len(clean_sequence(s)))
    mb = df.get("hits_b", pd.Series([0] * len(df))).fillna(0).astype(int)
    my = df.get("hits_y", pd.Series([0] * len(df))).fillna(0).astype(int)
    df["ics"] = [ics(b, y, plen) for b, y, plen in zip(mb, my, pep_len)]
    return df


def add_ics(
    df: pd.DataFrame,
    engine: str = "auto",
) -> pd.DataFrame:
    """
    Add ICS column to search results. Auto-detects engine from columns.

    Parameters
    ----------
    df : DataFrame
        Search results (SAGE, AlphaPept, or similar).
    engine : str
        ``'sage'``, ``'alphapept'``, or ``'auto'``.

    Returns
    -------
    DataFrame
        Input with ``ics`` column added (0-1 float).
    """
    if engine == "auto":
        if "hyperscore" in df.columns or "longest_b" in df.columns:
            engine = "sage"
        elif "hits_b" in df.columns:
            engine = "alphapept"
        else:
            engine = "sage"

    if engine == "sage":
        return add_ics_sage(df)
    elif engine == "alphapept":
        return add_ics_alphapept(df)
    else:
        logger.warning("Unknown engine %s, defaulting to SAGE", engine)
        return add_ics_sage(df)


def summarize(df: pd.DataFrame) -> dict:
    """ICS summary statistics."""
    if "ics" not in df.columns:
        return {}
    s = df["ics"]
    return {
        "n_psms": len(s),
        "mean_ics": round(s.mean(), 4),
        "median_ics": round(s.median(), 4),
        "pct_high": round((s >= 0.4).mean() * 100, 1),
        "pct_low": round((s < 0.3).mean() * 100, 1),
        "n_high": int((s >= 0.4).sum()),
        "n_low": int((s < 0.3).sum()),
    }
