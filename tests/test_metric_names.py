"""A metric name must not lie about which stage it came from.

THE DEFECT THIS LOCKS OUT
-------------------------
`3_persample/*/sample_lake.stats.json` holds the per-sample LAKE -- every candidate a
sample's own de novo peptides support, BEFORE parsimony. `extract_tables.py` wrote it
into `funnel.tsv` as **`median_db_proteins`**, so panels 2a-6a labelled it "per-sample
database (median)" and PREDICT's panel A printed **1,865,847** where the database is
**8,147**. Figure 1a printed the right number for the same quantity in the same cohort:
a 229-fold contradiction inside one paper. On SIHUMIx the same mistake is 310-fold.

Nothing caught it, for two compounding reasons:

1. the name was plausible. `median_db_proteins` reads like the database, and every
   consumer that trusted the name inherited the error;
2. the self-check that should have caught it, `chk("per-sample median db", ...)`, was
   called with NO EXPECTED VALUE. A check with nothing to compare against cannot fail.

The two stages now have names that state their stage, and this test is what keeps it
that way -- it runs in the ordinary suite, so the defect cannot return quietly:

    stage3_persample  median_persample_lake_proteins   the search space BEFORE parsimony
    stage4_parsimony  median_persample_db_proteins     the database SAGE searched

Peter, 2026-08-04: "make sure this protein miscounting is fixed forever."
"""

from __future__ import annotations

from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
FIGURES = REPO / "publication" / "figures"

# The name is banned outright rather than deprecated. A wrong name kept beside a right
# one only guarantees something eventually reads the wrong one: the alias lived for a
# single afternoon, and in that time every figure's funnel.tsv was extracted carrying
# ONLY the wrong name and no stage-4 rung at all.
BANNED = "median_db_proteins"

HONEST = {
    "stage3": "median_persample_lake_proteins",
    "stage4": "median_persample_db_proteins",
}


def _figure_sources() -> list[Path]:
    if not FIGURES.is_dir():
        return []
    return [
        p
        for p in list(FIGURES.rglob("*.py")) + list(FIGURES.rglob("*.ipynb"))
        if "__pycache__" not in p.parts and ".ipynb_checkpoints" not in p.parts
    ]


@pytest.mark.skipif(not FIGURES.is_dir(), reason="publication/figures not in this checkout")
def test_banned_metric_name_absent_from_figure_code():
    """No figure code may read or write the stage-mislabelling name."""
    offenders = []
    for p in _figure_sources():
        try:
            text = p.read_text(errors="replace")
        except OSError:
            continue
        for i, line in enumerate(text.splitlines(), 1):
            # Flag USE, not mention. The name appears legitimately in prose that explains
            # why it was retired -- this file's own docstring does it, and so do the
            # comments in extract_tables.py and cohort_panel.py that record the 229-fold
            # contradiction. What must never reappear is the name as a STRING LITERAL,
            # because that is a lookup: `f.get(("stage3_persample", "median_db_proteins"))`.
            # Matching on quotes is precise, and it means documenting the defect stays
            # possible while committing it does not.
            # .ipynb stores source as JSON strings, so a lookup written in a cell as
            #     f.get(("stage3_persample", "median_db_proteins"))
            # appears on disk as \"median_db_proteins\" -- the escaped closing quote
            # defeats a naive `"name"` match, and a first version of this test passed the
            # notebooks clean while they still contained the lookup. Drop backslashes
            # before matching so the notebook and the .py are held to one rule.
            probe = line.replace("\\", "")
            if f'"{BANNED}"' in probe or f"'{BANNED}'" in probe:
                offenders.append(f"{p.relative_to(REPO)}:{i}: {line.strip()[:90]}")
    assert not offenders, (
        f"`{BANNED}` names the stage-3 LAKE while claiming to be the database. "
        "Use `median_persample_lake_proteins` (stage 3) or "
        "`median_persample_db_proteins` (stage 4).\n  " + "\n  ".join(offenders)
    )


@pytest.mark.skipif(not FIGURES.is_dir(), reason="publication/figures not in this checkout")
def test_extracted_funnels_use_honest_names():
    """Any funnel.tsv present must carry the honest names and not the banned one.

    Skips when no funnel has been extracted -- a fresh clone has none, and absence is
    not a failure. What is a failure is a funnel that exists and still lies.
    """
    funnels = sorted(FIGURES.glob("Fig*/tables/funnel.tsv"))
    if not funnels:
        pytest.skip("no funnel.tsv extracted in this checkout")
    bad, stale = [], []
    for f in funnels:
        text = f.read_text(errors="replace")
        if BANNED in text:
            bad.append(str(f.relative_to(REPO)))
        if HONEST["stage3"] not in text:
            stale.append(f"{f.relative_to(REPO)} (no {HONEST['stage3']})")
    assert not bad, (
        f"these funnel.tsv still carry `{BANNED}`, which names the lake as the database. "
        "Re-extract them:\n  " + "\n  ".join(bad)
    )
    assert not stale, (
        "these funnel.tsv predate the rename and cannot be read correctly. Re-extract:\n  "
        + "\n  ".join(stale)
    )


@pytest.mark.skipif(not FIGURES.is_dir(), reason="publication/figures not in this checkout")
def test_extractor_emits_both_stages():
    """The extractor must write stage 3 AND stage 4, or the confusion has nowhere to go.

    The original defect was not only a bad name: stage 4 was never extracted at all, so
    a figure wanting the database had nothing correct to read even if it asked.
    """
    src = (FIGURES / "extract_tables.py").read_text(errors="replace")
    for stage, metric in HONEST.items():
        assert metric in src, f"extract_tables.py no longer emits {stage}'s `{metric}`"
    assert f'"{BANNED}"' not in src and f"'{BANNED}'" not in src, (
        f"extract_tables.py emits `{BANNED}` again — the alias was removed deliberately"
    )
