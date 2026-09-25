"""Where the golden-reference data is, resolved in one place.

WHY THIS EXISTS. The suite reported 42 skips with the reason "Data not available",
which reads as "this machine does not have the data". For 40 of them that was false:
the data was on the machine and the tests were looking in the wrong place.

- `test_e2e.py` resolved `<repo>/projects/e2e_test`. Correct before the code moved
  into `FastaLake/canonical/`; afterwards the data sat one level up.
- `test_empirical.py` hard-coded `Peter/FastaLake/e2e_test/stage2` — missing the
  `projects/` component. That path has never existed, so those tests had *never run*
  while reporting a reason that blamed the machine.

A skip whose message blames the environment for a bug in the test is worse than a
failure, because nobody investigates it. So resolution happens once, here, and
`describe()` prints which candidate won.

Override with FL_JANKO_DIR when the data lives elsewhere.

The RoPE-era e2e/Masters lookups were removed with the tests that used them
(2026-08-03); see FastaLake/ARCHIVE_2026-08-03/rope_era_tests_deleted_2026-08-03/.
"""

from __future__ import annotations

import os
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]

# Optional study fixtures are supplied by the user, never discovered at a lab path.
_POOL_CANDIDATES = tuple(Path(p) for p in [os.environ.get("FASTALAKE_TEST_DATA")] if p)


def _first_existing(candidates, marker: str | None = None) -> Path | None:
    """Return the first candidate that exists (and contains `marker` if given)."""
    for c in candidates:
        if c is None:
            continue
        c = Path(c)
        if c.is_dir() and (marker is None or (c / marker).exists()):
            return c
    return None


def janko_root() -> Path | None:
    cands = [os.environ.get("FL_JANKO_DIR")]
    cands += [p / "PeTr_Janko" for p in _POOL_CANDIDATES]
    return _first_existing(cands, marker="alphanovo_output")


def describe() -> str:
    def one(name, p):
        return f"  {name:<20} {p if p else 'NOT FOUND — dependent tests will skip'}"

    return "\n".join(
        [
            "golden-reference data resolution:",
            one("janko_root", janko_root()),
        ]
    )


if __name__ == "__main__":
    print(describe())
