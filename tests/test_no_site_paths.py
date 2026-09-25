"""Site-specific configuration must not be committed.

`tools/check_no_site_paths.py` has existed since 2026-08-02 and works. It was wired into
nothing, so it gated nothing, and the debt grew to 1,126 tracked files while a passing
test suite said the repository was healthy.

This is the wiring. It runs the scanner's own selftest first -- a scanner that cannot fail
has not passed -- and then the scan.

EXPECTED TO BE RED UNTIL THE DEBT IS PAID
-----------------------------------------
On 2026-08-04 the scan reports **290 files**, down from 1,126 once provenance records are
exempted (see below). That red is CORRECT and is the target of step 2 in
publication/notes/PLAN_2026-08-05.md. It is marked xfail with `strict=False` so the suite
stays green while the count comes down, and so that the day it reaches zero the test starts
PASSING and says so, rather than being quietly deleted.

Do not silence this by widening RECORD_PATTERNS. Records are exempt because rewriting them
would falsify provenance; code is not exempt because code is supposed to be portable.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCANNER = REPO / "tools" / "check_no_site_paths.py"


@pytest.mark.skipif(not SCANNER.is_file(), reason="scanner not in this checkout")
def test_scanner_selftest_passes():
    """The scanner must be able to FAIL on a planted defect and PASS on clean input.

    This runs before the scan itself on purpose. A green scan from a blind scanner is
    worse than a red one from a working scanner, because it is indistinguishable from
    having no problem.
    """
    r = subprocess.run(
        [sys.executable, str(SCANNER), "--selftest"],
        capture_output=True,
        text=True,
        cwd=REPO,
    )
    assert r.returncode == 0, (
        f"scanner selftest FAILED — it cannot be trusted:\n{r.stdout}\n{r.stderr}"
    )
    assert "selftest: PASS" in r.stdout


@pytest.mark.skipif(not SCANNER.is_file(), reason="scanner not in this checkout")
def test_no_site_specific_paths_in_tracked_files():
    r = subprocess.run([sys.executable, str(SCANNER)], capture_output=True, text=True, cwd=REPO)
    assert r.returncode == 0, (
        "tracked files carry site-specific configuration.\n"
        "Fix by moving the value into fastalake.site.sh and reading it in the script — "
        "never by deleting the path.\n"
        "See docs/get-started/environment.md.\n"
        "Full list: python3 tools/check_no_site_paths.py --list\n"
        f"{r.stdout[-2000:]}"
    )


@pytest.mark.skipif(not SCANNER.is_file(), reason="scanner not in this checkout")
def test_site_config_mechanism_exists():
    """The scanner tells you to move values into site config. That must be real.

    The environment guide and site template used to be missing, so the scanner's
    own remedy was a dead link and the 290 files had nowhere to go.
    """
    for rel in ("examples/site/fastalake.site.example.sh", "docs/get-started/environment.md"):
        p = REPO / rel
        assert p.is_file(), f"{rel} is referenced by the scanner's fix instructions but is missing"
        assert p.stat().st_size > 500, f"{rel} exists but is a stub"


def test_real_site_config_is_not_committed():
    """fastalake.site.sh holds real site values and must never be tracked."""
    r = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "fastalake.site.sh"],
        capture_output=True,
        text=True,
    )
    assert not r.stdout.strip(), (
        "fastalake.site.sh is TRACKED. It holds this site's real paths and accounts; "
        "only examples/site/fastalake.site.example.sh belongs in git."
    )
