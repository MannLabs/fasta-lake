"""The generated reference pages must match the code they describe."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TOOL = REPO / "tools" / "render_reference.py"


def test_generated_reference_is_current():
    r = subprocess.run(
        [sys.executable, str(TOOL), "--check"], capture_output=True, text=True, cwd=REPO
    )
    assert r.returncode == 0, (
        "docs/reference/generated/ is stale relative to the code.\n"
        "Regenerate with: python tools/render_reference.py --write\n" + r.stdout + r.stderr
    )
