"""The declared Python floor, the classifiers and the CI matrix agree."""

import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_python_floor_matches_classifiers_and_ci():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    assert project["requires-python"] == ">=3.11"
    versions = {
        c.rsplit(" :: ", 1)[1]
        for c in project["classifiers"]
        if re.fullmatch(r"Programming Language :: Python :: 3\.\d+", c)
    }
    assert min(versions, key=lambda v: int(v.split(".")[1])) == "3.11"
    ci = (ROOT / ".github/workflows/ci.yml").read_text()
    matrices = re.findall(r"python: \[([^\]]+)\]", ci)
    assert matrices
    for matrix in matrices:
        tested = [v.strip(" '\"") for v in matrix.split(",")]
        assert all(int(v.split(".")[1]) >= 11 for v in tested), matrix
