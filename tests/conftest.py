"""Shared test fixtures for FASTA Lake."""

from __future__ import annotations

import warnings
from pathlib import Path

import pytest

# Test modules import alphapepttools at collection time (pytest.importorskip),
# where pytest records import-time warnings without applying the configured
# filters. Import it once here, with the two known third-party import-time
# messages silenced (mudata's deprecated decorator arguments, alphabase's
# optional progress bar), so the cached module raises nothing later. Every
# other warning is an error (pyproject filterwarnings).
_KNOWN_IMPORT_WARNINGS = (
    "The decorator_name argument is deprecated",  # mudata/__init__.py
    "The docstring_style argument is deprecated",  # mudata/__init__.py
    "Dependency 'progressbar' not installed",  # alphabase data downloader
)
with warnings.catch_warnings(record=True) as _import_warnings:
    # Record rather than filter: mudata installs its own warning filter at
    # import time, which defeats any filter set before it. Recording captures
    # the warnings here so pytest's collection recorder never sees them.
    warnings.simplefilter("always")
    try:
        import alphapepttools  # noqa: F401
    except Exception:  # noqa: BLE001 - a broken optional backend must not abort the whole suite
        pass  # the modules that need it use pytest.importorskip and report it there
for _w in _import_warnings:
    if not str(_w.message).startswith(_KNOWN_IMPORT_WARNINGS):
        # Anything new must stay audible: re-emit it for the normal policy.
        warnings.warn_explicit(_w.message, _w.category, _w.filename, _w.lineno)

# Paths
TESTS_DIR = Path(__file__).parent
DATA_DIR = TESTS_DIR / "data"
REPO_ROOT = TESTS_DIR.parent
E2E_DIR = REPO_ROOT / "e2e_test"


@pytest.fixture
def tiny_fasta() -> Path:
    """Path to the tiny test FASTA with 20 proteins."""
    return DATA_DIR / "tiny.fasta"


@pytest.fixture
def tiny_predictions() -> Path:
    """Path to the tiny test predictions CSV with 15 peptides."""
    return DATA_DIR / "tiny_predictions.csv"


@pytest.fixture
def tiny_taxonomy() -> Path:
    """Taxonomy map for the tiny FASTA ids that carry no species (assembly, smORF)."""
    return DATA_DIR / "tiny_taxonomy.tsv"


@pytest.fixture
def tmp_output(tmp_path: Path) -> Path:
    """Temporary output directory for test results."""
    out = tmp_path / "output"
    out.mkdir()
    return out


@pytest.fixture
def sample_config_dict() -> dict:
    """Minimal valid config dictionary."""
    return {
        "database": {
            "sources": [{"tag": "REF", "path": str(DATA_DIR / "tiny.fasta")}],
        },
        "peptide_evidence": {
            "predictions_dir": str(DATA_DIR),
        },
        "output": {
            "directory": "/tmp/fastalake_test",
        },
    }


@pytest.fixture
def full_config_dict(sample_config_dict: dict) -> dict:
    """Full config dictionary with all sections."""
    return {
        **sample_config_dict,
        "inference": {"strategy": "species_budget", "min_peptides": 1},
        "clustering": {"enabled": True, "identity": 0.97, "coverage": 0.8},
        "rust": {"fasta_extractor": None, "lake_builder": None, "threads": None},
    }
