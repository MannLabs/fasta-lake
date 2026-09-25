"""Search options must preserve database and acquisition identity."""

import json
from pathlib import Path

import pytest

from fasta_lake.search import DIGESTIONS, SEARCH_MODES, sage_config


def test_original_default_search_is_unchanged():
    original = json.loads((Path(__file__).parent / "data/sage_default_release.json").read_text())
    assert sage_config("targets.fasta", "sample.mzML", "search") == original


@pytest.mark.parametrize("digestion", DIGESTIONS)
@pytest.mark.parametrize("search_mode", SEARCH_MODES)
def test_search_choices_preserve_targets_and_acquisition(digestion, search_mode):
    cfg = sage_config(
        "my targets.fasta",
        "my sample.mzML",
        "new output",
        digestion=digestion,
        search_mode=search_mode,
        min_length=9,
        max_length=30,
    )
    assert cfg["database"]["fasta"] == "my targets.fasta"
    assert cfg["database"]["generate_decoys"] is True
    assert cfg["database"]["decoy_tag"] == "rev_"
    assert cfg["mzml_paths"] == ["my sample.mzML"]
    assert cfg["output_directory"] == "new output"
    assert cfg["quant"]["lfq"] is True
    assert cfg["database"]["enzyme"]["min_len"] == 9
    assert cfg["database"]["enzyme"]["max_len"] == 30
    # No nested state can leak from one acquisition or mode to the next.
    cfg["database"]["enzyme"]["cleave_at"] = "DE"
    assert sage_config("a", "b", "c")["database"]["enzyme"]["cleave_at"] == "KR"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"digestion": "non-tryptic"},
        {"search_mode": "wide_window"},
        {"min_length": 0},
        {"min_length": 30, "max_length": 9},
        {"min_length": 5.5},
        {"min_length": True},
    ],
)
def test_invalid_settings_fail_before_execution(kwargs):
    with pytest.raises(ValueError):
        sage_config("a", "b", "c", **kwargs)
