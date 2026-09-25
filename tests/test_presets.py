"""Tests for v5 parameter presets."""

import pytest

from fasta_lake.presets import (
    DEFAULT,
    LENIENT,
    PRESETS,
    STRINGENT,
    apply_preset,
    describe_preset,
    get_preset,
    preset_sage_config,
)
from fasta_lake.run_config import RunConfig


class TestPresets:
    def test_all_presets_defined(self):
        assert set(PRESETS.keys()) == {"stringent", "default", "lenient", "proteomics"}

    def test_get_preset(self):
        assert get_preset("default").name == "default"

    def test_unknown_preset_raises(self):
        with pytest.raises(ValueError, match="Unknown preset"):
            get_preset("turbo")

    def test_stringent_tighter_than_default(self):
        # Higher min_length, lower q_threshold
        assert STRINGENT.search_min_length >= DEFAULT.search_min_length
        assert STRINGENT.search_fdr <= DEFAULT.search_fdr
        assert STRINGENT.search_precursor_ppm <= DEFAULT.search_precursor_ppm

    def test_lenient_looser_than_default(self):
        # Lower min_length, higher q_threshold, kmer rescue on
        assert LENIENT.search_min_length <= DEFAULT.search_min_length
        assert LENIENT.search_fdr >= DEFAULT.search_fdr
        assert LENIENT.kmer_rescue is True

    def test_default_min_len_is_9(self):
        assert DEFAULT.denovo_min_length == 9
        assert DEFAULT.search_min_length == 9

    def test_apply_preset_to_runconfig(self):
        cfg = RunConfig()
        apply_preset(cfg, "stringent")
        assert cfg.search_min_length == 10
        assert cfg.search_fdr == 0.001
        assert cfg.kmer_rescue is False

        apply_preset(cfg, "lenient")
        assert cfg.search_min_length == 7
        assert cfg.search_fdr == 0.05
        assert cfg.kmer_rescue is True

    def test_preset_sage_config(self):
        base = {
            "database": {
                "enzyme": {"min_len": 7, "max_len": 50, "missed_cleavages": 2},
                "max_variable_mods": 2,
            },
            "precursor_tol": {"ppm": [-10, 10]},
        }
        out = preset_sage_config(base, "stringent")
        assert out["database"]["enzyme"]["min_len"] == 10
        assert out["database"]["enzyme"]["missed_cleavages"] == 1
        assert out["database"]["max_variable_mods"] == 1
        assert out["precursor_tol"]["ppm"] == [-5.0, 5.0]
        # Base unmodified
        assert base["database"]["enzyme"]["min_len"] == 7

    def test_describe_preset_runs(self):
        for name in PRESETS:
            desc = describe_preset(name)
            assert name.upper() in desc
            assert "Stage 1-3" in desc
            assert "Stage 4" in desc
