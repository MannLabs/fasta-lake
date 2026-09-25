"""Check accounting conventions used by the HSTOOL comparison."""

import csv
import importlib.util
from pathlib import Path

import numpy as np
import pytest

SCRIPT = (
    Path(__file__).resolve().parents[1] / "publication/analyses/"
    "simplified_grouping_2026-09-10/hstool_full_comparison.py"
)
spec = importlib.util.spec_from_file_location("hstool_comparison", SCRIPT)
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


def test_coverage_unions_overlap_and_repeated_occurrences():
    # ABCD covers 1..4 and 6..9; BCD overlaps those intervals; X is unsupported.
    assert comparison.coverage("ABCDZABCDZ", {"ABCD", "BCD", "X"}) == (80.0, 1)


def test_coverage_treats_il_as_equivalent_without_counting_unsupported_peptide():
    assert comparison.coverage("MILKAAA", {"MLLK", "WWWW"}) == (400 / 7, 1)
    with pytest.raises(ValueError, match="Empty peptide"):
        comparison.coverage("AAA", {""})
    with pytest.raises(ValueError, match="empty sequence"):
        comparison.coverage("", {"AAA"})


@pytest.mark.parametrize(
    "q,expected",
    [
        (0, True),
        (0.01, True),
        (0.011, False),
        (-0.1, False),
        (float("nan"), False),
        (float("inf"), False),
    ],
)
def test_q_gate_rejects_out_of_range_or_nonfinite_values(q, expected):
    assert comparison.passing(q) is expected


def test_undefined_correlations_are_missing():
    assert comparison.correlation(np.array([1, 2]), np.array([1, 2])) is None
    assert comparison.correlation(np.ones(3), np.array([1, 2, 3])) is None
    assert comparison.correlation(np.array([1, 2, 3]), np.array([3, 2, 1])) == -1.0


def test_member_annotation_key_preserves_accession_bars():
    assert comparison.metadata_key("DB|UHGG|P2;P1;P1") == "P1;P2"
    assert comparison.metadata_key("sp|P12345|NAME") == "sp|P12345|NAME"


def test_unresolved_species_is_not_a_leader_switch():
    names = ["a"]
    design = {"a": {"group": "TechnicalReplicates"}}
    functions = [
        dict(
            method="original",
            sample="a",
            ko="K1",
            leading_species="Species A",
            normalised_function_intensity=1.0,
        ),
        dict(
            method="simplified",
            sample="a",
            ko="K1",
            leading_species="",
            normalised_function_intensity=1.0,
        ),
    ]
    rows, profiles = comparison.functional_comparison(functions, names, design, {"K1"})
    assert rows[0]["resolution_status_changed"]
    assert not rows[0]["both_resolved"]
    assert rows[0]["changed_when_both_resolved"] is None
    assert profiles[0]["pearson_log2_KO_profile"] is None


def test_detection_counts_missing_without_inventing_zero_intensity(tmp_path):
    design = {tier: {"group": tier} for tier in comparison.TIERS}
    names = list(design)
    data = {"group_matrix": {m: {s: {"G1": 10.0} for s in names} for m in comparison.METHODS}}
    comparison.export_matrices_and_detection(data, names, design, {"G1", "G2"}, tmp_path)
    rows = list(comparison.read_rows(tmp_path / "detection_completeness.tsv"))
    assert all(int(r["groups_detected_in_all"]) == 1 for r in rows)
    with (tmp_path / "simplified_study_group_matrix.tsv").open() as f:
        matrix = list(csv.reader(f, delimiter="\t"))
    assert matrix[2] == ["G2", "", "", "", "", ""]


def test_matched_correlations_exclude_rows_present_only_in_one_method():
    matched_spec = importlib.util.spec_from_file_location(
        "matched_reproducibility", SCRIPT.with_name("matched_reproducibility.py")
    )
    matched = importlib.util.module_from_spec(matched_spec)
    matched_spec.loader.exec_module(matched)
    pd = matched.pd
    original = pd.DataFrame({"a": [1.0, 2.0, 4.0], "b": [1.0, 2.0, 4.0]}, index=["G1", "G2", "G3"])
    direct = pd.DataFrame(
        {"a": [1.0, 2.0, 4.0, 100.0], "b": [4.0, 2.0, 1.0, 100.0]}, index=["G1", "G2", "G3", "G4"]
    )
    result = matched.common_pair(original, direct, "a", "b")
    assert result == {
        "common_positive_rows": 3,
        "original": 1.0,
        "simplified": -1.0,
        "change": -2.0,
    }
