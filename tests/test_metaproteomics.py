"""Experimental downstream helpers preserve units, design and signal accounting."""

import itertools

import numpy as np
import pandas as pd
import pytest

from fasta_lake.metaproteomics import (
    UNANNOTATED,
    aggregate_terms,
    bray_curtis,
    clr_complete,
    diversity_summary,
    permanova,
    relative_signal,
)

scipy_distance = pytest.importorskip("scipy.spatial.distance")
pdist, squareform = scipy_distance.pdist, scipy_distance.squareform


def test_multifunction_signal_is_split_once_and_unknown_signal_is_retained():
    matrix = pd.DataFrame(
        [[12, 6, np.nan], [3, np.nan, np.nan], [5, 4, np.nan]],
        index=["a", "b", "c"],
        columns=["s1", "s2", "empty"],
    )
    profile, coverage = aggregate_terms(matrix, {"a": ["K1", "K2", "K1"], "b": ["K1"]})
    np.testing.assert_array_equal(profile["s1"], [9, 6, 5])
    np.testing.assert_array_equal(profile["s2"], [3, 3, 4])
    assert profile["empty"].isna().all()
    assert coverage.loc["s1", "annotation_fraction"] == 0.75
    assert coverage.loc["s2", "annotation_fraction"] == 0.6
    assert np.isnan(coverage.loc["empty", "annotation_fraction"])
    np.testing.assert_allclose(relative_signal(profile)["s1"], [0.45, 0.30, 0.25])


@pytest.mark.parametrize(
    "annotations", [{"foreign": ["K1"]}, {"a": "K1"}, {"a": [""]}, {"a": [UNANNOTATED]}]
)
def test_ambiguous_or_foreign_annotation_input_fails(annotations):
    with pytest.raises(ValueError):
        aggregate_terms(pd.DataFrame([[1]], index=["a"], columns=["s"]), annotations)


def test_diversity_is_feature_based_and_empty_samples_stay_undefined():
    profile = pd.DataFrame([[1, 4, np.nan], [1, 0, np.nan]], columns=["equal", "single", "empty"])
    result = diversity_summary(profile)
    assert result.loc["equal", "observed_richness"] == 2
    assert result.loc["equal", "shannon"] == pytest.approx(np.log(2))
    assert result.loc["equal", "effective_shannon"] == pytest.approx(2)
    assert result.loc["single", "shannon"] == 0
    assert result.loc["empty"].isna().all()


def test_bray_curtis_distinguishes_composition_from_signal_and_excludes_empty():
    profile = pd.DataFrame([[1, 10, np.nan], [3, 30, np.nan]], columns=["a", "b", "empty"])
    normalized, excluded = bray_curtis(profile)
    assert excluded == ["empty"] and normalized.loc["a", "b"] == 0
    raw, _ = bray_curtis(profile, normalize=False)
    assert raw.loc["a", "b"] == pytest.approx(36 / 44)


def test_clr_has_no_hidden_pseudocount_and_is_scale_invariant():
    profile = pd.DataFrame(
        [[1, 10], [3, 30], [0, 2], [4, np.nan]],
        index=["a", "b", "zero", "missing"],
        columns=["s1", "s2"],
    )
    clr, excluded = clr_complete(profile)
    assert excluded == ["zero", "missing"]
    np.testing.assert_allclose(clr.mean(axis=1), 0, atol=1e-15)
    np.testing.assert_allclose(clr.loc["s1"], clr.loc["s2"])


@pytest.mark.parametrize("values", [[[-1, 2]], [[np.inf, 2]], [[1, -np.inf]]])
def test_invalid_signal_is_rejected(values):
    with pytest.raises(ValueError):
        relative_signal(pd.DataFrame(values))


def distance_fixture():
    samples = list("abcdef")
    x = np.array([0, 1, 2, 9, 10, 12.0])
    distance = pd.DataFrame(squareform(pdist(x[:, None])), index=samples, columns=samples)
    labels = pd.Series(["A"] * 3 + ["B"] * 3, index=samples)
    return x, distance, labels


def test_permanova_matches_independent_univariate_anova_and_seeded_permutations():
    x, distances, labels = distance_fixture()
    result = permanova(distances, labels.iloc[::-1], permutations=49, seed=44)
    total = np.sum((x - x.mean()) ** 2)
    within = sum(np.sum((v - v.mean()) ** 2) for v in (x[:3], x[3:]))
    expected = (total - within) / (within / 4)
    assert result["pseudo_f"] == pytest.approx(expected)
    assert result["r_squared"] == pytest.approx((total - within) / total)
    generator = np.random.default_rng(44)
    exceed = 0
    for _ in range(49):
        permuted = generator.permutation(labels.to_numpy())
        w = sum(np.sum((x[permuted == g] - x[permuted == g].mean()) ** 2) for g in ("A", "B"))
        exceed += ((total - w) / (w / 4)) >= expected - result["comparison_tolerance"]
    assert result["p_value"] == (1 + exceed) / 50


def test_blocked_null_contains_only_within_block_exchanges():
    _, distances, labels = distance_fixture()
    blocks = pd.Series([1, 2, 3, 1, 2, 3], index=labels.index)
    result = permanova(distances, labels, blocks=blocks, permutations=63, seed=2)
    allowed = []
    for switches in itertools.product([False, True], repeat=3):
        permuted = labels.copy()
        for i, swap in enumerate(switches):
            if swap:
                permuted.iloc[i], permuted.iloc[i + 3] = labels.iloc[i + 3], labels.iloc[i]
        allowed.append(permanova(distances, permuted, permutations=1)["pseudo_f"])
    assert all(
        any(np.isclose(value, candidate) for candidate in allowed)
        for value in result["null_statistics"]
    )


def test_numba_preserves_the_entire_seeded_null_and_result():
    pytest.importorskip("numba")
    _, distances, labels = distance_fixture()
    blocks = pd.Series([1, 2, 3, 1, 2, 3], index=labels.index)
    reference = permanova(distances, labels, blocks=blocks, permutations=99, seed=19)
    compiled = permanova(
        distances, labels, blocks=blocks, permutations=99, seed=19, backend="numba"
    )
    np.testing.assert_allclose(
        reference["null_statistics"], compiled["null_statistics"], rtol=1e-12
    )
    assert compiled["p_value"] == reference["p_value"]
    assert compiled["r_squared"] == pytest.approx(reference["r_squared"], abs=1e-12)


def test_nonexchangeable_blocks_and_foreign_metadata_are_rejected():
    _, distances, labels = distance_fixture()
    with pytest.raises(ValueError, match="No group labels can move"):
        permanova(distances, labels, blocks=labels)
    with pytest.raises(ValueError, match="exactly"):
        permanova(distances, labels.rename(index={"a": "foreign"}))


def test_zero_dispersion_is_not_a_misleading_p_value():
    _, distances, labels = distance_fixture()
    with pytest.raises(ValueError, match="nonzero"):
        permanova(distances * 0, labels)
