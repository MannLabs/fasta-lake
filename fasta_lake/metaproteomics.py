"""Experimental metaproteomics transforms for declared measured-signal profiles.

These helpers do not infer taxa, resolve shared-protein origin or turn intensity
into biomass/flux. They accept explicit annotations and experimental design.


Array/table backends: https://github.com/numpy/numpy and
https://github.com/pandas-dev/pandas. Bray-Curtis uses SciPy pdist
(https://github.com/scipy/scipy). Optional JIT uses Numba
(https://github.com/numba/numba); profile accounting and permutation design
are FastaLake code, not copied AlphaPeptTools implementations.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

UNANNOTATED = "__unannotated__"


def _matrix(frame):
    """Validate a nonempty feature-by-acquisition frame and return numeric values.

    Identifiers must be unique and nonmissing. NaN remains missing; negative
    values and infinities are rejected.
    """
    if (
        not isinstance(frame, pd.DataFrame)
        or not frame.index.is_unique
        or not frame.columns.is_unique
    ):
        raise ValueError("Use a DataFrame with unique feature rows and acquisition columns")
    if frame.empty or frame.index.hasnans or frame.columns.hasnans:
        raise ValueError("The profile and its feature/acquisition IDs must be nonempty")
    values = frame.to_numpy(dtype=float)
    if np.isinf(values).any() or np.any(values[np.isfinite(values)] < 0):
        raise ValueError("Profile values must be nonnegative finite signal or missing")
    return values


def aggregate_terms(matrix, annotations):
    """Split observed group signal equally across distinct declared functional terms.

    matrix has reporting groups in rows and acquisitions in columns. annotations
    maps group IDs to term iterables at ONE level (e.g. KO). Missing/empty term
    lists retain signal in __unannotated__; unknown group keys fail. This is an
    explicit attribution convention, not evidence for equal enzyme contribution.
    Missing group cells contribute no observed signal; entirely empty acquisitions
    remain NaN. Return the term matrix and an acquisition-level coverage table.
    """
    values = _matrix(matrix)
    if not set(annotations) <= set(matrix.index):
        raise ValueError("Annotation contains groups absent from the study matrix")
    mapping = {}
    for group in matrix.index:
        terms = annotations.get(group, ())
        if isinstance(terms, str):
            raise ValueError("Supply a list/set of terms per group, not an unsplit string")
        terms = sorted(set(terms))
        if any(not isinstance(t, str) or not t.strip() or t == UNANNOTATED for t in terms):
            raise ValueError("Annotation terms must be nonempty strings without the reserved key")
        mapping[group] = terms or [UNANNOTATED]
    terms = sorted({term for labels in mapping.values() for term in labels})
    output = pd.DataFrame(0.0, index=terms, columns=matrix.columns)
    for group, signal in zip(matrix.index, np.nan_to_num(values, nan=0.0)):
        labels = mapping[group]
        output.loc[labels] += signal / len(labels)
    measured = np.isfinite(values) & (values > 0)
    totals = np.nansum(values, axis=0)
    empty = totals == 0
    output.loc[:, empty] = np.nan
    unannotated = (
        output.loc[UNANNOTATED].to_numpy()
        if UNANNOTATED in output.index
        else np.zeros(len(matrix.columns))
    )
    coverage = pd.DataFrame(
        {
            "observed_signal": totals,
            "measured_groups": measured.sum(axis=0),
            "unannotated_signal": unannotated,
            "annotation_fraction": np.divide(
                totals - unannotated, totals, out=np.full(len(totals), np.nan), where=~empty
            ),
        },
        index=matrix.columns,
    )
    if not np.allclose(output.sum(axis=0), totals, rtol=1e-12, atol=1e-12):
        raise ValueError("Functional aggregation failed signal conservation")
    return output, coverage


def relative_signal(profile):
    """Close observed signal to one per nonempty acquisition; preserve empty columns.

    Unannotated signal remains in the denominator when supplied in the profile.
    These fractions describe measured signal, not absolute abundance or biomass.
    """
    _matrix(profile)
    total = profile.sum(axis=0)
    return profile.div(total.where(total > 0), axis=1)


def diversity_summary(profile):
    """Return observed richness, Shannon entropy and effective Shannon diversity.

    The unit is the supplied feature (KO, EC, reporting group, or explicit taxon),
    not automatically species. An empty acquisition has undefined diversity.
    """
    fractions = relative_signal(profile).to_numpy()
    positive = fractions > 0
    logs = np.zeros_like(fractions)
    np.log(fractions, out=logs, where=positive)
    shannon = -np.nansum(fractions * logs, axis=0)
    result = pd.DataFrame(
        {
            "observed_richness": positive.sum(axis=0).astype(float),
            "shannon": shannon,
            "effective_shannon": np.exp(shannon),
            "simpson_concentration": np.nansum(fractions**2, axis=0),
        },
        index=profile.columns,
    )
    result.loc[profile.sum(axis=0) <= 0] = np.nan
    return result


def bray_curtis(profile, *, normalize=True):
    """Return a labelled Bray–Curtis distance matrix and explicitly excluded samples.

    Default normalization compares composition; normalize=False compares observed
    signal magnitudes as well. Unmeasured feature cells contribute zero observed
    signal. Entirely empty acquisitions are excluded rather than assigned distances.
    """
    from scipy.spatial.distance import pdist, squareform

    _matrix(profile)
    included = profile.sum(axis=0) > 0
    selected = profile.loc[:, included]
    if selected.shape[1] < 2:
        raise ValueError("Bray–Curtis needs at least two nonempty acquisitions")
    selected = relative_signal(selected) if normalize else selected
    distances = squareform(pdist(selected.fillna(0).T.to_numpy(), metric="braycurtis"))
    return pd.DataFrame(distances, index=selected.columns, columns=selected.columns), list(
        profile.columns[~included]
    )


def clr_complete(profile):
    """Return sample × feature natural-log CLR on the complete positive core.

    No pseudocount or imputation is applied. All acquisitions must have some
    observed signal. Return transformed data and the excluded feature roster;
    this core can be much narrower than the full metaproteome.
    """
    _matrix(profile)
    if (profile.sum(axis=0) <= 0).any():
        raise ValueError("Declare nonempty acquisitions before CLR")
    complete = (profile > 0).all(axis=1) & profile.notna().all(axis=1)
    if complete.sum() < 2:
        raise ValueError("CLR needs at least two features positive in every selected acquisition")
    values = np.log(profile.loc[complete].T)
    return values.sub(values.mean(axis=1), axis=0), list(profile.index[~complete])


def _within_group_ss(squared, labels, counts):
    """Sum within-group squared distances, dividing each pair by its group size."""
    total = 0.0
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            if labels[i] == labels[j]:
                total += squared[i, j] / counts[labels[i]]
    return total


_compiled_within = None


def permanova(distances, labels, *, blocks=None, permutations=999, seed=0, backend="numpy"):
    """One-way experimental PERMANOVA with optional within-block label permutations.

    distances is a square DataFrame; labels and optional blocks are Series indexed
    by its exact acquisition roster. Unblocked permutations require independent,
    exchangeable units. Within-donor blocks test within-donor effects only: donor-
    constant conditions cannot be tested by shuffling within donors. Groups need
    at least two observations. This is an association test, not causality; assess
    dispersion differences separately. Numba accelerates only the numeric sum,
    with fastmath disabled and the same Python-generated seeded permutations.
    """
    global _compiled_within
    if backend not in {"numpy", "numba"}:
        raise ValueError("Choose numpy or numba backend")
    if type(permutations) is not int or permutations < 1:
        raise ValueError("Use a positive integer permutation count")
    if (
        not isinstance(distances, pd.DataFrame)
        or not distances.index.is_unique
        or not distances.index.equals(distances.columns)
    ):
        raise ValueError("Distance axes must be the same unique acquisition roster")
    d = distances.to_numpy(dtype=float)
    if (
        not np.isfinite(d).all()
        or np.any(d < 0)
        or not np.allclose(d, d.T, atol=1e-12)
        or not np.allclose(np.diag(d), 0, atol=1e-12)
    ):
        raise ValueError("Distances must be finite, nonnegative, symmetric, with zero diagonal")

    def align(series, name):
        """Align complete acquisition labels to the distance matrix and encode groups."""
        if (
            not isinstance(series, pd.Series)
            or not series.index.is_unique
            or set(series.index) != set(distances.index)
            or series.isna().any()
        ):
            raise ValueError(f"{name} must label exactly the distance-matrix acquisitions")
        return pd.factorize(series.reindex(distances.index), sort=True)[0]

    group = align(labels, "Group labels")
    counts = np.bincount(group)
    n, k = len(group), len(counts)
    if k < 2 or np.any(counts < 2):
        raise ValueError("PERMANOVA needs at least two groups with two observations each")
    block_indices = None
    if blocks is not None:
        block = align(blocks, "Block labels")
        block_indices = [np.flatnonzero(block == b) for b in np.unique(block)]
        if not any(len(np.unique(group[ix])) > 1 for ix in block_indices):
            raise ValueError("No group labels can move within the declared blocks")
    squared = d * d
    total = math.fsum(squared[np.triu_indices(n, 1)]) / n
    if total <= 0:
        raise ValueError("PERMANOVA requires nonzero total dissimilarity")
    if backend == "numba":
        if _compiled_within is None:
            from numba import njit

            _compiled_within = njit(cache=False, fastmath=False)(_within_group_ss)

        def within(values):
            """Compute within-group dispersion for one observed or permuted label vector."""
            return _compiled_within(squared, values, counts)
    else:

        def within(values):
            """Compute within-group dispersion for one observed or permuted label vector."""
            return sum(
                squared[np.ix_(values == g, values == g)].sum() / (2 * counts[g]) for g in range(k)
            )

    def statistic(values):
        """Compute the PERMANOVA pseudo-F statistic for the supplied group labels."""
        w = within(values)
        return ((total - w) / (k - 1)) / (w / (n - k)) if w > 0 else math.inf

    observed = statistic(group)
    if not math.isfinite(observed):
        raise ValueError("PERMANOVA is undefined here because within-group dispersion is zero")
    generator = np.random.default_rng(seed)
    null = np.empty(permutations)
    for i in range(permutations):
        permuted = group.copy()
        if block_indices is None:
            permuted = generator.permutation(group)
        else:
            for ix in block_indices:
                permuted[ix] = generator.permutation(group[ix])
        null[i] = statistic(permuted)
    tolerance = 1e-12 * max(1.0, abs(observed))
    return {
        "pseudo_f": observed,
        "r_squared": (total - within(group)) / total,
        "p_value": (1 + int(np.sum(null >= observed - tolerance))) / (permutations + 1),
        "permutations": permutations,
        "seed": seed,
        "blocked": blocks is not None,
        "samples": n,
        "groups": k,
        "backend": backend,
        "comparison_tolerance": tolerance,
        "null_statistics": null,
    }
