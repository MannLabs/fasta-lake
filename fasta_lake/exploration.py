"""Optional experimental functional profiles from a fixed study matrix.

Aggregation needs an explicit group-to-term map at one declared level. It cannot
resolve the origin of a shared peptide. Unannotated signal remains visible in
coverage and signal tables, and is excluded from diversity and ordination.


PCA backend: https://github.com/MannLabs/alphapepttools (0.4.0).
FastaLake implements the explicit profile accounting and design checks around
that backend; optional numerical transforms are in metaproteomics.py.
"""

from __future__ import annotations

import csv
import importlib.metadata
import json
from pathlib import Path


def read_group_terms(path):
    """Read a long table with study_group and term columns; deduplicate exact pairs."""
    result = {}
    with Path(path).open() as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        if rows.fieldnames != ["study_group", "term"]:
            raise ValueError("Annotation table needs exactly study_group and term columns")
        for row in rows:
            if None in row or not row["study_group"] or row["term"] is None:
                raise ValueError("Malformed study-group annotation row")
            terms = result.setdefault(row["study_group"], set())
            if row["term"].strip():
                terms.add(row["term"].strip())
    return result


def explore_study(
    matrix_path,
    annotation_path,
    output,
    *,
    level,
    metadata=None,
    group_column=None,
    block_column=None,
    permutations=999,
    seed=0,
    backend="numpy",
):
    """Export functional signal, coverage, diversity, Bray–Curtis and CLR PCA.

    matrix_path is the sum or directLFQ study_group_matrix.tsv chosen by the
    caller. Its values are never requantified. Diversity and composition refer
    only to annotated signal; coverage reports how much of the total it represents.
    Metadata must contain exactly the matrix acquisitions in a sample column.
    A one-way PERMANOVA runs only when group_column is explicitly supplied.
    block_column restricts permutations to repeated-measure units; the caller
    remains responsible for exchangeability and dispersion interpretation.
    """
    import anndata as ad
    import numpy as np
    import pandas as pd

    from fasta_lake.downstream import require_analysis
    from fasta_lake.metaproteomics import (
        UNANNOTATED,
        aggregate_terms,
        bray_curtis,
        clr_complete,
        diversity_summary,
        permanova,
        relative_signal,
    )
    from fasta_lake.multi_omics.exact import sha256
    from fasta_lake.quantification import _baseline

    require_analysis()
    import alphapepttools as apt

    if not isinstance(level, str) or not level.strip():
        raise ValueError("Declare one annotation level, such as KO, EC or genus")
    if block_column and not group_column:
        raise ValueError("A block column requires an explicit group column")
    if group_column and metadata is None:
        raise ValueError("Group comparisons require sample metadata")
    out = Path(output)
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"Exploration output exists: {out}; choose a new directory")
    paths = [Path(matrix_path).resolve(strict=True), Path(annotation_path).resolve(strict=True)]
    if metadata:
        paths.append(Path(metadata).resolve(strict=True))
    identities = {str(p): sha256(p) for p in paths}
    _baseline(paths[0])
    matrix = pd.read_csv(
        paths[0],
        sep="\t",
        index_col=0,
        dtype={"study_group": str},
        keep_default_na=False,
        na_values=[""],
    )
    profile, coverage = aggregate_terms(matrix, read_group_terms(paths[1]))
    if metadata:
        with Path(metadata).open() as stream:
            header = next(csv.reader(stream, delimiter="\t"), [])
        if not header or header[0] != "sample" or len(header) != len(set(header)):
            raise ValueError("Metadata needs a unique sample first column and distinct columns")
        obs = pd.read_csv(metadata, sep="\t", index_col=0, dtype=str, keep_default_na=False)
        if not obs.index.is_unique or set(obs.index) != set(matrix.columns):
            raise ValueError("Metadata must contain exactly the matrix acquisitions")
        obs = obs.reindex(matrix.columns)
        for name in (group_column, block_column):
            if name and (name not in obs or obs[name].eq("").any()):
                raise ValueError(f"Metadata column is absent or incomplete: {name}")
    else:
        obs = pd.DataFrame(index=matrix.columns)
    fractions = relative_signal(profile)
    annotated = profile.drop(index=UNANNOTATED, errors="ignore")
    included = annotated.sum(axis=0) > 0
    ordination = annotated.loc[:, included]
    status = {"pca": "EXCLUDED: fewer than three annotated acquisitions", "bray_curtis": "EXCLUDED"}
    distances, clr, scores, loadings, variance = None, None, None, None, None
    excluded_features = list(annotated.index)
    diversity = None
    if not annotated.empty:
        diversity = diversity_summary(annotated)
    if included.sum() >= 2:
        distances, _ = bray_curtis(ordination)
        status["bray_curtis"] = "PASS"
    if included.sum() >= 3:
        complete = (ordination > 0).all(axis=1) & ordination.notna().all(axis=1)
        if complete.sum() >= 2:
            clr, excluded_features = clr_complete(ordination)
            varying = clr.var(axis=0, ddof=0) > 1e-15
            if varying.sum() >= 2:
                data = ad.AnnData(
                    clr.loc[:, varying].to_numpy(),
                    obs=obs.loc[clr.index].copy(),
                    var=pd.DataFrame(index=clr.columns[varying]),
                )
                n_comps = min(3, data.n_obs - 1, data.n_vars)
                apt.tl.pca(
                    data,
                    n_comps=n_comps,
                    dim_space="obs",
                    svd_solver="full",
                    zero_center=True,
                    random_state=seed,
                )
                columns = [f"PC{i + 1}" for i in range(n_comps)]
                scores = pd.DataFrame(data.obsm["X_pca_obs"], index=data.obs_names, columns=columns)
                loadings = pd.DataFrame(
                    data.varm["PCs_pca_obs"], index=data.var_names, columns=columns
                )
                variance = pd.DataFrame(data.uns["variance_pca_obs"], index=columns)
                status["pca"] = "PASS"
            else:
                status["pca"] = "EXCLUDED: fewer than two varying complete CLR features"
        else:
            status["pca"] = "EXCLUDED: fewer than two complete positive annotated features"
    test = None
    if group_column:
        if distances is None:
            raise ValueError("PERMANOVA requires at least two annotated acquisitions")
        selected_metadata = obs.loc[distances.index]
        test = permanova(
            distances,
            selected_metadata[group_column],
            blocks=selected_metadata[block_column] if block_column else None,
            permutations=permutations,
            seed=seed,
            backend=backend,
        )
    provenance = {
        "status": "PASS",
        "experimental": True,
        "annotation_level": level,
        "inputs": identities,
        "implementation_sha256": sha256(__file__),
        "versions": {
            p: importlib.metadata.version(p)
            for p in ("alphapepttools", "anndata", "numpy", "pandas", "scipy")
        },
        "attribution": "Equal split across distinct supplied terms per reporting group",
        "coverage_denominator": "All measured study-group signal, including unannotated",
        "diversity_and_distance_denominator": "Annotated signal only; no unannotated pseudo-taxon",
        "clr": "Complete positive annotated core; natural log; no pseudocount or imputation",
        "pca": "Feature-centred CLR; no unit-variance scaling; AlphaPeptTools",
        "excluded_acquisitions": list(matrix.columns[~included]),
        "clr_excluded_features": excluded_features,
        "outputs": status,
        "group_column": group_column,
        "block_column": block_column,
    }
    out.mkdir(parents=True)
    for name, frame, index in [
        ("term_signal", profile, "term"),
        ("relative_signal", fractions, "term"),
        ("annotation_coverage", coverage, "sample"),
        ("diversity", diversity, "sample"),
        ("bray_curtis", distances, "sample"),
        ("clr_values", clr, "sample"),
        ("pca_scores", scores, "sample"),
        ("pca_loadings", loadings, "term"),
        ("pca_variance", variance, "component"),
    ]:
        if frame is not None:
            frame.to_csv(out / (name + ".tsv"), sep="\t", index_label=index)
    annotated_data = ad.AnnData(
        profile.T.to_numpy(), obs=obs.copy(), var=pd.DataFrame(index=profile.index)
    )
    annotated_data.uns["fastalake"] = json.dumps(provenance, sort_keys=True)
    annotated_data.write_h5ad(out / "functional_profile.h5ad", compression="gzip")
    if test:
        null = test.pop("null_statistics")
        pd.DataFrame({"permutation": np.arange(1, len(null) + 1), "pseudo_f": null}).to_csv(
            out / "permutation_statistics.tsv", sep="\t", index=False
        )
        (out / "PERMANOVA.json").write_text(json.dumps(test, indent=2) + "\n")
    _plot(coverage, diversity, scores, out, level, variance)
    if {str(p): sha256(p) for p in paths} != identities:
        raise ValueError("An input changed during functional exploration")
    provenance["files"] = {p.name: sha256(p) for p in sorted(out.iterdir()) if p.is_file()}
    (out / "COMPLETE.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return provenance


def _plot(coverage, diversity, scores, output, level, variance):
    """Render coverage beside annotated-feature diversity; label experimental scope."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    axes[0].bar(range(len(coverage)), coverage.annotation_fraction * 100, color="#2878a8")
    axes[0].set(ylabel="Annotated measured signal (%)", ylim=(0, 100), title="Annotation coverage")
    if diversity is not None:
        axes[1].bar(range(len(diversity)), diversity.effective_shannon, color="#dc9342")
    axes[1].set(ylabel=f"Effective Shannon diversity ({level})", title="Annotated features only")
    for ax in axes:
        ax.set_xticks(range(len(coverage)), coverage.index, rotation=90, fontsize=7)
    fig.suptitle("Experimental functional exploration")
    for ext in ("png", "pdf"):
        fig.savefig(output / ("Functional_QC." + ext), dpi=180)
    plt.close(fig)
    if scores is not None:
        fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
        ax.scatter(scores.PC1, scores.PC2, color="#2878a8")
        ax.margins(x=0.18, y=0.12)
        for sample, row in scores.iterrows():
            ax.annotate(
                sample, (row.PC1, row.PC2), fontsize=7, xytext=(3, 3), textcoords="offset points"
            )
        ax.set(
            xlabel=f"PC1 ({100 * variance.loc['PC1', 'variance_ratio']:.1f}%)",
            ylabel=f"PC2 ({100 * variance.loc['PC2', 'variance_ratio']:.1f}%)",
            title=f"Experimental {level} CLR PCA\nComplete positive core",
        )
        for ext in ("png", "pdf"):
            fig.savefig(output / ("PCA_clr." + ext), dpi=180)
        plt.close(fig)
