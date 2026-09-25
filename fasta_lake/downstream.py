"""Explicit QC and PCA of fixed study-group quantities with AlphaPeptTools.

Uses AlphaPeptTools 0.4.0: https://github.com/MannLabs/alphapepttools
(QC metrics, log transform and PCA). FastaLake supplies identity checks,
complete-case eligibility and reporting; it does not reimplement that backend.
"""

from __future__ import annotations

import importlib.metadata
import inspect
import json
import shutil
import tempfile
import warnings
from pathlib import Path

from fasta_lake.gzip_io import PANDAS_GZIP
from fasta_lake.quantification import _baseline, stamp

ALPHAPEPTTOOLS_VERSION = "0.4.0"


def require_analysis():
    """Load the pinned QC backend and return AlphaPeptTools, AnnData, NumPy and pandas.

    Raise RuntimeError before workflow output is created when the optional
    analysis environment is missing or its AlphaPeptTools version differs.
    """
    try:
        version = importlib.metadata.version("alphapepttools")
        if version != ALPHAPEPTTOOLS_VERSION:
            raise RuntimeError(f"Expected alphapepttools {ALPHAPEPTTOOLS_VERSION}; found {version}")
        with warnings.catch_warnings():
            # mudata, imported by alphapepttools, deprecates two arguments of its
            # own namespace-registration decorator at import time. Not ours to
            # fix and not acted on here; everything else stays audible.
            for message in ("The decorator_name argument", "The docstring_style argument"):
                warnings.filterwarnings("ignore", message, DeprecationWarning)
            import alphapepttools as apt
        import anndata as ad
        import numpy as np
        import pandas as pd
    except ImportError as error:
        raise RuntimeError(
            "Use Python >=3.11 and install: pip install fasta-lake[analysis]"
        ) from error
    return apt, ad, np, pd


def _names(path, expected, label):
    """Read unique nonblank IDs and require them to belong to the expected roster."""
    names = Path(path).read_text().splitlines()
    if not names or any(not n for n in names) or len(names) != len(set(names)):
        raise ValueError(f"{label} must contain one unique, nonblank identifier per line")
    if set(names) - set(expected):
        raise ValueError(
            f"Unknown identifiers in {label}: {sorted(set(names) - set(expected))[:5]}"
        )
    return names


def analyze_study(
    groups,
    output,
    *,
    quantification=None,
    metadata=None,
    replicate_column=None,
    pca_samples=None,
    pca_features=None,
):
    """Summarize study coverage and plot eligible PCA without imputing quantities.

    Parameters
    ----------
    groups : path-like
        Completed grouping bundle defining the acquisition and peptide roster.
    output : path-like
        New QC directory; existing output is refused.
    quantification : path-like or None
        Selected-method bundle. None uses the grouping bundle's sum matrix.
    metadata : path-like or None
        Acquisition metadata for joins and optional replicate summaries.
    replicate_column : str or None
        Metadata column defining the requested replicate groups.
    pca_samples, pca_features : path-like or None
        Optional files listing unique acquisition or feature IDs to include.

    Returns
    -------
    dict
        QC summary written with tables, plots and the analysis object.

    Notes
    -----
    PCA needs at least three nonempty acquisitions and two complete varying
    groups in the selected population. The report records eligibility and
    exclusions. Input quantities and missing cells are preserved. These are
    descriptive analyses, not independent biological replication.
    """
    apt, ad, np, pd = require_analysis()
    source = Path(groups).resolve()
    selected = Path(quantification).resolve() if quantification else source
    out = Path(output).absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"Output exists: {out}; choose a new directory")
    out = out.resolve()
    if out in source.parents or out in selected.parents:
        raise ValueError("Output must not contain an input directory")
    paths = [
        source / n for n in ["summary.json", "study_group_matrix.tsv", "peptide_to_study_group.tsv"]
    ]
    paths += [
        source / n
        for n in ["peptide_evidence.tsv.gz", "peptide_ambiguity.tsv"]
        if (source / n).exists()
    ]
    paths += [
        source / n
        for n in ["feature_assignments.tsv.gz", "acquisition_summary.tsv", "input_manifest.json"]
        if (source / n).exists()
    ]
    if selected != source:
        paths += [selected / "summary.json", selected / "study_group_matrix.tsv"]
    if selected != source:
        paths += [
            selected / n
            for n in ["peptide_support.tsv.gz", "ion_mapping.tsv"]
            if (selected / n).exists()
        ]
    paths += [Path(p).resolve() for p in [metadata, pca_samples, pca_features] if p]
    paths = list(dict.fromkeys(paths))
    inputs = [{"path": str(p), **stamp(p)} for p in paths]
    base_summary = json.loads(paths[0].read_text())
    if base_summary.get("schema_version") != 1 or not base_summary.get(
        "representative_sequence_compatibility_checked"
    ):
        raise ValueError("Expected a completed sequence-checked grouping bundle")
    samples, baseline_counts = _baseline(source / "study_group_matrix.tsv")
    if samples != base_summary.get("acquisition_names"):
        raise ValueError("Grouping summary and matrix rosters differ")
    base = pd.read_csv(
        source / "study_group_matrix.tsv", sep="\t", index_col=0, keep_default_na=False, dtype=str
    )
    method = "sum"
    if selected != source:
        record = json.loads((selected / "summary.json").read_text())
        method = record.get("method")
        if method not in {"sum", "directlfq"} or record.get("acquisition_names") != samples:
            raise ValueError("Invalid selected quantification summary or acquisition roster")
        for item in record["inputs"]:
            p = source / item["file"]
            if stamp(p) != item:
                raise ValueError(f"Selected quantification does not match source study: {p.name}")
        expected = next(i for i in record["outputs"] if i["file"] == "study_group_matrix.tsv")
        if stamp(selected / expected["file"]) != expected:
            raise ValueError("Selected quantification matrix checksum differs from its summary")
    matrix = selected / "study_group_matrix.tsv"
    names, _ = _baseline(matrix)
    frame = pd.read_csv(matrix, sep="\t", index_col=0, keep_default_na=False, dtype=str)
    if names != samples or not frame.index.equals(base.index):
        raise ValueError("Selected and baseline matrix axes must be identical")
    linear = frame.replace("", "0").astype(float).T
    linear.index.name, linear.columns.name = "sample", "study_group"
    positive = linear.to_numpy() > 0
    if np.any(positive & ~(base.replace("", "0").astype(float).T.to_numpy() > 0)):
        raise ValueError("Selected matrix has a new positive group/sample quantity")
    linear = linear.where(linear > 0)
    obs = pd.DataFrame(index=linear.index)
    if metadata:
        meta = pd.read_csv(metadata, sep="\t", dtype=str, keep_default_na=False)
        if "sample" not in meta or meta["sample"].duplicated().any():
            raise ValueError("Metadata needs a unique sample column")
        if set(meta["sample"]) != set(samples):
            raise ValueError("Metadata must contain exactly the matrix sample roster")
        obs = meta.set_index("sample").reindex(samples)
        reserved = {
            "total_sample_intensity",
            "num_features_detected",
            "fraction_detected_features",
            "baseline_sum_groups",
            "pca_included",
        }
        if reserved.intersection(obs.columns):
            raise ValueError("Metadata uses reserved QC column names")
    if replicate_column and (not metadata or replicate_column not in obs):
        raise ValueError("A replicate column requires the corresponding metadata column")
    data = ad.AnnData(X=linear.to_numpy(), obs=obs, var=pd.DataFrame(index=linear.columns))
    # Detection is evaluated on positive linear intensities, never on log2 values.
    with warnings.catch_warnings():
        if not data.n_vars:
            # No study groups at all is a legitimate empty study: AlphaPeptTools
            # averages over zero features, which numpy reports as "Mean of empty
            # slice" and "invalid value encountered in divide". The NaN metrics
            # that result are the correct statement for an empty study.
            warnings.filterwarnings("ignore", "Mean of empty slice", RuntimeWarning)
            warnings.filterwarnings("ignore", "invalid value encountered", RuntimeWarning)
        apt.metrics.calculate_qc_metrics(data)
    data.obs["baseline_sum_groups"] = baseline_counts
    if not np.array_equal(data.obs.num_features_detected, positive.sum(axis=1)):
        raise ValueError("AlphaPeptTools detection counts disagree with positive input cells")
    data.layers["log2"] = data.X.copy()
    apt.pp.nanlog(data, base=2, layer="log2", verbosity=0)
    if not np.array_equal(np.isfinite(data.layers["log2"]), positive):
        raise ValueError("Log transform changed the measured/missing mask")
    # Keep empty acquisitions in the exported object, with explicit PCA exclusion.
    requested = samples
    if pca_samples:
        requested = _names(pca_samples, samples, "PCA samples")
    sample_mask = data.obs_names.isin(requested) & (positive.sum(axis=1) > 0)
    log_values = data.layers["log2"][sample_mask]
    complete = (
        np.isfinite(log_values).all(axis=0) if len(log_values) else np.zeros(data.n_vars, bool)
    )
    eligible = np.zeros(data.n_vars, dtype=bool)
    if len(log_values) >= 2:
        eligible[complete] = np.ptp(log_values[:, complete], axis=0) > 0
    if pca_features:
        chosen = _names(pca_features, data.var_names, "PCA features")
        chosen_mask = data.var_names.isin(chosen)
        if np.any(chosen_mask & ~complete):
            raise ValueError("Requested PCA features must be measured in every included sample")
        eligible &= chosen_mask
    data.obs["pca_included"] = sample_mask
    data.var["pca_complete"] = complete
    data.var["pca_eligible"] = eligible
    pca_status = "SKIPPED: need at least three nonempty samples and two complete varying groups"
    if sample_mask.sum() >= 3 and eligible.sum() >= 2:
        subset = data[sample_mask, eligible].copy()
        n_comps = min(10, int(sample_mask.sum()) - 1, int(eligible.sum()))
        apt.tl.pca(
            subset,
            layer="log2",
            dim_space="obs",
            n_comps=n_comps,
            svd_solver="full",
            zero_center=True,
            random_state=0,
        )
        coordinates = np.full((data.n_obs, n_comps), np.nan)
        coordinates[sample_mask] = subset.obsm["X_pca_obs"]
        loadings = np.full((data.n_vars, n_comps), np.nan)
        loadings[eligible] = subset.varm["PCs_pca_obs"]
        data.obsm["X_pca"] = coordinates
        data.varm["PCs"] = loadings
        data.uns["pca_variance"] = subset.uns["variance_pca_obs"]
        pca_status = "PASS"
    # A sample/feature comparison is meaningful only with the exact mask recorded.
    cv = []
    if replicate_column:
        for label in sorted(set(data.obs[replicate_column]) - {""}):
            mask = data.obs[replicate_column].eq(label).to_numpy()
            sub = data[mask].copy()
            apt.metrics.coefficient_of_variation(sub, min_valid=2, key_added="cv_ddof0")
            n = np.isfinite(sub.X).sum(axis=0)
            population = sub.var["cv_ddof0"].to_numpy() * 100
            correction = np.sqrt(np.divide(n, n - 1, out=np.full(len(n), np.nan), where=n > 1))
            cv.append(
                pd.DataFrame(
                    {
                        "study_group": data.var_names,
                        "replicate_set": label,
                        "set_samples": int(mask.sum()),
                        "measured_samples": n,
                        "complete": n == mask.sum(),
                        "cv_population_pct": population,
                        "cv_sample_pct": population * correction,
                    }
                )
            )
    cv = (
        pd.concat(cv, ignore_index=True)
        if cv
        else pd.DataFrame(
            columns=[
                "study_group",
                "replicate_set",
                "set_samples",
                "measured_samples",
                "complete",
                "cv_population_pct",
                "cv_sample_pct",
            ]
        )
    )
    versions = {
        p: importlib.metadata.version(p)
        for p in [
            "alphapepttools",
            "anndata",
            "numpy",
            "pandas",
            "scanpy",
            "scikit-learn",
            "matplotlib",
        ]
    }
    provenance = dict(
        method=method,
        inputs=inputs,
        versions=versions,
        implementation=stamp(Path(__file__)),
        intensity_plot_implementation=stamp(Path(__file__).with_name("qc_plots.py")),
        backend_implementation=[
            stamp(Path(p))
            for p in sorted(
                {
                    inspect.getsourcefile(f)
                    for f in [
                        apt.metrics.calculate_qc_metrics,
                        apt.metrics.coefficient_of_variation,
                        apt.pp.nanlog,
                        apt.tl.pca,
                    ]
                }
            )
        ],
        transformations="Positive linear X; log2 layer; no normalization or imputation",
        feature_identity="Study reporting group; representative does not prove unique origin",
        pca=dict(
            status=pca_status,
            samples=int(sample_mask.sum()),
            groups=int(eligible.sum()),
            feature_centered=True,
            unit_variance_scaled=False,
        ),
        cv=dict(
            replicate_column=replicate_column,
            min_valid=2,
            upstream_ddof=0,
            sample_cv_correction="sqrt(n/(n-1)); n is observed values",
        ),
        peptide_evidence="Original peptide/ion ledgers remain linked by checksums and paths",
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{out.name}-", dir=out.parent) as directory:
        temp = Path(directory)
        data.obs.to_csv(temp / "sample_qc.tsv", sep="\t")
        data.var.to_csv(temp / "group_qc.tsv", sep="\t")
        cv.to_csv(temp / "replicate_cv.tsv.gz", sep="\t", index=False, compression=PANDAS_GZIP)
        shutil.copyfile(source / "peptide_to_study_group.tsv", temp / "peptide_to_study_group.tsv")
        if "X_pca" in data.obsm:
            columns = [f"PC{i + 1}" for i in range(data.obsm["X_pca"].shape[1])]
            pd.DataFrame(data.obsm["X_pca"], index=data.obs_names, columns=columns).to_csv(
                temp / "pca_scores.tsv", sep="\t", index_label="sample"
            )
            pd.DataFrame(data.varm["PCs"], index=data.var_names, columns=columns).to_csv(
                temp / "pca_loadings.tsv", sep="\t", index_label="study_group"
            )
            pd.DataFrame(data.uns["pca_variance"], index=columns).to_csv(
                temp / "pca_variance.tsv", sep="\t", index_label="component"
            )
        _plot(data, temp, method)
        from fasta_lake.qc_plots import plot_intensity_distributions

        provenance["intensity_distributions"] = plot_intensity_distributions(data, temp, method)
        data.uns["fastalake"] = json.dumps(provenance, sort_keys=True)
        data.write_h5ad(temp / "study.h5ad", compression="gzip")
        if inputs != [{"path": str(p), **stamp(p)} for p in paths]:
            raise ValueError("An input changed during downstream analysis")
        provenance.update(
            status="PASS",
            samples=data.n_obs,
            groups=data.n_vars,
            positive_cells=int(positive.sum()),
            outputs=[stamp(p) for p in sorted(temp.iterdir())],
        )
        (temp / "summary.json").write_text(json.dumps(provenance, indent=2) + "\n")
        temp.rename(out)
    return provenance


def _plot(data, out, method):
    """Save descriptive QC panels and the eligible PCA result using local plot styles."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    # Optional plotting backends/notebooks can change global fonts. Keep the
    # saved report readable and self-contained regardless of those styles.
    with plt.rc_context(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    ):
        fig, axes = plt.subplots(1, 3, figsize=(11, 4.0), constrained_layout=True)
        x = np.arange(1, data.n_obs + 1)
        axes[0].bar(x, data.obs.num_features_detected, color="#009E73")
        axes[0].set(xlabel="Acquisition row (sample_qc.tsv)", ylabel="Quantified study groups")
        if data.n_obs <= 12 and all(len(str(name)) <= 18 for name in data.obs_names):
            axes[0].set_xticks(x, data.obs_names, rotation=45 if data.n_obs > 4 else 0)
            axes[0].set_xlabel("Acquisition")
        else:
            axes[0].xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        axes[1].hist(
            data.var.fraction_detected_samples, bins=np.linspace(0, 1, 21), color="#56B4E9"
        )
        axes[1].set(xlabel="Fraction of acquisitions with quantity", ylabel="Study groups")

        def draw_pca(axis):
            """Draw available finite PCA scores, or show why no PCA coordinates were saved."""
            if "X_pca" in data.obsm:
                xy = data.obsm["X_pca"]
                if data.n_obs <= 20:
                    handles = []
                    colors = plt.get_cmap("tab20").colors
                    short_names = all(len(str(name)) <= 18 for name in data.obs_names)
                    for i, row in enumerate(xy):
                        if np.isfinite(row).all():
                            handles.append(
                                axis.scatter(
                                    row[0],
                                    row[1],
                                    color=colors[i % 20],
                                    s=25,
                                    label=str(data.obs_names[i]) if short_names else str(i + 1),
                                )
                            )
                    axis.figure.legend(
                        handles=handles,
                        loc="outside lower center",
                        ncol=min(6, len(handles)),
                        title="PCA acquisitions"
                        if short_names
                        else "Acquisition row (sample_qc.tsv)",
                        frameon=False,
                        fontsize=8,
                    )
                else:
                    axis.scatter(xy[:, 0], xy[:, 1], color="#009E73", s=20)
                ratio = data.uns["pca_variance"]["variance_ratio"]
                axis.set(
                    xlabel=f"PC1 ({ratio[0]:.1%})",
                    ylabel=f"PC2 ({ratio[1]:.1%})",
                    title=f"{int(data.var.pca_eligible.sum()):,} complete varying groups",
                )
            else:
                axis.text(
                    0.5, 0.5, "PCA skipped: insufficient\ncomplete varying evidence", ha="center"
                )
                axis.set_axis_off()

        draw_pca(axes[2])
        fig.suptitle(f"AlphaPeptTools QC: {method}; log2 PCA without imputation", fontsize=11)
        for ext in ["pdf", "png"]:
            fig.savefig(out / f"QC_overview.{ext}", dpi=180, bbox_inches="tight", pad_inches=0.12)
        plt.close(fig)
        if "X_pca" in data.obsm:
            fig, axis = plt.subplots(figsize=(6, 5), constrained_layout=True)
            draw_pca(axis)
            fig.suptitle(f"{method}: log2 PCA without imputation", fontsize=11)
            for ext in ("pdf", "png"):
                fig.savefig(out / f"PCA.{ext}", dpi=180, bbox_inches="tight", pad_inches=0.12)
            plt.close(fig)
