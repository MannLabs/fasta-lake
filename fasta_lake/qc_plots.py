"""Readable sample QC without changing quantities or filling missing values.

Plotting uses Matplotlib: https://github.com/matplotlib/matplotlib.
FastaLake computes the displayed quantiles and missingness from the unchanged
AnnData matrix prepared with AlphaPeptTools; see downstream.py for its citation.
"""

from pathlib import Path


def intensity_distribution_summary(data):
    """Summarise measured positive reporting-group intensities per acquisition.

    Quantiles use log2 intensities. Nonpositive and nonfinite entries are absent
    measurements, never zero-valued observations in the displayed distribution.
    An empty acquisition remains a row with missing quantiles.
    """
    import numpy as np
    import pandas as pd

    rows = []
    for sample, values in zip(data.obs_names, data.X):
        values = np.asarray(values, dtype=float)
        measured = values[np.isfinite(values) & (values > 0)]
        quantiles = (
            np.quantile(np.log2(measured), [0, 0.25, 0.5, 0.75, 1])
            if len(measured)
            else [np.nan] * 5
        )
        rows.append(
            dict(
                sample=str(sample),
                measured_groups=len(measured),
                missing_groups=data.n_vars - len(measured),
                total_groups=data.n_vars,
                missing_fraction=(data.n_vars - len(measured)) / data.n_vars
                if data.n_vars
                else np.nan,
                **dict(
                    zip(["log2_min", "log2_q25", "log2_median", "log2_q75", "log2_max"], quantiles)
                ),
            )
        )
    return pd.DataFrame(rows).set_index("sample")


def plot_intensity_distributions(data, output, method, *, per_page=24):
    """Save paginated violins, missingness bars and their exact summary table.

    A constant or single-observation sample is drawn as a point rather than an
    undefined density. All acquired samples are shown in their original order.
    Violin widths describe within-sample distributions, not sample yield.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.backends.backend_pdf import PdfPages

    if type(per_page) is not int or per_page < 1:
        raise ValueError("per_page must be a positive integer")
    out = Path(output)
    summary = intensity_distribution_summary(data)
    summary.to_csv(out / "intensity_distribution.tsv", sep="\t")
    with (
        plt.rc_context(
            {
                "font.family": "DejaVu Sans",
                "font.size": 10,
                "axes.spines.top": False,
                "axes.spines.right": False,
                "figure.facecolor": "white",
                "savefig.facecolor": "white",
            }
        ),
        PdfPages(out / "Protein_intensity_violins.pdf") as pdf,
    ):
        for start in range(0, max(1, data.n_obs), per_page):
            end = min(start + per_page, data.n_obs)
            fig, (axis, missing) = plt.subplots(
                2,
                1,
                figsize=(max(8, (end - start) * 0.42), 7.2),
                gridspec_kw={"height_ratios": [4, 1]},
                sharex=True,
            )
            for position, i in enumerate(range(start, end), 1):
                values = np.asarray(data.X[i], dtype=float)
                logs = np.log2(values[np.isfinite(values) & (values > 0)])
                if len(logs) > 1 and np.ptp(logs) > 0:
                    parts = axis.violinplot(
                        logs, positions=[position], widths=0.8, showextrema=False
                    )
                    for body in parts["bodies"]:
                        body.set_facecolor("#397f91")
                        body.set_edgecolor("#256175")
                        body.set_alpha(0.65)
                    q1, med, q3 = np.quantile(logs, [0.25, 0.5, 0.75])
                    axis.vlines(position, q1, q3, color="#1e4656", linewidth=3)
                    axis.scatter(position, med, color="white", edgecolor="#1e4656", s=22, zorder=4)
                elif len(logs):
                    axis.scatter(position, logs[0], color="#1e4656", s=25)
                else:
                    axis.text(
                        position,
                        0.12,
                        "No\nquantity",
                        ha="center",
                        fontsize=8,
                        transform=axis.get_xaxis_transform(),
                        color="#64748b",
                    )
                axis.text(
                    position,
                    1.02,
                    f"{len(logs):,}",
                    ha="center",
                    fontsize=8,
                    transform=axis.get_xaxis_transform(),
                )
            x = np.arange(1, end - start + 1)
            missing.bar(x, summary.iloc[start:end].missing_fraction * 100, color="#b9c2c9")
            missing.set(ylim=(0, 100), ylabel="Missing\n(% groups)")
            missing.set_yticks([0, 50, 100])
            names = list(data.obs_names[start:end])
            short = all(len(str(n)) <= 20 for n in names)
            missing.set_xticks(
                x,
                names if short else list(range(start + 1, end + 1)),
                rotation=50 if short and len(names) > 4 else 0,
                ha="right" if short and len(names) > 4 else "center",
            )
            missing.set_xlabel(
                "Acquisition" if short else "Acquisition row (intensity_distribution.tsv)"
            )
            axis.set_ylabel("Protein reporting-group intensity (log₂ a.u.)")
            axis.grid(axis="y", color="#e6ebef", linewidth=0.6)
            axis.set_axisbelow(True)
            axis.set_xlim(0.3, max(1, end - start) + 0.7)
            fig.suptitle(
                f"Protein intensity distributions · {method}", x=0.1, ha="left", fontsize=15
            )
            fig.text(
                0.1,
                0.918,
                "Numbers: quantified groups   ·   Point: median   ·   Line: middle 50%",
                color="#475569",
                fontsize=9,
            )
            fig.text(
                0.1,
                0.025,
                "Measured positive quantities only; no normalization or imputation.\n"
                "Missingness uses all matrix groups. "
                "Group quantities do not resolve shared-peptide origin.",
                fontsize=8,
                color="#475569",
            )
            fig.subplots_adjust(left=0.1, right=0.98, top=0.85, bottom=0.23, hspace=0.08)
            pdf.savefig(fig)
            suffix = "" if start == 0 else f"_page{start // per_page + 1}"
            fig.savefig(out / f"Protein_intensity_violins{suffix}.png", dpi=180)
            plt.close(fig)
    return {
        "method": method,
        "samples": data.n_obs,
        "per_page": per_page,
        "table": "intensity_distribution.tsv",
        "units": "log2 arbitrary intensity units",
        "missingness": "Nonpositive/nonfinite entries excluded; original roster retained",
    }
