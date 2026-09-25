"""Static community-context figures with explicit measures, references and citations.

Unipept assignments: https://github.com/unipept/unipept-cli;
Mesuere et al. (2016), https://doi.org/10.1093/bioinformatics/btw039.
These static summaries are FastaLake plots built with Matplotlib
(https://github.com/matplotlib/matplotlib); the interactive taxonomy views
separately reuse Unipept Visualizations.
"""

from pathlib import Path

COLORS = {
    "Host": "#a65f69",
    "Bacteria": "#398795",
    "Archaea": "#d29b43",
    "Other assigned taxa": "#8575ae",
    "Unresolved": "#b7c1c9",
}


def plot_community_context(payload, output, *, per_page=24):
    """Plot both count and summed-signal fractions; keep empty acquisitions visible.

    Taxonomic fractions and contaminant-reference overlap answer different
    questions. The latter is drawn only when an explicit reference was supplied.
    """
    try:
        import matplotlib
    except ImportError:
        return dict(
            status="SKIPPED",
            reason="Static exports require matplotlib; interactive HTML is available",
        )
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.patches import Patch

    if type(per_page) is not int or per_page < 1:
        raise ValueError("per_page must be a positive integer")
    out = Path(output)
    kinds = ["Holobiont_balance"]
    if payload["contaminant_reference_supplied"]:
        kinds.append("Contaminant_reference_overlap")
    with plt.rc_context(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    ):
        for kind in kinds:
            with PdfPages(out / (kind + ".pdf")) as pdf:
                for start in range(0, max(1, len(payload["samples"])), per_page):
                    subset = payload["samples"][start : start + per_page]
                    fig, axes = plt.subplots(
                        1, 2, figsize=(12, max(4.8, len(subset) * 0.34 + 2.8)), sharey=True
                    )
                    names = [x["sample"] for x in subset]
                    labels = (
                        names
                        if all(len(n) <= 22 for n in names)
                        else [str(i) for i in range(start + 1, start + len(names) + 1)]
                    )
                    y = np.arange(len(names))
                    for axis, measure, title in zip(
                        axes,
                        ["peptides", "signal"],
                        ["Accepted peptide counts", "Summed LFQ-feature intensity"],
                    ):
                        total = np.array([x["totals"][measure] for x in subset])
                        offset = np.zeros(len(subset))
                        categories = (
                            list(COLORS)
                            if kind == "Holobiont_balance"
                            else ["Reference overlap", "Remaining signal"]
                        )
                        for category in categories:
                            if kind == "Holobiont_balance":
                                values = np.array([x["balance"][category][measure] for x in subset])
                            else:
                                values = np.array([x["overlap"][measure] for x in subset])
                                if category == "Remaining signal":
                                    values = total - values
                            percent = (
                                np.divide(values, total, out=np.zeros(len(total)), where=total > 0)
                                * 100
                            )
                            colour = (
                                COLORS[category]
                                if kind == "Holobiont_balance"
                                else "#cb8548"
                                if category == "Reference overlap"
                                else "#d6dee3"
                            )
                            axis.barh(
                                y, percent, left=offset, color=colour, height=0.66, label=category
                            )
                            for i, value in enumerate(percent):
                                if value >= 9:
                                    axis.text(
                                        offset[i] + value / 2,
                                        i,
                                        f"{value:.0f}%",
                                        ha="center",
                                        va="center",
                                        fontsize=9,
                                        color="#152f3e",
                                    )
                            offset += percent
                        for i, value in enumerate(total):
                            if value == 0:
                                axis.text(
                                    2,
                                    i,
                                    "No quantity" if measure == "signal" else "No peptides",
                                    va="center",
                                    color="#64748b",
                                    fontsize=9,
                                )
                        axis.set(xlim=(0, 100), xlabel="Share of all evidence (%)", title=title)
                        axis.set_yticks(y, labels)
                        axis.grid(axis="x", alpha=0.15)
                        axis.set_axisbelow(True)
                    axes[0].invert_yaxis()
                    axes[0].set_ylabel(
                        "Acquisition" if labels == names else "Acquisition row (source table)"
                    )
                    heading = (
                        "Holobiont balance"
                        if kind == "Holobiont_balance"
                        else "Contaminant-reference overlap"
                    )
                    fig.suptitle(heading, x=0.08, ha="left", fontsize=17)
                    if kind == "Holobiont_balance":
                        handles = [Patch(color=c, label=n) for n, c in COLORS.items()]
                        note = (
                            f"Host: NCBI taxon {payload['host_taxon_id']}. "
                            "Unresolved evidence stays in the denominator."
                        )
                    else:
                        handles = [
                            Patch(color="#cb8548", label="Reference overlap"),
                            Patch(color="#d6dee3", label="Remaining evidence"),
                        ]
                        note = (
                            "Exact canonical peptide overlap; "
                            "shared sequences do not establish laboratory contamination."
                        )
                    fig.legend(
                        handles=handles,
                        loc="upper center",
                        bbox_to_anchor=(0.53, 0.91),
                        ncol=len(handles),
                        frameon=False,
                        fontsize=9,
                    )
                    fig.text(0.08, 0.11, note, fontsize=9, color="#475569")
                    fig.text(
                        0.08,
                        0.065,
                        "Reference: " + payload["database"] + ". Signal fractions are not biomass.",
                        fontsize=8,
                        color="#475569",
                    )
                    fig.text(
                        0.08,
                        0.025,
                        "Taxonomy: Unipept — Mesuere et al. (2016), "
                        "doi:10.1093/bioinformatics/btw039. "
                        "Plot: FastaLake.",
                        fontsize=8,
                        color="#475569",
                    )
                    fig.subplots_adjust(top=0.78, bottom=0.25, left=0.15, right=0.97, wspace=0.08)
                    pdf.savefig(fig)
                    suffix = "" if start == 0 else "_page" + str(start // per_page + 1)
                    for ext in ["png", "svg"]:
                        fig.savefig(out / (kind + suffix + "." + ext), dpi=180)
                    plt.close(fig)
    return dict(status="PASS", samples=len(payload["samples"]), figures=kinds, per_page=per_page)
