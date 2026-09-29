#!/usr/bin/env python3
"""Build inspectable teaching frames and GIFs in a new directory.

Run from the source root: python tools/animate_concepts.py --out concept_frames_01
Requires matplotlib and Pillow (the docs environment). All values are synthetic.
Study selection and aggregation are executed by the production Python functions.
Razor is an explanatory trace of the fixed-count Rust hash-acc comparator; its
source hash is recorded. It is not a substitute for running the Rust engine.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from fasta_lake.study import group_searches, study_dictionary  # noqa: E402

INK, BLUE, TEAL, ORANGE = "#183c4a", "#0072b2", "#007e70", "#b34b12"
MUTED, LINE, BG = "#526772", "#d5dfe4", "#f7fafb"
PEPTIDES = dict(p1="PEPTIDEAK", p2="AASDFGHYK", p3="VVAGDAGLK", p4="TVQALADAK", p5="GATPQDLVK")
GRAPH = dict(A={"p1", "p2", "p3"}, B={"p3", "p4"}, C={"p4", "p5"}, D={"p5"})


def fnv1a(accession):
    value = 14695981039346656037
    for byte in accession.encode():
        value = ((value ^ byte) * 1099511628211) & ((1 << 64) - 1)
    return value


def razor_trace():
    graph = dict(A={"uA", "s1", "s2"}, B={"uB", "s1"}, C={"s2", "s3"}, D={"s3"})
    memberships = {pep: {p for p in graph if pep in graph[p]} for pep in set.union(*graph.values())}
    unique = {p: sum(len(memberships[q]) == 1 for q in qs) for p, qs in graph.items()}
    assignments = {
        q: min(ps, key=lambda p: (-unique[p], -len(graph[p]), fnv1a(p), p))
        for q, ps in memberships.items()
    }
    assert assignments == dict(uA="A", uB="B", s1="A", s2="A", s3="C")
    return graph, unique, assignments


def production_example():
    """Execute a real grouping/quantification fixture; publish no temporary paths."""
    dictionary = study_dictionary(GRAPH, set(GRAPH))
    assert dictionary["selected"] == [("A", 3), ("C", 2)]
    assert dictionary["peptide_to_group"] == dict(p1="A", p2="A", p3="A", p4="C", p5="C")
    with tempfile.TemporaryDirectory(prefix="fastalake-teaching-") as directory:
        base = Path(directory)
        fasta = base / "targets.fasta"
        fasta.write_text(
            "".join(
                f">{p}\n" + "GG".join(PEPTIDES[q] for q in sorted(qs)) + "\n"
                for p, qs in GRAPH.items()
            )
        )
        specifications = {
            "S1": (
                ["p1", "p3", "p4"],
                [("p1", 2, 10), ("p1", 3, 20), ("p3", 2, 5), ("p4", 2, 40), ("p5", 2, 1000)],
            ),
            "S2": (["p2", "p5"], [("p2", 2, 70), ("p5", 2, 0)]),
        }
        for name, (accepted, features) in specifications.items():
            sample = base / "search" / name
            sample.mkdir(parents=True)
            (sample / "config.json").write_text(
                json.dumps({"database": {"fasta": str(fasta)}, "mzml_paths": [name + ".mzML"]})
            )
            members = {q: ";".join(p for p in GRAPH if q in GRAPH[p]) for q in PEPTIDES}
            rows = ["peptide\tproteins\tlabel\tpeptide_q\tfilename\trank"]
            for q in sorted(set(accepted) | {f[0] for f in features}):
                score = "0.001" if q in accepted else "0.20"
                rows.append(f"{PEPTIDES[q]}\t{members[q]}\t1\t{score}\t{name}.mzML\t1")
            (sample / "results.sage.tsv").write_text("\n".join(rows) + "\n")
            rows = [f"peptide\tcharge\tproteins\tq_value\tscore\tspectral_angle\t{name}.mzML"]
            rows += [
                f"{PEPTIDES[q]}\t{charge}\t{members[q]}\t0.5\t1\t0.9\t{value}"
                for q, charge, value in features
            ]
            (sample / "lfq.tsv").write_text("\n".join(rows) + "\n")
        result = group_searches(base / "search", base / "grouped")
        with (base / "grouped/study_group_matrix.tsv").open() as handle:
            matrix = list(csv.DictReader(handle, delimiter="\t"))
        assert result["intensity_in"] == result["intensity_out"] == 145
        assert len(matrix) == 2
        by_group = {r["study_group"]: r for r in matrix}
        assert float(by_group["A"]["S1"]) == 35 and float(by_group["A"]["S2"]) == 70
        assert float(by_group["C"]["S1"]) == 40 and by_group["C"]["S2"] == ""
    return {
        "selected": dictionary["selected"],
        "peptide_to_group": dictionary["peptide_to_group"],
        "matrix": matrix,
        "accepted_intensity": 145,
        "rejected_S1_intensity": 1000,
        "scope": "Synthetic production-function fixture; no empirical or FDR claim",
    }


def canvas(title, stage, subtitle):
    fig = plt.figure(figsize=(12, 7), dpi=110, facecolor="white")
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    ax.text(0.05, 0.945, stage.upper(), color=TEAL, size=11, weight="bold")
    ax.text(0.05, 0.88, title, color=INK, size=23, weight="bold")
    ax.text(0.05, 0.825, subtitle, color=MUTED, size=12)
    ax.plot([0.05, 0.95], [0.785, 0.785], color=LINE)
    ax.text(0.05, 0.035, "FASTALAKE  /  Synthetic worked example", color=MUTED, size=10)
    return fig, ax


def note(ax, heading, body, takeaway):
    ax.add_patch(
        FancyBboxPatch(
            (0.64, 0.23), 0.30, 0.49, boxstyle="round,pad=0.018", facecolor=BG, edgecolor=LINE
        )
    )
    ax.text(0.665, 0.675, heading, size=15, weight="bold", color=INK, va="top")
    ax.text(0.665, 0.60, body, size=13, color=INK, va="top", linespacing=1.65)
    ax.text(0.05, 0.115, takeaway, size=13, color=TEAL, weight="bold")


def graph_panel(ax, graph, peptides, active=None, selected=(), assignments=None, badges=None):
    xs = {p: 0.21 + i * 0.105 for i, p in enumerate(graph)}
    ys = {q: 0.61 - i * 0.078 for i, q in enumerate(peptides)}
    for p, x in xs.items():
        color = TEAL if p in selected else INK
        ax.text(x, 0.735, p, ha="center", size=18, color=color, weight="bold")
        if badges:
            ax.text(x, 0.68, badges[p], ha="center", size=10, color=MUTED)
    for q, y in ys.items():
        ax.text(
            0.065, y, q, size=14, color=INK, va="center", weight="bold" if q == active else "normal"
        )
        for p in graph:
            if q not in graph[p]:
                continue
            assigned = assignments is not None and assignments.get(q) == p
            ax.scatter(
                xs[p],
                y,
                s=170,
                facecolor=TEAL if assigned else "white",
                edgecolor=TEAL if assigned else MUTED,
                linewidth=1.6,
                zorder=3,
            )
            if assigned:
                ax.text(xs[p], y, "✓", ha="center", va="center", color="white", size=10)
        if q == active:
            ax.plot([0.14, 0.57], [y, y], color=BLUE, alpha=0.22, linewidth=14, zorder=1)
    ax.text(0.065, 0.17, "○ compatible     ● assigned", color=MUTED, size=11)


def write_story(out, name, frames):
    directory = out / name
    directory.mkdir()
    images, records = [], []
    for i, (fig, caption) in enumerate(frames):
        stem = f"step_{i:02}"
        fig.savefig(directory / (stem + ".svg"), metadata={"Date": None})
        svg = directory / (stem + ".svg")
        svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
        fig.savefig(directory / (stem + ".png"), dpi=110)
        plt.close(fig)
        with Image.open(directory / (stem + ".png")) as im:
            images.append(im.convert("RGB").copy())
        records.append({"image": stem + ".svg", "caption": caption})
    images[0].save(
        directory / (name + ".gif"),
        save_all=True,
        append_images=images[1:],
        duration=[3500] * (len(images) - 1) + [5500],
        loop=0,
        optimize=False,
    )
    (directory / (name + ".json")).write_text(json.dumps(records, indent=2) + "\n")
    images[-1].save(directory / "summary.png")


def build(out):
    out.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "svg.fonttype": "none",
            "svg.hashsalt": "fastalake-concepts-v1",
        }
    )
    checked = production_example()
    graph, unique, assigned = razor_trace()
    frames = []
    explanations = [
        (
            None,
            "Count the evidence",
            "U = unique peptides\nT = all mapped peptides\n\nA: 1 unique, 3 tota"
            "l\nB: 1 unique, 2 total\nC: 0 unique, 2 total\nD: 0 unique, 1 total",
            "Count once on the full graph. These counts stay fixed.",
        ),
        (
            "s1",
            "A wins shared s1",
            "A and B each have\none unique peptide.\n\nA has more total\nevidence: 3 versus 2.",
            "Compare unique evidence first, then total evidence.",
        ),
        (
            "s2",
            "A wins shared s2",
            "A has unique support;\nC has none.\n\nThe first criterion\nalready resolves s2.",
            "Earlier assignments do not change the evidence counts.",
        ),
        (
            "s3",
            "C wins shared s3",
            "C and D both have\nzero unique peptides.\n\nC wins on total\nevidence: 2 versus 1.",
            "A protein supported only by shared peptides can survive.",
        ),
        (
            "uA",
            "Unique peptides",
            "uA maps only to A.\nuB maps only to B.\n\nEach goes to its sole\ncompatible protein.",
            "Illustration order only: hash-acc assignments are order-independent.",
        ),
        (
            None,
            "Keep A, B and C",
            "Each receives at least\none assignment.\nD receives none.\n\nThis se"
            "lects complete\nsearchable sequences.",
            "Selection supports candidates; it does not identify a biological source.",
        ),
    ]
    # Shared rows precede unique rows, as in Rust's descending-degeneracy order.
    seen = {}
    for i, (active, heading, body, takeaway) in enumerate(explanations):
        if active:
            seen[active] = assigned[active]
        if i >= 4:
            seen.update(uA="A", uB="B")
        fig, ax = canvas(
            "Selection parsimony: assign predicted peptides",
            f"Before search · {i + 1}/6",
            "Fixed evidence counts → one assignment per peptide → selected FASTA",
        )
        graph_panel(
            ax,
            graph,
            ["uA", "uB", "s1", "s2", "s3"],
            active,
            set(seen.values()),
            seen,
            {p: f"U {unique[p]} / T {len(graph[p])}" for p in graph},
        )
        note(ax, heading, body, takeaway)
        frames.append((fig, heading + ". " + body.replace("\n", " ") + " " + takeaway))
    write_story(out, "razor", frames)
    frames = []
    stages = [
        (
            [],
            {},
            "Start with accepted peptides",
            "Choose the protein that\nexplains the most\nunassigned peptides.\n\n"
            "Initial gains:\nA: 3; B: 2; C: 2; D: 1.",
            "Here the graph pools accepted target peptides from both acquisitions.",
        ),
        (
            ["A"],
            dict(p1="A", p2="A", p3="A"),
            "Select A: gain 3",
            "A covers p1, p2 and p3.\n\nRecalculate the gains:\nB: 1; C: 2; D: 1"
            ".\n\nAlready covered p3\nno longer counts for B.",
            "Study grouping recalculates gains after every selection.",
        ),
        (
            ["A", "C"],
            checked["peptide_to_group"],
            "Select C: gain 2",
            "C covers p4 and p5.\nAll five peptides now\nhave an assignment.\n\nB"
            " and D add no\nuncovered evidence.",
            "Each peptide stays with the first selected representative explaining it.",
        ),
        (
            ["A", "C"],
            checked["peptide_to_group"],
            "Two reporting groups",
            "A: p1, p2, p3\nC: p4, p5\n\nEqual gains would use\nlexical accession"
            " order.\nThe heuristic need not\nfind a minimum cover.",
            "A group label does not prove a unique protein or species origin.",
        ),
    ]
    for i, (selected, mapping, heading, body, takeaway) in enumerate(stages):
        fig, ax = canvas(
            "Study parsimony: cover accepted evidence",
            f"After search · {i + 1}/4",
            "Observed study graph → greedy coverage → one peptide-to-group dictionary",
        )
        covered = set(mapping)
        graph_panel(
            ax,
            GRAPH,
            list(PEPTIDES),
            selected=selected,
            assignments=mapping,
            badges={
                p: ("selected" if p in selected else f"gain {len(qs - covered)}")
                for p, qs in GRAPH.items()
            },
        )
        note(ax, heading, body, takeaway)
        frames.append((fig, heading + ". " + body.replace("\n", " ") + " " + takeaway))
    write_story(out, "study", frames)
    frames = []
    stories = [
        (
            "Keep the local acceptance gate",
            "S1 p5 fails peptide q.\nIts intensity of 1,000\nis excluded even t"
            "hough\np5 is accepted in S2.\n\nS2 p5 has no positive\nquantity to c"
            "ontribute.",
            "Peptide acceptance is acquisition-specific; this is not an MS1 q-value filter.",
        ),
        (
            "Assign each included feature",
            "Both charge states of\np1 follow p1 to A.\np3 also follows A.\np4 f"
            "ollows C.\n\nEach included feature\ncontributes once.",
            "Two features for p1 still count as one canonical peptide.",
        ),
        (
            "Sum within each acquisition",
            "S1 A: 10 + 20 + 5 = 35\nS1 C: 40\nS2 A: 70\nS2 C: no positive value"
            "\n\nInput = output = 145.",
            "Sum conservation checks accounting, not biological accuracy.",
        ),
        (
            "Read the resulting matrix",
            "Rows: reporting groups\nColumns: acquisitions\n\nBlank means unquan"
            "tified.\nIt does not establish\nbiological absence.",
            "This illustrates sum quantification. directLFQ estimates quantities differently.",
        ),
    ]
    for i, (heading, body, takeaway) in enumerate(stories):
        fig, ax = canvas(
            "Aggregation: evidence becomes a matrix",
            f"Quantification · {i + 1}/4",
            "Use the same study dictionary: p1, p2, p3 → A; p4, p5 → C",
        )
        if i < 3:
            rows = [
                ["S1", "p1 / 2+", "10", "A"],
                ["S1", "p1 / 3+", "20", "A"],
                ["S1", "p3 / 2+", "5", "A"],
                ["S1", "p4 / 2+", "40", "C"],
                ["S1", "p5 / 2+", "1,000", "excluded"],
                ["S2", "p2 / 2+", "70", "A"],
                ["S2", "p5 / 2+", "0", "no quantity"],
            ]
            table = ax.table(
                cellText=rows,
                colLabels=["Run", "Peptide / charge", "Intensity", "Destination"],
                bbox=[0.05, 0.24, 0.55, 0.49],
                colWidths=[0.08, 0.19, 0.11, 0.17],
            )
            for r in (5, 7):
                for c in range(4):
                    table[(r, c)].set_text_props(color=ORANGE)
        else:
            table = ax.table(
                cellText=[["A", "35", "70"], ["C", "40", ""]],
                colLabels=["Study group", "S1", "S2"],
                bbox=[0.065, 0.37, 0.51, 0.30],
            )
            ax.text(
                0.065,
                0.28,
                "S1: 4 features / 3 quantified peptides\nS2: 1 feature / 1 quantified peptide",
                size=12,
                color=MUTED,
            )
        table.auto_set_font_size(False)
        table.set_fontsize(11 if i < 3 else 17)
        for (r, _), cell in table.get_celld().items():
            cell.set_edgecolor(LINE)
            cell.set_facecolor(BG if r == 0 else "white")
            if r == 0:
                cell.set_text_props(weight="bold", color=INK)
        note(ax, heading, body, takeaway)
        frames.append((fig, heading + ". " + body.replace("\n", " ") + " " + takeaway))
    write_story(out, "aggregation", frames)
    checked["source_sha256"] = {
        str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (ROOT / "fasta_lake/study.py", ROOT / "rust/parsimony_engine/src/main.rs")
    }
    checked["razor_assignments"] = assigned
    (out / "CHECK.json").write_text(json.dumps(checked, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "stories": 3, "frames": 14, "accepted_intensity": 145}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    build(parser.parse_args().out)
