#!/usr/bin/env python3
"""Draw a vertical workflow and a stepwise GIF from one shared set of labels."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from matplotlib import get_data_path
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
STAGES = [
    ("Protein reservoir", "Complete protein sequences from your chosen sources"),
    ("Study-supported sequence space", "Exact lookup of retained queries pooled across the study"),
    ("Acquisition-specific candidate set", "Exact lookup of this acquisition's retained queries"),
    ("Selected acquisition-specific database", "Selection parsimony retains complete proteins"),
    ("Sage database search", "Search the original spectra against the final target FASTA"),
    ("Study-level protein groups", "Assign accepted peptides to fixed reporting representatives"),
    ("Quantification and QC", "Sum or directLFQ; coverage, missingness and eligible PCA"),
]
WIDTH, HEIGHT = 860, 1330
INK, MUTED, GREEN = "#173934", "#54645f", "#28786a"
PALE, BORDER, PAPER = "#e9f4ef", "#c4d3cd", "#fbfcf9"
FONTS = Path(get_data_path()) / "fonts" / "ttf"


def font(size, bold=False):
    """Load a bundled font so the artwork does not depend on system font discovery."""
    return ImageFont.truetype(
        str(FONTS / ("DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf")), size
    )


def fit(draw, text, maximum, start, bold=False):
    """Choose a font size that keeps a single label inside its allotted width."""
    size = start
    while draw.textlength(text, font=font(size, bold)) > maximum and size > 12:
        size -= 1
    return font(size, bold)


def frame(active=None):
    """Draw all stages, emphasizing one while leaving the complete route readable."""
    canvas = Image.new("RGB", (WIDTH, HEIGHT), PAPER)
    d = ImageDraw.Draw(canvas)
    d.text((64, 45), "FastaLake", fill=INK, font=font(42, True))
    d.text((64, 105), "Protein databases shaped by your sample.", fill=MUTED, font=font(23))
    d.text((64, 153), "FOLLOW THE WORKFLOW", fill=GREEN, font=font(15, True))
    ys = [204, 341, 478, 615, 844, 981, 1118]
    for i, ((title, description), y) in enumerate(zip(STAGES, ys)):
        if i:
            top = ys[i - 1] + 102
            if i == 4:
                top = 806
            d.line((WIDTH // 2, top + 3, WIDTH // 2, y - 13), fill=BORDER, width=3)
            d.polygon(
                [(WIDTH // 2 - 7, y - 20), (WIDTH // 2 + 7, y - 20), (WIDTH // 2, y - 11)],
                fill=BORDER,
            )
        focus = active is None or i == active
        fill = PALE if focus else "#ffffff"
        outline = GREEN if focus else BORDER
        d.rounded_rectangle(
            (64, y, WIDTH - 64, y + 102),
            radius=14,
            fill=fill,
            outline=outline,
            width=3 if focus else 1,
        )
        d.ellipse((82, y + 20, 122, y + 60), fill=GREEN if focus else "#e4eae6")
        d.text((96, y + 26), str(i + 1), fill="white" if focus else INK, font=font(20, True))
        d.text((142, y + 19), title, fill=INK, font=fit(d, title, WIDTH - 230, 25, True))
        d.text((142, y + 61), description, fill=MUTED, font=fit(d, description, WIDTH - 230, 18))
    d.line((WIDTH // 2, 720, WIDTH // 2, 739), fill=BORDER, width=3)
    d.rounded_rectangle((133, 740, 727, 805), radius=10, fill="#fff3dc", outline="#d5b46d")
    d.text(
        (155, 750), "Optional: add matched DNA / RNA sequences", fill="#77551b", font=font(20, True)
    )
    d.text((155, 779), "After selection; keep the de novo baseline.", fill="#77551b", font=font(17))
    label = "Complete route" if active is None else f"Step {active + 1} of {len(STAGES)}"
    d.text((64, 1263), label, fill=GREEN, font=font(17, True))
    d.text((WIDTH - 356, 1263), "Predictions and mzML prepared upstream", fill=MUTED, font=font(13))
    return canvas


def main():
    """Write a static PNG, a looping GIF and a label/timing record."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "docs" / "assets")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    static = frame()
    static.save(args.out / "workflow.png")
    frames = [frame(i) for i in range(len(STAGES))] + [static]
    palette = static.quantize(colors=96)
    frames = [f.quantize(palette=palette) for f in frames]
    durations = [1800] * len(STAGES) + [3000]
    frames[0].save(
        args.out / "workflow.gif",
        save_all=True,
        append_images=frames[1:],
        duration=durations,
        loop=0,
        optimize=False,
        disposal=2,
    )
    record = {
        "generator": "tools/render_workflow.py",
        "dimensions": [WIDTH, HEIGHT],
        "frame_ms": durations,
        "stages": STAGES,
        "optional": "Matched DNA/RNA exact additions after selection and before search",
    }
    (args.out / "workflow.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"Wrote {len(frames)} frames and static artwork to {args.out}")


if __name__ == "__main__":
    main()
