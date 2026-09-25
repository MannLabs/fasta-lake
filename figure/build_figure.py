#!/usr/bin/env python3
"""Build the README SVGs, standalone HTML and source-label ledger from one data file.

Run with Python >=3.10 and the versions in requirements-figure.txt. The generator
writes SVG markup directly; no diagram framework, browser, network or JavaScript
is needed to build it. Google Fonts originals and OFL licences are in fonts/.
"""
from __future__ import annotations

import ast
import base64
import hashlib
import io
import json
from html import escape
from pathlib import Path

from fontTools import subset
from fontTools.ttLib import TTFont

HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE / "pipeline_data.json").read_text())
SOURCES = json.loads((HERE / "code_sources.json").read_text())
assert DATA["source_commit"] == SOURCES["commit"]
PALETTES = {
    "light": dict(bg="#fcfdfc", ink="#243348", secondary="#4d5d70", muted="#607184",
                  line="#bac8d0", rule="#e2e8eb", input="#ffffff", border="#dce4e8",
                  study="#ebf5f0", study_line="#d4e8de", study_accent="#34856c",
                  acquisition="#eef2fd", acquisition_line="#dce3f4", acquisition_accent="#617ed0",
                  results="#fbf2e6", results_line="#efe2cf", results_accent="#b18040"),
    "dark": dict(bg="#111923", ink="#e8edf5", secondary="#becbd9", muted="#a1b2c4",
                 line="#52657b", rule="#2b3b4c", input="#192432", border="#35465a",
                 study="#1c342f", study_line="#2b5145", study_accent="#85c6a9",
                 acquisition="#202e48", acquisition_line="#354568", acquisition_accent="#9bb3fa",
                 results="#352e26", results_line="#554636", results_accent="#dfb475"),
}
USED = set()
BOUNDS = []


def label(key):
    """Resolve a visible label and require a source for it."""
    value = DATA["labels"][key]
    assert value["sources"], key
    for source in value["sources"]:
        assert source in SOURCES["sources"], (key, source)
    USED.add(key)
    return value["text"]


def font_css():
    """Subset both variable fonts to the actual figure/page vocabulary and embed WOFF2."""
    chars = "".join(v["text"] for v in DATA["labels"].values())
    chars += "".join(s["detail"] for s in DATA["stages"])
    chars += "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz →−/()[]:,.;!?—–×≤"
    rules = []
    records = []
    for family, filename in [("Plus Jakarta Sans", "PlusJakartaSans[wght].ttf"),
                             ("JetBrains Mono", "JetBrainsMono[wght].ttf")]:
        original = HERE / "fonts" / filename
        font = TTFont(original, recalcTimestamp=False)
        options = subset.Options()
        options.recalc_timestamp = False
        options.layout_features = ["*"]
        sub = subset.Subsetter(options=options)
        sub.populate(text=chars)
        sub.subset(font)
        font.flavor = "woff2"
        output = io.BytesIO()
        font.save(output, reorderTables=True)
        raw = output.getvalue()
        target = filename.replace("[wght].ttf", "-subset.woff2")
        (HERE / "fonts" / target).write_bytes(raw)
        encoded = base64.b64encode(raw).decode("ascii")
        rules.append(f"@font-face{{font-family:'{family}';font-style:normal;font-weight:100 900;"
                     f"font-display:block;src:url(data:font/woff2;base64,{encoded}) format('woff2')}}")
        records.append(dict(file=target, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
    return "\n".join(rules), records


FONT_CSS, FONT_RECORDS = font_css()


class Drawing:
    """Small deterministic SVG writer with explicit geometry and label bounds."""
    def __init__(self, theme):
        """Start an empty drawing with the selected light or dark palette."""
        self.theme = theme
        self.c = PALETTES[theme]
        self.parts = []

    def add(self, markup):
        """Append trusted SVG markup to the drawing in paint order."""
        self.parts.append(markup)

    def rect(self, x, y, w, h, fill, stroke=None, radius=16, attrs=""):
        """Draw a rounded rectangle, optionally with a border and extra SVG attributes."""
        border = f' stroke="{stroke}" stroke-width="1.2"' if stroke else ""
        self.add(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{radius}" '
                 f'fill="{fill}"{border} {attrs}/>')

    def path(self, d, stroke, width=1.8, fill="none", attrs=""):
        """Draw an SVG path with consistent rounded stroke joins."""
        self.add(f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{width}" '
                 f'stroke-linecap="round" stroke-linejoin="round" {attrs}/>')

    def text(self, key, x, y, size=18, weight=400, color=None, mono=False, max_width=None):
        """Resolve a cited label, escape its text and record its allocated bounds."""
        text = label(key)
        assert size * 880 / DATA["width"] >= 11, (key, size)
        family = "JetBrains Mono" if mono else "Plus Jakarta Sans"
        color = color or self.c["ink"]
        self.add(f'<text data-label="{key}" x="{x}" y="{y}" fill="{color}" '
                 f'font-family="{family}" font-size="{size}" font-weight="{weight}"'
                 f' font-variant-ligatures="none">{escape(text)}</text>')
        if max_width:
            BOUNDS.append(dict(theme=self.theme, key=key, x=x, y=y, max_width=max_width, size=size))

    def chevron(self, x, y, direction="right", color=None):
        """Draw a rightward or downward arrowhead at the given coordinates."""
        d = f"M{x-4} {y-6} L{x+2} {y} L{x-4} {y+6}"
        if direction == "down":
            d = f"M{x-6} {y-4} L{x} {y+2} L{x+6} {y-4}"
        self.path(d, color or self.c["line"], 2.0)

    def arrow(self, d, x, y, down=False):
        """Draw a connecting path and its terminal arrowhead."""
        self.path(d, self.c["line"], 1.6)
        self.chevron(x, y, "down" if down else "right")

    def card(self, stage, x, y, w, h):
        """Open a stage group with accessible detail text and a coloured card."""
        group = stage["group"]
        self.add(f'<g id="{self.theme}-{stage["id"]}" class="stage" '
                 f'aria-label="{escape(" ".join(label(k) for k in stage["title"]))}">')
        self.add(f'<title>{escape(stage["detail"])}</title>')
        self.rect(x, y, w, h, self.c[group], self.c[group+"_line"], 18)
        self.rect(x+1, y+23, 4, h-46, self.c[group+"_accent"], radius=2)

    def close_card(self):
        """Close the current stage group."""
        self.add('</g>')

    def icon(self, kind, x, y, color):
        """Draw the named stage icon from vector primitives."""
        self.add(f'<g transform="translate({x} {y})" opacity="0.85" aria-hidden="true">')
        if kind == "lookup":
            self.path("M1 2 H19 M1 9 H14 M1 16 H11 M25 22 L31 28", color, 1.7)
            self.add(f'<circle cx="22" cy="18" r="7" stroke="{color}" stroke-width="1.7" fill="none"/>')
        elif kind == "select":
            self.path("M0 4 H18 M0 12 H18 M0 20 H18 M23 12 L27 16 L34 7", color, 1.7)
        elif kind == "spectra":
            self.path("M0 26 H34 M4 26 V16 M10 26 V6 M17 26 V12 M24 26 V0 M30 26 V18", color, 1.7)
        elif kind == "matrix":
            for row in range(3):
                for col in range(4):
                    self.rect(col*8, row*8, 5, 5, color, radius=1.3)
        elif kind == "qc":
            self.path("M0 0 V27 H34 M6 20 L13 11 L21 16 L31 5", color, 1.7)
        self.add('</g>')

    def render(self):
        """Return the complete SVG for this theme and record visible label bounds."""
        c = self.c
        self.add(f'<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="1080" '
                 f'viewBox="0 0 1200 1080" role="img" aria-labelledby="{self.theme}-title {self.theme}-desc">')
        self.add(f'<title id="{self.theme}-title">{escape(label("accessible_title"))}</title>')
        self.add(f'<desc id="{self.theme}-desc">{escape(label("accessible_description"))}</desc>')
        self.add(f'<defs><style>{FONT_CSS}\ntext{{text-rendering:geometricPrecision}}</style></defs>')
        self.rect(0, 0, 1200, 1080, c["bg"], radius=24)
        # A simple lake mark, drawn as vectors rather than an external logo asset.
        self.path("M48 50 C57 43 65 57 74 50 S91 57 100 50 M48 61 C57 54 65 68 74 61 S91 68 100 61 M48 72 C57 65 65 79 74 72 S91 79 100 72", c["study_accent"], 2.5)
        self.text("brand", 118, 76, 34, 750, max_width=350)
        self.text("subtitle", 48, 115, 19, 400, c["secondary"], max_width=850)
        self.rect(925, 46, 227, 38, c["input"], c["border"], 19)
        self.add(f'<circle cx="944" cy="65" r="3.5" fill="{c["study_accent"]}"/>')
        self.text("badge", 957, 71, 16, 500, c["secondary"], True, 184)
        self.path("M48 140 H1152", c["rule"], 1)

        self.text("study", 48, 174, 16, 700, c["study_accent"], max_width=310)
        self.rect(48, 196, 294, 64, c["input"], c["border"], 13)
        self.text("reservoir", 66, 222, 20, 650, max_width=257)
        self.text("reservoir_type", 66, 246, 16, 400, c["secondary"], True, 257)
        self.rect(48, 274, 294, 80, c["input"], c["border"], 13)
        self.text("predictions", 66, 302, 20, 650, max_width=220)
        self.text("csv", 295, 301, 16, 500, c["muted"], True, 32)
        self.text("predictions_gloss", 66, 332, 15.5, 400, c["secondary"], max_width=258)
        self.arrow("M350 228 H377", 379, 228)
        self.arrow("M350 314 H377", 379, 314)
        st = DATA["stages"][0]
        self.card(st, 392, 196, 760, 158)
        self.text("pool_title", 417, 232, 26, 700, max_width=640)
        self.text("pool_body", 417, 267, 18.5, 400, c["secondary"], max_width=700)
        self.text("pool_il", 417, 293, 18.5, 400, c["secondary"], max_width=700)
        self.text("pool_filter", 417, 329, 16, 500, c["study_accent"], True, 700)
        self.icon("lookup", 1099, 215, c["study_accent"])
        self.close_card()

        # First elbow feeds only acquisition candidate retrieval. The original
        # eligible predictions enter selection separately, below the lookup filter.
        self.arrow("M544 354 V373 Q544 383 534 383 H354 Q344 383 344 393 V438", 344, 440, True)
        self.text("acquisition", 48, 422, 16, 700, c["acquisition_accent"], max_width=290)
        self.rect(428, 398, 344, 34, c["input"], c["border"], 10)
        self.text("all_predictions", 444, 420, 16, 400, c["secondary"], True, 312)
        self.arrow("M600 432 V441", 600, 443, True)
        self.rect(808, 398, 344, 34, c["input"], c["border"], 10)
        self.text("spectra", 830, 420, 16, 400, c["secondary"], True, 300)
        self.arrow("M980 432 V441", 980, 443, True)

        for i, x in enumerate([48, 428, 808]):
            stage = DATA["stages"][i+1]
            self.card(stage, x, 452, 344, 218)
            for j, key in enumerate(stage["title"]):
                self.text(key, x+24, 489+j*31, 25, 700, max_width=288)
            for j, key in enumerate(stage["body"]):
                self.text(key, x+24, 557+j*26, 18, 400, c["secondary"], max_width=296)
            if i == 2:
                self.text(stage["tool"][0], x+24, 619, 16, 500, c["acquisition_accent"], True, 296)
                self.text(stage["tool"][1], x+24, 644, 16, 500, c["secondary"], True, 296)
            else:
                self.text(stage["tool"][0], x+24, 640, 16, 500, c["acquisition_accent"], True, 296)
            self.icon(["lookup", "select", "spectra"][i], x+285, 476, c["acquisition_accent"])
            self.close_card()
        self.chevron(412, 562)
        self.chevron(792, 562)

        # Summaries are an auxiliary study output; grouping reads the search
        # results directly, not these member-set matrices.
        self.arrow("M980 670 V696 Q980 706 970 706 H356 Q344 706 344 718 V782", 344, 784, True)
        self.path("M980 706 V721", c["line"], 1.6)
        self.rect(836, 724, 316, 38, c["input"], c["border"], 10)
        self.icon("matrix", 850, 732, c["muted"])
        self.text("summary_output", 893, 749, 17, 500, c["secondary"], max_width=244)
        self.text("results", 48, 762, 16, 700, c["results_accent"], max_width=280)
        group = DATA["stages"][4]
        self.card(group, 48, 798, 596, 198)
        self.text("group_title", 72, 839, 26, 700, max_width=520)
        self.text("group_body_1", 72, 884, 18.5, 400, c["secondary"], max_width=546)
        self.text("group_body_2", 72, 912, 18.5, 400, c["secondary"], max_width=546)
        self.text("group_tool", 72, 968, 16, 500, c["results_accent"], True, 546)
        self.close_card()
        self.chevron(664, 897)
        qc = DATA["stages"][5]
        self.card(qc, 684, 798, 468, 198)
        self.text("qc_title", 708, 839, 26, 700, max_width=368)
        self.text("qc_body_1", 708, 884, 18.5, 400, c["secondary"], max_width=419)
        self.text("qc_body_2", 708, 912, 18.5, 400, c["secondary"], max_width=419)
        self.text("qc_tool", 708, 968, 16, 500, c["results_accent"], True, 419)
        self.icon("qc", 1100, 816, c["results_accent"])
        self.close_card()
        self.path("M48 1022 H1152", c["rule"], 1)
        self.text("acquisition_gloss", 48, 1053, 16, 400, c["muted"], max_width=550)
        self.text("footer", 843, 1053, 16, 400, c["muted"], max_width=309)
        self.add('</svg>')
        return '\n'.join(self.parts)


def source_link(source):
    """Return a human-readable citation and URL pinned to the audited code."""
    s = SOURCES["sources"][source]
    url = SOURCES["repository"] + "/blob/" + SOURCES["commit"] + "/" + s["path"]
    url += f'#L{s["first_line"]}-L{s["last_line"]}'
    title = f'{s["path"]}::{s["symbol"]} L{s["first_line"]}–{s["last_line"]}'
    return title, url


def implementation_ref(symbol):
    """Cite figure-generation claims to this builder, separately from pipeline claims."""
    tree = ast.parse(Path(__file__).read_text())
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == symbol)
    return f'[build_figure.py::{symbol} L{node.lineno}–{node.end_lineno}](build_figure.py#L{node.lineno}-L{node.end_lineno})'


def build_html(svgs):
    """Combine both SVG themes with expandable stage explanations in one page."""
    details = []
    for stage in DATA["stages"]:
        title = ' '.join(label(k) for k in stage["title"])
        details.append(f'<details class="detail {stage["group"]}"><summary>{escape(title)}</summary>'
                       f'<p>{escape(stage["detail"])}</p></details>')
    sections = []
    for key in ["scope", "confidence", "optional", "sources"]:
        sections.append(f'<section><h2>{escape(label(key+"_heading"))}</h2>'
                        f'<p>{escape(label(key+"_detail"))}</p></section>')
    _, url = source_link("record")
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light dark"><title>{escape(label("html_title"))}</title>
<style>{FONT_CSS}
:root{{color-scheme:light dark;--bg:#fcfdfc;--ink:#243348;--secondary:#4d5d70;--line:#dce4e8;--surface:#fff;--teal:#34856c}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font-family:'Plus Jakarta Sans',sans-serif}}
main{{max-width:1600px;margin:auto}}figure{{margin:0}}.theme svg{{display:block;width:100%;height:auto}}.dark{{display:none}}
.notes{{max-width:1104px;margin:16px auto 72px;padding:0 32px;font-size:16px;line-height:1.7}}
h1,h2{{font-weight:650;letter-spacing:-.025em}}h1{{font-size:26px;margin:28px 0 20px}}h2{{font-size:19px;margin:28px 0 8px}}p{{color:var(--secondary);margin:8px 0 16px}}
.details-grid{{display:grid;grid-template-columns:1fr 1fr;gap:12px}}.detail{{padding:16px 20px;border:1px solid var(--line);border-radius:12px;background:var(--surface)}}
summary{{cursor:pointer;font-weight:600;line-height:1.5}}details p{{font-size:15px;margin-bottom:0}}a{{color:var(--teal);text-underline-offset:3px}}
@media(prefers-color-scheme:dark){{:root{{--bg:#111923;--ink:#e8edf5;--secondary:#becbd9;--line:#35465a;--surface:#192432;--teal:#85c6a9}}.light{{display:none}}.dark{{display:block}}}}
@media(max-width:680px){{.details-grid{{grid-template-columns:1fr}}.notes{{padding:0 22px}}}}
@media print{{.notes{{display:none}}body{{background:#fff}}.light{{display:block}}.dark{{display:none}}}}
</style></head><body><main>
<figure aria-label="{escape(label('accessible_title'))}"><div class="theme light">{svgs['light']}</div><div class="theme dark">{svgs['dark']}</div></figure>
<article class="notes"><h1>{escape(label('details_heading'))}</h1><div class="details-grid">{''.join(details)}</div>
{''.join(sections)}<p><a href="{url}">{escape(label('source_link'))}</a></p></article>
</main></body></html>'''


def main():
    """Write both SVGs, HTML, the README snippet, source ledger and build receipts."""
    svgs = {theme: Drawing(theme).render() for theme in PALETTES}
    for theme, svg in svgs.items():
        (HERE / f'fastalake-pipeline-{theme}.svg').write_text(svg+'\n')
    (HERE / 'fastalake-pipeline.html').write_text(build_html(svgs)+'\n')
    snippet = f'''<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/figure/fastalake-pipeline-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="docs/figure/fastalake-pipeline-light.svg">
  <img src="docs/figure/fastalake-pipeline-light.svg" alt="{escape(label('accessible_title'))}: {escape(label('accessible_description'))}" width="1200">
</picture>
'''
    (HERE / 'README_SNIPPET.md').write_text(snippet)
    lines = ['## Final figure labels and source lines', '',
             'Generated from `pipeline_data.json` and `code_sources.json`; the code commit is '+DATA['source_commit']+'.', '',
             '| Literal label | Exact source |', '|---|---|']
    for key in sorted(USED):
        entry = DATA['labels'][key]
        refs = '<br>'.join(f'[{title}]({url})' for title,url in map(source_link,entry['sources']))
        if key == 'sources_detail':
            refs += '<br>' + implementation_ref('font_css') + '<br>' + implementation_ref('main')
        if key in ('footer', 'sources_detail'):
            refs += '<br>[code_sources.json::commit L3](code_sources.json#L3)'
        lines.append('| '+entry['text'].replace('|',r'\|')+' | '+refs+' |')
    lines += ['', '### HTML stage details and SVG tooltips', '', '| Literal detail | Exact source |', '|---|---|']
    for stage in DATA['stages']:
        refs = '<br>'.join(f'[{title}]({url})' for title,url in map(source_link,stage['sources']))
        lines.append('| '+stage['detail'].replace('|',r'\|')+' | '+refs+' |')
    ledger = '\n'.join(lines)+'\n'
    (HERE / 'LABEL_SOURCES.md').write_text(ledger)
    inventory = HERE / 'PIPELINE_INVENTORY.md'
    if inventory.exists():
        text = inventory.read_text().split('## Final figure labels and source lines')[0].rstrip()
        inventory.write_text(text+'\n\n'+ledger)
    (HERE / 'LAYOUT_BOUNDS.json').write_text(json.dumps(BOUNDS,indent=2)+'\n')
    record = dict(source_commit=DATA['source_commit'], fonts=FONT_RECORDS,
                  literal_labels=len(USED), minimum_font_px_at_880=min(b['size'] for b in BOUNDS)*880/1200,
                  outputs={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in [
                      'fastalake-pipeline-light.svg','fastalake-pipeline-dark.svg','fastalake-pipeline.html',
                      'README_SNIPPET.md','LABEL_SOURCES.md']})
    (HERE / 'BUILD_RECEIPT.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))


if __name__ == '__main__':
    main()
