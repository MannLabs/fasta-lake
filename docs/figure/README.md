# FastaLake pipeline figure

[Light PNG](renders/svg-light-1600-browser-light.png) · [Dark PNG](renders/svg-dark-1600-browser-dark.png) · [One-page PDF](fastalake-pipeline.pdf) · [Standalone HTML](fastalake-pipeline.html)

The figure follows MannLabs/fasta-lake main at `e072ca2a9b73dcd16f85ebc73bec9f2d9f9e2773`. Current v11 selection, search and reporting definitions were checked. The README uses both SVG themes; the site also offers the standalone HTML and PDF.

## Assets

- `fastalake-pipeline-light.svg` and `fastalake-pipeline-dark.svg`: self-contained vectors, embedded subset fonts, no scripts or external asset requests.
- `README_SNIPPET.md`: the `<picture>` block for the repository README.
- `fastalake-pipeline.html`: the same figure with native expandable stage explanations; follows the browser’s colour preference.
- `pipeline_data.json`: the shared source of stage labels and explanations.
- `build_figure.py`: handwritten SVG/HTML generation and a generated label-to-source ledger.
- `code_sources.json`, `PIPELINE_INVENTORY.md` and `LABEL_SOURCES.md`: pinned code provenance, full inventory and exact source lines for each final label.
- `fonts/`: original fonts, generated WOFF2 subsets, pinned source manifest and licences.

The `renders/` folder contains twelve browser PNGs: each fixed SVG theme at 880/1600 px under light/dark browser preferences, and the HTML at both widths/preferences. The smallest diagram text is 11.37 px at 880 px. All twelve final PNGs were inspected.

## Rebuild and render

Use Python 3.10 or newer in a separate virtual environment:

```sh
python -m pip install -r docs/figure/requirements-figure.txt
python docs/figure/build_figure.py
python -m playwright install chromium
python docs/figure/render_figure.py
```

An existing browser can be selected with `--chromium /path/to/chrome`. The builder uses only the bundled font files; the renderer blocks HTTP(S) requests. It also exports the one-page light PDF.

After changing the artwork, inspect all twelve new PNGs and update `VISUAL_REVIEW.json` with their hashes and findings. Then run:

```sh
python docs/figure/verify_figure.py
```

This requires `xmllint` on PATH and Git history containing the recorded source commit (use a full clone, or fetch that commit before verification). It checks XML, embedded fonts, scripts/hrefs, literal labels and citations, reproducible outputs, browser geometry and the saved visual review. The current receipts are [`VERIFICATION.json`](VERIFICATION.json), [`RENDER_RECEIPT.json`](RENDER_RECEIPT.json) and [`VISUAL_REVIEW.json`](VISUAL_REVIEW.json). The verifier reads the pinned code revision from Git and checks its hashes. Manuscript snapshots and private review records are not included.
