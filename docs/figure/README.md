# Workflow illustration

The illustration follows reference preparation, peptide lookup, protein selection,
searching and study-level reporting.

[Open the illustration](fastalake-pipeline.html) ·
[Light SVG](fastalake-pipeline-light.svg) ·
[Dark SVG](fastalake-pipeline-dark.svg) ·
[Download PDF](fastalake-pipeline.pdf)

## Update the artwork

The source labels, generator and fonts are in
[`docs/figure/`](https://github.com/MannLabs/fasta-lake/tree/main/docs/figure).
Use a separate Python environment:

```bash
python -m pip install -r docs/figure/requirements-figure.txt
python docs/figure/build_figure.py
python -m playwright install chromium
python docs/figure/render_figure.py
```

Inspect the light and dark renders at both sizes before committing the SVG,
HTML and PDF. Generated render records stay local. Bundled font licences are
kept with the fonts.
