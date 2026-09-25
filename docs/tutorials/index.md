# Tutorials

The tutorials are Jupyter notebooks. Each one opens the outputs of a bundled example, so you can read them here with their stored results, or run them yourself on the same example after a Docker check or a native run. They are executed in continuous integration against freshly produced example outputs, so the code on this site is code that runs.

| Notebook | What you learn | Example it reads |
|---|---|---|
| [Your first result](01_first_result.ipynb) | What a completed run contains, how peptides map to study groups, what a blank cell means | CAMPI |
| [QC, replicates and PCA](02_qc_and_pca.ipynb) | Coverage, completeness, replicate CV and when PCA is eligible | Synthetic three-acquisition study |
| [Community and function](03_community_and_function.ipynb) | Unipept taxonomy, eggNOG annotation and the peptide to protein evidence network | CAMPI community |
| [Matched DNA and RNA](04_molecular_additions.ipynb) | What specimen-matched candidates add, and the ledger that explains every inclusion | CAMPI with synthetic molecular values |

## Run a tutorial yourself

Produce the example outputs once (each example refuses to overwrite an existing folder):

```bash
python tools/run_tutorials.py --layout       # the example commands, in order
export FASTALAKE_TUTORIAL_DATA=$PWD/tutorial_data
python tools/run_tutorials.py --check        # executes every notebook against those outputs
```

Then open any notebook in `docs/tutorials/` with Jupyter. The first cell reads `FASTALAKE_TUTORIAL_DATA`; point it at your own run directory to apply the same analysis to your study.
