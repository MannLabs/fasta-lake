"""Verify the installed optional backends on the measured CAMPI example."""

import json
import sys
from pathlib import Path

import anndata
import numpy as np
import pandas as pd

groups, quantities, analysis = map(Path, sys.argv[1:])
source = json.loads((groups / "summary.json").read_text())
quant = json.loads((quantities / "summary.json").read_text())
record = json.loads((analysis / "summary.json").read_text())
assert quant["method"] == "directlfq" and quant["directlfq_version"] == "0.3.3"
assert quant["acquisition_names"] == source["acquisition_names"]
matrix = pd.read_csv(quantities / "study_group_matrix.tsv", sep="\t", index_col=0)
objects = list(analysis.glob("*.h5ad"))
assert len(objects) == 1, objects
adata = anndata.read_h5ad(objects[0])
assert list(adata.obs_names) == source["acquisition_names"]
np.testing.assert_allclose(adata.X, matrix.T.to_numpy(), equal_nan=True)
assert (analysis / "QC_overview.png").is_file()
assert record["versions"]["alphapepttools"] == "0.4.0"
assert len(source["acquisition_names"]) == 2
assert record["pca"]["status"].startswith("SKIPPED:")
assert "X_pca" not in adata.obsm
print("PASS: directLFQ quantities, sample identity, AlphaPeptTools QC; two-sample PCA skipped.")
