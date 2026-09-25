"""Synthetic installed-workflow check, not measured biological evidence."""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import anndata
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    fasta = out / "synthetic.fasta"
    fasta.write_text(">A\nMACDEK\n>B\nMQQQQK\n>C\nMSSSSK\n")
    # The third protein has a missing observation; it must stay missing and leave PCA.
    for sample, values in zip(("S1", "S2", "S3"), ((100, 20, 10), (150, 21, 0), (300, 8, 50))):
        folder = out / "search" / sample
        folder.mkdir(parents=True)
        mzml = out / (sample + ".mzML")
        mzml.write_text("Synthetic format fixture: no measured spectra and no search performed.\n")
        (folder / "config.json").write_text(
            json.dumps({"database": {"fasta": str(fasta)}, "mzml_paths": [str(mzml)]})
        )
        psms = ["peptide\tproteins\tlabel\tpeptide_q\tfilename\n"]
        lfq = ["peptide\tcharge\tproteins\t" + sample + ".mzML\n"]
        for (peptide, protein), value in zip(
            (("ACDEK", "A"), ("QQQQK", "B"), ("SSSSK", "C")), values
        ):
            psms.append(f"{peptide}\t{protein}\t1\t0.001\t{sample}.mzML\n")
            lfq.append(f"{peptide}\t2\t{protein}\t{value}\n")
        (folder / "results.sage.tsv").write_text("".join(psms))
        (folder / "lfq.tsv").write_text("".join(lfq))
    from fasta_lake.study import group_searches

    group_searches(out / "search", out / "study")
    for method in ("sum", "directlfq"):

        def command(arguments):
            result = subprocess.run(
                [sys.executable, "-m", "fasta_lake.cli", *arguments], text=True, capture_output=True
            )
            (out / (method + "_" + arguments[0] + ".log")).write_text(result.stdout + result.stderr)
            if result.returncode:
                raise RuntimeError(result.stderr)

        command(
            [
                "quantify-study",
                "--groups",
                str(out / "study"),
                "--method",
                method,
                "--out",
                str(out / method),
                "--threads",
                "1",
            ]
        )
        data = anndata.read_h5ad(out / method / "qc" / "study.h5ad")
        record = json.loads((out / method / "qc" / "summary.json").read_text())
        assert record["pca"]["status"] == "PASS", record["pca"]
        for name in ("QC_overview.png", "QC_overview.pdf", "PCA.png", "PCA.pdf"):
            assert (out / method / "qc" / name).stat().st_size > 1000, name
        assert list(data.obs_names) == ["S1", "S2", "S3"]
        assert np.isnan(data["S2", "C"].X).all()
        eligible = np.asarray(data.var["pca_eligible"], dtype=bool)
        assert not eligible[list(data.var_names).index("C")]
        scores = np.asarray(data.obsm["X_pca"])
        x = data.layers["log2"][:, eligible]
        np.testing.assert_allclose(
            (x - x.mean(axis=0)) @ data.varm["PCs"][eligible], scores, atol=1e-5, rtol=1e-5
        )
        assert np.isfinite(scores).all()
    receipt = {
        "status": "PASS",
        "scope": "synthetic three-acquisition workflow, no biological evidence",
        "methods": ["sum", "directlfq"],
        "pca_fits": 2,
        "automatic_qc": True,
        "missing_observation_preserved": True,
    }
    (out / "VALIDATION.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
