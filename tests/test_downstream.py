"""Real AlphaPeptTools checks of sample identity, missingness and analysis definitions."""

import json

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner
from test_simplified_study import P1, P2, P3, search

from fasta_lake.cli import cli
from fasta_lake.downstream import analyze_study
from fasta_lake.quantification import quantify_study
from fasta_lake.study import group_searches

pytest.importorskip("alphapepttools")
ad = pytest.importorskip("anndata")


@pytest.fixture
def study(tmp_path):
    # Positive sub-unit intensities must remain detected despite negative log2 values.
    profiles = {"NA": [0.5, 40, 5], "S2": [1, 10, 0], "S3": [2, 30, 20], "empty": [0, 0, 0]}
    for sample, values in profiles.items():
        search(
            tmp_path,
            sample,
            {"A": P1, "B": P2, "C": P3},
            [(p, g, 1, 0.001) for p, g in zip([P1, P2, P3], "ABC")],
            [(p, 2, g, v) for p, g, v in zip([P1, P2, P3], "ABC", values)],
        )
    group_searches(tmp_path / "search", tmp_path / "groups")
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("sample\treplicates\tsubject\nS3\tR\tX\nempty\t\tX\nNA\tR\tX\nS2\tR\tX\n")
    return tmp_path / "groups", metadata


def test_real_qc_pca_cv_and_h5ad_preserve_identity(study, tmp_path):
    groups, metadata = study
    out = tmp_path / "analysis"
    original = (groups / "study_group_matrix.tsv").read_bytes()
    result = analyze_study(groups, out, metadata=metadata, replicate_column="replicates")
    data = ad.read_h5ad(out / "study.h5ad")
    assert data.obs_names.tolist() == ["NA", "S2", "S3", "empty"]
    assert data.obs.loc["NA", "replicates"] == "R"
    assert data.var_names.tolist() == ["A", "B", "C"]
    np.testing.assert_array_equal(data.obs.num_features_detected, [3, 2, 3, 0])
    assert data.X[0, 0] == 0.5 and data.layers["log2"][0, 0] == -1
    assert np.array_equal(np.isfinite(data.X), np.isfinite(data.layers["log2"]))
    assert data.obs.pca_included.tolist() == [True, True, True, False]
    assert data.var.pca_eligible.tolist() == [True, True, False]
    assert result["pca"]["status"] == "PASS"
    # Independent centered matrix reconstruction, insensitive to arbitrary PCA signs.
    x = data.layers["log2"][:3, :2]
    reconstructed = data.obsm["X_pca"][:3] @ data.varm["PCs"][:2].T
    np.testing.assert_allclose(reconstructed, x - x.mean(axis=0), rtol=1e-6, atol=1e-6)
    cv = pd.read_csv(out / "replicate_cv.tsv.gz", sep="\t").set_index("study_group")
    assert cv.loc["C", "measured_samples"] == 2 and not cv.loc["C", "complete"]
    assert cv.loc["C", "cv_sample_pct"] == pytest.approx(np.std([5, 20], ddof=1) / 12.5 * 100)
    assert cv.loc["C", "cv_population_pct"] == pytest.approx(60)
    assert (groups / "study_group_matrix.tsv").read_bytes() == original
    assert json.loads(data.uns["fastalake"])["method"] == "sum"
    assert (out / "QC_overview.pdf").stat().st_size > 1000


@pytest.mark.parametrize("defect,match", [("missing", "exactly"), ("duplicate", "unique")])
def test_bad_metadata_is_rejected(study, tmp_path, defect, match):
    groups, metadata = study
    s = metadata.read_text()
    metadata.write_text(s.replace("empty\t\tX\n", "") if defect == "missing" else s + "S3\tR\tX\n")
    with pytest.raises(ValueError, match=match):
        analyze_study(groups, tmp_path / "out", metadata=metadata)
    assert not (tmp_path / "out").exists()


def test_pca_feature_roster_cannot_impute_sparse_groups(study, tmp_path):
    groups, _ = study
    mask = tmp_path / "features.txt"
    mask.write_text("A\nC\n")
    with pytest.raises(ValueError, match="measured in every"):
        analyze_study(groups, tmp_path / "out", pca_features=mask)


def test_selected_source_checksum_mismatch_fails(study, tmp_path):
    groups, _ = study
    selected = tmp_path / "selected"
    quantify_study(groups, selected, method="sum")
    p = groups / "peptide_to_study_group.tsv"
    p.write_text(p.read_text() + "CORRUPTED\tA\n")
    with pytest.raises(ValueError, match="does not match"):
        analyze_study(groups, tmp_path / "out", quantification=selected)


def test_cli_reads_real_directlfq_output(study, tmp_path):
    pytest.importorskip("directlfq")
    groups, metadata = study
    selected = tmp_path / "lfq"
    quantify_study(groups, selected, method="directlfq")
    result = CliRunner().invoke(
        cli,
        [
            "analyze-study",
            "--groups",
            str(groups),
            "--quantification",
            str(selected),
            "--out",
            str(tmp_path / "out"),
            "--metadata",
            str(metadata),
        ],
    )
    assert result.exit_code == 0, result.output
    observed = ad.read_h5ad(tmp_path / "out/study.h5ad")
    matrix = pd.read_csv(
        selected / "study_group_matrix.tsv",
        sep="\t",
        index_col=0,
        keep_default_na=False,
        na_values=[""],
    ).astype(float)
    np.testing.assert_allclose(observed.X, matrix.to_numpy().T, equal_nan=True)
    assert json.loads((tmp_path / "out/summary.json").read_text())["method"] == "directlfq"


def test_all_empty_keeps_roster_and_skips_pca(tmp_path):
    search(tmp_path, "empty", {"A": P1}, [(P1, "A", 1, 0.001)], [])
    group_searches(tmp_path / "search", tmp_path / "groups")
    result = analyze_study(tmp_path / "groups", tmp_path / "out")
    assert result["groups"] == 0 and result["samples"] == 1
    assert result["pca"]["status"].startswith("SKIPPED")
    assert ad.read_h5ad(tmp_path / "out/study.h5ad").obs_names.tolist() == ["empty"]


@pytest.mark.parametrize("method", ["sum", "directlfq"])
def test_quantify_cli_writes_qc_and_pca_by_default(study, tmp_path, method):
    if method == "directlfq":
        pytest.importorskip("directlfq")
    groups, _ = study
    out = tmp_path / "automatic"
    result = CliRunner().invoke(
        cli, ["quantify-study", "--groups", str(groups), "--out", str(out), "--method", method]
    )
    assert result.exit_code == 0, result.output
    record = json.loads((out / "qc/summary.json").read_text())
    assert record["method"] == method
    assert record["pca"]["status"] == "PASS"
    for name in ("QC_overview.pdf", "QC_overview.png", "PCA.pdf", "PCA.png"):
        assert (out / "qc" / name).stat().st_size > 1000
    data = ad.read_h5ad(out / "qc/study.h5ad")
    assert data.obs_names.tolist() == ["NA", "S2", "S3", "empty"]
    assert np.isnan(data["S2", "C"].X).all()
    assert not data.var.loc["C", "pca_eligible"]
    assert json.loads((out / "WORKFLOW.json").read_text())["qc"]["status"] == "PASS"


def test_grouping_tool_writes_qc_by_default(study, tmp_path):
    import subprocess
    import sys
    from pathlib import Path

    out = tmp_path / "automatic_groups"
    tool = Path(__file__).resolve().parents[1] / "tools/group_study.py"
    result = subprocess.run(
        [sys.executable, str(tool), "--search", str(tmp_path / "search"), "--out", str(out)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    record = json.loads((out / "WORKFLOW.json").read_text())
    assert record["qc"]["pca"]["status"] == "PASS"
    assert (out / "qc/PCA.png").stat().st_size > 1000


def test_default_qc_dependency_failure_precedes_quantification(study, tmp_path, monkeypatch):
    from fasta_lake import downstream

    def unavailable():
        raise RuntimeError("analysis dependency unavailable")

    monkeypatch.setattr(downstream, "require_analysis", unavailable)
    out = tmp_path / "must_not_exist"
    result = CliRunner().invoke(
        cli, ["quantify-study", "--groups", str(study[0]), "--out", str(out)]
    )
    assert result.exit_code != 0
    assert "analysis dependency unavailable" in result.output
    assert not out.exists()
