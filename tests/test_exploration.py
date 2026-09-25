"""The public functional route must preserve denominators and experimental design."""

import json

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from fasta_lake.exploration import explore_study, read_group_terms
from fasta_lake.exploration_cli import explore_cmd

pytest.importorskip("alphapepttools")
pytest.importorskip("anndata")


@pytest.fixture
def inputs(tmp_path):
    names = [f"donor{d}_{condition}" for d in range(3) for condition in ("A", "B")]
    matrix = pd.DataFrame(
        [[4, 7, 5, 8, 6, 9], [9, 5, 11, 7, 10, 8], [3, 6, 4, 5, 5, 7], [100, 20, 80, 40, 90, 30]],
        index=["g1", "g2", "g3", "unknown"],
        columns=names,
    )
    matrix.to_csv(tmp_path / "matrix.tsv", sep="\t", index_label="study_group")
    (tmp_path / "terms.tsv").write_text("study_group\tterm\ng1\tK1\ng2\tK2\ng3\tK3\n")
    metadata = pd.DataFrame(
        {
            "sample": names,
            "condition": ["A", "B"] * 3,
            "donor": [f"donor{d}" for d in range(3) for _ in range(2)],
        }
    )
    metadata.to_csv(tmp_path / "metadata.tsv", sep="\t", index=False)
    return tmp_path, matrix


def test_installed_cli_functional_route_and_conditional_denominator(inputs):
    root, matrix = inputs
    out = root / "result"
    args = [
        "--matrix",
        str(root / "matrix.tsv"),
        "--annotations",
        str(root / "terms.tsv"),
        "--out",
        str(out),
        "--level",
        "KO",
        "--metadata",
        str(root / "metadata.tsv"),
        "--group-column",
        "condition",
        "--block-column",
        "donor",
        "--permutations",
        "19",
    ]
    result = CliRunner().invoke(explore_cmd, args)
    assert result.exit_code == 0, result.output
    receipt = json.loads((out / "COMPLETE.json").read_text())
    assert receipt["experimental"] and receipt["outputs"]["pca"] == "PASS"
    profile = pd.read_csv(out / "term_signal.tsv", sep="\t", index_col=0)
    assert np.allclose(profile.sum(axis=0), matrix.sum(axis=0))
    diversity = pd.read_csv(out / "diversity.tsv", sep="\t", index_col=0)
    assert (diversity.observed_richness == 3).all()  # unknown is never a pseudo-taxon
    coverage = pd.read_csv(out / "annotation_coverage.tsv", sep="\t", index_col=0)
    assert np.allclose(
        coverage.annotation_fraction,
        matrix.loc[["g1", "g2", "g3"]].sum(axis=0) / matrix.sum(axis=0),
    )
    clr = pd.read_csv(out / "clr_values.tsv", sep="\t", index_col=0)
    assert np.allclose(clr.sum(axis=1), 0, atol=1e-12)
    assert "__unannotated__" not in clr.columns
    assert (out / "Functional_QC.png").stat().st_size > 1000
    assert (out / "PCA_clr.pdf").stat().st_size > 1000
    test = json.loads((out / "PERMANOVA.json").read_text())
    assert test["blocked"] and test["permutations"] == 19


def test_all_unannotated_has_coverage_and_explicit_pca_exclusion(inputs):
    root, _ = inputs
    (root / "terms.tsv").write_text("study_group\tterm\n")
    result = explore_study(root / "matrix.tsv", root / "terms.tsv", root / "out", level="KO")
    assert result["outputs"]["pca"].startswith("EXCLUDED")
    assert not (root / "out/diversity.tsv").exists()
    coverage = pd.read_csv(root / "out/annotation_coverage.tsv", sep="\t")
    assert (coverage.annotation_fraction == 0).all()
    assert not (root / "out/PERMANOVA.json").exists()


def test_foreign_metadata_rejected_before_output(inputs):
    root, _ = inputs
    p = root / "metadata.tsv"
    p.write_text(p.read_text().replace("donor0_A", "foreign"))
    with pytest.raises(ValueError, match="exactly the matrix acquisitions"):
        explore_study(root / "matrix.tsv", root / "terms.tsv", root / "out", level="KO", metadata=p)
    assert not (root / "out").exists()


def test_unknown_group_cannot_be_silently_ignored(inputs):
    root, _ = inputs
    (root / "terms.tsv").write_text("study_group\tterm\nforeign\tK1\n")
    with pytest.raises(ValueError, match="absent from the study"):
        explore_study(root / "matrix.tsv", root / "terms.tsv", root / "out", level="KO")


def test_duplicate_term_pairs_are_deduplicated(inputs):
    root, _ = inputs
    (root / "terms.tsv").write_text("study_group\tterm\ng1\tK1\ng1\tK1\n")
    assert read_group_terms(root / "terms.tsv") == {"g1": {"K1"}}
