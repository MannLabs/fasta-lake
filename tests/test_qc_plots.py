"""Check the measured distribution and missingness shown to the reader."""

import numpy as np
import pytest

from fasta_lake.qc_plots import intensity_distribution_summary, plot_intensity_distributions

ad = pytest.importorskip("anndata")
pd = pytest.importorskip("pandas")


def test_distribution_retains_subunit_values_empty_samples_and_constant_samples(tmp_path):
    data = ad.AnnData(
        np.array([[0.5, 1, 2], [4, 4, 4], [0, np.nan, 0]]),
        obs=pd.DataFrame(index=["negative log2", "constant", "empty"]),
    )
    result = intensity_distribution_summary(data)
    assert result.index.tolist() == list(data.obs_names)
    assert result.loc["negative log2", "log2_min"] == -1
    assert result.loc["negative log2", "log2_median"] == 0
    assert result.loc["constant", "log2_median"] == 2
    assert result.loc["empty", "missing_fraction"] == 1
    assert np.isnan(result.loc["empty", "log2_median"])
    before = data.X.copy()
    plot_intensity_distributions(data, tmp_path, "sum", per_page=2)
    np.testing.assert_array_equal(data.X, before)
    assert (tmp_path / "Protein_intensity_violins_page2.png").stat().st_size > 1000
    saved = pd.read_csv(tmp_path / "intensity_distribution.tsv", sep="\t", index_col=0)
    pd.testing.assert_frame_equal(saved, result)


def test_single_observation_is_not_a_density(tmp_path):
    data = ad.AnnData(np.array([[2.0]]), obs=pd.DataFrame(index=["single"]))
    plot_intensity_distributions(data, tmp_path, "directlfq")
    assert (tmp_path / "Protein_intensity_violins.pdf").stat().st_size > 1000


def test_zero_groups_keeps_empty_acquisition(tmp_path):
    data = ad.AnnData(np.empty((1, 0)), obs=pd.DataFrame(index=["empty"]))
    result = intensity_distribution_summary(data)
    assert result.loc["empty", "total_groups"] == 0
    assert np.isnan(result.loc["empty", "missing_fraction"])
    plot_intensity_distributions(data, tmp_path, "sum")
