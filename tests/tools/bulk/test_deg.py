"""Tests for bulk differential expression."""

import numpy as np
import pytest
from anndata import AnnData

from scLucid.tools.bulk import run_bulk_de


def test_run_bulk_de_welch():
    np.random.seed(0)
    n_genes = 50
    # Condition A: low expression
    a = np.random.poisson(10, size=(3, n_genes)).astype(float)
    # Condition B: high expression
    b = np.random.poisson(50, size=(3, n_genes)).astype(float)
    counts = np.vstack([a, b])
    obs = {"condition": ["A"] * 3 + ["B"] * 3}
    adata = AnnData(X=counts, obs=obs)

    result = run_bulk_de(
        adata,
        condition_col="condition",
        condition1="A",
        condition2="B",
        method="welch",
    )
    assert not result.empty
    assert "log2fc" in result.columns
    assert "pvals_adj" in result.columns
    assert result["valid_for_publication_inference"].all()
    assert result["inference_level"].iloc[0] == "sample_level"
    assert "bulk" in adata.uns["sclucid"]["tools"]
    assert adata.uns["sclucid"]["tools"]["bulk"]["de"]["method"] == "welch"


def test_run_bulk_de_insufficient_replicates_raises():
    counts = np.random.poisson(10, size=(2, 10)).astype(float)
    obs = {"condition": ["A", "B"]}
    adata = AnnData(X=counts, obs=obs)
    with pytest.raises(ValueError):
        run_bulk_de(
            adata,
            condition_col="condition",
            condition1="A",
            condition2="B",
            method="welch",
        )


def test_run_bulk_de_descriptive_fallback():
    counts = np.random.poisson(10, size=(2, 10)).astype(float)
    obs = {"condition": ["A", "B"]}
    adata = AnnData(X=counts, obs=obs)
    result = run_bulk_de(
        adata,
        condition_col="condition",
        condition1="A",
        condition2="B",
        method="welch",
        fallback_to_descriptive=True,
    )
    assert not result.empty
    assert not result["valid_for_publication_inference"].iloc[0]
    assert result["inference_level"].iloc[0] == "descriptive_sample_level"


def test_run_bulk_de_ignores_unrequested_condition_levels():
    a = np.array([[10.0, 90.0, 90.0], [12.0, 88.0, 90.0]])
    b = np.array([[90.0, 10.0, 90.0], [88.0, 12.0, 90.0]])
    c = np.array([[1.0, 1000.0, 1.0], [1.0, 1000.0, 1.0]])
    counts = np.vstack([a, b, c])
    obs = {"condition": ["A"] * 2 + ["B"] * 2 + ["C"] * 2}
    adata = AnnData(X=counts, obs=obs)

    result = run_bulk_de(
        adata,
        condition_col="condition",
        condition1="A",
        condition2="B",
        method="welch",
        min_counts_per_gene=0,
    )

    assert not result.empty
    assert result["n_samples_condition1"].iloc[0] == 2
    assert result["n_samples_condition2"].iloc[0] == 2

    b_cpm = b / b.sum(axis=1, keepdims=True) * 1e6
    expected_b_gene0_mean = np.log2(b_cpm[:, 0] + 1.0).mean()
    observed_b_gene0_mean = result.loc[
        result["gene"] == "0", "mean_logcpm_condition2"
    ].iloc[0]
    assert np.isclose(observed_b_gene0_mean, expected_b_gene0_mean)
