"""Tests for bulk differential expression."""

import numpy as np
import pytest
from anndata import AnnData

from scLucid.tools.bulk import run_bulk_de, BulkDEConfig


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
    assert result["valid_for_publication_inference"].iloc[0] == False
    assert result["inference_level"].iloc[0] == "descriptive_sample_level"
