"""Tests for bulk diagnostics."""

import numpy as np
import pytest
from anndata import AnnData

from scLucid.tools.bulk import diagnose_bulk_data_quality


def test_diagnose_bulk_data_quality_passes_with_replicates():
    n_samples, n_genes = 6, 100
    counts = np.random.poisson(100, size=(n_samples, n_genes)).astype(float)
    obs = {
        "condition": ["A", "A", "A", "B", "B", "B"],
    }
    adata = AnnData(X=counts, obs=obs)
    result = diagnose_bulk_data_quality(adata, condition_col="condition")
    assert result["replicate_requirement_met"]
    assert result["n_samples"] == n_samples
    assert result["recommended_method"] == "welch"


def test_diagnose_bulk_data_quality_warns_on_insufficient_replicates():
    counts = np.random.poisson(100, size=(3, 100)).astype(float)
    obs = {"condition": ["A", "B", "B"]}
    adata = AnnData(X=counts, obs=obs)
    result = diagnose_bulk_data_quality(adata, condition_col="condition")
    assert not result["replicate_requirement_met"]
    assert result["recommended_method"] == "descriptive"
    assert any("replicate" in w for w in result["warnings"])
