"""Tests for bulk normalization."""

import numpy as np
import pytest
from anndata import AnnData

from scLucid.tools.bulk import normalize_bulk_counts, estimate_size_factors_median_ratio


def test_normalize_bulk_cpm():
    counts = np.array([[1, 2, 0], [4, 0, 1]], dtype=float)
    adata = AnnData(X=counts)
    result = normalize_bulk_counts(adata)
    np.testing.assert_allclose(result.sum(axis=1).values, [1e6, 1e6], rtol=1e-6)


def test_normalize_bulk_rejects_negative():
    counts = np.array([[1, -1, 0]], dtype=float)
    adata = AnnData(X=counts)
    with pytest.raises(ValueError):
        normalize_bulk_counts(adata)


def test_estimate_size_factors_median_ratio():
    counts = np.array([[10, 20, 30], [40, 50, 60]], dtype=float)
    adata = AnnData(X=counts)
    sf = estimate_size_factors_median_ratio(adata)
    assert sf.shape == (2,)
    assert np.all(sf > 0)
