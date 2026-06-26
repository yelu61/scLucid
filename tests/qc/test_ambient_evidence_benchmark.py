"""Tests for ambient evidence benchmark runner."""

from pathlib import Path

import numpy as np
import pytest
from anndata import AnnData

from validation.qc.run_ambient_evidence_benchmark import _available, _ensure_metrics


def test_ensure_metrics_adds_barcode_metrics():
    X = np.random.poisson(2, (20, 10)).astype(float)
    adata = AnnData(X)
    adata.layers["counts"] = X.copy()
    _ensure_metrics(adata)
    assert "barcode_total_counts" in adata.obs
    assert "n_genes_by_counts" in adata.obs


def test_available_reports_backend_status():
    # These are optional dependencies; the function should just return bools.
    assert isinstance(_available("cellbender"), bool)
    assert isinstance(_available("soupx"), bool)
    assert isinstance(_available("scar"), bool)


@pytest.mark.parametrize("backend", ["cellbender", "soupx", "scar"])
def test_unknown_backend_returns_false(backend):
    assert isinstance(_available(backend), bool)
