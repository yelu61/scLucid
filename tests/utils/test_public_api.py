"""Public API checks for scLucid.utils."""

import numpy as np
import pytest
from anndata import AnnData

import scLucid.utils as utils


@pytest.mark.unit
def test_utils_exports_resolve():
    for symbol in utils.__all__:
        assert hasattr(utils, symbol), f"scLucid.utils missing exported symbol: {symbol}"


@pytest.mark.unit
def test_storage_roundtrip_via_public_api():
    adata = AnnData(np.random.randn(6, 4))
    utils.save_result(adata, "qc", "unit_test_payload", {"ok": True})
    loaded = utils.load_result(adata, "qc", "unit_test_payload")
    assert loaded["ok"] is True
