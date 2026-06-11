"""Tests for bulk deconvolution backward compatibility."""

import numpy as np
import pandas as pd
from anndata import AnnData

from scLucid.tools.bulk import deconvolve_bulk
from scLucid.analysis.bulk import deconvolve_bulk as legacy_deconvolve_bulk


def _make_reference(n_cells=20, n_genes=200):
    counts = np.random.poisson(5, size=(n_cells, n_genes)).astype(float)
    obs = {
        "cell_type": ["T"] * 10 + ["B"] * 10,
        "sampleID": ["S1"] * 20,
    }
    var_names = [f"G{i}" for i in range(n_genes)]
    import pandas as pd
    return AnnData(
        X=counts,
        obs=obs,
        var=pd.DataFrame({"gene": var_names}, index=var_names),
    )


def test_deconvolve_bulk_legacy_namespace():
    ref = _make_reference()
    bulk = pd.DataFrame(
        np.random.poisson(100, size=(200, 3)).astype(float),
        index=[f"G{i}" for i in range(200)],
        columns=["S1", "S2", "S3"],
    )
    out = legacy_deconvolve_bulk(ref, bulk, cell_type_key="cell_type", method="DWLS")
    assert "bulk_deconvolution" in out.uns["sclucid"]["tools"]


def test_deconvolve_bulk_analysis_namespace():
    ref = _make_reference()
    bulk = pd.DataFrame(
        np.random.poisson(100, size=(200, 3)).astype(float),
        index=[f"G{i}" for i in range(200)],
        columns=["S1", "S2", "S3"],
    )
    out = deconvolve_bulk(ref, bulk, cell_type_key="cell_type", method="DWLS")
    assert "deconvolution" in out.uns["sclucid"]["analysis"]["bulk"]
    assert "proportions" in out.uns["sclucid"]["analysis"]["bulk"]["deconvolution"]
