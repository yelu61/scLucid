"""Tests for tumor-focused bulk utilities."""

import numpy as np
import pandas as pd
from anndata import AnnData

from scLucid.tools.bulk import (
    associate_tme_with_response,
    bulk_immune_landscape,
    deconvolve_tumor_tme,
    estimate_tumor_purity_from_bulk,
)


def _make_reference(n_cells=30, n_genes=200):
    counts = np.random.poisson(5, size=(n_cells, n_genes)).astype(float)
    obs = {
        "cell_type": ["T_cell"] * 10 + ["B_cell"] * 10 + ["fibroblast"] * 10,
        "sampleID": ["S1"] * 30,
    }
    var_names = [f"G{i}" for i in range(n_genes)]
    return AnnData(
        X=counts,
        obs=obs,
        var=pd.DataFrame({"gene": var_names}, index=var_names),
    )


def test_deconvolve_tumor_tme():
    ref = _make_reference()
    bulk = pd.DataFrame(
        np.random.poisson(100, size=(200, 4)).astype(float),
        index=[f"G{i}" for i in range(200)],
        columns=["A1", "A2", "B1", "B2"],
    )
    out = deconvolve_tumor_tme(
        ref, bulk, cell_type_key="cell_type", method="DWLS", key_added="tumor_tme"
    )
    assert "tumor_tme" in out.uns["sclucid"]["tools"]
    record = out.uns["sclucid"]["tools"]["tumor_tme"]
    assert "compartment_proportions" in record
    assert all(c in record["compartment_proportions"].columns for c in ["malignant", "immune", "stromal", "other"])
    assert "result_warning" in record


def test_estimate_tumor_purity_from_bulk_expression():
    bulk = pd.DataFrame(
        np.random.poisson(50, size=(20, 5)).astype(float),
        index=[f"G{i}" for i in range(20)],
        columns=["S1", "S2", "S3", "S4", "S5"],
    )
    result = estimate_tumor_purity_from_bulk(bulk, method="expression")
    assert not result.empty
    assert "tumor_purity" in result.columns
    assert not result["valid_for_publication_inference"].any()


def test_bulk_immune_landscape_mean_expression():
    bulk = pd.DataFrame(
        np.random.poisson(50, size=(50, 4)).astype(float),
        index=[f"G{i}" for i in range(50)],
        columns=["S1", "S2", "S3", "S4"],
    )
    result = bulk_immune_landscape(bulk, method="mean_expression")
    assert not result.empty
    assert any("T_cells" in c for c in result.columns)
    assert not result["valid_for_publication_inference"].any()


def test_associate_tme_with_response_mannwhitney():
    np.random.seed(0)
    proportions = pd.DataFrame(
        {
            "immune": np.concatenate([np.random.uniform(0.1, 0.3, 5), np.random.uniform(0.4, 0.6, 5)]),
            "stromal": np.random.uniform(0.1, 0.4, 10),
        },
        index=[f"S{i}" for i in range(10)],
    )
    metadata = pd.DataFrame(
        {"response": ["NR"] * 5 + ["R"] * 5},
        index=[f"S{i}" for i in range(10)],
    )
    result = associate_tme_with_response(
        proportions, metadata, response_col="response", method="mannwhitney"
    )
    assert not result.empty
    assert "pval" in result.columns
    assert "pvals_adj" in result.columns
