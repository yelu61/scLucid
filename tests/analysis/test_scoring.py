"""Tests for scoring helpers."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scLucid.analysis import run_module_scoring_workflow


@pytest.fixture
def scoring_adata():
    x = np.array(
        [
            [5, 4, 0, 0, 1, 1],
            [4, 5, 0, 0, 1, 1],
            [0, 1, 5, 4, 2, 2],
            [0, 0, 4, 5, 2, 2],
        ],
        dtype=float,
    )
    adata = AnnData(x)
    adata.var_names = ["CD3D", "CD3E", "NKG7", "GNLY", "ACTB", "GAPDH"]
    adata.obs_names = [f"cell_{idx}" for idx in range(adata.n_obs)]
    adata.obs["celltype"] = pd.Categorical(["T", "T", "NK", "NK"])
    adata.obs["sample_id"] = pd.Categorical(["S1", "S1", "S2", "S2"])
    adata.obs["condition"] = pd.Categorical(["ctrl", "ctrl", "tx", "tx"])
    adata.raw = adata.copy()
    return adata


@pytest.mark.unit
def test_run_module_scoring_workflow_returns_summary_tables(scoring_adata):
    modules = {
        "T_core": ["CD3D", "CD3E", "MISSING1"],
        "NK_core": ["NKG7", "GNLY"],
    }

    scored, results = run_module_scoring_workflow(
        scoring_adata,
        modules,
        groupby="celltype",
        sample_col="sample_id",
        condition_col="condition",
        use_raw=True,
    )

    assert "T_core_score" in scored.obs.columns
    assert "NK_core_score" in scored.obs.columns
    assert set(results) >= {
        "module_summary",
        "group_mean_scores",
        "sample_mean_scores",
        "condition_mean_scores",
    }
    assert results["module_summary"]["scored"].all()
    assert set(results["group_mean_scores"]["celltype"].astype(str)) == {"T", "NK"}


@pytest.mark.unit
def test_run_module_scoring_workflow_tracks_unscored_modules(scoring_adata):
    modules = {
        "valid": ["CD3D", "CD3E"],
        "invalid": ["DOES_NOT_EXIST"],
    }

    _, results = run_module_scoring_workflow(
        scoring_adata,
        modules,
        use_raw=True,
        min_genes_required=2,
    )

    summary = results["module_summary"].set_index("module")
    assert bool(summary.loc["valid", "scored"]) is True
    assert bool(summary.loc["invalid", "scored"]) is False
