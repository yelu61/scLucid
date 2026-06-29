"""Tests for scoring helpers."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scLucid.analysis import (
    calculate_signature_matrix,
    plot_delta_heatmap,
    plot_score_violin_with_stats,
    run_module_scoring_workflow,
    score_by_gene_sets,
)


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


def test_score_by_gene_sets_accepts_legacy_log1p_norm_alias(scoring_adata):
    scoring_adata.raw = None
    scoring_adata.layers["log1p_norm"] = np.log1p(scoring_adata.X.copy())

    result = score_by_gene_sets(
        scoring_adata,
        {"T_core": ["CD3D", "CD3E"]},
        layer="normalized",
        use_raw=False,
    )

    assert "T_core_score" in result.obs.columns


def test_calculate_signature_matrix_falls_back_when_raw_missing(scoring_adata):
    scoring_adata.raw = None
    matrix = calculate_signature_matrix(
        scoring_adata,
        {"T_core": ["CD3D", "CD3E"], "NK_core": ["NKG7", "GNLY"]},
        groupby="celltype",
        use_raw=True,
        z_score=False,
    )

    assert set(matrix.index) == {"T_core", "NK_core"}
    assert set(matrix.columns.astype(str)) == {"T", "NK"}


def test_calculate_signature_matrix_constant_zscore_has_no_nan():
    adata = AnnData(np.ones((4, 2), dtype=float))
    adata.var_names = ["GeneA", "GeneB"]
    adata.obs["celltype"] = pd.Categorical(["A", "A", "B", "B"])

    matrix = calculate_signature_matrix(
        adata,
        {"constant": ["GeneA", "GeneB"]},
        groupby="celltype",
        use_raw=True,
        z_score=True,
    )

    assert not matrix.isna().any().any()


def test_plot_delta_heatmap_reports_no_common_groups(scoring_adata):
    with pytest.raises(ValueError, match="No common"):
        plot_delta_heatmap(
            scoring_adata,
            {"T_core": ["CD3D", "CD3E"]},
            groupby="celltype",
            compare_group="condition",
            ref_group="ctrl",
            target_group="tx",
            use_raw=True,
        )


def test_plot_score_violin_with_stats_validates_input(scoring_adata):
    scoring_adata.obs["toy_score"] = [1.0, 1.1, 2.0, 2.1]

    with pytest.raises(KeyError, match="missing_score"):
        plot_score_violin_with_stats(
            scoring_adata,
            score_key="missing_score",
            groupby="condition",
            group1="ctrl",
            group2="tx",
        )

    with pytest.raises(ValueError, match="at least 2 observations"):
        plot_score_violin_with_stats(
            scoring_adata,
            score_key="toy_score",
            groupby="sample_id",
            group1="S1",
            group2="missing",
        )
