"""Tests for cell-type proportion statistics."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from pydantic import ValidationError

from scLucid.analysis.config import ProportionConfig as PublicProportionConfig
from scLucid.analysis.proportion.config import ProportionConfig
from scLucid.analysis.proportion.pseudobulk import celltype_proportion_analysis
from scLucid.analysis.proportion.stats import composition_transform, run_statistical_test


def test_proportion_config_accepts_kruskal():
    cfg = ProportionConfig(
        celltype_col="cell_type",
        sample_col="sample",
        condition_col="condition",
        test_method="kruskal",
    )

    assert cfg.test_method == "kruskal"


def test_proportion_config_defaults_to_clr_sample_level():
    cfg = ProportionConfig(
        celltype_col="cell_type",
        sample_col="sample",
        condition_col="condition",
    )

    assert cfg.test_method == "clr-t-test"
    assert cfg.composition_transform == "clr"


def test_proportion_config_rejects_unimplemented_fisher():
    with pytest.raises(ValidationError):
        ProportionConfig(
            celltype_col="cell_type",
            sample_col="sample",
            condition_col="condition",
            test_method="fisher",
        )


def test_public_proportion_config_rejects_unimplemented_fisher():
    with pytest.raises(ValidationError):
        PublicProportionConfig(
            celltype_col="cell_type",
            sample_col="sample",
            condition_col="condition",
            test_method="fisher",
        )


def test_deseq2_proportion_backend_returns_standard_columns():
    pytest.importorskip("pydeseq2")
    count_df = pd.DataFrame(
        {
            "T": [40, 42, 38, 180, 176, 184],
            "B": [180, 176, 184, 40, 42, 38],
        },
        index=["a1", "a2", "a3", "b1", "b2", "b3"],
    )
    sample_to_cond = pd.Series(["A", "A", "A", "B", "B", "B"], index=count_df.index)

    result = run_statistical_test(
        count_df,
        condition_col="condition",
        test_method="deseq2",
        sample_to_cond=sample_to_cond,
        multiple_testing_correction=None,
    )

    assert set(result["cell_type"]) == {"T", "B"}
    assert {"log2fc", "statistic", "pval", "padj", "direction", "method"}.issubset(
        result.columns
    )
    assert set(result["method"]) == {"deseq2"}
    assert result[result["cell_type"] == "T"]["log2fc"].iloc[0] > 0


def test_two_group_statistics_use_condition2_minus_condition1_direction():
    count_df = pd.DataFrame(
        {"T": [0.2, 0.3, 0.8, 0.9]},
        index=["s1", "s2", "s3", "s4"],
    )
    sample_to_cond = pd.Series(["A", "A", "B", "B"], index=count_df.index)

    result = run_statistical_test(
        count_df,
        condition_col="condition",
        test_method="t-test",
        sample_to_cond=sample_to_cond,
        multiple_testing_correction=None,
    )

    row = result.iloc[0]
    assert row["condition1"] == "A"
    assert row["condition2"] == "B"
    assert row["direction"] == "B - A"
    assert row["mean_diff"] == pytest.approx(0.6)
    assert row["inference_level"] == "exploratory_legacy_proportion"
    assert bool(row["valid_for_publication_inference"]) is False
    assert "compositional" in row["compositional_data_warning"]


def test_clr_ttest_reports_compositional_effect_ci_and_fdr():
    prop_df = pd.DataFrame(
        {
            "T": [0.7, 0.65, 0.2, 0.25],
            "B": [0.2, 0.25, 0.6, 0.55],
            "NK": [0.1, 0.1, 0.2, 0.2],
        },
        index=["a1", "a2", "b1", "b2"],
    )
    sample_to_cond = pd.Series(["A", "A", "B", "B"], index=prop_df.index)

    clr = composition_transform(prop_df, pseudocount=0.01)
    result = run_statistical_test(
        prop_df,
        condition_col="condition",
        test_method="clr-t-test",
        sample_to_cond=sample_to_cond,
        composition_pseudocount=0.01,
    )

    assert clr.mean(axis=1).abs().max() < 1e-10
    assert {"effect_size", "ci_lower", "ci_upper", "padj", "inference_level"}.issubset(
        result.columns
    )
    assert set(result["method"]) == {"clr-t-test"}
    assert set(result["inference_level"]) == {"sample_level"}


def test_pseudobulk_single_replicate_proportion_marked_descriptive():
    adata = AnnData(X=np.ones((6, 1)))
    adata.obs["sample"] = ["s1", "s1", "s1", "s2", "s2", "s2"]
    adata.obs["condition"] = ["A", "A", "A", "B", "B", "B"]
    adata.obs["cell_type"] = ["T", "T", "B", "T", "B", "B"]

    cfg = ProportionConfig(
        celltype_col="cell_type",
        sample_col="sample",
        condition_col="condition",
        test_method="clr-t-test",
        auto_configure=True,
        plot_types=[],
    )
    _, stat_df = celltype_proportion_analysis(adata, cfg)

    assert not stat_df.empty
    assert set(stat_df["inference_level"]) == {"descriptive_sample_level"}
    assert not stat_df["valid_for_publication_inference"].any()
    assert stat_df["pval"].isna().all()


def test_chi_square_returns_global_pvalue_and_celltype_contributions():
    count_df = pd.DataFrame(
        {"T": [10, 12, 2, 3], "B": [2, 1, 10, 12]},
        index=["s1", "s2", "s3", "s4"],
    )
    sample_to_cond = pd.Series(["A", "A", "B", "B"], index=count_df.index)

    result = run_statistical_test(
        count_df,
        condition_col="condition",
        test_method="chi-square",
        sample_to_cond=sample_to_cond,
        multiple_testing_correction=None,
    )

    assert set(result["cell_type"]) == {"T", "B"}
    assert {"overall_pval", "method_note", "observed_A", "expected_A"}.issubset(result.columns)
    assert result["statistic"].gt(0).all()


def test_run_statistical_test_kruskal():
    count_df = pd.DataFrame(
        {
            "T": [0.1, 0.2, 0.8, 0.9, 0.4, 0.5],
            "B": [0.6, 0.5, 0.2, 0.1, 0.3, 0.2],
        },
        index=["s1", "s2", "s3", "s4", "s5", "s6"],
    )
    sample_to_cond = pd.Series(
        ["A", "A", "B", "B", "C", "C"],
        index=count_df.index,
    )

    result = run_statistical_test(
        count_df,
        condition_col="condition",
        test_method="kruskal",
        sample_to_cond=sample_to_cond,
    )

    assert set(result["cell_type"]) == {"T", "B"}
    assert {"statistic", "pval", "padj"}.issubset(result.columns)


def test_pseudobulk_kruskal_uses_sample_level_metadata():
    adata = AnnData(X=np.ones((6, 1)))
    adata.obs["sample"] = ["s1", "s1", "s2", "s2", "s3", "s3"]
    adata.obs["condition"] = ["A", "A", "B", "B", "C", "C"]
    adata.obs["cell_type"] = ["T", "B", "T", "T", "B", "B"]

    cfg = ProportionConfig(
        celltype_col="cell_type",
        sample_col="sample",
        condition_col="condition",
        test_method="kruskal",
        auto_configure=True,
        plot_types=[],
    )

    prop_df, stat_df = celltype_proportion_analysis(adata, cfg)

    assert list(prop_df.index) == ["s1", "s2", "s3"]
    assert not stat_df.empty
    assert set(stat_df["cell_type"]) == {"T", "B"}


def test_pseudobulk_two_group_plots_use_sample_level_metadata():
    adata = AnnData(X=np.ones((8, 1)))
    adata.obs["sample"] = ["s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4"]
    adata.obs["condition"] = ["A", "A", "A", "A", "B", "B", "B", "B"]
    adata.obs["patient"] = ["p1", "p1", "p2", "p2", "p1", "p1", "p2", "p2"]
    adata.obs["cell_type"] = ["T", "B", "T", "T", "B", "B", "T", "B"]

    cfg = ProportionConfig(
        celltype_col="cell_type",
        sample_col="sample",
        condition_col="condition",
        pairing_col="patient",
        test_method="t-test",
        auto_configure=False,
        plot_types=["box", "diff", "shift", "paired_shift", "composition_pca", "clr_heatmap"],
    )

    prop_df, stat_df = celltype_proportion_analysis(adata, cfg)

    assert list(prop_df.index) == ["s1", "s2", "s3", "s4"]
    assert "effect_size_cohens_d" in stat_df.columns
