"""Tests for cell-type proportion statistics."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from pydantic import ValidationError

from scLucid.analysis.config import ProportionConfig as PublicProportionConfig
from scLucid.analysis.proportion.config import ProportionConfig
from scLucid.analysis.proportion.pseudobulk import (
    _auto_configure_analysis,
    celltype_proportion_analysis,
)
from scLucid.analysis.proportion.stats import (
    _run_ancom_like_clr_test,
    composition_transform,
    run_statistical_test,
)


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


def test_proportion_auto_configure_false_preserves_method_and_empty_plot_plan():
    adata = AnnData(X=np.ones((24, 1)))
    adata.obs["sample"] = np.repeat([f"s{i}" for i in range(6)], 4)
    adata.obs["condition"] = adata.obs["sample"].map(
        {"s0": "A", "s1": "A", "s2": "A", "s3": "B", "s4": "B", "s5": "B"}
    )
    config = ProportionConfig(
        celltype_col="cell_type",
        sample_col="sample",
        condition_col="condition",
        auto_configure=False,
        test_method="clr-t-test",
        plot_types=[],
    )

    resolved = _auto_configure_analysis(adata, config)

    assert resolved.test_method == "clr-t-test"
    assert resolved.plot_types == []


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
        legacy_exploratory=True,
    )

    row = result.iloc[0]
    assert row["condition1"] == "A"
    assert row["condition2"] == "B"
    assert row["direction"] == "B - A"
    assert row["mean_diff"] == pytest.approx(0.6)
    assert row["inference_level"] == "exploratory_legacy_proportion"
    assert row["claim_level"] == "exploratory_hypothesis_generation"
    assert row["model_type"] == "raw_proportion_legacy_test"
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
    assert set(result["claim_level"]) == {"sample_level_clr_compositional_inference"}
    assert set(result["model_type"]) == {"sample_level_clr_test"}


def test_composition_transform_closes_count_input():
    count_df = pd.DataFrame({"T": [10, 5], "B": [30, 15], "NK": [10, 30]})

    clr = composition_transform(count_df, pseudocount=0.01)
    expected = composition_transform(
        count_df.div(count_df.sum(axis=1), axis=0),
        pseudocount=0.01,
    )

    pd.testing.assert_frame_equal(clr, expected)
    assert clr.mean(axis=1).abs().max() < 1e-10


def test_composition_transform_keeps_closed_proportions():
    prop_df = pd.DataFrame({"T": [0.2, 0.5], "B": [0.3, 0.25], "NK": [0.5, 0.25]})

    clr = composition_transform(prop_df, pseudocount=0.01)
    manual = np.log(prop_df + 0.01)
    manual = manual.sub(manual.mean(axis=1), axis=0)

    pd.testing.assert_frame_equal(clr, manual)


def test_composition_transform_converts_percent_input():
    percent_df = pd.DataFrame({"T": [20.0, 50.0], "B": [30.0, 25.0], "NK": [50.0, 25.0]})
    prop_df = percent_df / 100.0

    clr = composition_transform(percent_df, pseudocount=0.01)
    expected = composition_transform(prop_df, pseudocount=0.01)

    pd.testing.assert_frame_equal(clr, expected)


def test_composition_transform_closes_subcomposition_input():
    sub_df = pd.DataFrame({"T": [0.2, 0.1], "B": [0.3, 0.2], "NK": [0.1, 0.2]})

    clr = composition_transform(sub_df, pseudocount=0.01)
    expected = composition_transform(sub_df.div(sub_df.sum(axis=1), axis=0), pseudocount=0.01)

    pd.testing.assert_frame_equal(clr, expected)


def test_composition_transform_rejects_negative_values():
    with pytest.raises(ValueError, match="non-negative"):
        composition_transform(pd.DataFrame({"T": [0.5], "B": [-0.1]}))


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
    assert set(stat_df["claim_level"]) == {"descriptive_effect_size_only"}
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
    assert set(result["claim_level"]) == {"exploratory_global_composition_screen"}
    assert set(result["model_type"]) == {"global_chi_square_contingency_contribution"}
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
        legacy_exploratory=True,
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
        legacy_exploratory=True,
        plot_types=[],
    )

    prop_df, stat_df = celltype_proportion_analysis(adata, cfg)

    assert list(prop_df.index) == ["s1", "s2", "s3"]
    assert not stat_df.empty
    assert set(stat_df["cell_type"]) == {"T", "B"}


def test_recommend_sccoda_reference_selects_stable_abundant_cell_type():
    from scLucid.analysis.proportion.sccoda import recommend_sccoda_reference

    count_df = pd.DataFrame(
        {
            "Stable": [500, 520, 480, 510],
            "Rare": [10, 15, 8, 12],
            "Variable": [100, 300, 50, 250],
        },
        index=["s1", "s2", "s3", "s4"],
    )

    ref, diag = recommend_sccoda_reference(count_df)

    assert ref == "Stable"
    assert {"mean_frac", "cv"}.issubset(diag.columns)
    assert diag.loc["Stable", "mean_frac"] > diag.loc["Rare", "mean_frac"]


def test_ancom_like_clr_test_returns_expected_columns_and_tags():
    count_df = pd.DataFrame(
        {
            "T": [70, 65, 20, 25],
            "B": [20, 25, 60, 55],
            "NK": [10, 10, 20, 20],
        },
        index=["a1", "a2", "b1", "b2"],
    )
    metadata_df = pd.DataFrame(
        {
            "sample": ["a1", "a2", "b1", "b2"],
            "condition": ["A", "A", "B", "B"],
        }
    )

    result = _run_ancom_like_clr_test(
        count_df,
        metadata_df,
        condition_col="condition",
        sample_col="sample",
    )

    assert set(result["cell_type"]) == {"T", "B", "NK"}
    assert {"pval", "padj", "w_statistic", "significant", "method"}.issubset(result.columns)
    assert set(result["method"]) == {"ancom-like-clr"}
    assert set(result["inference_level"]) == {"sample_level"}
    assert set(result["claim_level"]) == {"ancom_like_clr_heuristic"}
    assert "ANCOM-BC" in result["proportion_review_note"].iloc[0]


def test_raw_proportion_tests_require_legacy_exploratory_flag():
    count_df = pd.DataFrame(
        {"T": [0.2, 0.3, 0.8, 0.9]},
        index=["s1", "s2", "s3", "s4"],
    )
    sample_to_cond = pd.Series(["A", "A", "B", "B"], index=count_df.index)

    with pytest.raises(ValueError, match="legacy_exploratory=True"):
        run_statistical_test(
            count_df,
            condition_col="condition",
            test_method="t-test",
            sample_to_cond=sample_to_cond,
        )

    result = run_statistical_test(
        count_df,
        condition_col="condition",
        test_method="t-test",
        sample_to_cond=sample_to_cond,
        legacy_exploratory=True,
    )
    assert not result.empty
    assert set(result["inference_level"]) == {"exploratory_legacy_proportion"}
    assert "compositional" in result["compositional_data_warning"].iloc[0]


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
        legacy_exploratory=True,
        plot_types=["box", "diff", "shift", "paired_shift", "composition_pca", "clr_heatmap"],
    )

    prop_df, stat_df = celltype_proportion_analysis(adata, cfg)

    assert list(prop_df.index) == ["s1", "s2", "s3", "s4"]
    assert "effect_size_cohens_d" in stat_df.columns


def test_proportion_rejects_sample_condition_conflicts():
    """A sample may not be silently assigned the first of several conditions."""
    adata = AnnData(X=np.ones((8, 1)))
    adata.obs["sample"] = ["s1"] * 4 + ["s2"] * 4
    adata.obs["condition"] = ["A", "A", "B", "B"] + ["B"] * 4
    adata.obs["cell_type"] = ["T", "B"] * 4
    cfg = ProportionConfig(
        celltype_col="cell_type",
        sample_col="sample",
        condition_col="condition",
        plot_types=[],
    )

    with pytest.raises(ValueError, match="conflicting sample mappings"):
        celltype_proportion_analysis(adata, cfg)


def test_proportion_paired_design_counts_complete_experimental_units():
    rows = []
    for donor_index, donor in enumerate(["p1", "p2", "p3"]):
        for condition in ["pre", "post"]:
            sample = f"{donor}_{condition}"
            t_cells = 2 + donor_index + (3 if condition == "post" else 0)
            b_cells = 8 - t_cells
            for cell_type in ["T"] * t_cells + ["B"] * b_cells:
                rows.append(
                    {
                        "sample": sample,
                        "condition": condition,
                        "patient": donor,
                        "cell_type": cell_type,
                    }
                )
    adata = AnnData(X=np.ones((len(rows), 1)), obs=pd.DataFrame(rows))
    cfg = ProportionConfig(
        celltype_col="cell_type",
        sample_col="sample",
        condition_col="condition",
        experimental_unit_col="patient",
        pairing_col="patient",
        test_method="clr-paired-t-test",
        auto_configure=False,
        plot_types=[],
    )

    _, stat_df = celltype_proportion_analysis(adata, cfg)

    design = adata.uns["sclucid"]["proportion"]["design"]
    assert design["experimental_unit_col"] == "patient"
    assert design["n_complete_pairs"] == 3
    assert design["replicate_basis"] == "complete_pairs"
    assert set(stat_df["n_complete_pairs"]) == {3}
    assert stat_df["valid_for_publication_inference"].all()


def test_proportion_rejects_pair_column_that_does_not_identify_repeated_units():
    rows = []
    for donor in ["p1", "p2", "p3"]:
        for condition in ["pre", "post"]:
            for cell_type in ["T", "T", "B", "B"]:
                rows.append(
                    {
                        "sample": f"{donor}_{condition}",
                        "condition": condition,
                        "patient": donor,
                        "shared_pair": "all_patients",
                        "cell_type": cell_type,
                    }
                )
    adata = AnnData(X=np.ones((len(rows), 1)), obs=pd.DataFrame(rows))
    config = ProportionConfig(
        celltype_col="cell_type",
        sample_col="sample",
        condition_col="condition",
        experimental_unit_col="patient",
        pairing_col="shared_pair",
        test_method="clr-paired-t-test",
        auto_configure=False,
        plot_types=[],
    )

    with pytest.raises(ValueError, match="must uniquely identify"):
        celltype_proportion_analysis(adata, config)


def test_proportion_does_not_promote_legacy_raw_test_to_formal_inference():
    adata = AnnData(X=np.ones((16, 1)))
    adata.obs["sample"] = np.repeat(["a1", "a2", "b1", "b2"], 4)
    adata.obs["condition"] = adata.obs["sample"].map(
        {"a1": "A", "a2": "A", "b1": "B", "b2": "B"}
    )
    adata.obs["cell_type"] = ["T", "T", "B", "B"] * 4
    cfg = ProportionConfig(
        celltype_col="cell_type",
        sample_col="sample",
        condition_col="condition",
        test_method="t-test",
        auto_configure=False,
        legacy_exploratory=True,
        plot_types=[],
    )

    _, stat_df = celltype_proportion_analysis(adata, cfg)

    assert set(stat_df["inference_level"]) == {"exploratory_legacy_proportion"}
    assert not stat_df["valid_for_publication_inference"].any()


def test_proportion_rejects_condition_batch_confounding():
    adata = AnnData(X=np.ones((16, 1)))
    adata.obs["sample"] = np.repeat(["a1", "a2", "b1", "b2"], 4)
    adata.obs["condition"] = adata.obs["sample"].map(
        {"a1": "A", "a2": "A", "b1": "B", "b2": "B"}
    )
    adata.obs["batch"] = adata.obs["condition"].map({"A": "batch1", "B": "batch2"})
    adata.obs["cell_type"] = ["T", "T", "B", "B"] * 4
    cfg = ProportionConfig(
        celltype_col="cell_type",
        sample_col="sample",
        condition_col="condition",
        batch_col="batch",
        test_method="clr-ols",
        auto_configure=False,
        plot_types=[],
    )

    with pytest.raises(ValueError, match="rank deficient"):
        celltype_proportion_analysis(adata, cfg)


def test_proportion_rejects_missing_values_in_configured_covariate():
    adata = AnnData(X=np.ones((16, 1)))
    adata.obs["sample"] = np.repeat(["a1", "a2", "b1", "b2"], 4)
    adata.obs["condition"] = adata.obs["sample"].map(
        {"a1": "A", "a2": "A", "b1": "B", "b2": "B"}
    )
    adata.obs["batch"] = adata.obs["sample"].map(
        {"a1": "batch1", "a2": None, "b1": "batch1", "b2": "batch2"}
    )
    adata.obs["cell_type"] = ["T", "T", "B", "B"] * 4
    config = ProportionConfig(
        celltype_col="cell_type",
        sample_col="sample",
        condition_col="condition",
        batch_col="batch",
        test_method="clr-ols",
        auto_configure=False,
        plot_types=[],
    )

    with pytest.raises(ValueError, match="complete non-empty metadata in 'batch'"):
        celltype_proportion_analysis(adata, config)


def test_ancom_like_clr_allows_unpaired_two_condition_design():
    adata = AnnData(X=np.ones((24, 1)))
    adata.obs["sample"] = np.repeat(["a1", "a2", "a3", "b1", "b2", "b3"], 4)
    adata.obs["condition"] = adata.obs["sample"].map(
        {"a1": "A", "a2": "A", "a3": "A", "b1": "B", "b2": "B", "b3": "B"}
    )
    adata.obs["cell_type"] = ["T", "T", "B", "B"] * 6
    config = ProportionConfig(
        celltype_col="cell_type",
        sample_col="sample",
        condition_col="condition",
        test_method="ancom-like-clr",
        auto_configure=False,
        plot_types=[],
    )

    _, results = celltype_proportion_analysis(adata, config)

    assert not results.empty
    assert set(results["method"]) == {"ancom-like-clr"}
    assert not results["valid_for_publication_inference"].any()
