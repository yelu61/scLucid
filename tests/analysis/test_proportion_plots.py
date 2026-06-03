"""Tests for proportion plotting helpers."""

import pandas as pd
import pytest

from scLucid.analysis.proportion import (
    plot_batch_effect,
    plot_box_summary,
    plot_celltype_alluvial,
    plot_celltype_variability,
    plot_composition,
    plot_composition_pca,
    plot_composition_transform_heatmap,
    plot_diff_stats,
    plot_grouped_celltype_counts,
    plot_grouped_proportion_bar,
    plot_individual_boxplots,
    plot_paired_proportion_shifts,
    plot_proportion_heatmap,
    plot_proportion_shifts,
    plot_proportion_with_ci,
    transform_composition,
)


@pytest.mark.unit
def test_plot_grouped_celltype_counts_returns_figure():
    count_df = pd.DataFrame(
        {
            "group": ["A", "A", "B", "B"],
            "cell_type": ["T", "NK", "T", "NK"],
            "count": [10, 5, 7, 8],
        }
    )
    fig = plot_grouped_celltype_counts(count_df, annotate=True)
    assert fig is not None
    assert len(fig.axes) == 1


@pytest.mark.unit
def test_plot_grouped_proportion_bar_returns_figure():
    group_props = pd.DataFrame(
        {"T": [0.6, 0.4], "NK": [0.4, 0.6]},
        index=["A", "B"],
    )
    fig = plot_grouped_proportion_bar(group_props, group_order=["B", "A"])
    assert fig is not None
    assert len(fig.axes) == 1


@pytest.mark.unit
def test_plot_celltype_alluvial_returns_figure():
    group_props = pd.DataFrame(
        {"T": [0.7, 0.5, 0.4], "NK": [0.3, 0.5, 0.6]},
        index=["A", "B", "C"],
    )
    fig = plot_celltype_alluvial(group_props)
    assert fig is not None
    assert len(fig.axes) == 1


def _sample_prop_df():
    return pd.DataFrame(
        {
            "T": [0.7, 0.6, 0.3, 0.2],
            "NK": [0.2, 0.3, 0.5, 0.6],
            "B": [0.1, 0.1, 0.2, 0.2],
        },
        index=["s1", "s2", "s3", "s4"],
    )


def _sample_condition():
    return pd.Series(["A", "A", "B", "B"], index=["s1", "s2", "s3", "s4"], name="condition")


def _sample_stats_df():
    return pd.DataFrame(
        {
            "cell_type": ["T", "NK", "B"],
            "mean_diff": [0.4, -0.3, -0.1],
            "pval": [0.01, 0.2, 0.04],
            "padj": [0.03, 0.2, 0.06],
            "effect_size_cohens_d": [2.0, -1.2, -0.5],
        }
    )


@pytest.mark.unit
def test_plot_box_summary_uses_condition_groups():
    fig = plot_box_summary(_sample_prop_df(), _sample_condition())

    assert fig is not None
    assert len(fig.axes) == 1
    assert fig.axes[0].get_xlabel() == "Cell Type"


@pytest.mark.unit
def test_new_proportion_plot_helpers_return_figures():
    prop_df = _sample_prop_df()
    condition = _sample_condition()
    stat_df = _sample_stats_df()
    shift_df = prop_df.copy()
    shift_df["condition"] = condition

    figs = [
        plot_composition(prop_df, condition),
        plot_diff_stats(prop_df, stat_df, condition),
        plot_individual_boxplots(prop_df, condition, stat_df),
        plot_proportion_shifts(shift_df, "condition", "A", "B"),
        plot_proportion_with_ci(prop_df, condition),
        plot_celltype_variability(prop_df),
    ]

    assert all(fig is not None for fig in figs)


@pytest.mark.unit
def test_compositional_transform_and_heatmap_return_expected_shape():
    prop_df = _sample_prop_df()

    clr = transform_composition(prop_df, method="clr")
    fig = plot_composition_transform_heatmap(prop_df, transform="clr")

    assert clr.shape == prop_df.shape
    assert clr.mean(axis=1).abs().max() < 1e-10
    assert fig is not None


@pytest.mark.unit
def test_composition_pca_and_paired_shift_return_figures():
    prop_df = _sample_prop_df()
    condition = _sample_condition()
    pair = pd.Series(["p1", "p2", "p1", "p2"], index=prop_df.index, name="patient")

    fig_pca = plot_composition_pca(prop_df, condition, top_loadings=2)
    fig_paired = plot_paired_proportion_shifts(
        prop_df,
        condition=condition,
        pair=pair,
        condition1="A",
        condition2="B",
    )

    assert fig_pca is not None
    assert fig_paired is not None


@pytest.mark.unit
def test_heatmap_clustering_and_batch_method_validation():
    prop_df = _sample_prop_df()
    batch = pd.Series(["x", "x", "y", "y"], index=prop_df.index, name="batch")

    fig_clustered = plot_proportion_heatmap(
        prop_df,
        cluster_samples=True,
        cluster_celltypes=True,
    )

    assert fig_clustered is not None
    with pytest.raises(ValueError):
        plot_batch_effect(prop_df, batch=batch, method="umap")
