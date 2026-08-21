"""Tests for notebook-derived stable plotting APIs."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import scLucid.analysis as analysis
from scLucid.analysis.differential_expression import (
    plot_categorized_gene_heatmap,
    plot_grouped_marker_dotplot,
)
from scLucid.analysis.proportion import (
    plot_composition_shift_bubble,
    plot_composition_shift_effect,
    summarize_composition_shift,
)


def _small_adata() -> AnnData:
    adata = AnnData(
        np.array(
            [
                [5.0, 1.0, 0.0, 2.0],
                [4.0, 2.0, 0.0, 1.0],
                [1.0, 5.0, 3.0, 0.0],
                [0.0, 4.0, 4.0, 1.0],
                [3.0, 0.0, 1.0, 5.0],
                [1.0, 3.0, 5.0, 0.0],
            ]
        )
    )
    adata.obs_names = [f"cell_{idx}" for idx in range(adata.n_obs)]
    adata.var_names = ["GeneA", "GeneB", "GeneC", "GeneD"]
    adata.obs["cell_type"] = pd.Categorical(["T", "T", "B", "B", "T", "B"])
    adata.obs["condition"] = pd.Categorical(["ctrl", "case", "ctrl", "case", "ctrl", "case"])
    adata.raw = adata.copy()
    return adata


@pytest.mark.unit
def test_composition_shift_summary_and_plots_are_public_api():
    prop_df = pd.DataFrame(
        {
            "sample": ["s1", "s2", "s3", "s4"],
            "condition": ["ctrl", "ctrl", "case", "case"],
            "T": [0.2, 0.3, 0.5, 0.6],
            "B": [0.8, 0.7, 0.5, 0.4],
        }
    )
    stat_df = pd.DataFrame({"cell_type": ["T", "B"], "padj": [0.04, 0.20]})

    summary = summarize_composition_shift(
        prop_df,
        condition_col="condition",
        condition1="ctrl",
        condition2="case",
        stat_df=stat_df,
    )

    assert analysis.summarize_composition_shift is summarize_composition_shift
    assert {"cell_type", "delta", "abs_delta", "log2_fc", "padj"}.issubset(summary.columns)
    assert summary.loc[summary["cell_type"] == "T", "delta"].iloc[0] == pytest.approx(0.3)

    bubble = plot_composition_shift_bubble(summary, q_col="padj")
    effect = plot_composition_shift_effect(summary)

    assert bubble is not None
    assert effect is not None
    plt.close(bubble)
    plt.close(effect)


@pytest.mark.unit
def test_grouped_marker_dotplot_filters_missing_genes_and_returns_scanpy_plot():
    adata = _small_adata()

    dotplot = plot_grouped_marker_dotplot(
        adata,
        {"immune": ["GeneA", "GeneB", "MissingGene"]},
        celltype_col="cell_type",
        condition_col="condition",
        celltype_order=["T", "B"],
        condition_order=["ctrl", "case"],
    )

    assert analysis.plot_grouped_marker_dotplot is plot_grouped_marker_dotplot
    assert dotplot is not None
    assert "_sclucid_celltype_condition" not in adata.obs
    plt.close("all")


@pytest.mark.unit
def test_categorized_gene_heatmap_returns_mean_expression_table():
    adata = _small_adata()

    ax, mean_df = plot_categorized_gene_heatmap(
        adata,
        {"lineage": ["GeneA", "GeneC"], "state": ["GeneD", "MissingGene"]},
        groupby=["cell_type", "condition"],
        title="Functional marker programs",
    )

    assert analysis.plot_categorized_gene_heatmap is plot_categorized_gene_heatmap
    assert list(mean_df.index) == ["GeneA", "GeneC", "GeneD"]
    assert {"T | ctrl", "B | case"}.issubset(mean_df.columns)
    assert ax.get_title() == "Functional marker programs"
    plt.close(ax.figure)
