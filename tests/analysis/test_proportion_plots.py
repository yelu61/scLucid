"""Tests for proportion plotting helpers."""

import pandas as pd
import pytest

from scLucid.analysis.proportion import (
    plot_celltype_alluvial,
    plot_grouped_celltype_counts,
    plot_grouped_proportion_bar,
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
