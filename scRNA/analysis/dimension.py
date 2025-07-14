"""
Dimensionality reduction and visualization functions for single-cell RNA-seq data.

This module provides functions for visualizing marker gene expression and
cell type compositions across clusters.
"""

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from adjustText import adjust_text


def plot_embedding(
    adata: sc.AnnData,
    color_by: str,
    basis: str = "umap",
    title: Optional[str] = None,
    show_labels: bool = True,
    **kwargs,
) -> None:
    """

    Custom wrapper for `sc.pl.embedding` with improved labeling.

    Args:
        adata: AnnData object.
        color_by: Key in `adata.obs` to color by.
        basis: Basis for the embedding (e.g., 'umap').
        title: Plot title.
        show_labels: Whether to show cluster/group labels on the plot.
        **kwargs: Additional arguments passed to `sc.pl.embedding`.
    """
    fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 7)))

    sc.pl.embedding(
        adata,
        basis=basis,
        color=color_by,
        ax=ax,
        show=False,
        legend_loc="on data" if show_labels else "right margin",
        **kwargs,
    )

    if show_labels:
        ax.legend_.remove()  # Remove default legend to replace with adjust_text

        texts = []
        for label in adata.obs[color_by].cat.categories:
            x, y = np.median(
                adata[adata.obs[color_by] == label].obsm[f"X_{basis}"], axis=0
            )
            texts.append(ax.text(x, y, label, fontsize=10, fontweight="bold"))

        adjust_text(
            texts, ax=ax, arrowprops=dict(arrowstyle="-", color="black", lw=0.5)
        )

    ax.set_title(title if title else f"{basis.upper()} colored by {color_by}")
    plt.tight_layout()
    plt.show()


def plot_marker_heatmap(
    adata: sc.AnnData,
    markers_df: pd.DataFrame,
    groupby: str,
    n_genes: int = 5,
    **kwargs,
):
    """
    Plots a heatmap of top marker genes from a markers DataFrame.

    Args:
        adata: AnnData object.
        markers_df: DataFrame from `find_markers` or `filter_markers`.
        groupby: Key in `adata.obs` used for grouping.
        n_genes: Number of top genes to show per group.
        **kwargs: Additional arguments passed to `sc.pl.heatmap`.
    """
    top_markers = markers_df.groupby("group").head(n_genes)
    marker_dict = top_markers.groupby("group")["names"].apply(list).to_dict()

    sc.pl.heatmap(
        adata, marker_dict, groupby=groupby, dendrogram=True, swap_axes=True, **kwargs
    )


def plot_composition(
    adata: sc.AnnData, group_key: str, stack_key: str, normalize: bool = True, **kwargs
):
    """
    Generates a stacked bar plot showing the composition of groups.

    Args:
        adata: AnnData object.
        group_key: Key in `adata.obs` for the x-axis groups (e.g., 'leiden').
        stack_key: Key in `adata.obs` for the stacked bars (e.g., 'cell_type_annotated').
        normalize: If True, show proportions; otherwise, show raw counts.
    """
    composition = pd.crosstab(
        adata.obs[group_key],
        adata.obs[stack_key],
        normalize="index" if normalize else "all",
    )

    composition.plot.bar(stacked=True, figsize=kwargs.pop("figsize", (12, 7)), **kwargs)
    plt.ylabel("Proportion" if normalize else "Count")
    plt.legend(title=stack_key, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.show()
