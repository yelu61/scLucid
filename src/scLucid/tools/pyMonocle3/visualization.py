"""
Visualization functions for pyMonocle3 (R-free)
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
from typing import Optional, List, Tuple, Union
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import logging

from .core import CellDataSet

log = logging.getLogger(__name__)


def plot_cells(
    cds: CellDataSet,
    reduction_method: str = "UMAP",
    color_cells_by: Optional[str] = None,
    group_cells_by: Optional[str] = None,
    genes: Optional[List[str]] = None,
    show_trajectory_graph: bool = True,
    principal_graph: Optional[dict] = None,
    alpha: Optional[float] = None,
    min_expr: float = 0.0,
    cell_size: float = 0.75,
    cell_stroke_size: float = 0.15,
    normalize: bool = True,
    trajectory_graph_segment_size: float = 0.5,
    label_cell_groups: bool = True,
    label_groups_by_cluster: bool = True,
    label_branch_points: bool = True,
    label_roots: bool = True,
    label_leaves: bool = True,
    graph_label_size: int = 2,
    cell_group_label_size: int = 3,
    show_group_labels: bool = True,
    group_label_font_size: int = 10,
    labels_per_group: int = 1,
    label_groups: bool = False,
    figsize: Tuple[int, int] = (10, 8),
    cmap: str = "viridis",
    **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot cells in reduced dimensional space

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    reduction_method : str
        Which reduction to plot ("UMAP", "tSNE", "PCA")
    color_cells_by : str, optional
        Column to color cells by
    group_cells_by : str, optional
        Column to group cells by
    genes : list, optional
        Genes to plot expression
    show_trajectory_graph : bool
        Show trajectory graph overlay
    principal_graph : dict, optional
        Pre-computed principal graph
    alpha : float, optional
        Point transparency
    min_expr : float
        Minimum expression for plotting
    cell_size : float
        Size of cells
    cell_stroke_size : float
        Width of cell borders
    normalize : bool
        Normalize expression
    trajectory_graph_segment_size : float
        Width of graph edges
    label_cell_groups : bool
        Label cell groups
    label_branch_points : bool
        Label branch points
    label_roots : bool
        Label root cells
    label_leaves : bool
        Label leaf cells
    graph_label_size : int
        Size of graph labels
    cell_group_label_size : int
        Size of cell group labels
    show_group_labels : bool
        Show group labels
    group_label_font_size : int
        Font size for group labels
    labels_per_group : int
        Number of labels per group
    label_groups : bool
        Whether to label groups
    figsize : tuple
        Figure size
    cmap : str
        Colormap for expression
    **kwargs
        Additional arguments for scatter

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """
    if reduction_method not in cds.reducedDims:
        raise ValueError(
            f"{reduction_method} not found in reducedDims. "
            f"Available: {list(cds.reducedDims.keys())}"
        )

    coords = cds.reducedDims[reduction_method]

    fig, ax = plt.subplots(figsize=figsize)

    # Determine coloring
    if genes is not None:
        # Plot gene expression
        if isinstance(genes, str):
            genes = [genes]

        for gene in genes:
            if gene not in cds.gene_metadata.index:
                log.warning(f"Gene {gene} not found")
                continue

            gene_idx = cds.gene_metadata.index.get_loc(gene)
            expr = cds.expression_data[gene_idx, :]
            if sp.issparse(expr):
                expr = expr.toarray().flatten()

            scatter = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=expr,
                s=cell_size * 100,
                cmap=cmap,
                alpha=alpha or 0.6,
                edgecolors='black',
                linewidths=cell_stroke_size,
                **kwargs
            )
            plt.colorbar(scatter, label=f'{gene} expression')
            ax.set_title(f'Gene expression: {gene}')

    elif color_cells_by is not None:
        # Plot by metadata
        if color_cells_by not in cds.cell_metadata.columns:
            raise ValueError(f"Column '{color_cells_by}' not found")

        values = cds.cell_metadata[color_cells_by]

        if values.dtype in ['object', 'category']:
            # Categorical
            categories = values.unique()
            colors = plt.cm.tab20(np.linspace(0, 1, len(categories)))

            for i, cat in enumerate(categories):
                mask = values == cat
                ax.scatter(
                    coords[mask, 0],
                    coords[mask, 1],
                    c=[colors[i]],
                    s=cell_size * 100,
                    label=cat,
                    alpha=alpha or 0.6,
                    edgecolors='black',
                    linewidths=cell_stroke_size,
                    **kwargs
                )

            if show_group_labels:
                ax.legend(loc='best', fontsize=group_label_font_size)
        else:
            # Continuous
            scatter = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=values,
                s=cell_size * 100,
                cmap=cmap,
                alpha=alpha or 0.6,
                edgecolors='black',
                linewidths=cell_stroke_size,
                **kwargs
            )
            plt.colorbar(scatter, label=color_cells_by)

        ax.set_title(f'Colored by {color_cells_by}')

    else:
        # Simple scatter
        ax.scatter(
            coords[:, 0],
            coords[:, 1],
            s=cell_size * 100,
            alpha=alpha or 0.6,
            edgecolors='black',
            linewidths=cell_stroke_size,
            **kwargs
        )

    # Add trajectory graph
    if show_trajectory_graph and cds.principal_graph is not None:
        graph = cds.principal_graph

        for edge in graph['edge_list']:
            i, j, _ = edge
            x = [coords[i, 0], coords[j, 0]]
            y = [coords[i, 1], coords[j, 1]]
            ax.plot(x, y, 'k-', linewidth=trajectory_graph_segment_size, alpha=0.5)

    ax.set_xlabel(f'{reduction_method}_1')
    ax.set_ylabel(f'{reduction_method}_2')

    plt.tight_layout()
    return fig, ax


def plot_genes_by_group(
    cds: CellDataSet,
    markers: pd.DataFrame,
    group_cells_by: str = "cluster",
    ordering_type: str = "cluster_row_col",
    max_width: int = 12,
    norm_method: str = "size_only",
    scale_max: float = 3,
    scale_min: float = -3,
    figsize: Tuple[int, int] = (12, 8),
    cmap: str = "RdBu_r",
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot heatmap of marker genes by group

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    markers : pd.DataFrame
        Marker genes from top_markers()
    group_cells_by : str
        Column to group cells by
    ordering_type : str
        How to order genes and groups
    max_width : int
        Maximum figure width
    norm_method : str
        Normalization method
    scale_max : float
        Maximum value for color scale
    scale_min : float
        Minimum value for color scale
    figsize : tuple
        Figure size
    cmap : str
        Colormap

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """
    import seaborn as sns

    # Get top markers per group
    top_markers = markers.groupby('cell_group').head(10)['gene'].unique()

    # Get expression matrix
    expr = cds.expression_data
    if sp.issparse(expr):
        expr = expr.toarray()

    gene_mask = cds.gene_metadata.index.isin(top_markers)
    expr_subset = expr[gene_mask, :]

    # Aggregate by group
    groups = cds.cell_metadata[group_cells_by].unique()
    agg_expr = []

    for group in groups:
        group_cells = cds.cell_metadata[group_cells_by] == group
        group_expr = expr_subset[:, group_cells].mean(axis=1)
        agg_expr.append(group_expr)

    heatmap_data = pd.DataFrame(
        np.array(agg_expr).T,
        index=cds.gene_metadata.index[gene_mask],
        columns=groups
    )

    # Normalize
    heatmap_data = heatmap_data.apply(
        lambda x: (x - x.mean()) / (x.std() + 1e-10),
        axis=1
    )
    heatmap_data = heatmap_data.clip(scale_min, scale_max)

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(heatmap_data, cmap=cmap, center=0, ax=ax, xticklabels=True, yticklabels=True)
    ax.set_xlabel('Cell Group')
    ax.set_ylabel('Gene')
    ax.set_title('Marker Gene Expression by Group')

    plt.tight_layout()
    return fig, ax


def plot_pseudotime_heatmap(
    cds: CellDataSet,
    genes: List[str],
    num_bins: int = 100,
    max_col: float = 2.5,
    min_col: float = -2.5,
    trend_formula: str = "~ sm.ns(pseudotime, df=3)",
    cell_attributes: Optional[List[str]] = None,
    normalize: bool = True,
    return_heatmap: bool = False,
    figsize: Tuple[int, int] = (10, 10),
    cmap: str = "viridis",
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot heatmap of gene expression over pseudotime

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet with pseudotime
    genes : list
        Genes to plot
    num_bins : int
        Number of pseudotime bins
    max_col : float
        Maximum value for color scale
    min_col : float
        Minimum value for color scale
    trend_formula : str
        Formula for fitting trends
    cell_attributes : list, optional
        Additional cell attributes to show
    normalize : bool
        Normalize expression
    return_heatmap : bool
        Return heatmap data
    figsize : tuple
        Figure size
    cmap : str
        Colormap

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """
    import seaborn as sns

    if 'pseudotime' not in cds.cell_metadata.columns:
        raise ValueError("No pseudotime found. Run order_cells first.")

    pseudotime = cds.cell_metadata['pseudotime'].values

    # Get expression
    expr = cds.expression_data
    if sp.issparse(expr):
        expr = expr.toarray()

    gene_mask = cds.gene_metadata.index.isin(genes)
    expr_subset = expr[gene_mask, :]

    # Bin cells by pseudotime
    bins = np.linspace(pseudotime.min(), pseudotime.max(), num_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    binned_expr = []
    for i in range(num_bins):
        mask = (pseudotime >= bins[i]) & (pseudotime < bins[i + 1])
        if mask.sum() > 0:
            binned_expr.append(expr_subset[:, mask].mean(axis=1))
        else:
            binned_expr.append(np.zeros(expr_subset.shape[0]))

    heatmap_data = pd.DataFrame(
        np.array(binned_expr).T,
        index=cds.gene_metadata.index[gene_mask],
        columns=bin_centers
    )

    # Normalize by gene
    if normalize:
        heatmap_data = heatmap_data.apply(
            lambda x: (x - x.mean()) / (x.std() + 1e-10),
            axis=1
        )
        heatmap_data = heatmap_data.clip(min_col, max_col)

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        heatmap_data,
        cmap=cmap,
        xticklabels=10,
        yticklabels=True,
        ax=ax
    )
    ax.set_xlabel('Pseudotime')
    ax.set_ylabel('Gene')
    ax.set_title('Gene Expression along Pseudotime')

    plt.tight_layout()

    if return_heatmap:
        return fig, ax, heatmap_data
    return fig, ax


def plot_trajectory(
    cds: CellDataSet,
    reduction_method: str = "UMAP",
    color_by: str = "pseudotime",
    show_graph: bool = True,
    show_cell_names: bool = False,
    cell_size: float = 0.5,
    edge_width: float = 1.0,
    figsize: Tuple[int, int] = (10, 8),
    cmap: str = "viridis",
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot trajectory with cells colored by pseudotime or other variable

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    reduction_method : str
        Reduction to use
    color_by : str
        Column to color by (default: pseudotime)
    show_graph : bool
        Show principal graph
    show_cell_names : bool
        Show cell names
    cell_size : float
        Size of cells
    edge_width : float
        Width of graph edges
    figsize : tuple
        Figure size
    cmap : str
        Colormap

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """
    if reduction_method not in cds.reducedDims:
        raise ValueError(f"{reduction_method} not found in reducedDims")

    coords = cds.reducedDims[reduction_method]

    fig, ax = plt.subplots(figsize=figsize)

    # Get coloring
    if color_by in cds.cell_metadata.columns:
        values = cds.cell_metadata[color_by]
    elif color_by == "pseudotime":
        raise ValueError(f"'{color_by}' not found. Run order_cells first.")
    else:
        values = np.arange(len(coords))

    # Plot cells
    scatter = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=values,
        s=cell_size * 100,
        cmap=cmap,
        alpha=0.6,
        edgecolors='black',
        linewidths=0.1,
    )
    plt.colorbar(scatter, label=color_by)

    # Plot graph
    if show_graph and cds.principal_graph is not None:
        for edge in cds.principal_graph['edge_list']:
            i, j, _ = edge
            ax.plot(
                [coords[i, 0], coords[j, 0]],
                [coords[i, 1], coords[j, 1]],
                'k-',
                linewidth=edge_width,
                alpha=0.5
            )

    ax.set_xlabel(f'{reduction_method}_1')
    ax.set_ylabel(f'{reduction_method}_2')
    ax.set_title(f'Trajectory colored by {color_by}')

    plt.tight_layout()
    return fig, ax
