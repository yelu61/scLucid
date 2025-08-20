"""
Dimensionality reduction and visualization functions for single-cell RNA-seq data.

This module provides comprehensive visualization tools for single-cell RNA-seq data,
including functions for exploring dimensionality reduction embeddings, visualizing marker
gene expression, examining cell type compositions across clusters, and analyzing
pseudotime trajectories. The functions are designed to create publication-quality
figures with customizable parameters.

Key features:
- Customizable embedding plots (UMAP, t-SNE, PCA) with automatic labeling
- Marker gene expression visualization across different cell groups
- Cell type composition analysis across clusters or conditions
- Functional enrichment visualization
- Trajectory and pseudotime visualization
- Multi-modal data visualization
"""

import logging
import os
from typing import Dict, List, Literal, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import seaborn as sns
from adjustText import adjust_text
from matplotlib.patches import Patch

# Configure logging
log = logging.getLogger(__name__)

__all__ = [
    "plot_embedding",
    "plot_marker_expression",
    "plot_marker_heatmap",
    "plot_composition",
    "plot_enrichment",
    "plot_pseudotime",
    "plot_spatial",
    "plot_volcano",
    "plot_feature_correlation",
    "plot_multi_modality",
    "plot_ridge",
    "plot_coexpression"
]


def plot_embedding(
    adata: sc.AnnData,
    color_by: Union[str, List[str]],
    basis: str = "umap",
    title: Optional[str] = None,
    show_labels: bool = True,
    palette: Optional[Union[str, Dict[str, str]]] = None,
    size: float = 10,
    alpha: float = 0.8,
    ncols: int = 1,
    figsize: Optional[Tuple[float, float]] = None,
    legend_loc: str = "on data",
    label_size: int = 10,
    save: Optional[str] = None,
    dpi: int = 300,
    ax: Optional[plt.Axes] = None,
    show: bool = True,
    adjust_label_positions: bool = True,
    **kwargs,
) -> Union[plt.Figure, List[plt.Figure]]:
    """
    Enhanced wrapper for scanpy's embedding plot with improved labeling and customization.

    This function creates publication-quality embedding plots (UMAP, t-SNE, PCA) with
    customizable colors, labels, and styling. It supports automatic label positioning
    to avoid overlaps and can generate multiple plots for different colorings.

    Args:
        adata: AnnData object containing the embedding coordinates
        color_by: Key(s) in adata.obs for coloring cells (categorical or continuous)
        basis: Dimensionality reduction to use ('umap', 'tsne', 'pca', etc.)
        title: Plot title (will be generated automatically if None)
        show_labels: Whether to show cluster/group labels on the plot
        palette: Color palette for categorical variables (name or dict mapping categories to colors)
        size: Point size for the scatter plot
        alpha: Transparency level for the scatter plot
        ncols: Number of columns if plotting multiple features
        figsize: Figure size (width, height) in inches
        legend_loc: Location of the legend ('on data', 'right margin', 'best', etc.)
        label_size: Font size for cluster labels
        save: Path to save the figure (if not None)
        dpi: Resolution for saved figure
        ax: Matplotlib axes to plot on (creates new figure if None)
        show: Whether to display the plot
        adjust_label_positions: Whether to use adjustText to position labels optimally
        **kwargs: Additional arguments passed to sc.pl.embedding

    Returns:
        Figure object(s) for further customization

    Examples:
        >>> # Basic UMAP colored by cluster
        >>> plot_embedding(adata, 'leiden')
        >>>
        >>> # t-SNE with custom colors and no labels
        >>> plot_embedding(
        ...     adata,
        ...     'cell_type',
        ...     basis='tsne',
        ...     show_labels=False,
        ...     palette={'T cells': 'green', 'B cells': 'blue'}
        ... )
        >>>
        >>> # Multiple plots showing different features
        >>> plot_embedding(
        ...     adata,
        ...     ['leiden', 'cell_type', 'CD3D'],
        ...     ncols=2,
        ...     save='embeddings.png'
        ... )
    """
    # Handle multiple color_by as a list
    if isinstance(color_by, list):
        if len(color_by) == 1:
            return plot_embedding(
                adata,
                color_by[0],
                basis,
                title,
                show_labels,
                palette,
                size,
                alpha,
                ncols,
                figsize,
                legend_loc,
                label_size,
                save,
                dpi,
                ax,
                show,
                adjust_label_positions,
                **kwargs,
            )

        # Create a grid of plots
        n_plots = len(color_by)
        n_rows = int(np.ceil(n_plots / ncols))

        # Calculate figsize if not provided
        if figsize is None:
            figsize = (5 * ncols, 5 * n_rows)

        # Create figure and axes
        fig, axes = plt.subplots(n_rows, ncols, figsize=figsize, squeeze=False)
        axes = axes.flatten()

        # Create each plot
        for i, color in enumerate(color_by):
            if i < n_plots:
                plot_embedding(
                    adata,
                    color,
                    basis,
                    title=None,
                    show_labels=show_labels,
                    palette=palette,
                    size=size,
                    alpha=alpha,
                    ax=axes[i],
                    legend_loc=legend_loc,
                    label_size=label_size,
                    show=False,
                    adjust_label_positions=adjust_label_positions,
                    **kwargs,
                )
                axes[i].set_title(color)
            else:
                axes[i].axis("off")

        plt.tight_layout()

        if save:
            plt.savefig(save, dpi=dpi, bbox_inches="tight")

        if show:
            plt.show()

        return fig

    # Check if the embedding exists
    embed_key = f"X_{basis}"
    if embed_key not in adata.obsm:
        log.error(f"Embedding '{embed_key}' not found in adata.obsm")
        raise ValueError(f"Embedding '{embed_key}' not found in adata.obsm")

    # Check if color_by exists
    if color_by not in adata.obs and color_by not in adata.var_names:
        log.error(f"'{color_by}' not found in adata.obs or adata.var_names")
        raise ValueError(f"'{color_by}' not found in adata.obs or adata.var_names")

    # Create figure if ax not provided
    if ax is None:
        if figsize is None:
            figsize = (8, 7)
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # Determine if color_by is categorical
    categorical = False
    if color_by in adata.obs:
        if pd.api.types.is_categorical_dtype(adata.obs[color_by]):
            categorical = True
            categories = adata.obs[color_by].cat.categories
        elif pd.api.types.is_bool_dtype(adata.obs[color_by]):
            # Convert boolean to categorical for better visualization
            categorical = True
            adata.obs[f"{color_by}_cat"] = (
                adata.obs[color_by]
                .map({True: "True", False: "False"})
                .astype("category")
            )
            color_by = f"{color_by}_cat"
            categories = adata.obs[color_by].cat.categories

    # Plot the embedding
    sc.pl.embedding(
        adata,
        basis=basis,
        color=color_by,
        ax=ax,
        show=False,
        legend_loc=legend_loc
        if not (categorical and show_labels and adjust_label_positions)
        else "none",
        size=size,
        alpha=alpha,
        palette=palette,
        **kwargs,
    )

    # Add labels for categorical data if requested
    if categorical and show_labels:
        if adjust_label_positions:
            # Remove the default legend if we're going to use adjustText
            if hasattr(ax, "legend_") and ax.legend_ is not None:
                ax.legend_.remove()

            # Add labels at the center of each category
            texts = []
            for label in categories:
                mask = adata.obs[color_by] == label
                if np.sum(mask) > 0:  # Only add label if category has cells
                    x, y = np.median(adata[mask].obsm[embed_key], axis=0)
                    # Get color if palette is provided
                    if isinstance(palette, dict) and label in palette:
                        color = palette[label]
                    else:
                        # Try to get color from the scatter plot
                        sc_artists = [
                            c
                            for c in ax.collections
                            if isinstance(c, plt.PathCollection)
                        ]
                        if sc_artists and hasattr(sc_artists[0], "get_facecolor"):
                            # Find index of this category in the plot
                            try:
                                idx = np.where(categories == label)[0][0]
                                if idx < len(sc_artists):
                                    colors = sc_artists[idx].get_facecolor()
                                    if len(colors) > 0:
                                        color = colors[0]
                                    else:
                                        color = "black"
                                else:
                                    color = "black"
                            except:
                                color = "black"
                        else:
                            color = "black"

                    texts.append(
                        ax.text(
                            x,
                            y,
                            label,
                            fontsize=label_size,
                            fontweight="bold",
                            ha="center",
                            va="center",
                            color="black",
                            bbox=dict(
                                facecolor="white",
                                alpha=0.7,
                                edgecolor="none",
                                boxstyle="round,pad=0.2",
                            ),
                        )
                    )

            # Adjust label positions to avoid overlaps
            if texts:
                try:
                    adjust_text(
                        texts,
                        ax=ax,
                        arrowprops=dict(arrowstyle="-", color="black", lw=0.5),
                        expand_points=(1.5, 1.5),
                    )
                except Exception as e:
                    log.warning(f"Error adjusting text positions: {str(e)}")
        else:
            # Keep the default legend
            pass

    # Set title
    if title:
        ax.set_title(title)
    else:
        ax.set_title(f"{basis.upper()} colored by {color_by}")

    # Finalize plot
    plt.tight_layout()

    # Save if requested
    if save:
        plt.savefig(save, dpi=dpi, bbox_inches="tight")

    # Show if requested
    if show:
        plt.show()

    return fig


def plot_marker_expression(
    adata: sc.AnnData,
    markers: Union[str, List[str]],
    basis: str = "umap",
    ncols: int = 4,
    figsize: Optional[Tuple[float, float]] = None,
    cmap: str = "viridis",
    size: float = 10,
    alpha: float = 0.8,
    use_raw: bool = False,
    layer: Optional[str] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """
    Visualize expression of marker genes on an embedding plot.

    This function creates a grid of embedding plots, each showing the expression
    of a different marker gene. It's useful for visualizing the spatial distribution
    of gene expression across cell populations.

    Args:
        adata: AnnData object containing the data
        markers: Gene name(s) to visualize
        basis: Dimensionality reduction to use ('umap', 'tsne', 'pca', etc.)
        ncols: Number of columns in the plot grid
        figsize: Figure size (width, height) in inches
        cmap: Colormap for gene expression
        size: Point size for the scatter plot
        alpha: Transparency level for the scatter plot
        use_raw: Whether to use adata.raw for expression values
        layer: Layer in adata.layers to use for expression values
        vmin: Minimum value for color scale (default: 2nd percentile)
        vmax: Maximum value for color scale (default: 98th percentile)
        save: Path to save the figure (if not None)
        dpi: Resolution for saved figure
        show: Whether to display the plot
        **kwargs: Additional arguments passed to sc.pl.embedding

    Returns:
        Figure object for further customization

    Examples:
        >>> # Visualize a single marker gene
        >>> plot_marker_expression(adata, 'CD3D')
        >>>
        >>> # Visualize multiple T cell markers
        >>> plot_marker_expression(
        ...     adata,
        ...     ['CD3D', 'CD4', 'CD8A', 'FOXP3'],
        ...     cmap='Reds',
        ...     layer='log1p_norm'
        ... )
        >>>
        >>> # Use custom color scaling
        >>> plot_marker_expression(
        ...     adata,
        ...     ['IL6', 'TNF', 'IFNG'],
        ...     vmin=0,
        ...     vmax=3,
        ...     save='cytokine_expression.png'
        ... )
    """
    # Convert single marker to list
    if isinstance(markers, str):
        markers = [markers]

    # Check if markers exist
    missing_markers = [m for m in markers if m not in adata.var_names]
    if missing_markers:
        if len(missing_markers) == len(markers):
            log.error("None of the specified markers found in adata.var_names")
            raise ValueError("None of the specified markers found in adata.var_names")
        else:
            log.warning(
                f"Some markers not found in adata.var_names: {', '.join(missing_markers)}"
            )
            markers = [m for m in markers if m not in missing_markers]

    # Determine grid layout
    n_markers = len(markers)
    n_rows = int(np.ceil(n_markers / ncols))

    # Calculate figsize if not provided
    if figsize is None:
        figsize = (4 * min(ncols, n_markers), 4 * n_rows)

    # Create figure and axes
    fig, axes = plt.subplots(
        n_rows, min(ncols, n_markers), figsize=figsize, squeeze=False
    )
    axes = axes.flatten()

    # Plot each marker
    for i, marker in enumerate(markers):
        if i < len(axes):
            try:
                sc.pl.embedding(
                    adata,
                    basis=basis,
                    color=marker,
                    ax=axes[i],
                    show=False,
                    cmap=cmap,
                    size=size,
                    alpha=alpha,
                    use_raw=use_raw,
                    layer=layer,
                    vmin=vmin,
                    vmax=vmax,
                    **kwargs,
                )

                # Remove colorbar title which is redundant
                if axes[i].get_title() == marker:
                    axes[i].set_title(marker, fontweight="bold")

            except Exception as e:
                log.warning(f"Error plotting marker {marker}: {str(e)}")
                axes[i].text(
                    0.5,
                    0.5,
                    f"Error plotting {marker}",
                    ha="center",
                    va="center",
                    transform=axes[i].transAxes,
                )
                axes[i].set_title(marker)

    # Hide empty axes
    for i in range(n_markers, len(axes)):
        axes[i].axis("off")

    plt.tight_layout()

    # Save if requested
    if save:
        plt.savefig(save, dpi=dpi, bbox_inches="tight")

    # Show if requested
    if show:
        plt.show()

    return fig


def plot_marker_heatmap(
    adata: sc.AnnData,
    markers_df: Optional[pd.DataFrame] = None,
    groupby: str = None,
    markers: Optional[Union[List[str], Dict[str, List[str]]]] = None,
    n_genes: int = 5,
    layer: Optional[str] = None,
    use_raw: bool = False,
    standard_scale: Optional[Literal["var", "group"]] = "var",
    cmap: str = "viridis",
    colorbar_title: str = "Expression",
    dendrogram: bool = True,
    swap_axes: bool = True,
    figsize: Optional[Tuple[float, float]] = None,
    show_gene_labels: bool = True,
    gene_label_size: int = 8,
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """
    Create a heatmap visualization of marker gene expression across cell groups.

    This function generates a heatmap showing the expression of marker genes across
    different cell groups. It supports multiple input formats for specifying markers
    and offers customization options for the visualization.

    Args:
        adata: AnnData object containing the data
        markers_df: DataFrame from find_markers() or filter_markers() (optional)
        groupby: Key in adata.obs for grouping cells
        markers: Dictionary mapping groups to marker lists, or a list of marker genes
        n_genes: Number of top genes to show per group when using markers_df
        layer: Layer in adata.layers to use for expression values
        use_raw: Whether to use adata.raw for expression values
        standard_scale: Whether to scale data by 'var' (genes) or 'group'
        cmap: Colormap for the heatmap
        colorbar_title: Title for the colorbar
        dendrogram: Whether to show dendrograms
        swap_axes: Whether to swap the x and y axes
        figsize: Figure size (width, height) in inches
        show_gene_labels: Whether to show gene labels
        gene_label_size: Font size for gene labels
        save: Path to save the figure (if not None)
        dpi: Resolution for saved figure
        show: Whether to display the plot
        **kwargs: Additional arguments passed to sc.pl.heatmap

    Returns:
        Figure object for further customization

    Examples:
        >>> # Using marker dataframe from find_markers()
        >>> markers_df = find_markers(adata, groupby='leiden')
        >>> plot_marker_heatmap(adata, markers_df=markers_df, groupby='leiden')
        >>>
        >>> # Using a dictionary of markers
        >>> markers_dict = {
        ...     'T cells': ['CD3D', 'CD4', 'CD8A'],
        ...     'B cells': ['MS4A1', 'CD79A'],
        ...     'NK cells': ['NCAM1', 'NKG7']
        ... }
        >>> plot_marker_heatmap(adata, markers=markers_dict, groupby='cell_type')
        >>>
        >>> # Using a simple list of markers
        >>> plot_marker_heatmap(
        ...     adata,
        ...     markers=['CD3D', 'CD4', 'CD8A', 'MS4A1', 'CD79A', 'NCAM1'],
        ...     groupby='cell_type',
        ...     standard_scale=None,
        ...     cmap='Reds'
        ... )
    """
    log.info("Creating marker gene heatmap visualization")

    # Validate and process input markers
    if markers_df is not None:
        # Case 1: Using a dataframe from find_markers()
        if "group" not in markers_df.columns or "names" not in markers_df.columns:
            log.error("markers_df must contain 'group' and 'names' columns")
            raise ValueError("markers_df must contain 'group' and 'names' columns")

        # Extract top markers per group
        log.info(f"Using top {n_genes} markers per group from markers_df")
        marker_dict = {}
        for group in markers_df["group"].unique():
            group_markers = markers_df[markers_df["group"] == group]
            if "logfoldchanges" in group_markers.columns:
                group_markers = group_markers.sort_values(
                    "logfoldchanges", ascending=False
                )
            marker_dict[group] = group_markers["names"].head(n_genes).tolist()

    elif markers is not None:
        if isinstance(markers, dict):
            # Case 2: Dictionary mapping groups to marker lists
            marker_dict = markers
            log.info(f"Using provided marker dictionary with {len(marker_dict)} groups")
        elif isinstance(markers, (list, tuple)):
            # Case 3: Simple list of markers
            if groupby is None:
                log.error("groupby must be specified when markers is a list")
                raise ValueError("groupby must be specified when markers is a list")

            marker_dict = {"Markers": markers}
            log.info(f"Using provided list of {len(markers)} markers")
        else:
            log.error("markers must be a dictionary or list")
            raise TypeError("markers must be a dictionary or list")
    else:
        log.error("Either markers_df or markers must be provided")
        raise ValueError("Either markers_df or markers must be provided")

    # Validate groupby
    if groupby is not None and groupby not in adata.obs.columns:
        log.error(f"groupby '{groupby}' not found in adata.obs")
        raise ValueError(f"groupby '{groupby}' not found in adata.obs")

    # Validate that markers exist in the dataset
    all_markers = [gene for genes in marker_dict.values() for gene in genes]
    missing_markers = [gene for gene in all_markers if gene not in adata.var_names]

    if missing_markers:
        log.warning(f"{len(missing_markers)} markers not found in dataset")
        if len(missing_markers) < 10:
            log.warning(f"Missing markers: {', '.join(missing_markers)}")
        else:
            log.warning(
                f"First 10 missing markers: {', '.join(missing_markers[:10])}..."
            )

        # Remove missing markers from each group
        for group in marker_dict:
            marker_dict[group] = [
                gene for gene in marker_dict[group] if gene in adata.var_names
            ]

    # Calculate figsize if not provided
    if figsize is None:
        n_groups = len(adata.obs[groupby].cat.categories) if groupby else 1
        n_genes = len(all_markers)

        if swap_axes:
            figsize = (max(6, min(12, n_groups * 0.6)), max(5, min(16, n_genes * 0.25)))
        else:
            figsize = (max(6, min(16, n_genes * 0.4)), max(5, min(12, n_groups * 0.4)))

    # Create the heatmap
    try:
        # Store the current figure size to restore it later
        original_figsize = plt.rcParams["figure.figsize"]
        plt.rcParams["figure.figsize"] = figsize

        # Generate the heatmap
        ax = sc.pl.heatmap(
            adata,
            marker_dict,
            groupby=groupby,
            layer=layer,
            use_raw=use_raw,
            standard_scale=standard_scale,
            cmap=cmap,
            dendrogram=dendrogram,
            swap_axes=swap_axes,
            show_gene_labels=show_gene_labels,
            gene_symbols=None,
            show=False,
            **kwargs,
        )

        # Get the figure from the axes
        fig = ax["heatmap_ax"].figure

        # Customize gene labels if needed
        if show_gene_labels and gene_label_size != 8:  # 8 is scanpy's default
            if swap_axes:
                for text in ax["heatmap_ax"].get_xticklabels():
                    text.set_fontsize(gene_label_size)
            else:
                for text in ax["heatmap_ax"].get_yticklabels():
                    text.set_fontsize(gene_label_size)

        # Customize colorbar if needed
        if "colorbar_ax" in ax and ax["colorbar_ax"] is not None:
            ax["colorbar_ax"].set_title(colorbar_title)

        # Restore original figure size
        plt.rcParams["figure.figsize"] = original_figsize

        # Save if requested
        if save:
            plt.savefig(save, dpi=dpi, bbox_inches="tight")

        # Show if requested
        if show:
            plt.show()

        return fig

    except Exception as e:
        log.error(f"Error creating heatmap: {str(e)}")
        # Restore original figure size
        plt.rcParams["figure.figsize"] = plt.rcParams.get(
            "figure.figsize", original_figsize
        )
        raise RuntimeError(f"Heatmap creation failed: {str(e)}")


def plot_composition(
    adata: sc.AnnData,
    group_key: str,
    stack_key: str,
    normalize: bool = True,
    kind: Literal["stacked", "percent", "nested"] = "stacked",
    palette: Optional[Union[str, Dict[str, str]]] = None,
    figsize: Optional[Tuple[float, float]] = None,
    rotation: float = 0,
    show_counts: bool = False,
    legend_loc: str = "right",
    title: Optional[str] = None,
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """
    Generate a composition plot showing the distribution of cell types across groups.

    This function creates a visualization of how cell types (or other categorical variables)
    are distributed across different groups, such as clusters, conditions, or time points.

    Args:
        adata: AnnData object containing the data
        group_key: Key in adata.obs for the x-axis groups (e.g., 'leiden', 'condition')
        stack_key: Key in adata.obs for the stacked bars (e.g., 'cell_type')
        normalize: If True, show proportions; otherwise, show raw counts
        kind: Type of plot to create:
            - 'stacked': Standard stacked bar chart
            - 'percent': 100% stacked bar chart
            - 'nested': Nested pie charts
        palette: Color palette for stack_key categories
        figsize: Figure size (width, height) in inches
        rotation: Rotation angle for x-tick labels
        show_counts: Whether to display count values on bars
        legend_loc: Location of the legend ('right', 'bottom', etc.)
        title: Plot title
        save: Path to save the figure (if not None)
        dpi: Resolution for saved figure
        show: Whether to display the plot
        **kwargs: Additional arguments passed to plotting functions

    Returns:
        Figure object for further customization

    Examples:
        >>> # Standard stacked bar chart of cell types across clusters
        >>> plot_composition(adata, 'leiden', 'cell_type')
        >>>
        >>> # 100% stacked bars showing condition composition
        >>> plot_composition(
        ...     adata,
        ...     'sample',
        ...     'cell_type',
        ...     kind='percent',
        ...     rotation=45,
        ...     palette='Set3'
        ... )
        >>>
        >>> # Nested pie charts of cell states across conditions
        >>> plot_composition(
        ...     adata,
        ...     'condition',
        ...     'cell_state',
        ...     kind='nested',
        ...     normalize=False,
        ...     show_counts=True,
        ...     title='Cell state distribution by condition'
        ... )
    """
    log.info(f"Creating composition plot of '{stack_key}' across '{group_key}'")

    # Check that keys exist
    for key in [group_key, stack_key]:
        if key not in adata.obs.columns:
            log.error(f"Column '{key}' not found in adata.obs")
            raise ValueError(f"Column '{key}' not found in adata.obs")

    # Create contingency table
    if normalize and kind != "percent":
        # Normalize by row (each group sums to 1)
        composition = pd.crosstab(
            adata.obs[group_key], adata.obs[stack_key], normalize="index"
        )
    else:
        # Raw counts
        composition = pd.crosstab(adata.obs[group_key], adata.obs[stack_key])

    # Check if we have valid data
    if composition.empty:
        log.error("No valid data for composition plot")
        raise ValueError("No valid data for composition plot")

    # Calculate figsize if not provided
    if figsize is None:
        if kind == "nested":
            # For nested pie charts
            n_groups = len(composition.index)
            figsize = (max(6, min(12, n_groups * 2)), 6)
        else:
            # For bar charts
            n_groups = len(composition.index)
            n_stacks = len(composition.columns)
            figsize = (max(6, min(15, n_groups * 0.5)), max(5, min(10, n_stacks * 0.3)))

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Generate the plot based on kind
    if kind == "stacked" or kind == "percent":
        # Stacked bar charts
        if kind == "percent" and not normalize:
            # Convert to percentages (each column will sum to 100%)
            composition = composition.div(composition.sum(axis=1), axis=0) * 100

        composition.plot.bar(
            stacked=True,
            ax=ax,
            colormap=palette if isinstance(palette, str) else None,
            color=palette if isinstance(palette, dict) else None,
            **kwargs,
        )

        # Set labels
        if normalize and kind != "percent":
            ax.set_ylabel("Proportion")
        elif kind == "percent":
            ax.set_ylabel("Percentage")
        else:
            ax.set_ylabel("Count")

        # Add count labels if requested
        if show_counts:
            for i, (idx, row) in enumerate(composition.iterrows()):
                total = row.sum()
                bottom = 0
                for j, val in enumerate(row):
                    if val > 0:  # Only label non-zero values
                        # Calculate position
                        if normalize and kind != "percent":
                            height = val
                            label = f"{val:.2f}"
                        elif kind == "percent":
                            height = val / 100 if val <= 1 else val
                            label = f"{val:.1f}%"
                        else:
                            height = val
                            label = str(int(val))

                        # Add label to bar
                        if height > 0.05:  # Only label bars with enough height
                            ax.text(
                                i,
                                bottom + height / 2,
                                label,
                                ha="center",
                                va="center",
                                fontsize=8,
                                fontweight="bold",
                                color="white",
                            )

                        bottom += height

    elif kind == "nested":
        # Nested pie charts
        n_groups = len(composition.index)
        n_stacks = len(composition.columns)

        # Calculate grid layout
        if n_groups <= 3:
            n_cols = n_groups
            n_rows = 1
        else:
            n_cols = 3
            n_rows = int(np.ceil(n_groups / 3))

        # Create new figure with subplots
        plt.close(fig)
        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
        axes = axes.flatten()

        # Generate a color palette if not provided
        if palette is None:
            palette = sns.color_palette("husl", n_stacks)
            colors = {cat: palette[i] for i, cat in enumerate(composition.columns)}
        elif isinstance(palette, str):
            palette = sns.color_palette(palette, n_stacks)
            colors = {cat: palette[i] for i, cat in enumerate(composition.columns)}
        else:
            colors = palette

        # Create a pie chart for each group
        for i, (group, row) in enumerate(composition.iterrows()):
            if i < len(axes):
                # Filter out zero values
                non_zero = row[row > 0]

                # Create pie chart
                wedges, texts, autotexts = axes[i].pie(
                    non_zero,
                    autopct=lambda pct: f"{pct:.1f}%" if pct > 5 else "",
                    startangle=90,
                    colors=[colors.get(cat, "gray") for cat in non_zero.index],
                )

                # Customize text
                for autotext in autotexts:
                    autotext.set_fontsize(8)
                    autotext.set_fontweight("bold")

                axes[i].set_title(group)

                # Add count if requested
                if show_counts:
                    total = row.sum()
                    axes[i].text(
                        0,
                        0,
                        f"n={int(total)}",
                        ha="center",
                        va="center",
                        fontweight="bold",
                        fontsize=10,
                    )

        # Hide empty axes
        for i in range(n_groups, len(axes)):
            axes[i].axis("off")

        # Create a single legend for all pie charts
        handles = [
            Patch(color=colors.get(cat, "gray"), label=cat)
            for cat in composition.columns
        ]
        fig.legend(
            handles=handles,
            loc="center right" if legend_loc == "right" else "lower center",
            bbox_to_anchor=(1.05, 0.5) if legend_loc == "right" else (0.5, 0),
            title=stack_key,
        )

        plt.tight_layout()

    else:
        log.error(f"Unknown plot kind: {kind}")
        raise ValueError(f"Unknown plot kind: {kind}")

    # Rotate x-tick labels if needed
    if kind != "nested" and rotation != 0:
        plt.xticks(rotation=rotation, ha="right" if rotation > 30 else "center")

    # Add legend
    if kind != "nested":
        if legend_loc == "right":
            plt.legend(title=stack_key, bbox_to_anchor=(1.05, 1), loc="upper left")
        else:
            plt.legend(title=stack_key, loc=legend_loc)

    # Add title
    if title:
        plt.title(title)
    else:
        plt.title(f"Composition of {stack_key} across {group_key}")

    # Finalize layout
    plt.tight_layout()

    # Save if requested
    if save:
        plt.savefig(save, dpi=dpi, bbox_inches="tight")

    # Show if requested
    if show:
        plt.show()

    return fig


def plot_enrichment(
    adata: sc.AnnData,
    cluster: Union[str, List[str]],
    enrich_key: str = "enrichment",
    n_top_terms: int = 10,
    figsize: Optional[Tuple[float, float]] = None,
    color: str = "steelblue",
    title_fontsize: int = 14,
    label_fontsize: int = 10,
    max_label_length: int = 50,
    save_dir: Optional[str] = None,
    show: bool = True,
    **kwargs,
) -> Union[plt.Figure, List[plt.Figure]]:
    """
    Visualize functional enrichment analysis results for clusters.

    This function creates bar plots showing the most significantly enriched functional
    terms (e.g., GO terms, pathways) for specific clusters or cell types.

    Args:
        adata: AnnData object containing enrichment results
        cluster: Name(s) of the cluster(s) to visualize
        enrich_key: Key in adata.uns storing enrichment analysis results
        n_top_terms: Number of top terms to display
        figsize: Figure size (width, height) in inches
        color: Bar color
        title_fontsize: Font size for plot title
        label_fontsize: Font size for term labels
        max_label_length: Maximum length for term labels before truncation
        save_dir: Directory to save plot files (if not None)
        show: Whether to display the image
        **kwargs: Additional arguments passed to plotting functions

    Returns:
        Figure object(s) for further customization

    Examples:
        >>> # Visualize enrichment for a single cluster
        >>> plot_enrichment(adata, 'Cluster_1')
        >>>
        >>> # Visualize enrichment for multiple clusters and save results
        >>> plot_enrichment(
        ...     adata,
        ...     ['T cells', 'B cells', 'NK cells'],
        ...     n_top_terms=15,
        ...     save_dir='./enrichment_plots/'
        ... )
        >>>
        >>> # Customize appearance
        >>> plot_enrichment(
        ...     adata,
        ...     'Macrophages',
        ...     color='firebrick',
        ...     max_label_length=30,
        ...     title_fontsize=16
        ... )
    """
    # Handle multiple clusters
    if isinstance(cluster, list):
        figures = []
        for c in cluster:
            try:
                fig = plot_enrichment(
                    adata,
                    c,
                    enrich_key,
                    n_top_terms,
                    figsize,
                    color,
                    title_fontsize,
                    label_fontsize,
                    max_label_length,
                    save_dir,
                    show,
                    **kwargs,
                )
                figures.append(fig)
            except Exception as e:
                log.warning(f"Error plotting enrichment for cluster '{c}': {str(e)}")
        return figures

    # Check if enrichment results exist
    if enrich_key not in adata.uns or cluster not in adata.uns[enrich_key]:
        log.error(
            f"Enrichment results for cluster '{cluster}' not found in `adata.uns['{enrich_key}']`"
        )
        raise KeyError(
            f"Enrichment results for cluster '{cluster}' not found in `adata.uns['{enrich_key}']`. "
            "Please run `run_enrichment()` first."
        )

    # Get enrichment results
    results_df = adata.uns[enrich_key][cluster]

    if results_df.empty:
        log.warning(f"No enrichment results to plot for cluster '{cluster}'")
        return None

    # Get top terms
    top_terms = results_df.head(n_top_terms)

    # Truncate long term names
    if max_label_length > 0:
        top_terms = top_terms.copy()
        top_terms["Term_Short"] = top_terms["Term"].apply(
            lambda x: x
            if len(x) <= max_label_length
            else x[: max_label_length - 3] + "..."
        )
    else:
        top_terms["Term_Short"] = top_terms["Term"]

    # Calculate figsize if not provided
    if figsize is None:
        figsize = (10, max(5, min(15, n_top_terms * 0.4)))

    # Create plot
    fig, ax = plt.subplots(figsize=figsize)

    # Plot bars
    y_pos = np.arange(len(top_terms))
    bars = ax.barh(
        y_pos, -np.log10(top_terms["Adjusted P-value"]), color=color, **kwargs
    )

    # Add term labels
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_terms["Term_Short"], fontsize=label_fontsize)

    # Add p-value labels
    for i, bar in enumerate(bars):
        width = bar.get_width()
        p_val = top_terms["Adjusted P-value"].iloc[i]
        if p_val < 0.001:
            p_label = "p < 0.001"
        else:
            p_label = f"p = {p_val:.3f}"

        ax.text(
            width + 0.1,
            bar.get_y() + bar.get_height() / 2,
            p_label,
            va="center",
            fontsize=8,
            color="grey",
        )

    # Add counts to bars
    for i, bar in enumerate(bars):
        gene_count = (
            f"{top_terms['Genes_found'].iloc[i]}/{top_terms['Genes_in_Term'].iloc[i]}"
        )
        ax.text(
            0.1,
            bar.get_y() + bar.get_height() / 2,
            gene_count,
            va="center",
            fontsize=8,
            color="white",
            fontweight="bold",
        )

    # Customize plot
    ax.set_title(f"Enriched Terms for {cluster}", fontsize=title_fontsize)
    ax.set_xlabel("-log10(Adjusted P-value)", fontsize=12)
    ax.invert_yaxis()  # Terms with smallest p-value at the top
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.grid(axis="x", linestyle="--", alpha=0.6)

    plt.tight_layout()

    # Save if requested
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        safe_cluster = cluster.replace("/", "_").replace(" ", "_")
        save_path = os.path.join(save_dir, f"{safe_cluster}_enrichment.png")
        plt.savefig(save_path, dpi=300, bbox_inches="tight", facecolor="white")
        log.info(f"Saved enrichment plot to {save_path}")

    # Show if requested
    if show:
        plt.show()
    else:
        plt.close()

    return fig


def plot_pseudotime(
    adata: sc.AnnData,
    pseudotime_key: str,
    embedding_key: str = "X_umap",
    color_key: Optional[str] = None,
    lineage_key: Optional[str] = None,
    cluster_key: Optional[str] = None,
    cmap: str = "viridis",
    lineage_colors: Optional[Union[str, List[str]]] = None,
    background: bool = True,
    background_alpha: float = 0.1,
    point_size: float = 10,
    figsize: Optional[Tuple[float, float]] = None,
    ncols: int = 3,
    colorbar: bool = True,
    arrow: bool = False,
    arrow_color: str = "black",
    arrow_size: int = 100,
    arrow_density: int = 30,
    title: Optional[str] = None,
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """
    Visualize pseudotime and trajectory analysis results on embeddings.

    This function creates visualizations of pseudotemporal ordering and cell trajectories,
    supporting various trajectory inference methods (Slingshot, PAGA, Monocle, etc.).

    Args:
        adata: AnnData object containing the data
        pseudotime_key: Key in adata.obs containing pseudotime values
        embedding_key: Key in adata.obsm for the embedding coordinates (e.g., 'X_umap')
        color_key: Additional observation to color points by (optional)
        lineage_key: Key in adata.obs containing lineage/branch assignments
        cluster_key: Key in adata.obs for cluster annotations (for reference)
        cmap: Colormap for pseudotime values
        lineage_colors: Colors for different lineages/branches
        background: Whether to show all cells in the background
        background_alpha: Transparency for background cells
        point_size: Size of the scatter points
        figsize: Figure size (width, height) in inches
        ncols: Number of columns for multi-panel plots
        colorbar: Whether to show the colorbar
        arrow: Whether to show velocity arrows (requires velocity data)
        arrow_color: Color for velocity arrows
        arrow_size: Size of arrow heads
        arrow_density: Density of arrows on the plot
        title: Plot title
        save: Path to save the figure (if not None)
        dpi: Resolution for saved figure
        show: Whether to display the plot
        **kwargs: Additional arguments passed to plotting functions

    Returns:
        Figure object for further customization

    Examples:
        >>> # Basic pseudotime visualization
        >>> plot_pseudotime(adata, 'slingshot_pseudotime')
        >>>
        >>> # Visualize multiple lineages
        >>> plot_pseudotime(
        ...     adata,
        ...     'dpt_pseudotime',
        ...     lineage_key='dpt_branches',
        ...     cluster_key='leiden'
        ... )
        >>>
        >>> # Add RNA velocity arrows to trajectory
        >>> plot_pseudotime(
        ...     adata,
        ...     'monocle3_pseudotime',
        ...     arrow=True,
        ...     arrow_density=50,
        ...     background=True
        ... )
    """
    log.info(f"Creating pseudotime visualization for '{pseudotime_key}'")

    # Check required keys
    if pseudotime_key not in adata.obs:
        log.error(f"Pseudotime key '{pseudotime_key}' not found in adata.obs")
        raise ValueError(f"Pseudotime key '{pseudotime_key}' not found in adata.obs")

    if embedding_key not in adata.obsm:
        log.error(f"Embedding key '{embedding_key}' not found in adata.obsm")
        raise ValueError(f"Embedding key '{embedding_key}' not found in adata.obsm")

    # Extract embedding coordinates
    X_embed = adata.obsm[embedding_key]

    # Extract pseudotime values
    pseudotime = adata.obs[pseudotime_key]

    # Check for missing values
    if pseudotime.isna().any():
        log.warning(
            f"Pseudotime contains {pseudotime.isna().sum()} missing values. These will be excluded."
        )
        has_pseudotime = ~pseudotime.isna()
    else:
        has_pseudotime = np.ones(len(pseudotime), dtype=bool)

    # Handle lineage information if provided
    if lineage_key is not None:
        if lineage_key not in adata.obs:
            log.error(f"Lineage key '{lineage_key}' not found in adata.obs")
            raise ValueError(f"Lineage key '{lineage_key}' not found in adata.obs")

        lineage = adata.obs[lineage_key].astype(str)

        # Handle missing lineage values
        if lineage.isna().any():
            log.warning(
                f"Lineage contains {lineage.isna().sum()} missing values. These will be excluded."
            )
            has_lineage = ~lineage.isna()
            has_data = has_pseudotime & has_lineage
        else:
            has_data = has_pseudotime

        # Get unique lineages (excluding NaN)
        unique_lineages = np.unique(lineage[~lineage.isna()])
        n_lineages = len(unique_lineages)

        # Generate colors for lineages if not provided
        if lineage_colors is None:
            lineage_colors = sns.color_palette("tab10", n_lineages)
        elif isinstance(lineage_colors, str):
            lineage_colors = sns.color_palette(lineage_colors, n_lineages)

        # Create a dictionary mapping lineages to colors
        lineage_color_dict = {
            l: lineage_colors[i] for i, l in enumerate(unique_lineages)
        }

        # Determine if we need a multi-panel plot
        multi_panel = True
    else:
        has_data = has_pseudotime
        multi_panel = False

    # Calculate figsize if not provided
    if figsize is None:
        if multi_panel:
            n_panels = n_lineages + 1  # All lineages + combined
            n_cols = min(ncols, n_panels)
            n_rows = int(np.ceil(n_panels / n_cols))
            figsize = (5 * n_cols, 5 * n_rows)
        else:
            figsize = (10, 8)

    # Create plots
    if multi_panel:
        # Multi-panel plot for lineages
        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
        axes = axes.flatten()

        # First panel: Combined view
        ax = axes[0]

        # Show all cells in background if requested
        if background:
            ax.scatter(
                X_embed[:, 0],
                X_embed[:, 1],
                s=point_size * 0.5,
                c="lightgray",
                alpha=background_alpha,
            )

        # Color points by pseudotime and lineage
        for lin in unique_lineages:
            mask = (lineage == lin) & has_data
            if np.sum(mask) > 0:
                sc = ax.scatter(
                    X_embed[mask, 0],
                    X_embed[mask, 1],
                    s=point_size,
                    c=pseudotime[mask],
                    cmap=cmap,
                    edgecolor=lineage_color_dict[lin],
                    linewidth=0.5,
                    alpha=0.8,
                )

        ax.set_title("All Lineages")

        # Add colorbar
        if colorbar:
            cbar = plt.colorbar(sc, ax=ax)
            cbar.set_label("Pseudotime")

        # Add arrows if requested
        if arrow:
            _add_velocity_arrows(
                adata, ax, density=arrow_density, size=arrow_size, color=arrow_color
            )

        # Individual panels for each lineage
        for i, lin in enumerate(unique_lineages):
            if i + 1 < len(axes):
                ax = axes[i + 1]

                # Show all cells in background if requested
                if background:
                    ax.scatter(
                        X_embed[:, 0],
                        X_embed[:, 1],
                        s=point_size * 0.5,
                        c="lightgray",
                        alpha=background_alpha,
                    )

                # Select cells in this lineage
                mask = (lineage == lin) & has_data

                if np.sum(mask) > 0:
                    sc = ax.scatter(
                        X_embed[mask, 0],
                        X_embed[mask, 1],
                        s=point_size,
                        c=pseudotime[mask],
                        cmap=cmap,
                        alpha=0.8,
                    )

                    ax.set_title(f"Lineage {lin}")

                    # Add colorbar
                    if colorbar:
                        cbar = plt.colorbar(sc, ax=ax)
                        cbar.set_label("Pseudotime")

                    # Add arrows if requested
                    if arrow:
                        _add_velocity_arrows(
                            adata,
                            ax,
                            mask=mask,
                            density=arrow_density,
                            size=arrow_size,
                            color=arrow_color,
                        )
                else:
                    ax.text(
                        0.5,
                        0.5,
                        f"No cells in lineage {lin}",
                        ha="center",
                        va="center",
                        transform=ax.transAxes,
                    )

        # Hide empty axes
        for i in range(n_panels, len(axes)):
            axes[i].axis("off")

    else:
        # Single panel plot
        fig, ax = plt.subplots(figsize=figsize)

        # Show all cells in background if requested
        if background and not np.all(has_data):
            ax.scatter(
                X_embed[~has_data, 0],
                X_embed[~has_data, 1],
                s=point_size * 0.5,
                c="lightgray",
                alpha=background_alpha,
            )

        # Color points by pseudotime
        sc = ax.scatter(
            X_embed[has_data, 0],
            X_embed[has_data, 1],
            s=point_size,
            c=pseudotime[has_data],
            cmap=cmap,
            alpha=0.8,
        )

        # Add colorbar
        if colorbar:
            cbar = plt.colorbar(sc, ax=ax)
            cbar.set_label("Pseudotime")

        # Add arrows if requested
        if arrow:
            _add_velocity_arrows(
                adata,
                ax,
                mask=has_data,
                density=arrow_density,
                size=arrow_size,
                color=arrow_color,
            )

    # Set axis labels for all axes
    for ax in fig.axes:
        if not ax.get_subplotspec().is_last_row():
            ax.set_xlabel("")
        else:
            ax.set_xlabel(f"{embedding_key.split('_')[1].upper()}1")

        if not ax.get_subplotspec().is_first_col():
            ax.set_ylabel("")
        else:
            ax.set_ylabel(f"{embedding_key.split('_')[1].upper()}2")

    # Set title
    if title:
        fig.suptitle(title, fontsize=16)

    plt.tight_layout()

    # Save if requested
    if save:
        plt.savefig(save, dpi=dpi, bbox_inches="tight")
        log.info(f"Saved pseudotime plot to {save}")

    # Show if requested
    if show:
        plt.show()

    return fig


def _add_velocity_arrows(
    adata: sc.AnnData,
    ax: plt.Axes,
    basis: str = "umap",
    mask: Optional[np.ndarray] = None,
    density: int = 30,
    size: float = 120,
    color: str = "black",
    alpha: float = 0.8,
    scale: Optional[float] = None,
    **kwargs,
) -> None:
    """
    Helper function to add velocity arrows to trajectory plots.

    Args:
        adata: AnnData object with velocity information
        ax: Matplotlib axes to add arrows to
        basis: Basis for embedding and velocity projection
        mask: Boolean mask to select subset of cells
        density: Density of arrows
        size: Size of arrowheads
        color: Color of arrows
        alpha: Transparency of arrows
        scale: Scaling factor for arrows
        **kwargs: Additional arguments for plt.quiver
    """
    vel_keys = [
        f"{basis}_velocity",  # Scanpy/scVelo
        f"velocity_{basis}",  # Alternative format
        "velocity",  # Generic
    ]

    # Check if velocity data exists
    vel_key = None
    for key in vel_keys:
        if key in adata.obsm:
            vel_key = key
            break

    if vel_key is None:
        log.warning(
            "No velocity data found in adata.obsm. Velocity arrows will not be shown."
        )
        return

    embed_key = f"X_{basis}"
    if embed_key not in adata.obsm:
        log.warning(
            f"Embedding '{embed_key}' not found in adata.obsm. Velocity arrows will not be shown."
        )
        return

    # Extract coordinates and velocities
    X_emb = adata.obsm[embed_key]
    V_emb = adata.obsm[vel_key]

    # Apply mask if provided
    if mask is not None:
        X_emb = X_emb[mask]
        V_emb = V_emb[mask]

    # Sample points for arrows based on density
    n_cells = X_emb.shape[0]
    n_arrows = min(n_cells, density * density)

    # Grid-based sampling for more uniform coverage
    grid_size = int(np.sqrt(n_arrows))

    # Get min/max coordinates
    xmin, xmax = X_emb[:, 0].min(), X_emb[:, 0].max()
    ymin, ymax = X_emb[:, 1].min(), X_emb[:, 1].max()

    # Create grid
    grid_x = np.linspace(xmin, xmax, grid_size)
    grid_y = np.linspace(ymin, ymax, grid_size)

    # Find nearest cells to grid points
    arrows_x = []
    arrows_y = []
    arrows_u = []
    arrows_v = []

    from scipy.spatial import cKDTree

    tree = cKDTree(X_emb)

    for ix in grid_x:
        for iy in grid_y:
            # Find closest cell to this grid point
            dist, idx = tree.query([ix, iy], k=1)

            # Only add arrow if there's a cell near this grid point
            if dist < (xmax - xmin) / grid_size:
                arrows_x.append(X_emb[idx, 0])
                arrows_y.append(X_emb[idx, 1])
                arrows_u.append(V_emb[idx, 0])
                arrows_v.append(V_emb[idx, 1])

    # Normalize velocity vectors
    if len(arrows_u) > 0:
        # Calculate vector magnitudes
        magnitudes = np.sqrt(np.array(arrows_u) ** 2 + np.array(arrows_v) ** 2)

        # Normalize vectors
        arrows_u = np.array(arrows_u) / (magnitudes + 1e-10)
        arrows_v = np.array(arrows_v) / (magnitudes + 1e-10)

        # Apply scaling
        if scale is None:
            scale = (xmax - xmin) / 50  # Automatic scaling based on plot size

        arrows_u *= scale
        arrows_v *= scale

        # Add arrows to plot
        ax.quiver(
            arrows_x,
            arrows_y,
            arrows_u,
            arrows_v,
            color=color,
            alpha=alpha,
            scale=1,
            scale_units="xy",
            width=0.001,
            headwidth=size / 10,
            headlength=size / 7,
            headaxislength=size / 7,
            **kwargs,
        )


def plot_spatial(
    adata: sc.AnnData,
    color_by: Union[str, List[str]],
    use_raw: bool = False,
    layer: Optional[str] = None,
    basis: str = "spatial",
    library_id: Optional[str] = None,
    spot_size: Optional[float] = None,
    cmap: str = "viridis",
    palette: Optional[Union[str, Dict[str, str]]] = None,
    alpha: float = 1.0,
    img: bool = True,
    img_alpha: float = 1.0,
    crop: bool = True,
    colorbar: bool = True,
    figsize: Optional[Tuple[float, float]] = None,
    ncols: int = 2,
    legend_loc: str = "right margin",
    title: Optional[str] = None,
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """
    Visualize spatial transcriptomics data with enhanced customization.

    This function creates publication-quality visualizations of spatial transcriptomics
    data, supporting both gene expression and categorical variables.

    Args:
        adata: AnnData object with spatial coordinates
        color_by: Key(s) in adata.obs or adata.var_names for coloring spots
        use_raw: Whether to use adata.raw for gene expression
        layer: Layer in adata.layers to use for gene expression
        basis: Key in adata.obsm containing spatial coordinates
        library_id: Library ID for adata.uns['spatial']
        spot_size: Size of spots (if None, determined automatically)
        cmap: Colormap for continuous variables
        palette: Color palette for categorical variables
        alpha: Transparency of spots
        img: Whether to show the tissue image in the background
        img_alpha: Transparency of the background image
        crop: Whether to crop the plot to the tissue region
        colorbar: Whether to show colorbar for continuous variables
        figsize: Figure size (width, height) in inches
        ncols: Number of columns for multiple plots
        legend_loc: Location of the legend
        title: Plot title
        save: Path to save the figure (if not None)
        dpi: Resolution for saved figure
        show: Whether to display the plot
        **kwargs: Additional arguments passed to sc.pl.spatial

    Returns:
        Figure object for further customization

    Examples:
        >>> # Visualize expression of a single gene
        >>> plot_spatial(adata, 'MS4A1')
        >>>
        >>> # Visualize multiple features
        >>> plot_spatial(
        ...     adata,
        ...     ['MS4A1', 'CD3D', 'leiden'],
        ...     spot_size=1.2,
        ...     ncols=3,
        ...     img_alpha=0.7
        ... )
        >>>
        >>> # Customize colors for a categorical variable
        >>> plot_spatial(
        ...     adata,
        ...     'cluster',
        ...     palette={'Tumor': 'red', 'Immune': 'blue', 'Stroma': 'green'},
        ...     save='spatial_clusters.png'
        ... )
    """
    log.info(f"Creating spatial plot for '{color_by}'")

    # Validate input
    if f"X_{basis}" not in adata.obsm:
        log.error(f"Spatial coordinates not found in adata.obsm['X_{basis}']")
        raise ValueError(f"Spatial coordinates not found in adata.obsm['X_{basis}']")

    # Convert single color to list
    if isinstance(color_by, str):
        color_by = [color_by]

    # Check if all color_by items exist
    for color in color_by:
        if color not in adata.obs.columns and color not in adata.var_names:
            log.error(f"'{color}' not found in adata.obs.columns or adata.var_names")
            raise ValueError(
                f"'{color}' not found in adata.obs.columns or adata.var_names"
            )

    # Determine library_id if not provided
    if img and library_id is None:
        if "spatial" in adata.uns:
            library_ids = list(adata.uns["spatial"].keys())
            if len(library_ids) == 1:
                library_id = library_ids[0]
                log.info(f"Using library_id: {library_id}")
            elif len(library_ids) > 1:
                log.warning(
                    f"Multiple library_ids found: {library_ids}. Using the first one: {library_ids[0]}"
                )
                library_id = library_ids[0]
            else:
                log.warning("No library_ids found in adata.uns['spatial']")
                img = False
        else:
            log.warning("No spatial image data found in adata.uns['spatial']")
            img = False

    # Calculate appropriate figsize if not provided
    if figsize is None:
        n_colors = len(color_by)
        n_cols = min(ncols, n_colors)
        n_rows = int(np.ceil(n_colors / n_cols))
        figsize = (6 * n_cols, 6 * n_rows)

    try:
        # Store original settings
        orig_figsize = plt.rcParams["figure.figsize"]
        plt.rcParams["figure.figsize"] = figsize

        # Create spatial plot
        sc.pl.spatial(
            adata,
            color=color_by,
            use_raw=use_raw,
            layer=layer,
            library_id=library_id,
            spot_size=spot_size,
            img=img,
            img_alpha=img_alpha,
            alpha=alpha,
            crop=crop,
            show=False,
            cmap=cmap,
            palette=palette,
            colorbar=colorbar,
            ncols=ncols,
            legend_loc=legend_loc,
            title=title,
            **kwargs,
        )

        # Get current figure
        fig = plt.gcf()

        # Restore original settings
        plt.rcParams["figure.figsize"] = orig_figsize

        # Save if requested
        if save:
            plt.savefig(save, dpi=dpi, bbox_inches="tight")
            log.info(f"Saved spatial plot to {save}")

        # Show if requested
        if show:
            plt.show()
        else:
            plt.close()

        return fig

    except Exception as e:
        # Restore original settings in case of error
        plt.rcParams["figure.figsize"] = plt.rcParams.get(
            "figure.figsize", orig_figsize
        )
        log.error(f"Error creating spatial plot: {str(e)}")
        raise RuntimeError(f"Spatial plot creation failed: {str(e)}")


def plot_volcano(
    df: pd.DataFrame,
    x: str = "logfoldchanges",
    y: str = "pvals_adj",
    gene_col: str = "names",
    group_col: Optional[str] = "group",
    group: Optional[str] = None,
    highlight_genes: Optional[List[str]] = None,
    min_log2fc: float = 1.0,
    max_pval: float = 0.05,
    figsize: Optional[Tuple[float, float]] = None,
    cmap: str = "coolwarm",
    point_size: float = 10,
    max_overlaps: int = 25,
    label_top_n: int = 10,
    title: Optional[str] = None,
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """
    Create a publication-quality volcano plot from differential expression results.

    This function visualizes differential expression results with customizable highlighting
    of significant genes and optional gene labels.

    Args:
        df: DataFrame containing differential expression results
        x: Column name for log2 fold change values
        y: Column name for p-values
        gene_col: Column name for gene names
        group_col: Column name for group information
        group: Specific group to plot (if None, plot all)
        highlight_genes: List of genes to highlight
        min_log2fc: Minimum absolute log2 fold change for significance
        max_pval: Maximum adjusted p-value for significance
        figsize: Figure size (width, height) in inches
        cmap: Colormap for fold change values
        point_size: Size of scatter points
        max_overlaps: Maximum overlaps for text adjustment
        label_top_n: Number of top genes to label
        title: Plot title
        save: Path to save the figure (if not None)
        dpi: Resolution for saved figure
        show: Whether to display the plot
        **kwargs: Additional arguments passed to plotting functions

    Returns:
        Figure object for further customization

    Examples:
        >>> # Basic volcano plot from find_markers() results
        >>> markers_df = find_markers(adata, groupby='leiden')
        >>> plot_volcano(markers_df)
        >>>
        >>> # Volcano plot for a specific group with custom thresholds
        >>> plot_volcano(
        ...     markers_df,
        ...     group='CD8+ T cells',
        ...     min_log2fc=2.0,
        ...     max_pval=0.01,
        ...     label_top_n=15
        ... )
        >>>
        >>> # Highlight specific genes of interest
        >>> plot_volcano(
        ...     markers_df,
        ...     highlight_genes=['IFNG', 'IL2', 'TNF', 'GZMB'],
        ...     save='cytokine_volcano.png'
        ... )
    """
    log.info("Creating volcano plot")

    # Filter for specific group if requested
    if group is not None and group_col in df.columns:
        df = df[df[group_col] == group].copy()
        log.info(f"Filtered for group: {group}")

    # Check if required columns exist
    for col in [x, y, gene_col]:
        if col not in df.columns:
            log.error(f"Required column '{col}' not found in DataFrame")
            raise ValueError(f"Required column '{col}' not found in DataFrame")

    # Filter out rows with missing values
    df = df.dropna(subset=[x, y, gene_col])

    if df.empty:
        log.error("No valid data for volcano plot after filtering")
        raise ValueError("No valid data for volcano plot after filtering")

    # Transform p-values to -log10 scale
    df["neg_log10_pval"] = -np.log10(df[y].clip(lower=1e-300))  # Clip to avoid infinity

    # Determine significance
    df["significant"] = (df[y] < max_pval) & (abs(df[x]) > min_log2fc)
    df["regulation"] = "Not significant"
    df.loc[(df[y] < max_pval) & (df[x] > min_log2fc), "regulation"] = "Up"
    df.loc[(df[y] < max_pval) & (df[x] < -min_log2fc), "regulation"] = "Down"

    # Determine point colors
    color_dict = {"Up": "red", "Down": "blue", "Not significant": "gray"}

    # Calculate figsize if not provided
    if figsize is None:
        figsize = (10, 8)

    # Create plot
    fig, ax = plt.subplots(figsize=figsize)

    # Get color values based on fold change for colormap
    if cmap is not None:
        colors = None
        cmap_obj = plt.get_cmap(cmap)

        # Normalize fold changes to [-1, 1] for colormap
        vmin, vmax = (
            -max(abs(df[x].min()), abs(df[x].max())),
            max(abs(df[x].min()), abs(df[x].max())),
        )
        norm = plt.Normalize(vmin, vmax)

        # Create scatter plot with colormap
        sc = ax.scatter(
            df[x],
            df["neg_log10_pval"],
            s=point_size,
            c=df[x],
            cmap=cmap_obj,
            norm=norm,
            alpha=0.7,
            **kwargs,
        )

        # Add colorbar
        cbar = plt.colorbar(sc, ax=ax)
        cbar.set_label("Log2 Fold Change")

    else:
        # Create scatter plot with categorical colors
        for regulation, color in color_dict.items():
            subset = df[df["regulation"] == regulation]
            ax.scatter(
                subset[x],
                subset["neg_log10_pval"],
                s=point_size,
                c=color,
                alpha=0.7,
                label=regulation,
                **kwargs,
            )

    # Add horizontal line for p-value threshold
    ax.axhline(y=-np.log10(max_pval), color="gray", linestyle="--", alpha=0.5)
    ax.text(
        df[x].max() * 0.98,
        -np.log10(max_pval) * 1.02,
        f"p-value = {max_pval}",
        ha="right",
        va="bottom",
        color="gray",
        fontsize=8,
    )

    # Add vertical lines for fold change thresholds
    ax.axvline(x=min_log2fc, color="gray", linestyle="--", alpha=0.5)
    ax.axvline(x=-min_log2fc, color="gray", linestyle="--", alpha=0.5)

    # Add labels
    genes_to_label = set()

    # 1. Add user-specified genes
    if highlight_genes is not None:
        genes_to_label.update(highlight_genes)

    # 2. Add top significant genes
    if label_top_n > 0:
        # Get top up-regulated genes
        top_up = (
            df[df["regulation"] == "Up"]
            .nlargest(label_top_n // 2, x)[gene_col]
            .tolist()
        )
        genes_to_label.update(top_up)

        # Get top down-regulated genes
        top_down = (
            df[df["regulation"] == "Down"]
            .nsmallest(label_top_n // 2, x)[gene_col]
            .tolist()
        )
        genes_to_label.update(top_down)

        # Get top significant genes by p-value
        top_sig = (
            df[df["significant"]]
            .nlargest(label_top_n // 2, "neg_log10_pval")[gene_col]
            .tolist()
        )
        genes_to_label.update(top_sig)

    # Add text labels
    if genes_to_label:
        from adjustText import adjust_text

        # Filter for genes that actually exist in the data
        genes_to_label = [g for g in genes_to_label if g in df[gene_col].values]

        texts = []
        for gene in genes_to_label:
            gene_data = df[df[gene_col] == gene].iloc[0]
            texts.append(
                ax.text(
                    gene_data[x],
                    gene_data["neg_log10_pval"],
                    gene,
                    ha="center",
                    va="center",
                    fontsize=9,
                    fontweight="bold",
                )
            )

        # Adjust text positions to avoid overlaps
        try:
            adjust_text(
                texts,
                arrowprops=dict(arrowstyle="-", color="black", lw=0.5),
                expand_points=(1.5, 1.5),
                force_text=(0.5, 0.5),
                max_overlaps=max_overlaps,
            )
        except Exception as e:
            log.warning(f"Error adjusting text positions: {str(e)}")

    # Add legend for categorical colors
    if cmap is None:
        ax.legend(title="Regulation")

    # Set axis labels and title
    ax.set_xlabel("Log2 Fold Change", fontsize=12)
    ax.set_ylabel("-log10(Adjusted p-value)", fontsize=12)

    if title:
        ax.set_title(title, fontsize=14)
    elif group is not None:
        ax.set_title(f"Volcano Plot - {group}", fontsize=14)
    else:
        ax.set_title("Volcano Plot", fontsize=14)

    # Add summary text
    n_up = sum(df["regulation"] == "Up")
    n_down = sum(df["regulation"] == "Down")
    n_total = len(df)

    ax.text(
        0.02,
        0.98,
        f"Total: {n_total}\nUp: {n_up} ({n_up / n_total:.1%})\nDown: {n_down} ({n_down / n_total:.1%})",
        transform=ax.transAxes,
        va="top",
        ha="left",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7),
    )

    # Finalize plot
    plt.tight_layout()

    # Save if requested
    if save:
        plt.savefig(save, dpi=dpi, bbox_inches="tight")
        log.info(f"Saved volcano plot to {save}")

    # Show if requested
    if show:
        plt.show()

    return fig


def plot_feature_correlation(
    adata: sc.AnnData,
    features: List[str],
    method: Literal["pearson", "spearman"] = "spearman",
    layer: Optional[str] = None,
    use_raw: bool = False,
    cluster_features: bool = True,
    cmap: str = "coolwarm",
    vmin: Optional[float] = -1,
    vmax: Optional[float] = 1,
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """
    Create a correlation heatmap for selected features.

    This function visualizes the correlation between genes, proteins, or other features
    in the dataset, supporting both Pearson and Spearman correlation methods.

    Args:
        adata: AnnData object containing the data
        features: List of features (genes/proteins) to correlate
        method: Correlation method ('pearson' or 'spearman')
        layer: Layer in adata.layers to use for expression values
        use_raw: Whether to use adata.raw for expression values
        cluster_features: Whether to cluster features by similarity
        cmap: Colormap for the heatmap
        vmin: Minimum value for color scale
        vmax: Maximum value for color scale
        figsize: Figure size (width, height) in inches
        title: Plot title
        save: Path to save the figure (if not None)
        dpi: Resolution for saved figure
        show: Whether to display the plot
        **kwargs: Additional arguments passed to sns.heatmap

    Returns:
        Figure object for further customization

    Examples:
        >>> # Correlate T cell markers
        >>> plot_feature_correlation(
        ...     adata,
        ...     ['CD3D', 'CD4', 'CD8A', 'FOXP3', 'IL2RA', 'CTLA4']
        ... )
        >>>
        >>> # Correlate genes across conditions using Pearson correlation
        >>> plot_feature_correlation(
        ...     adata,
        ...     ['IFNG', 'TNF', 'IL2', 'GZMB', 'PRF1', 'TNFRSF9'],
        ...     method='pearson',
        ...     title='Cytotoxic gene correlation',
        ...     save='cytotoxic_correlation.png'
        ... )
    """
    log.info(f"Creating feature correlation heatmap for {len(features)} features")

    # Check if features exist
    if use_raw and adata.raw is not None:
        var_names = adata.raw.var_names
    else:
        var_names = adata.var_names

    missing_features = [f for f in features if f not in var_names]
    if missing_features:
        log.warning(f"{len(missing_features)} features not found in dataset")
        if len(missing_features) < 10:
            log.warning(f"Missing features: {', '.join(missing_features)}")
        else:
            log.warning(
                f"First 10 missing features: {', '.join(missing_features[:10])}..."
            )

        # Remove missing features
        features = [f for f in features if f not in missing_features]

        if not features:
            log.error("No valid features remaining for correlation analysis")
            raise ValueError("No valid features remaining for correlation analysis")

    # Extract data for correlation
    if use_raw and adata.raw is not None:
        log.info("Using adata.raw for expression values")
        X = adata.raw[:, features].X
    elif layer is not None:
        log.info(f"Using layer '{layer}' for expression values")
        X = adata[:, features].layers[layer]
    else:
        log.info("Using adata.X for expression values")
        X = adata[:, features].X

    # Convert to dense array if sparse
    if scipy.sparse.issparse(X):
        X = X.toarray()

    # Create dataframe for correlation
    df = pd.DataFrame(X, columns=features)

    # Calculate correlation
    log.info(f"Calculating {method} correlation")
    corr = df.corr(method=method)

    # Calculate figsize if not provided
    if figsize is None:
        n_features = len(features)
        figsize = (max(6, min(15, n_features * 0.5)), max(5, min(15, n_features * 0.5)))

    # Create plot
    fig, ax = plt.subplots(figsize=figsize)

    # Cluster features if requested
    if cluster_features and len(features) > 2:
        try:
            from scipy.cluster import hierarchy
            from scipy.spatial import distance

            # Calculate linkage
            corr_linkage = hierarchy.linkage(distance.pdist(corr), method="average")

            # Get clustered indices
            idx = hierarchy.dendrogram(corr_linkage, no_plot=True)["leaves"]

            # Reorder correlation matrix
            corr = corr.iloc[idx, idx]

            log.info("Features clustered by correlation similarity")
        except Exception as e:
            log.warning(f"Error clustering features: {str(e)}")

    # Create heatmap
    sns.heatmap(
        corr,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        center=0,
        square=True,
        linewidths=0.5,
        annot=len(features) <= 20,  # Only show annotations for small matrices
        fmt=".2f" if len(features) <= 20 else None,
        cbar_kws={"shrink": 0.8},
        ax=ax,
        **kwargs,
    )

    # Add title
    if title:
        ax.set_title(title, fontsize=14)
    else:
        ax.set_title(f"{method.capitalize()} correlation between features", fontsize=14)

    # Rotate axis labels if there are many features
    if len(features) > 10:
        plt.xticks(rotation=90)
        plt.yticks(rotation=0)

    plt.tight_layout()

    # Save if requested
    if save:
        plt.savefig(save, dpi=dpi, bbox_inches="tight")
        log.info(f"Saved correlation heatmap to {save}")

    # Show if requested
    if show:
        plt.show()

    return fig


def plot_multi_modality(
    adata: sc.AnnData,
    mod1: str,
    mod2: str,
    features: List[Tuple[str, str]],
    basis: str = "umap",
    ncols: int = 2,
    figsize: Optional[Tuple[float, float]] = None,
    palette: Optional[str] = None,
    cmap: str = "viridis",
    size: float = 10,
    alpha: float = 0.8,
    show_correlation: bool = True,
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """
    Visualize and compare features across multiple modalities.

    This function creates paired visualizations of features from different modalities
    (e.g., RNA and protein) to facilitate comparisons of their expression patterns.

    Args:
        adata: AnnData object with multi-modal data
        mod1: Name of the first modality (e.g., 'rna', 'protein')
        mod2: Name of the second modality (e.g., 'protein', 'rna')
        features: List of tuples with paired features [(mod1_feature, mod2_feature), ...]
        basis: Dimensionality reduction to use ('umap', 'tsne', 'pca', etc.)
        ncols: Number of columns in the grid
        figsize: Figure size (width, height) in inches
        palette: Color palette for categorical variables
        cmap: Colormap for continuous variables
        size: Point size for scatter plots
        alpha: Transparency for scatter plots
        show_correlation: Whether to display correlation coefficients
        save: Path to save the figure (if not None)
        dpi: Resolution for saved figure
        show: Whether to display the plot
        **kwargs: Additional arguments passed to sc.pl.embedding

    Returns:
        Figure object for further customization

    Examples:
        >>> # Compare RNA and protein expression for paired markers
        >>> plot_multi_modality(
        ...     adata,
        ...     'rna',
        ...     'protein',
        ...     [('CD3E', 'CD3'), ('CD4', 'CD4'), ('CD8A', 'CD8')]
        ... )
        >>>
        >>> # Compare genes across modalities with correlation statistics
        >>> plot_multi_modality(
        ...     adata,
        ...     'spliced',
        ...     'unspliced',
        ...     [('IFNG', 'IFNG'), ('TNF', 'TNF'), ('IL2', 'IL2')],
        ...     basis='tsne',
        ...     show_correlation=True,
        ...     save='spliced_unspliced_comparison.png'
        ... )
    """
    log.info(f"Creating multi-modal visualization comparing {mod1} and {mod2}")

    # Validate input
    if f"X_{basis}" not in adata.obsm:
        log.error(f"Embedding '{basis}' not found in adata.obsm")
        raise ValueError(f"Embedding '{basis}' not found in adata.obsm")

    # Get embedding coordinates
    X_embed = adata.obsm[f"X_{basis}"]

    # Check if modalities exist
    if not hasattr(adata, "mod") and mod1 != "rna" and mod2 != "rna":
        log.error("Multi-modal data not found in adata.mod")
        raise ValueError("Multi-modal data not found in adata.mod")

    # Define helper to get expression for a feature from a modality
    def get_expression(modality, feature):
        if modality == "rna":
            # Check if the feature exists in the RNA data
            if feature in adata.var_names:
                return (
                    adata[:, feature].X.toarray().flatten()
                    if scipy.sparse.issparse(adata[:, feature].X)
                    else adata[:, feature].X.flatten()
                )
            else:
                log.error(f"Feature '{feature}' not found in RNA data")
                raise ValueError(f"Feature '{feature}' not found in RNA data")
        else:
            # Check if the feature exists in the specified modality
            if hasattr(adata, "mod") and modality in adata.mod:
                if feature in adata.mod[modality].var_names:
                    return (
                        adata.mod[modality][:, feature].X.toarray().flatten()
                        if scipy.sparse.issparse(adata.mod[modality][:, feature].X)
                        else adata.mod[modality][:, feature].X.flatten()
                    )
                else:
                    log.error(f"Feature '{feature}' not found in {modality} data")
                    raise ValueError(
                        f"Feature '{feature}' not found in {modality} data"
                    )
            else:
                log.error(f"Modality '{modality}' not found in adata.mod")
                raise ValueError(f"Modality '{modality}' not found in adata.mod")

    # Calculate figsize if not provided
    if figsize is None:
        n_pairs = len(features)
        n_cols = min(ncols, n_pairs)
        n_rows = int(np.ceil(n_pairs / n_cols)) * 2  # Each pair gets 2 rows
        figsize = (6 * n_cols, 5 * n_rows)

    # Create figure
    fig = plt.figure(figsize=figsize)

    # Create grid layout
    n_pairs = len(features)
    n_cols = min(ncols, n_pairs)
    n_rows = int(np.ceil(n_pairs / n_cols)) * 2  # Each pair gets 2 rows

    # Create subplots
    for i, (feat1, feat2) in enumerate(features):
        # Calculate grid positions for this pair
        row_idx = (i // n_cols) * 2
        col_idx = i % n_cols

        # Get expression values
        try:
            expr1 = get_expression(mod1, feat1)
            expr2 = get_expression(mod2, feat2)

            # Calculate correlation
            corr = np.corrcoef(expr1, expr2)[0, 1]

            # Create subplot for first modality
            ax1 = plt.subplot(n_rows, n_cols, row_idx * n_cols + col_idx + 1)
            sc1 = ax1.scatter(
                X_embed[:, 0],
                X_embed[:, 1],
                s=size,
                c=expr1,
                cmap=cmap,
                alpha=alpha,
                **kwargs,
            )
            ax1.set_title(f"{mod1.upper()}: {feat1}")
            ax1.set_xlabel(f"{basis.upper()}1")
            ax1.set_ylabel(f"{basis.upper()}2")
            plt.colorbar(sc1, ax=ax1, label=f"{mod1.upper()} Expression")

            # Create subplot for second modality
            ax2 = plt.subplot(n_rows, n_cols, (row_idx + 1) * n_cols + col_idx + 1)
            sc2 = ax2.scatter(
                X_embed[:, 0],
                X_embed[:, 1],
                s=size,
                c=expr2,
                cmap=cmap,
                alpha=alpha,
                **kwargs,
            )
            ax2.set_title(f"{mod2.upper()}: {feat2}")
            ax2.set_xlabel(f"{basis.upper()}1")
            ax2.set_ylabel(f"{basis.upper()}2")
            plt.colorbar(sc2, ax=ax2, label=f"{mod2.upper()} Expression")

            # Add correlation text if requested
            if show_correlation:
                corr_text = f"Correlation: {corr:.3f}"
                # Add correlation text between the plots
                fig.text(
                    (col_idx + 0.5) / n_cols,
                    1 - (row_idx + 1) / n_rows,
                    corr_text,
                    ha="center",
                    va="center",
                    fontsize=10,
                    fontweight="bold",
                    bbox=dict(facecolor="white", alpha=0.7, boxstyle="round"),
                )

        except Exception as e:
            log.warning(f"Error plotting feature pair ({feat1}, {feat2}): {str(e)}")
            # Create empty subplot with error message
            ax = plt.subplot(n_rows, n_cols, row_idx * n_cols + col_idx + 1)
            ax.text(
                0.5,
                0.5,
                f"Error: {str(e)}",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.axis("off")

            ax = plt.subplot(n_rows, n_cols, (row_idx + 1) * n_cols + col_idx + 1)
            ax.axis("off")

    # Add overall title
    plt.suptitle(
        f"Comparison of {mod1.upper()} and {mod2.upper()} features", fontsize=16, y=0.98
    )

    plt.tight_layout()

    # Save if requested
    if save:
        plt.savefig(save, dpi=dpi, bbox_inches="tight")
        log.info(f"Saved multi-modal visualization to {save}")

    # Show if requested
    if show:
        plt.show()

    return fig


def plot_ridge(
    adata: sc.AnnData,
    features: List[str],
    groupby: str,
    layer: Optional[str] = None,
    use_raw: bool = False,
    figsize: Optional[Tuple[float, float]] = None,
    cmap: str = "viridis",
) -> plt.Figure:
    """
    Create a ridge plot to visualize feature distributions across groups.
    """
    log.info(f"Creating ridge plot for {len(features)} features across '{groupby}'")

    # Extract data
    df_list = []
    for feature in features:
        if use_raw and adata.raw is not None:
            expr = adata.raw[:, feature].X.toarray().flatten()
        elif layer:
            expr = adata[:, feature].layers[layer].toarray().flatten()
        else:
            expr = adata[:, feature].X.toarray().flatten()

        temp_df = pd.DataFrame(
            {"expression": expr, "group": adata.obs[groupby].values, "feature": feature}
        )
        df_list.append(temp_df)

    plot_df = pd.concat(df_list)

    # Create plot
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    pal = sns.cubehelix_palette(10, rot=-0.25, light=0.7)

    g = sns.FacetGrid(
        plot_df, row="group", hue="group", aspect=10, height=0.75, palette=pal
    )
    g.map(
        sns.kdeplot,
        "expression",
        bw_adjust=0.5,
        clip_on=False,
        fill=True,
        alpha=1,
        linewidth=1.5,
    )
    g.map(sns.kdeplot, "expression", clip_on=False, color="w", lw=2, bw_adjust=0.5)
    g.map(plt.axhline, y=0, lw=2, clip_on=False)

    def label(x, color, label):
        ax = plt.gca()
        ax.text(
            0,
            0.2,
            label,
            fontweight="bold",
            color=color,
            ha="left",
            va="center",
            transform=ax.transAxes,
        )

    g.map(label, "expression")
    g.fig.subplots_adjust(hspace=-0.25)
    g.set_titles("")
    g.set(yticks=[], ylabel="")
    g.despine(bottom=True, left=True)

    plt.suptitle(f"Expression distribution across {groupby}", y=0.98)

    return g.fig


def plot_coexpression(
    adata: sc.AnnData,
    x_gene: str,
    y_gene: str,
    color_by: Optional[str] = None,
    layer: Optional[str] = None,
    use_raw: bool = False,
    figsize: Optional[Tuple[float, float]] = (7, 6),
) -> plt.Figure:
    """
    Create a scatter plot to visualize co-expression of two genes.
    """
    log.info(f"Plotting co-expression of '{x_gene}' and '{y_gene}'")

    # Extract data
    if use_raw and adata.raw is not None:
        data_source = adata.raw
    else:
        data_source = adata

    df = sc.get.obs_df(data_source, keys=[x_gene, y_gene], layer=layer)

    if color_by:
        if color_by in adata.obs.columns:
            df[color_by] = adata.obs[color_by].values
        else:
            df[color_by] = sc.get.obs_df(data_source, keys=[color_by], layer=layer)[
                color_by
            ].values

    # Create plot
    fig, ax = plt.subplots(figsize=figsize)
    sns.scatterplot(
        data=df,
        x=x_gene,
        y=y_gene,
        hue=color_by,
        s=10,
        alpha=0.7,
        ax=ax,
        edgecolor=None,
    )

    ax.set_title(f"Co-expression of {x_gene} and {y_gene}")
    plt.tight_layout()

    return fig


def plot_differential_abundance(
    diff_abundance_df: pd.DataFrame,
    pval_threshold: float = 0.05,
    figsize: tuple = (8, 6)
):
    """
    Visualize differential abundance results using a volcano-style plot.
    """
    df = diff_abundance_df.copy()
    df['-log10(pvalue)'] = -np.log10(df['pvalue'])
    df['significant'] = df['pvalue'] < pval_threshold
    
    plt.figure(figsize=figsize)
    sns.scatterplot(
        data=df,
        x='log2fc_abundance',
        y='-log10(pvalue)',
        hue='significant',
        palette={True: 'red', False: 'grey'},
        size='mean_abundance_group1',
        sizes=(50, 500),
        alpha=0.7
    )
    
    # Add labels for significant points
    for i, row in df[df['significant']].iterrows():
        plt.text(row['log2fc_abundance'], row['-log10(pvalue)'], row['cell_type'], fontsize=9)
        
    plt.title("Differential Cell Type Abundance")
    plt.xlabel("Log2 Fold Change in Abundance")
    plt.ylabel("-log10(p-value)")
    plt.axhline(-np.log10(pval_threshold), ls='--', color='grey')
    plt.axvline(0, ls='--', color='grey')
    plt.legend(title='Significant')
    plt.tight_layout()
    plt.show()