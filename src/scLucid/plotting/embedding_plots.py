"""
Plotting functions for single-cell RNA-seq data.
"""

import logging
from typing import Dict, List, Literal, Optional, Tuple, Union

import matplotlib.patheffects as PathEffects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from .plotting_utils import _get_palette_map, _subset_adata

# Try importing adjustText softly
try:
    from adjustText import adjust_text
except ImportError:
    adjust_text = None

log = logging.getLogger(__name__)


def plot_embedding(
    adata: sc.AnnData,
    color_by: Union[str, List[str]],
    basis: str = "umap",
    subset: Optional[pd.Series] = None,
    title: Optional[str] = None,
    show_labels: bool = True,
    palette: Optional[Union[str, Dict[str, str]]] = None,
    size: float = 12,
    alpha: float = 0.8,
    ncols: int = 3,
    figsize: Optional[Tuple[float, float]] = None,
    legend_loc: Literal[
        "right margin", "left margin", "top margin", "bottom margin", "none"
    ] = "right margin",
    label_size: int = 10,
    save: Optional[str] = None,
    dpi: int = 300,
    ax: Optional[plt.Axes] = None,
    show: bool = True,
    legend_style: Literal["on_data", "legend", "both", "none"] = "on_data",
    rasterized: bool = True,
    **kwargs,
) -> plt.Figure:
    """
    Enhanced wrapper for scanpy's embedding plot with subsetting.

    This function provides advanced features like on-data labeling, automatic
    grid plotting for multiple features, and robust color mapping. It now also
    supports on-the-fly subsetting of cells.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix.
    color_by : Union[str, List[str]]
        The key in `adata.obs` or a gene in `adata.var_names` to color points by.
        If a list, creates a grid of plots.
    basis : str, default 'umap'
        The embedding to use (e.g., 'umap', 'tsne').
    subset : Optional[pd.Series], default None
        A boolean Series to subset the data before plotting.
        Example: `adata.obs['cell_type'].isin(['B_cell', 'T_cell'])`.
    title : Optional[str], default None
        Title for the plot.
    show_labels : bool, default True
        If True and `color_by` is categorical, plot labels on the data centroids.
    palette : Optional[Union[str, Dict[str, str]]], default None
        Colors to use for categorical data.
    size : float, default 12
        Point size.
    alpha : float, default 0.8
        Point transparency.
    ncols : int, default 3
        Number of columns for the grid if `color_by` is a list.
    figsize : Optional[Tuple[float, float]], default None
        Figure size.
    legend_loc : str, default 'right margin'
        Location of the legend.
    label_size : int, default 10
        Font size for on-data labels.
    save : Optional[str], default None
        Path to save the figure.
    dpi : int, default 300
        DPI for saving the figure.
    ax : Optional[plt.Axes], default None
        A matplotlib axes object to plot on.
    show : bool, default True
        Whether to display the plot.
    legend_style : Literal["on_data", "legend", "both", "none"], default 'on_data'
        How to display the legend.
    rasterized : bool, default True
        Rasterize the scatter points for smaller file sizes.
    **kwargs
        Additional keyword arguments passed to `sc.pl.embedding`.

    Returns:
    -------
    plt.Figure
        The matplotlib Figure object.
    """
    # --- 1. Handle Subsetting ---
    adata_to_plot = _subset_adata(adata, subset)

    # --- 2. Recursive call for list handling ---
    if isinstance(color_by, list):
        if len(color_by) == 1:
            return plot_embedding(
                adata_to_plot,
                color_by=color_by[0],
                basis=basis,
                title=title,
                show_labels=show_labels,
                palette=palette,
                size=size,
                alpha=alpha,
                ncols=ncols,
                figsize=figsize,
                legend_loc=legend_loc,
                label_size=label_size,
                save=save,
                dpi=dpi,
                ax=ax,
                show=show,
                legend_style=legend_style,
                rasterized=rasterized,
                subset=None,
                **kwargs,
            )
        n_plots = len(color_by)
        n_rows = int(np.ceil(n_plots / ncols))
        if figsize is None:
            figsize = (4 * ncols, 3.5 * n_rows)
        fig, axes = plt.subplots(n_rows, ncols, figsize=figsize, squeeze=False)
        axes = axes.flatten()
        for i, color in enumerate(color_by):
            plot_embedding(
                adata_to_plot,
                color,
                basis,
                ax=axes[i],
                show=False,
                save=None,
                show_labels=show_labels,
                palette=palette,
                size=size,
                alpha=alpha,
                subset=None,
                **kwargs,
            )
        for i in range(n_plots, len(axes)):
            axes[i].axis("off")
        plt.tight_layout()
        if save:
            plt.savefig(save, dpi=dpi, bbox_inches="tight")
        if show:
            plt.show()
        return fig

    # --- 3. Standard Plotting Logic ---
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize or (6, 5))
    else:
        fig = ax.figure

    # Use adata_to_plot for all subsequent operations
    sc.pl.embedding(
        adata_to_plot,
        basis=basis,
        color=color_by,
        ax=ax,
        show=False,
        size=size,
        alpha=alpha,
        palette=palette,
        legend_loc="none" if legend_style == "on_data" else legend_loc,
        **kwargs,
    )

    if rasterized:
        [collection.set_rasterized(True) for collection in ax.collections]

    # --- 4. Label Logic (Only for categorical data) ---
    is_categorical = (
        pd.api.types.is_categorical_dtype(adata_to_plot.obs[color_by])
        if color_by in adata_to_plot.obs
        else False
    )
    if is_categorical and show_labels and legend_style in ["on_data", "both"]:
        color_map = _get_palette_map(adata_to_plot, color_by, palette)
        categories = adata_to_plot.obs[color_by].cat.categories
        embed_key = f"X_{basis}"
        texts = []
        for label in categories:
            mask = adata_to_plot.obs[color_by] == label
            if not np.any(mask):
                continue
            coords = adata_to_plot.obsm[embed_key][mask]
            x, y = np.median(coords, axis=0)
            bg_color = color_map.get(label, "white")
            txt = ax.text(
                x,
                y,
                label,
                fontsize=label_size,
                fontweight="bold",
                ha="center",
                va="center",
                color="black",
            )
            txt.set_path_effects(
                [PathEffects.withStroke(linewidth=3, foreground=bg_color, alpha=0.8)]
            )
            texts.append(txt)
        if adjust_text and texts:
            try:
                adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))
            except Exception as e:
                log.debug(f"adjust_text failed: {e}")

    ax.set_title(title if title else f"{basis.upper()}: {color_by}")
    if ax.get_legend() and legend_style == "on_data":
        ax.get_legend().remove()
    if save:
        plt.savefig(save, dpi=dpi, bbox_inches="tight")
        log.info(f"Saved embedding plot to {save}")
    if show:
        plt.show()

    return fig


def plot_faceted_embedding(
    adata: sc.AnnData,
    split_by: str,
    color_by: str,
    basis: str = "umap",
    subset: Optional[pd.Series] = None,
    col_wrap: int = 4,
    point_size: float = 10,
    alpha: float = 0.8,
    palette: Optional[Union[str, Dict[str, str]]] = None,
    figsize: Tuple[float, float] = (12, 10),
    main_title: Optional[str] = None,
    frameon: bool = True,
    show_grid: bool = False,
    show_ticks: bool = False,
    aspect: Optional[float] = None,
    legend_loc: Literal[
        "right margin", "left margin", "top margin", "bottom margin", "none"
    ] = "right margin",
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """
    Creates faceted scatter plots of embeddings, colored by a categorical variable.

    This function generates a grid of plots (facets), where each plot corresponds
    to a unique category in the `split_by` column of `adata.obs`. Within each
    facet, cells are plotted on the specified embedding (e.g., UMAP, t-SNE)
    and colored according to their category in the `color_by` column.

    This function is designed for coloring by categorical data. For coloring
    by gene expression or other continuous features, use `plot_faceted_feature`.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix.
    split_by : str
        The column in `adata.obs` that defines the facets (e.g., 'sample', 'batch').
    color_by : str
        The categorical column in `adata.obs` to use for coloring points (e.g., 'cell_type').
    basis : str, default 'umap'
        The embedding to use. The function looks for `X_{basis}` in `adata.obsm`.
    col_wrap : int, default 4
        The maximum number of facets per row.
    point_size : float, default 10
        Size of the scatter plot points.
    alpha : float, default 0.8
        Point transparency.
    palette : Optional[Union[str, Dict[str, str]]], default None
        Colors to use for the `color_by` categories. Can be a named seaborn
        palette (e.g., 'tab20'), a matplotlib colormap, or a dictionary mapping
        categories to colors. If None, it will try to use colors stored in
        `adata.uns[f'{color_by}_colors']` or generate a default palette.
    figsize : Tuple[float, float], default (12, 10)
        Width and height of the figure in inches.
    main_title : Optional[str], default None
        A main title for the entire figure.
    frameon : bool, default True
        Whether to draw a frame around each subplot.
    show_grid : bool, default False
        Whether to draw a grid on each subplot.
    show_ticks : bool, default False
        Whether to show axis ticks and labels on each subplot.
    aspect : Optional[float], default None
        The aspect ratio of the subplots (e.g., 1 for a square plot). If None,
        the aspect ratio is determined automatically.
    legend_loc : Literal["right margin", "top margin", "none"], default "right margin"
        Position of the legend. "right margin" places it to the right of the
        grid, "top margin" places it above, and "none" hides it.
    save : Optional[str], default None
        Path to save the figure. File extension determines the format (e.g., '.png', '.pdf').
    dpi : int, default 300
        The resolution in dots per inch for the saved figure.
    show : bool, default True
        Whether to display the plot.
    **kwargs
        Additional keyword arguments passed to `seaborn.scatterplot`.

    Returns:
    -------
    plt.Figure
        The matplotlib Figure object containing the plot.
    """
    adata_to_plot = _subset_adata(adata, subset)

    # --- 1. Input Validation ---
    embed_key = f"X_{basis}" if not basis.startswith("X_") else basis
    if embed_key not in adata_to_plot.obsm:
        raise ValueError(f"Embedding '{embed_key}' not found in adata.obsm.")
    if split_by not in adata_to_plot.obs:
        raise ValueError(f"'{split_by}' not found in adata.obs.")
    is_categorical = color_by in adata_to_plot.obs and pd.api.types.is_categorical_dtype(
        adata_to_plot.obs[color_by]
    )
    if not is_categorical:
        raise ValueError(f"'{color_by}' not found as a categorical column in adata.obs.")

    # --- 2. Data Preparation ---
    coords = adata_to_plot.obsm[embed_key][:, :2]
    plot_df = pd.DataFrame(coords, columns=["Dim1", "Dim2"])
    if isinstance(adata_to_plot.obs[split_by].dtype, pd.CategoricalDtype):
        categories_order = adata_to_plot.obs[split_by].cat.categories.tolist()
    else:
        categories_order = sorted(adata_to_plot.obs[split_by].unique())
    plot_df[split_by] = pd.Categorical(
        adata_to_plot.obs[split_by].values, categories=categories_order, ordered=True
    )
    plot_df[color_by] = adata_to_plot.obs[color_by].values

    # --- 3. Palette Handling ---
    final_palette = (
        palette if palette is not None else _get_palette_map(adata_to_plot, color_by, None)
    )

    # --- 4. Plotting with FacetGrid ---
    g = sns.FacetGrid(plot_df, col=split_by, col_wrap=col_wrap, sharex=True, sharey=True)
    g.map_dataframe(
        sns.scatterplot,
        x="Dim1",
        y="Dim2",
        hue=color_by,
        palette=final_palette,
        s=point_size,
        alpha=alpha,
        linewidth=0,
        **kwargs,
    )

    # --- 5. Styling and Annotations ---
    g.set_axis_labels(f"{basis}_1", f"{basis}_2")
    g.set_titles("{col_name}")

    if legend_loc == "right margin":
        g.add_legend(title=color_by, bbox_to_anchor=(1.01, 0.5), loc="center left")
    elif legend_loc == "top margin":
        g.add_legend(title=color_by, bbox_to_anchor=(0.5, 1.02), loc="lower center", ncol=3)

    # Main title
    if main_title:
        g.fig.suptitle(
            main_title,
            y=1.02 + (0.02 if legend_loc == "top margin" else 0),
            fontsize=16,
        )

    # Axis and frame styling
    for ax in g.axes.flat:
        ax.grid(show_grid)
        ax.set_frame_on(frameon)
        if not show_ticks:
            ax.set_xticks([])
            ax.set_yticks([])
        # Remove individual legends if a global one is added
        if legend_loc != "none" and ax.get_legend():
            ax.get_legend().remove()

    if aspect:
        for ax in g.axes.flat:
            ax.set_aspect(aspect)

    # --- 6. Final Adjustments and Output ---
    # Adjust layout to make space for title and legend
    g.fig.set_size_inches(figsize)
    g.fig.tight_layout(
        rect=[
            0,
            0,
            0.95 if legend_loc == "right margin" else 1,
            0.96 if main_title else 1,
        ]
    )

    if save:
        g.fig.savefig(save, dpi=dpi, bbox_inches="tight")
        log.info(f"Saved faceted embedding plot to {save}")

    if show:
        plt.show()

    return g.fig


# =============================================================================
# Marker Visualization
# =============================================================================
