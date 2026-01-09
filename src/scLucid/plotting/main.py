"""
Dimensionality reduction and visualization functions for single-cell RNA-seq data.
"""

import logging
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.patheffects as PathEffects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.spatial import distance

# Try importing adjustText softly
try:
    from adjustText import adjust_text
except ImportError:
    adjust_text = None

# Configure logging
log = logging.getLogger(__name__)

__all__ = [
    "plot_embedding",
    "plot_faceted_embedding",
    "plot_dotplot",
    "plot_stacked_violin",
    "plot_split_violin_with_stats",
    "plot_marker_expression",
    "plot_faceted_feature",
    "plot_marker_heatmap",
    "plot_volcano",
    "plot_ridge",
    "plot_feature_correlation",
    "plot_coexpression",
    "plot_differential_abundance",
]


# =============================================================================
# Helper Functions
# =============================================================================


def _subset_adata(adata: sc.AnnData, subset: Optional[pd.Series]) -> sc.AnnData:
    """Helper to safely subset AnnData."""
    if subset is not None:
        if (
            not isinstance(subset, pd.Series)
            or not subset.dtype == bool
            or len(subset) != adata.n_obs
        ):
            raise ValueError(
                "`subset` must be a boolean pandas Series of length adata.n_obs."
            )
        log.info(f"Subsetting data to {subset.sum()} cells.")
        return adata[subset].copy()
    return adata


def _combine_groupby(adata: sc.AnnData, groupby_main: str, groupby_sub: str) -> str:
    """Helper to create a combined groupby column and return its name."""
    if groupby_main not in adata.obs or groupby_sub not in adata.obs:
        raise ValueError("groupby_main and groupby_sub must be columns in adata.obs.")

    combined_col_name = f"{groupby_main}_{groupby_sub}"
    adata.obs[combined_col_name] = (
        adata.obs[groupby_main].astype(str) + "_" + adata.obs[groupby_sub].astype(str)
    )

    main_cats = adata.obs[groupby_main].cat.categories
    sub_cats = adata.obs[groupby_sub].cat.categories
    combined_order = [f"{m}_{s}" for m in main_cats for s in sub_cats]

    adata.obs[combined_col_name] = pd.Categorical(
        adata.obs[combined_col_name], categories=combined_order, ordered=True
    )
    return combined_col_name


def _get_palette_map(
    adata: sc.AnnData, key: str, palette: Optional[Union[str, Dict]] = None
) -> Dict[str, Any]:
    """
    Robustly resolve color mapping for a categorical column.
    Priority: User Dict > adata.uns > Scanpy default generation.
    """
    if not pd.api.types.is_categorical_dtype(adata.obs[key]):
        return {}

    categories = adata.obs[key].cat.categories

    # 1. User provided dictionary
    if isinstance(palette, dict):
        # Fill missing keys with gray
        return {cat: palette.get(cat, "#cccccc") for cat in categories}

    # 2. Check adata.uns (Scanpy convention: key_colors)
    uns_key = f"{key}_colors"
    if uns_key in adata.uns:
        colors = adata.uns[uns_key]
        if len(colors) >= len(categories):
            return dict(zip(categories, colors[: len(categories)]))

    # 3. Generate new palette (Seaborn/Matplotlib)
    # If user provided a string palette name (e.g., 'tab20'), use it
    palette_name = palette if isinstance(palette, str) else "tab20"
    if len(categories) > 20 and palette_name == "tab20":
        palette_name = "husl"  # Fallback for many categories

    generated_colors = sns.color_palette(palette_name, n_colors=len(categories))
    return dict(zip(categories, generated_colors))


def _sort_genes_within_categories(
    agg_df: pd.DataFrame, marker_dict: Dict[str, List[str]]
) -> List[str]:
    """
    Cluster genes within categories based on their expression across groups.

    Parameters
    ----------
    agg_df : pd.DataFrame
        Aggregated expression data (Groups × Genes)
    marker_dict : Dict[str, List[str]]
        Gene categories

    Returns
    -------
    List[str]
        Sorted gene names
    """
    sorted_genes = []

    for category, genes in marker_dict.items():
        valid_genes = [g for g in genes if g in agg_df.columns]

        if len(valid_genes) < 3:
            sorted_genes.extend(valid_genes)
            continue

        # Subset dataframe (genes are columns here)
        sub_df = agg_df[valid_genes]
        sub_df_T = sub_df.T  # genes x cells

        sub_df_T = sub_df_T.loc[sub_df_T.var(axis=1) > 0]
        if len(sub_df_T) < 2:
            sorted_genes.extend(valid_genes)
            continue
        try:
            Z = hierarchy.linkage(
                distance.pdist(sub_df_T.values, metric="correlation"), method="average"
            )
            leaves = hierarchy.dendrogram(Z, no_plot=True)["leaves"]
            sorted_genes.extend(sub_df_T.index[leaves].tolist())
        except Exception as e:
            log.warning(f"Clustering failed for '{category}': {e}")
            sorted_genes.extend(valid_genes)

    return sorted_genes


# =============================================================================
# Embedding Visualization Functions
# =============================================================================


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

    Returns
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
                legend_style,
                rasterized,
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
                adjust_text(
                    texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5)
                )
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

    Returns
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
    is_categorical = (
        color_by in adata_to_plot.obs
        and pd.api.types.is_categorical_dtype(adata_to_plot.obs[color_by])
    )
    if not is_categorical:
        raise ValueError(
            f"'{color_by}' not found as a categorical column in adata.obs."
        )

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
        palette
        if palette is not None
        else _get_palette_map(adata_to_plot, color_by, None)
    )

    # --- 4. Plotting with FacetGrid ---
    g = sns.FacetGrid(
        plot_df, col=split_by, col_wrap=col_wrap, sharex=True, sharey=True
    )
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
        g.add_legend(
            title=color_by, bbox_to_anchor=(0.5, 1.02), loc="lower center", ncol=3
        )

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


def plot_dotplot(
    adata: sc.AnnData,
    var_names: Union[str, List[str], Dict[str, List[str]]],
    groupby: Optional[str] = None,
    groupby_main: Optional[str] = None,
    groupby_sub: Optional[str] = None,
    subset: Optional[pd.Series] = None,
    auto_order_categories: bool = True,
    use_raw: bool = None,
    layer: Optional[str] = None,
    standard_scale: Optional[Literal["var", "group"]] = None,
    cmap: Literal["viridis", "Reds", "RdBu_r"] = "RdBu_r",
    dot_min: float = 0.0,
    dot_max: Optional[float] = None,
    title: Optional[str] = None,
    swap_axes: bool = False,
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
    **kwargs,
) -> sc.pl.DotPlot:
    """
    Enhanced dotplot function with combined grouping and subsetting capabilities.

    This function is a powerful wrapper around `sc.pl.dotplot`, adding features
    for combined grouping (e.g., celltype_condition), on-the-fly data subsetting,
    and intelligent category ordering.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix.
    var_names : Union[str, List[str], Dict[str, List[str]]]
        Genes or variables to plot. Can be a single gene, a list of genes, or a
        dictionary of gene lists for grouped plotting.
    groupby : Optional[str], default None
        The key in `adata.obs` to group the data by. Ignored if `groupby_main`
        and `groupby_sub` are provided.
    groupby_main, groupby_sub : Optional[str], default None
        Two keys in `adata.obs` to create a combined grouping. For example,
        `groupby_main='celltype'` and `groupby_sub='condition'` will create
        groups like 'Fibroblast_Control', 'Fibroblast_Stimulated', etc.
    subset : Optional[pd.Series], default None
        A boolean Series (e.g., `adata.obs['celltype'].isin(['B', 'T'])`) to
        subset the data before plotting. This is more efficient than subsetting
        `adata` beforehand.
    auto_order_categories : bool, default True
        If True, and `categories_order` is not provided in `kwargs`, it will
        automatically use the categorical order from `adata.obs[groupby]`.
    use_raw : bool, default None
        Use `adata.raw` for plotting. If None, uses `adata.raw` if it exists.
    layer : Optional[str], default None
        The layer in `adata.layers` to use for plotting.
    standard_scale : Optional[Literal["var", "group"]], default None
        Whether to standardize the data (z-score). 'var' standardizes across
        genes (columns), 'group' across groups (rows).
    cmap : str, default 'Reds'
        Colormap for the dot colors.
    dot_min, dot_max : float, Optional[float]
        Size range for the dots.
    title : Optional[str], default None
        Main title for the plot.
    swap_axes : bool, default False
        If True, swaps the x and y axes.
    save : Optional[str], default None
        Path to save the figure.
    dpi : int, default 300
        Resolution for the saved figure.
    show : bool, default True
        Whether to display the plot.
    **kwargs
        Additional keyword arguments passed to `sc.pl.dotplot`.

    Returns
    -------
    sc.pl.DotPlot
        The DotPlot object from scanpy, which can be used for further
        customization (e.g., `dp.add_totals()`).
    """
    # --- 1. Input Validation ---
    if not (groupby or (groupby_main and groupby_sub)):
        raise ValueError(
            "Provide either 'groupby' or both 'groupby_main' and 'groupby_sub'."
        )

    # --- 2. Handle Subsetting and Grouping using Helper Functions ---
    adata_to_plot = _subset_adata(adata, subset)

    if groupby_main and groupby_sub:
        groupby = _combine_groupby(adata_to_plot, groupby_main, groupby_sub)

    # --- 3. Handle Category Ordering ---
    if auto_order_categories and "categories_order" not in kwargs:
        if pd.api.types.is_categorical_dtype(adata_to_plot.obs[groupby]):
            kwargs["categories_order"] = adata_to_plot.obs[
                groupby
            ].cat.categories.tolist()

    # --- 4. Call sc.pl.dotplot with show=False to get the object ---
    dp = sc.pl.dotplot(
        adata_to_plot,
        var_names=var_names,
        groupby=groupby,
        use_raw=use_raw,
        layer=layer,
        standard_scale=standard_scale,
        cmap=cmap,
        dot_min=dot_min,
        dot_max=dot_max,
        title=title,
        swap_axes=swap_axes,
        show=False,  
        **kwargs,
    )

    # --- 5. Save, Show, and Return ---
    fig = dp["mainplot_ax"].figure

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
        log.info(f"Saved dotplot to {save}")

    if show:
        plt.show()

    return dp


def plot_stacked_violin(
    adata: sc.AnnData,
    var_names: Union[str, List[str]],
    groupby: Optional[str] = None,
    groupby_main: Optional[str] = None,
    groupby_sub: Optional[str] = None,
    subset: Optional[pd.Series] = None,
    use_raw: bool = None,
    layer: Optional[str] = None,
    title: Optional[str] = None,
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """
    Creates a stacked violin plot to visualize gene expression across groups.

    This function is a wrapper around `sc.pl.stacked_violin`, enhanced with
    subsetting and combined grouping capabilities.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix.
    var_names : Union[str, List[str]]
        Genes or variables to plot.
    groupby : Optional[str], default None
        The key in `adata.obs` to group the data by. Ignored if `groupby_main`
        and `groupby_sub` are provided.
    groupby_main, groupby_sub : Optional[str], default None
        Two keys in `adata.obs` to create a combined grouping. For example,
        `groupby_main='celltype'` and `groupby_sub='condition'` will create
        groups like 'Fibroblast_Control', 'Fibroblast_Stimulated', etc.
    subset : Optional[pd.Series], default None
        A boolean Series to subset the data before plotting.
    use_raw : bool, default None
        Use `adata.raw` for plotting. If None, uses `adata.raw` if it exists.
    layer : Optional[str], default None
        The layer in `adata.layers` to use for plotting.
    title : Optional[str], default None
        Title for the plot.
    save : Optional[str], default None
        Path to save the figure.
    dpi : int, default 300
        Resolution for the saved figure.
    show : bool, default True
        Whether to display the plot.
    **kwargs
        Additional keyword arguments passed to `sc.pl.stacked_violin`.

    Returns
    -------
    plt.Figure
        The matplotlib Figure object.

    Examples
    --------
    >>> # Standard stacked violin plot
    >>> plot_stacked_violin(adata, var_names=['LYZ', 'S100A8'], groupby='cell_type')

    >>> # Stacked violin plot with combined grouping
    >>> plot_stacked_violin(adata, var_names=['CD3D', 'MS4A1'], groupby_main='cell_type', groupby_sub='condition')
    """
    # --- 1. Input Validation and Data Preparation ---
    if not (groupby or (groupby_main and groupby_sub)):
        raise ValueError(
            "Provide either 'groupby' or both 'groupby_main' and 'groupby_sub'."
        )

    adata_to_plot = _subset_adata(adata, subset)

    if groupby_main and groupby_sub:
        groupby = _combine_groupby(adata_to_plot, groupby_main, groupby_sub)

    # --- 2. Call sc.pl.stacked_violin with show=False to get the object ---
    sc.pl.stacked_violin(
        adata_to_plot,
        var_names=var_names,
        groupby=groupby,
        use_raw=use_raw,
        layer=layer,
        title=title,
        show=False, 
        **kwargs,
    )

    # --- 3. Save, Show, and Return ---
    fig = plt.gcf()

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
        log.info(f"Saved stacked violin plot to {save}")

    if show:
        plt.show()

    return fig


def plot_split_violin_with_stats(
    adata: sc.AnnData,
    genes: Union[str, List[str]],
    celltype_col: str,
    condition_col: str,
    subset: Optional[pd.Series] = None,
    test: Literal["Mann-Whitney", "t-test"] = "Mann-Whitney",
    palette: str = "viridis",
    use_raw: bool = True,
    layer: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """
    Creates a split violin plot with statistical annotations.

    This function visualizes the expression of genes across different cell types,
    split by a condition (e.g., control vs. stimulated), and automatically
    adds statistical significance tests.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix.
    genes : Union[str, List[str]]
        Gene(s) to plot.
    celltype_col : str
        The column in `adata.obs` that defines the cell types (x-axis).
    condition_col : str
        The column in `adata.obs` that defines the condition for splitting (color).
    subset : Optional[pd.Series], default None
        A boolean Series to subset the data before plotting.
    test : Literal["Mann-Whitney", "t-test"], default 'Mann-Whitney'
        The statistical test to use.
    palette : str, default 'viridis'
        Colormap for the violins.
    use_raw : bool, default True
        Whether to use `adata.raw` for gene expression data.
    layer : Optional[str], default None
        The layer in `adata.layers` to use for plotting.
    figsize : Optional[Tuple[float, float]], default None
        Figure size. If None, it's calculated automatically.
    save : Optional[str], default None
        Path to save the figure.
    dpi : int, default 300
        Resolution for the saved figure.
    show : bool, default True
        Whether to display the plot.
    **kwargs
        Additional keyword arguments passed to `seaborn.violinplot`.

    Returns
    -------
    plt.Figure
        The matplotlib Figure object.

    Raises
    ------
    ImportError
        If the `statannotations` library is not installed.
    """
    try:
        from statannotations.Annotator import Annotator
    except ImportError:
        log.error("statannotations not installed.")
        raise ImportError("Please install statannotations: pip install statannotations")

    # --- 1. Handle Subsetting and Data Preparation ---
    adata_to_plot = _subset_adata(adata, subset)

    if isinstance(genes, str):
        genes = [genes]

    # Use sc.get.obs_df which handles raw/layer logic robustly
    try:
        df = sc.get.obs_df(
            adata_to_plot,
            keys=[celltype_col, condition_col] + genes,
            use_raw=use_raw,
            layer=layer,
        )
    except KeyError as e:
        raise ValueError(f"Could not find data for plotting: {e}")

    # --- 2. Plotting ---
    if figsize is None:
        figsize = (4 * len(genes), 5)

    fig, axes = plt.subplots(1, len(genes), figsize=figsize, squeeze=False)
    axes = axes.flatten()

    groups = df[condition_col].unique()
    if len(groups) != 2:
        raise ValueError(
            f"'{condition_col}' must have exactly 2 groups, found {len(groups)}"
        )

    for i, gene in enumerate(genes):
        ax = axes[i]
        sns.violinplot(
            data=df,
            x=celltype_col,
            y=gene,
            hue=condition_col,
            split=True,
            inner="quartile",
            palette=palette,
            ax=ax,
            cut=0,
            **kwargs,
        )
        ax.set_title(gene)
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=45)

        # Stats
        box_pairs = [
            ((ct, groups[0]), (ct, groups[1]))
            for ct in df[celltype_col].unique()
        ]
        try:
            annot = Annotator(
                ax, box_pairs, data=df, x=celltype_col, y=gene, hue=condition_col
            )
            annot.configure(test=test, text_format="star", loc="inside", verbose=0)
            annot.apply_and_annotate()
        except Exception as e:
            log.warning(f"Stats failed for {gene}: {e}")

        # Legend cleanup
        if i < len(genes) - 1 and ax.get_legend():
            ax.get_legend().remove()

    plt.tight_layout()

    # --- 3. Save, Show, and Return ---
    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
        log.info(f"Saved split violin plot to {save}")

    if show:
        plt.show()

    return fig


def plot_marker_expression(
    adata: sc.AnnData,
    markers: Union[str, List[str]],
    subset: Optional[pd.Series] = None,
    basis: str = "umap",
    ncols: int = 4,
    use_raw: bool = False,
    layer: Optional[str] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    cmap: str = "viridis",
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """
    Creates a grid of embedding plots, each colored by the expression of a single gene.

    This is a convenient wrapper for visualizing the spatial distribution of multiple
    marker genes on an embedding like UMAP or t-SNE.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix.
    markers : Union[str, List[str]]
        Gene(s) to plot.
    subset : Optional[pd.Series], default None
        A boolean Series to subset the data before plotting.
    basis : str, default 'umap'
        The embedding to use (e.g., 'umap', 'tsne').
    ncols : int, default 4
        Number of plots per row.
    use_raw : bool, default False
        Whether to use `adata.raw` for gene expression data.
    layer : Optional[str], default None
        The layer in `adata.layers` to use for gene expression data.
    vmin, vmax : Optional[float], default None
        The data range that the colormap covers. If None, it's determined
        automatically for each plot. Setting them manually ensures a consistent
        color scale across all genes for better comparison.
    cmap : str, default 'viridis'
        The colormap to use for all plots.
    save : Optional[str], default None
        Path to save the figure.
    dpi : int, default 300
        Resolution for the saved figure.
    show : bool, default True
        Whether to display the plot.
    **kwargs
        Additional keyword arguments passed to `sc.pl.embedding`.

    Returns
    -------
    plt.Figure
        The matplotlib Figure object.

    Examples
    --------
    >>> # Plot expression of two genes
    >>> plot_marker_expression(adata, markers=['CD3D', 'MS4A1'])
    """
    # --- 1. Handle Subsetting and Input Validation ---
    adata_to_plot = _subset_adata(adata, subset)

    if isinstance(markers, str):
        markers = [markers]

    # Clean markers
    valid_markers = []
    source = adata_to_plot.raw if (use_raw and adata_to_plot.raw) else adata_to_plot
    for m in markers:
        if m in source.var_names:
            valid_markers.append(m)
        else:
            log.warning(f"Marker {m} not found.")

    if not valid_markers:
        raise ValueError("No valid markers.")

    # --- 2. Plotting ---
    n_plots = len(valid_markers)
    n_rows = int(np.ceil(n_plots / ncols))

    fig, axes = plt.subplots(
        n_rows, ncols, figsize=(4 * ncols, 3.5 * n_rows), squeeze=False
    )
    axes = axes.flatten()

    for i, m in enumerate(valid_markers):
        sc.pl.embedding(
            adata_to_plot,
            basis=basis,
            color=m,
            ax=axes[i],
            show=False,
            use_raw=use_raw,
            layer=layer,
            legend_loc="none",
            vmin=vmin, # Apply consistent color scale
            vmax=vmax, # Apply consistent color scale
            cmap=cmap,   # Apply consistent colormap
            **kwargs,
        )
        axes[i].set_title(m)

    for i in range(n_plots, len(axes)):
        axes[i].axis("off")

    plt.tight_layout()

    # --- 3. Save, Show, and Return ---
    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
        log.info(f"Saved marker expression plot to {save}")

    if show:
        plt.show()

    return fig


def plot_faceted_feature(
    adata: sc.AnnData,
    feature: str,
    split_by: str,
    basis: str = "umap",
    col_wrap: int = 2,
    cmap: Literal["viridis", "Reds", "RdBu_r"] = "RdBu_r",
    point_size: float = 2,
    alpha: float = 0.8,
    use_raw: bool = True,
    layer: Optional[str] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    percentile_range: Tuple[float, float] = (1, 99),
    panel_order: Optional[List[str]] = None,
    figsize: Optional[Tuple[float, float]] = None,
    main_title: Optional[str] = None,
    colorbar_loc: Literal["right", "bottom"] = "right",
    frameon: bool = True,
    show_ticks: bool = False,
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """
    Creates faceted plots of a continuous feature (e.g., gene expression).

    This function generates a grid of plots, where each plot shows the expression
    of a single feature (like a gene) for a specific category in `split_by`.
    It uses `gridspec` for robust layout and a shared colorbar.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix.
    feature : str
        The name of the feature (gene) to plot.
    split_by : str
        The column in `adata.obs` that defines the facets.
    basis : str, default 'umap'
        The embedding to use. The function looks for `X_{basis}` in `adata.obsm`.
    col_wrap : int, default 4
        The maximum number of facets per row.
    cmap : str, default 'viridis'
        The colormap to use for the feature.
    point_size : float, default 2
        Size of the scatter plot points.
    alpha : float, default 0.8
        Point transparency.
    use_raw : bool, default True
        Whether to use `adata.raw` for feature extraction.
    layer : Optional[str], default None
        The layer in `adata.layers` to use for feature extraction.
    vmin, vmax : Optional[float], default None
        The data range that the colormap covers. If None, it's determined
        automatically from `percentile_range`.
    percentile_range : Tuple[float, float], default (1, 99)
        The lower and upper percentiles used to determine `vmin` and `vmax`
        if they are not explicitly provided.
    panel_order : Optional[List[str]], default None
        The order in which to display the facets. If None, uses the default
        categorical order of `split_by`.
    figsize : Optional[Tuple[float, float]], default None
        The figure size. If None, it's calculated automatically.
    main_title : Optional[str], default None
        A main title for the entire figure.
    colorbar_loc : {'right', 'bottom'}, default 'right'
        Location of the shared colorbar.
    frameon : bool, default True
        Whether to draw a frame around each subplot.
    show_ticks : bool, default False
        Whether to show axis ticks and labels.
    save : Optional[str], default None
        Path to save the figure.
    dpi : int, default 300
        The resolution in dots per inch for the saved figure.
    show : bool, default True
        Whether to display the plot.
    **kwargs
        Additional keyword arguments passed to `ax.scatter`.

    Returns
    -------
    plt.Figure
        The matplotlib Figure object.
    """
    # --- 1. Input Validation and Data Preparation ---
    embed_key = f"X_{basis}" if not basis.startswith("X_") else basis
    if embed_key not in adata.obsm:
        raise ValueError(f"Embedding '{embed_key}' not found in adata.obsm.")
    if split_by not in adata.obs:
        raise ValueError(f"'{split_by}' not found in adata.obs.")

    # Determine panel order
    if panel_order is None:
        if pd.api.types.is_categorical_dtype(adata.obs[split_by]):
            groups = adata.obs[split_by].cat.categories.tolist()
        else:
            groups = sorted(adata.obs[split_by].unique())
    else:
        groups = panel_order

    n_groups = len(groups)
    if n_groups == 0:
        raise ValueError(f"No valid groups found in '{split_by}'.")

    log.info(
        f"Plotting feature '{feature}' for {n_groups} groups: {', '.join(map(str, groups))}"
    )

    # Get coordinates
    coords = adata.obsm[embed_key][:, :2]
    plot_df = pd.DataFrame(coords, columns=["Dim1", "Dim2"])
    plot_df["group"] = pd.Categorical(
        adata.obs[split_by], categories=groups, ordered=True
    )

    # Get feature expression
    try:
        expr = sc.get.obs_df(adata, keys=[feature], use_raw=use_raw, layer=layer)[
            feature
        ].values
    except KeyError:
        raise ValueError(
            f"Feature '{feature}' not found. Check `use_raw` and `layer` parameters."
        )
    plot_df["expression"] = expr

    # Determine color range
    if vmin is None or vmax is None:
        v_min, v_max = np.percentile(plot_df["expression"], percentile_range)
        vmin = v_min if vmin is None else vmin
        vmax = v_max if vmax is None else vmax

    # --- 2. Layout with GridSpec ---
    ncols = min(n_groups, col_wrap)
    nrows = int(np.ceil(n_groups / ncols))

    if figsize is None:
        # Estimate figure size
        subplot_w, subplot_h = 4.5, 4
        fig_w = ncols * subplot_w
        fig_h = nrows * subplot_h
        if colorbar_loc == "right":
            fig_w += 1  # Add space for colorbar
        else:  # bottom
            fig_h += 0.8
        figsize = (fig_w, fig_h)

    fig = plt.figure(figsize=figsize)

    # Define grid layout
    if colorbar_loc == "right":
        gs = gridspec.GridSpec(
            nrows, ncols + 1, width_ratios=[1] * ncols + [0.05], wspace=0.3, hspace=0.4
        )
        cbar_gs_indices = (slice(None), -1)
    else:  # bottom
        gs = gridspec.GridSpec(
            nrows + 1, ncols, height_ratios=[1] * nrows + [0.05], wspace=0.3, hspace=0.3
        )
        cbar_gs_indices = (-1, slice(None))

    # --- 3. Plotting on each facet ---
    scatter_collection = []
    for i, group_name in enumerate(groups):
        row, col = i // ncols, i % ncols
        ax = fig.add_subplot(gs[row, col])

        group_df = plot_df[plot_df["group"] == group_name]

        sc = ax.scatter(
            group_df["Dim1"],
            group_df["Dim2"],
            c=group_df["expression"],
            cmap=cmap,
            s=point_size,
            alpha=alpha,
            vmin=vmin,
            vmax=vmax,
            edgecolor="none",
            **kwargs,
        )
        scatter_collection.append(sc)

        ax.set_title(group_name)
        ax.set_xlabel(f"{basis}_1" if show_ticks else "")
        ax.set_ylabel(f"{basis}_2" if show_ticks else "")
        ax.grid(False)
        ax.set_frame_on(frameon)
        if not show_ticks:
            ax.set_xticks([])
            ax.set_yticks([])

    # --- 4. Add Shared Colorbar ---
    cax = fig.add_subplot(gs[cbar_gs_indices])
    fig.colorbar(
        scatter_collection[0],
        cax=cax,
        orientation="vertical" if colorbar_loc == "right" else "horizontal",
    )
    cax.set_title(feature, fontsize=10)

    # --- 5. Final Adjustments ---
    if main_title:
        fig.suptitle(
            main_title, y=0.98 if colorbar_loc == "right" else 0.97, fontsize=16
        )

    # Adjust layout to prevent overlap
    top_space = 0.96 if main_title else 1.0
    rect = [0, 0, 0.92 if colorbar_loc == "right" else 1.0, top_space]
    fig.tight_layout(rect=rect)

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
        log.info(f"Saved faceted feature plot to {save}")

    if show:
        plt.show()

    return fig


def plot_marker_heatmap(
    adata: sc.AnnData,
    markers: Optional[Union[List[str], Dict[str, List[str]]]] = None,
    markers_df: Optional[pd.DataFrame] = None,
    groupby: str = None,
    groupby_main: Optional[str] = None,
    groupby_sub: Optional[str] = None,
    subset: Optional[pd.Series] = None,
    n_genes: int = 5,
    layer: Optional[str] = None,
    use_raw: bool = False,
    standard_scale: Literal["var", "group", None] = "var",
    cmap: Literal["viridis", "Reds", "RdBu_r"] = "RdBu_r",
    show_gene_categories: bool = True,
    cluster_genes: bool = False,  # whether to re-cluster genes within categories
    figsize: Optional[Tuple[float, float]] = None,
    save: Optional[str] = None,
    show: bool = True,
    **kwargs,
) -> Tuple[plt.Figure, pd.DataFrame]:
    """
    Creates a heatmap of marker genes across different groups.

    This function is highly efficient, only extracting the required genes from the
    matrix. It supports combined grouping, on-the-fly subsetting, and
    intelligent gene ordering within categories.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix.
    markers : Optional[Union[List[str], Dict[str, List[str]]], default None
        Genes to plot. Can be a list of genes or a dictionary where keys are
        category names and values are lists of genes.
    markers_df : Optional[pd.DataFrame], default None
        A DataFrame (e.g., from `sc.tl.rank_genes_groups`) with 'group',
        'names', and optionally 'logfoldchanges' columns.
    groupby : str, default None
        The key in `adata.obs` to group the data by. Ignored if
        `groupby_main` and `groupby_sub` are provided.
    groupby_main, groupby_sub : Optional[str], default None
        Two keys in `adata.obs` to create a combined grouping. For example,
        `groupby_main='celltype'` and `groupby_sub='condition'` will create
        groups like 'Fibroblast_Control', 'Fibroblast_Stimulated', etc.
    subset : Optional[pd.Series], default None
        A boolean Series to subset the data before plotting.
    n_genes : int, default 5
        Number of top genes to take from each group if `markers_df` is provided.
    layer : Optional[str], default None
        The layer in `adata.layers` to use for plotting.
    use_raw : bool, default False
        Whether to use `adata.raw` for gene expression data.
    standard_scale : Literal["var", "group", None], default 'var'
        How to standardize the data. 'var' (z-score across genes) is standard.
    cmap : str, default 'viridis'
        The colormap to use for the heatmap.
    show_gene_categories : bool, default True
        Whether to add a color bar on the left to indicate gene categories.
    cluster_genes : bool, default False
        Whether to re-cluster genes within each category.
    figsize : Optional[Tuple[float, float]], default None
        Figure size. If None, it's calculated automatically.
    save : Optional[str], default None
        Path to save the figure.
    show : bool, default True
        Whether to display the plot.
    **kwargs
        Additional keyword arguments passed to `sns.clustermap`.

    Returns
    -------
    Tuple[plt.Figure, pd.DataFrame]
        The matplotlib Figure object and the DataFrame used for plotting (genes x groups).
    """
    adata_to_plot = _subset_adata(adata, subset)

    # 0. Validate Inputs
    if groupby_main and groupby_sub:
        groupby = _combine_groupby(adata_to_plot, groupby_main, groupby_sub)
    elif groupby is None:
        raise ValueError(
            "Either 'groupby' or both 'groupby_main' and 'groupby_sub' must be provided."
        )

    # 1. Parse Markers
    marker_dict = {}

    if markers_df is not None:
        # Expects standard scanpy rank_genes_groups output format
        if "logfoldchanges" in markers_df.columns:
            for group in markers_df["group"].unique():
                grp_df = markers_df[markers_df["group"] == group].sort_values(
                    "logfoldchanges", ascending=False
                )
                marker_dict[group] = grp_df["names"].head(n_genes).tolist()
        else:
            # Fallback
            for group in markers_df["group"].unique():
                marker_dict[group] = (
                    markers_df[markers_df["group"] == group]["names"]
                    .head(n_genes)
                    .tolist()
                )

    elif isinstance(markers, dict):
        marker_dict = markers
    elif isinstance(markers, list):
        marker_dict = {"Markers": markers}
    else:
        raise ValueError("Provide markers or markers_df")

    # 2. Validate Genes
    source = adata_to_plot.raw if (use_raw and adata_to_plot.raw) else adata_to_plot
    valid_genes_set = set()
    for cat, glist in marker_dict.items():
        valid_genes = [g for g in glist if g in source.var_names]
        marker_dict[cat] = valid_genes
        valid_genes_set.update(valid_genes)

    all_genes = list(valid_genes_set)
    if not all_genes:
        raise ValueError("No valid genes found in adata.")

    # 3. Efficient Data Extraction (Subset ONLY needed genes)
    if use_raw and adata_to_plot.raw:
        X_subset = adata_to_plot.raw[:, all_genes].X
    elif layer:
        X_subset = adata_to_plot[:, all_genes].layers[layer]
    else:
        X_subset = adata_to_plot[:, all_genes].X

    if scipy.sparse.issparse(X_subset):
        X_subset = X_subset.toarray()

    # Create DataFrame: Cells x Genes
    expr_df = pd.DataFrame(X_subset, columns=all_genes, index=adata_to_plot.obs.index)

    # 4. Handle Grouping & Aggregation
    groups = adata_to_plot.obs[groupby]
    expr_df["__group__"] = groups.values
    agg_df = expr_df.groupby("__group__").mean()

    if cluster_genes:
        final_gene_order = _sort_genes_within_categories(agg_df, marker_dict)
    else:
        final_gene_order = []
        seen = set()
        for cat, glist in marker_dict.items():
            for g in glist:
                if g not in seen and g in agg_df.columns:
                    final_gene_order.append(g)
                    seen.add(g)

    # 5. Scale (z-score)
    # axis=0: zscore across groups (standard_scale='group')
    # axis=1: zscore across genes (standard_scale='var') - Standard for heatmaps
    if standard_scale == "var":
        # Z-score columns to see how a gene varies across groups
        agg_df = (agg_df - agg_df.mean()) / agg_df.std()
    elif standard_scale == "group":
        # Z-score rows to see how a group varies across genes
        agg_df = agg_df.sub(agg_df.mean(axis=1), axis=0).div(agg_df.std(axis=1), axis=0)

    # 6. Organize Rows (Genes) & Columns (Groups) for plotting
    plot_df = agg_df.T

    # Reorder genes based on category + optional clustering
    final_gene_order = []
    row_colors = []

    # Setup Colors
    if show_gene_categories:
        cat_palette = sns.color_palette("Set2", n_colors=len(marker_dict))
        cat_color_map = dict(zip(marker_dict.keys(), cat_palette))

    for cat, genes in marker_dict.items():
        # Get subset of plot_df
        current_genes = [g for g in genes if g in plot_df.index]
        if not current_genes:
            continue

        if cluster_genes:
            sub = plot_df.loc[current_genes]
            if len(current_genes) > 2:
                try:
                    Z = hierarchy.linkage(sub, method="average", metric="correlation")
                    leaves = hierarchy.leaves_list(Z)
                    current_genes = sub.index[leaves].tolist()
                except:
                    pass  # Keep list order on failure

        final_gene_order.extend(current_genes)
        if show_gene_categories:
            row_colors.extend([cat_color_map[cat]] * len(current_genes))

    # Remove duplicates while keeping order
    seen = set()
    final_gene_order = [x for x in final_gene_order if not (x in seen or seen.add(x))]

    # Apply ordering
    plot_df = plot_df.loc[final_gene_order]

    # Prepare colors Series
    row_colors_series = None
    if show_gene_categories:
        # Re-align colors to the unique gene list (first occurrence wins for color)
        # This is a bit tricky if genes are shared. Assuming markers are mostly unique.
        gene_to_cat = {}
        for c, gs in marker_dict.items():
            for g in gs:
                if g not in gene_to_cat:
                    gene_to_cat[g] = cat_color_map[c]
        row_colors_series = pd.Series(
            [gene_to_cat[g] for g in plot_df.index],
            index=plot_df.index,
            name="Category",
        )

    # 7. Plotting
    if figsize is None:
        figsize = (
            max(6, len(plot_df.columns) * 0.5 + 2),
            max(8, len(plot_df.index) * 0.2),
        )

    g = sns.clustermap(
        plot_df,
        row_cluster=False,  # We handled gene ordering manually
        col_cluster=kwargs.get("col_cluster", True),
        row_colors=row_colors_series,
        cmap=cmap,
        figsize=figsize,
        center=0 if standard_scale else None,
        vmin=-2,
        vmax=2,  # reasonable defaults for z-score
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        **{k: v for k, v in kwargs.items() if k not in ["row_cluster", "col_cluster"]},
    )

    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.tick_params(axis="y", labelsize=8)

    # Add Category Legend manually if needed
    if show_gene_categories:
        handles = [mpatches.Patch(color=c, label=l) for l, c in cat_color_map.items()]
        g.fig.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(0.02, 0.75),
            frameon=False,
            title="Category",
        )

    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")
    if show:
        plt.show()

    return g.fig, plot_df


def plot_ranked_genes(
    adata: sc.AnnData,
    group: str,
    n_genes: int = 10,
    groupby: str = "leiden",
    key: str = "rank_genes_groups",
    figsize: Tuple[float, float] = (8, 5),
    save: Optional[str] = None,
    show: bool = True,
) -> plt.Figure:
    """
    Plot top-ranked genes for a specific group as a horizontal bar chart.

    Parameters
    ----------
    adata : sc.AnnData
        AnnData with rank_genes_groups results
    group : str
        Group name to plot
    n_genes : int
        Number of top genes to show
    groupby : str
        Grouping key
    key : str
        Key in adata.uns for results

    Examples
    --------
    >>> sc.tl.rank_genes_groups(adata, 'leiden')
    >>> plot_ranked_genes(adata, group='0', n_genes=15)
    """
    if key not in adata.uns:
        raise ValueError("Run sc.tl.rank_genes_groups() first")

    # Extract data
    result = adata.uns[key]
    groups = result["names"].dtype.names

    if group not in groups:
        raise ValueError(f"Group '{group}' not found. Available: {groups}")

    genes = result["names"][group][:n_genes]
    scores = result["scores"][group][:n_genes]

    # Plot
    fig, ax = plt.subplots(figsize=figsize)

    y_pos = np.arange(len(genes))
    ax.barh(y_pos, scores, align="center", color="steelblue")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(genes)
    ax.invert_yaxis()
    ax.set_xlabel("Score")
    ax.set_title(f"Top {n_genes} Marker Genes: {group}")

    plt.tight_layout()

    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")
    if show:
        plt.show()

    return fig


# =============================================================================
# Other Visualizations (Volcano, Correlation, etc.)
# =============================================================================


def plot_volcano(
    df: pd.DataFrame,
    x: str = "logfoldchanges",
    y: str = "pvals_adj",
    gene_col: str = "names",
    highlight_genes: Optional[List[str]] = None,
    min_log2fc: float = 1.0,
    max_pval: float = 0.05,
    save: Optional[str] = None,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """Simplified but robust volcano plot."""
    req_cols = [x, y, gene_col]
    if not all(c in df.columns for c in req_cols):
        raise ValueError(f"DataFrame missing columns. Required: {req_cols}")

    df = df[[x, y, gene_col]].dropna().copy()
    if df.empty:
        raise ValueError("No valid data after filtering NaN values")

    df = df[df[y] > 0].copy()

    if df.empty:
        raise ValueError("No valid positive p-values found")

    df["nlog10p"] = -np.log10(df[y].clip(lower=1e-300))

    df["color"] = "grey"
    df.loc[(df[y] < max_pval) & (df[x] > min_log2fc), "color"] = "red"
    df.loc[(df[y] < max_pval) & (df[x] < -min_log2fc), "color"] = "blue"

    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot non-sig first to be at back
    ax.scatter(
        df[df.color == "grey"][x],
        df[df.color == "grey"]["nlog10p"],
        c="grey",
        s=5,
        alpha=0.5,
        rasterized=True,
    )
    # Plot sig
    ax.scatter(
        df[df.color != "grey"][x],
        df[df.color != "grey"]["nlog10p"],
        c=df[df.color != "grey"]["color"],
        s=15,
        alpha=0.7,
        rasterized=True,
    )

    # Thresholds
    ax.axhline(-np.log10(max_pval), ls="--", c="k", lw=0.5)
    ax.axvline(min_log2fc, ls="--", c="k", lw=0.5)
    ax.axvline(-min_log2fc, ls="--", c="k", lw=0.5)

    # Labels
    texts = []
    if highlight_genes:
        subset = df[df[gene_col].isin(highlight_genes)]
        for _, row in subset.iterrows():
            texts.append(
                ax.text(
                    row[x], row["nlog10p"], row[gene_col], fontweight="bold", fontsize=9
                )
            )

    # Add top significant if space permits and no specific highlights
    if not highlight_genes:
        top_genes = df.sort_values(y).head(10)
        for _, row in top_genes.iterrows():
            texts.append(ax.text(row[x], row["nlog10p"], row[gene_col], fontsize=8))

    if adjust_text and texts:
        adjust_text(
            texts, ax=ax, arrowprops=dict(arrowstyle="-", color="black", lw=0.5)
        )

    ax.set_xlabel("Log2 Fold Change")
    ax.set_ylabel("-Log10 Adj. P-value")

    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    return fig


def plot_feature_correlation(
    adata: sc.AnnData,
    features: List[str],
    subset: Optional[pd.Series] = None,
    method: str = "spearman",
    layer: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """Robust feature correlation heatmap."""
    adata_to_plot = _subset_adata(adata, subset)
    # Filter features
    valid = [f for f in features if f in adata_to_plot.var_names]
    if len(valid) < 2:
        raise ValueError("Need at least 2 valid features.")

    # Extract ONLY these features
    if layer:
        X = adata_to_plot[:, valid].layers[layer]
    else:
        X = adata_to_plot[:, valid].X

    if scipy.sparse.issparse(X):
        X = X.toarray()

    df = pd.DataFrame(X, columns=valid)
    corr = df.corr(method=method)

    g = sns.clustermap(
        corr,
        center=0,
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
        figsize=(max(6, len(valid) * 0.4), max(6, len(valid) * 0.4)),
        annot=len(valid) < 15,
        fmt=".2f",
        **kwargs,
    )

    if save:
        plt.savefig(save, bbox_inches="tight")
    if show:
        plt.show()
    return g.fig


def plot_ridge(
    adata: sc.AnnData,
    features: List[str],
    groupby: Optional[str] = None,
    groupby_main: Optional[str] = None,
    groupby_sub: Optional[str] = None,
    subset: Optional[pd.Series] = None,
    layer: Optional[str] = None,
    use_raw: bool = False,
    palette: str = "mako",
    figsize: Optional[Tuple[float, float]] = None,
    main_title: Optional[str] = None,
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """
    Creates a ridge plot to visualize feature distributions across groups.

    This function generates a grid of density plots, where each plot shows the
    distribution of a single feature for a specific category in `groupby`.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data object.
    features : List[str]
        List of features (genes) to plot.
    groupby : Optional[str], default None
        The key in `adata.obs` to group cells by. Ignored if `groupby_main`
        and `groupby_sub` are provided.
    groupby_main, groupby_sub : Optional[str], default None
        Two keys in `adata.obs` to create a combined grouping.
    subset : Optional[pd.Series], default None
        A boolean Series to subset the data before plotting.
    layer : Optional[str], default None
        The layer in `adata.layers` to use for expression values.
    use_raw : bool, default False
        Whether to use `adata.raw` for expression values.
    palette : str, default 'mako'
        The seaborn palette to use for the groups.
    figsize : Optional[Tuple[float, float]], default None
        Figure size. If None, it's calculated automatically.
    main_title : Optional[str], default None
        A main title for the entire figure.
    save : Optional[str], default None
        Path to save the figure.
    dpi : int, default 300
        Resolution for the saved figure.
    show : bool, default True
        Whether to display the plot.
    **kwargs
        Additional keyword arguments passed to `seaborn.kdeplot`.

    Returns
    -------
    plt.Figure
        The matplotlib Figure object.

    Examples
    --------
    >>> # Basic ridge plot
    >>> plot_ridge(adata, features=['CD3D', 'CD4'], groupby='cell_type')

    >>> # Ridge plot with combined grouping
    >>> plot_ridge(adata, features=['LYZ', 'S100A8'], groupby_main='cell_type', groupby_sub='condition')
    """
    # --- 1. Input Validation and Data Preparation ---
    if not (groupby or (groupby_main and groupby_sub)):
        raise ValueError(
            "Provide either 'groupby' or both 'groupby_main' and 'groupby_sub'."
        )

    adata_to_plot = _subset_adata(adata, subset)

    if groupby_main and groupby_sub:
        groupby = _combine_groupby(adata_to_plot, groupby_main, groupby_sub)

    # --- 2. Data Extraction ---
    df_list = []
    for feature in features:
        try:
            # Use sc.get.obs_df which handles raw/layer logic robustly
            expr = sc.get.obs_df(
                adata_to_plot, keys=[feature], use_raw=use_raw, layer=layer
            )[feature].values
        except KeyError:
            raise ValueError(
                f"Feature '{feature}' not found. Check `use_raw` and `layer` parameters."
            )
        temp_df = pd.DataFrame(
            {"expression": expr, "group": adata_to_plot.obs[groupby].values}
        )
        df_list.append(temp_df)

    plot_df = pd.concat(df_list)
    plot_df['feature'] = np.repeat(features, [len(adata_to_plot.obs[groupby].unique())] * len(features))

    # --- 3. Plotting with FacetGrid ---
    n_groups = len(plot_df['group'].unique())
    n_features = len(features)
    
    if figsize is None:
        # Estimate figure size
        subplot_h = 0.8
        subplot_w = 4
        fig_w = n_features * subplot_w
        fig_h = n_groups * subplot_h
        figsize = (fig_w, fig_h)

    # Set a temporary style for ridge plots
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    g = sns.FacetGrid(
        plot_df, row="group", col="feature", sharex=True, sharey=False, height=subplot_h, aspect=subplot_w/subplot_h, palette=palette
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

    g.set_titles("{col_name}")
    g.set_axis_labels("Expression", "")
    g.set(xticks=[], yticks=[])

    # --- 4. Final Adjustments ---
    if main_title:
        g.fig.suptitle(main_title, y=0.98, fontsize=16)

    # Adjust layout to prevent overlap
    g.fig.tight_layout(rect=[0, 0, 1, 0.96])

    if save:
        g.fig.savefig(save, dpi=dpi, bbox_inches="tight")
        log.info(f"Saved ridge plot to {save}")

    if show:
        plt.show()

    return g.fig


def plot_coexpression(
    adata: sc.AnnData,
    x_gene: str,
    y_gene: str,
    color_by: Optional[str] = None,
    subset: Optional[pd.Series] = None,
    layer: Optional[str] = None,
    use_raw: bool = False,
    kind: Literal["scatter", "hexbin"] = "scatter",
    gridsize: int = 50,
    figsize: Tuple[float, float] = (7, 6),
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """
    Creates a scatter or hexbin plot to visualize co-expression of two genes.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data object.
    x_gene : str
        Gene name for x-axis.
    y_gene : str
        Gene name for y-axis.
    color_by : Optional[str], default None
        Column in `adata.obs` or a gene name to color points by.
    subset : Optional[pd.Series], default None
        A boolean Series to subset the data before plotting.
    layer : Optional[str], default None
        The layer in `adata.layers` to use for expression values.
    use_raw : bool, default False
        Whether to use `adata.raw` for expression values.
    kind : Literal["scatter", "hexbin"], default 'scatter'
        The kind of plot to draw. 'scatter' is a standard scatter plot,
        'hexbin' is a 2D histogram for dense data.
    gridsize : int, default 50
        The number of hexagons in the x-direction for the 'hexbin' plot.
        Ignored if `kind='scatter'`.
    figsize : Tuple[float, float], default (7, 6)
        Figure size.
    save : Optional[str], default None
        Path to save the figure.
    dpi : int, default 300
        The resolution in dots per inch for the saved figure.
    show : bool, default True
        Whether to display the plot.
    **kwargs
        Additional keyword arguments passed to `ax.scatter` or `ax.hexbin`.

    Returns
    -------
    plt.Figure
        The matplotlib Figure object.

    Examples
    --------
    >>> # Basic co-expression plot
    >>> plot_coexpression(adata, x_gene='CD3D', y_gene='CD8A')

    >>> # Color by cell type
    >>> plot_coexpression(adata, 'CD3D', 'CD8A', color_by='cell_type')

    >>> # Use hexbin for dense data
    >>> plot_coexpression(adata, 'EPCAM', 'KRT8', kind='hexbin')
    """
    # --- 1. Handle Subsetting and Data Preparation ---
    adata_to_plot = _subset_adata(adata, subset)

    try:
        # Use sc.get.obs_df which handles raw/layer logic robustly
        df = sc.get.obs_df(
            adata_to_plot, keys=[x_gene, y_gene], use_raw=use_raw, layer=layer
        )
    except KeyError as e:
        raise ValueError(f"Could not find data for plotting: {e}")

    if color_by:
        if color_by in adata_to_plot.obs:
            df[color_by] = adata_to_plot.obs[color_by].values
        else:
            try:
                df[color_by] = sc.get.obs_df(
                    adata_to_plot, keys=[color_by], use_raw=use_raw, layer=layer
                )[color_by].values
            except KeyError:
                raise ValueError(f"Color feature '{color_by}' not found.")

    # --- 2. Create plot ---
    fig, ax = plt.subplots(figsize=figsize)

    if kind == "scatter":
        sns.scatterplot(
            data=df,
            x=x_gene,
            y=y_gene,
            hue=color_by,
            s=10,
            alpha=0.7,
            ax=ax,
            edgecolor=None,
            **kwargs,
        )
    elif kind == "hexbin":
        # Hexbin does not support 'hue' directly, so we plot groups sequentially
        if color_by:
            groups = df[color_by].unique()
            cmap = sns.color_palette("tab10", n_colors=len(groups))
            for i, group in enumerate(groups):
                group_df = df[df[color_by] == group]
                ax.hexbin(
                    group_df[x_gene],
                    group_df[y_gene],
                    gridsize=gridsize,
                    cmap=sns.light_palette(cmap[i], as_cmap=True),
                    edgecolor="none",
                    label=group,
                    **kwargs,
                )
            ax.legend()
        else:
            ax.hexbin(
                df[x_gene],
                df[y_gene],
                gridsize=gridsize,
                cmap="viridis",
                edgecolor="none",
                **kwargs,
            )

    ax.set_title(f"Co-expression of {x_gene} and {y_gene}")
    ax.set_xlabel(f"{x_gene} Expression")
    ax.set_ylabel(f"{y_gene} Expression")
    plt.tight_layout()

    # --- 3. Save, Show, and Return ---
    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
        log.info(f"Saved co-expression plot to {save}")

    if show:
        plt.show()

    return fig


def plot_differential_abundance(
    diff_abundance_df: pd.DataFrame,
    pval_threshold: float = 0.05,
    figsize: Tuple[float, float] = (8, 6),
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
) -> plt.Figure:
    """
    Visualize differential abundance results using a volcano-style plot.

    Parameters
    ----------
    diff_abundance_df : pd.DataFrame
        DataFrame with columns: 'cell_type', 'log2fc_abundance', 'pvalue', 'mean_abundance_group1'
    pval_threshold : float, default=0.05
        P-value threshold for significance
    figsize : Tuple[float, float], default=(8, 6)
        Figure size (width, height) in inches
    save : Optional[str], default=None
        Path to save the figure
    dpi : int, default=300
        Resolution for saved figure
    show : bool, default=True
        Whether to display the plot

    Returns
    -------
    plt.Figure
        Figure object

    Examples
    --------
    >>> # Visualize differential abundance
    >>> plot_differential_abundance(diff_abundance_df)
    """
    df = diff_abundance_df.copy()
    df["-log10(pvalue)"] = -np.log10(df["pvalue"])
    df["significant"] = df["pvalue"] < pval_threshold

    fig, ax = plt.subplots(figsize=figsize)

    sns.scatterplot(
        data=df,
        x="log2fc_abundance",
        y="-log10(pvalue)",
        hue="significant",
        palette={True: "red", False: "grey"},
        size="mean_abundance_group1",
        sizes=(50, 500),
        alpha=0.7,
        ax=ax,
    )

    # Add labels for significant points
    for i, row in df[df["significant"]].iterrows():
        ax.text(
            row["log2fc_abundance"],
            row["-log10(pvalue)"],
            row["cell_type"],
            fontsize=9,
        )

    ax.set_title("Differential Cell Type Abundance")
    ax.set_xlabel("Log2 Fold Change in Abundance")
    ax.set_ylabel("-log10(p-value)")
    ax.axhline(-np.log10(pval_threshold), ls="--", color="grey")
    ax.axvline(0, ls="--", color="grey")
    ax.legend(title="Significant")

    plt.tight_layout()

    if save:
        plt.savefig(save, dpi=dpi, bbox_inches="tight")
        log.info(f"Saved differential abundance plot to {save}")

    if show:
        plt.show()

    return fig
