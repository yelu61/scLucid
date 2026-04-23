"""
Plotting functions for single-cell RNA-seq data.
"""

import logging
from typing import Dict, List, Literal, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

# Try importing adjustText softly
try:
    from adjustText import adjust_text
except ImportError:
    adjust_text = None

log = logging.getLogger(__name__)


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

    Returns:
    -------
    sc.pl.DotPlot
        The DotPlot object from scanpy, which can be used for further
        customization (e.g., `dp.add_totals()`).
    """
    # --- 1. Input Validation ---
    if not (groupby or (groupby_main and groupby_sub)):
        raise ValueError("Provide either 'groupby' or both 'groupby_main' and 'groupby_sub'.")

    # --- 2. Handle Subsetting and Grouping using Helper Functions ---
    adata_to_plot = _subset_adata(adata, subset)

    if groupby_main and groupby_sub:
        groupby = _combine_groupby(adata_to_plot, groupby_main, groupby_sub)

    # --- 3. Handle Category Ordering ---
    if auto_order_categories and "categories_order" not in kwargs:
        if pd.api.types.is_categorical_dtype(adata_to_plot.obs[groupby]):
            kwargs["categories_order"] = adata_to_plot.obs[groupby].cat.categories.tolist()

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

    Returns:
    -------
    plt.Figure
        The matplotlib Figure object.

    Examples:
    --------
    >>> # Standard stacked violin plot
    >>> plot_stacked_violin(adata, var_names=['LYZ', 'S100A8'], groupby='cell_type')

    >>> # Stacked violin plot with combined grouping
    >>> plot_stacked_violin(adata, var_names=['CD3D', 'MS4A1'], groupby_main='cell_type', groupby_sub='condition')
    """
    # --- 1. Input Validation and Data Preparation ---
    if not (groupby or (groupby_main and groupby_sub)):
        raise ValueError("Provide either 'groupby' or both 'groupby_main' and 'groupby_sub'.")

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

    Returns:
    -------
    plt.Figure
        The matplotlib Figure object.

    Raises:
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
        raise ValueError(f"'{condition_col}' must have exactly 2 groups, found {len(groups)}")

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
        box_pairs = [((ct, groups[0]), (ct, groups[1])) for ct in df[celltype_col].unique()]
        try:
            annot = Annotator(ax, box_pairs, data=df, x=celltype_col, y=gene, hue=condition_col)
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

    Returns:
    -------
    plt.Figure
        The matplotlib Figure object.

    Examples:
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

    fig, axes = plt.subplots(n_rows, ncols, figsize=(4 * ncols, 3.5 * n_rows), squeeze=False)
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
            vmin=vmin,  # Apply consistent color scale
            vmax=vmax,  # Apply consistent color scale
            cmap=cmap,  # Apply consistent colormap
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
