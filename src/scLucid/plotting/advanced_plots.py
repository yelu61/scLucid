"""
Plotting functions for single-cell RNA-seq data.
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

log = logging.getLogger(__name__)


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
