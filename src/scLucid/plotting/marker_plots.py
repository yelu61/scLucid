"""
Plotting functions for single-cell RNA-seq data.
"""

import logging
from typing import Dict, List, Literal, Optional, Tuple, Union

import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import seaborn as sns
from scipy.cluster import hierarchy

# Try importing adjustText softly
try:
    from adjustText import adjust_text
except ImportError:
    adjust_text = None

log = logging.getLogger(__name__)


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

    Returns:
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

    log.info(f"Plotting feature '{feature}' for {n_groups} groups: {', '.join(map(str, groups))}")

    # Get coordinates
    coords = adata.obsm[embed_key][:, :2]
    plot_df = pd.DataFrame(coords, columns=["Dim1", "Dim2"])
    plot_df["group"] = pd.Categorical(adata.obs[split_by], categories=groups, ordered=True)

    # Get feature expression
    try:
        expr = sc.get.obs_df(adata, keys=[feature], use_raw=use_raw, layer=layer)[feature].values
    except KeyError:
        raise ValueError(f"Feature '{feature}' not found. Check `use_raw` and `layer` parameters.")
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
        fig.suptitle(main_title, y=0.98 if colorbar_loc == "right" else 0.97, fontsize=16)

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

    Returns:
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
                    markers_df[markers_df["group"] == group]["names"].head(n_genes).tolist()
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

    Examples:
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
                ax.text(row[x], row["nlog10p"], row[gene_col], fontweight="bold", fontsize=9)
            )

    # Add top significant if space permits and no specific highlights
    if not highlight_genes:
        top_genes = df.sort_values(y).head(10)
        for _, row in top_genes.iterrows():
            texts.append(ax.text(row[x], row["nlog10p"], row[gene_col], fontsize=8))

    if adjust_text and texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="black", lw=0.5))

    ax.set_xlabel("Log2 Fold Change")
    ax.set_ylabel("-Log10 Adj. P-value")

    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    return fig
