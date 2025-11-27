"""
Dimensionality reduction and visualization functions for single-cell RNA-seq data.
Revised for robustness and efficiency.
"""

import logging
from typing import Dict, List, Literal, Optional, Tuple, Union, Any

import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
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

def _process_cmap(cmap_arg):
    """Helper to convert color list to Colormap."""
    if isinstance(cmap_arg, list):
        return mcolors.LinearSegmentedColormap.from_list("custom_cmap", cmap_arg)
    return cmap_arg


def _get_palette_map(
    adata: sc.AnnData,
    key: str,
    palette: Optional[Union[str, Dict]] = None
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
        return {cat: palette.get(cat, '#cccccc') for cat in categories}

    # 2. Check adata.uns (Scanpy convention: key_colors)
    uns_key = f"{key}_colors"
    if uns_key in adata.uns:
        colors = adata.uns[uns_key]
        if len(colors) >= len(categories):
            return dict(zip(categories, colors[:len(categories)]))

    # 3. Generate new palette (Seaborn/Matplotlib)
    # If user provided a string palette name (e.g., 'tab20'), use it
    palette_name = palette if isinstance(palette, str) else 'tab20'
    if len(categories) > 20 and palette_name == 'tab20':
        palette_name = 'husl' # Fallback for many categories
    
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
        sub_df_T = sub_df.T # genes x cells
        
        sub_df_T = sub_df_T.loc[sub_df_T.var(axis=1) > 0]
        if len(sub_df_T) < 2:
                sorted_genes.extend(valid_genes)
                continue
        try:
            Z = hierarchy.linkage(
                distance.pdist(sub_df_T.values, metric='correlation'),
                method='average'
            )
            leaves = hierarchy.dendrogram(Z, no_plot=True)['leaves']
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
    title: Optional[str] = None,
    show_labels: bool = True,
    palette: Optional[Union[str, Dict[str, str]]] = None,
    size: float = 12,
    alpha: float = 0.8,
    ncols: int = 3,
    figsize: Optional[Tuple[float, float]] = None,
    legend_loc: str = "right margin",
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
    Enhanced wrapper for scanpy's embedding plot.
    """
    # Recursive call for list handling
    if isinstance(color_by, list):
        if len(color_by) == 1:
            return plot_embedding(adata, color_by[0], basis, title, show_labels, palette, size, alpha, 
                                  ncols, figsize, legend_loc, label_size, save, dpi, ax, show, legend_style, rasterized, **kwargs)
        
        n_plots = len(color_by)
        n_rows = int(np.ceil(n_plots / ncols))
        if figsize is None:
            figsize = (4 * ncols, 3.5 * n_rows)
            
        fig, axes = plt.subplots(n_rows, ncols, figsize=figsize, squeeze=False)
        axes = axes.flatten()
        
        for i, color in enumerate(color_by):
            plot_embedding(adata, color, basis, ax=axes[i], show=False, save=None, 
                           show_labels=show_labels, palette=palette, size=size, alpha=alpha, **kwargs)
            
        for i in range(n_plots, len(axes)):
            axes[i].axis('off')
            
        plt.tight_layout()
        if save:
            plt.savefig(save, dpi=dpi, bbox_inches='tight')
        if show:
            plt.show()
        return fig

    # Standard Plotting Logic
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize or (6, 5))
    else:
        fig = ax.figure

    # Basic plot
    sc.pl.embedding(
        adata,
        basis=basis,
        color=color_by,
        ax=ax,
        show=False,
        size=size,
        alpha=alpha,
        palette=palette,
        legend_loc="none" if legend_style == "on_data" else legend_loc,
        **kwargs
    )
    
    if rasterized:
        for collection in ax.collections:
            collection.set_rasterized(True)

    # Label Logic (Only for categorical data)
    is_categorical = pd.api.types.is_categorical_dtype(adata.obs[color_by]) if color_by in adata.obs else False
    
    if is_categorical and show_labels and legend_style in ["on_data", "both"]:
        # Get color mapping for halo
        color_map = _get_palette_map(adata, color_by, palette)
        categories = adata.obs[color_by].cat.categories
        
        embed_key = f"X_{basis}"
        texts = []
        
        for label in categories:
            mask = adata.obs[color_by] == label
            if not np.any(mask):
                continue
                
            # Calculate centroid
            coords = adata.obsm[embed_key][mask]
            x, y = np.median(coords, axis=0)
            
            # Halo color
            bg_color = color_map.get(label, 'white')
            
            txt = ax.text(
                x, y, label,
                fontsize=label_size,
                fontweight="bold",
                ha="center", va="center",
                color="black"
            )
            txt.set_path_effects([
                PathEffects.withStroke(linewidth=3, foreground=bg_color, alpha=0.8)
            ])
            texts.append(txt)
            
        if adjust_text and texts:
            try:
                adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
            except Exception as e:
                log.debug(f"adjust_text failed: {e}")

    # Title and cleanup
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
    col_wrap: int = 4,
    point_size: float = 10,  # Increased default size
    alpha: float = 0.8,
    palette: Optional[Union[str, Dict[str, str]]] = None,
    figsize: Tuple[float, float] = (12, 10),
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
    **kwargs,
) -> plt.Figure:
    """Faceted embedding plots."""
    embed_key = f"X_{basis}" if not basis.startswith("X_") else basis
    
    if embed_key not in adata.obsm:
        raise ValueError(f"{embed_key} not found in obsm")
        
    # Prepare data
    coords = adata.obsm[embed_key][:, :2]
    plot_df = pd.DataFrame(coords, columns=["Dim1", "Dim2"])
    plot_df[split_by] = adata.obs[split_by].values
    plot_df[color_by] = adata.obs[color_by].values
    
    # Get robust palette
    final_palette = palette
    if final_palette is None and pd.api.types.is_categorical_dtype(adata.obs[color_by]):
        final_palette = _get_palette_map(adata, color_by, None)

    g = sns.FacetGrid(plot_df, col=split_by, col_wrap=col_wrap, sharex=True, sharey=True)
    
    g.map_dataframe(
        sns.scatterplot,
        x="Dim1", y="Dim2",
        hue=color_by,
        palette=final_palette,
        s=point_size,
        alpha=alpha,
        linewidth=0,
        **kwargs
    )
    
    g.add_legend(title=color_by, bbox_to_anchor=(1.01, 0.5), loc='center left')
    g.set_axis_labels(f"{basis}_1", f"{basis}_2")
    g.set_titles("{col_name}")
    
    for ax in g.axes.flat:
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        
    g.fig.set_size_inches(figsize)
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
        
    return g.fig


# =============================================================================
# Marker Visualization
# =============================================================================

def plot_dotplot(adata, markers, groupby, **kwargs):
    """Wrapper for scanpy dotplot."""
    # Ensure kwargs doesn't contain show/save which we handle
    show = kwargs.pop('show', True)
    save = kwargs.pop('save', None)
    
    dp = sc.pl.dotplot(adata, markers, groupby=groupby, show=False, **kwargs)
    
    # sc.pl.dotplot returns a DotPlot object, not Figure
    fig = dp['mainplot_ax'].figure if isinstance(dp, dict) else plt.gcf()
    
    if save:
        plt.savefig(save, bbox_inches='tight', dpi=300)
    if show:
        plt.show()
    return fig


def plot_stacked_violin(adata, markers, groupby, **kwargs):
    """Wrapper for scanpy stacked_violin."""
    show = kwargs.pop('show', True)
    save = kwargs.pop('save', None)
    
    sc.pl.stacked_violin(adata, markers, groupby=groupby, show=False, **kwargs)
    fig = plt.gcf()
    
    if save:
        plt.savefig(save, bbox_inches='tight', dpi=300)
    if show:
        plt.show()
    return fig


def plot_split_violin_with_stats(
    adata: sc.AnnData,
    genes: Union[str, List[str]],
    celltype_col: str,
    condition_col: str,
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
    """Split violin plot with statistical annotation."""
    try:
        from statannotations.Annotator import Annotator
    except ImportError:
        log.error("statannotations not installed.")
        raise ImportError("Please install statannotations: pip install statannotations")

    if isinstance(genes, str):
        genes = [genes]
        
    # Data extraction logic revised for safety
    if use_raw and adata.raw is not None:
        source = adata.raw
    else:
        source = adata
        
    # Check availability
    available_genes = [g for g in genes if g in source.var_names]
    if not available_genes:
        raise ValueError("No valid genes found.")

    df = sc.get.obs_df(adata, keys=[celltype_col, condition_col] + available_genes, use_raw=use_raw, layer=layer)
    
    if figsize is None:
        figsize = (4 * len(available_genes), 5)
        
    fig, axes = plt.subplots(1, len(available_genes), figsize=figsize, squeeze=False)
    axes = axes.flatten()
    
    groups = df[condition_col].unique()
    if len(groups) != 2:
        raise ValueError(f"'{condition_col}' must have exactly 2 groups, found {len(groups)}")
        
    for i, gene in enumerate(available_genes):
        ax = axes[i]
        sns.violinplot(
            data=df, x=celltype_col, y=gene, hue=condition_col,
            split=True, inner="quartile", palette=palette, ax=ax, cut=0, **kwargs
        )
        ax.set_title(gene)
        ax.set_xlabel("")
        ax.tick_params(axis='x', rotation=45)
        
        # Stats
        box_pairs = [((ct, groups[0]), (ct, groups[1])) for ct in df[celltype_col].unique()]
        try:
            annot = Annotator(ax, box_pairs, data=df, x=celltype_col, y=gene, hue=condition_col)
            annot.configure(test=test, text_format="star", loc="inside", verbose=0)
            annot.apply_and_annotate()
        except Exception as e:
            log.warning(f"Stats failed for {gene}: {e}")
            
        # Legend cleanup
        if i < len(available_genes) - 1 and ax.get_legend():
            ax.get_legend().remove()
            
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=dpi)
    if show:
        plt.show()
    return fig


def plot_marker_expression(
    adata: sc.AnnData,
    markers: Union[str, List[str]],
    basis: str = "umap",
    ncols: int = 4,
    use_raw: bool = False,
    layer: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = True,
    **kwargs
) -> plt.Figure:
    """
    Grid wrapper for sc.pl.embedding to show gene expression.
    """
    if isinstance(markers, str):
        markers = [markers]
        
    # Clean markers
    valid_markers = []
    source = adata.raw if (use_raw and adata.raw) else adata
    for m in markers:
        if m in source.var_names:
            valid_markers.append(m)
        else:
            log.warning(f"Marker {m} not found.")
            
    if not valid_markers:
        raise ValueError("No valid markers.")
        
    n_plots = len(valid_markers)
    n_rows = int(np.ceil(n_plots / ncols))
    
    fig, axes = plt.subplots(n_rows, ncols, figsize=(4*ncols, 3.5*n_rows), squeeze=False)
    axes = axes.flatten()
    
    for i, m in enumerate(valid_markers):
        sc.pl.embedding(
            adata, basis=basis, color=m, ax=axes[i], show=False,
            use_raw=use_raw, layer=layer, legend_loc='none', **kwargs
        )
        axes[i].set_title(m)
        
    for i in range(n_plots, len(axes)):
        axes[i].axis('off')
        
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=300)
    if show:
        plt.show()
    return fig


def plot_faceted_feature(
    adata: sc.AnnData,
    feature: str,
    split_by: str,
    basis: str = "umap",
    col_wrap: int = 4,
    cmap: str = "viridis",
    point_size: float = 2,
    alpha: float = 0.8,
    use_raw: bool = True,
    layer: Optional[str] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    sort_order: bool = True,
    figsize: Tuple[float, float] = (12, 10),
    save: Optional[str] = None,
    show: bool = True,
) -> plt.Figure:
    """Faceted feature plot with improved colorbar and memory handling."""
    embed_key = f"X_{basis}" if not basis.startswith("X_") else basis
    
    # 1. Get Expression Data safely
    if feature in adata.obs:
        vals = adata.obs[feature].values
    else:
        # Use sc.get.obs_df which handles raw/layer logic nicely
        try:
            vals = sc.get.obs_df(adata, keys=[feature], use_raw=use_raw, layer=layer)[feature].values
        except KeyError:
            raise ValueError(f"Feature {feature} not found.")

    # 2. Build Dataframe
    coords = adata.obsm[embed_key][:, :2]
    plot_df = pd.DataFrame({
        "Dim1": coords[:, 0],
        "Dim2": coords[:, 1],
        "group": adata.obs[split_by].values,
        "expression": vals
    })

    if sort_order:
        plot_df = plot_df.sort_values("expression")
        
    # Auto-scale
    if vmin is None: vmin = np.percentile(plot_df['expression'], 1) if sort_order else 0
    if vmax is None: vmax = np.percentile(plot_df['expression'], 99)
    
    # 3. Plot
    g = sns.FacetGrid(plot_df, col="group", col_wrap=col_wrap, height=3.5)
    
    cmap_obj = plt.get_cmap(cmap)
    
    def scatter_func(x, y, c, **kws):
        plt.scatter(x, y, c=c, cmap=cmap_obj, s=point_size, alpha=alpha, vmin=vmin, vmax=vmax, edgecolor='none')
        
    g.map(scatter_func, "Dim1", "Dim2", "expression")
    
    # 4. Styling & Colorbar
    g.set_titles("{col_name}")
    g.set_axis_labels("", "")
    for ax in g.axes.flat:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal')

    # Safely add colorbar to the figure, not inside axes
    g.fig.subplots_adjust(right=0.92)
    n_cols = min(col_wrap, len(plot_df['group'].unique()))
    if n_cols <= 2:
        cbar_x = 0.93
    elif n_cols <= 3:
        cbar_x = 0.94
    else:
        cbar_x = 0.95
    cax = g.fig.add_axes([cbar_x, 0.3, 0.015, 0.4])
    
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=cmap_obj, norm=norm)
    g.fig.colorbar(sm, cax=cax, label=feature)
    
    g.fig.set_size_inches(figsize)
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
        
    return g.fig


def plot_marker_heatmap(
    adata: sc.AnnData,
    markers: Optional[Union[List[str], Dict[str, List[str]]]] = None,
    markers_df: Optional[pd.DataFrame] = None,
    groupby: str = None,
    n_genes: int = 5,
    layer: Optional[str] = None,
    use_raw: bool = False,
    standard_scale: Literal["var", "group", None] = "var",
    cmap: str = "viridis",
    show_gene_categories: bool = True,
    cluster_genes: bool = False,  # whether to re-cluster genes within categories
    figsize: Optional[Tuple[float, float]] = None,
    save: Optional[str] = None,
    show: bool = True,
    **kwargs
) -> Tuple[plt.Figure, pd.DataFrame]:
    """
    Refactored marker heatmap for better memory usage and logic.
    """
    if groupby is None:
        raise ValueError("groupby is required")
        
    # 1. Parse Markers
    marker_dict = {}
    
    if markers_df is not None:
        # Expects standard scanpy rank_genes_groups output format
        if 'logfoldchanges' in markers_df.columns:
             for group in markers_df['group'].unique():
                 grp_df = markers_df[markers_df['group'] == group].sort_values('logfoldchanges', ascending=False)
                 marker_dict[group] = grp_df['names'].head(n_genes).tolist()
        else:
             # Fallback
             for group in markers_df['group'].unique():
                 marker_dict[group] = markers_df[markers_df['group'] == group]['names'].head(n_genes).tolist()
    
    elif isinstance(markers, dict):
        marker_dict = markers
    elif isinstance(markers, list):
        marker_dict = {'Markers': markers}
    else:
        raise ValueError("Provide markers or markers_df")

    # 2. Validate Genes
    source = adata.raw if (use_raw and adata.raw) else adata
    valid_genes_set = set()
    for cat, glist in marker_dict.items():
        valid_genes = [g for g in glist if g in source.var_names]
        marker_dict[cat] = valid_genes
        valid_genes_set.update(valid_genes)
        
    all_genes = list(valid_genes_set)
    if not all_genes:
        raise ValueError("No valid genes found in adata.")

    # 3. Efficient Data Extraction (Subset ONLY needed genes)
    # This avoids converting the whole matrix to dense
    if use_raw and adata.raw:
        X_subset = adata.raw[:, all_genes].X
    elif layer:
        X_subset = adata[:, all_genes].layers[layer]
    else:
        X_subset = adata[:, all_genes].X
        
    if scipy.sparse.issparse(X_subset):
        X_subset = X_subset.toarray()
        
    # Create DataFrame: Cells x Genes
    expr_df = pd.DataFrame(X_subset, columns=all_genes, index=adata.obs.index)
    
    # 4. Handle Grouping & Aggregation
    if isinstance(groupby, str):
        groups = adata.obs[groupby]
    else:
        # Multi-index grouping
        groups = adata.obs[groupby].apply(lambda x: '_'.join(x.astype(str)), axis=1)
        
    expr_df['__group__'] = groups.values
    agg_df = expr_df.groupby('__group__').mean()
    
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
    if standard_scale == 'var':
        from scipy.stats import zscore
        # We need to apply zscore to the aggregated data (Groups x Genes)
        # Usually we want to see how a gene varies across groups -> zscore columns
        #agg_df = agg_df.apply(zscore, axis=0) 
        agg_df = (agg_df - agg_df.mean()) / agg_df.std()
    elif standard_scale == 'group':
        # See how a group varies across genes (rare)
        #agg_df = agg_df.apply(zscore, axis=1)
        agg_df = agg_df.sub(agg_df.mean(axis=1), axis=0).div(agg_df.std(axis=1), axis=0)


    # 6. Organize Rows (Genes) & Columns (Groups) for plotting
    # We want Genes as Rows, Groups as Columns
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
        if not current_genes: continue
            
        if cluster_genes:
            # Cluster just these genes based on their expression across groups
            # Use original expr_df for correlation or the aggregated one? 
            # Aggregated is faster and usually sufficient for marker heatmaps
            sub = plot_df.loc[current_genes]
            if len(current_genes) > 2:
                try:
                    Z = hierarchy.linkage(sub, method='average', metric='correlation')
                    leaves = hierarchy.leaves_list(Z)
                    current_genes = sub.index[leaves].tolist()
                except:
                    pass # Keep list order on failure
        
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
                if g not in gene_to_cat: gene_to_cat[g] = cat_color_map[c]
        row_colors_series = pd.Series([gene_to_cat[g] for g in plot_df.index], index=plot_df.index, name="Category")

    # 7. Plotting
    if figsize is None:
        figsize = (max(6, len(plot_df.columns)*0.5 + 2), max(8, len(plot_df.index)*0.2))

    g = sns.clustermap(
        plot_df,
        row_cluster=False, # We handled gene ordering manually
        col_cluster=kwargs.get('col_cluster', True),
        row_colors=row_colors_series,
        cmap=cmap,
        figsize=figsize,
        center=0 if standard_scale else None,
        vmin=-2, vmax=2, # reasonable defaults for z-score
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        **{k:v for k,v in kwargs.items() if k not in ['row_cluster', 'col_cluster']}
    )
    
    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.tick_params(axis='y', labelsize=8)
    
    # Add Category Legend manually if needed
    if show_gene_categories:
        handles = [mpatches.Patch(color=c, label=l) for l, c in cat_color_map.items()]
        g.fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.02, 0.75), frameon=False, title="Category")

    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
        
    return g.fig, plot_df


def plot_ranked_genes(
    adata: sc.AnnData,
    group: str,
    n_genes: int = 10,
    groupby: str = 'leiden',
    key: str = 'rank_genes_groups',
    figsize: Tuple[float, float] = (8, 5),
    save: Optional[str] = None,
    show: bool = True
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
        raise ValueError(f"Run sc.tl.rank_genes_groups() first")
    
    # Extract data
    result = adata.uns[key]
    groups = result['names'].dtype.names
    
    if group not in groups:
        raise ValueError(f"Group '{group}' not found. Available: {groups}")
    
    genes = result['names'][group][:n_genes]
    scores = result['scores'][group][:n_genes]
    
    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    
    y_pos = np.arange(len(genes))
    ax.barh(y_pos, scores, align='center', color='steelblue')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(genes)
    ax.invert_yaxis()
    ax.set_xlabel('Score')
    ax.set_title(f'Top {n_genes} Marker Genes: {group}')
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
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
    **kwargs
) -> plt.Figure:
    """Simplified but robust volcano plot."""
    # ... (Logic largely similar to original but ensure column existence)
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
    
    df['color'] = 'grey'
    df.loc[(df[y] < max_pval) & (df[x] > min_log2fc), 'color'] = 'red'
    df.loc[(df[y] < max_pval) & (df[x] < -min_log2fc), 'color'] = 'blue'
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot non-sig first to be at back
    ax.scatter(df[df.color=='grey'][x], df[df.color=='grey']['nlog10p'], c='grey', s=5, alpha=0.5, rasterized=True)
    # Plot sig
    ax.scatter(df[df.color!='grey'][x], df[df.color!='grey']['nlog10p'], c=df[df.color!='grey']['color'], s=15, alpha=0.7, rasterized=True)
    
    # Thresholds
    ax.axhline(-np.log10(max_pval), ls='--', c='k', lw=0.5)
    ax.axvline(min_log2fc, ls='--', c='k', lw=0.5)
    ax.axvline(-min_log2fc, ls='--', c='k', lw=0.5)
    
    # Labels
    texts = []
    if highlight_genes:
        subset = df[df[gene_col].isin(highlight_genes)]
        for _, row in subset.iterrows():
            texts.append(ax.text(row[x], row['nlog10p'], row[gene_col], fontweight='bold', fontsize=9))
    
    # Add top significant if space permits and no specific highlights
    if not highlight_genes:
        top_genes = df.sort_values(y).head(10)
        for _, row in top_genes.iterrows():
            texts.append(ax.text(row[x], row['nlog10p'], row[gene_col], fontsize=8))

    if adjust_text and texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))
        
    ax.set_xlabel("Log2 Fold Change")
    ax.set_ylabel("-Log10 Adj. P-value")
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    return fig


def plot_feature_correlation(
    adata: sc.AnnData,
    features: List[str],
    method: str = "spearman",
    layer: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = True,
    **kwargs
) -> plt.Figure:
    """Robust feature correlation heatmap."""
    # Filter features
    valid = [f for f in features if f in adata.var_names]
    if len(valid) < 2:
        raise ValueError("Need at least 2 valid features.")
        
    # Extract ONLY these features
    if layer:
        X = adata[:, valid].layers[layer]
    else:
        X = adata[:, valid].X
        
    if scipy.sparse.issparse(X):
        X = X.toarray()
        
    df = pd.DataFrame(X, columns=valid)
    corr = df.corr(method=method)
    
    g = sns.clustermap(
        corr, 
        center=0, cmap="coolwarm", vmin=-1, vmax=1,
        figsize=(max(6, len(valid)*0.4), max(6, len(valid)*0.4)),
        annot=len(valid)<15, fmt=".2f",
        **kwargs
    )
    
    if save:
        plt.savefig(save, bbox_inches='tight')
    if show:
        plt.show()
    return g.fig


def plot_ridge(
    adata: sc.AnnData,
    features: List[str],
    groupby: str,
    layer: Optional[str] = None,
    use_raw: bool = False,
    figsize: Optional[Tuple[float, float]] = None,
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
) -> plt.Figure:
    """
    Create a ridge plot to visualize feature distributions across groups.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data object
    features : List[str]
        List of features (genes) to plot
    groupby : str
        Column in adata.obs to group cells by
    layer : Optional[str], default=None
        Layer in adata.layers to use for expression values
    use_raw : bool, default=False
        Whether to use adata.raw for expression values
    figsize : Optional[Tuple[float, float]], default=None
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
    >>> plot_ridge(
    ...     adata,
    ...     ['CD3D', 'CD4', 'CD8A'],
    ...     groupby='celltype'
    ... )
    """
    log.info(f"Creating ridge plot for {len(features)} features across '{groupby}'")

    # Extract data
    df_list = []
    for feature in features:
        if use_raw and adata.raw is not None:
            expr = adata.raw[:, feature].X
        elif layer:
            expr = adata[:, feature].layers[layer]
        else:
            expr = adata[:, feature].X

        if scipy.sparse.issparse(expr):
            expr = expr.toarray().flatten()
        else:
            expr = expr.flatten()

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

    if save:
        plt.savefig(save, dpi=dpi, bbox_inches="tight")
        log.info(f"Saved ridge plot to {save}")

    if show:
        plt.show()

    return g.fig


def plot_coexpression(
    adata: sc.AnnData,
    x_gene: str,
    y_gene: str,
    color_by: Optional[str] = None,
    layer: Optional[str] = None,
    use_raw: bool = False,
    figsize: Tuple[float, float] = (7, 6),
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
) -> plt.Figure:
    """
    Create a scatter plot to visualize co-expression of two genes.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data object
    x_gene : str
        Gene name for x-axis
    y_gene : str
        Gene name for y-axis
    color_by : Optional[str], default=None
        Column in adata.obs or gene name to color points by
    layer : Optional[str], default=None
        Layer in adata.layers to use for expression values
    use_raw : bool, default=False
        Whether to use adata.raw for expression values
    figsize : Tuple[float, float], default=(7, 6)
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
    >>> # Basic co-expression plot
    >>> plot_coexpression(adata, 'CD3D', 'CD8A')
    >>>
    >>> # Color by cell type
    >>> plot_coexpression(
    ...     adata,
    ...     'CD3D',
    ...     'CD8A',
    ...     color_by='celltype'
    ... )
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

    if save:
        plt.savefig(save, dpi=dpi, bbox_inches="tight")
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
