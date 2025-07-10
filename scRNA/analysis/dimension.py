"""
Dimensionality reduction and visualization functions for single-cell RNA-seq data.

This module provides functions for visualizing marker gene expression and
cell type compositions across clusters.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from adjustText import adjust_text
from typing import Optional, List, Dict, Tuple, Literal, Union
from matplotlib.colors import LinearSegmentedColormap

from .manager import Manager


def plot_marker_expression(
    adata,
    marker_config: str,
    cell_types: Optional[List[str]] = None,
    basis: str = 'umap',
    n_markers: int = 3,
    figsize: Tuple[float, float] = (12, 10),
    ncols: int = 3,
    color_map: str = 'viridis',
    use_raw: bool = False,
    layer: Optional[str] = None,
    show_legend: bool = True,
    save: Optional[str] = None,
    **kwargs
) -> None:
    """
    Visualize expression of marker genes for selected cell types.
    
    Args:
        adata: AnnData object
        marker_config: Path to marker configuration file
        cell_types: List of cell types to show markers for. If None, use all major cell types
        basis: Basis for embedding visualization (e.g., 'umap', 'tsne')
        n_markers: Number of top markers to show for each cell type
        figsize: Figure size
        ncols: Number of columns in the grid
        color_map: Colormap for gene expression
        use_raw: Whether to use raw data for expression values
        layer: Layer to use for expression values
        show_legend: Whether to show legend
        save: Path to save the figure
        **kwargs: Additional arguments to pass to sc.pl.embedding
    """
    # Load marker config
    mgr = Manager(marker_config)
    mgr.intersect_with(adata)
    
    # Get cell types
    if cell_types is None:
        cell_types = [
            cell_type for cell_type, cell in mgr.CELLS.items()
            if cell.level == 'major' and len(cell.markers) >= n_markers
        ]
    
    # Filter cell types with too few markers
    filtered_cell_types = []
    for cell_type in cell_types:
        if cell_type not in mgr.CELLS:
            print(f"Warning: Cell type '{cell_type}' not found in marker config")
            continue
        if len(mgr.CELLS[cell_type].markers) < n_markers:
            print(f"Warning: Cell type '{cell_type}' has fewer than {n_markers} markers")
            continue
        filtered_cell_types.append(cell_type)
    
    cell_types = filtered_cell_types
    
    if not cell_types:
        print("No valid cell types with sufficient markers found")
        return
    
    # Calculate number of plots
    n_plots = len(cell_types) * n_markers
    nrows = int(np.ceil(n_plots / ncols))
    
    # Create figure
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    if nrows * ncols == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # Plot markers
    plot_idx = 0
    for i, cell_type in enumerate(cell_types):
        # Get markers for this cell type
        markers = mgr.CELLS[cell_type].markers[:n_markers]
        cell_color = mgr.CELLS[cell_type].color or 'black'
        
        for j, marker in enumerate(markers):
            if plot_idx < len(axes):
                if marker in adata.var_names:
                    # Plot expression
                    sc.pl.embedding(
                        adata,
                        basis=basis,
                        color=marker,
                        use_raw=use_raw,
                        layer=layer,
                        ax=axes[plot_idx],
                        show=False,
                        cmap=color_map,
                        title=f"{cell_type}: {marker}",
                        **kwargs
                    )
                    
                    # Add colored border to indicate cell type
                    for spine in axes[plot_idx].spines.values():
                        spine.set_edgecolor(cell_color)
                        spine.set_linewidth(2)
                else:
                    axes[plot_idx].text(
                        0.5, 0.5, 
                        f"Gene '{marker}' not found",
                        ha='center', va='center'
                    )
                    axes[plot_idx].set_title(f"{cell_type}: {marker}")
                    axes[plot_idx].axis('off')
            
            plot_idx += 1
    
    # Hide unused axes
    for i in range(plot_idx, len(axes)):
        axes[i].axis('off')
    
    # Add legend for cell types
    if show_legend:
        handles = []
        for cell_type in cell_types:
            color = mgr.CELLS[cell_type].color or 'black'
            handles.append(mpatches.Patch(color=color, label=cell_type))
        
        fig.legend(
            handles=handles,
            loc='center right',
            bbox_to_anchor=(1.0, 0.5),
            title="Cell Types"
        )
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    plt.show()


def plot_cell_type_composition(
    adata,
    cluster_key: str,
    annotation_key: str,
    marker_config: Optional[str] = None,
    normalize: bool = True,
    plot_type: Literal['bar', 'heatmap', 'pie'] = 'bar',
    figsize: Optional[Tuple[float, float]] = None,
    color_dict: Optional[Dict[str, str]] = None,
    save: Optional[str] = None,
    sort_clusters: bool = False,
    min_fraction: float = 0.02
) -> None:
    """
    Visualize cell type composition across clusters.
    
    Args:
        adata: AnnData object
        cluster_key: Key in adata.obs for cluster information
        annotation_key: Key in adata.obs for cell type annotation
        marker_config: Path to marker configuration file (for colors)
        normalize: Whether to normalize counts to fractions
        plot_type: Type of plot to create
        figsize: Figure size
        color_dict: Dictionary mapping cell types to colors
        save: Path to save the figure
        sort_clusters: Whether to sort clusters by composition similarity
        min_fraction: Minimum fraction to include in the plot (smaller fractions are grouped as "Other")
    """
    # Check if keys exist in adata.obs
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")
    if annotation_key not in adata.obs.columns:
        raise ValueError(f"Annotation key '{annotation_key}' not found in adata.obs")
    
    # Get clusters and cell types
    clusters = adata.obs[cluster_key].cat.categories.tolist()
    cell_types = adata.obs[annotation_key].cat.categories.tolist()
    
    # Get color dictionary from marker config or adata.uns
    if color_dict is None:
        if marker_config is not None:
            mgr = Manager(marker_config)
            color_dict = {
                cell_type: mgr.CELLS[cell_type].color
                for cell_type in cell_types
                if cell_type in mgr.CELLS and mgr.CELLS[cell_type].color
            }
        elif f"{annotation_key}_colors" in adata.uns:
            color_dict = adata.uns[f"{annotation_key}_colors"]
        else:
            color_dict = {}
    
    # Calculate composition
    composition = pd.crosstab(
        adata.obs[cluster_key], 
        adata.obs[annotation_key]
    )
    
    # Normalize if requested
    if normalize:
        composition = composition.div(composition.sum(axis=1), axis=0)
    
    # Sort clusters by composition similarity if requested
    if sort_clusters and len(clusters) > 1:
        from scipy.cluster.hierarchy import linkage, dendrogram
        
        # Hierarchical clustering of clusters based on composition
        Z = linkage(composition.values, method='ward')
        d = dendrogram(Z, no_plot=True)
        
        # Reorder clusters based on dendrogram
        cluster_order = [clusters[i] for i in d['leaves']]
        composition = composition.loc[cluster_order]
    
    # Group small fractions as "Other"
    if normalize and min_fraction > 0:
        small_fractions = (composition < min_fraction).all(axis=0)
        other_types = small_fractions[small_fractions].index.tolist()
        
        if other_types:
            # Create a copy to avoid SettingWithCopyWarning
            composition_with_other = composition.copy()
            
            # Sum small fractions into "Other" category
            composition_with_other['Other'] = composition[other_types].sum(axis=1)
            
            # Drop original small fraction columns
            composition_with_other = composition_with_other.drop(columns=other_types)
            
            composition = composition_with_other
    
    # Set figure size
    if figsize is None:
        if plot_type == 'heatmap':
            figsize = (max(8, len(clusters) * 0.3), max(6, len(composition.columns) * 0.3))
        elif plot_type == 'pie':
            n_pies = len(clusters)
            ncols = min(4, n_pies)
            nrows = int(np.ceil(n_pies / ncols))
            figsize = (ncols * 4, nrows * 4)
        else:  # bar
            figsize = (max(8, len(clusters) * 0.5), 6)
    
    # Create plot
    if plot_type == 'bar':
        # Bar plot
        ax = composition.plot(
            kind='bar', 
            stacked=True, 
            figsize=figsize,
            color=[color_dict.get(ct, None) for ct in composition.columns]
        )
        
        ax.set_xlabel('Cluster')
        ax.set_ylabel('Fraction' if normalize else 'Count')
        ax.set_title('Cell Type Composition by Cluster')
        ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Add counts or percentages on bars
        if normalize:
            for i, (_, row) in enumerate(composition.iterrows()):
                cumsum = 0
                for j, val in enumerate(row):
                    if val >= 0.05:  # Only label if fraction is at least 5%
                        ax.text(
                            i, 
                            cumsum + val/2, 
                            f"{val:.1%}", 
                            ha='center', 
                            va='center',
                            fontsize=8,
                            fontweight='bold',
                            color='white' if val > 0.3 else 'black'
                        )
                    cumsum += val
    
    elif plot_type == 'heatmap':
        # Heatmap
        plt.figure(figsize=figsize)
        sns.heatmap(
            composition.T,
            annot=True,
            fmt='.0%' if normalize else '.0f',
            cmap='YlGnBu',
            linewidths=0.5,
            cbar_kws={'label': 'Fraction' if normalize else 'Count'}
        )
        
        plt.xlabel('Cluster')
        plt.ylabel('Cell Type')
        plt.title('Cell Type Composition by Cluster')
    
    elif plot_type == 'pie':
        # Multiple pie charts
        fig, axes = plt.subplots(
            int(np.ceil(len(clusters) / 4)),
            min(4, len(clusters)),
            figsize=figsize
        )
        
        if len(clusters) == 1:
            axes = np.array([axes])
        
        axes = axes.flatten()
        
        for i, cluster in enumerate(composition.index):
            if i < len(axes):
                # Get data for this cluster
                data = composition.loc[cluster]
                
                # Drop zeros
                data = data[data > 0]
                
                # Create pie chart
                wedges, texts, autotexts = axes[i].pie(
                    data,
                    labels=None,
                    autopct='%1.1f%%',
                    colors=[color_dict.get(ct, None) for ct in data.index],
                    startangle=90
                )
                
                # Customize text
                for autotext in autotexts:
                    autotext.set_fontsize(8)
                    autotext.set_fontweight('bold')
                
                axes[i].set_title(f'Cluster {cluster}')
        
        # Hide unused axes
        for i in range(len(clusters), len(axes)):
            axes[i].axis('off')
        
        # Add legend
        handles = [
            mpatches.Patch(color=color_dict.get(ct, None), label=ct)
            for ct in composition.columns
        ]
        
        fig.legend(
            handles=handles,
            loc='center right',
            bbox_to_anchor=(1.0, 0.5),
            title="Cell Types"
        )
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    plt.show()


def umap_with_annotated_clusters(
    adata,
    cluster_key: str,
    annotation_key: Optional[str] = None,
    basis: str = 'umap',
    figsize: Tuple[float, float] = (10, 8),
    text_size: int = 10,
    point_size: int = None,
    cluster_centers: Optional[pd.DataFrame] = None,
    color_by: Literal['cluster', 'annotation'] = 'cluster',
    palette: Optional[str] = None,
    show_legend: bool = True,
    save: Optional[str] = None,
    **kwargs
) -> None:
    """
    Create a UMAP or other embedding with annotated cluster labels.
    
    Args:
        adata: AnnData object
        cluster_key: Key in adata.obs for cluster information
        annotation_key: Key in adata.obs for cell type annotation
        basis: Basis for embedding visualization (e.g., 'umap', 'tsne')
        figsize: Figure size
        text_size: Size of text labels
        point_size: Size of points
        cluster_centers: DataFrame with cluster centers (if None, computed from data)
        color_by: Whether to color by cluster or annotation
        palette: Color palette to use
        show_legend: Whether to show legend
        save: Path to save the figure
        **kwargs: Additional arguments to pass to sc.pl.embedding
    """
    # Check if keys exist
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")
    
    if annotation_key is not None and annotation_key not in adata.obs.columns:
        print(f"Warning: Annotation key '{annotation_key}' not found in adata.obs. Using cluster IDs only.")
        annotation_key = None
    
    # Determine color_key based on color_by parameter
    color_key = annotation_key if color_by == 'annotation' and annotation_key is not None else cluster_key
    
    # Check if embedding exists
    embedding_key = f"X_{basis}"
    if embedding_key not in adata.obsm:
        raise ValueError(f"Embedding '{embedding_key}' not found in adata.obsm")
    
    # Create figure
    plt.figure(figsize=figsize)
    
    # Plot embedding
    sc.pl.embedding(
        adata,
        basis=basis,
        color=color_key,
        palette=palette,
        size=point_size,
        legend=False,  # We'll add our own legend
        show=False,
        **kwargs
    )
    
    # Calculate cluster centers if not provided
    if cluster_centers is None:
        cluster_centers = pd.DataFrame(
            index=adata.obs[cluster_key].cat.categories,
            columns=['x', 'y']
        )
        
        for cluster in cluster_centers.index:
            mask = adata.obs[cluster_key] == cluster
            x_vals = adata.obsm[embedding_key][mask, 0]
            y_vals = adata.obsm[embedding_key][mask, 1]
            
            cluster_centers.loc[cluster, 'x'] = np.median(x_vals)
            cluster_centers.loc[cluster, 'y'] = np.median(y_vals)
    
    # Prepare cluster labels
    texts = []
    for cluster in cluster_centers.index:
        x, y = cluster_centers.loc[cluster, ['x', 'y']]
        
        # Create label (cluster ID or annotation)
        if annotation_key is not None and color_by == 'annotation':
            # Use most common annotation for this cluster
            cluster_mask = adata.obs[cluster_key] == cluster
            annotations = adata.obs.loc[cluster_mask, annotation_key]
            most_common = annotations.value_counts().idxmax()
            label = f"{cluster}: {most_common}"
        else:
            label = str(cluster)
        
        # Add text
        texts.append(plt.text(
            x, y, label,
            fontsize=text_size,
            ha='center',
            va='center',
            fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1)
        ))
    
    # Adjust text positions to avoid overlap
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))
    
    # Add legend if requested
    if show_legend:
        if color_key in adata.obs:
            handles = []
            categories = adata.obs[color_key].cat.categories
            
            # Get colors
            if f"{color_key}_colors" in adata.uns:
                colors = adata.uns[f"{color_key}_colors"]
                if isinstance(colors, dict):
                    # Dictionary mapping categories to colors
                    colors = [colors.get(cat, 'gray') for cat in categories]
                elif len(colors) < len(categories):
                    # Not enough colors, use default colormap
                    cmap = plt.get_cmap(palette or 'tab20')
                    colors = [cmap(i / len(categories)) for i in range(len(categories))]
            else:
                # Use default colormap
                cmap = plt.get_cmap(palette or 'tab20')
                colors = [cmap(i / len(categories)) for i in range(len(categories))]
            
            # Create legend handles
            for i, category in enumerate(categories):
                if i < len(colors):
                    color = colors[i]
                    handles.append(mpatches.Patch(color=color, label=category))
            
            plt.legend(
                handles=handles,
                loc='center right',
                bbox_to_anchor=(1.2, 0.5),
                title=color_key
            )
    
    plt.title(f"{basis.upper()} with Cluster Labels")
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    plt.show()


def plot_top_genes_by_cluster(
    adata,
    cluster_key: str,
    gene_key: str = 'rank_genes_groups',
    n_genes: int = 5,
    figsize: Tuple[float, float] = (12, 10),
    cmap: str = 'viridis',
    layer: Optional[str] = None,
    save: Optional[str] = None,
) -> None:
    """
    Visualize top marker genes for each cluster on the embedding.
    
    Args:
        adata: AnnData object
        cluster_key: Key in adata.obs for cluster information
        gene_key: Key in adata.uns for rank_genes_groups results
        n_genes: Number of top genes to show per cluster
        figsize: Figure size
        cmap: Colormap for gene expression
        layer: Layer to use for expression values
        save: Path to save the figure
    """
    # Check if keys exist
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")
    
    if gene_key not in adata.uns:
        raise ValueError(f"Gene key '{gene_key}' not found in adata.uns")
    
    # Get clusters
    clusters = adata.obs[cluster_key].cat.categories.tolist()
    
    # Get top genes for each cluster
    top_genes = {}
    for i, cluster in enumerate(clusters):
        if 'names' in adata.uns[gene_key]:
            cluster_genes = adata.uns[gene_key]['names'][str(i)][:n_genes]
            top_genes[cluster] = cluster_genes
    
    # Calculate number of rows and columns
    n_clusters = len(top_genes)
    n_cols = min(3, n_clusters)
    n_rows = int(np.ceil(n_clusters / n_cols))
    
    # Create figure
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_rows * n_cols == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # Plot top genes for each cluster
    for i, (cluster, genes) in enumerate(top_genes.items()):
        if i < len(axes):
            # Subset data for this cluster
            cluster_data = adata[adata.obs[cluster_key] == cluster]
            
            # Create expression matrix for these genes
            if layer is not None and layer in cluster_data.layers:
                X = cluster_data.layers[layer]
            else:
                X = cluster_data.X
            
            # Calculate mean expression for these genes
            gene_indices = [cluster_data.var_names.get_loc(gene) for gene in genes if gene in cluster_data.var_names]
            
            if isinstance(X, np.ndarray):
                gene_expr = np.mean(X[:, gene_indices], axis=1)
            else:  # sparse matrix
                gene_expr = np.array(X[:, gene_indices].mean(axis=1)).flatten()
            
            # Plot UMAP with expression
            sc.pl.umap(
                cluster_data,
                color=genes,
                ncols=1,
                ax=axes[i],
                show=False,
                title=f"Cluster {cluster}",
                cmap=cmap,
                size=50
            )
    
    # Hide unused axes
    for i in range(len(top_genes), len(axes)):
        axes[i].axis('off')
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    plt.show()