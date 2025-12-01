"""
PyMonocle3: Python implementation of Monocle3 for single-cell RNA-seq analysis
Author: Python Port of cole-trapnell-lab/monocle3
License: MIT
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse import csr_matrix, csc_matrix, issparse
from scipy.spatial.distance import pdist, squareform
from scipy.stats import ranksums
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, Union, List, Dict, Tuple
import warnings
import h5py
import os
import pickle
from dataclasses import dataclass, field
from pathlib import Path
import json
from datetime import datetime
import umap
import igraph as ig
from collections import defaultdict


# ============================================================================
# Core Data Structures
# ============================================================================

@dataclass
class CellDataSet:
    """
    Main class for storing single-cell data (CDS object)
    
    Attributes:
        expression_data: Gene expression matrix (genes x cells)
        cell_metadata: Cell annotations (DataFrame)
        gene_metadata: Gene/feature annotations (DataFrame)
        reducedDims: Dictionary of dimensionality reduction results
        clusters: Cell cluster assignments
        partitions: Cell partition assignments
        principal_graph: Principal graph for trajectory analysis
        matrix_class: Type of matrix storage ('memory' or 'BPCells')
    """
    expression_data: Union[np.ndarray, sp.spmatrix]
    cell_metadata: pd.DataFrame
    gene_metadata: pd.DataFrame
    reducedDims: Dict[str, np.ndarray] = field(default_factory=dict)
    clusters: Optional[pd.Series] = None
    partitions: Optional[pd.Series] = None
    principal_graph: Optional[Dict] = None
    matrix_class: str = 'memory'
    _counts_row_order: Optional[str] = None
    _matrix_path: Optional[str] = None
    
    def __post_init__(self):
        """Validate and initialize CDS object"""
        # Ensure expression data is sparse for efficiency
        if not issparse(self.expression_data):
            if isinstance(self.expression_data, np.ndarray):
                self.expression_data = csr_matrix(self.expression_data)
        
        # Validate dimensions
        n_genes, n_cells = self.expression_data.shape
        assert len(self.cell_metadata) == n_cells, "Cell metadata length mismatch"
        assert len(self.gene_metadata) == n_genes, "Gene metadata length mismatch"
    
    @property
    def counts(self):
        """Get counts matrix"""
        return self.expression_data
    
    @counts.setter
    def counts(self, value, bpcells_warn=True):
        """Set counts matrix with warning for BPCells"""
        if self.matrix_class == 'BPCells' and bpcells_warn:
            warnings.warn(
                "Setting counts matrix for BPCells storage. "
                "Remember to update counts_row_order matrix using set_cds_row_order_matrix()"
            )
        self.expression_data = value
    
    @property
    def n_cells(self):
        return self.expression_data.shape[1]
    
    @property
    def n_genes(self):
        return self.expression_data.shape[0]


# ============================================================================
# Data Loading Functions
# ============================================================================

def new_cell_data_set(
    expression_data: Union[np.ndarray, sp.spmatrix],
    cell_metadata: pd.DataFrame,
    gene_metadata: pd.DataFrame,
    matrix_control: Optional[Dict] = None
) -> CellDataSet:
    """
    Create a new CellDataSet object
    
    Parameters:
        expression_data: Gene expression matrix (genes x cells)
        cell_metadata: Cell annotations
        gene_metadata: Gene annotations
        matrix_control: Dictionary with matrix storage options
            - matrix_class: 'memory' (default) or 'BPCells'
            - matrix_path: Path for BPCells storage
    
    Returns:
        CellDataSet object
    """
    matrix_control = matrix_control or {}
    matrix_class = matrix_control.get('matrix_class', 'memory')
    matrix_path = matrix_control.get('matrix_path', None)
    
    cds = CellDataSet(
        expression_data=expression_data,
        cell_metadata=cell_metadata,
        gene_metadata=gene_metadata,
        matrix_class=matrix_class,
        _matrix_path=matrix_path
    )
    
    return cds


def load_mm_data(
    mat_path: str,
    feature_anno_path: str,
    cell_anno_path: str,
    matrix_control: Optional[Dict] = None,
    sep: str = '\t'
) -> CellDataSet:
    """
    Load data from MatrixMarket format
    
    Parameters:
        mat_path: Path to .mtx file
        feature_anno_path: Path to feature/gene annotations
        cell_anno_path: Path to cell annotations
        matrix_control: Matrix storage options
        sep: Separator for annotation files
    
    Returns:
        CellDataSet object
    """
    from scipy.io import mmread
    
    # Load expression matrix
    expression_data = mmread(mat_path).T.tocsr()  # Transpose to genes x cells
    
    # Load annotations
    gene_metadata = pd.read_csv(feature_anno_path, sep=sep)
    cell_metadata = pd.read_csv(cell_anno_path, sep=sep)
    
    return new_cell_data_set(expression_data, cell_metadata, gene_metadata, matrix_control)


def load_cellranger_data(
    data_dir: str,
    genome: Optional[str] = None,
    matrix_control: Optional[Dict] = None
) -> CellDataSet:
    """
    Load data from CellRanger output directory
    
    Parameters:
        data_dir: Path to CellRanger output directory
        genome: Genome name (if multiple genomes present)
        matrix_control: Matrix storage options
    
    Returns:
        CellDataSet object
    """
    from scipy.io import mmread
    
    data_path = Path(data_dir)
    
    # Handle different CellRanger versions
    if (data_path / 'filtered_feature_bc_matrix').exists():
        base_path = data_path / 'filtered_feature_bc_matrix'
        if genome:
            base_path = base_path / genome
    else:
        base_path = data_path
    
    # Load files
    matrix_file = base_path / 'matrix.mtx.gz' if (base_path / 'matrix.mtx.gz').exists() else base_path / 'matrix.mtx'
    features_file = base_path / 'features.tsv.gz' if (base_path / 'features.tsv.gz').exists() else base_path / 'genes.tsv'
    barcodes_file = base_path / 'barcodes.tsv.gz' if (base_path / 'barcodes.tsv.gz').exists() else base_path / 'barcodes.tsv'
    
    expression_data = mmread(matrix_file).T.tocsr()
    gene_metadata = pd.read_csv(features_file, sep='\t', header=None, names=['gene_id', 'gene_short_name', 'feature_type'])
    cell_metadata = pd.read_csv(barcodes_file, sep='\t', header=None, names=['barcode'])
    cell_metadata.index = cell_metadata['barcode']
    
    return new_cell_data_set(expression_data, cell_metadata, gene_metadata, matrix_control)


# ============================================================================
# Preprocessing Functions
# ============================================================================

def preprocess_cds(
    cds: CellDataSet,
    method: str = 'PCA',
    num_dim: int = 50,
    norm_method: str = 'log',
    use_genes: Optional[List[str]] = None,
    pseudo_count: float = 1.0,
    scaling: bool = True,
    verbose: bool = True
) -> CellDataSet:
    """
    Preprocess CDS: normalize and dimensionality reduction
    
    Parameters:
        cds: CellDataSet object
        method: Dimensionality reduction method ('PCA', 'LSI')
        num_dim: Number of dimensions to compute
        norm_method: Normalization method ('log', 'size_only')
        use_genes: Subset of genes to use
        pseudo_count: Pseudocount for log normalization
        scaling: Whether to scale features
        verbose: Print progress messages
    
    Returns:
        Updated CellDataSet object
    """
    if verbose:
        print(f"Preprocessing data using {method}...")
    
    # Get expression data
    expr_data = cds.expression_data.copy()
    
    # Subset genes if specified
    if use_genes is not None:
        gene_idx = cds.gene_metadata.index.isin(use_genes)
        expr_data = expr_data[gene_idx, :]
    
    # Normalize
    if norm_method == 'log':
        # Size factor normalization
        size_factors = np.array(expr_data.sum(axis=0)).flatten()
        size_factors = size_factors / np.median(size_factors)
        
        # Normalize and log transform
        if issparse(expr_data):
            expr_data = expr_data.multiply(1 / size_factors)
            expr_data.data = np.log1p(expr_data.data)
        else:
            expr_data = expr_data / size_factors
            expr_data = np.log1p(expr_data)
    
    # Convert to dense for PCA
    if issparse(expr_data):
        expr_data_dense = expr_data.toarray()
    else:
        expr_data_dense = expr_data
    
    # Perform dimensionality reduction
    if method == 'PCA':
        if scaling:
            # Center and scale
            means = np.mean(expr_data_dense, axis=1, keepdims=True)
            stds = np.std(expr_data_dense, axis=1, keepdims=True)
            stds[stds == 0] = 1
            expr_data_dense = (expr_data_dense - means) / stds
        
        pca = PCA(n_components=num_dim)
        reduced_dims = pca.fit_transform(expr_data_dense.T)
        
        cds.reducedDims['PCA'] = reduced_dims
        
        if verbose:
            var_explained = pca.explained_variance_ratio_
            print(f"PC1-10 explain {var_explained[:10].sum():.2%} of variance")
    
    elif method == 'LSI':
        # LSI (Latent Semantic Indexing) for ATAC-seq data
        # TF-IDF transformation
        idf = np.log1p(expr_data_dense.shape[1] / (np.sum(expr_data_dense > 0, axis=1) + 1))
        tf_idf = expr_data_dense * idf[:, np.newaxis]
        
        # SVD
        from sklearn.decomposition import TruncatedSVD
        svd = TruncatedSVD(n_components=num_dim)
        reduced_dims = svd.fit_transform(tf_idf.T)
        
        cds.reducedDims['LSI'] = reduced_dims
    
    return cds


def detect_genes(
    cds: CellDataSet,
    min_expr: float = 0.1
) -> CellDataSet:
    """
    Detect expressed genes in each cell
    
    Parameters:
        cds: CellDataSet object
        min_expr: Minimum expression threshold
    
    Returns:
        Updated CellDataSet with num_genes_expressed in cell_metadata
    """
    expr_data = cds.expression_data
    
    if issparse(expr_data):
        genes_detected = (expr_data > min_expr).sum(axis=0)
        genes_detected = np.array(genes_detected).flatten()
    else:
        genes_detected = np.sum(expr_data > min_expr, axis=0)
    
    cds.cell_metadata['num_genes_expressed'] = genes_detected
    
    return cds


# ============================================================================
# Dimensionality Reduction
# ============================================================================

def reduce_dimension(
    cds: CellDataSet,
    max_components: int = 2,
    reduction_method: str = 'UMAP',
    preprocess_method: str = 'PCA',
    verbose: bool = True,
    **kwargs
) -> CellDataSet:
    """
    Reduce dimensions for visualization (UMAP, tSNE)
    
    Parameters:
        cds: CellDataSet object
        max_components: Number of components (usually 2 for visualization)
        reduction_method: 'UMAP' or 'tSNE'
        preprocess_method: Which preprocessing to use as input
        verbose: Print progress
        **kwargs: Additional arguments for UMAP/tSNE
    
    Returns:
        Updated CellDataSet object
    """
    if verbose:
        print(f"Reducing dimensions using {reduction_method}...")
    
    # Get preprocessed data
    if preprocess_method not in cds.reducedDims:
        raise ValueError(f"Must run preprocess_cds with method={preprocess_method} first")
    
    input_data = cds.reducedDims[preprocess_method]
    
    if reduction_method == 'UMAP':
        # Set default UMAP parameters
        umap_params = {
            'n_components': max_components,
            'n_neighbors': 15,
            'min_dist': 0.1,
            'metric': 'euclidean',
            'random_state': 42
        }
        umap_params.update(kwargs)
        
        reducer = umap.UMAP(**umap_params)
        embedding = reducer.fit_transform(input_data)
        
        cds.reducedDims['UMAP'] = embedding
    
    elif reduction_method == 'tSNE':
        from sklearn.manifold import TSNE
        
        tsne_params = {
            'n_components': max_components,
            'perplexity': 30,
            'random_state': 42
        }
        tsne_params.update(kwargs)
        
        reducer = TSNE(**tsne_params)
        embedding = reducer.fit_transform(input_data)
        
        cds.reducedDims['tSNE'] = embedding
    
    return cds


# ============================================================================
# Clustering
# ============================================================================

def cluster_cells(
    cds: CellDataSet,
    reduction_method: str = 'UMAP',
    k: int = 20,
    cluster_method: str = 'louvain',
    resolution: Optional[float] = None,
    random_seed: int = 42,
    verbose: bool = True
) -> CellDataSet:
    """
    Cluster cells using graph-based clustering
    
    Parameters:
        cds: CellDataSet object
        reduction_method: Which reduction to use for clustering
        k: Number of nearest neighbors
        cluster_method: 'louvain' or 'leiden'
        resolution: Resolution parameter for clustering
        random_seed: Random seed
        verbose: Print progress
    
    Returns:
        Updated CellDataSet with clusters
    """
    if verbose:
        print(f"Clustering cells using {cluster_method}...")
    
    # Get reduced dimensions
    if reduction_method not in cds.reducedDims:
        raise ValueError(f"Reduction method {reduction_method} not found. Run reduce_dimension first.")
    
    reduced_data = cds.reducedDims[reduction_method]
    
    # Build k-NN graph
    nbrs = NearestNeighbors(n_neighbors=k, metric='euclidean').fit(reduced_data)
    distances, indices = nbrs.kneighbors(reduced_data)
    
    # Create igraph
    edges = []
    weights = []
    for i in range(len(indices)):
        for j, neighbor in enumerate(indices[i]):
            if i < neighbor:  # Avoid duplicate edges
                edges.append((i, neighbor))
                weights.append(1 / (1 + distances[i, j]))
    
    g = ig.Graph(n=len(reduced_data), edges=edges)
    g.es['weight'] = weights
    
    # Cluster
    np.random.seed(random_seed)
    if cluster_method == 'louvain':
        if resolution is None:
            resolution = 1.0
        clusters = g.community_multilevel(weights='weight', return_levels=False)
    elif cluster_method == 'leiden':
        if resolution is None:
            resolution = 1.0
        clusters = g.community_leiden(weights='weight', resolution_parameter=resolution)
    else:
        raise ValueError(f"Unknown cluster method: {cluster_method}")
    
    # Store cluster assignments
    cluster_labels = np.array(clusters.membership) + 1  # 1-indexed
    cds.clusters = pd.Series(cluster_labels, index=cds.cell_metadata.index, name='cluster')
    cds.cell_metadata['cluster'] = cluster_labels
    
    if verbose:
        print(f"Found {len(np.unique(cluster_labels))} clusters")
    
    return cds


def partition_cells(
    cds: CellDataSet,
    reduction_method: str = 'UMAP',
    knn: int = 25,
    verbose: bool = True
) -> CellDataSet:
    """
    Partition cells into disconnected components
    
    Parameters:
        cds: CellDataSet object
        reduction_method: Which reduction to use
        knn: Number of nearest neighbors
        verbose: Print progress
    
    Returns:
        Updated CellDataSet with partitions
    """
    if verbose:
        print("Partitioning cells...")
    
    reduced_data = cds.reducedDims[reduction_method]
    
    # Build k-NN graph
    nbrs = NearestNeighbors(n_neighbors=knn, metric='euclidean').fit(reduced_data)
    distances, indices = nbrs.kneighbors(reduced_data)
    
    # Create graph
    edges = []
    for i in range(len(indices)):
        for neighbor in indices[i]:
            edges.append((i, neighbor))
    
    g = ig.Graph(n=len(reduced_data), edges=edges)
    
    # Find connected components
    components = g.components()
    partition_labels = np.array(components.membership) + 1
    
    cds.partitions = pd.Series(partition_labels, index=cds.cell_metadata.index, name='partition')
    cds.cell_metadata['partition'] = partition_labels
    
    if verbose:
        print(f"Found {len(np.unique(partition_labels))} partitions")
    
    return cds


# ============================================================================
# Trajectory Analysis
# ============================================================================

def learn_graph(
    cds: CellDataSet,
    use_partition: bool = True,
    close_loop: bool = True,
    learn_graph_control: Optional[Dict] = None,
    verbose: bool = True
) -> CellDataSet:
    """
    Learn principal graph for trajectory inference
    
    Parameters:
        cds: CellDataSet object
        use_partition: Whether to learn separate graphs per partition
        close_loop: Whether to close loops in the graph
        learn_graph_control: Additional control parameters
        verbose: Print progress
    
    Returns:
        Updated CellDataSet with principal graph
    """
    if verbose:
        print("Learning principal graph...")
    
    learn_graph_control = learn_graph_control or {}
    
    # Use UMAP coordinates
    if 'UMAP' not in cds.reducedDims:
        raise ValueError("Must run reduce_dimension with UMAP first")
    
    coords = cds.reducedDims['UMAP']
    
    # Simple implementation using minimum spanning tree
    from scipy.spatial.distance import pdist, squareform
    from scipy.sparse.csgraph import minimum_spanning_tree
    
    if use_partition and cds.partitions is not None:
        # Learn graph for each partition
        graphs = {}
        for partition in np.unique(cds.partitions):
            mask = cds.partitions == partition
            partition_coords = coords[mask]
            
            # Compute distance matrix
            dist_matrix = squareform(pdist(partition_coords))
            
            # Get MST
            mst = minimum_spanning_tree(dist_matrix)
            
            graphs[partition] = {
                'coords': partition_coords,
                'mst': mst,
                'cell_indices': np.where(mask)[0]
            }
        
        cds.principal_graph = graphs
    else:
        # Single graph
        dist_matrix = squareform(pdist(coords))
        mst = minimum_spanning_tree(dist_matrix)
        
        cds.principal_graph = {
            'coords': coords,
            'mst': mst,
            'cell_indices': np.arange(len(coords))
        }
    
    return cds


def order_cells(
    cds: CellDataSet,
    root_cells: Optional[List[str]] = None,
    reduction_method: str = 'UMAP',
    verbose: bool = True
) -> CellDataSet:
    """
    Order cells in pseudotime along trajectory
    
    Parameters:
        cds: CellDataSet object
        root_cells: Cell IDs to use as trajectory root
        reduction_method: Which reduction to use
        verbose: Print progress
    
    Returns:
        Updated CellDataSet with pseudotime
    """
    if verbose:
        print("Ordering cells in pseudotime...")
    
    if cds.principal_graph is None:
        raise ValueError("Must run learn_graph first")
    
    coords = cds.reducedDims[reduction_method]
    
    # Simple pseudotime: distance from root along graph
    if root_cells is not None:
        root_idx = [cds.cell_metadata.index.get_loc(c) for c in root_cells]
        root_coord = coords[root_idx].mean(axis=0)
    else:
        # Use cell closest to origin
        root_idx = np.argmin(np.linalg.norm(coords, axis=1))
        root_coord = coords[root_idx]
    
    # Compute distances from root
    distances = np.linalg.norm(coords - root_coord, axis=1)
    
    cds.cell_metadata['pseudotime'] = distances
    
    if verbose:
        print(f"Pseudotime range: {distances.min():.2f} to {distances.max():.2f}")
    
    return cds


# ============================================================================
# Differential Expression
# ============================================================================

def graph_test(
    cds: CellDataSet,
    neighbor_graph: str = 'principal_graph',
    reduction_method: str = 'UMAP',
    k: int = 25,
    method: str = 'Moran_I',
    alternative: str = 'greater',
    verbose: bool = True
) -> pd.DataFrame:
    """
    Test genes for differential expression along trajectory
    
    Parameters:
        cds: CellDataSet object
        neighbor_graph: Which graph to use
        reduction_method: Reduction method
        k: Number of neighbors
        method: Test method ('Moran_I')
        alternative: Alternative hypothesis
        verbose: Print progress
    
    Returns:
        DataFrame with test results
    """
    if verbose:
        print(f"Testing genes using {method}...")
    
    expr_data = cds.expression_data
    coords = cds.reducedDims[reduction_method]
    
    # Build neighbor graph
    nbrs = NearestNeighbors(n_neighbors=k, metric='euclidean').fit(coords)
    distances, indices = nbrs.kneighbors(coords)
    
    # Compute weights (inverse distance)
    weights = 1 / (distances + 1e-10)
    weights = weights / weights.sum(axis=1, keepdims=True)
    
    results = []
    
    # Test each gene
    n_genes = cds.n_genes
    for i in range(n_genes):
        if i % 1000 == 0 and verbose:
            print(f"Testing gene {i}/{n_genes}")
        
        gene_expr = expr_data[i, :].toarray().flatten() if issparse(expr_data) else expr_data[i, :]
        
        # Moran's I statistic
        n = len(gene_expr)
        mean_expr = gene_expr.mean()
        
        numerator = 0
        denominator = 0
        
        for j in range(n):
            for l, neighbor in enumerate(indices[j]):
                w = weights[j, l]
                numerator += w * (gene_expr[j] - mean_expr) * (gene_expr[neighbor] - mean_expr)
            
            denominator += (gene_expr[j] - mean_expr) ** 2
        
        W = weights.sum()
        morans_i = (n / W) * (numerator / denominator) if denominator > 0 else 0
        
        # Simple p-value (would need proper null distribution)
        # For now, use normal approximation
        from scipy.stats import norm
        expected_i = -1 / (n - 1)
        std_i = 1 / np.sqrt(n)  # Simplified
        z_score = (morans_i - expected_i) / std_i
        
        if alternative == 'greater':
            p_value = 1 - norm.cdf(z_score)
        else:
            p_value = norm.cdf(z_score)
        
        results.append({
            'gene_id': cds.gene_metadata.index[i],
            'gene_short_name': cds.gene_metadata.iloc[i].get('gene_short_name', ''),
            'morans_i': morans_i,
            'morans_test_statistic': z_score,
            'p_value': p_value
        })
    
    results_df = pd.DataFrame(results)
    
    # Multiple testing correction (Benjamini-Hochberg)
    from statsmodels.stats.multitest import multipletests
    _, results_df['q_value'], _, _ = multipletests(results_df['p_value'], method='fdr_bh')
    
    results_df = results_df.sort_values('q_value')
    
    return results_df


def top_markers(
    cds: CellDataSet,
    group_cells_by: str = 'cluster',
    genes_to_test_per_group: int = 25,
    reduction_method: str = 'UMAP',
    verbose: bool = True
) -> pd.DataFrame:
    """
    Find top marker genes for each group
    
    Parameters:
        cds: CellDataSet object
        group_cells_by: Column in cell_metadata to group by
        genes_to_test_per_group: Number of top genes per group
        reduction_method: Reduction method
        verbose: Print progress
    
    Returns:
        DataFrame with top markers
    """
    if verbose:
        print(f"Finding top markers for {group_cells_by}...")
    
    if group_cells_by not in cds.cell_metadata.columns:
        raise ValueError(f"{group_cells_by} not found in cell_metadata")
    
    expr_data = cds.expression_data
    groups = cds.cell_metadata[group_cells_by]
    
    all_markers = []
    
    for group in np.unique(groups):
        group_mask = groups == group
        other_mask = ~group_mask
        
        # Compute mean expression in vs out of group
        if issparse(expr_data):
            group_mean = np.array(expr_data[:, group_mask].mean(axis=1)).flatten()
            other_mean = np.array(expr_data[:, other_mask].mean(axis=1)).flatten()
        else:
            group_mean = expr_data[:, group_mask].mean(axis=1)
            other_mean = expr_data[:, other_mask].mean(axis=1)
        
        # Log fold change
        log_fc = np.log2((group_mean + 1) / (other_mean + 1))
        
        # Rank and select top genes
        top_indices = np.argsort(-log_fc)[:genes_to_test_per_group]
        
        for idx in top_indices:
            # Wilcoxon test
            if issparse(expr_data):
                group_expr = expr_data[idx, group_mask].toarray().flatten()
                other_expr = expr_data[idx, other_mask].toarray().flatten()
            else:
                group_expr = expr_data[idx, group_mask]
                other_expr = expr_data[idx, other_mask]
            
            stat, p_value = ranksums(group_expr, other_expr)
            
            all_markers.append({
                'gene_id': cds.gene_metadata.index[idx],
                'gene_short_name': cds.gene_metadata.iloc[idx].get('gene_short_name', ''),
                'cell_group': group,
                'mean_expression': group_mean[idx],
                'log2_fold_change': log_fc[idx],
                'p_value': p_value
            })
    
    markers_df = pd.DataFrame(all_markers)
    
    # FDR correction
    from statsmodels.stats.multitest import multipletests
    _, markers_df['q_value'], _, _ = multipletests(markers_df['p_value'], method='fdr_bh')
    
    markers_df = markers_df.sort_values(['cell_group', 'q_value'])
    
    return markers_df


# ============================================================================
# Visualization Functions
# ============================================================================

def plot_cells(
    cds: CellDataSet,
    x: int = 1,
    y: int = 2,
    reduction_method: str = 'UMAP',
    color_cells_by: Optional[str] = None,
    genes: Optional[List[str]] = None,
    show_trajectory_graph: bool = False,
    label_cell_groups: bool = True,
    cell_size: float = 0.5,
    alpha: float = 0.5,
    figsize: Tuple[int, int] = (10, 8)
) -> plt.Figure:
    """
    Plot cells in reduced dimensions
    
    Parameters:
        cds: CellDataSet object
        x, y: Which dimensions to plot
        reduction_method: Which reduction to use
        color_cells_by: Column to color cells by
        genes: Genes to plot expression
        show_trajectory_graph: Show principal graph
        label_cell_groups: Label groups
        cell_size: Size of points
        alpha: Transparency
        figsize: Figure size
    
    Returns:
        Matplotlib figure
    """
    if reduction_method not in cds.reducedDims:
        raise ValueError(f"Reduction method {reduction_method} not found")
    
    coords = cds.reducedDims[reduction_method]
    
    # Handle multiple genes
    n_genes = len(genes) if genes else 1
    n_plots = n_genes if genes else 1
    
    fig, axes = plt.subplots(1, n_plots, figsize=(figsize[0] * n_plots, figsize[1]))
    if n_plots == 1:
        axes = [axes]
    
    for plot_idx in range(n_plots):
        ax = axes[plot_idx]
        
        # Determine colors
        if genes:
            gene = genes[plot_idx]
            gene_idx = cds.gene_metadata.index.get_loc(gene) if gene in cds.gene_metadata.index else None
            
            if gene_idx is not None:
                if issparse(cds.expression_data):
                    colors = cds.expression_data[gene_idx, :].toarray().flatten()
                else:
                    colors = cds.expression_data[gene_idx, :]
                cmap = 'viridis'
                label = f'{gene} expression'
            else:
                colors = 'gray'
                cmap = None
                label = f'{gene} (not found)'
        
        elif color_cells_by and color_cells_by in cds.cell_metadata.columns:
            color_col = cds.cell_metadata[color_cells_by]
            
            if pd.api.types.is_numeric_dtype(color_col):
                colors = color_col.values
                cmap = 'viridis'
            else:
                # Categorical
                unique_vals = color_col.unique()
                color_map = {val: i for i, val in enumerate(unique_vals)}
                colors = color_col.map(color_map).values
                cmap = 'tab20'
            
            label = color_cells_by
        else:
            colors = 'steelblue'
            cmap = None
            label = None
        
        # Plot
        scatter = ax.scatter(
            coords[:, x-1], coords[:, y-1],
            c=colors, s=cell_size, alpha=alpha,
            cmap=cmap
        )
        
        # Add trajectory graph
        if show_trajectory_graph and cds.principal_graph is not None:
            if isinstance(cds.principal_graph, dict) and 'mst' in cds.principal_graph:
                mst = cds.principal_graph['mst'].tocoo()
                for i, j, v in zip(mst.row, mst.col, mst.data):
                    ax.plot(
                        [coords[i, x-1], coords[j, x-1]],
                        [coords[i, y-1], coords[j, y-1]],
                        'k-', alpha=0.3, linewidth=1
                    )
        
        # Labels
        ax.set_xlabel(f'{reduction_method}_{x}')
        ax.set_ylabel(f'{reduction_method}_{y}')
        
        if label:
            ax.set_title(label)
        
        if cmap and isinstance(colors, np.ndarray):
            plt.colorbar(scatter, ax=ax)
        
        # Label groups
        if label_cell_groups and color_cells_by and color_cells_by in cds.cell_metadata.columns:
            if not pd.api.types.is_numeric_dtype(cds.cell_metadata[color_cells_by]):
                for group in cds.cell_metadata[color_cells_by].unique():
                    mask = cds.cell_metadata[color_cells_by] == group
                    center = coords[mask].mean(axis=0)
                    ax.text(center[x-1], center[y-1], str(group), 
                           fontsize=12, fontweight='bold',
                           ha='center', va='center')
    
    plt.tight_layout()
    return fig


def plot_genes_by_group(
    cds: CellDataSet,
    genes: List[str],
    group_cells_by: str = 'cluster',
    ordering_type: str = 'cluster_row_col',
    max_size: int = 3,
    normalize: bool = True,
    figsize: Tuple[int, int] = (12, 8)
) -> plt.Figure:
    """
    Plot heatmap of gene expression by group
    
    Parameters:
        cds: CellDataSet object
        genes: List of genes to plot
        group_cells_by: How to group cells
        ordering_type: How to order rows/columns
        max_size: Maximum dot size
        normalize: Z-score normalize expression
        figsize: Figure size
    
    Returns:
        Matplotlib figure
    """
    # Get expression data for genes
    gene_indices = [cds.gene_metadata.index.get_loc(g) for g in genes if g in cds.gene_metadata.index]
    
    if not gene_indices:
        raise ValueError("No genes found")
    
    if issparse(cds.expression_data):
        expr_matrix = cds.expression_data[gene_indices, :].toarray()
    else:
        expr_matrix = cds.expression_data[gene_indices, :]
    
    # Group cells
    groups = cds.cell_metadata[group_cells_by]
    unique_groups = sorted(groups.unique())
    
    # Compute mean expression per group
    mean_expr = np.zeros((len(gene_indices), len(unique_groups)))
    pct_expr = np.zeros((len(gene_indices), len(unique_groups)))
    
    for i, group in enumerate(unique_groups):
        mask = groups == group
        mean_expr[:, i] = expr_matrix[:, mask].mean(axis=1)
        pct_expr[:, i] = (expr_matrix[:, mask] > 0).mean(axis=1) * 100
    
    # Normalize if requested
    if normalize:
        mean_expr = (mean_expr - mean_expr.mean(axis=1, keepdims=True)) / (mean_expr.std(axis=1, keepdims=True) + 1e-10)
    
    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create dot plot
    for i in range(len(gene_indices)):
        for j in range(len(unique_groups)):
            size = (pct_expr[i, j] / 100) * max_size * 100
            color = mean_expr[i, j]
            ax.scatter(j, i, s=size, c=[color], cmap='RdBu_r', 
                      vmin=-2, vmax=2, edgecolors='black', linewidth=0.5)
    
    # Labels
    ax.set_xticks(range(len(unique_groups)))
    ax.set_xticklabels(unique_groups, rotation=45, ha='right')
    ax.set_yticks(range(len(gene_indices)))
    ax.set_yticklabels([genes[i] for i in range(len(gene_indices))])
    ax.set_xlabel(group_cells_by)
    ax.set_ylabel('Genes')
    
    plt.colorbar(plt.cm.ScalarMappable(cmap='RdBu_r', norm=plt.Normalize(vmin=-2, vmax=2)), 
                 ax=ax, label='Mean expression (z-score)')
    
    plt.tight_layout()
    return fig


def plot_pseudotime_heatmap(
    cds: CellDataSet,
    genes: List[str],
    num_bins: int = 100,
    scale_max: float = 3,
    figsize: Tuple[int, int] = (10, 8)
) -> plt.Figure:
    """
    Plot gene expression along pseudotime
    
    Parameters:
        cds: CellDataSet object
        genes: Genes to plot
        num_bins: Number of pseudotime bins
        scale_max: Maximum z-score for colormap
        figsize: Figure size
    
    Returns:
        Matplotlib figure
    """
    if 'pseudotime' not in cds.cell_metadata.columns:
        raise ValueError("Must run order_cells first")
    
    # Get gene indices
    gene_indices = [cds.gene_metadata.index.get_loc(g) for g in genes if g in cds.gene_metadata.index]
    
    if not gene_indices:
        raise ValueError("No genes found")
    
    # Get expression
    if issparse(cds.expression_data):
        expr_matrix = cds.expression_data[gene_indices, :].toarray()
    else:
        expr_matrix = cds.expression_data[gene_indices, :]
    
    # Get pseudotime and sort cells
    pseudotime = cds.cell_metadata['pseudotime'].values
    sort_idx = np.argsort(pseudotime)
    
    expr_matrix = expr_matrix[:, sort_idx]
    pseudotime = pseudotime[sort_idx]
    
    # Bin cells by pseudotime
    bins = np.linspace(pseudotime.min(), pseudotime.max(), num_bins + 1)
    bin_indices = np.digitize(pseudotime, bins) - 1
    
    # Compute mean expression per bin
    binned_expr = np.zeros((len(gene_indices), num_bins))
    for i in range(num_bins):
        mask = bin_indices == i
        if mask.sum() > 0:
            binned_expr[:, i] = expr_matrix[:, mask].mean(axis=1)
    
    # Z-score normalize
    binned_expr = (binned_expr - binned_expr.mean(axis=1, keepdims=True)) / (binned_expr.std(axis=1, keepdims=True) + 1e-10)
    binned_expr = np.clip(binned_expr, -scale_max, scale_max)
    
    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(binned_expr, aspect='auto', cmap='RdBu_r', 
                   vmin=-scale_max, vmax=scale_max, interpolation='nearest')
    
    ax.set_xlabel('Pseudotime')
    ax.set_ylabel('Genes')
    ax.set_yticks(range(len(gene_indices)))
    ax.set_yticklabels([genes[i] for i in range(len(gene_indices))])
    
    plt.colorbar(im, ax=ax, label='Expression (z-score)')
    plt.tight_layout()
    
    return fig


# ============================================================================
# Save/Load Functions
# ============================================================================

def save_monocle_objects(
    cds: CellDataSet,
    directory: str,
    compress: bool = True
) -> None:
    """
    Save CellDataSet to disk
    
    Parameters:
        cds: CellDataSet object
        directory: Directory to save to
        compress: Whether to compress
    """
    os.makedirs(directory, exist_ok=True)
    
    # Save expression data
    if cds.matrix_class == 'BPCells':
        # Would save BPCells matrix here
        warnings.warn("BPCells saving not fully implemented in Python version")
    
    # Save as scipy sparse matrix
    sp.save_npz(os.path.join(directory, 'expression_data.npz'), 
                cds.expression_data.tocsr())
    
    # Save metadata
    cds.cell_metadata.to_csv(os.path.join(directory, 'cell_metadata.csv'))
    cds.gene_metadata.to_csv(os.path.join(directory, 'gene_metadata.csv'))
    
    # Save reduced dims
    for name, data in cds.reducedDims.items():
        np.save(os.path.join(directory, f'reducedDims_{name}.npy'), data)
    
    # Save clusters/partitions
    if cds.clusters is not None:
        cds.clusters.to_csv(os.path.join(directory, 'clusters.csv'))
    
    if cds.partitions is not None:
        cds.partitions.to_csv(os.path.join(directory, 'partitions.csv'))
    
    # Save principal graph
    if cds.principal_graph is not None:
        with open(os.path.join(directory, 'principal_graph.pkl'), 'wb') as f:
            pickle.dump(cds.principal_graph, f)
    
    # Save metadata
    metadata = {
        'matrix_class': cds.matrix_class,
        'n_cells': cds.n_cells,
        'n_genes': cds.n_genes,
        'saved_date': datetime.now().isoformat()
    }
    
    with open(os.path.join(directory, 'metadata.json'), 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"Saved CellDataSet to {directory}")


def load_monocle_objects(directory: str) -> CellDataSet:
    """
    Load CellDataSet from disk
    
    Parameters:
        directory: Directory to load from
    
    Returns:
        CellDataSet object
    """
    # Load expression data
    expression_data = sp.load_npz(os.path.join(directory, 'expression_data.npz'))
    
    # Load metadata
    cell_metadata = pd.read_csv(os.path.join(directory, 'cell_metadata.csv'), index_col=0)
    gene_metadata = pd.read_csv(os.path.join(directory, 'gene_metadata.csv'), index_col=0)
    
    # Create CDS
    cds = CellDataSet(
        expression_data=expression_data,
        cell_metadata=cell_metadata,
        gene_metadata=gene_metadata
    )
    
    # Load reduced dims
    for file in os.listdir(directory):
        if file.startswith('reducedDims_') and file.endswith('.npy'):
            name = file.replace('reducedDims_', '').replace('.npy', '')
            cds.reducedDims[name] = np.load(os.path.join(directory, file))
    
    # Load clusters/partitions
    if os.path.exists(os.path.join(directory, 'clusters.csv')):
        cds.clusters = pd.read_csv(os.path.join(directory, 'clusters.csv'), index_col=0, squeeze=True)
        cds.cell_metadata['cluster'] = cds.clusters
    
    if os.path.exists(os.path.join(directory, 'partitions.csv')):
        cds.partitions = pd.read_csv(os.path.join(directory, 'partitions.csv'), index_col=0, squeeze=True)
        cds.cell_metadata['partition'] = cds.partitions
    
    # Load principal graph
    if os.path.exists(os.path.join(directory, 'principal_graph.pkl')):
        with open(os.path.join(directory, 'principal_graph.pkl'), 'rb') as f:
            cds.principal_graph = pickle.load(f)
    
    print(f"Loaded CellDataSet from {directory}")
    return cds


# ============================================================================
# Example Usage
# ============================================================================

def example_analysis():
    """
    Example analysis workflow
    """
    print("PyMonocle3 Example Analysis")
    print("=" * 50)
    
    # 1. Generate example data
    print("\n1. Generating example data...")
    np.random.seed(42)
    n_genes = 2000
    n_cells = 1000
    
    # Simulate expression data with trajectory structure
    pseudotime = np.random.rand(n_cells)
    expression_data = np.random.poisson(
        5 * np.exp(-2 * np.abs(np.random.randn(n_genes, 1) - pseudotime))
    )
    expression_data = csr_matrix(expression_data)
    
    gene_metadata = pd.DataFrame({
        'gene_id': [f'Gene_{i}' for i in range(n_genes)],
        'gene_short_name': [f'GENE{i}' for i in range(n_genes)]
    })
    gene_metadata.index = gene_metadata['gene_id']
    
    cell_metadata = pd.DataFrame({
        'cell_id': [f'Cell_{i}' for i in range(n_cells)],
        'batch': np.random.choice(['A', 'B'], n_cells)
    })
    cell_metadata.index = cell_metadata['cell_id']
    
    # 2. Create CDS
    print("\n2. Creating CellDataSet...")
    cds = new_cell_data_set(expression_data, cell_metadata, gene_metadata)
    print(f"Created CDS with {cds.n_genes} genes and {cds.n_cells} cells")
    
    # 3. Preprocess
    print("\n3. Preprocessing...")
    cds = preprocess_cds(cds, num_dim=30)
    
    # 4. Reduce dimensions
    print("\n4. Reducing dimensions with UMAP...")
    cds = reduce_dimension(cds, max_components=2, reduction_method='UMAP')
    
    # 5. Cluster cells
    print("\n5. Clustering cells...")
    cds = cluster_cells(cds, resolution=1.0)
    
    # 6. Learn trajectory
    print("\n6. Learning trajectory...")
    cds = learn_graph(cds)
    
    # 7. Order cells
    print("\n7. Ordering cells in pseudotime...")
    cds = order_cells(cds)
    
    # 8. Find markers
    print("\n8. Finding marker genes...")
    markers = top_markers(cds, genes_to_test_per_group=10)
    print(f"Found {len(markers)} marker genes")
    print("\nTop markers:")
    print(markers.head(10)[['gene_short_name', 'cell_group', 'log2_fold_change', 'q_value']])
    
    # 9. Visualize
    print("\n9. Creating visualizations...")
    
    # Plot by cluster
    fig1 = plot_cells(cds, color_cells_by='cluster', label_cell_groups=True)
    plt.savefig('cells_by_cluster.png', dpi=300, bbox_inches='tight')
    print("Saved: cells_by_cluster.png")
    
    # Plot by pseudotime
    fig2 = plot_cells(cds, color_cells_by='pseudotime', show_trajectory_graph=True)
    plt.savefig('cells_by_pseudotime.png', dpi=300, bbox_inches='tight')
    print("Saved: cells_by_pseudotime.png")
    
    # Plot top marker genes
    top_genes = markers.head(6)['gene_short_name'].tolist()
    fig3 = plot_genes_by_group(cds, top_genes, group_cells_by='cluster')
    plt.savefig('marker_genes_heatmap.png', dpi=300, bbox_inches='tight')
    print("Saved: marker_genes_heatmap.png")
    
    # 10. Save
    print("\n10. Saving CellDataSet...")
    save_monocle_objects(cds, 'example_cds')
    
    # 11. Load
    print("\n11. Testing load...")
    cds_loaded = load_monocle_objects('example_cds')
    
    print("\n" + "=" * 50)
    print("Analysis complete!")
    print(f"Final CDS: {cds_loaded.n_genes} genes, {cds_loaded.n_cells} cells")
    print(f"Clusters: {len(cds_loaded.clusters.unique())}")
    print(f"Pseudotime range: {cds_loaded.cell_metadata['pseudotime'].min():.2f} - "
          f"{cds_loaded.cell_metadata['pseudotime'].max():.2f}")
    
    return cds


if __name__ == "__main__":
    # Run example analysis
    cds = example_analysis()
    
    print("\n" + "=" * 50)
    print("PyMonocle3 is ready to use!")
    print("=" * 50)