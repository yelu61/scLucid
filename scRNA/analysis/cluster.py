"""
Clustering functions for single-cell RNA-seq data.

This module provides functions for marker-guided clustering,
optimal resolution selection, and cluster evaluation.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import entropy
from sklearn import metrics
from typing import Optional, List, Tuple, Dict, Literal, Union
import matplotlib.pyplot as plt

from .manager import Manager


def marker_guided_clustering(
    adata,
    marker_config: str,
    resolution_range: Tuple[float, float, int] = (0.1, 2.0, 10),
    metric: Literal["marker_separation", "silhouette", "calinski_harabasz"] = "marker_separation",
    clustering_method: Literal["leiden", "louvain"] = "leiden",
    neighbors_key: Optional[str] = None,
    key_added: Optional[str] = None,
    plot: bool = True,
    min_cluster_genes: int = 3,
    random_state: int = 42,
) -> sc.AnnData:
    """
    Perform clustering with resolution optimization guided by marker genes.
    
    Args:
        adata: AnnData object
        marker_config: Path to marker configuration file
        resolution_range: (start, end, steps) for resolution search
        metric: Metric to evaluate clustering quality
        clustering_method: Method to use for clustering ("leiden" or "louvain")
        neighbors_key: Key to use for neighbors graph. If None, uses "neighbors"
        key_added: Key under which to add the clustering result
        plot: Whether to plot the evaluation metrics for different resolutions
        min_cluster_genes: Minimum number of marker genes required to evaluate a cell type
        random_state: Random seed for clustering
        
    Returns:
        AnnData with optimized clustering results
    """
    # Prepare marker genes
    mgr = Manager(marker_config)
    mgr.intersect_with(adata)
    
    # Get marker genes dictionary filtered by min_cluster_genes
    marker_genes = {
        cell_type: cell.markers 
        for cell_type, cell in mgr.CELLS.items() 
        if len(cell.markers) >= min_cluster_genes
    }
    
    # Generate resolutions to try
    start, end, steps = resolution_range
    resolutions = np.linspace(start, end, steps)
    
    # Set key_added
    if key_added is None:
        key_added = clustering_method
    
    # Create dictionary to store evaluation results
    eval_results = []
    
    # Run clustering at different resolutions
    for res in resolutions:
        # Perform clustering
        if clustering_method == "leiden":
            sc.tl.leiden(adata, resolution=res, key_added=f"{key_added}_{res:.2f}", 
                         neighbors_key=neighbors_key, random_state=random_state)
            cluster_key = f"{key_added}_{res:.2f}"
        else:
            sc.tl.louvain(adata, resolution=res, key_added=f"{key_added}_{res:.2f}", 
                          neighbors_key=neighbors_key, random_state=random_state)
            cluster_key = f"{key_added}_{res:.2f}"
        
        # Evaluate clustering
        if metric == "marker_separation":
            score = evaluate_marker_separation(adata, cluster_key, marker_genes)
        elif metric == "silhouette":
            score = evaluate_silhouette(adata, cluster_key)
        elif metric == "calinski_harabasz":
            score = evaluate_calinski_harabasz(adata, cluster_key)
        else:
            raise ValueError(f"Unknown metric: {metric}")
        
        # Store results
        n_clusters = len(adata.obs[cluster_key].unique())
        eval_results.append({
            'resolution': res,
            'n_clusters': n_clusters,
            'score': score
        })
        
        print(f"Resolution {res:.2f}: {n_clusters} clusters, {metric} score = {score:.4f}")
    
    # Convert to DataFrame
    eval_df = pd.DataFrame(eval_results)
    
    # Find optimal resolution
    if metric in ["marker_separation", "silhouette", "calinski_harabasz"]:
        # Higher is better
        optimal_idx = eval_df['score'].idxmax()
    else:
        # This shouldn't happen with current metrics, but for future-proofing
        raise ValueError(f"Unknown optimization direction for metric: {metric}")
    
    optimal_res = eval_df.loc[optimal_idx, 'resolution']
    optimal_score = eval_df.loc[optimal_idx, 'score']
    optimal_n_clusters = eval_df.loc[optimal_idx, 'n_clusters']
    
    print(f"\nOptimal resolution: {optimal_res:.2f} with {optimal_n_clusters} clusters")
    print(f"Optimal {metric} score: {optimal_score:.4f}")
    
    # Use the optimal resolution for the final clustering
    if clustering_method == "leiden":
        sc.tl.leiden(adata, resolution=optimal_res, key_added=key_added, 
                     neighbors_key=neighbors_key, random_state=random_state)
    else:
        sc.tl.louvain(adata, resolution=optimal_res, key_added=key_added, 
                      neighbors_key=neighbors_key, random_state=random_state)
    
    # Store evaluation results
    adata.uns[f"{key_added}_evaluation"] = {
        'metric': metric,
        'results': eval_df.to_dict('records'),
        'optimal_resolution': optimal_res,
        'optimal_score': optimal_score,
        'optimal_n_clusters': optimal_n_clusters
    }
    
    # Plot evaluation metrics
    if plot:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Plot score vs resolution
        axes[0].plot(eval_df['resolution'], eval_df['score'], 'o-')
        axes[0].set_xlabel('Resolution')
        axes[0].set_ylabel(metric.replace('_', ' ').title())
        axes[0].set_title(f"{metric.replace('_', ' ').title()} vs Resolution")
        axes[0].axvline(x=optimal_res, color='r', linestyle='--')
        axes[0].grid(True)
        
        # Plot number of clusters vs resolution
        axes[1].plot(eval_df['resolution'], eval_df['n_clusters'], 'o-')
        axes[1].set_xlabel('Resolution')
        axes[1].set_ylabel('Number of Clusters')
        axes[1].set_title('Number of Clusters vs Resolution')
        axes[1].axvline(x=optimal_res, color='r', linestyle='--')
        axes[1].grid(True)
        
        plt.tight_layout()
        plt.show()
    
    return adata


def evaluate_resolution(
    adata,
    resolution_range: Tuple[float, float, int] = (0.1, 2.0, 10),
    clustering_method: Literal["leiden", "louvain"] = "leiden",
    metric: Literal["silhouette", "calinski_harabasz", "davies_bouldin"] = "silhouette",
    representation: str = "X_pca",
    neighbors_key: Optional[str] = None,
    plot: bool = True,
    random_state: int = 42,
) -> Dict:
    """
    Evaluate different clustering resolutions using internal metrics.
    
    Args:
        adata: AnnData object
        resolution_range: (start, end, steps) for resolution search
        clustering_method: Method to use for clustering ("leiden" or "louvain")
        metric: Metric to evaluate clustering quality
        representation: Representation to use for metric calculation
        neighbors_key: Key to use for neighbors graph. If None, uses "neighbors"
        plot: Whether to plot the evaluation metrics
        random_state: Random seed for clustering
        
    Returns:
        Dictionary containing evaluation results and optimal resolution
    """
    # Generate resolutions to try
    start, end, steps = resolution_range
    resolutions = np.linspace(start, end, steps)
    
    # Create dictionary to store evaluation results
    eval_results = []
    
    # Run clustering at different resolutions
    for res in resolutions:
        # Perform clustering
        if clustering_method == "leiden":
            sc.tl.leiden(adata, resolution=res, key_added=f"leiden_res{res:.2f}", 
                         neighbors_key=neighbors_key, random_state=random_state)
            cluster_key = f"leiden_res{res:.2f}"
        else:
            sc.tl.louvain(adata, resolution=res, key_added=f"louvain_res{res:.2f}", 
                          neighbors_key=neighbors_key, random_state=random_state)
            cluster_key = f"louvain_res{res:.2f}"
        
        # Get cluster labels
        labels = adata.obs[cluster_key].cat.codes.values
        n_clusters = len(np.unique(labels))
        
        # Skip if only one cluster (some metrics will fail)
        if n_clusters <= 1:
            eval_results.append({
                'resolution': res,
                'n_clusters': n_clusters,
                'score': np.nan
            })
            continue
        
        # Get representation for metric calculation
        X = adata.obsm[representation]
        
        # Calculate metric
        if metric == "silhouette":
            score = metrics.silhouette_score(X, labels)
        elif metric == "calinski_harabasz":
            score = metrics.calinski_harabasz_score(X, labels)
        elif metric == "davies_bouldin":
            score = -metrics.davies_bouldin_score(X, labels)  # Negate so higher is better
        else:
            raise ValueError(f"Unknown metric: {metric}")
        
        # Store results
        eval_results.append({
            'resolution': res,
            'n_clusters': n_clusters,
            'score': score
        })
        
        print(f"Resolution {res:.2f}: {n_clusters} clusters, {metric} score = {score:.4f}")
    
    # Convert to DataFrame
    eval_df = pd.DataFrame(eval_results)
    
    # Find optimal resolution (ignoring NaN values)
    if metric in ["silhouette", "calinski_harabasz"] or metric == "davies_bouldin":
        # Higher is better (davies_bouldin is already negated)
        valid_df = eval_df.dropna(subset=['score'])
        if len(valid_df) > 0:
            optimal_idx = valid_df['score'].idxmax()
            optimal_res = eval_df.loc[optimal_idx, 'resolution']
            optimal_score = eval_df.loc[optimal_idx, 'score']
            optimal_n_clusters = eval_df.loc[optimal_idx, 'n_clusters']
        else:
            optimal_res = resolutions[0]
            optimal_score = np.nan
            optimal_n_clusters = 0
    else:
        # This shouldn't happen with current metrics, but for future-proofing
        raise ValueError(f"Unknown optimization direction for metric: {metric}")
    
    print(f"\nOptimal resolution: {optimal_res:.2f} with {optimal_n_clusters} clusters")
    print(f"Optimal {metric} score: {optimal_score:.4f}")
    
    # Plot evaluation metrics
    if plot:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Plot score vs resolution
        axes[0].plot(eval_df['resolution'], eval_df['score'], 'o-')
        axes[0].set_xlabel('Resolution')
        axes[0].set_ylabel(metric.replace('_', ' ').title())
        axes[0].set_title(f"{metric.replace('_', ' ').title()} vs Resolution")
        axes[0].axvline(x=optimal_res, color='r', linestyle='--')
        axes[0].grid(True)
        
        # Plot number of clusters vs resolution
        axes[1].plot(eval_df['resolution'], eval_df['n_clusters'], 'o-')
        axes[1].set_xlabel('Resolution')
        axes[1].set_ylabel('Number of Clusters')
        axes[1].set_title('Number of Clusters vs Resolution')
        axes[1].axvline(x=optimal_res, color='r', linestyle='--')
        axes[1].grid(True)
        
        plt.tight_layout()
        plt.show()
    
    # Return evaluation results
    return {
        'metric': metric,
        'results': eval_df.to_dict('records'),
        'optimal_resolution': optimal_res,
        'optimal_score': optimal_score,
        'optimal_n_clusters': optimal_n_clusters
    }


def evaluate_marker_separation(
    adata,
    cluster_key: str,
    marker_genes: Dict[str, List[str]],
) -> float:
    """
    Evaluate clustering by how well it separates marker genes.
    
    Calculates a score based on how exclusively marker genes are
    expressed in different clusters.
    
    Args:
        adata: AnnData object
        cluster_key: Key in adata.obs containing cluster assignments
        marker_genes: Dictionary mapping cell types to lists of marker genes
        
    Returns:
        Marker separation score (higher is better)
    """
    # Get expression matrix (normalized or log-normalized)
    if "normalized" in adata.layers:
        X = adata.layers["normalized"]
    elif "log1p_norm" in adata.layers:
        X = adata.layers["log1p_norm"]
    else:
        X = adata.X
    
    # Get clusters
    clusters = adata.obs[cluster_key].cat.categories
    n_clusters = len(clusters)
    
    # Initialize scores
    cluster_scores = []
    
    # For each marker set, calculate score
    for cell_type, markers in marker_genes.items():
        # Get indices of marker genes
        marker_indices = [adata.var_names.get_loc(gene) for gene in markers if gene in adata.var_names]
        
        if not marker_indices:
            continue
        
        # Calculate mean expression of markers in each cluster
        cluster_means = np.zeros(n_clusters)
        
        for i, cluster in enumerate(clusters):
            mask = adata.obs[cluster_key] == cluster
            if isinstance(X, np.ndarray):
                cluster_means[i] = np.mean(X[mask][:, marker_indices])
            else:  # sparse matrix
                cluster_means[i] = np.mean(X[mask][:, marker_indices].mean())
        
        # Skip if all means are zero
        if np.sum(cluster_means) == 0:
            continue
        
        # Normalize to get a probability distribution
        cluster_probs = cluster_means / np.sum(cluster_means)
        
        # Calculate entropy (lower entropy means better separation)
        marker_entropy = entropy(cluster_probs)
        
        # Convert to a score where higher is better
        # Max entropy is log(n_clusters), so we normalize to [0, 1] and invert
        max_entropy = np.log(n_clusters)
        if max_entropy > 0:
            marker_score = 1 - (marker_entropy / max_entropy)
        else:
            marker_score = 0
        
        cluster_scores.append(marker_score)
    
    # Return average score across all marker sets
    if cluster_scores:
        return np.mean(cluster_scores)
    else:
        return 0.0


def evaluate_silhouette(
    adata,
    cluster_key: str,
    representation: str = "X_pca",
) -> float:
    """
    Evaluate clustering using silhouette score.
    
    Args:
        adata: AnnData object
        cluster_key: Key in adata.obs containing cluster assignments
        representation: Representation to use for distance calculation
        
    Returns:
        Silhouette score (higher is better)
    """
    # Get cluster labels
    labels = adata.obs[cluster_key].cat.codes.values
    
    # Check if we have more than one cluster
    if len(np.unique(labels)) <= 1:
        return 0.0
    
    # Get representation
    X = adata.obsm[representation]
    
    # Calculate silhouette score
    return metrics.silhouette_score(X, labels)


def evaluate_calinski_harabasz(
    adata,
    cluster_key: str,
    representation: str = "X_pca",
) -> float:
    """
    Evaluate clustering using Calinski-Harabasz index.
    
    Args:
        adata: AnnData object
        cluster_key: Key in adata.obs containing cluster assignments
        representation: Representation to use for distance calculation
        
    Returns:
        Calinski-Harabasz index (higher is better)
    """
    # Get cluster labels
    labels = adata.obs[cluster_key].cat.codes.values
    
    # Check if we have more than one cluster
    if len(np.unique(labels)) <= 1:
        return 0.0
    
    # Get representation
    X = adata.obsm[representation]
    
    # Calculate Calinski-Harabasz index
    return metrics.calinski_harabasz_score(X, labels)


def merge_clusters(
    adata,
    cluster_key: str,
    marker_config: str,
    similarity_threshold: float = 0.8,
    method: Literal["marker_overlap", "expression_correlation"] = "marker_overlap",
    inplace: bool = True,
    key_added: Optional[str] = None,
) -> sc.AnnData:
    """
    Merge similar clusters based on marker overlap or expression correlation.
    
    Args:
        adata: AnnData object
        cluster_key: Key in adata.obs containing cluster assignments
        marker_config: Path to marker configuration file
        similarity_threshold: Threshold for merging clusters
        method: Method to calculate similarity between clusters
        inplace: Whether to modify adata inplace
        key_added: Key under which to add the merged clustering result
        
    Returns:
        AnnData with merged clusters
    """
    if not inplace:
        adata = adata.copy()
    
    if key_added is None:
        key_added = f"{cluster_key}_merged"
    
    # Get clusters
    clusters = adata.obs[cluster_key].cat.categories
    n_clusters = len(clusters)
    
    # Calculate similarity matrix
    if method == "marker_overlap":
        # Find top marker genes for each cluster
        sc.tl.rank_genes_groups(adata, groupby=cluster_key, method='wilcoxon')
        
        # Get top markers for each cluster
        top_markers = {}
        for i, cluster in enumerate(clusters):
            genes = [gene for gene in adata.uns['rank_genes_groups']['names'][str(i)]]
            scores = [score for score in adata.uns['rank_genes_groups']['scores'][str(i)]]
            
            # Keep only positive markers with score > 0
            pos_markers = [(gene, score) for gene, score in zip(genes, scores) if score > 0]
            top_markers[cluster] = [gene for gene, _ in pos_markers[:50]]  # Top 50 positive markers
        
        # Calculate Jaccard similarity
        similarity_matrix = np.zeros((n_clusters, n_clusters))
        for i, c1 in enumerate(clusters):
            for j, c2 in enumerate(clusters):
                if i == j:
                    similarity_matrix[i, j] = 1.0
                else:
                    set1 = set(top_markers[c1])
                    set2 = set(top_markers[c2])
                    if not set1 or not set2:
                        similarity_matrix[i, j] = 0.0
                    else:
                        similarity_matrix[i, j] = len(set1.intersection(set2)) / len(set1.union(set2))
    
    elif method == "expression_correlation":
        # Get expression matrix
        if "normalized" in adata.layers:
            X = adata.layers["normalized"]
        elif "log1p_norm" in adata.layers:
            X = adata.layers["log1p_norm"]
        else:
            X = adata.X
            
        # Calculate mean expression profile for each cluster
        mean_profiles = np.zeros((n_clusters, adata.n_vars))
        for i, cluster in enumerate(clusters):
            mask = adata.obs[cluster_key] == cluster
            if isinstance(X, np.ndarray):
                mean_profiles[i] = np.mean(X[mask], axis=0)
            else:  # sparse matrix
                mean_profiles[i] = np.array(X[mask].mean(axis=0)).flatten()
        
        # Calculate correlation matrix
        similarity_matrix = np.corrcoef(mean_profiles)
    
    else:
        raise ValueError(f"Unknown method: {method}")
    
    # Create cluster graph based on similarity
    import networkx as nx
    G = nx.Graph()
    
    # Add nodes
    for i, cluster in enumerate(clusters):
        G.add_node(cluster)
    
    # Add edges for clusters with similarity above threshold
    for i, c1 in enumerate(clusters):
        for j, c2 in enumerate(clusters):
            if i < j and similarity_matrix[i, j] >= similarity_threshold:
                G.add_edge(c1, c2, weight=similarity_matrix[i, j])
    
    # Find connected components (clusters to merge)
    connected_components = list(nx.connected_components(G))
    
    # Create mapping from old to new clusters
    cluster_mapping = {}
    for i, component in enumerate(connected_components):
        for cluster in component:
            cluster_mapping[cluster] = f"Cluster_{i+1}"
    
    # Apply mapping to create new cluster assignments
    adata.obs[key_added] = adata.obs[cluster_key].map(cluster_mapping).astype('category')
    
    # Print merge results
    print(f"Original clusters: {n_clusters}")
    print(f"Merged clusters: {len(connected_components)}")
    
    # Add merge info to adata.uns
    adata.uns[f"{key_added}_info"] = {
        'original_key': cluster_key,
        'similarity_threshold': similarity_threshold,
        'method': method,
        'original_clusters': n_clusters,
        'merged_clusters': len(connected_components),
        'mapping': cluster_mapping
    }
    
    return adata