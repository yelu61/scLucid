"""
Cell type annotation functions for single-cell RNA-seq data.

This module provides functions to annotate cell clusters with cell type labels
based on marker gene expression.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from sklearn.metrics.pairwise import cosine_similarity
from typing import Optional, Literal, Dict, List, Union

from .manager import Manager


def score_cell_types(
    adata, 
    marker_config: str, 
    key_added: str = "cell_type_scores",
    min_genes: int = 3,
    method: Literal["scanpy", "average"] = "scanpy",
) -> sc.AnnData:
    """
    Score cells based on cell type marker genes.
    
    Args:
        adata: AnnData object
        marker_config: Path to marker configuration file
        key_added: Key under which to add scores to adata.uns
        min_genes: Minimum number of marker genes required to compute a score
        method: Method used to compute scores ("scanpy" or "average")
        
    Returns:
        AnnData with cell type scores added to obs
    """
    mgr = Manager(marker_config)
    mgr.intersect_with(adata)
    
    # Get marker genes present in the dataset
    cell_type_markers = {}
    for cell_type, cell in mgr.CELLS.items():
        markers = cell.markers
        if len(markers) >= min_genes:
            cell_type_markers[cell_type] = markers
        else:
            print(f"Warning: Cell type '{cell_type}' has fewer than {min_genes} markers in the dataset ({len(markers)}).")
    
    # Score cells
    if method == "scanpy":
        for cell_type, markers in cell_type_markers.items():
            sc.tl.score_genes(adata, markers, score_name=f"{cell_type}_score")
    elif method == "average":
        # Get expression matrix
        if "normalized" in adata.layers:
            X = adata.layers["normalized"]
        elif "log1p_norm" in adata.layers:
            X = adata.layers["log1p_norm"]
        else:
            X = adata.X
            
        # Calculate average expression scores
        for cell_type, markers in cell_type_markers.items():
            marker_indices = [adata.var_names.get_loc(gene) for gene in markers if gene in adata.var_names]
            if len(marker_indices) > 0:
                if isinstance(X, np.ndarray):
                    scores = np.mean(X[:, marker_indices], axis=1)
                else:  # sparse matrix
                    scores = np.array(X[:, marker_indices].mean(axis=1)).flatten()
                adata.obs[f"{cell_type}_score"] = scores
    
    # Collect all score columns
    score_cols = [col for col in adata.obs.columns if col.endswith("_score")]
    adata.uns[key_added] = {
        "cell_types": [col.replace("_score", "") for col in score_cols],
        "method": method,
        "min_genes": min_genes
    }
    
    return adata


def annotate_clusters(
    adata,
    cluster_key: str,
    marker_config: str,
    method: Literal["correlation", "max_score", "marker_enrichment"] = "correlation",
    use_raw: bool = False,
    threshold: float = 0.05,
    key_added: Optional[str] = None,
) -> sc.AnnData:
    """
    Annotate clusters using marker genes.
    
    Args:
        adata: AnnData object
        cluster_key: Key in adata.obs containing cluster assignments
        marker_config: Path to marker configuration file
        method: Method for annotation ("correlation", "max_score", or "marker_enrichment")
        use_raw: Whether to use raw counts for marker enrichment analysis
        threshold: P-value threshold for marker enrichment significance
        key_added: Key under which to add annotations (defaults to f"{cluster_key}_annotation")
        
    Returns:
        AnnData with cluster annotations added to obs
    """
    if method == "correlation":
        return _annotate_by_correlation(adata, cluster_key, marker_config, key_added)
    elif method == "max_score":
        return _annotate_by_max_score(adata, cluster_key, marker_config, key_added)
    elif method == "marker_enrichment":
        return _annotate_by_marker_enrichment(adata, cluster_key, marker_config, use_raw, threshold, key_added)
    else:
        raise ValueError(f"Unknown method: {method}")


def _annotate_by_correlation(
    adata,
    cluster_key: str,
    marker_config: str,
    key_added: Optional[str] = None,
) -> sc.AnnData:
    """
    Annotate clusters by correlation with marker gene expression.
    
    Computes the correlation between each cluster's average gene expression profile
    and the binary marker profile of each cell type.
    
    Args:
        adata: AnnData object
        cluster_key: Key in adata.obs containing cluster assignments
        marker_config: Path to marker configuration file
        key_added: Key under which to add annotations
        
    Returns:
        AnnData with cluster annotations added to obs
    """
    if key_added is None:
        key_added = f"{cluster_key}_annotation"
    
    mgr = Manager(marker_config)
    mgr.intersect_with(adata)
    
    # Get expression matrix
    if "normalized" in adata.layers:
        X = adata.layers["normalized"]
    elif "log1p_norm" in adata.layers:
        X = adata.layers["log1p_norm"]
    else:
        X = adata.X
        
    # Compute average expression for each cluster
    clusters = adata.obs[cluster_key].cat.categories
    cluster_means = np.zeros((len(clusters), adata.n_vars))
    
    for i, cluster in enumerate(clusters):
        mask = adata.obs[cluster_key] == cluster
        if isinstance(X, np.ndarray):
            cluster_means[i] = np.mean(X[mask], axis=0)
        else:  # sparse matrix
            cluster_means[i] = np.array(X[mask].mean(axis=0)).flatten()
    
    # Create binary marker profiles for each cell type
    cell_types = []
    marker_profiles = []
    
    for cell_type, cell in mgr.CELLS.items():
        markers = cell.markers
        if len(markers) >= 3:  # At least 3 markers
            cell_types.append(cell_type)
            profile = np.zeros(adata.n_vars)
            for marker in markers:
                if marker in adata.var_names:
                    idx = adata.var_names.get_loc(marker)
                    profile[idx] = 1
            marker_profiles.append(profile)
    
    marker_profiles = np.vstack(marker_profiles)
    
    # Compute correlations
    correlations = cosine_similarity(cluster_means, marker_profiles)
    
    # Assign the most correlated cell type to each cluster
    annotations = {}
    correlation_scores = {}
    
    for i, cluster in enumerate(clusters):
        if len(cell_types) > 0:
            best_idx = np.argmax(correlations[i])
            best_cell_type = cell_types[best_idx]
            correlation_score = correlations[i, best_idx]
            
            annotations[cluster] = best_cell_type
            correlation_scores[cluster] = correlation_score
        else:
            annotations[cluster] = "Unknown"
            correlation_scores[cluster] = 0.0
    
    # Add annotations to adata
    annotation_series = pd.Series(
        adata.obs[cluster_key].map(annotations).values,
        index=adata.obs.index,
        name=key_added
    )
    adata.obs[key_added] = annotation_series
    adata.uns[f"{key_added}_correlation"] = correlation_scores
    
    # Add color mapping
    adata.uns[f"{key_added}_colors"] = {
        cell_type: mgr[cell_type].color 
        for cell_type in set(annotations.values())
        if cell_type in mgr.CELLS and mgr[cell_type].color
    }
    
    return adata


def _annotate_by_max_score(
    adata,
    cluster_key: str,
    marker_config: str,
    key_added: Optional[str] = None,
) -> sc.AnnData:
    """
    Annotate clusters by maximum marker gene score.
    
    Uses the cell type score calculated by score_cell_types to assign
    the cell type with the highest average score to each cluster.
    
    Args:
        adata: AnnData object
        cluster_key: Key in adata.obs containing cluster assignments
        marker_config: Path to marker configuration file
        key_added: Key under which to add annotations
        
    Returns:
        AnnData with cluster annotations added to obs
    """
    if key_added is None:
        key_added = f"{cluster_key}_annotation"
    
    # First, score cell types if not already done
    if not any(col.endswith("_score") for col in adata.obs.columns):
        adata = score_cell_types(adata, marker_config)
    
    # Get cell type names from score columns
    cell_types = [col.replace("_score", "") for col in adata.obs.columns if col.endswith("_score")]
    
    if not cell_types:
        raise ValueError("No cell type scores found. Run score_cell_types first.")
    
    # Calculate average score per cluster
    clusters = adata.obs[cluster_key].cat.categories
    cluster_scores = {}
    
    for cluster in clusters:
        mask = adata.obs[cluster_key] == cluster
        cluster_scores[cluster] = {
            cell_type: adata.obs.loc[mask, f"{cell_type}_score"].mean()
            for cell_type in cell_types
        }
    
    # Assign the cell type with highest score to each cluster
    annotations = {}
    max_scores = {}
    
    for cluster, scores in cluster_scores.items():
        if scores:
            best_cell_type = max(scores.items(), key=lambda x: x[1])[0]
            max_score = scores[best_cell_type]
            
            annotations[cluster] = best_cell_type
            max_scores[cluster] = max_score
        else:
            annotations[cluster] = "Unknown"
            max_scores[cluster] = 0.0
    
    # Add annotations to adata
    annotation_series = pd.Series(
        adata.obs[cluster_key].map(annotations).values,
        index=adata.obs.index,
        name=key_added
    )
    adata.obs[key_added] = annotation_series
    adata.uns[f"{key_added}_scores"] = max_scores
    
    # Add color mapping from manager
    mgr = Manager(marker_config)
    adata.uns[f"{key_added}_colors"] = {
        cell_type: mgr[cell_type].color 
        for cell_type in set(annotations.values())
        if cell_type in mgr.CELLS and mgr[cell_type].color
    }
    
    return adata


def _annotate_by_marker_enrichment(
    adata,
    cluster_key: str,
    marker_config: str,
    use_raw: bool = False,
    threshold: float = 0.05,
    key_added: Optional[str] = None,
) -> sc.AnnData:
    """
    Annotate clusters by marker gene enrichment analysis.
    
    Performs statistical tests to identify which cell type markers are
    significantly enriched in each cluster.
    
    Args:
        adata: AnnData object
        cluster_key: Key in adata.obs containing cluster assignments
        marker_config: Path to marker configuration file
        use_raw: Whether to use raw counts for enrichment analysis
        threshold: P-value threshold for significance
        key_added: Key under which to add annotations
        
    Returns:
        AnnData with cluster annotations added to obs
    """
    if key_added is None:
        key_added = f"{cluster_key}_annotation"
    
    mgr = Manager(marker_config)
    mgr.intersect_with(adata)
    
    # Get expression data
    if use_raw and adata.raw is not None:
        X = adata.raw.X
        var_names = adata.raw.var_names
    else:
        X = adata.X
        var_names = adata.var_names
    
    # Prepare for enrichment analysis
    clusters = adata.obs[cluster_key].cat.categories
    cell_types = [name for name, cell in mgr.CELLS.items() if len(cell.markers) >= 3]
    
    # Store enrichment results
    enrichment_results = {}
    
    # For each cluster, calculate marker enrichment for each cell type
    for cluster in clusters:
        cluster_mask = adata.obs[cluster_key] == cluster
        enrichment_results[cluster] = {}
        
        for cell_type in cell_types:
            markers = mgr[cell_type].markers
            marker_indices = [var_names.get_loc(gene) for gene in markers if gene in var_names]
            
            if len(marker_indices) < 3:
                continue
                
            # Calculate mean expression in cluster vs outside
            if isinstance(X, np.ndarray):
                in_cluster_expr = np.mean(X[cluster_mask][:, marker_indices], axis=0)
                out_cluster_expr = np.mean(X[~cluster_mask][:, marker_indices], axis=0)
            else:  # sparse matrix
                in_cluster_expr = np.array(X[cluster_mask][:, marker_indices].mean(axis=0)).flatten()
                out_cluster_expr = np.array(X[~cluster_mask][:, marker_indices].mean(axis=0)).flatten()
            
            # Perform statistical test (Mann-Whitney U test)
            try:
                statistic, pvalue = stats.mannwhitneyu(in_cluster_expr, out_cluster_expr, alternative='greater')
                enrichment_results[cluster][cell_type] = {
                    'pvalue': pvalue,
                    'statistic': statistic,
                    'fold_change': np.mean(in_cluster_expr) / np.mean(out_cluster_expr) if np.mean(out_cluster_expr) > 0 else float('inf'),
                    'mean_in_cluster': np.mean(in_cluster_expr),
                    'mean_out_cluster': np.mean(out_cluster_expr)
                }
            except ValueError:
                # Handle ties or other errors
                enrichment_results[cluster][cell_type] = {
                    'pvalue': 1.0,
                    'statistic': 0.0,
                    'fold_change': 1.0,
                    'mean_in_cluster': np.mean(in_cluster_expr),
                    'mean_out_cluster': np.mean(out_cluster_expr)
                }
    
    # Assign cell types based on lowest p-value
    annotations = {}
    pvalues = {}
    
    for cluster in clusters:
        if cluster in enrichment_results and enrichment_results[cluster]:
            # Filter by threshold and find best match
            significant = {
                cell_type: result for cell_type, result in enrichment_results[cluster].items()
                if result['pvalue'] < threshold
            }
            
            if significant:
                best_cell_type = min(significant.items(), key=lambda x: x[1]['pvalue'])[0]
                annotations[cluster] = best_cell_type
                pvalues[cluster] = enrichment_results[cluster][best_cell_type]['pvalue']
            else:
                annotations[cluster] = "Unknown"
                pvalues[cluster] = 1.0
        else:
            annotations[cluster] = "Unknown"
            pvalues[cluster] = 1.0
    
    # Add annotations to adata
    annotation_series = pd.Series(
        adata.obs[cluster_key].map(annotations).values,
        index=adata.obs.index,
        name=key_added
    )
    adata.obs[key_added] = annotation_series
    adata.uns[f"{key_added}_pvalues"] = pvalues
    adata.uns[f"{key_added}_enrichment"] = enrichment_results
    
    # Add color mapping
    adata.uns[f"{key_added}_colors"] = {
        cell_type: mgr[cell_type].color 
        for cell_type in set(annotations.values())
        if cell_type in mgr.CELLS and mgr[cell_type].color
    }
    
    return adata


def extract_high_confidence_cells(
    adata,
    annotation_key: str,
    score_threshold: float = 0.5,
    percentile: float = 75,
    min_cells: int = 10,
    inplace: bool = False,
) -> sc.AnnData:
    """
    Extract cells with high confidence annotations based on cell type scores.
    
    Args:
        adata: AnnData object
        annotation_key: Key in adata.obs containing cell type annotations
        score_threshold: Absolute threshold for cell type scores
        percentile: Percentile threshold within each cell type
        min_cells: Minimum number of cells to keep for each cell type
        inplace: Whether to modify adata inplace
        
    Returns:
        AnnData with high confidence cells or a boolean mask
    """
    if not inplace:
        adata = adata.copy()
    
    # Get unique cell types
    cell_types = adata.obs[annotation_key].unique()
    
    # Create mask for high confidence cells
    high_confidence = np.zeros(adata.n_obs, dtype=bool)
    
    for cell_type in cell_types:
        # Skip unknown/unassigned cells
        if cell_type == "Unknown":
            continue
            
        # Get cells of this type
        cell_mask = adata.obs[annotation_key] == cell_type
        
        # Get score column
        score_col = f"{cell_type}_score"
        if score_col not in adata.obs.columns:
            print(f"Warning: Score column '{score_col}' not found for cell type '{cell_type}'")
            continue
            
        # Get scores for this cell type
        scores = adata.obs.loc[cell_mask, score_col]
        
        # Apply absolute and percentile thresholds
        abs_threshold = score_threshold
        perc_threshold = np.percentile(scores, percentile) if len(scores) > 0 else 0
        threshold = max(abs_threshold, perc_threshold)
        
        # Mark high confidence cells
        high_conf_mask = (adata.obs[score_col] >= threshold) & cell_mask
        
        # Ensure minimum number of cells
        if sum(high_conf_mask) < min_cells and sum(cell_mask) >= min_cells:
            # Take top min_cells by score
            top_indices = adata.obs.loc[cell_mask, score_col].nlargest(min_cells).index
            high_conf_mask = adata.obs.index.isin(top_indices)
        
        # Update overall mask
        high_confidence = high_confidence | high_conf_mask
    
    if inplace:
        adata.obs["high_confidence"] = high_confidence
        return adata
    else:
        return adata[high_confidence]