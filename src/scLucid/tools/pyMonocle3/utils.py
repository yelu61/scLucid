"""
Utility functions for pyMonocle3 (R-free)
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
from typing import Optional, List, Dict, Union, Tuple, Any
import logging

log = logging.getLogger(__name__)


def detect_sparse_type(matrix) -> str:
    """Detect sparse matrix type"""
    if sp.issparse(matrix):
        return type(matrix).__name__
    return "dense"


def convert_to_dense(matrix) -> np.ndarray:
    """Convert matrix to dense numpy array"""
    if sp.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix)


def convert_to_sparse(matrix, format: str = 'csr') -> sp.spmatrix:
    """Convert matrix to sparse format"""
    if sp.issparse(matrix):
        if format == 'csr':
            return matrix.tocsr()
        elif format == 'csc':
            return matrix.tocsc()
        elif format == 'coo':
            return matrix.tocoo()
    return sp.csr_matrix(matrix) if format == 'csr' else sp.csc_matrix(matrix)


def calculate_size_factors(
    counts: np.ndarray,
    method: str = "mean-geometric-mean-total"
) -> np.ndarray:
    """
    Calculate size factors for normalization

    Parameters
    ----------
    counts : np.ndarray
        Count matrix (genes x cells)
    method : str
        Method for size factor calculation

    Returns
    -------
    np.ndarray
        Size factors for each cell
    """
    # Calculate total counts per cell
    total_counts = np.sum(counts, axis=0)

    if method == "mean-geometric-mean-total":
        # Calculate geometric mean of totals
        log_totals = np.log(total_counts + 1)
        geometric_mean = np.exp(np.mean(log_totals)) - 1
        size_factors = total_counts / geometric_mean
    elif method == "median":
        median_total = np.median(total_counts)
        size_factors = total_counts / median_total
    elif method == "upper-quartile":
        upper_q = np.percentile(total_counts, 75)
        size_factors = total_counts / upper_q
    else:
        raise ValueError(f"Unknown method: {method}")

    return size_factors


def normalize_expression(
    counts: np.ndarray,
    size_factors: Optional[np.ndarray] = None,
    method: str = "log",
    pseudo_count: float = 1.0
) -> np.ndarray:
    """
    Normalize expression data

    Parameters
    ----------
    counts : np.ndarray
        Count matrix (genes x cells)
    size_factors : np.ndarray, optional
        Pre-calculated size factors
    method : str
        Normalization method
    pseudo_count : float
        Pseudo-count for log transformation

    Returns
    -------
    np.ndarray
        Normalized expression
    """
    if size_factors is None:
        size_factors = calculate_size_factors(counts)

    # Size factor normalization
    norm_counts = counts / size_factors[np.newaxis, :]

    if method == "log":
        return np.log1p(norm_counts + pseudo_count)
    elif method == "sqrt":
        return np.sqrt(norm_counts + pseudo_count)
    elif method == "none":
        return norm_counts
    else:
        raise ValueError(f"Unknown method: {method}")


def select_highly_variable_genes(
    expression: np.ndarray,
    gene_names: Optional[List[str]] = None,
    n_top_genes: int = 2000,
    min_mean: float = 0.0125,
    max_mean: float = 3.0,
    min_disp: float = 0.5,
) -> pd.DataFrame:
    """
    Select highly variable genes

    Parameters
    ----------
    expression : np.ndarray
        Expression matrix (genes x cells)
    gene_names : list, optional
        Gene names
    n_top_genes : int
        Number of top genes to select
    min_mean : float
        Minimum mean expression
    max_mean : float
        Maximum mean expression
    min_disp : float
        Minimum dispersion

    Returns
    -------
    pd.DataFrame
        HVG statistics
    """
    if gene_names is None:
        gene_names = [f"gene_{i}" for i in range(expression.shape[0])]

    # Calculate mean and variance
    mean_expr = np.mean(expression, axis=1)
    var_expr = np.var(expression, axis=1)

    # Calculate dispersion
    dispersion = var_expr / (mean_expr + 1e-10)

    # Filter by mean
    mean_filter = (mean_expr >= min_mean) & (mean_expr <= max_mean)

    # Create results DataFrame
    hvg_stats = pd.DataFrame({
        'gene': gene_names,
        'means': mean_expr,
        'variances': var_expr,
        'dispersions': dispersion,
        'highly_variable': mean_filter & (dispersion >= min_disp),
    })

    # Rank by dispersion and select top genes
    hvg_stats['dispersion_rank'] = hvg_stats['dispersions'].rank(ascending=False)
    hvg_stats.loc[
        hvg_stats['dispersion_rank'] <= n_top_genes, 'highly_variable'
    ] = True

    return hvg_stats.sort_values('dispersions', ascending=False)


def calculate_gene_loadings(
    pca_model,
    gene_names: List[str],
    n_components: int = 10
) -> pd.DataFrame:
    """
    Calculate gene loadings from PCA model

    Parameters
    ----------
    pca_model : sklearn.decomposition.PCA
        Fitted PCA model
    gene_names : list
        Gene names
    n_components : int
        Number of components

    Returns
    -------
    pd.DataFrame
        Gene loadings
    """
    n_components = min(n_components, len(pca_model.components_))

    loadings = pd.DataFrame(
        pca_model.components_[:n_components].T,
        index=gene_names,
        columns=[f"PC{i+1}" for i in range(n_components)]
    )

    return loadings


def aggregate_expression_by_group(
    expression: np.ndarray,
    groups: pd.Series,
    group_names: Optional[List] = None,
    aggregation: str = "mean"
) -> pd.DataFrame:
    """
    Aggregate expression by group

    Parameters
    ----------
    expression : np.ndarray
        Expression matrix (genes x cells)
    groups : pd.Series
        Group assignments for each cell
    group_names : list, optional
        Names for groups
    aggregation : str
        Aggregation method ("mean", "sum", "median")

    Returns
    -------
    pd.DataFrame
        Aggregated expression
    """
    if group_names is None:
        group_names = sorted(groups.unique())

    results = []
    for group in group_names:
        mask = groups == group
        if aggregation == "mean":
            agg = np.mean(expression[:, mask], axis=1)
        elif aggregation == "sum":
            agg = np.sum(expression[:, mask], axis=1)
        elif aggregation == "median":
            agg = np.median(expression[:, mask], axis=1)
        else:
            raise ValueError(f"Unknown aggregation: {aggregation}")

        results.append(agg)

    return pd.DataFrame(
        np.array(results).T,
        columns=group_names
    )


def calculate_correlation_matrix(
    data: np.ndarray,
    method: str = "pearson"
) -> np.ndarray:
    """
    Calculate correlation matrix

    Parameters
    ----------
    data : np.ndarray
        Data matrix
    method : str
        Correlation method ("pearson" or "spearman")

    Returns
    -------
    np.ndarray
        Correlation matrix
    """
    if method == "pearson":
        return np.corrcoef(data)
    elif method == "spearman":
        from scipy.stats import spearmanr
        corr, _ = spearmanr(data.T)
        return corr
    else:
        raise ValueError(f"Unknown method: {method}")


def subsample_cells(
    cds,
    n_cells: Optional[int] = None,
    frac: Optional[float] = None,
    random_state: int = 42
):
    """
    Subsample cells from CellDataSet

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    n_cells : int, optional
        Number of cells to keep
    frac : float, optional
        Fraction of cells to keep
    random_state : int
        Random seed

    Returns
    -------
    CellDataSet
        Subsampled CellDataSet
    """
    np.random.seed(random_state)

    if n_cells is not None:
        n_keep = min(n_cells, cds.n_cells)
    elif frac is not None:
        n_keep = int(cds.n_cells * frac)
    else:
        raise ValueError("Either n_cells or frac must be provided")

    keep_idx = np.random.choice(cds.n_cells, n_keep, replace=False)

    from .core import CellDataSet

    return CellDataSet(
        expression_data=cds.expression_data[:, keep_idx],
        cell_metadata=cds.cell_metadata.iloc[keep_idx],
        gene_metadata=cds.gene_metadata,
        reducedDims={k: v[keep_idx] for k, v in cds.reducedDims.items()},
        clusters=cds.clusters.iloc[keep_idx] if cds.clusters is not None else None,
        partitions=cds.partitions.iloc[keep_idx] if cds.partitions is not None else None,
        principal_graph=None,  # Graph needs to be recomputed
    )


def merge_datasets(
    cds_list: List,
    merge_groups: Optional[List[str]] = None,
) -> Any:
    """
    Merge multiple CellDataSets

    Parameters
    ----------
    cds_list : list
        List of CellDataSets
    merge_groups : list, optional
        Group labels for each dataset

    Returns
    -------
    CellDataSet
        Merged CellDataSet
    """
    from .core import CellDataSet

    if len(cds_list) == 0:
        raise ValueError("No datasets to merge")

    if merge_groups is None:
        merge_groups = [f"dataset_{i}" for i in range(len(cds_list))]

    # Concatenate expression
    expr_list = [cds.expression_data for cds in cds_list]
    merged_expr = sp.hstack(expr_list) if sp.issparse(expr_list[0]) else np.hstack(expr_list)

    # Concatenate metadata
    meta_list = []
    for cds, group in zip(cds_list, merge_groups):
        meta = cds.cell_metadata.copy()
        meta['dataset'] = group
        meta_list.append(meta)
    merged_meta = pd.concat(meta_list, ignore_index=True)

    # Use gene metadata from first dataset
    merged_gene_meta = cds_list[0].gene_metadata.copy()

    cds_merged = CellDataSet(
        expression_data=merged_expr,
        cell_metadata=merged_meta,
        gene_metadata=merged_gene_meta,
    )

    log.info(f"Merged {len(cds_list)} datasets: {cds_merged.n_cells} cells")

    return cds_merged


def validate_cds(cds) -> Tuple[bool, str]:
    """
    Validate CellDataSet integrity

    Parameters
    ----------
    cds : CellDataSet
        CellDataSet to validate

    Returns
    -------
    tuple
        (is_valid, error_message)
    """
    # Check dimensions
    n_genes, n_cells = cds.expression_data.shape

    if len(cds.cell_metadata) != n_cells:
        return False, f"Cell metadata mismatch: {len(cds.cell_metadata)} vs {n_cells}"

    if len(cds.gene_metadata) != n_genes:
        return False, f"Gene metadata mismatch: {len(cds.gene_metadata)} vs {n_genes}"

    # Check reduced dimensions
    for key, reduction in cds.reducedDims.items():
        if reduction.shape[0] != n_cells:
            return False, f"Reduction '{key}' shape mismatch: {reduction.shape[0]} vs {n_cells}"

    # Check clusters
    if cds.clusters is not None and len(cds.clusters) != n_cells:
        return False, f"Clusters length mismatch: {len(cds.clusters)} vs {n_cells}"

    return True, "Valid"


def estimate_memory_usage(cds) -> Dict[str, float]:
    """
    Estimate memory usage of CellDataSet

    Parameters
    ----------
    cds : CellDataSet
        CellDataSet to analyze

    Returns
    -------
    dict
        Memory usage in MB for each component
    """
    def get_size(obj):
        """Get size of object in MB"""
        if isinstance(obj, np.ndarray):
            return obj.nbytes / 1024 / 1024
        elif sp.issparse(obj):
            return (obj.data.nbytes + obj.indices.nbytes + obj.indptr.nbytes) / 1024 / 1024
        elif isinstance(obj, pd.DataFrame):
            return obj.memory_usage(deep=True).sum() / 1024 / 1024
        else:
            return 0

    usage = {
        'expression': get_size(cds.expression_data),
        'cell_metadata': get_size(cds.cell_metadata),
        'gene_metadata': get_size(cds.gene_metadata),
    }

    for key, reduction in cds.reducedDims.items():
        usage[f'reduction_{key}'] = get_size(reduction)

    usage['total'] = sum(usage.values())

    return usage
