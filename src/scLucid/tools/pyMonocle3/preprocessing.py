"""
Preprocessing functions for pyMonocle3 (R-free)
"""

import logging
from typing import List, Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.preprocessing import scale

from .core import CellDataSet

log = logging.getLogger(__name__)


def detect_genes(
    cds: CellDataSet,
    min_expr: float = 0.1,
    min_cells: int = 10,
) -> CellDataSet:
    """
    Detect expressed genes and calculate gene statistics

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    min_expr : float
        Minimum expression threshold
    min_cells : int
        Minimum number of cells expressing the gene

    Returns:
    -------
    CellDataSet
        Updated CellDataSet with gene statistics
    """
    expr = cds.expression_data

    # Calculate gene statistics
    if sp.issparse(expr):
        expr_dense = expr.toarray()
    else:
        expr_dense = expr

    # Mean expression per gene
    mean_expr = np.mean(expr_dense, axis=1)

    # Proportion of cells expressing the gene
    prop_expr = np.mean(expr_dense > min_expr, axis=1)

    # Number of cells expressing the gene
    num_cells = np.sum(expr_dense > min_expr, axis=1)

    # Update gene metadata
    cds.gene_metadata["mean_expression"] = mean_expr
    cds.gene_metadata["prop_expressed"] = prop_expr
    cds.gene_metadata["num_cells_expressed"] = num_cells

    # Mark valid genes
    cds.gene_metadata["use_for_ordering"] = num_cells >= min_cells

    log.info(f"Detected {np.sum(num_cells >= min_cells)} genes expressed in >= {min_cells} cells")

    return cds


def estimate_size_factors(
    cds: CellDataSet,
    locfunc: callable = np.median,
    round_exprs: bool = True,
    method: str = "mean-geometric-mean-total",
) -> CellDataSet:
    """
    Estimate size factors for normalization

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    locfunc : callable
        Location function for size factor estimation
    round_exprs : bool
        Whether to round expression values
    method : str
        Normalization method

    Returns:
    -------
    CellDataSet
        CellDataSet with size factors added to cell_metadata
    """
    expr = cds.expression_data

    if sp.issparse(expr):
        expr = expr.toarray()

    if round_exprs:
        expr = np.round(expr)

    # Calculate total counts per cell
    total_counts = np.sum(expr, axis=0)

    if method == "mean-geometric-mean-total":
        # Calculate geometric mean of totals
        geometric_mean = np.exp(np.mean(np.log(total_counts + 1))) - 1
        size_factors = total_counts / geometric_mean
    elif method == "median":
        median_total = locfunc(total_counts)
        size_factors = total_counts / median_total
    else:
        raise ValueError(f"Unknown method: {method}")

    cds.cell_metadata["Size_Factor"] = size_factors

    log.info(
        f"Size factors estimated: mean={np.mean(size_factors):.3f}, std={np.std(size_factors):.3f}"
    )

    return cds


def preprocess_cds(
    cds: CellDataSet,
    num_dim: int = 50,
    norm_method: str = "log",
    pseudo_count: float = 1.0,
    scaling: bool = True,
    method: str = "PCA",
    use_genes: Optional[List[str]] = None,
) -> CellDataSet:
    """
    Preprocess CellDataSet: normalize, reduce dimensions

    This is the main preprocessing function that combines normalization,
    scaling, and dimensionality reduction.

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    num_dim : int
        Number of dimensions for reduction
    norm_method : str
        Normalization method ("log" or "size_only")
    pseudo_count : float
        Pseudo-count for log normalization
    scaling : bool
        Whether to scale genes to unit variance
    method : str
        Reduction method ("PCA" or "LSI")
    use_genes : list, optional
        Specific genes to use

    Returns:
    -------
    CellDataSet
        Preprocessed CellDataSet with reduced dimensions
    """
    log.info("Starting preprocessing...")

    # Step 1: Detect genes
    cds = detect_genes(cds)

    # Step 2: Estimate size factors
    cds = estimate_size_factors(cds)

    # Step 3: Normalize expression
    expr = cds.expression_data
    if sp.issparse(expr):
        expr = expr.toarray()

    size_factors = cds.cell_metadata["Size_Factor"].values

    if norm_method == "log":
        # Size factor normalization + log transform
        norm_expr = expr / size_factors[np.newaxis, :]
        norm_expr = np.log1p(norm_expr + pseudo_count)
    elif norm_method == "size_only":
        norm_expr = expr / size_factors[np.newaxis, :]
    else:
        raise ValueError(f"Unknown norm_method: {norm_method}")

    # Step 4: Select genes for ordering
    if use_genes is not None:
        gene_mask = cds.gene_metadata.index.isin(use_genes)
    else:
        gene_mask = cds.gene_metadata["use_for_ordering"].values

    norm_expr_subset = norm_expr[gene_mask, :]

    # Step 5: Scale if requested
    if scaling:
        norm_expr_subset = scale(norm_expr_subset, axis=1, with_mean=True, with_std=True)
        # Replace NaN with 0
        norm_expr_subset = np.nan_to_num(norm_expr_subset, nan=0.0)

    # Step 6: Dimensionality reduction
    if method == "PCA":
        from sklearn.decomposition import PCA

        # Transpose for PCA (samples x features)
        pca = PCA(n_components=min(num_dim, norm_expr_subset.shape[0], norm_expr_subset.shape[1]))
        reduced = pca.fit_transform(norm_expr_subset.T)

        cds.reducedDims["PCA"] = reduced

        # Store variance explained
        cds.preprocessing_params["pca_variance_explained"] = pca.explained_variance_ratio_

        log.info(f"PCA completed: {reduced.shape[1]} components")

    elif method == "LSI":
        # Latent Semantic Indexing (common for ATAC-seq)
        from sklearn.decomposition import TruncatedSVD

        svd = TruncatedSVD(n_components=min(num_dim, norm_expr_subset.shape[0] - 1))

        # For LSI, use TF-IDF like transformation
        tfidf = norm_expr_subset * np.log1p(
            norm_expr_subset.shape[1] / (1 + np.sum(norm_expr_subset > 0, axis=1)[:, np.newaxis])
        )

        reduced = svd.fit_transform(tfidf.T)
        cds.reducedDims["LSI"] = reduced

        log.info(f"LSI completed: {reduced.shape[1]} components")

    else:
        raise ValueError(f"Unknown method: {method}")

    # Store preprocessing parameters
    cds.preprocessing_params.update(
        {
            "num_dim": num_dim,
            "norm_method": norm_method,
            "pseudo_count": pseudo_count,
            "scaling": scaling,
            "method": method,
        }
    )

    log.info("Preprocessing completed")

    return cds


def align_cds(
    cds_list: List[CellDataSet],
    method: str = "mutual_nearest_neighbor",
    k: int = 20,
) -> CellDataSet:
    """
    Align multiple CellDataSets (batch correction)

    Parameters
    ----------
    cds_list : list
        List of CellDataSets to align
    method : str
        Alignment method
    k : int
        Number of nearest neighbors

    Returns:
    -------
    CellDataSet
        Aligned CellDataSet
    """
    if len(cds_list) < 2:
        raise ValueError("Need at least 2 CellDataSets to align")

    if method == "mutual_nearest_neighbor":
        # Simple MNN alignment
        from sklearn.neighbors import NearestNeighbors

        # Get PCA representations
        refs = [cds.reducedDims["PCA"] for cds in cds_list]

        # Find mutual nearest neighbors
        # This is a simplified version - full MNN is more complex
        aligned = refs[0].copy()

        for i, ref in enumerate(refs[1:], 1):
            # Find nearest neighbors between datasets
            nn = NearestNeighbors(n_neighbors=k)
            nn.fit(ref)
            distances, indices = nn.kneighbors(aligned)

            # Compute correction vector
            correction = np.mean(ref[indices] - aligned[:, np.newaxis, :], axis=1)

            # Apply correction
            aligned = aligned + correction

        # Create merged CellDataSet
        merged_expr = sp.hstack([cds.expression_data for cds in cds_list])
        merged_meta = pd.concat([cds.cell_metadata for cds in cds_list], ignore_index=True)
        merged_meta["batch"] = np.concatenate([[i] * cds.n_cells for i, cds in enumerate(cds_list)])
        merged_gene_meta = cds_list[0].gene_metadata.copy()

        cds_aligned = CellDataSet(
            expression_data=merged_expr,
            cell_metadata=merged_meta,
            gene_metadata=merged_gene_meta,
        )
        cds_aligned.reducedDims["PCA_aligned"] = aligned

        log.info(f"Alignment completed: {cds_aligned.n_cells} cells from {len(cds_list)} batches")

        return cds_aligned

    else:
        raise ValueError(f"Unknown alignment method: {method}")
