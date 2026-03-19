"""
Dimensionality reduction for pyMonocle3 (R-free)
"""

import numpy as np
import scipy.sparse as sp
from typing import Optional, Union, List
import logging

from .core import CellDataSet

log = logging.getLogger(__name__)


def reduce_dimension(
    cds: CellDataSet,
    reduction_method: str = "UMAP",
    preprocess_method: str = "PCA",
    umap_metric: str = "cosine",
    umap_min_dist: float = 0.1,
    umap_n_neighbors: int = 15,
    umap_fast_sgd: bool = True,
    cores: int = 1,
    verbose: bool = False,
) -> CellDataSet:
    """
    Perform dimensionality reduction on CellDataSet

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    reduction_method : str
        Reduction method ("UMAP", "tSNE", or "PCA")
    preprocess_method : str
        Preprocessing reduction to use as input
    umap_metric : str
        Distance metric for UMAP
    umap_min_dist : float
        Minimum distance for UMAP
    umap_n_neighbors : int
        Number of neighbors for UMAP
    umap_fast_sgd : bool
        Use fast SGD for UMAP
    cores : int
        Number of cores to use
    verbose : bool
        Verbose output

    Returns
    -------
    CellDataSet
        CellDataSet with reduced dimensions added
    """
    # Check if preprocessing reduction exists
    if preprocess_method not in cds.reducedDims:
        raise ValueError(
            f"{preprocess_method} not found in reducedDims. "
            f"Run preprocess_cds first."
        )

    reduced_data = cds.reducedDims[preprocess_method]

    if reduction_method == "UMAP":
        try:
            import umap
        except ImportError:
            raise ImportError(
                "UMAP is required for dimensionality reduction. "
                "Install with: pip install umap-learn"
            )

        reducer = umap.UMAP(
            n_neighbors=umap_n_neighbors,
            min_dist=umap_min_dist,
            metric=umap_metric,
            n_components=2,
            random_state=42,
        )

        embedding = reducer.fit_transform(reduced_data)
        cds.reducedDims['UMAP'] = embedding

        log.info(f"UMAP completed: {embedding.shape}")

    elif reduction_method == "tSNE":
        from sklearn.manifold import TSNE

        tsne = TSNE(
            n_components=2,
            perplexity=30,
            learning_rate='auto',
            init='pca',
            random_state=42,
            n_jobs=cores,
            verbose=verbose,
        )

        embedding = tsne.fit_transform(reduced_data)
        cds.reducedDims['tSNE'] = embedding

        log.info(f"tSNE completed: {embedding.shape}")

    elif reduction_method == "PCA":
        # Already done in preprocessing, just pass through
        cds.reducedDims['PCA_viz'] = reduced_data[:, :2] if reduced_data.shape[1] > 2 else reduced_data

        log.info(f"PCA visualization dimensions ready")

    else:
        raise ValueError(f"Unknown reduction_method: {reduction_method}")

    return cds


def run_pca(
    cds: CellDataSet,
    n_components: int = 50,
    use_genes: Optional[List[str]] = None,
    layer: Optional[str] = None,
) -> CellDataSet:
    """
    Run PCA on CellDataSet

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    n_components : int
        Number of principal components
    use_genes : list, optional
        Specific genes to use
    layer : str, optional
        Expression layer to use

    Returns
    -------
    CellDataSet
        CellDataSet with PCA results
    """
    from sklearn.decomposition import PCA

    # Get expression data
    expr = cds.expression_data
    if sp.issparse(expr):
        expr = expr.toarray()

    # Subset genes if specified
    if use_genes is not None:
        gene_mask = cds.gene_metadata.index.isin(use_genes)
        expr = expr[gene_mask, :]

    # Transpose for PCA (cells x genes)
    expr_t = expr.T

    # Run PCA
    n_components = min(n_components, expr_t.shape[0], expr_t.shape[1])
    pca = PCA(n_components=n_components, random_state=42)

    reduced = pca.fit_transform(expr_t)

    cds.reducedDims['PCA'] = reduced
    cds.preprocessing_params['pca_variance_explained'] = pca.explained_variance_ratio_

    log.info(f"PCA: {reduced.shape[1]} components, "
             f"explained variance: {np.sum(pca.explained_variance_ratio_):.2%}")

    return cds


def run_umap(
    cds: CellDataSet,
    reduction_key: str = "PCA",
    n_components: int = 2,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    metric: str = "euclidean",
    random_state: int = 42,
) -> CellDataSet:
    """
    Run UMAP on existing reduction

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    reduction_key : str
        Key of reduction to use as input
    n_components : int
        Number of output dimensions
    n_neighbors : int
        Number of neighbors
    min_dist : float
        Minimum distance between points
    metric : str
        Distance metric
    random_state : int
        Random seed

    Returns
    -------
    CellDataSet
        CellDataSet with UMAP results
    """
    try:
        import umap
    except ImportError:
        raise ImportError(
            "UMAP is required. Install with: pip install umap-learn"
        )

    if reduction_key not in cds.reducedDims:
        raise ValueError(f"Reduction '{reduction_key}' not found. "
                        f"Available: {list(cds.reducedDims.keys())}")

    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=n_components,
        metric=metric,
        random_state=random_state,
    )

    embedding = reducer.fit_transform(cds.reducedDims[reduction_key])

    cds.reducedDims[f'UMAP_{n_components}D'] = embedding

    log.info(f"UMAP {n_components}D completed: {embedding.shape}")

    return cds
