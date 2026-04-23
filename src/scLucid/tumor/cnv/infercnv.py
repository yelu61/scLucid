"""
Infer CNV from single-cell RNA-seq data.

This module implements CNV inference methods based on gene expression patterns.
It identifies tumor cells by detecting chromosomal amplifications and deletions.
"""

import logging
from typing import List, Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.ndimage import gaussian_filter1d

log = logging.getLogger(__name__)


class CNVAnalyzer:
    """
    Analyze copy number variations from scRNA-seq data.

    This class provides methods to infer CNV status from gene expression
    and identify malignant cells based on CNV patterns.

    Parameters
    ----------
    gene_order : pd.DataFrame
        DataFrame with gene positions (chromosome, start, end)
    window_size : int
        Size of smoothing window for CNV inference

    Attributes:
    ----------
    cnv_matrix_ : np.ndarray
        Inferred CNV matrix (cells x genomic_windows)
    tumor_scores_ : pd.Series
        Tumor malignancy scores per cell
    """

    def __init__(
        self,
        gene_order: Optional[pd.DataFrame] = None,
        window_size: int = 100,
    ):
        self.gene_order = gene_order
        self.window_size = window_size
        self.cnv_matrix_: Optional[np.ndarray] = None
        self.tumor_scores_: Optional[pd.Series] = None

    def fit(
        self,
        adata: AnnData,
        reference_cells: Optional[Union[str, List[str]]] = None,
        reference_key: str = "cell_type",
    ) -> "CNVAnalyzer":
        """
        Infer CNV from expression data.

        Parameters
        ----------
        adata : AnnData
            Single-cell expression data
        reference_cells : str or list, optional
            Cell types to use as reference (normal cells)
        reference_key : str
            Column in adata.obs containing cell type labels

        Returns:
        -------
        CNVAnalyzer
            Fitted analyzer
        """
        # Get expression matrix
        X = adata.X if not hasattr(adata.X, "toarray") else adata.X.toarray()

        # Get reference cells
        if reference_cells is not None:
            if isinstance(reference_cells, str):
                ref_mask = adata.obs[reference_key] == reference_cells
            else:
                ref_mask = adata.obs[reference_key].isin(reference_cells)
            ref_expr = X[ref_mask].mean(axis=0)
        else:
            # Use all cells as reference (less optimal)
            ref_expr = X.mean(axis=0)

        # Calculate relative expression (log2 ratio)
        cnv_ratio = np.log2((X + 1) / (ref_expr + 1))

        # Smooth along chromosomes if gene order provided
        if self.gene_order is not None:
            cnv_ratio = self._smooth_by_chromosome(cnv_ratio, adata.var_names)

        self.cnv_matrix_ = cnv_ratio

        # Calculate tumor scores
        self.tumor_scores_ = self._calculate_tumor_scores(cnv_ratio)

        return self

    def _smooth_by_chromosome(
        self,
        cnv_ratio: np.ndarray,
        gene_names: pd.Index,
    ) -> np.ndarray:
        """Smooth CNV signals along chromosomes."""
        smoothed = np.zeros_like(cnv_ratio)

        for chrom in self.gene_order["chromosome"].unique():
            chrom_genes = self.gene_order[self.gene_order["chromosome"] == chrom].index

            # Find indices of chromosomal genes
            mask = gene_names.isin(chrom_genes)
            if mask.sum() > 0:
                chrom_data = cnv_ratio[:, mask]
                # Apply Gaussian smoothing
                smoothed[:, mask] = gaussian_filter1d(chrom_data, sigma=self.window_size, axis=1)

        return smoothed

    def _calculate_tumor_scores(self, cnv_ratio: np.ndarray) -> pd.Series:
        """Calculate tumor malignancy scores from CNV patterns."""
        # Score based on deviation from normal
        score = np.abs(cnv_ratio).mean(axis=1)
        return pd.Series(score)

    def predict_tumor_cells(
        self,
        threshold: float = 0.5,
    ) -> pd.Series:
        """
        Predict tumor cells based on CNV scores.

        Parameters
        ----------
        threshold : float
            Threshold for tumor classification

        Returns:
        -------
        pd.Series
            Boolean series indicating tumor cells
        """
        if self.tumor_scores_ is None:
            raise ValueError("Must call fit() first")

        return self.tumor_scores_ > threshold


def infer_cnv(
    adata: AnnData,
    reference_cells: Optional[Union[str, List[str]]] = None,
    reference_key: str = "cell_type",
    gene_order: Optional[pd.DataFrame] = None,
    key_added: str = "cnv",
    copy: bool = False,
) -> AnnData:
    """
    Infer CNV from single-cell expression data.

    Parameters
    ----------
    adata : AnnData
        Single-cell expression data
    reference_cells : str or list, optional
        Cell types to use as normal reference
    reference_key : str
        Column in adata.obs with cell type labels
    gene_order : pd.DataFrame, optional
        Gene position information
    key_added : str
        Key for storing results in adata.obsm
    copy : bool
        Return a copy of adata

    Returns:
    -------
    AnnData
        Annotated data with CNV information
    """
    if copy:
        adata = adata.copy()

    analyzer = CNVAnalyzer(gene_order=gene_order)
    analyzer.fit(adata, reference_cells, reference_key)

    # Store results
    adata.obsm[f"X_{key_added}"] = analyzer.cnv_matrix_
    adata.obs[f"{key_added}_score"] = analyzer.tumor_scores_

    log.info(f"CNV inference complete. Results stored in obsm['X_{key_added}']")

    return adata


def find_tumor_cells(
    adata: AnnData,
    method: str = "cnv_score",
    threshold: float = 0.5,
    key: str = "cnv",
) -> pd.Series:
    """
    Identify tumor cells based on CNV analysis.

    Parameters
    ----------
    adata : AnnData
        Data with CNV information
    method : str
        Method for tumor identification ("cnv_score", "threshold", "clustering")
    threshold : float
        Threshold for classification
    key : str
        Key for CNV data in adata

    Returns:
    -------
    pd.Series
        Boolean series indicating tumor cells
    """
    if method == "cnv_score":
        scores = adata.obs[f"{key}_score"]
        return scores > threshold

    elif method == "clustering":
        from sklearn.cluster import KMeans

        cnv_data = adata.obsm[f"X_{key}"]
        kmeans = KMeans(n_clusters=2, random_state=42).fit(cnv_data)

        # Assign tumor label to cluster with higher CNV burden
        cluster_means = [np.abs(cnv_data[kmeans.labels_ == i]).mean() for i in range(2)]
        tumor_cluster = np.argmax(cluster_means)

        return pd.Series(kmeans.labels_ == tumor_cluster, index=adata.obs_names)

    else:
        raise ValueError(f"Unknown method: {method}")


def identify_clones(
    adata: AnnData,
    cnv_key: str = "cnv",
    n_clusters: int = 5,
    method: str = "hierarchical",
) -> pd.Series:
    """
    Identify tumor clones based on CNV patterns.

    Parameters
    ----------
    adata : AnnData
        Data with CNV information
    cnv_key : str
        Key for CNV data
    n_clusters : int
        Number of expected clones
    method : str
        Clustering method ("hierarchical", "kmeans", "leiden")

    Returns:
    -------
    pd.Series
        Clone assignments for each cell
    """
    cnv_data = adata.obsm[f"X_{cnv_key}"]

    if method == "hierarchical":
        from sklearn.cluster import AgglomerativeClustering

        clusterer = AgglomerativeClustering(n_clusters=n_clusters)
    elif method == "kmeans":
        from sklearn.cluster import KMeans

        clusterer = KMeans(n_clusters=n_clusters, random_state=42)
    elif method == "leiden":
        # Use scanpy's leiden implementation
        import scanpy as sc

        adata_cnv = AnnData(X=cnv_data)
        sc.pp.neighbors(adata_cnv, n_neighbors=15)
        sc.tl.leiden(adata_cnv, resolution=1.0)
        return adata_cnv.obs["leiden"]
    else:
        raise ValueError(f"Unknown method: {method}")

    labels = clusterer.fit_predict(cnv_data)
    return pd.Series(labels, index=adata.obs_names, name="clone")


def calculate_cnv_score(
    adata: AnnData,
    cnv_key: str = "cnv",
    method: str = "mean_absolute",
) -> pd.Series:
    """
    Calculate overall CNV burden score.

    Parameters
    ----------
    adata : AnnData
        Data with CNV information
    cnv_key : str
        Key for CNV data
    method : str
        Scoring method ("mean_absolute", "variance", "gini")

    Returns:
    -------
    pd.Series
        CNV scores per cell
    """
    cnv_data = adata.obsm[f"X_{cnv_key}"]

    if method == "mean_absolute":
        scores = np.abs(cnv_data).mean(axis=1)
    elif method == "variance":
        scores = np.var(cnv_data, axis=1)
    elif method == "gini":
        scores = np.array([_gini_coefficient(np.abs(x)) for x in cnv_data])
    else:
        raise ValueError(f"Unknown method: {method}")

    return pd.Series(scores, index=adata.obs_names, name="cnv_burden")


def _gini_coefficient(x: np.ndarray) -> float:
    """Calculate Gini coefficient for inequality."""
    sorted_x = np.sort(x)
    n = len(x)
    cumsum = np.cumsum(sorted_x)
    return (n + 1 - 2 * np.sum(cumsum) / cumsum[-1]) / n
