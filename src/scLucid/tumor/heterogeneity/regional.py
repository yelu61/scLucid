"""
Regional heterogeneity and spatial pattern analysis.

This module provides tools for analyzing regional heterogeneity
in tumor samples with or without spatial information.
"""

import logging
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData

log = logging.getLogger(__name__)


class RegionalAnalyzer:
    """
    Analyze regional heterogeneity in tumors.

    Parameters
    ----------
    spatial_key : str
        Key for spatial coordinates in adata.obsm

    Attributes:
    ----------
    regions_ : pd.DataFrame
        Regional heterogeneity metrics
    """

    def __init__(
        self,
        spatial_key: str = "spatial",
    ):
        self.spatial_key = spatial_key
        self.regions_: Optional[pd.DataFrame] = None

    def analyze_regional_heterogeneity(
        self,
        adata: AnnData,
        clone_key: str = "clone_id",
        region_key: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Analyze heterogeneity across tumor regions.

        Parameters
        ----------
        adata : AnnData
            Expression data
        clone_key : str
            Column containing clone IDs
        region_key : str, optional
            Column defining regions. If None, will infer from spatial coordinates

        Returns:
        -------
        pd.DataFrame
            Regional heterogeneity metrics
        """
        if region_key is None and self.spatial_key in adata.obsm:
            # Infer regions from spatial coordinates
            regions = self._cluster_spatial_regions(adata)
        elif region_key is not None:
            regions = adata.obs[region_key]
        else:
            raise ValueError("Either region_key must be provided or spatial coordinates must exist")

        results = []

        for region in regions.unique():
            mask = regions == region
            region_adata = adata[mask]

            # Calculate clonal composition
            clone_counts = region_adata.obs[clone_key].value_counts()
            clone_props = clone_counts / clone_counts.sum()

            result = {
                "region": region,
                "n_cells": region_adata.n_obs,
                "n_clones": len(clone_counts),
                "shannon_diversity": -np.sum(clone_props * np.log(clone_props + 1e-10)),
                "simpson_diversity": 1 - np.sum(clone_props**2),
                "dominant_clone": clone_counts.index[0],
                "dominant_clone_freq": clone_props.iloc[0],
            }

            results.append(result)

        self.regions_ = pd.DataFrame(results)
        return self.regions_

    def _cluster_spatial_regions(
        self,
        adata: AnnData,
        n_regions: int = 4,
    ) -> pd.Series:
        """Cluster spatial coordinates into regions."""
        from sklearn.cluster import KMeans

        coords = adata.obsm[self.spatial_key]

        kmeans = KMeans(n_clusters=n_regions, random_state=42, n_init=10)
        labels = kmeans.fit_predict(coords)

        return pd.Series(labels, index=adata.obs_names, name="region")

    def identify_spatial_patterns(
        self,
        adata: AnnData,
        clone_key: str = "clone_id",
        method: str = "moran",
    ) -> pd.DataFrame:
        """
        Identify spatial patterns of clones.

        Parameters
        ----------
        adata : AnnData
            Expression data with spatial coordinates
        clone_key : str
            Column containing clone IDs
        method : str
            Spatial autocorrelation method ("moran", "geary")

        Returns:
        -------
        pd.DataFrame
            Spatial pattern metrics
        """
        if self.spatial_key not in adata.obsm:
            raise ValueError(f"Spatial coordinates not found in obsm['{self.spatial_key}']")

        results = []

        # For each clone, calculate spatial autocorrelation
        for clone in adata.obs[clone_key].unique():
            # Binary indicator for clone membership
            is_clone = (adata.obs[clone_key] == clone).astype(int).values

            if method == "moran":
                score = self._calculate_moran_i(adata, is_clone)
            elif method == "geary":
                score = self._calculate_geary_c(adata, is_clone)
            else:
                raise ValueError(f"Unknown method: {method}")

            results.append(
                {
                    "clone": clone,
                    "spatial_score": score,
                    "n_cells": is_clone.sum(),
                    "method": method,
                }
            )

        return pd.DataFrame(results)

    def _calculate_moran_i(
        self,
        adata: AnnData,
        values: np.ndarray,
        k_neighbors: int = 10,
    ) -> float:
        """
        Calculate Moran's I spatial autocorrelation.
        """
        from sklearn.neighbors import NearestNeighbors

        coords = adata.obsm[self.spatial_key]

        # Build nearest neighbor graph
        nbrs = NearestNeighbors(n_neighbors=k_neighbors + 1).fit(coords)
        distances, indices = nbrs.kneighbors(coords)

        # Calculate Moran's I
        n = len(values)
        z = values - values.mean()

        # Spatial weights (inverse distance)
        W = np.zeros((n, n))
        for i in range(n):
            for j_idx, j in enumerate(indices[i][1:]):  # Skip self
                W[i, j] = 1.0 / (distances[i][j_idx + 1] + 1e-6)

        W_sum = W.sum()

        numerator = 0
        for i in range(n):
            for j in range(n):
                numerator += W[i, j] * z[i] * z[j]

        denominator = (z**2).sum()

        if denominator == 0:
            return 0.0

        I = (n / W_sum) * (numerator / denominator)
        return I

    def _calculate_geary_c(
        self,
        adata: AnnData,
        values: np.ndarray,
        k_neighbors: int = 10,
    ) -> float:
        """
        Calculate Geary's C spatial autocorrelation.
        """
        from sklearn.neighbors import NearestNeighbors

        coords = adata.obsm[self.spatial_key]

        # Build nearest neighbor graph
        nbrs = NearestNeighbors(n_neighbors=k_neighbors + 1).fit(coords)
        distances, indices = nbrs.kneighbors(coords)

        n = len(values)

        # Spatial weights
        W = np.zeros((n, n))
        for i in range(n):
            for j_idx, j in enumerate(indices[i][1:]):
                W[i, j] = 1.0 / (distances[i][j_idx + 1] + 1e-6)

        W_sum = W.sum()

        numerator = 0
        for i in range(n):
            for j in range(n):
                numerator += W[i, j] * (values[i] - values[j]) ** 2

        denominator = 2 * W_sum * ((values - values.mean()) ** 2).sum() / (n - 1)

        if denominator == 0:
            return 0.0

        C = ((n - 1) / (2 * W_sum)) * (numerator / denominator)
        return C


def analyze_regional_heterogeneity(
    adata: AnnData,
    clone_key: str = "clone_id",
    region_key: Optional[str] = None,
    spatial_key: str = "spatial",
) -> pd.DataFrame:
    """
    Analyze heterogeneity across tumor regions.

    Parameters
    ----------
    adata : AnnData
        Expression data
    clone_key : str
        Column containing clone IDs
    region_key : str, optional
        Column defining regions
    spatial_key : str
        Key for spatial coordinates

    Returns:
    -------
    pd.DataFrame
        Regional heterogeneity metrics
    """
    analyzer = RegionalAnalyzer(spatial_key=spatial_key)
    return analyzer.analyze_regional_heterogeneity(adata, clone_key, region_key)


def identify_spatial_patterns(
    adata: AnnData,
    clone_key: str = "clone_id",
    method: str = "moran",
    spatial_key: str = "spatial",
) -> pd.DataFrame:
    """
    Identify spatial patterns of clones.

    Parameters
    ----------
    adata : AnnData
        Expression data with spatial coordinates
    clone_key : str
        Column containing clone IDs
    method : str
        Spatial autocorrelation method
    spatial_key : str
        Key for spatial coordinates

    Returns:
    -------
    pd.DataFrame
        Spatial pattern metrics
    """
    analyzer = RegionalAnalyzer(spatial_key=spatial_key)
    return analyzer.identify_spatial_patterns(adata, clone_key, method)


def calculate_regional_expression_differences(
    adata: AnnData,
    region_key: str,
    method: str = "wilcoxon",
) -> pd.DataFrame:
    """
    Calculate expression differences between regions.

    Parameters
    ----------
    adata : AnnData
        Expression data
    region_key : str
        Column defining regions
    method : str
        Statistical test method

    Returns:
    -------
    pd.DataFrame
        Differential expression results
    """
    from scipy import stats

    regions = adata.obs[region_key].unique()

    if len(regions) != 2:
        raise ValueError("This function requires exactly 2 regions")

    region1_mask = adata.obs[region_key] == regions[0]
    region2_mask = adata.obs[region_key] == regions[1]

    results = []

    for gene in adata.var_names:
        expr1 = adata[region1_mask, gene].X
        expr2 = adata[region2_mask, gene].X

        if hasattr(expr1, "toarray"):
            expr1 = expr1.toarray().flatten()
        if hasattr(expr2, "toarray"):
            expr2 = expr2.toarray().flatten()

        if method == "wilcoxon":
            try:
                stat, pval = stats.ranksums(expr1, expr2)
            except:
                continue
        elif method == "t-test":
            stat, pval = stats.ttest_ind(expr1, expr2)
        else:
            raise ValueError(f"Unknown method: {method}")

        results.append(
            {
                "gene": gene,
                "region1_mean": np.mean(expr1),
                "region2_mean": np.mean(expr2),
                "log2fc": np.log2((np.mean(expr1) + 1e-6) / (np.mean(expr2) + 1e-6)),
                "statistic": stat,
                "pvalue": pval,
            }
        )

    results_df = pd.DataFrame(results)

    # Add FDR correction
    from scipy.stats import false_discovery_rate

    if len(results_df) > 0:
        results_df["fdr"] = false_discovery_rate(results_df["pvalue"].fillna(1))[1]

    return results_df.sort_values("pvalue")
