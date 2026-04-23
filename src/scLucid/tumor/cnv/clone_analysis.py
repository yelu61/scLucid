"""
Clone analysis and phylogenetics from CNV data.

This module provides tools for analyzing clonal structure
and evolution from copy number variation data.
"""

import logging
from typing import Dict, Optional

import numpy as np
import pandas as pd
from anndata import AnnData

log = logging.getLogger(__name__)


class CloneAnalyzer:
    """
    Analyze clonal structure from CNV data.

    Parameters
    ----------
    n_clusters : int
        Number of clones to identify
    method : str
        Clustering method ("kmeans", "hierarchical", "leiden")

    Attributes:
    ----------
    clone_ids_ : pd.Series
        Clone assignments for each cell
    clone_profiles_ : pd.DataFrame
        Average CNV profile per clone
    """

    def __init__(
        self,
        n_clusters: int = 5,
        method: str = "kmeans",
    ):
        self.n_clusters = n_clusters
        self.method = method
        self.clone_ids_: Optional[pd.Series] = None
        self.clone_profiles_: Optional[pd.DataFrame] = None

    def fit(self, adata: AnnData, cnv_key: str = "cnv") -> "CloneAnalyzer":
        """
        Identify clones from CNV data.

        Parameters
        ----------
        adata : AnnData
            Expression data with CNV
        cnv_key : str
            Key for CNV scores in adata.obsm

        Returns:
        -------
        CloneAnalyzer
            Fitted analyzer
        """
        if cnv_key in adata.obsm:
            cnv_data = adata.obsm[cnv_key]
        elif f"{cnv_key}_score" in adata.obs.columns:
            cnv_data = adata.obs[[f"{cnv_key}_score"]].values
        else:
            raise ValueError(f"CNV data not found with key {cnv_key}")

        # Cluster cells
        if self.method == "kmeans":
            from sklearn.cluster import KMeans

            clusterer = KMeans(n_clusters=self.n_clusters, random_state=42, n_init=10)
            labels = clusterer.fit_predict(cnv_data)
            self.centroids_ = clusterer.cluster_centers_

        elif self.method == "hierarchical":
            from sklearn.cluster import AgglomerativeClustering

            clusterer = AgglomerativeClustering(n_clusters=self.n_clusters)
            labels = clusterer.fit_predict(cnv_data)

        else:
            raise ValueError(f"Unknown method: {self.method}")

        self.clone_ids_ = pd.Series(labels, index=adata.obs_names, name="clone_id")

        # Calculate clone profiles
        self.clone_profiles_ = self._calculate_clone_profiles(cnv_data, labels)

        return self

    def _calculate_clone_profiles(
        self,
        cnv_data: np.ndarray,
        labels: np.ndarray,
    ) -> pd.DataFrame:
        """Calculate average CNV profile per clone."""
        profiles = {}

        for clone_id in np.unique(labels):
            mask = labels == clone_id
            profiles[f"clone_{clone_id}"] = cnv_data[mask].mean(axis=0)

        return pd.DataFrame(profiles)

    def calculate_clonal_diversity(self) -> Dict[str, float]:
        """
        Calculate clonal diversity metrics.

        Returns:
        -------
        dict
            Diversity metrics
        """
        if self.clone_ids_ is None:
            raise ValueError("CloneAnalyzer not fitted yet")

        counts = self.clone_ids_.value_counts()
        proportions = counts / counts.sum()

        # Shannon diversity
        shannon = -np.sum(proportions * np.log(proportions + 1e-10))

        # Simpson diversity
        simpson = 1 - np.sum(proportions**2)

        # Richness
        richness = len(proportions)

        return {
            "shannon_diversity": shannon,
            "simpson_diversity": simpson,
            "richness": richness,
            "dominant_clone_freq": proportions.max(),
        }

    def compare_clones(self, adata: AnnData) -> pd.DataFrame:
        """
        Compare clones by expression differences.

        Parameters
        ----------
        adata : AnnData
            Expression data

        Returns:
        -------
        pd.DataFrame
            Clone comparison results
        """
        if self.clone_ids_ is None:
            raise ValueError("CloneAnalyzer not fitted yet")

        adata.obs["clone_id"] = self.clone_ids_

        results = []
        clones = self.clone_ids_.unique()

        for i, clone1 in enumerate(clones):
            for clone2 in clones[i + 1 :]:
                mask1 = self.clone_ids_ == clone1
                mask2 = self.clone_ids_ == clone2

                # Compare key markers
                for gene in adata.var_names[:100]:  # Limit to first 100 for speed
                    expr1 = adata[mask1, gene].X.mean()
                    expr2 = adata[mask2, gene].X.mean()

                    if hasattr(expr1, "toarray"):
                        expr1 = expr1.toarray().flatten()[0]
                    if hasattr(expr2, "toarray"):
                        expr2 = expr2.toarray().flatten()[0]

                    results.append(
                        {
                            "clone1": clone1,
                            "clone2": clone2,
                            "gene": gene,
                            "mean1": expr1,
                            "mean2": expr2,
                            "diff": expr1 - expr2,
                        }
                    )

        return pd.DataFrame(results)


def identify_clones(
    adata: AnnData,
    n_clusters: int = 5,
    cnv_key: str = "cnv",
    key_added: str = "clone_id",
) -> pd.Series:
    """
    Identify tumor clones from CNV data.

    Parameters
    ----------
    adata : AnnData
        Expression data with CNV
    n_clusters : int
        Number of clones
    cnv_key : str
        Key for CNV data
    key_added : str
        Key for storing clone IDs

    Returns:
    -------
    pd.Series
        Clone assignments
    """
    analyzer = CloneAnalyzer(n_clusters=n_clusters)
    analyzer.fit(adata, cnv_key)

    adata.obs[key_added] = analyzer.clone_ids_

    log.info(f"Identified {n_clusters} clones")

    return analyzer.clone_ids_


def calculate_clonal_diversity(
    adata: AnnData,
    clone_key: str = "clone_id",
) -> Dict[str, float]:
    """
    Calculate clonal diversity metrics.

    Parameters
    ----------
    adata : AnnData
        Expression data
    clone_key : str
        Column with clone IDs

    Returns:
    -------
    dict
        Diversity metrics
    """
    counts = adata.obs[clone_key].value_counts()
    proportions = counts / counts.sum()

    shannon = -np.sum(proportions * np.log(proportions + 1e-10))
    simpson = 1 - np.sum(proportions**2)

    return {
        "shannon_diversity": shannon,
        "simpson_diversity": simpson,
        "richness": len(proportions),
        "dominant_clone_freq": proportions.max(),
    }


def infer_clonal_phylogeny(
    adata: AnnData,
    clone_key: str = "clone_id",
    method: str = "nj",
) -> Dict:
    """
    Infer phylogenetic tree of clones.

    Parameters
    ----------
    adata : AnnData
        Expression data
    clone_key : str
        Column with clone IDs
    method : str
        Tree building method

    Returns:
    -------
    dict
        Phylogenetic tree
    """
    from ..evolution.phylogeny import build_phylogenetic_tree

    # Build tree from clone-level data
    tree = build_phylogenetic_tree(adata, clone_key=clone_key, method=method)

    return tree
