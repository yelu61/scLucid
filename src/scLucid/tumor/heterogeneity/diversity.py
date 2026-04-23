"""
Tumor diversity indices and heterogeneity quantification.

This module provides tools for calculating various diversity indices
to quantify intratumoral heterogeneity.
"""

import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial.distance import pdist

log = logging.getLogger(__name__)


def shannon_diversity_index(proportions: np.ndarray) -> float:
    """
    Calculate Shannon diversity index.

    H' = -sum(p_i * log(p_i))

    Parameters
    ----------
    proportions : np.ndarray
        Proportion of each clone/type

    Returns:
    -------
    float
        Shannon diversity index
    """
    proportions = proportions[proportions > 0]  # Remove zeros
    if len(proportions) == 0:
        return 0.0
    return -np.sum(proportions * np.log(proportions))


def simpson_diversity_index(proportions: np.ndarray) -> float:
    """
    Calculate Simpson's diversity index.

    D = 1 - sum(p_i^2)

    Parameters
    ----------
    proportions : np.ndarray
        Proportion of each clone/type

    Returns:
    -------
    float
        Simpson's diversity index
    """
    return 1 - np.sum(proportions**2)


def inverse_simpson_index(proportions: np.ndarray) -> float:
    """
    Calculate inverse Simpson's index.

    1/D = 1 / sum(p_i^2)

    Parameters
    ----------
    proportions : np.ndarray
        Proportion of each clone/type

    Returns:
    -------
    float
        Inverse Simpson's index
    """
    proportions = proportions[proportions > 0]
    if len(proportions) == 0:
        return 0.0
    return 1 / np.sum(proportions**2)


def gini_simpson_index(proportions: np.ndarray) -> float:
    """
    Calculate Gini-Simpson index.

    GS = 1 - sum(p_i^2)

    Parameters
    ----------
    proportions : np.ndarray
        Proportion of each clone/type

    Returns:
    -------
    float
        Gini-Simpson index
    """
    return 1 - np.sum(proportions**2)


def berger_parker_index(proportions: np.ndarray) -> float:
    """
    Calculate Berger-Parker dominance index.

    D = max(p_i)

    Parameters
    ----------
    proportions : np.ndarray
        Proportion of each clone/type

    Returns:
    -------
    float
        Berger-Parker index
    """
    return np.max(proportions)


def fisher_alpha(counts: np.ndarray) -> float:
    """
    Estimate Fisher's alpha diversity.

    Uses the relationship S = alpha * ln(1 + N/alpha)
    where S is species richness and N is total individuals.

    Parameters
    ----------
    counts : np.ndarray
        Count of each clone/type

    Returns:
    -------
    float
        Fisher's alpha
    """
    from scipy.optimize import fsolve

    S = np.sum(counts > 0)
    N = np.sum(counts)

    if S <= 1 or N <= 1:
        return 0.0

    def equation(alpha):
        return alpha * np.log(1 + N / alpha) - S

    try:
        alpha = fsolve(equation, S)[0]
        return max(0, alpha)
    except:
        return 0.0


class DiversityAnalyzer:
    """
    Analyze diversity and heterogeneity in tumor samples.

    Parameters
    ----------
    metrics : list
        List of diversity metrics to calculate

    Attributes:
    ----------
    diversity_scores_ : pd.DataFrame
        Diversity scores per sample
    """

    def __init__(
        self,
        metrics: Optional[List[str]] = None,
    ):
        self.metrics = metrics or [
            "shannon",
            "simpson",
            "inverse_simpson",
            "gini_simpson",
            "berger_parker",
            "fisher_alpha",
        ]
        self.diversity_scores_: Optional[pd.DataFrame] = None

    def calculate_diversity_indices(
        self,
        adata: AnnData,
        groupby: str,
        sample_key: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Calculate diversity indices for tumor samples.

        Parameters
        ----------
        adata : AnnData
            Expression data
        groupby : str
            Column defining groups (e.g., clones, cell types)
        sample_key : str, optional
            Column defining samples

        Returns:
        -------
        pd.DataFrame
            Diversity indices
        """
        results = []

        if sample_key is None:
            # Single sample analysis
            proportions = self._get_proportions(adata, groupby)
            result = self._calculate_metrics(proportions)
            result["sample"] = "all"
            results.append(result)
        else:
            # Multi-sample analysis
            for sample in adata.obs[sample_key].unique():
                mask = adata.obs[sample_key] == sample
                sample_adata = adata[mask]
                proportions = self._get_proportions(sample_adata, groupby)
                result = self._calculate_metrics(proportions)
                result["sample"] = sample
                results.append(result)

        self.diversity_scores_ = pd.DataFrame(results)
        return self.diversity_scores_

    def _get_proportions(
        self,
        adata: AnnData,
        groupby: str,
    ) -> np.ndarray:
        """Get proportions of each group."""
        counts = adata.obs[groupby].value_counts()
        proportions = counts.values / counts.sum()
        return proportions

    def _calculate_metrics(self, proportions: np.ndarray) -> Dict:
        """Calculate all diversity metrics."""
        result = {}
        counts = (proportions * 1000).astype(int)  # Approximate counts

        if "shannon" in self.metrics:
            result["shannon"] = shannon_diversity_index(proportions)
        if "simpson" in self.metrics:
            result["simpson"] = simpson_diversity_index(proportions)
        if "inverse_simpson" in self.metrics:
            result["inverse_simpson"] = inverse_simpson_index(proportions)
        if "gini_simpson" in self.metrics:
            result["gini_simpson"] = gini_simpson_index(proportions)
        if "berger_parker" in self.metrics:
            result["berger_parker"] = berger_parker_index(proportions)
        if "fisher_alpha" in self.metrics:
            result["fisher_alpha"] = fisher_alpha(counts)

        result["richness"] = np.sum(proportions > 0)
        result["evenness"] = result.get("shannon", 0) / np.log(result["richness"] + 1e-6)

        return result

    def estimate_intratumoral_heterogeneity(
        self,
        adata: AnnData,
        clone_key: str = "clone_id",
        sample_key: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Estimate intratumoral heterogeneity using multiple metrics.

        Parameters
        ----------
        adata : AnnData
            Expression data
        clone_key : str
            Column containing clone IDs
        sample_key : str, optional
            Column containing sample IDs

        Returns:
        -------
        pd.DataFrame
            Heterogeneity estimates
        """
        results = []

        samples = [None] if sample_key is None else adata.obs[sample_key].unique()

        for sample in samples:
            if sample is None:
                sample_adata = adata
            else:
                mask = adata.obs[sample_key] == sample
                sample_adata = adata[mask]

            # Clonal diversity
            clone_proportions = self._get_proportions(sample_adata, clone_key)

            # Calculate ITH metrics
            result = {
                "sample": sample or "all",
                "n_clones": len(clone_proportions),
                "clonal_diversity_shannon": shannon_diversity_index(clone_proportions),
                "clonal_diversity_simpson": simpson_diversity_index(clone_proportions),
                "dominant_clone_freq": np.max(clone_proportions),
                "purity": 1 - simpson_diversity_index(clone_proportions),
            }

            results.append(result)

        return pd.DataFrame(results)


def calculate_diversity_indices(
    adata: AnnData,
    groupby: str,
    metrics: Optional[List[str]] = None,
    sample_key: Optional[str] = None,
) -> pd.DataFrame:
    """
    Calculate diversity indices for tumor samples.

    Parameters
    ----------
    adata : AnnData
        Expression data
    groupby : str
        Column defining groups
    metrics : list, optional
        List of metrics to calculate
    sample_key : str, optional
        Column defining samples

    Returns:
    -------
    pd.DataFrame
        Diversity indices
    """
    analyzer = DiversityAnalyzer(metrics=metrics)
    return analyzer.calculate_diversity_indices(adata, groupby, sample_key)


def estimate_intratumoral_heterogeneity(
    adata: AnnData,
    clone_key: str = "clone_id",
    sample_key: Optional[str] = None,
) -> pd.DataFrame:
    """
    Estimate intratumoral heterogeneity.

    Parameters
    ----------
    adata : AnnData
        Expression data
    clone_key : str
        Column containing clone IDs
    sample_key : str, optional
        Column containing sample IDs

    Returns:
    -------
    pd.DataFrame
        Heterogeneity estimates
    """
    analyzer = DiversityAnalyzer()
    return analyzer.estimate_intratumoral_heterogeneity(adata, clone_key, sample_key)


def calculate_transcriptional_diversity(
    adata: AnnData,
    sample_key: Optional[str] = None,
    n_pcs: int = 20,
) -> pd.DataFrame:
    """
    Calculate transcriptional diversity within samples.

    Uses PCA variance and cell-to-cell distance distribution.

    Parameters
    ----------
    adata : AnnData
        Expression data
    sample_key : str, optional
        Column defining samples
    n_pcs : int
        Number of principal components

    Returns:
    -------
    pd.DataFrame
        Transcriptional diversity metrics
    """
    results = []

    # Ensure PCA exists
    if "X_pca" not in adata.obsm:
        from scanpy.preprocessing import pca

        pca(adata, n_comps=n_pcs)

    samples = [None] if sample_key is None else adata.obs[sample_key].unique()

    for sample in samples:
        if sample is None:
            sample_adata = adata
        else:
            mask = adata.obs[sample_key] == sample
            sample_adata = adata[mask]

        if sample_adata.n_obs < 10:
            continue

        # Get PCA coordinates
        X = sample_adata.obsm["X_pca"][:, :n_pcs]

        # Calculate pairwise distances
        distances = pdist(X, metric="euclidean")

        result = {
            "sample": sample or "all",
            "n_cells": sample_adata.n_obs,
            "mean_pairwise_distance": np.mean(distances),
            "std_pairwise_distance": np.std(distances),
            "median_pairwise_distance": np.median(distances),
            "max_pairwise_distance": np.max(distances),
        }

        results.append(result)

    return pd.DataFrame(results)
