"""
Cancer stemness analysis.

This module provides tools for quantifying cancer stemness
and identifying cancer stem cell populations.
"""

import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from anndata import AnnData

log = logging.getLogger(__name__)


# Cancer stem cell signatures
STEMNESS_SIGNATURES = {
    "core_stemness": [
        "PROM1",
        "CD44",
        "ALDH1A1",
        "ALDH1A3",
        "NANOG",
        "SOX2",
        "POU5F1",
        "MYC",
        "KLF4",
        "LIN28A",
        "SALL4",
        "DPPA4",
    ],
    "embryonic_stemness": [
        "NANOG",
        "SOX2",
        "POU5F1",
        "LIN28A",
        "SALL4",
        "DPPA2",
        "DPPA4",
        "ZFP42",
        "TDGF1",
        "GDF3",
        "UTF1",
        "LEFTY1",
    ],
    "tissue_stemness": [
        "LGR5",
        "ASCL2",
        "SMOC2",
        "OLFM4",
        "AXIN2",
        "EPHB2",
        "TNFRSF19",
        "SOX9",
        "CD44",
        "ITGB1",
        "EPCAM",
        "PROM1",
    ],
    "epithelial_plasticity": [
        "CDH1",
        "CDH2",
        "VIM",
        "CLDN3",
        "CLDN4",
        "CLDN7",
        "EPCAM",
        "MUC1",
        "CD24",
        "CD44",
        "ALCAM",
        "MET",
    ],
    "quiescence": [
        "CDKN1B",
        "CDKN1C",
        "TGFB1",
        "BTG1",
        "BTG2",
        "ID1",
        "ID2",
        "ID3",
        "CCND1",
        "MYC",
        "MTOR",
        "AMPK",
    ],
}


class StemnessAnalyzer:
    """
    Analyze cancer stemness in tumor cells.

    Parameters
    ----------
    signatures : dict
        Stemness gene signatures
    method : str
        Scoring method

    Attributes:
    ----------
    stemness_scores_ : pd.DataFrame
        Stemness scores per cell
    """

    def __init__(
        self,
        signatures: Optional[Dict[str, List[str]]] = None,
        method: str = "mean",
    ):
        self.signatures = signatures or STEMNESS_SIGNATURES
        self.method = method
        self.stemness_scores_: Optional[pd.DataFrame] = None

    def calculate_stemness_score(self, adata: AnnData) -> pd.DataFrame:
        """
        Calculate stemness scores for each cell.

        Parameters
        ----------
        adata : AnnData
            Expression data

        Returns:
        -------
        pd.DataFrame
            Stemness scores
        """
        scores = {}

        for sig_name, genes in self.signatures.items():
            available = [g for g in genes if g in adata.var_names]

            if len(available) == 0:
                log.warning(f"No genes found for {sig_name}")
                continue

            expr = adata[:, available].X.mean(axis=1)
            if hasattr(expr, "toarray"):
                expr = expr.toarray().flatten()

            scores[sig_name] = expr

        self.stemness_scores_ = pd.DataFrame(scores, index=adata.obs_names)

        # Calculate composite score
        self.stemness_scores_["composite_stemness"] = self.stemness_scores_.mean(axis=1)

        return self.stemness_scores_

    def identify_cancer_stem_cells(
        self,
        adata: AnnData,
        threshold: Optional[float] = None,
        top_n_percent: float = 10.0,
    ) -> pd.Series:
        """
        Identify cancer stem cell population.

        Parameters
        ----------
        adata : AnnData
            Expression data
        threshold : float, optional
            Stemness threshold
        top_n_percent : float
            Percentage of top cells to consider as CSCs

        Returns:
        -------
        pd.Series
            Boolean indicating CSCs
        """
        if self.stemness_scores_ is None:
            self.calculate_stemness_score(adata)

        composite = self.stemness_scores_["composite_stemness"]

        if threshold is None:
            threshold = np.percentile(composite, 100 - top_n_percent)

        is_csc = composite > threshold

        return pd.Series(is_csc, index=adata.obs_names, name="is_cancer_stem_cell")

    def analyze_stemness_dynamics(
        self,
        adata: AnnData,
        pseudotime_key: str = "pseudotime",
    ) -> pd.DataFrame:
        """
        Analyze stemness dynamics along pseudotime.

        Parameters
        ----------
        adata : AnnData
            Expression data with pseudotime
        pseudotime_key : str
            Column containing pseudotime

        Returns:
        -------
        pd.DataFrame
            Stemness dynamics
        """
        if pseudotime_key not in adata.obs.columns:
            raise ValueError(f"Pseudotime column {pseudotime_key} not found")

        if self.stemness_scores_ is None:
            self.calculate_stemness_score(adata)

        pseudotime = adata.obs[pseudotime_key]

        # Bin pseudotime
        bins = np.linspace(0, 1, 11)
        bin_centers = (bins[:-1] + bins[1:]) / 2

        results = []
        for i, (start, end) in enumerate(zip(bins[:-1], bins[1:])):
            mask = (pseudotime >= start) & (pseudotime < end)

            if mask.sum() < 5:
                continue

            result = {
                "pseudotime_bin": bin_centers[i],
                "n_cells": mask.sum(),
            }

            for sig in self.stemness_scores_.columns:
                result[sig] = self.stemness_scores_.loc[mask, sig].mean()

            results.append(result)

        return pd.DataFrame(results)


def calculate_stemness_score(
    adata: AnnData,
    method: str = "mean",
    key_added: str = "stemness",
) -> pd.DataFrame:
    """
    Calculate stemness scores for tumor cells.

    Parameters
    ----------
    adata : AnnData
        Expression data
    method : str
        Scoring method
    key_added : str
        Key prefix for storing results

    Returns:
    -------
    pd.DataFrame
        Stemness scores
    """
    analyzer = StemnessAnalyzer(method=method)
    scores = analyzer.calculate_stemness_score(adata)

    # Store in adata
    for col in scores.columns:
        adata.obs[f"{key_added}_{col}"] = scores[col]

    log.info(f"Calculated stemness scores for {len(scores)} cells")

    return scores


def identify_cancer_stem_cells(
    adata: AnnData,
    threshold: Optional[float] = None,
    top_n_percent: float = 10.0,
    key_added: str = "is_cancer_stem_cell",
) -> pd.Series:
    """
    Identify cancer stem cell population.

    Parameters
    ----------
    adata : AnnData
        Expression data
    threshold : float, optional
        Stemness threshold
    top_n_percent : float
        Percentage of top cells to consider as CSCs
    key_added : str
        Key for storing results

    Returns:
    -------
    pd.Series
        Boolean indicating CSCs
    """
    analyzer = StemnessAnalyzer()
    is_csc = analyzer.identify_cancer_stem_cells(adata, threshold, top_n_percent)

    adata.obs[key_added] = is_csc

    log.info(f"Identified {is_csc.sum()} cancer stem cells ({is_csc.sum()/len(adata)*100:.1f}%)")

    return is_csc


def compare_stemness_between_groups(
    adata: AnnData,
    groupby: str,
    method: str = "mann-whitney",
) -> pd.DataFrame:
    """
    Compare stemness between cell groups.

    Parameters
    ----------
    adata : AnnData
        Expression data
    groupby : str
        Column to group cells
    method : str
        Statistical test method

    Returns:
    -------
    pd.DataFrame
        Comparison results
    """
    from scipy import stats

    analyzer = StemnessAnalyzer()
    scores = analyzer.calculate_stemness_score(adata)

    groups = adata.obs[groupby].unique()

    results = []

    for sig in scores.columns:
        group_scores = []
        for group in groups:
            mask = adata.obs[groupby] == group
            group_scores.append(scores.loc[mask, sig].values)

        if len(groups) == 2:
            if method == "mann-whitney":
                stat, pval = stats.mannwhitneyu(group_scores[0], group_scores[1])
            elif method == "t-test":
                stat, pval = stats.ttest_ind(group_scores[0], group_scores[1])
            else:
                raise ValueError(f"Unknown method: {method}")

            results.append(
                {
                    "signature": sig,
                    "group1": groups[0],
                    "group2": groups[1],
                    "mean1": np.mean(group_scores[0]),
                    "mean2": np.mean(group_scores[1]),
                    "statistic": stat,
                    "pvalue": pval,
                }
            )

    return pd.DataFrame(results)
