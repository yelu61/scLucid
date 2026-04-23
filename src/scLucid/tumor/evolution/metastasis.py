"""
Metastasis tracking and dissemination analysis.

This module provides tools for predicting metastasis risk
and analyzing tumor cell dissemination patterns.
"""

import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from anndata import AnnData

log = logging.getLogger(__name__)

# Metastasis-associated gene signatures
METASTASIS_SIGNATURES = {
    "emt": [
        "CDH2",
        "VIM",
        "SNAI1",
        "SNAI2",
        "TWIST1",
        "ZEB1",
        "ZEB2",
        "FN1",
        "MMP2",
        "MMP3",
        "MMP9",
        "TGFB1",
        "TGFB2",
    ],
    "invasion": [
        "MMP1",
        "MMP2",
        "MMP3",
        "MMP7",
        "MMP9",
        "MMP13",
        "PLAU",
        "PLAUR",
        "CTSD",
        "CTSB",
        "HEBP2",
    ],
    "angiogenesis": [
        "VEGFA",
        "VEGFB",
        "VEGFC",
        "ANGPT1",
        "ANGPT2",
        "PDGFA",
        "FGF2",
        "HGF",
        "TGFA",
        "IGF1",
        "IGF2",
    ],
    "stemness": [
        "PROM1",
        "CD44",
        "ALDH1A1",
        "ALDH1A3",
        "NANOG",
        "SOX2",
        "POU5F1",
        "KLF4",
        "MYC",
        "BMI1",
        "EZH2",
    ],
    "proliferation": [
        "MKI67",
        "PCNA",
        "CCNB1",
        "CCND1",
        "CDK1",
        "CDK4",
        "CDK6",
        "E2F1",
        "TOP2A",
        "AURKA",
        "AURKB",
    ],
}

# Organ-specific metastasis signatures
ORGAN_TROPISM = {
    "bone": ["ITGAV", "ITGB3", "CXCR4", "CTGF", "PTHRP", "MMP1"],
    "lung": ["CXCR4", "CXCL12", "ANGPTL4", "S100A4", "S100A9", "COX2"],
    "liver": ["CXCR4", "CXCL12", "IGF1", "IGF1R", "PIGF", "CEACAM5"],
    "brain": ["COX2", "HBEGF", "ANGPTL4", "PTGS2", "S100A4", "L1CAM"],
    "lymph_node": ["CCR7", "CXCR4", "CXCL13", "CCL21", "LTB", "CD40"],
}


class MetastasisTracker:
    """
    Track and analyze metastasis potential.

    Parameters
    ----------
    signatures : dict
        Metastasis-associated gene signatures

    Attributes:
    ----------
    risk_scores_ : pd.DataFrame
        Metastasis risk scores per cell
    """

    def __init__(
        self,
        signatures: Optional[Dict[str, List[str]]] = None,
    ):
        self.signatures = signatures or METASTASIS_SIGNATURES
        self.risk_scores_: Optional[pd.DataFrame] = None

    def predict_metastasis_risk(
        self,
        adata: AnnData,
        organ: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Predict metastasis risk for each cell.

        Parameters
        ----------
        adata : AnnData
            Expression data
        organ : str, optional
            Target organ for tropism prediction

        Returns:
        -------
        pd.DataFrame
            Risk scores per cell
        """
        scores = {}

        # Calculate signature scores
        for sig_name, genes in self.signatures.items():
            available = [g for g in genes if g in adata.var_names]

            if len(available) == 0:
                continue

            expr = adata[:, available].X.mean(axis=1)
            if hasattr(expr, "toarray"):
                expr = expr.toarray().flatten()

            scores[sig_name] = expr

        # Organ-specific tropism
        if organ and organ in ORGAN_TROPISM:
            tropism_genes = ORGAN_TROPISM[organ]
            available = [g for g in tropism_genes if g in adata.var_names]

            if len(available) > 0:
                expr = adata[:, available].X.mean(axis=1)
                if hasattr(expr, "toarray"):
                    expr = expr.toarray().flatten()
                scores[f"{organ}_tropism"] = expr

        self.risk_scores_ = pd.DataFrame(scores, index=adata.obs_names)

        # Calculate composite risk score
        self.risk_scores_["composite_risk"] = self.risk_scores_.mean(axis=1)

        return self.risk_scores_

    def analyze_dissemination(
        self,
        adata: AnnData,
        primary_key: str = "is_primary",
        metastatic_key: str = "is_metastatic",
    ) -> pd.DataFrame:
        """
        Analyze dissemination from primary to metastatic sites.

        Parameters
        ----------
        adata : AnnData
            Expression data
        primary_key : str
            Column indicating primary tumor cells
        metastatic_key : str
            Column indicating metastatic cells

        Returns:
        -------
        pd.DataFrame
            Dissemination analysis results
        """
        results = []

        # Check if keys exist
        has_primary = primary_key in adata.obs.columns
        has_metastatic = metastatic_key in adata.obs.columns

        if not (has_primary or has_metastatic):
            log.warning("Neither primary nor metastatic labels found")
            return pd.DataFrame()

        # Compare signature scores
        for sig_name in self.signatures.keys():
            if sig_name not in self.risk_scores_.columns:
                continue

            scores = self.risk_scores_[sig_name]

            if has_primary and has_metastatic:
                primary_mask = adata.obs[primary_key]
                meta_mask = adata.obs[metastatic_key]

                primary_mean = scores[primary_mask].mean()
                meta_mean = scores[meta_mask].mean()

                result = {
                    "signature": sig_name,
                    "primary_mean": primary_mean,
                    "metastatic_mean": meta_mean,
                    "fold_change": meta_mean / (primary_mean + 1e-6),
                }

                # Statistical test
                from scipy import stats

                _, pval = stats.mannwhitneyu(
                    scores[primary_mask], scores[meta_mask], alternative="two-sided"
                )
                result["pvalue"] = pval

                results.append(result)

        return pd.DataFrame(results)

    def identify_seeding_clones(
        self,
        adata: AnnData,
        clone_key: str = "clone_id",
        site_key: str = "site",
        primary_value: str = "primary",
    ) -> pd.DataFrame:
        """
        Identify clones that likely seeded metastases.

        Parameters
        ----------
        adata : AnnData
            Expression data
        clone_key : str
            Column containing clone IDs
        site_key : str
            Column indicating sample site
        primary_value : str
            Value indicating primary tumor

        Returns:
        -------
        pd.DataFrame
            Potential seeding clones
        """
        results = []

        if site_key not in adata.obs.columns:
            raise ValueError(f"Site column {site_key} not found")

        sites = adata.obs[site_key].unique()
        metastatic_sites = [s for s in sites if s != primary_value]

        for clone in adata.obs[clone_key].unique():
            clone_mask = adata.obs[clone_key] == clone
            clone_adata = adata[clone_mask]

            # Check if clone appears in both primary and metastasis
            clone_sites = clone_adata.obs[site_key].unique()

            in_primary = primary_value in clone_sites
            in_metastasis = any(s in clone_sites for s in metastatic_sites)

            if in_primary and in_metastasis:
                # Calculate frequency in each site
                site_freqs = clone_adata.obs[site_key].value_counts(normalize=True)

                result = {
                    "clone": clone,
                    "n_cells": clone_adata.n_obs,
                    "primary_freq": site_freqs.get(primary_value, 0),
                    "seeding_potential": site_freqs.get(primary_value, 0)
                    * len([s for s in metastatic_sites if s in site_freqs]),
                }

                for site in metastatic_sites:
                    result[f"{site}_freq"] = site_freqs.get(site, 0)

                results.append(result)

        results_df = pd.DataFrame(results)
        if len(results_df) > 0:
            results_df = results_df.sort_values("seeding_potential", ascending=False)

        return results_df


def predict_metastasis_risk(
    adata: AnnData,
    organ: Optional[str] = None,
    key_added: str = "metastasis_risk",
) -> pd.DataFrame:
    """
    Predict metastasis risk for tumor cells.

    Parameters
    ----------
    adata : AnnData
        Expression data
    organ : str, optional
        Target organ for tropism
    key_added : str
        Key prefix for storing results

    Returns:
    -------
    pd.DataFrame
        Risk scores
    """
    tracker = MetastasisTracker()
    scores = tracker.predict_metastasis_risk(adata, organ)

    # Store in adata
    for col in scores.columns:
        adata.obs[f"{key_added}_{col}"] = scores[col]

    log.info(f"Calculated metastasis risk for {len(scores)} cells")

    return scores


def analyze_dissemination(
    adata: AnnData,
    primary_key: str = "is_primary",
    metastatic_key: str = "is_metastatic",
) -> pd.DataFrame:
    """
    Analyze dissemination from primary to metastasis.

    Parameters
    ----------
    adata : AnnData
        Expression data
    primary_key : str
        Primary tumor indicator
    metastatic_key : str
        Metastatic indicator

    Returns:
    -------
    pd.DataFrame
        Dissemination analysis
    """
    tracker = MetastasisTracker()
    tracker.predict_metastasis_risk(adata)

    return tracker.analyze_dissemination(adata, primary_key, metastatic_key)


def compare_primary_vs_metastasis(
    adata: AnnData,
    site_key: str = "site",
    primary_value: str = "primary",
    method: str = "wilcoxon",
) -> pd.DataFrame:
    """
    Compare gene expression between primary and metastatic tumors.

    Parameters
    ----------
    adata : AnnData
        Expression data
    site_key : str
        Column indicating site
    primary_value : str
        Value for primary tumor
    method : str
        Statistical test method

    Returns:
    -------
    pd.DataFrame
        Differential expression results
    """
    from scipy import stats

    sites = adata.obs[site_key].unique()
    metastatic_sites = [s for s in sites if s != primary_value]

    if len(metastatic_sites) == 0:
        raise ValueError("No metastatic sites found")

    primary_mask = adata.obs[site_key] == primary_value
    meta_mask = adata.obs[site_key].isin(metastatic_sites)

    results = []

    for gene in adata.var_names:
        expr_primary = adata[primary_mask, gene].X
        expr_meta = adata[meta_mask, gene].X

        if hasattr(expr_primary, "toarray"):
            expr_primary = expr_primary.toarray().flatten()
        if hasattr(expr_meta, "toarray"):
            expr_meta = expr_meta.toarray().flatten()

        if method == "wilcoxon":
            try:
                stat, pval = stats.ranksums(expr_primary, expr_meta)
            except:
                continue
        elif method == "t-test":
            stat, pval = stats.ttest_ind(expr_primary, expr_meta)
        else:
            raise ValueError(f"Unknown method: {method}")

        results.append(
            {
                "gene": gene,
                "primary_mean": np.mean(expr_primary),
                "metastatic_mean": np.mean(expr_meta),
                "log2fc": np.log2((np.mean(expr_meta) + 1e-6) / (np.mean(expr_primary) + 1e-6)),
                "statistic": stat,
                "pvalue": pval,
            }
        )

    return pd.DataFrame(results).sort_values("pvalue")
