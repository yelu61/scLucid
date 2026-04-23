"""
Cancer hallmark signatures and scoring utilities.

This module provides curated gene signatures for cancer hallmarks
and utilities for calculating signature scores.

Note: Gene signature data is now loaded from resources/ directory via GeneSetManager.
This module provides backward-compatible wrapper functions and scoring classes.
"""

import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import zscore

log = logging.getLogger(__name__)


class HallmarkCalculator:
    """
    Calculate cancer hallmark signature scores.

    Parameters
    ----------
    signatures : dict
        Dictionary of hallmark signatures
    method : str
        Scoring method ("mean", "ssgsea", "zscore")

    Attributes:
    ----------
    scores_ : pd.DataFrame
        Hallmark scores per cell
    """

    def __init__(
        self,
        signatures: Optional[Dict[str, List[str]]] = None,
        method: str = "mean",
    ):
        # Load from resources if not provided
        if signatures is None:
            from ...utils.manager import GeneSetManager

            gsm = GeneSetManager(species="human")
            try:
                self.signatures = gsm.load_geneset("cancer_hallmarks")
            except FileNotFoundError:
                log.warning("Could not load hallmarks from resources, using empty set")
                self.signatures = {}
        else:
            self.signatures = signatures
        self.method = method
        self.scores_: Optional[pd.DataFrame] = None

    def fit(self, adata: AnnData) -> "HallmarkCalculator":
        """
        Calculate hallmark scores.

        Parameters
        ----------
        adata : AnnData
            Expression data

        Returns:
        -------
        HallmarkCalculator
            Fitted calculator
        """
        scores = {}

        for hallmark, genes in self.signatures.items():
            # Find available genes
            available = [g for g in genes if g in adata.var_names]

            if len(available) == 0:
                log.warning(f"No genes found for {hallmark}")
                scores[hallmark] = np.zeros(adata.n_obs)
                continue

            # Calculate score
            if self.method == "mean":
                expr = adata[:, available].X.mean(axis=1)
                if hasattr(expr, "toarray"):
                    expr = expr.toarray().flatten()

            elif self.method == "zscore":
                expr = adata[:, available].X.mean(axis=1)
                if hasattr(expr, "toarray"):
                    expr = expr.toarray().flatten()
                expr = zscore(expr)

            else:
                raise ValueError(f"Unknown method: {self.method}")

            scores[hallmark] = expr

        self.scores_ = pd.DataFrame(scores, index=adata.obs_names)
        return self


def load_hallmark_signatures(
    custom_signatures: Optional[Dict[str, List[str]]] = None,
) -> Dict[str, List[str]]:
    """
    Load cancer hallmark signatures.

    Parameters
    ----------
    custom_signatures : dict, optional
        Custom signatures to add or override

    Returns:
    -------
    dict
        Hallmark signature gene sets

    Note:
    ----
    This function now loads from resources/genesets_cancer_hallmarks.json
    via GeneSetManager.
    """
    from ...utils.manager import GeneSetManager

    gsm = GeneSetManager(species="human")

    try:
        signatures = gsm.load_geneset("cancer_hallmarks")
    except FileNotFoundError:
        log.warning("Could not load hallmarks from resources, returning empty dict")
        signatures = {}

    if custom_signatures:
        signatures.update(custom_signatures)

    return signatures


def calculate_signature_scores(
    adata: AnnData,
    signatures: Optional[Dict[str, List[str]]] = None,
    method: str = "mean",
    key_added: str = "hallmark",
) -> pd.DataFrame:
    """
    Calculate signature scores for all hallmarks.

    Parameters
    ----------
    adata : AnnData
        Expression data
    signatures : dict, optional
        Custom signatures (if None, loads from resources)
    method : str
        Scoring method
    key_added : str
        Key for storing results

    Returns:
    -------
    pd.DataFrame
        Signature scores per cell
    """
    calculator = HallmarkCalculator(signatures=signatures, method=method)
    calculator.fit(adata)

    # Store in adata
    for hallmark in calculator.scores_.columns:
        adata.obs[f"{key_added}_{hallmark}"] = calculator.scores_[hallmark]

    log.info(f"Calculated {len(calculator.scores_.columns)} signature scores")

    return calculator.scores_


def get_signature_summary(
    adata: AnnData,
    signature_prefix: str = "hallmark",
    groupby: Optional[str] = None,
) -> pd.DataFrame:
    """
    Get summary statistics for hallmark scores.

    Parameters
    ----------
    adata : AnnData
        Data with hallmark scores
    signature_prefix : str
        Prefix for hallmark score columns
    groupby : str, optional
        Column to group by

    Returns:
    -------
    pd.DataFrame
        Summary statistics
    """
    hallmark_cols = [c for c in adata.obs.columns if c.startswith(signature_prefix)]

    if groupby is not None:
        return adata.obs.groupby(groupby)[hallmark_cols].mean()

    return adata.obs[hallmark_cols].describe()


# Deprecated: HALLMARK_SIGNATURES is now loaded from resources
# Use load_hallmark_signatures() instead
HALLMARK_SIGNATURES = load_hallmark_signatures()
