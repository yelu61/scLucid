"""
Malignant cell classification.

This module provides tools for distinguishing malignant cells
from normal cells based on expression patterns and CNVs.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union
from anndata import AnnData
import logging

log = logging.getLogger(__name__)


class MalignancyClassifier:
    """
    Classify cells as malignant or normal.

    Parameters
    ----------
    method : str
        Classification method ("threshold", "cnv", "ml")
    cnv_threshold : float
        Threshold for CNV-based classification

    Attributes
    ----------
    is_malignant_ : pd.Series
        Classification results
    """

    def __init__(
        self,
        method: str = "cnv",
        cnv_threshold: float = 0.1,
    ):
        self.method = method
        self.cnv_threshold = cnv_threshold
        self.is_malignant_: Optional[pd.Series] = None

    def fit(
        self,
        adata: AnnData,
        reference_adata: Optional[AnnData] = None,
    ) -> "MalignancyClassifier":
        """
        Classify cells based on malignancy.

        Parameters
        ----------
        adata : AnnData
            Expression data
        reference_adata : AnnData, optional
            Reference normal cells

        Returns
        -------
        MalignancyClassifier
            Fitted classifier
        """
        if self.method == "cnv":
            self.is_malignant_ = self._classify_by_cnv(adata)
        elif self.method == "threshold":
            self.is_malignant_ = self._classify_by_threshold(adata, reference_adata)
        elif self.method == "ml":
            self.is_malignant_ = self._classify_by_ml(adata, reference_adata)
        else:
            raise ValueError(f"Unknown method: {self.method}")

        return self

    def _classify_by_cnv(self, adata: AnnData) -> pd.Series:
        """Classify based on CNV burden."""
        if "cnv_score" in adata.obs.columns:
            cnv_score = adata.obs["cnv_score"]
        elif "cnv" in adata.obsm:
            cnv_score = np.abs(adata.obsm["cnv"]).mean(axis=1)
        else:
            # Estimate from expression
            cnv_score = self._estimate_cnv_from_expression(adata)

        is_malignant = cnv_score > self.cnv_threshold

        return pd.Series(is_malignant, index=adata.obs_names, name="is_malignant")

    def _estimate_cnv_from_expression(self, adata: AnnData) -> np.ndarray:
        """Estimate CNV burden from expression data."""
        # Calculate expression variance across chromosome regions
        # Simplified: use overall expression variance as proxy
        expr = adata.X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray()

        # Calculate deviation from median
        median_expr = np.median(expr, axis=0)
        deviation = np.abs(expr - median_expr).mean(axis=1)

        return deviation

    def _classify_by_threshold(
        self,
        adata: AnnData,
        reference_adata: Optional[AnnData],
    ) -> pd.Series:
        """Classify based on expression threshold."""
        # Use proliferation markers as proxy
        proliferation_markers = ["MKI67", "PCNA", "CCNB1", "TOP2A"]
        available = [g for g in proliferation_markers if g in adata.var_names]

        if len(available) == 0:
            # Fall back to high expression genes
            mean_expr = np.array(adata.X.mean(axis=0)).flatten()
            top_genes = np.argsort(mean_expr)[-10:]
            scores = adata[:, top_genes].X.mean(axis=1)
        else:
            scores = adata[:, available].X.mean(axis=1)

        if hasattr(scores, 'toarray'):
            scores = scores.toarray().flatten()

        # Determine threshold
        if reference_adata is not None:
            ref_scores = reference_adata[:, available].X.mean(axis=1)
            if hasattr(ref_scores, 'toarray'):
                ref_scores = ref_scores.toarray().flatten()
            threshold = np.percentile(ref_scores, 95)
        else:
            threshold = np.percentile(scores, 75)

        is_malignant = scores > threshold

        return pd.Series(is_malignant, index=adata.obs_names, name="is_malignant")

    def _classify_by_ml(
        self,
        adata: AnnData,
        reference_adata: Optional[AnnData],
    ) -> pd.Series:
        """Classify using machine learning."""
        if reference_adata is None:
            log.warning("No reference provided, using unsupervised classification")
            return self._classify_by_cnv(adata)

        from sklearn.ensemble import RandomForestClassifier

        # Prepare training data
        X_ref = reference_adata.X
        if hasattr(X_ref, 'toarray'):
            X_ref = X_ref.toarray()

        X_tumor = adata.X
        if hasattr(X_tumor, 'toarray'):
            X_tumor = X_tumor.toarray()

        X_train = np.vstack([X_ref, X_tumor[:100]])  # Use some tumor cells
        y_train = np.array([0] * len(X_ref) + [1] * 100)

        # Train classifier
        clf = RandomForestClassifier(n_estimators=100, random_state=42)
        clf.fit(X_train, y_train)

        # Predict
        predictions = clf.predict(X_tumor)

        return pd.Series(predictions == 1, index=adata.obs_names, name="is_malignant")


def classify_malignant_cells(
    adata: AnnData,
    method: str = "cnv",
    reference_adata: Optional[AnnData] = None,
    key_added: str = "is_malignant",
) -> pd.Series:
    """
    Classify cells as malignant or normal.

    Parameters
    ----------
    adata : AnnData
        Expression data
    method : str
        Classification method
    reference_adata : AnnData, optional
        Reference normal cells
    key_added : str
        Key for storing results

    Returns
    -------
    pd.Series
        Classification results
    """
    classifier = MalignancyClassifier(method=method)
    classifier.fit(adata, reference_adata)

    adata.obs[key_added] = classifier.is_malignant_

    n_malignant = classifier.is_malignant_.sum()
    log.info(f"Classified {n_malignant} cells as malignant ({n_malignant/len(adata)*100:.1f}%)")

    return classifier.is_malignant_


def score_malignancy_potential(
    adata: AnnData,
    malignant_key: str = "is_malignant",
) -> pd.Series:
    """
    Calculate malignancy potential score for each cell.

    Parameters
    ----------
    adata : AnnData
        Expression data
    malignant_key : str
        Column with malignant classification

    Returns
    -------
    pd.Series
        Malignancy scores
    """
    # Combine multiple indicators
    scores = pd.Series(0.0, index=adata.obs_names)

    # CNV score
    if "cnv_score" in adata.obs.columns:
        scores += adata.obs["cnv_score"]

    # Proliferation score
    proliferation_markers = ["MKI67", "PCNA", "CCNB1", "TOP2A"]
    available = [g for g in proliferation_markers if g in adata.var_names]
    if len(available) > 0:
        prolif = adata[:, available].X.mean(axis=1)
        if hasattr(prolif, 'toarray'):
            prolif = prolif.toarray().flatten()
        scores += (prolif - prolif.min()) / (prolif.max() - prolif.min() + 1e-6)

    # Normalize to [0, 1]
    scores = (scores - scores.min()) / (scores.max() - scores.min() + 1e-6)

    adata.obs["malignancy_score"] = scores

    return scores
