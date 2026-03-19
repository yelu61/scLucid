"""
Therapy response prediction and patient stratification.

This module provides tools for predicting therapy response
and stratifying patients based on molecular features.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union
from anndata import AnnData
import logging
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score

log = logging.getLogger(__name__)

# Therapy response signatures
RESPONSE_SIGNATURES = {
    "chemotherapy_sensitive": [
        "BCL2L11", "BAD", "BAK1", "BAX", "CASP8", "FAS", "TNFRSF10A",
        "CHEK1", "CHEK2", "TP53", "ATM", "ATR",
    ],
    "chemotherapy_resistant": [
        "ABCB1", "ABCC1", "ABCC2", "ABCG2", "GSTP1", "MGMT",
        "BCL2", "BCL2L1", "MCL1", "BIRC5", "XIAP",
    ],
    "immunotherapy_responsive": [
        "CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG",
        "CXCL9", "CXCL10", "CXCL11", "TBX21", "EOMES",
        "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT",
    ],
    "immunotherapy_nonresponsive": [
        "FOXP3", "IL10", "TGFB1", "VEGFA", "ARG1", "IDO1",
        "PDL1", "PDL2", "SIGLEC15", "VTCN1",
    ],
    "targeted_therapy_sensitive": [
        "EGFR", "ERBB2", "BRAF", "ALK", "ROS1", "RET",
        "KIT", "PDGFRA", "MET", "NTRK1", "NTRK2", "NTRK3",
    ],
    "targeted_therapy_resistant": [
        "NRAS", "KRAS", "PIK3CA", "PTEN", "MTOR", "AKT1",
        "EGFR_T790M", "EGFR_C797S", "BRAF_V600E",
    ],
}

# Drug response biomarkers
DRUG_BIOMARKERS = {
    "trastuzumab": {
        "predictive": {"ERBB2_amplification": "positive", "PIK3CA_mutation": "negative"},
        "genes": ["ERBB2", "PIK3CA", "PTEN", "HER3"],
    },
    "imatinib": {
        "predictive": {"BCR_ABL_fusion": "positive", "ABL_mutation": "negative"},
        "genes": ["BCR", "ABL1", "SRC", "LYN"],
    },
    "gefitinib": {
        "predictive": {"EGFR_mutation": "positive", "KRAS_mutation": "negative"},
        "genes": ["EGFR", "KRAS", "MET", "HER2"],
    },
    "vemurafenib": {
        "predictive": {"BRAF_V600E": "positive", "NRAS_mutation": "negative"},
        "genes": ["BRAF", "NRAS", "MAPK1", "MAPK3"],
    },
    "pembrolizumab": {
        "predictive": {"PDL1_expression": "high", "TMB": "high", "MSI": "positive"},
        "genes": ["CD274", "PDCD1", "CTLA4", "IFNG"],
    },
    "nivolumab": {
        "predictive": {"PDL1_expression": "high", "TMB": "high"},
        "genes": ["CD274", "PDCD1", "LAG3", "IDO1"],
    },
    "olaparib": {
        "predictive": {"BRCA1_mutation": "positive", "BRCA2_mutation": "positive", "HRD": "positive"},
        "genes": ["BRCA1", "BRCA2", "RAD51", "PARP1"],
    },
}


class ResponsePredictor:
    """
    Predict therapy response from single-cell data.

    Parameters
    ----------
    method : str
        Prediction method ("signature", "ml", "hybrid")
    signatures : dict
        Response signatures

    Attributes
    ----------
    predictions_ : pd.DataFrame
        Response predictions per cell
    """

    def __init__(
        self,
        method: str = "signature",
        signatures: Optional[Dict] = None,
    ):
        self.method = method
        self.signatures = signatures or RESPONSE_SIGNATURES
        self.predictions_: Optional[pd.DataFrame] = None
        self.model_: Optional[object] = None

    def fit(self, adata: AnnData, response_labels: Optional[pd.Series] = None) -> "ResponsePredictor":
        """
        Fit response predictor.

        Parameters
        ----------
        adata : AnnData
            Expression data
        response_labels : pd.Series, optional
            Known response labels for ML training

        Returns
        -------
        ResponsePredictor
            Fitted predictor
        """
        if self.method == "signature":
            self._fit_signature(adata)
        elif self.method == "ml":
            if response_labels is None:
                raise ValueError("ML method requires response_labels")
            self._fit_ml(adata, response_labels)
        elif self.method == "hybrid":
            self._fit_hybrid(adata, response_labels)

        return self

    def _fit_signature(self, adata: AnnData):
        """Fit using signature-based approach."""
        scores = {}

        for signature_name, genes in self.signatures.items():
            available = [g for g in genes if g in adata.var_names]

            if len(available) == 0:
                continue

            expr = adata[:, available].X.mean(axis=1)
            if hasattr(expr, 'toarray'):
                expr = expr.toarray().flatten()

            scores[signature_name] = expr

        self.predictions_ = pd.DataFrame(scores, index=adata.obs_names)

    def _fit_ml(self, adata: AnnData, response_labels: pd.Series):
        """Fit using machine learning."""
        # Prepare features
        X = adata.X
        if hasattr(X, 'toarray'):
            X = X.toarray()

        y = response_labels.loc[adata.obs_names]

        # Scale features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # Train model
        model = RandomForestClassifier(n_estimators=100, random_state=42)
        model.fit(X_scaled, y)

        self.model_ = model
        self.scaler_ = scaler

        # Get predictions
        probs = model.predict_proba(X_scaled)
        self.predictions_ = pd.DataFrame(
            probs,
            columns=[f"prob_{c}" for c in model.classes_],
            index=adata.obs_names
        )

    def _fit_hybrid(self, adata: AnnData, response_labels: Optional[pd.Series] = None):
        """Fit using hybrid approach."""
        # Start with signature scores
        self._fit_signature(adata)

        # Add ML if labels available
        if response_labels is not None:
            self._fit_ml(adata, response_labels)

    def predict_therapy_response(
        self,
        adata: AnnData,
        therapy_type: str = "chemotherapy",
    ) -> pd.DataFrame:
        """
        Predict response to a specific therapy type.

        Parameters
        ----------
        adata : AnnData
            Expression data
        therapy_type : str
            Type of therapy ("chemotherapy", "immunotherapy", "targeted")

        Returns
        -------
        pd.DataFrame
            Response predictions
        """
        if self.predictions_ is None:
            self.fit(adata)

        predictions = pd.DataFrame(index=adata.obs_names)

        if therapy_type == "chemotherapy":
            if "chemotherapy_sensitive" in self.predictions_.columns:
                sensitive = self.predictions_["chemotherapy_sensitive"]
                resistant = self.predictions_.get("chemotherapy_resistant", 0)
                predictions["response_score"] = sensitive - resistant
                predictions["predicted_response"] = predictions["response_score"] > 0

        elif therapy_type == "immunotherapy":
            if "immunotherapy_responsive" in self.predictions_.columns:
                responsive = self.predictions_["immunotherapy_responsive"]
                nonresponsive = self.predictions_.get("immunotherapy_nonresponsive", 0)
                predictions["response_score"] = responsive - nonresponsive
                predictions["predicted_response"] = predictions["response_score"] > 0

        elif therapy_type == "targeted":
            if "targeted_therapy_sensitive" in self.predictions_.columns:
                sensitive = self.predictions_["targeted_therapy_sensitive"]
                resistant = self.predictions_.get("targeted_therapy_resistant", 0)
                predictions["response_score"] = sensitive - resistant
                predictions["predicted_response"] = predictions["response_score"] > 0

        return predictions

    def stratify_patients(
        self,
        adata: AnnData,
        n_strata: int = 3,
        by: str = "resistance",
    ) -> pd.Series:
        """
        Stratify patients into response groups.

        Parameters
        ----------
        adata : AnnData
            Expression data
        n_strata : int
            Number of strata (2 or 3)
        by : str
            Stratification basis ("resistance", "sensitivity")

        Returns
        -------
        pd.Series
            Stratum assignments
        """
        if self.predictions_ is None:
            self.fit(adata)

        # Calculate composite score
        if by == "resistance":
            score_cols = [c for c in self.predictions_.columns if "resistant" in c]
        else:
            score_cols = [c for c in self.predictions_.columns if "sensitive" in c or "responsive" in c]

        if len(score_cols) == 0:
            score_cols = self.predictions_.columns

        composite = self.predictions_[score_cols].mean(axis=1)

        # Stratify
        if n_strata == 2:
            labels = ["low", "high"]
        elif n_strata == 3:
            labels = ["low", "intermediate", "high"]
        else:
            labels = [f"stratum_{i}" for i in range(n_strata)]

        strata = pd.qcut(composite, q=n_strata, labels=labels)

        return strata


def predict_therapy_response(
    adata: AnnData,
    therapy_type: str = "chemotherapy",
    method: str = "signature",
    key_added: str = "therapy_response",
) -> pd.DataFrame:
    """
    Predict therapy response for all cells.

    Parameters
    ----------
    adata : AnnData
        Expression data
    therapy_type : str
        Type of therapy
    method : str
        Prediction method
    key_added : str
        Key prefix for storing results

    Returns
    -------
    pd.DataFrame
        Response predictions
    """
    predictor = ResponsePredictor(method=method)
    predictions = predictor.predict_therapy_response(adata, therapy_type)

    # Store in adata
    for col in predictions.columns:
        adata.obs[f"{key_added}_{col}"] = predictions[col]

    log.info(f"Predicted {therapy_type} response for {len(predictions)} cells")

    return predictions


def stratify_patients(
    adata: AnnData,
    n_strata: int = 3,
    by: str = "resistance",
    key_added: str = "stratum",
) -> pd.Series:
    """
    Stratify patients into response groups.

    Parameters
    ----------
    adata : AnnData
        Expression data
    n_strata : int
        Number of strata
    by : str
        Stratification basis
    key_added : str
        Key for storing results

    Returns
    -------
    pd.Series
        Stratum assignments
    """
    predictor = ResponsePredictor()
    predictor.fit(adata)

    strata = predictor.stratify_patients(adata, n_strata=n_strata, by=by)

    adata.obs[key_added] = strata

    return strata


def evaluate_biomarker(
    adata: AnnData,
    drug: str,
    expression_threshold: float = 0.5,
) -> pd.DataFrame:
    """
    Evaluate biomarker expression for a specific drug.

    Parameters
    ----------
    adata : AnnData
        Expression data
    drug : str
        Drug name
    expression_threshold : float
        Threshold for positive expression

    Returns
    -------
    pd.DataFrame
        Biomarker evaluation results
    """
    if drug not in DRUG_BIOMARKERS:
        log.warning(f"No biomarkers defined for {drug}")
        return pd.DataFrame()

    biomarkers = DRUG_BIOMARKERS[drug]
    genes = biomarkers["genes"]

    results = pd.DataFrame(index=adata.obs_names)

    for gene in genes:
        if gene in adata.var_names:
            expr = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, 'toarray') else adata[:, gene].X
            results[f"{gene}_expression"] = expr
            results[f"{gene}_positive"] = expr > expression_threshold

    # Calculate composite biomarker score
    available_genes = [g for g in genes if g in adata.var_names]
    if len(available_genes) > 0:
        expr_matrix = adata[:, available_genes].X
        if hasattr(expr_matrix, 'toarray'):
            expr_matrix = expr_matrix.toarray()
        results["biomarker_score"] = expr_matrix.mean(axis=1)

    return results
