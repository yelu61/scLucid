"""
Malignancy scoring and characterization.

This module provides methods to score cells based on their
malignant characteristics including proliferation and
metastatic potential.
"""

import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from anndata import AnnData

log = logging.getLogger(__name__)

# Cancer-type-specific marker gene sets for malignancy scoring.
# These complement the universal markers and are used when a cancer_type is specified.
CANCER_MARKER_SETS: Dict[str, Dict[str, List[str]]] = {
    "pancreatic": {
        "proliferation": ["MKI67", "PCNA", "TOP2A", "AURKA", "CCNB1", "CCNE1", "AURKB"],
        "oncogenes": ["KRAS", "MYC", "EGFR", "ERBB2", "MET", "ARNTL"],
        "tumor_suppressors": ["TP53", "CDKN2A", "SMAD4", "ARID1A", "PTEN", "BRCA2"],
    },
    "lung": {
        "proliferation": ["MKI67", "PCNA", "TOP2A", "AURKA", "CCNB1", "MCM2"],
        "oncogenes": ["EGFR", "KRAS", "BRAF", "MET", "ALK", "ROS1", "MYC"],
        "tumor_suppressors": ["TP53", "STK11", "KEAP1", "PTEN", "RB1", "CDKN2A"],
    },
    "colorectal": {
        "proliferation": ["MKI67", "PCNA", "TOP2A", "AURKA", "CCNB1", "TYMS"],
        "oncogenes": ["KRAS", "NRAS", "BRAF", "PIK3CA", "MYC", "EGFR"],
        "tumor_suppressors": ["TP53", "APC", "SMAD4", "PTEN", "FBXW7"],
    },
    "breast": {
        "proliferation": ["MKI67", "PCNA", "TOP2A", "AURKA", "CCNB1", "CCNE1"],
        "oncogenes": ["ERBB2", "MYC", "PIK3CA", "KRAS", "EGFR"],
        "tumor_suppressors": ["TP53", "BRCA1", "BRCA2", "PTEN", "RB1", "CDH1"],
    },
    "liver": {
        "proliferation": ["MKI67", "PCNA", "TOP2A", "AURKA", "CCNB1", "STMN1"],
        "oncogenes": ["MYC", "CTNNB1", "TP53", "AXIN1", "ARID1A", "TERT"],
        "tumor_suppressors": ["TP53", "RB1", "PTEN", "CDKN2A", "AXIN1"],
    },
}


def _resolve_cancer_type(cancer_type: Optional[str]) -> Optional[str]:
    """Normalize cancer type string to canonical key."""
    if cancer_type is None:
        return None
    ct = cancer_type.lower().strip()
    # Common aliases
    aliases = {
        "pdac": "pancreatic",
        "pancreatic cancer": "pancreatic",
        "pancreas": "pancreatic",
        "paad": "pancreatic",
        "luad": "lung",
        "lusc": "lung",
        "nsclc": "lung",
        "lung cancer": "lung",
        "coad": "colorectal",
        "read": "colorectal",
        "crc": "colorectal",
        "colorectal cancer": "colorectal",
        "brca": "breast",
        "breast cancer": "breast",
        "hcc": "liver",
        "lihc": "liver",
        "liver cancer": "liver",
        "hepatocellular": "liver",
    }
    return aliases.get(ct, ct)


class MalignancyScorer:
    """
    Score cells based on malignant characteristics.

    Combines multiple features including proliferation markers,
    CNV burden, and oncogene expression to calculate malignancy score.

    Parameters
    ----------
    proliferation_genes : list
        List of proliferation marker genes
    oncogene_genes : list
        List of oncogene genes
    tumor_suppressor_genes : list
        List of tumor suppressor genes
    cancer_type : str, optional
        Cancer type to load cancer-specific marker sets (e.g. "pancreatic", "lung",
        "colorectal", "breast", "liver"). Falls back to universal markers if unknown.
    """

    def __init__(
        self,
        proliferation_genes: Optional[List[str]] = None,
        oncogene_genes: Optional[List[str]] = None,
        tumor_suppressor_genes: Optional[List[str]] = None,
        cancer_type: Optional[str] = None,
    ):
        canonical = _resolve_cancer_type(cancer_type)
        cancer_markers = CANCER_MARKER_SETS.get(canonical) if canonical else None

        if cancer_markers is not None:
            log.info(f"Using cancer-specific marker set for '{canonical}' malignancy scoring.")

        self.proliferation_genes = proliferation_genes or (
            cancer_markers["proliferation"] if cancer_markers else [
                "MKI67", "PCNA", "TOP2A", "AURKA", "CCNB1",
            ]
        )
        self.oncogene_genes = oncogene_genes or (
            cancer_markers["oncogenes"] if cancer_markers else [
                "MYC", "KRAS", "EGFR", "BRAF", "PIK3CA",
            ]
        )
        self.tumor_suppressor_genes = tumor_suppressor_genes or (
            cancer_markers["tumor_suppressors"] if cancer_markers else [
                "TP53", "PTEN", "RB1", "CDKN2A", "APC",
            ]
        )
        self.scores_: Optional[pd.Series] = None

    def fit(self, adata: AnnData) -> "MalignancyScorer":
        """
        Calculate malignancy scores.

        Parameters
        ----------
        adata : AnnData
            Single-cell expression data

        Returns:
        -------
        MalignancyScorer
            Fitted scorer
        """
        scores = np.zeros(adata.n_obs)

        # Proliferation score
        prolif_score = self._calculate_gene_set_score(adata, self.proliferation_genes)
        scores += prolif_score * 0.4

        # Oncogene score
        oncogene_score = self._calculate_gene_set_score(adata, self.oncogene_genes)
        scores += oncogene_score * 0.35

        # Tumor suppressor loss (inverted)
        ts_score = self._calculate_gene_set_score(adata, self.tumor_suppressor_genes)
        scores += (1 - ts_score) * 0.25

        self.scores_ = pd.Series(scores, index=adata.obs_names)

        return self

    def _calculate_gene_set_score(
        self,
        adata: AnnData,
        gene_set: List[str],
    ) -> np.ndarray:
        """Calculate average expression score for a gene set."""
        available_genes = [g for g in gene_set if g in adata.var_names]

        if len(available_genes) == 0:
            return np.zeros(adata.n_obs)

        expr = adata[:, available_genes].X.mean(axis=1)
        if hasattr(expr, "toarray"):
            expr = expr.toarray().flatten()

        # Normalize to 0-1
        expr = (expr - expr.min()) / (expr.max() - expr.min() + 1e-10)

        return expr


def score_malignancy(
    adata: AnnData,
    proliferation_genes: Optional[List[str]] = None,
    oncogene_genes: Optional[List[str]] = None,
    key_added: str = "malignancy",
    copy: bool = False,
    cancer_type: Optional[str] = None,
) -> AnnData:
    """
    Calculate malignancy scores for all cells.

    Parameters
    ----------
    adata : AnnData
        Single-cell expression data
    proliferation_genes : list, optional
        Custom proliferation markers
    oncogene_genes : list, optional
        Custom oncogenes
    key_added : str
        Key for storing scores
    copy : bool
        Return a copy of adata
    cancer_type : str, optional
        Cancer type to load cancer-specific marker sets (e.g. "pancreatic", "lung",
        "colorectal", "breast", "liver", "pdac").

    Returns:
    -------
    AnnData
        Annotated data with malignancy scores
    """
    if copy:
        adata = adata.copy()

    scorer = MalignancyScorer(
        proliferation_genes=proliferation_genes,
        oncogene_genes=oncogene_genes,
        cancer_type=cancer_type,
    )
    scorer.fit(adata)

    adata.obs[key_added] = scorer.scores_

    log.info(f"Malignancy scoring complete. Scores in obs['{key_added}']")

    return adata


def calculate_proliferation_index(
    adata: AnnData,
    gene_set: str = "classic",
    custom_genes: Optional[List[str]] = None,
) -> pd.Series:
    """
    Calculate proliferation index based on cell cycle genes.

    Parameters
    ----------
    adata : AnnData
        Expression data
    gene_set : str
        Predefined gene set ("classic", "seurat", "tracer", "custom")
    custom_genes : list, optional
        Custom gene list if gene_set="custom"

    Returns:
    -------
    pd.Series
        Proliferation index per cell
    """
    gene_sets = {
        "classic": ["MKI67", "PCNA", "TOP2A", "AURKA", "CCNB1", "CDK1"],
        "seurat": ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4"],
        "tracer": ["MKI67", "PCNA", "TOP2A", "TYMS", "RRM2"],
    }

    if gene_set == "custom":
        genes = custom_genes
    else:
        genes = gene_sets.get(gene_set, gene_sets["classic"])

    available_genes = [g for g in genes if g in adata.var_names]

    if len(available_genes) == 0:
        raise ValueError("No proliferation genes found in data")

    expr = adata[:, available_genes].X.mean(axis=1)
    if hasattr(expr, "toarray"):
        expr = expr.toarray().flatten()

    return pd.Series(expr, index=adata.obs_names, name="proliferation_index")


def estimate_metastatic_potential(
    adata: AnnData,
    emt_genes: Optional[List[str]] = None,
    invasion_genes: Optional[List[str]] = None,
) -> pd.Series:
    """
    Estimate metastatic potential based on EMT and invasion signatures.

    Parameters
    ----------
    adata : AnnData
        Expression data
    emt_genes : list, optional
        EMT marker genes
    invasion_genes : list, optional
        Invasion-related genes

    Returns:
    -------
    pd.Series
        Metastatic potential scores
    """
    # Default EMT signature
    if emt_genes is None:
        emt_genes = ["VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "ZEB1", "TWIST1"]

    # Default invasion signature
    if invasion_genes is None:
        invasion_genes = ["MMP2", "MMP9", "MMP14", "SERPINB5", "SPP1"]

    # Calculate EMT score
    emt_available = [g for g in emt_genes if g in adata.var_names]
    if emt_available:
        emt_score = adata[:, emt_available].X.mean(axis=1)
        if hasattr(emt_score, "toarray"):
            emt_score = emt_score.toarray().flatten()
    else:
        emt_score = np.zeros(adata.n_obs)

    # Calculate invasion score
    invasion_available = [g for g in invasion_genes if g in adata.var_names]
    if invasion_available:
        invasion_score = adata[:, invasion_available].X.mean(axis=1)
        if hasattr(invasion_score, "toarray"):
            invasion_score = invasion_score.toarray().flatten()
    else:
        invasion_score = np.zeros(adata.n_obs)

    # Combine scores
    metastatic_score = 0.6 * emt_score + 0.4 * invasion_score

    return pd.Series(metastatic_score, index=adata.obs_names, name="metastatic_potential")
