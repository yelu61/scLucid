"""
CNV signature analysis and extraction.

This module provides tools for extracting and analyzing
copy number variation signatures similar to mutational signatures.
"""

import logging
from typing import Optional

import pandas as pd
from anndata import AnnData

log = logging.getLogger(__name__)


# Known CNV signatures from literature
REFERENCE_CNV_SIGNATURES = {
    "CX1": {
        "description": "Extensive duplication",
        "pattern": [1, 1, 2, 2, 2, 1, 1, 1, 1, 1],
    },
    "CX2": {
        "description": "Chromosome 1q gain, 16q loss",
        "pattern": [1, 1, 1, 1, 1, 1, 1, 0, 1, 0],
    },
    "CX3": {
        "description": "Extensive LOH",
        "pattern": [1, 0, 1, 1, 0, 0, 1, 0, 1, 0],
    },
    "CX4": {
        "description": "Chromosome 8 gain",
        "pattern": [1, 1, 1, 1, 1, 1, 1, 2, 1, 1],
    },
    "CX5": {
        "description": "Amplification",
        "pattern": [1, 1, 1, 3, 1, 1, 1, 1, 1, 1],
    },
    "CX6": {
        "description": "Chromosome 20 gain",
        "pattern": [1, 1, 1, 1, 1, 1, 1, 1, 1, 2],
    },
    "CX7": {
        "description": "Flat",
        "pattern": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    },
    "CX8": {
        "description": "Chromosome 17 amplification",
        "pattern": [1, 1, 1, 1, 1, 1, 2, 1, 1, 1],
    },
}


class CNVSigExtractor:
    """
    Extract CNV signatures from single-cell data.

    Parameters
    ----------
    n_components : int
        Number of signatures to extract
    method : str
        Decomposition method ("nmf", "pca")

    Attributes:
    ----------
    signatures_ : pd.DataFrame
        Extracted CNV signatures
    exposures_ : pd.DataFrame
        Signature exposure per cell
    """

    def __init__(
        self,
        n_components: int = 5,
        method: str = "nmf",
    ):
        self.n_components = n_components
        self.method = method
        self.signatures_: Optional[pd.DataFrame] = None
        self.exposures_: Optional[pd.DataFrame] = None

    def fit(self, adata: AnnData, cnv_key: str = "cnv") -> "CNVSigExtractor":
        """
        Extract CNV signatures.

        Parameters
        ----------
        adata : AnnData
            Expression data with CNV
        cnv_key : str
            Key for CNV matrix in adata.obsm

        Returns:
        -------
        CNVSigExtractor
            Fitted extractor
        """
        if cnv_key not in adata.obsm:
            raise ValueError(f"CNV data not found in obsm['{cnv_key}']")

        cnv_matrix = adata.obsm[cnv_key]

        if self.method == "nmf":
            from sklearn.decomposition import NMF

            model = NMF(n_components=self.n_components, random_state=42, max_iter=500)
            self.exposures_ = pd.DataFrame(
                model.fit_transform(cnv_matrix),
                index=adata.obs_names,
                columns=[f"CNV_Sig{i+1}" for i in range(self.n_components)],
            )
            self.signatures_ = pd.DataFrame(
                model.components_,
                index=[f"CNV_Sig{i+1}" for i in range(self.n_components)],
                columns=range(cnv_matrix.shape[1]) if hasattr(cnv_matrix, "shape") else None,
            )

        elif self.method == "pca":
            from sklearn.decomposition import PCA

            model = PCA(n_components=self.n_components)
            self.exposures_ = pd.DataFrame(
                model.fit_transform(cnv_matrix),
                index=adata.obs_names,
                columns=[f"PC{i+1}" for i in range(self.n_components)],
            )
            self.signatures_ = pd.DataFrame(
                model.components_, index=[f"PC{i+1}" for i in range(self.n_components)]
            )

        else:
            raise ValueError(f"Unknown method: {self.method}")

        return self

    def assign_signature(self, adata: AnnData) -> pd.Series:
        """
        Assign dominant signature to each cell.

        Parameters
        ----------
        adata : AnnData
            Expression data

        Returns:
        -------
        pd.Series
            Dominant signature per cell
        """
        if self.exposures_ is None:
            raise ValueError("CNVSigExtractor not fitted yet")

        dominant = self.exposures_.idxmax(axis=1)

        adata.obs["cnv_signature"] = dominant

        return dominant

    def compare_to_reference(self) -> pd.DataFrame:
        """
        Compare extracted signatures to reference signatures.

        Returns:
        -------
        pd.DataFrame
            Comparison results
        """
        if self.signatures_ is None:
            raise ValueError("CNVSigExtractor not fitted yet")

        results = []

        # Simplified comparison using correlation
        for sig_name, sig_data in REFERENCE_CNV_SIGNATURES.items():
            # This is a placeholder for actual signature comparison
            results.append(
                {
                    "reference": sig_name,
                    "description": sig_data["description"],
                    "best_match": "unknown",
                    "correlation": 0.0,
                }
            )

        return pd.DataFrame(results)


def extract_cnv_signatures(
    adata: AnnData,
    n_components: int = 5,
    cnv_key: str = "cnv",
    key_added: str = "cnv_sig",
) -> pd.DataFrame:
    """
    Extract CNV signatures from data.

    Parameters
    ----------
    adata : AnnData
        Expression data with CNV
    n_components : int
        Number of signatures
    cnv_key : str
        Key for CNV data
    key_added : str
        Key prefix for storing results

    Returns:
    -------
    pd.DataFrame
        Signature exposures
    """
    extractor = CNVSigExtractor(n_components=n_components)
    extractor.fit(adata, cnv_key)

    # Store exposures in adata
    for col in extractor.exposures_.columns:
        adata.obs[f"{key_added}_{col}"] = extractor.exposures_[col]

    log.info(f"Extracted {n_components} CNV signatures")

    return extractor.exposures_


def assign_cnv_signature(
    adata: AnnData,
    cnv_key: str = "cnv",
    key_added: str = "cnv_signature",
) -> pd.Series:
    """
    Assign CNV signature to each cell.

    Parameters
    ----------
    adata : AnnData
        Expression data
    cnv_key : str
        Key for CNV data
    key_added : str
        Key for storing results

    Returns:
    -------
    pd.Series
        Signature assignments
    """
    extractor = CNVSigExtractor()
    extractor.fit(adata, cnv_key)
    dominant = extractor.assign_signature(adata)

    adata.obs[key_added] = dominant

    return dominant
