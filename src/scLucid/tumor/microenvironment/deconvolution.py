"""
Tumor Microenvironment (TME) deconvolution and profiling.

This module provides tools to deconvolve and characterize
the tumor microenvironment composition.
"""

import numpy as np
import pandas as pd
from typing import Optional, Union, List, Dict
from anndata import AnnData
import logging

log = logging.getLogger(__name__)


class TMEProfiler:
    """
    Profile tumor microenvironment composition.

    This class estimates the proportions of different TME components
    including malignant cells, immune cells, and stromal cells.

    Parameters
    ----------
    cell_type_key : str
        Column in adata.obs containing cell type annotations
    """

    def __init__(self, cell_type_key: str = "cell_type"):
        self.cell_type_key = cell_type_key
        self.proportions_: Optional[pd.DataFrame] = None
        self.immune_score_: Optional[pd.Series] = None
        self.stromal_score_: Optional[pd.Series] = None

    def fit(self, adata: AnnData) -> "TMEProfiler":
        """
        Calculate TME composition.

        Parameters
        ----------
        adata : AnnData
            Single-cell data with cell type annotations

        Returns
        -------
        TMEProfiler
            Fitted profiler
        """
        cell_types = adata.obs[self.cell_type_key]

        # Calculate proportions
        self.proportions_ = cell_types.value_counts(normalize=True)

        # Calculate immune score
        immune_types = self._get_immune_types()
        self.immune_score_ = cell_types.isin(immune_types).mean()

        # Calculate stromal score
        stromal_types = self._get_stromal_types()
        self.stromal_score_ = cell_types.isin(stromal_types).mean()

        return self

    def _get_immune_types(self) -> List[str]:
        """Get list of immune cell type names."""
        return [
            "T_cell", "B_cell", "NK_cell", "Macrophage",
            "Monocyte", "Neutrophil", "DC", "Mast_cell"
        ]

    def _get_stromal_types(self) -> List[str]:
        """Get list of stromal cell type names."""
        return [
            "Fibroblast", "Endothelial", "Pericyte",
            "Stromal", "CAF"
        ]

    def get_immune_infiltration(self) -> pd.Series:
        """Get immune infiltration scores by cell type."""
        immune_types = self._get_immune_types()
        return self.proportions_[
            self.proportions_.index.isin(immune_types)
        ]

    def get_stromal_content(self) -> pd.Series:
        """Get stromal content by cell type."""
        stromal_types = self._get_stromal_types()
        return self.proportions_[
            self.proportions_.index.isin(stromal_types)
        ]


def deconvolve_tme(
    adata: AnnData,
    cell_type_key: str = "cell_type",
    key_added: str = "tme",
    copy: bool = False,
) -> AnnData:
    """
    Deconvolve tumor microenvironment composition.

    Parameters
    ----------
    adata : AnnData
        Single-cell data with cell type annotations
    cell_type_key : str
        Column with cell type labels
    key_added : str
        Key for storing results
    copy : bool
        Return a copy of adata

    Returns
    -------
    AnnData
        Annotated data with TME information
    """
    if copy:
        adata = adata.copy()

    profiler = TMEProfiler(cell_type_key=cell_type_key)
    profiler.fit(adata)

    # Store results
    adata.uns[f"{key_added}_proportions"] = profiler.proportions_
    adata.uns[f"{key_added}_immune_score"] = profiler.immune_score_
    adata.uns[f"{key_added}_stromal_score"] = profiler.stromal_score_

    log.info(f"TME profiling complete. Results stored in uns['{key_added}_*']")

    return adata


def estimate_stromal_content(
    adata: AnnData,
    cell_type_key: str = "cell_type",
    method: str = "proportion",
) -> pd.Series:
    """
    Estimate stromal content in each sample.

    Parameters
    ----------
    adata : AnnData
        Single-cell data
    cell_type_key : str
        Column with cell type labels
    method : str
        Method for estimation ("proportion", "score", "genes")

    Returns
    -------
    pd.Series
        Stromal content scores
    """
    if method == "proportion":
        stromal_types = ["Fibroblast", "Endothelial", "Pericyte", "CAF"]
        is_stromal = adata.obs[cell_type_key].isin(stromal_types)
        return pd.Series(is_stromal.astype(int), index=adata.obs_names)

    elif method == "genes":
        # Use stromal signature genes
        stromal_genes = ["ACTA2", "PDGFRB", "COL1A1", "COL3A1", "VIM"]
        available_genes = [g for g in stromal_genes if g in adata.var_names]

        if len(available_genes) == 0:
            raise ValueError("No stromal signature genes found")

        expr = adata[:, available_genes].X.mean(axis=1)
        if hasattr(expr, 'toarray'):
            expr = expr.toarray().flatten()

        return pd.Series(expr, index=adata.obs_names)

    else:
        raise ValueError(f"Unknown method: {method}")


def analyze_immune_infiltration(
    adata: AnnData,
    cell_type_key: str = "cell_type",
    groupby: Optional[str] = None,
) -> pd.DataFrame:
    """
    Analyze immune cell infiltration patterns.

    Parameters
    ----------
    adata : AnnData
        Single-cell data
    cell_type_key : str
        Column with cell type labels
    groupby : str, optional
        Column to group by (e.g., patient_id)

    Returns
    -------
    pd.DataFrame
        Immune infiltration summary
    """
    immune_types = [
        "T_cell", "CD4_T", "CD8_T", "Treg",
        "B_cell", "NK_cell", "Macrophage", "Monocyte",
        "Neutrophil", "DC", "Mast_cell"
    ]

    df = adata.obs.copy()
    is_immune = df[cell_type_key].isin(immune_types)
    df["is_immune"] = is_immune

    if groupby is not None:
        summary = df.groupby(groupby)[cell_type_key].value_counts(
            normalize=True
        ).unstack(fill_value=0)
    else:
        # Overall proportions
        summary = df[cell_type_key].value_counts(normalize=True).to_frame().T

    return summary
