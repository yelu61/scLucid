"""
Tumor Microenvironment (TME) deconvolution and profiling.

This module provides tools to deconvolve and characterize
the tumor microenvironment composition.
"""

import logging
from typing import Dict, List, Optional

import pandas as pd
from anndata import AnnData

log = logging.getLogger(__name__)


def _normalize_compartment_label(label: str) -> str:
    """Normalize user-provided cell labels for compartment mapping."""
    return str(label).lower().replace("-", " ").replace("_", " ").strip()


class TMEProfiler:
    """
    Profile tumor microenvironment composition.

    This class estimates the proportions of different TME components
    including malignant cells, immune cells, and stromal cells.

    Parameters
    ----------
    cell_type_key : str
        Column in adata.obs containing cell type annotations
    compartment_map : dict, optional
        Mapping from input cell type labels to canonical compartments
        ("immune", "stromal", "malignant", "other"). If None, uses an
        internal normalization map that covers common label variants.
    """

    # Canonical compartment map. Keys are lower-cased input labels.
    _DEFAULT_COMPARTMENT_MAP: Dict[str, str] = {
        # T cells
        "t_cell": "immune",
        "t cell": "immune",
        "t-cell": "immune",
        "tcell": "immune",
        "cd4_t": "immune",
        "cd4 t": "immune",
        "cd4": "immune",
        "cd8_t": "immune",
        "cd8 t": "immune",
        "cd8": "immune",
        "treg": "immune",
        "regulatory t": "immune",
        # B cells
        "b_cell": "immune",
        "b cell": "immune",
        "b-cell": "immune",
        "bcell": "immune",
        "plasma": "immune",
        "plasma_cell": "immune",
        "plasma cell": "immune",
        # NK / myeloid
        "nk": "immune",
        "nk_cell": "immune",
        "nk cell": "immune",
        "macrophage": "immune",
        "monocyte": "immune",
        "neutrophil": "immune",
        "dc": "immune",
        "dendritic": "immune",
        "mast_cell": "immune",
        "mast cell": "immune",
        "myeloid": "immune",
        # Stromal
        "fibroblast": "stromal",
        "endothelial": "stromal",
        "pericyte": "stromal",
        "stromal": "stromal",
        "caf": "stromal",
        "myofibroblast": "stromal",
        # Malignant / tumor labels
        "tumor": "malignant",
        "malignant": "malignant",
        "cancer": "malignant",
        "carcinoma": "malignant",
        "adenocarcinoma": "malignant",
        "tumor epithelial": "malignant",
        "malignant epithelial": "malignant",
        # Common non-malignant epithelial labels (treated as reference)
        "epithelial": "other",
        "normal_epithelial": "other",
        "normal epithelial": "other",
        "normal": "other",
    }

    def __init__(
        self,
        cell_type_key: str = "cell_type",
        compartment_map: Optional[Dict[str, str]] = None,
    ):
        self.cell_type_key = cell_type_key
        raw_map = compartment_map or self._DEFAULT_COMPARTMENT_MAP
        self.compartment_map = {
            _normalize_compartment_label(k): v for k, v in raw_map.items()
        }
        self.proportions_: Optional[pd.Series] = None
        self.immune_score_: Optional[float] = None
        self.stromal_score_: Optional[float] = None
        self.malignant_score_: Optional[float] = None
        self.normalized_labels_: Optional[pd.Series] = None
        self.unmapped_labels_: List[str] = []

    def fit(self, adata: AnnData) -> "TMEProfiler":
        """
        Calculate TME composition.

        Parameters
        ----------
        adata : AnnData
            Single-cell data with cell type annotations

        Returns:
        -------
        TMEProfiler
            Fitted profiler
        """
        cell_types = adata.obs[self.cell_type_key].astype(str)

        # Normalize labels using the compartment map
        normalized = cell_types.map(_normalize_compartment_label)
        compartments = normalized.map(self.compartment_map)
        self.unmapped_labels_ = sorted(
            normalized[compartments.isna()].unique().tolist()
        )
        if self.unmapped_labels_:
            log.debug(
                f"TMEProfiler: unmapped labels treated as malignant/other: "
                f"{self.unmapped_labels_}"
            )
        # Unmapped labels default to "other" so they do not inflate malignant
        # estimates; they remain listed for manual review.
        compartments = compartments.fillna("other")
        self.normalized_labels_ = compartments

        # Calculate proportions
        self.proportions_ = compartments.value_counts(normalize=True)

        # Calculate compartment scores
        self.immune_score_ = float((compartments == "immune").mean())
        self.stromal_score_ = float((compartments == "stromal").mean())
        self.malignant_score_ = float((compartments == "malignant").mean())

        return self

    def _get_immune_types(self) -> List[str]:
        """Get list of canonical immune compartment labels."""
        return ["immune"]

    def _get_stromal_types(self) -> List[str]:
        """Get list of canonical stromal compartment labels."""
        return ["stromal"]

    def get_immune_infiltration(self) -> pd.Series:
        """Get immune infiltration proportion."""
        if self.proportions_ is None:
            raise ValueError("Profiler has not been fitted.")
        value = float(self.proportions_.get("immune", 0.0))
        return pd.Series([value], index=["immune"])

    def get_stromal_content(self) -> pd.Series:
        """Get stromal content proportion."""
        if self.proportions_ is None:
            raise ValueError("Profiler has not been fitted.")
        value = float(self.proportions_.get("stromal", 0.0))
        return pd.Series([value], index=["stromal"])


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

    Returns:
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
    adata.uns[f"{key_added}_malignant_score"] = profiler.malignant_score_
    adata.uns[f"{key_added}_compartment_claim"] = (
        "annotation-derived TME composition; not a bulk deconvolution"
    )

    # Per-cell TME scores (useful for downstream analysis and visualization)
    compartments = profiler.normalized_labels_
    if compartments is not None:
        adata.obs[f"{key_added}_compartment"] = compartments.astype(str)
        adata.obs[f"{key_added}_is_immune"] = (compartments == "immune").astype(int)
        adata.obs[f"{key_added}_is_stromal"] = (compartments == "stromal").astype(int)
        adata.obs[f"{key_added}_is_malignant"] = (compartments == "malignant").astype(int)

    if profiler.unmapped_labels_:
        adata.uns[f"{key_added}_unmapped_labels"] = profiler.unmapped_labels_

    log.info(
        f"TME profiling complete. Results in uns['{key_added}_*'] and "
        f"obs['{key_added}_compartment']"
    )

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

    Returns:
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
        if hasattr(expr, "toarray"):
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

    Returns:
    -------
    pd.DataFrame
        Immune infiltration summary
    """
    immune_types = [
        "T_cell",
        "CD4_T",
        "CD8_T",
        "Treg",
        "B_cell",
        "NK_cell",
        "Macrophage",
        "Monocyte",
        "Neutrophil",
        "DC",
        "Mast_cell",
    ]

    df = adata.obs.copy()
    is_immune = df[cell_type_key].isin(immune_types)
    df["is_immune"] = is_immune

    if groupby is not None:
        summary = (
            df.groupby(groupby)[cell_type_key].value_counts(normalize=True).unstack(fill_value=0)
        )
    else:
        # Overall proportions
        summary = df[cell_type_key].value_counts(normalize=True).to_frame().T

    return summary
