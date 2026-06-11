"""Tumor-focused spatial transcriptomics utilities.

These functions interpret spatial data through the lens of tumor biology:
tumor-stroma boundaries, immune infiltration, spatial niches, and therapy
response signatures.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial import cKDTree

from .config import SpatialNeighborsConfig
from .neighbors import build_spatial_neighbors
from .zones import find_tissue_zones, TissueZonesConfig


def find_tumor_stroma_boundary(
    adata: AnnData,
    tumor_label: Union[str, List[str]],
    *,
    spatial_key: str = "spatial",
    label_key: str = "tumor_label",
    core_radius: Optional[float] = None,
    margin_radius: Optional[float] = None,
    key_added: str = "tumor_boundary",
) -> AnnData:
    """Label tumor core, boundary, margin, and stroma regions.

    Parameters
    ----------
    adata
        Spatial AnnData with ``obsm[spatial_key]`` and ``obs[label_key]``.
    tumor_label
        Label(s) considered tumor (e.g., "tumor", "malignant").
    spatial_key
        Key for coordinates in ``obsm``.
    label_key
        Key for region labels in ``obs``.
    core_radius
        Distance threshold (same units as coordinates) for tumor core. If
        None, estimated as the 25th percentile distance from tumor spots to
        nearest non-tumor spot.
    margin_radius
        Distance threshold for tumor margin. If None, estimated as the 75th
        percentile distance from tumor spots to nearest non-tumor spot.
    key_added
        Column name for boundary labels in ``adata.obs``.

    Returns
    -------
    AnnData
        ``adata`` with ``obs[key_added]`` and boundary statistics in ``uns``.
    """
    if spatial_key not in adata.obsm:
        raise KeyError(f"Spatial key '{spatial_key}' not found in adata.obsm")
    if label_key not in adata.obs.columns:
        raise KeyError(f"Label key '{label_key}' not found in adata.obs")

    coords = np.asarray(adata.obsm[spatial_key])[:, :2]
    labels = adata.obs[label_key].astype(str)
    tumor_labels = {tumor_label} if isinstance(tumor_label, str) else set(str(x) for x in tumor_label)
    is_tumor = labels.isin(tumor_labels).values

    if is_tumor.sum() == 0:
        raise ValueError(f"No spots found with tumor_label in {tumor_labels}")
    if (~is_tumor).sum() == 0:
        raise ValueError("All spots are labeled as tumor; cannot define boundary")

    tumor_coords = coords[is_tumor]
    non_tumor_coords = coords[~is_tumor]

    tree = cKDTree(non_tumor_coords)
    dist_to_non_tumor, _ = tree.query(tumor_coords, k=1)

    # Estimate radii if not provided
    if core_radius is None:
        core_radius = float(np.percentile(dist_to_non_tumor, 25))
    if margin_radius is None:
        margin_radius = float(np.percentile(dist_to_non_tumor, 75))

    boundary_labels = np.full(adata.n_obs, "stroma", dtype=object)
    boundary_labels[is_tumor] = "tumor"  # temporary

    tumor_indices = np.where(is_tumor)[0]
    for idx, d in zip(tumor_indices, dist_to_non_tumor):
        if d <= core_radius:
            boundary_labels[idx] = "core"
        elif d <= margin_radius:
            boundary_labels[idx] = "boundary"
        else:
            boundary_labels[idx] = "margin"

    adata.obs[key_added] = pd.Categorical(boundary_labels)

    stats = {
        "n_core": int((boundary_labels == "core").sum()),
        "n_boundary": int((boundary_labels == "boundary").sum()),
        "n_margin": int((boundary_labels == "margin").sum()),
        "n_stroma": int((boundary_labels == "stroma").sum()),
        "core_radius": core_radius,
        "margin_radius": margin_radius,
        "tumor_label": tumor_label,
        "label_key": label_key,
        "inference_level": "exploratory_spatial",
        "valid_for_publication_inference": False,
        "result_warning": (
            "Boundary labels are exploratory spatial annotations; validate "
            "with pathology or image-based annotations before publication inference."
        ),
    }

    adata.uns.setdefault("sclucid", {}).setdefault("tools", {}).setdefault("spatial", {})
    adata.uns["sclucid"]["tools"]["spatial"]["tumor_boundary"] = stats

    return adata


def compute_immune_infiltration_score(
    adata: AnnData,
    immune_markers: Optional[List[str]] = None,
    *,
    method: str = "mean_expression",
    layer: Optional[str] = None,
    key_added: str = "immune_infiltration_score",
) -> AnnData:
    """Compute a spatial immune infiltration score per spot/cell.

    Parameters
    ----------
    adata
        Spatial AnnData.
    immune_markers
        List of immune marker genes. Uses a built-in pan-immune set if None.
    method
        ``"mean_expression"`` only in this clean-room implementation.
    layer
        Expression layer. If None, uses ``adata.X``.
    key_added
        Column name in ``adata.obs``.

    Returns
    -------
    AnnData
        ``adata`` with ``obs[key_added]`` populated.
    """
    markers = immune_markers or [
        "PTPRC",
        "CD3E",
        "CD68",
        "CD79A",
        "CD8A",
        "CD4",
        "FOXP3",
    ]

    X = adata.layers[layer] if layer is not None else adata.X
    X = np.asarray(X.toarray() if hasattr(X, "toarray") else X)

    present = [g for g in markers if g in adata.var_names]
    if not present:
        score = np.zeros(adata.n_obs)
    else:
        gene_mask = adata.var_names.isin(present)
        score = X[:, gene_mask].mean(axis=1)

    adata.obs[key_added] = score
    adata.uns.setdefault("sclucid", {}).setdefault("tools", {}).setdefault("spatial", {})
    adata.uns["sclucid"]["tools"]["spatial"][key_added] = {
        "method": method,
        "markers_used": present,
        "n_markers_requested": len(markers),
        "inference_level": "exploratory_spatial",
        "valid_for_publication_inference": False,
        "result_warning": (
            "Spot-level immune scores are exploratory; they do not replace "
            "single-cell annotation or IHC quantification."
        ),
    }
    return adata


def analyze_spatial_niches(
    adata: AnnData,
    niche_config: Optional[TissueZonesConfig] = None,
    *,
    spatial_key: str = "spatial",
    key_added: str = "spatial_niches",
    n_components: int = 4,
) -> AnnData:
    """Identify tumor-relevant spatial niches.

    Wraps ``find_tissue_zones`` with tumor-focused defaults and renames zones
    to interpretable niche labels where possible.

    Parameters
    ----------
    adata
        Spatial AnnData.
    niche_config
        Optional ``TissueZonesConfig``. Defaults to expression input and
        ``n_components`` niches.
    spatial_key
        Key for spatial coordinates (currently unused but reserved).
    key_added
        Column name in ``adata.obs``.
    n_components
        Number of niches if ``niche_config`` is None.

    Returns
    -------
    AnnData
        ``adata`` with niche labels in ``obs[key_added]``.
    """
    config = niche_config or TissueZonesConfig(n_components=n_components, input="expression")
    find_tissue_zones(adata, config=config)

    # Rename the generic "tissue_zones" column to the requested key.
    zone_key = config.key_added
    adata.obs[key_added] = adata.obs[zone_key].astype(str).astype("category")

    # Compute simple marker statistics per niche to help interpretation.
    X = np.asarray(adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X)
    means = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    means["_niche"] = adata.obs[key_added].values
    niche_means = means.groupby("_niche").mean()

    adata.uns.setdefault("sclucid", {}).setdefault("tools", {}).setdefault("spatial", {})
    adata.uns["sclucid"]["tools"]["spatial"][key_added] = {
        "n_niches": int(adata.obs[key_added].nunique()),
        "niche_sizes": adata.obs[key_added].value_counts().to_dict(),
        "top_genes_per_niche": {
            str(zone): adata.var_names[np.argsort(-niche_means.loc[zone].values)[:5]].tolist()
            for zone in niche_means.index
        },
        "inference_level": "exploratory_spatial",
        "valid_for_publication_inference": False,
        "result_warning": (
            "Niche labels are data-driven and exploratory; map them to "
            "biological regions using marker inspection or pathology."
        ),
    }
    return adata


def spatial_ici_response_signature(
    adata: AnnData,
    signature_dict: Optional[Dict[str, List[str]]] = None,
    *,
    boundary_label: Optional[str] = None,
    boundary_key: str = "tumor_boundary",
    layer: Optional[str] = None,
    key_added: str = "ici_response_score",
) -> pd.DataFrame:
    """Score ICI response signatures spatially.

    Parameters
    ----------
    adata
        Spatial AnnData.
    signature_dict
        Custom signatures. If None, uses a small built-in set: cytotoxic,
        exhausted T, antigen presentation, and TLS-like.
    boundary_label
        If provided, restrict scoring to this boundary region (e.g., "margin").
    boundary_key
        Column with boundary labels.
    layer
        Expression layer.
    key_added
        Prefix for per-spot score columns.

    Returns
    -------
    pd.DataFrame
        Per-spot scores and per-region summary.
    """
    signatures = signature_dict or {
        "ICI_cytotoxic": ["GZMA", "GZMB", "PRF1", "GNLY"],
        "ICI_exhausted": ["PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT"],
        "ICI_apc": ["HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "CD74"],
        "ICI_tls": ["CD19", "CD20", "CXCL13", "CCL19", "CCL21"],
    }

    if boundary_label is not None:
        if boundary_key not in adata.obs.columns:
            raise KeyError(f"boundary_key '{boundary_key}' not found in adata.obs")
        mask = adata.obs[boundary_key] == boundary_label
        if mask.sum() == 0:
            raise ValueError(f"No spots with boundary_label='{boundary_label}'")
        subset = adata[mask].copy()
    else:
        subset = adata

    X = subset.layers[layer] if layer is not None else subset.X
    X = np.asarray(X.toarray() if hasattr(X, "toarray") else X)

    scores = {}
    for name, markers in signatures.items():
        present = [g for g in markers if g in subset.var_names]
        if not present:
            scores[f"{key_added}_{name}"] = np.nan
            continue
        gene_mask = subset.var_names.isin(present)
        scores[f"{key_added}_{name}"] = X[:, gene_mask].mean(axis=1)

    score_df = pd.DataFrame(scores, index=subset.obs_names)
    score_df["inference_level"] = "exploratory_spatial"
    score_df["valid_for_publication_inference"] = False
    score_df["method"] = "mean_expression"
    score_df["result_warning"] = (
        "Spatial ICI response signatures are exploratory enrichment scores; "
        "they do not establish clinical predictive validity."
    )

    if boundary_label is not None:
        score_df["region"] = boundary_label

    return score_df
