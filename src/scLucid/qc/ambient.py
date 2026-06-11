"""Python-native ambient RNA risk diagnostics.

This module deliberately diagnoses risk instead of correcting expression values.
Correction backends such as CellBender/scAR can be connected later without adding
R dependencies to the default QC path.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import scipy.sparse as sparse
from anndata import AnnData

from ..utils import sanitize_for_hdf5


def _matrix_from_layer(adata: AnnData, layer: Optional[str]):
    if layer is not None:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in adata.layers")
        return adata.layers[layer]
    return adata.X


def _row_sums(X) -> np.ndarray:
    return np.asarray(X.sum(axis=1)).ravel() if sparse.issparse(X) else np.asarray(X.sum(axis=1)).ravel()


def _col_sums(X) -> np.ndarray:
    return np.asarray(X.sum(axis=0)).ravel() if sparse.issparse(X) else np.asarray(X.sum(axis=0)).ravel()


def _ambient_risk_method_metadata() -> Dict[str, Any]:
    return {
        "calibration_status": "heuristic_unvalidated",
        "risk_score_weights": {
            "top_gene_dominance": 0.4,
            "low_count_enrichment": 0.4,
            "enriched_gene_breadth": 0.2,
        },
        "risk_score_note": (
            "Risk score is a Python-native diagnostic heuristic and is not calibrated "
            "as a cross-tissue ambient RNA contamination probability."
        ),
    }


def diagnose_ambient_rna(
    adata: AnnData,
    *,
    layer: Optional[str] = None,
    low_count_quantile: float = 0.1,
    top_n_genes: int = 20,
    dominance_threshold: float = 0.35,
    low_count_enrichment_threshold: float = 1.5,
    min_cells: int = 20,
) -> Dict[str, Any]:
    """Estimate ambient RNA contamination risk from count-distribution patterns.

    The diagnostic is Python-only and heuristic. It is intended to decide whether
    users should inspect or run an external correction backend, not to replace
    CellBender/scAR/SoupX-like correction.
    """
    if adata.n_obs == 0 or adata.n_vars == 0:
        return {
            "available": False,
            "diagnostic_only": True,
            **_ambient_risk_method_metadata(),
            "risk_level": "unknown",
            "risk_score": 0.0,
            "reason": "empty_adata",
            "method_note": "Ambient RNA diagnostic only; no expression correction was applied.",
        }

    X = _matrix_from_layer(adata, layer)
    totals = _row_sums(X)
    positive_totals = totals[totals > 0]
    if positive_totals.size < min_cells:
        return {
            "available": False,
            "diagnostic_only": True,
            **_ambient_risk_method_metadata(),
            "risk_level": "unknown",
            "risk_score": 0.0,
            "reason": "too_few_nonzero_cells",
            "n_nonzero_cells": int(positive_totals.size),
            "method_note": "Ambient RNA diagnostic only; no expression correction was applied.",
        }

    threshold = float(np.quantile(positive_totals, low_count_quantile))
    low_mask = totals <= threshold
    low_mask &= totals > 0
    if int(low_mask.sum()) < max(5, min_cells // 4):
        return {
            "available": False,
            "diagnostic_only": True,
            **_ambient_risk_method_metadata(),
            "risk_level": "unknown",
            "risk_score": 0.0,
            "reason": "too_few_low_count_cells",
            "low_count_threshold": threshold,
            "n_low_count_cells": int(low_mask.sum()),
            "method_note": "Ambient RNA diagnostic only; no expression correction was applied.",
        }

    gene_totals = _col_sums(X).astype(float)
    low_gene_totals = _col_sums(X[low_mask]).astype(float)
    total_counts = float(gene_totals.sum())
    low_total_counts = float(low_gene_totals.sum())
    if total_counts <= 0 or low_total_counts <= 0:
        return {
            "available": False,
            "diagnostic_only": True,
            **_ambient_risk_method_metadata(),
            "risk_level": "unknown",
            "risk_score": 0.0,
            "reason": "zero_library_after_subset",
            "method_note": "Ambient RNA diagnostic only; no expression correction was applied.",
        }

    top_n = int(min(max(1, top_n_genes), adata.n_vars))
    top_idx = np.argsort(gene_totals)[-top_n:][::-1]
    global_frac = gene_totals[top_idx] / total_counts
    low_frac = low_gene_totals[top_idx] / low_total_counts
    enrichment = np.divide(low_frac, global_frac, out=np.zeros_like(low_frac), where=global_frac > 0)

    top_gene_fraction = float(global_frac.sum())
    low_count_top_gene_fraction = float(low_frac.sum())
    median_enrichment = float(np.median(enrichment)) if enrichment.size else 0.0
    max_enrichment = float(np.max(enrichment)) if enrichment.size else 0.0
    n_enriched_top_genes = int((enrichment >= low_count_enrichment_threshold).sum())

    dominance_score = min(1.0, top_gene_fraction / dominance_threshold)
    enrichment_score = min(1.0, median_enrichment / low_count_enrichment_threshold)
    breadth_score = min(1.0, n_enriched_top_genes / max(1.0, top_n / 2.0))
    risk_score = float(0.4 * dominance_score + 0.4 * enrichment_score + 0.2 * breadth_score)
    if risk_score >= 0.7:
        risk_level = "high"
    elif risk_score >= 0.4:
        risk_level = "moderate"
    else:
        risk_level = "low"

    top_genes = [
        {
            "gene": str(adata.var_names[i]),
            "global_fraction": float(gene_totals[i] / total_counts),
            "low_count_fraction": float(low_gene_totals[i] / low_total_counts),
            "low_count_enrichment": float(enrichment[pos]),
        }
        for pos, i in enumerate(top_idx)
    ]

    correction_status = adata.uns.get("sclucid", {}).get("qc", {}).get(
        "ambient_correction_status", {}
    )
    if not correction_status:
        correction_status = {
            "corrected": False,
            "backend": None,
            "note": "No ambient RNA correction metadata detected in scLucid QC namespace.",
        }

    return {
        "available": True,
        "diagnostic_only": True,
        **_ambient_risk_method_metadata(),
        "method": "python_heuristic_low_count_enrichment",
        "method_note": "Ambient RNA diagnostic only; no expression correction was applied.",
        "layer": layer or "X",
        "risk_level": risk_level,
        "risk_score": risk_score,
        "low_count_quantile": float(low_count_quantile),
        "low_count_threshold": threshold,
        "n_low_count_cells": int(low_mask.sum()),
        "top_n_genes": top_n,
        "top_gene_fraction": top_gene_fraction,
        "low_count_top_gene_fraction": low_count_top_gene_fraction,
        "median_low_count_enrichment": median_enrichment,
        "max_low_count_enrichment": max_enrichment,
        "n_enriched_top_genes": n_enriched_top_genes,
        "top_genes": top_genes,
        "correction_status": correction_status,
        "recommendation": (
            "Inspect ambient RNA and consider Python backends such as CellBender or scAR."
            if risk_level in {"moderate", "high"}
            else "Ambient RNA risk appears low by this heuristic diagnostic."
        ),
    }


def record_ambient_correction_status(
    adata: AnnData,
    *,
    corrected: bool,
    backend: str,
    output_layer: Optional[str] = None,
    details: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Record externally performed ambient RNA correction metadata."""
    status = {
        "corrected": bool(corrected),
        "backend": backend,
        "output_layer": output_layer,
        "details": details or {},
    }
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {})[
        "ambient_correction_status"
    ] = status
    return status


def register_external_ambient_result(
    adata: AnnData,
    *,
    backend: str,
    source_path: Optional[str] = None,
    corrected_adata: Optional[AnnData] = None,
    corrected_layer: Optional[str] = None,
    output_layer: str = "ambient_corrected",
    details: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Register externally corrected counts from CellBender/scAR/SoupX-like tools.

    This function intentionally accepts already-loaded AnnData objects to avoid
    imposing a single file format. If ``corrected_adata`` is supplied, it must
    have matching cell and gene names and its selected matrix is copied into
    ``adata.layers[output_layer]``.
    """
    copied_matrix = False
    if corrected_adata is not None:
        if not corrected_adata.obs_names.equals(adata.obs_names):
            raise ValueError("corrected_adata.obs_names must match adata.obs_names")
        if not corrected_adata.var_names.equals(adata.var_names):
            raise ValueError("corrected_adata.var_names must match adata.var_names")
        corrected_X = _matrix_from_layer(corrected_adata, corrected_layer)
        adata.layers[output_layer] = corrected_X.copy()
        copied_matrix = True

    status = record_ambient_correction_status(
        adata,
        corrected=True,
        backend=backend,
        output_layer=output_layer if copied_matrix else corrected_layer,
        details={
            "source_path": str(Path(source_path)) if source_path else None,
            "matrix_copied": copied_matrix,
            "corrected_layer": corrected_layer,
            **(details or {}),
        },
    )
    return status


def diagnose_empty_droplets(
    adata: AnnData,
    *,
    layer: Optional[str] = None,
    cell_call_key: Optional[str] = None,
    empty_count_quantile: float = 0.1,
    top_n_genes: int = 20,
    min_barcodes: int = 100,
) -> Dict[str, Any]:
    """Diagnose empty-droplet/background profile from a raw 10x-like matrix.

    This is a diagnostic interface, not an EmptyDrops replacement. It reports
    whether raw low-count barcodes are available and summarizes the putative
    ambient/background profile for review.
    """
    if adata.n_obs < min_barcodes or adata.n_vars == 0:
        return {
            "available": False,
            "diagnostic_only": True,
            "risk_level": "unknown",
            "reason": "too_few_barcodes_for_empty_droplet_diagnostic",
            "n_barcodes": int(adata.n_obs),
            "method_note": "Empty-droplet diagnostic only; not an EmptyDrops replacement.",
        }

    X = _matrix_from_layer(adata, layer)
    totals = _row_sums(X).astype(float)
    positive = totals[totals > 0]
    if positive.size < min_barcodes:
        return {
            "available": False,
            "diagnostic_only": True,
            "risk_level": "unknown",
            "reason": "too_few_nonzero_barcodes",
            "n_nonzero_barcodes": int(positive.size),
            "method_note": "Empty-droplet diagnostic only; not an EmptyDrops replacement.",
        }

    if cell_call_key and cell_call_key in adata.obs.columns:
        calls = adata.obs[cell_call_key]
        empty_mask = calls.astype(str).str.lower().isin(
            {"empty", "background", "ambient", "false", "0"}
        )
        called_cell_mask = ~empty_mask
    else:
        threshold = float(np.quantile(positive, empty_count_quantile))
        empty_mask = (totals > 0) & (totals <= threshold)
        called_cell_mask = totals > threshold

    empty_mask = np.asarray(empty_mask, dtype=bool)
    called_cell_mask = np.asarray(called_cell_mask, dtype=bool)
    n_empty = int(np.asarray(empty_mask).sum())
    n_called = int(np.asarray(called_cell_mask).sum())
    if n_empty < 10:
        return {
            "available": False,
            "diagnostic_only": True,
            "risk_level": "unknown",
            "reason": "too_few_putative_empty_droplets",
            "n_putative_empty_droplets": n_empty,
            "method_note": "Empty-droplet diagnostic only; not an EmptyDrops replacement.",
        }

    empty_gene_totals = _col_sums(X[empty_mask]).astype(float)
    all_gene_totals = _col_sums(X).astype(float)
    empty_total = float(empty_gene_totals.sum())
    all_total = float(all_gene_totals.sum())
    if empty_total <= 0 or all_total <= 0:
        return {
            "available": False,
            "diagnostic_only": True,
            "risk_level": "unknown",
            "reason": "zero_background_counts",
            "method_note": "Empty-droplet diagnostic only; not an EmptyDrops replacement.",
        }

    top_n = min(int(top_n_genes), adata.n_vars)
    top_idx = np.argsort(empty_gene_totals)[-top_n:][::-1]
    background_fraction = empty_gene_totals[top_idx] / empty_total
    global_fraction = all_gene_totals[top_idx] / all_total
    enrichment = np.divide(
        background_fraction,
        global_fraction,
        out=np.zeros_like(background_fraction),
        where=global_fraction > 0,
    )
    top_background_fraction = float(background_fraction.sum())
    barcode_rank_gap = np.nan
    sorted_totals = np.sort(positive)[::-1]
    if sorted_totals.size >= 20:
        log_totals = np.log10(sorted_totals + 1)
        diffs = np.diff(log_totals)
        barcode_rank_gap = float(abs(diffs.min())) if diffs.size else np.nan

    risk_score = float(min(1.0, 0.6 * top_background_fraction / 0.5 + 0.4 * np.nan_to_num(barcode_rank_gap, nan=0.0)))
    risk_level = "high" if risk_score >= 0.7 else "moderate" if risk_score >= 0.4 else "low"

    result = {
        "available": True,
        "diagnostic_only": True,
        "method": "python_barcode_rank_background_profile",
        "method_note": "Empty-droplet diagnostic only; not an EmptyDrops replacement.",
        "layer": layer or "X",
        "risk_level": risk_level,
        "risk_score": risk_score,
        "n_barcodes": int(adata.n_obs),
        "n_called_cells": n_called,
        "n_putative_empty_droplets": n_empty,
        "empty_count_quantile": float(empty_count_quantile),
        "barcode_rank_gap": barcode_rank_gap,
        "top_background_fraction": top_background_fraction,
        "top_background_genes": [
            {
                "gene": str(adata.var_names[i]),
                "background_fraction": float(empty_gene_totals[i] / empty_total),
                "global_fraction": float(all_gene_totals[i] / all_total),
                "background_enrichment": float(enrichment[pos]),
            }
            for pos, i in enumerate(top_idx)
        ],
        "recommendation": (
            "Raw matrix contains substantial background signal; consider CellBender/scAR "
            "or a validated empty-droplet workflow before downstream interpretation."
            if risk_level in {"moderate", "high"}
            else "Empty-droplet background risk appears low by this barcode-rank diagnostic."
        ),
    }
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {})[
        "empty_droplet_summary"
    ] = sanitize_for_hdf5(result)
    return result
