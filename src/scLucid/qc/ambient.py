"""Python-native ambient RNA risk diagnostics and lightweight correction.

This module deliberately starts with diagnostics.  A zero-dependency linear
background-subtraction correction is provided for users who need a fast,
interpretable correction without installing CellBender/scAR.  External backends
(CellBender/scAR/SoupX) can still be connected via
:func:`register_external_ambient_result`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np
import scipy.sparse as sparse
from anndata import AnnData

from ..utils import sanitize_for_hdf5
from .artifacts import record_qc_artifact_contract

AMBIENT_CORRECTED_COUNTS_LAYER = "ambient_corrected_counts"


def _matrix_from_layer(adata: AnnData, layer: Optional[str]):
    if layer is not None:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in adata.layers")
        return adata.layers[layer]
    return adata.X


def _row_sums(X) -> np.ndarray:
    return (
        np.asarray(X.sum(axis=1)).ravel()
        if sparse.issparse(X)
        else np.asarray(X.sum(axis=1)).ravel()
    )


def _col_sums(X) -> np.ndarray:
    return (
        np.asarray(X.sum(axis=0)).ravel()
        if sparse.issparse(X)
        else np.asarray(X.sum(axis=0)).ravel()
    )


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


def infer_ambient_input_context(
    adata: AnnData,
    *,
    layer: Optional[str] = None,
    matrix_source: str = "auto",
    cell_call_key: Optional[str] = None,
) -> Dict[str, Any]:
    """Infer whether the matrix is raw-like or filtered-like for ambient QC.

    The context is used to avoid overstating what can be learned from a
    filtered-feature matrix. Raw matrices with empty/background barcodes can
    support empty-droplet/background modeling; filtered matrices should default
    to diagnostic/review or externally registered correction.
    """
    candidate_call_keys = [
        key
        for key in [
            cell_call_key,
            "likely_empty_droplet",
            "empty_droplet",
            "barcode_type",
            "cell_call",
        ]
        if key and key in adata.obs
    ]
    has_empty_labels = False
    for key in candidate_call_keys:
        labels = adata.obs[key].astype(str).str.lower()
        if labels.isin({"empty", "background", "ambient", "false", "0", "true"}).any():
            has_empty_labels = True
            break
    matrix_type = matrix_source
    if matrix_source == "auto":
        matrix_type = "raw_like" if has_empty_labels else "filtered_like"

    if matrix_type == "raw_like":
        suitable_backends = ["cellbender", "soupx", "scar", "linear_background_subtraction"]
        correction_recommendation = (
            "Raw-like matrix detected; model-based ambient correction can be considered "
            "before formal normalization."
        )
    elif matrix_type == "filtered_like":
        suitable_backends = ["external_filtered_matrix_correction", "diagnostic_only"]
        correction_recommendation = (
            "Filtered matrix detected; CellBender/empty-droplet modeling is limited. "
            "Prefer diagnostic review or externally registered DecontX/SoupX-like output."
        )
    else:
        suitable_backends = ["external_registration", "diagnostic_only"]
        correction_recommendation = (
            "Ambient matrix source is user-declared; record method details before using "
            "corrected counts downstream."
        )

    return {
        "schema_version": "ambient_input_context_v1",
        "matrix_source": matrix_source,
        "matrix_type": matrix_type,
        "layer": layer or "X",
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "cell_call_keys_detected": candidate_call_keys,
        "has_empty_droplet_labels": has_empty_labels,
        "counts_layer_present": "counts" in adata.layers,
        "suitable_backends": suitable_backends,
        "correction_recommendation": correction_recommendation,
    }


def build_ambient_layer_contract(
    adata: AnnData,
    *,
    input_context: Optional[Dict[str, Any]] = None,
    correction_summary: Optional[Dict[str, Any]] = None,
    output_layer: str = AMBIENT_CORRECTED_COUNTS_LAYER,
) -> Dict[str, Any]:
    """Build the ambient layer contract for downstream preprocessing."""
    correction_summary = correction_summary or {}
    corrected = bool(correction_summary.get("corrected", False))
    corrected_layer = correction_summary.get("output_layer") or output_layer
    corrected_layer_present = bool(corrected_layer in adata.layers)
    counts_layer = "counts" if "counts" in adata.layers else None
    recommended_layer = (
        corrected_layer
        if corrected and corrected_layer_present
        else counts_layer
        if counts_layer is not None
        else None
    )
    return {
        "schema_version": "ambient_layer_contract_v1",
        "input_context": input_context or infer_ambient_input_context(adata),
        "correction": correction_summary,
        "counts_layer": counts_layer,
        "corrected_layer": corrected_layer if corrected else None,
        "corrected_layer_present": corrected_layer_present,
        "recommended_preprocess_counts_layer": recommended_layer,
        "canonical_flow": (
            "counts -> ambient_corrected_counts(optional) -> normalized -> raw -> HVG -> scaled -> PCA -> graph"
        ),
        "review_required": bool(
            correction_summary.get("review_required", False)
            or (corrected and not corrected_layer_present)
        ),
        "risk_note": (
            correction_summary.get("risk_note")
            or "Ambient correction was not applied; use counts for preprocessing unless external corrected counts are registered."
        ),
    }


def record_ambient_layer_contract(
    adata: AnnData,
    *,
    input_context: Optional[Dict[str, Any]] = None,
    correction_summary: Optional[Dict[str, Any]] = None,
    output_layer: str = AMBIENT_CORRECTED_COUNTS_LAYER,
) -> Dict[str, Any]:
    """Store the ambient layer contract in the scLucid QC namespace."""
    contract = build_ambient_layer_contract(
        adata,
        input_context=input_context,
        correction_summary=correction_summary,
        output_layer=output_layer,
    )
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {})["ambient_layer_contract"] = (
        sanitize_for_hdf5(contract)
    )
    record_qc_artifact_contract(adata)
    return contract


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
    input_context = infer_ambient_input_context(adata, layer=layer)
    if adata.n_obs == 0 or adata.n_vars == 0:
        return {
            "available": False,
            "diagnostic_only": True,
            "input_context": input_context,
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
    enrichment = np.divide(
        low_frac, global_frac, out=np.zeros_like(low_frac), where=global_frac > 0
    )

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

    correction_status = (
        adata.uns.get("sclucid", {}).get("qc", {}).get("ambient_correction_status", {})
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
        "input_context": input_context,
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
            input_context["correction_recommendation"]
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
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {})["ambient_correction_status"] = status
    return status


def register_external_ambient_result(
    adata: AnnData,
    *,
    backend: str,
    source_path: Optional[str] = None,
    corrected_adata: Optional[AnnData] = None,
    corrected_layer: Optional[str] = None,
    obs_column_map: Optional[Dict[str, str]] = None,
    output_layer: str = AMBIENT_CORRECTED_COUNTS_LAYER,
    details: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Register externally corrected counts from CellBender/scAR/SoupX-like tools.

    This function intentionally accepts already-loaded AnnData objects to avoid
    imposing a single file format. If ``corrected_adata`` is supplied, it must
    have matching cell and gene names and its selected matrix is copied into
    ``adata.layers[output_layer]``.
    """
    copied_matrix = False
    copied_obs_columns: Dict[str, str] = {}
    if corrected_adata is not None:
        if not corrected_adata.obs_names.equals(adata.obs_names):
            raise ValueError("corrected_adata.obs_names must match adata.obs_names")
        if not corrected_adata.var_names.equals(adata.var_names):
            raise ValueError("corrected_adata.var_names must match adata.var_names")
        corrected_X = _matrix_from_layer(corrected_adata, corrected_layer)
        adata.layers[output_layer] = corrected_X.copy()
        copied_matrix = True
        for canonical, source in (obs_column_map or {}).items():
            if source not in corrected_adata.obs:
                raise KeyError(f"corrected_adata.obs column '{source}' not found")
            adata.obs[canonical] = corrected_adata.obs[source].reindex(adata.obs_names).to_numpy()
            copied_obs_columns[canonical] = source
    elif obs_column_map:
        for canonical, source in obs_column_map.items():
            if source not in adata.obs:
                raise KeyError(f"adata.obs column '{source}' not found")
            if canonical != source:
                adata.obs[canonical] = adata.obs[source].to_numpy()
            copied_obs_columns[canonical] = source

    status = record_ambient_correction_status(
        adata,
        corrected=True,
        backend=backend,
        output_layer=output_layer if copied_matrix else corrected_layer,
        details={
            "source_path": str(Path(source_path)) if source_path else None,
            "matrix_copied": copied_matrix,
            "corrected_layer": corrected_layer,
            "obs_column_map": copied_obs_columns,
            **(details or {}),
        },
    )
    record_ambient_layer_contract(
        adata,
        input_context=infer_ambient_input_context(
            adata,
            layer=output_layer if copied_matrix else corrected_layer,
            matrix_source="external",
        ),
        correction_summary=status,
        output_layer=output_layer,
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
    record: bool = False,
    key_added: str = "empty_droplet_summary",
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
        empty_mask = (
            calls.astype(str).str.lower().isin({"empty", "background", "ambient", "false", "0"})
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

    risk_score = float(
        min(
            1.0,
            0.6 * top_background_fraction / 0.5 + 0.4 * np.nan_to_num(barcode_rank_gap, nan=0.0),
        )
    )
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
    if record:
        adata.uns.setdefault("sclucid", {}).setdefault("qc", {})[key_added] = sanitize_for_hdf5(
            result
        )
    return result


def _estimate_cell_contamination(
    X,
    background_profile: np.ndarray,
    ambient_gene_mask: np.ndarray,
    robust_quantile: float = 0.5,
    chunk_size: int = 4096,
) -> Tuple[np.ndarray, np.ndarray]:
    """Estimate per-cell ambient contamination fraction ``rho_i``.

    For each cell, the observed counts in high-ambient genes are compared to
    the expected counts from the background profile.  The contamination
    fraction is the robust median ratio, clipped to ``[0, 1]``.

    Parameters
    ----------
    X
        Count matrix (sparse or dense), same orientation as ``background_profile``.
    background_profile
        Per-gene ambient expectation (length ``n_vars``).
    ambient_gene_mask
        Boolean mask of genes used for the estimate.
    robust_quantile
        Quantile used for the robust ratio summary per cell.

    Returns:
    -------
    rho : np.ndarray
        Per-cell contamination fraction, clipped to ``[0, 1]``.
    expected_counts : np.ndarray
        Per-cell total expected ambient counts (used for diagnostics).
    """
    X = _ensure_array(X)
    bg = np.asarray(background_profile, dtype=float)
    ambient_idx = np.flatnonzero(ambient_gene_mask)
    if len(ambient_idx) == 0:
        return np.zeros(X.shape[0]), np.zeros(X.shape[0])

    # Observed counts in ambient genes
    observed = np.asarray(X[:, ambient_idx].sum(axis=1)).ravel().astype(float)
    # Expected ambient counts per gene per cell (background profile scaled to cell library)
    cell_totals = np.asarray(X.sum(axis=1)).ravel().astype(float)
    bg_total = float(bg[ambient_idx].sum())
    if bg_total <= 0:
        return np.zeros(X.shape[0]), np.zeros(X.shape[0])

    ambient_weights = bg[ambient_idx] / bg_total
    expected_counts = cell_totals.copy()
    rho = np.zeros(X.shape[0], dtype=float)

    # Per-cell, per-gene ratio; summarise robustly in row chunks to avoid
    # materialising a full n_cells x n_ambient_genes dense matrix.
    q = robust_quantile * 100
    for start in range(0, X.shape[0], max(1, chunk_size)):
        stop = min(start + max(1, chunk_size), X.shape[0])
        if sparse.issparse(X):
            observed_chunk = X[start:stop, :][:, ambient_idx].toarray()
        else:
            observed_chunk = np.asarray(X[start:stop, :][:, ambient_idx], dtype=float)
        expected_chunk = cell_totals[start:stop, None] * ambient_weights[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio_chunk = observed_chunk / np.maximum(expected_chunk, 1e-12)
            ratio_chunk = np.where(np.isfinite(ratio_chunk), ratio_chunk, 0.0)
            rho[start:stop] = np.percentile(ratio_chunk, q, axis=1)

    rho = np.nan_to_num(rho, nan=0.0, posinf=1.0, neginf=0.0)
    rho = np.clip(rho, 0.0, 1.0)
    # If observed ambient counts are zero, set rho to zero
    rho[observed == 0] = 0.0
    return rho, expected_counts


def _residual_ambient_score(
    X,
    corrected_X,
    ambient_gene_mask: np.ndarray,
    low_count_mask: np.ndarray,
) -> float:
    """Compute a lightweight residual ambient-RNA score after correction.

    The score is the median, across low-count cells, of the fraction of counts
    that remain in the ambient-marker genes relative to the cell library.  A
    high score suggests that additional model-based correction may be warranted.
    """
    X = _ensure_array(X)
    corrected_X = _ensure_array(corrected_X)
    if low_count_mask.sum() == 0 or ambient_gene_mask.sum() == 0:
        return 0.0

    ambient_idx = np.flatnonzero(ambient_gene_mask)
    after = np.asarray(corrected_X[low_count_mask, :][:, ambient_idx].sum(axis=1)).ravel()
    totals = np.asarray(X[low_count_mask, :].sum(axis=1)).ravel()
    with np.errstate(divide="ignore", invalid="ignore"):
        residual_frac = np.divide(
            after,
            totals,
            out=np.zeros_like(after, dtype=float),
            where=totals > 0,
        )
    return float(np.median(residual_frac))


def _ensure_array(X):
    """Return a CSR sparse or dense array with a consistent interface."""
    if sparse.issparse(X):
        return X.tocsr(copy=False)
    return np.asarray(X)


def correct_ambient_rna_linear(
    adata: AnnData,
    *,
    layer: Optional[str] = None,
    output_layer: str = AMBIENT_CORRECTED_COUNTS_LAYER,
    empty_droplet_key: Optional[str] = None,
    empty_count_quantile: float = 0.1,
    ambient_enrichment_threshold: float = 1.5,
    min_background_counts: int = 10,
    robust_quantile: float = 0.5,
    max_removed_fraction_per_gene: float = 0.95,
    record: bool = True,
    key_added: str = "ambient_correction_summary",
) -> Dict[str, Any]:
    """Lightweight ambient RNA correction by linear background subtraction.

    This is a fast, zero-dependency correction that uses the empty-droplet
    profile estimated from the raw matrix.  It is intended as a practical
    alternative to CellBender/scAR when only a quick, interpretable correction
    is needed.  The correction is conservative: contamination fractions are
    clipped to ``[0, 1]`` and counts are never allowed to go negative.

    Parameters
    ----------
    adata
        AnnData object containing raw counts (ideally including empty droplets).
    layer
        Layer to correct.  Defaults to ``adata.X``.
    output_layer
        Layer name for the corrected matrix.
    empty_droplet_key
        Optional obs column with empty-droplet labels.  If provided, cells
        labelled as empty/background are used directly.  Otherwise empty
        droplets are defined by ``empty_count_quantile``.
    empty_count_quantile
        Quantile used to define putative empty droplets when
        ``empty_droplet_key`` is not provided.
    ambient_enrichment_threshold
        Minimum empty-droplet / global enrichment for a gene to be used in the
        contamination estimate.
    min_background_counts
        Minimum total counts in empty droplets required to proceed.
    robust_quantile
        Quantile used to summarise per-cell contamination ratios.
    max_removed_fraction_per_gene
        Maximum fraction of a gene's counts that may be removed from any single
        cell.  Prevents over-subtraction for highly expressed genes.
    record
        If True, store the correction summary in ``adata.uns['sclucid']['qc']``.
    key_added
        Key used when ``record=True``.

    Returns:
    -------
    dict
        Correction summary with ``corrected``, ``output_layer``, ``n_cells``,
        ``removed_counts``, ``mean_rho``, ``top_corrected_genes``,
        ``review_required``, ``risk_note``.
    """
    X = _matrix_from_layer(adata, layer)
    totals = _row_sums(X).astype(float)
    positive = totals[totals > 0]
    if len(positive) == 0:
        result = {
            "corrected": False,
            "output_layer": output_layer,
            "reason": "no_nonzero_barcodes",
            "review_required": True,
            "risk_note": "No nonzero barcodes found; cannot estimate ambient profile.",
        }
        if record:
            adata.uns.setdefault("sclucid", {}).setdefault("qc", {})[key_added] = sanitize_for_hdf5(
                result
            )
            record_ambient_correction_status(
                adata,
                corrected=False,
                backend="linear_background_subtraction",
                output_layer=output_layer,
                details=result,
            )
            record_ambient_layer_contract(
                adata,
                input_context=infer_ambient_input_context(adata, layer=layer),
                correction_summary=result,
                output_layer=output_layer,
            )
        return result

    # Define empty droplets
    if empty_droplet_key and empty_droplet_key in adata.obs.columns:
        calls = adata.obs[empty_droplet_key]
        empty_mask = (
            calls.astype(str).str.lower().isin({"empty", "background", "ambient", "false", "0"})
        )
    else:
        threshold = float(np.quantile(positive, empty_count_quantile))
        empty_mask = (totals > 0) & (totals <= threshold)

    empty_mask = np.asarray(empty_mask, dtype=bool)
    n_empty = int(empty_mask.sum())
    if n_empty < 10:
        result = {
            "corrected": False,
            "output_layer": output_layer,
            "reason": "too_few_empty_droplets",
            "n_putative_empty_droplets": n_empty,
            "review_required": True,
            "risk_note": "Too few empty droplets for reliable ambient correction.",
        }
        if record:
            adata.uns.setdefault("sclucid", {}).setdefault("qc", {})[key_added] = sanitize_for_hdf5(
                result
            )
            record_ambient_correction_status(
                adata,
                corrected=False,
                backend="linear_background_subtraction",
                output_layer=output_layer,
                details=result,
            )
            record_ambient_layer_contract(
                adata,
                input_context=infer_ambient_input_context(adata, layer=layer),
                correction_summary=result,
                output_layer=output_layer,
            )
        return result

    empty_gene_totals = _col_sums(X[empty_mask]).astype(float)
    all_gene_totals = _col_sums(X).astype(float)
    empty_total = float(empty_gene_totals.sum())
    if empty_total < min_background_counts:
        result = {
            "corrected": False,
            "output_layer": output_layer,
            "reason": "insufficient_background_counts",
            "empty_total_counts": empty_total,
            "review_required": True,
            "risk_note": "Background counts too low for reliable correction.",
        }
        if record:
            adata.uns.setdefault("sclucid", {}).setdefault("qc", {})[key_added] = sanitize_for_hdf5(
                result
            )
            record_ambient_correction_status(
                adata,
                corrected=False,
                backend="linear_background_subtraction",
                output_layer=output_layer,
                details=result,
            )
            record_ambient_layer_contract(
                adata,
                input_context=infer_ambient_input_context(adata, layer=layer),
                correction_summary=result,
                output_layer=output_layer,
            )
        return result

    background_fraction = empty_gene_totals / empty_total
    global_fraction = np.divide(
        all_gene_totals,
        float(all_gene_totals.sum()),
        out=np.zeros_like(all_gene_totals, dtype=float),
        where=all_gene_totals > 0,
    )
    enrichment = np.divide(
        background_fraction,
        global_fraction,
        out=np.zeros_like(background_fraction),
        where=global_fraction > 0,
    )

    # Ambient genes: enriched in empty droplets and present in cells
    ambient_gene_mask = enrichment >= ambient_enrichment_threshold
    if ambient_gene_mask.sum() < 3:
        # Relax threshold if too few genes pass
        top_n = max(3, min(20, adata.n_vars // 10))
        ambient_gene_mask = np.zeros_like(enrichment, dtype=bool)
        ambient_gene_mask[np.argsort(enrichment)[-top_n:]] = True

    rho, expected_counts = _estimate_cell_contamination(
        X, background_fraction, ambient_gene_mask, robust_quantile=robust_quantile
    )

    # Build correction: subtract expected ambient counts per gene per cell
    bg = np.asarray(background_fraction, dtype=float)
    bg_total = float(bg.sum())
    if bg_total <= 0:
        result = {
            "corrected": False,
            "output_layer": output_layer,
            "reason": "zero_background_profile",
            "review_required": True,
            "risk_note": "Background profile sums to zero; cannot correct.",
        }
        if record:
            adata.uns.setdefault("sclucid", {}).setdefault("qc", {})[key_added] = sanitize_for_hdf5(
                result
            )
            record_ambient_correction_status(
                adata,
                corrected=False,
                backend="linear_background_subtraction",
                output_layer=output_layer,
                details=result,
            )
            record_ambient_layer_contract(
                adata,
                input_context=infer_ambient_input_context(adata, layer=layer),
                correction_summary=result,
                output_layer=output_layer,
            )
        return result

    cell_totals = np.asarray(X.sum(axis=1)).ravel().astype(float)
    expected_cell_scale = rho * cell_totals / bg_total

    # Shrink per-gene subtraction so we never remove more than the configured
    # fraction of the original count in that cell. Process in chunks to avoid
    # materialising a full n_cells x n_genes dense matrix for large datasets.
    chunk_cells = max(1, min(adata.n_obs, int(50_000_000 // max(adata.n_vars, 1))))
    if sparse.issparse(X):
        corrected = X.astype(float).copy()
        for start in range(0, adata.n_obs, chunk_cells):
            stop = min(start + chunk_cells, adata.n_obs)
            chunk = corrected[start:stop, :]
            row, col = chunk.nonzero()
            original = chunk.data
            # row is relative to the chunk, so offset by start for expected_cell_scale
            subtract = expected_cell_scale[start + row] * bg[col]
            max_remove = original * max_removed_fraction_per_gene
            chunk.data = np.maximum(original - np.minimum(subtract, max_remove), 0.0)
            # Keep sparse structure integer-like by zeroing values below 1.
            chunk.data[chunk.data < 1.0] = 0.0
            chunk.eliminate_zeros()
    else:
        X_arr = np.asarray(X, dtype=float)
        corrected = np.empty_like(X_arr, dtype=float)
        chunk_size = max(1, int(20_000_000 // max(adata.n_vars, 1)))
        for start in range(0, adata.n_obs, chunk_size):
            stop = min(start + chunk_size, adata.n_obs)
            expected_chunk = expected_cell_scale[start:stop, None] * bg[None, :]
            max_remove = X_arr[start:stop, :] * max_removed_fraction_per_gene
            corrected[start:stop, :] = np.maximum(
                X_arr[start:stop, :] - np.minimum(expected_chunk, max_remove),
                0.0,
            )
        corrected[corrected < 1.0] = 0.0

    adata.layers[output_layer] = corrected

    removed_counts = float(np.asarray(X.sum() - corrected.sum()))
    mean_rho = float(np.mean(rho))

    # Top corrected genes by absolute removed counts
    gene_removed = np.asarray(X.sum(axis=0)).ravel() - np.asarray(corrected.sum(axis=0)).ravel()
    top_idx = np.argsort(-gene_removed)[:10]
    top_corrected_genes = [
        {
            "gene": str(adata.var_names[i]),
            "removed_counts": float(gene_removed[i]),
            "background_fraction": float(background_fraction[i]),
        }
        for i in top_idx
        if gene_removed[i] > 0
    ]

    review_required = mean_rho > 0.3 or removed_counts > float(X.sum()) * 0.2
    risk_note = (
        "High estimated contamination fraction; review corrected counts before downstream analysis."
        if review_required
        else (
            "Conservative linear fallback applied for review. Prefer validated "
            "model-based or external ambient correction when raw/background data are available."
        )
    )

    # Residual ambient score on low-count cells
    low_count_threshold = float(np.quantile(positive, empty_count_quantile))
    low_count_mask = (totals > 0) & (totals <= low_count_threshold)
    residual_score = _residual_ambient_score(X, corrected, ambient_gene_mask, low_count_mask)

    result = {
        "corrected": True,
        "output_layer": output_layer,
        "method": "linear_background_subtraction",
        "calibration_status": "conservative_fallback_not_model_based",
        "diagnostic_only": False,
        "correction_note": (
            "This Python-native correction is a conservative fallback and is not "
            "equivalent to CellBender, SoupX, DecontX, or scAR."
        ),
        "n_cells": int(adata.n_obs),
        "n_putative_empty_droplets": n_empty,
        "n_ambient_genes": int(ambient_gene_mask.sum()),
        "removed_counts": removed_counts,
        "removed_fraction": float(removed_counts / max(float(X.sum()), 1.0)),
        "mean_rho": mean_rho,
        "median_rho": float(np.median(rho)),
        "max_rho": float(np.max(rho)),
        "rho_q95": float(np.percentile(rho, 95)),
        "residual_ambient_score": residual_score,
        "top_corrected_genes": top_corrected_genes,
        "review_required": review_required,
        "risk_note": risk_note,
    }

    if record:
        adata.uns.setdefault("sclucid", {}).setdefault("qc", {})[key_added] = sanitize_for_hdf5(
            result
        )
        record_ambient_correction_status(
            adata,
            corrected=True,
            backend="linear_background_subtraction",
            output_layer=output_layer,
            details=result,
        )
        record_ambient_layer_contract(
            adata,
            input_context=infer_ambient_input_context(adata, layer=layer),
            correction_summary=result,
            output_layer=output_layer,
        )
    return result
