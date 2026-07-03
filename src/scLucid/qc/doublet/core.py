"""Doublet detection core — configuration helpers and simple public API.

Extracted for maintainability.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData

from ...utils import get_marker_manager
from ..config import DoubletConfig, MarkerConfig

log = logging.getLogger(__name__)

# Constants for column naming consistency across submodules
LINEAGE_SCORES_KEY = "lineage_module_scores"
HEURISTIC_SCORE_COL = "heuristic_confidence_score"
HEURISTIC_PRED_COL = "heuristic_predicted"
FINAL_PRED_COL = "predicted_doublet"
ALGORITHM_SCORE_COL = "algorithm_doublet_score"
ALGORITHM_PRED_COL = "algorithm_predicted_doublet"
COMBINED_SCORE_COL = "combined_doublet_score"
HETEROTYPIC_RISK_COL = "heterotypic_doublet_risk"
HOMOTYPIC_RISK_COL = "homotypic_doublet_risk"
EXPECTED_TOTAL_RATE_COL = "expected_total_doublet_rate"
EXPECTED_HETEROTYPIC_RATE_COL = "expected_heterotypic_doublet_rate"
EXPECTED_HOMOTYPIC_RATE_COL = "expected_homotypic_doublet_rate"

DOUBLET_OBS_COLUMNS = {
    "final": [FINAL_PRED_COL, COMBINED_SCORE_COL],
    "algorithm": [ALGORITHM_SCORE_COL, ALGORITHM_PRED_COL],
    "heuristic": [HEURISTIC_SCORE_COL, HEURISTIC_PRED_COL],
    "risk_decomposition": [
        HETEROTYPIC_RISK_COL,
        HOMOTYPIC_RISK_COL,
        EXPECTED_TOTAL_RATE_COL,
        EXPECTED_HETEROTYPIC_RATE_COL,
        EXPECTED_HOMOTYPIC_RATE_COL,
    ],
    "intermediate_patterns": [
        "scrublet_",
        "scanpy_scrublet_",
        "solo_",
        "doubletdetection_",
        "scdblfinder_",
    ],
}

__all__ = [
    "LINEAGE_SCORES_KEY",
    "HEURISTIC_SCORE_COL",
    "HEURISTIC_PRED_COL",
    "FINAL_PRED_COL",
    "ALGORITHM_SCORE_COL",
    "ALGORITHM_PRED_COL",
    "COMBINED_SCORE_COL",
    "HETEROTYPIC_RISK_COL",
    "HOMOTYPIC_RISK_COL",
    "EXPECTED_TOTAL_RATE_COL",
    "EXPECTED_HETEROTYPIC_RATE_COL",
    "EXPECTED_HOMOTYPIC_RATE_COL",
    "DOUBLET_OBS_COLUMNS",
    "generate_doublet_rates",
    "create_custom_marker_dict",
    "audit_doublets",
]


def _coerce_expected_rate(rate: Optional[float], *, default: float = 0.1) -> float:
    """Clamp expected-rate inputs to a stable [0, 1] float."""
    if rate is None:
        return float(default)
    return max(0.0, min(1.0, float(rate)))


def _expected_rate_series(
    adata: AnnData,
    expected_rate: Optional[Union[float, Dict[str, float]]],
    *,
    sample_key: str,
    default_rate: float = 0.1,
) -> pd.Series:
    """Return a per-cell expected-rate series with stable fallback behavior."""
    if isinstance(expected_rate, dict) and sample_key in adata.obs.columns:
        mapped = adata.obs[sample_key].map(expected_rate)
        fallback = (
            float(np.mean(list(expected_rate.values()))) if expected_rate else float(default_rate)
        )
        return pd.to_numeric(mapped, errors="coerce").fillna(fallback).astype(float)
    if expected_rate is None:
        return pd.Series(float(default_rate), index=adata.obs_names, dtype=float)
    return pd.Series(float(expected_rate), index=adata.obs_names, dtype=float)


def _expected_rate_topk_predictions(
    scores,
    expected_rate: float,
    *,
    eligible_mask=None,
    universe_size: Optional[int] = None,
) -> tuple[np.ndarray, float]:
    """Flag the top expected-rate eligible cells by score, robust to ties."""
    scores = np.asarray(scores, dtype=float).ravel()
    predicted = np.zeros(scores.shape[0], dtype=bool)
    if scores.size == 0:
        return predicted, float("nan")

    eligible = np.isfinite(scores)
    if eligible_mask is not None:
        eligible &= np.asarray(eligible_mask, dtype=bool).ravel()
    if not np.any(eligible):
        finite = np.isfinite(scores)
        fallback_threshold = float(np.nanmax(scores[finite])) if np.any(finite) else float("nan")
        return predicted, fallback_threshold

    rate = _coerce_expected_rate(expected_rate)
    base_size = int(universe_size) if universe_size is not None else scores.size
    n_expected = int(np.ceil(rate * max(base_size, 0)))
    if rate > 0 and n_expected == 0:
        n_expected = 1
    n_expected = min(n_expected, int(np.sum(eligible)))
    if n_expected <= 0:
        threshold = float(np.nanmax(scores[eligible]))
        return predicted, threshold

    eligible_indices = np.flatnonzero(eligible)
    eligible_scores = scores[eligible]
    if np.nanmax(eligible_scores) <= 0 and (
        np.nanmax(eligible_scores) - np.nanmin(eligible_scores) <= 1e-12
    ):
        return predicted, float(np.nanmax(eligible_scores))
    order = np.argsort(eligible_scores, kind="mergesort")
    selected = eligible_indices[order[-n_expected:]]
    predicted[selected] = True
    threshold = float(np.min(scores[selected]))
    return predicted, threshold


def _expected_rate_grouped_predictions(
    scores: pd.Series,
    *,
    expected_rate: Optional[Union[float, Dict[str, float]]],
    groups: Optional[pd.Series] = None,
    eligible_mask: Optional[pd.Series] = None,
    default_rate: float = 0.1,
) -> tuple[pd.Series, Union[float, Dict[str, float]]]:
    """Apply expected-rate thresholding globally or per group with shared semantics."""
    scores = pd.to_numeric(scores, errors="coerce").astype(float)
    if eligible_mask is not None:
        eligible_mask = pd.Series(eligible_mask, index=scores.index, dtype=bool)

    if isinstance(expected_rate, dict) and groups is not None:
        groups = pd.Series(groups, index=scores.index)
        predicted = pd.Series(False, index=scores.index, dtype=bool)
        thresholds: Dict[str, float] = {}
        fallback_rate = (
            float(np.mean(list(expected_rate.values()))) if expected_rate else float(default_rate)
        )
        for group_name, idx in groups.groupby(groups, observed=False).groups.items():
            rate = _coerce_expected_rate(expected_rate.get(group_name), default=fallback_rate)
            sample_scores = scores.loc[idx]
            sample_eligible = eligible_mask.loc[idx] if eligible_mask is not None else None
            sample_pred, threshold = _expected_rate_topk_predictions(
                sample_scores.to_numpy(),
                rate,
                eligible_mask=sample_eligible.to_numpy() if sample_eligible is not None else None,
                universe_size=sample_scores.shape[0],
            )
            predicted.loc[idx] = sample_pred
            thresholds[str(group_name)] = threshold
        return predicted, thresholds

    rate = (
        _coerce_expected_rate(float(expected_rate), default=default_rate)
        if expected_rate is not None and not isinstance(expected_rate, dict)
        else float(default_rate)
    )
    predicted, threshold = _expected_rate_topk_predictions(
        scores.to_numpy(),
        rate,
        eligible_mask=eligible_mask.to_numpy() if eligible_mask is not None else None,
        universe_size=scores.shape[0],
    )
    return pd.Series(predicted, index=scores.index, dtype=bool), threshold


def audit_doublets(
    adata: AnnData,
    prediction_cols: Optional[List[str]] = None,
    score_cols: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Audit doublet predictions and scores remaining in ``adata``.

    Useful after filtering to confirm that predicted doublets were removed, or to
    inspect residual doublet signal.

    Parameters
    ----------
    adata
        AnnData object (typically after QC filtering).
    prediction_cols
        Boolean ``.obs`` columns to summarize. Defaults to common scLucid doublet
        prediction columns.
    score_cols
        Numeric ``.obs`` columns to summarize. Defaults to common scLucid doublet
        score columns.

    Returns:
    -------
    Dict[str, Any]
        Summary with ``predictions`` and ``scores`` sub-dicts.

    Examples:
    --------
    >>> summary = scl.qc.audit_doublets(adata_filtered)
    >>> print(summary["predictions"]["predicted_doublet"])
    """
    if prediction_cols is None:
        prediction_cols = [
            "predicted_doublet",
            "scrublet_predicted",
            "doubletdetection_predicted",
            "heuristic_predicted",
        ]
    if score_cols is None:
        score_cols = [
            "algorithm_doublet_score",
            "doublet_score",
            "heterotypic_doublet_risk",
            "homotypic_doublet_risk",
        ]

    predictions = {}
    for col in prediction_cols:
        if col not in adata.obs.columns:
            continue
        n = int(adata.obs[col].sum()) if adata.obs[col].dtype == bool else None
        pct = (n / adata.n_obs * 100) if n is not None else None
        predictions[col] = {"count": n, "percent": pct}
        if n is not None:
            print(f"  {col}: {n} cells ({pct:.2f}%)")
            if col == "predicted_doublet" and n > 0:
                print(
                    f"    NOTE: {n} predicted doublets remain after filtering — "
                    "review FilterConfig criteria."
                )

    scores = {}
    for col in score_cols:
        if col not in adata.obs.columns:
            continue
        vals = pd.to_numeric(adata.obs[col], errors="coerce")
        scores[col] = {
            "median": float(vals.median()),
            "p95": float(vals.quantile(0.95)),
            "max": float(vals.max()),
        }
        print(
            f"  {col}: median={vals.median():.3f}, "
            f"p95={vals.quantile(0.95):.3f}, max={vals.max():.3f}"
        )

    return {"predictions": predictions, "scores": scores}


def _create_doublet_marker_config_from_manager(
    adata: AnnData, cfg: DoubletConfig
) -> Dict[str, MarkerConfig]:
    """
    Create MarkerConfig objects using the marker manager and the main DoubletConfig.

    Args:
        adata: The AnnData object to intersect markers with.
        cfg: The main DoubletConfig object, providing context like species, tissue,
             and default evaluation parameters.

    Returns:
        A dictionary mapping lineage names to MarkerConfig objects.
    """
    try:
        # --- ❗ Select the correct adata object upfront ❗ ---
        # Ensure the intersection is performed on the same data that will be used for scoring.
        adata_for_intersection = (
            adata.raw.to_adata() if cfg.default_use_raw and adata.raw else adata
        )
        log.info(
            f"Performing marker intersection on {'adata.raw' if cfg.default_use_raw and adata.raw else 'adata'}."
        )

        case_sensitive = cfg.marker_species.lower() == "mouse"
        manager = get_marker_manager(
            species=cfg.marker_species,
            tissue=cfg.marker_tissue,
            case_sensitive=case_sensitive,
        )
        # Intersect with the correctly chosen data object
        manager.intersect_with(adata_for_intersection)
        markers_dict = manager.get_doublet_lineage_markers()

        if not markers_dict:
            log.warning(
                "No lineage markers found. Ensure `doublet_lineage = true` is set in your TOML file for desired cell types."
            )
            return {}

    except Exception as e:
        log.error(
            f"Failed to load markers via marker_manager: {e}. Cannot create heuristic configs."
        )
        return {}

    marker_configs = {}
    for lineage, genes in markers_dict.items():
        if genes:
            marker_configs[lineage] = MarkerConfig(
                genes=genes,
                expression_threshold=cfg.default_expression_threshold,
                min_genes_required=cfg.default_min_genes_required,
                use_raw=cfg.default_use_raw,
            )
    log.info(f"Auto-generated {len(marker_configs)} marker configurations for doublet detection.")
    return marker_configs


def generate_doublet_rates(
    adata: AnnData,
    sample_key: str = "sampleID",
    chemistry: str = "v3",
    custom_rate: float = 0.008,
    custom_rate_model: Literal["scale", "fixed"] = "scale",
    max_rate: float = 0.20,
    min_rate: float = 0.001,
) -> Dict[str, float]:
    """
    Automatically generate expected doublet rates based on cell count per sample or platform.

    This function supports three modes:
    1.  Platform Models ('v2', 'v3', 'HT'): Applies 10x-style linear/non-linear scaling
        based on cell count.
    2.  Fixed Models ('BD'): Applies a single, fixed rate to all samples.
    3.  Custom Models ('custom'):
        - If `custom_rate_model='scale'`, applies 10x-style scaling using `custom_rate`
          as the factor per 1000 cells.
        - If `custom_rate_model='fixed'`, applies `custom_rate` as a single fixed rate
          to all samples (ignoring cell counts).

    Args:
        adata: AnnData object containing cell count information.
        sample_key: Column name in adata.obs used to distinguish samples.
        chemistry: The technology platform ('v2', 'v3', 'HT', 'BD', or 'custom').
        custom_rate: The rate to use when `chemistry='custom'`. Interpreted based on
                     `custom_rate_model`. Defaults to 0.008.
        custom_rate_model: Defines how to use `custom_rate` ('scale' or 'fixed').
                           Defaults to 'scale'.
        max_rate: (For 'scale' models) Maximum doublet rate cap.
        min_rate: (For 'scale' models) Minimum doublet rate floor.

    Returns:
        Dictionary mapping sample IDs to calculated doublet rates.

    Notes:
        The returned rates can substantially affect final labels when
        ``DoubletConfig.final_label_strategy='expected_rate_rank'`` because the
        merged evidence score is thresholded by expected rate within detection
        groups. Use ``final_label_strategy='algorithm_label'`` when the final
        call should follow the algorithm's binary prediction directly.
    """
    log.info("Automatically generating doublet rates based on sample chemistry and cell counts...")

    # Define platform-specific models: (model_type, rate_value)
    # 'scale' model_type uses the rate_value as a scaling factor per 1000 cells.
    # 'fixed' model_type uses the rate_value as the final, fixed rate for all samples.
    chemistry_models = {
        "v2": ("scale", 0.007),  # 10x v2 chemistry
        "v3": ("scale", 0.008),  # 10x v3 chemistry
        "HT": ("scale", 0.016),  # 10x High-throughput
        "BD": ("fixed", 0.025),  # BD Rhapsody
    }

    cell_counts = adata.obs[sample_key].value_counts()
    doublet_rates = {}

    model_type = None
    rate_value = None

    # --- 1. Determine Model Type and Rate ---
    if chemistry in chemistry_models:
        model_type, rate_value = chemistry_models[chemistry]
        log.info(
            f"Using known platform model '{chemistry}': type='{model_type}', base_rate={rate_value}"
        )
    elif chemistry == "custom":
        model_type = custom_rate_model
        rate_value = custom_rate
        log.info(f"Using 'custom' model: type='{model_type}', custom_rate={rate_value}")
    else:
        log.warning(f"Unknown chemistry '{chemistry}'. Falling back to default 'v3' model.")
        model_type, rate_value = chemistry_models["v3"]
        chemistry = "v3"  # Set chemistry for scaling logic

    # --- 2. Apply Model to Calculate Rates ---
    if model_type == "fixed":
        log.info(
            f"Applying fixed doublet rate of {rate_value:.4f} to all {len(cell_counts)} samples."
        )
        for sample, n_cells in cell_counts.items():
            doublet_rates[sample] = rate_value
            log.info(f"  - Sample '{sample}': {n_cells} cells -> Doublet rate: {rate_value:.4f}")

    elif model_type == "scale":
        log.info(f"Applying scaling model with base rate of {rate_value:.4f} per 1000 cells.")
        for sample, n_cells in cell_counts.items():
            # Standard linear scaling (10x Genomics documented baseline).
            # References:
            #   - 10x Genomics Support: ~0.8% per 1k cells (v2/v3), ~1.6% (HT).
            #   - McGinnis et al. 2019; Xi & Li 2021.
            n_k = n_cells / 1000.0
            rate = n_k * rate_value

            # Poisson saturation correction for high loading.
            # Under Poisson loading the multiplet fraction is
            #   p_multi = 1 - exp(-λ) - λ·exp(-λ)
            #   p_non_empty = 1 - exp(-λ)
            # For small λ the linear approximation rate ≈ λ/2 holds.
            # At high loading the exact Poisson formula gives a modest
            # upward correction (consistent with 10x observed curves).
            if n_cells > 5000:
                lambda_ = 2.0 * rate  # from rate ≈ λ/2
                if lambda_ < 20.0:
                    p_nonempty = 1.0 - math.exp(-lambda_)
                    if p_nonempty > 1e-12:
                        p_multiplet = 1.0 - math.exp(-lambda_) - lambda_ * math.exp(-lambda_)
                        rate = p_multiplet / p_nonempty

            # Apply rate constraints
            rate = max(min_rate, min(rate, max_rate))
            doublet_rates[sample] = rate
            log.info(f"  - Sample '{sample}': {n_cells} cells -> Doublet rate: {rate:.4f}")

    else:
        # This case should not be reachable if logic is sound
        raise ValueError(f"Internal error: Unrecognized model_type '{model_type}'")

    return doublet_rates


def create_custom_marker_dict(
    lineage_definitions: Dict[str, Dict], save_path: Optional[Union[str, Path]] = None
) -> Dict[str, MarkerConfig]:
    """
    Create custom marker dictionary from user-defined lineage specifications.

    This function allows users to define their own marker sets with custom
    parameters for specialized doublet detection scenarios.

    Args:
        lineage_definitions: Dictionary defining lineages and their parameters
        save_path: Optional path to save the configuration for future use

    Returns:
        Dictionary mapping lineage names to MarkerConfig objects

    Example:
        lineage_defs = {
            "T_cells": {
                "genes": ["CD3D", "CD3E", "CD8A"],
                "expression_threshold": 0.5,
                "min_genes_required": 1
            },
            "Epithelial": {
                "genes": r"^KRT[0-9]+",  # Regex pattern
                "expression_threshold": 1.0,
                "min_genes_required": 2
            }
        }
        marker_configs = create_custom_marker_dict(lineage_defs)
    """
    config_dict = {}

    for lineage, definition in lineage_definitions.items():
        # Validate required 'genes' field
        if "genes" not in definition:
            raise ValueError(f"Missing 'genes' field for lineage '{lineage}'")

        config_dict[lineage] = MarkerConfig(**definition)

    # Save configuration if requested
    if save_path:
        import json

        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)

        # Convert to serializable format
        serializable_dict = {}
        for lineage, config in config_dict.items():
            serializable_dict[lineage] = {
                "genes": config.genes,
                "expression_threshold": config.expression_threshold,
                "min_genes_required": config.min_genes_required,
                "use_raw": config.use_raw,
            }

        with open(save_path, "w") as f:
            json.dump(serializable_dict, f, indent=2)
        log.info(f"Marker configuration saved to {save_path}")

    return config_dict
