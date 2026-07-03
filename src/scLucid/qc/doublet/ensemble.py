"""Ensemble doublet detection pipeline and evidence profiling.

Extracted from core.py for maintainability.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData

from ..artifacts import record_qc_artifact_contract
from ..config import DoubletConfig
from .algorithms import (
    _run_doubletdetection,
    _run_scanpy_scrublet,
    _run_scdblfinder,
    _run_scrublet,
    _run_solo,
)
from .core import (
    ALGORITHM_PRED_COL,
    ALGORITHM_SCORE_COL,
    COMBINED_SCORE_COL,
    EXPECTED_HETEROTYPIC_RATE_COL,
    EXPECTED_HOMOTYPIC_RATE_COL,
    EXPECTED_TOTAL_RATE_COL,
    FINAL_PRED_COL,
    HETEROTYPIC_RISK_COL,
    HEURISTIC_PRED_COL,
    HEURISTIC_SCORE_COL,
    HOMOTYPIC_RISK_COL,
    LINEAGE_SCORES_KEY,
    _expected_rate_grouped_predictions,
    _expected_rate_series,
)
from .heuristic import _plot_doublet_summary, _run_heuristic

log = logging.getLogger(__name__)


__all__: list[str] = ["predict_doublets"]


# --- Helper Functions ---
def _score_distribution_stats(
    scores: pd.Series, groups: Optional[pd.Series] = None
) -> Dict[str, Any]:
    """Compute per-group score distribution statistics for audit purposes."""
    scores = pd.to_numeric(scores, errors="coerce").fillna(0.0).astype(float)
    if groups is None:
        return {
            "global": {
                "min": float(scores.min()),
                "max": float(scores.max()),
                "median": float(scores.median()),
                "p95": float(scores.quantile(0.95)),
                "n": int(scores.shape[0]),
            }
        }
    groups = pd.Series(groups, index=scores.index)
    stats: Dict[str, Any] = {}
    for group_name, idx in groups.groupby(groups, observed=False).groups.items():
        vals = scores.loc[idx]
        stats[str(group_name)] = {
            "min": float(vals.min()),
            "max": float(vals.max()),
            "median": float(vals.median()),
            "p95": float(vals.quantile(0.95)),
            "n": int(vals.shape[0]),
        }
    return stats


def _rank_normalize_within_group(
    values: pd.Series, groups: Optional[pd.Series] = None
) -> pd.Series:
    """Convert scores to percentile ranks within each detection group.

    Constant or single-value groups return zeros to avoid artificial signal.
    """
    values = pd.to_numeric(values, errors="coerce").fillna(0.0).astype(float)
    if groups is None:
        if values.max() <= values.min() or values.nunique() <= 1:
            return pd.Series(0.0, index=values.index)
        return values.rank(method="min", pct=True).fillna(0.0)

    groups = pd.Series(groups, index=values.index)
    normalized = pd.Series(0.0, index=values.index, dtype=float)
    for _, idx in groups.groupby(groups, observed=False).groups.items():
        group_vals = values.loc[idx]
        if group_vals.max() <= group_vals.min() or group_vals.nunique() <= 1:
            normalized.loc[idx] = 0.0
        else:
            normalized.loc[idx] = group_vals.rank(method="min", pct=True).fillna(0.0)
    return normalized


def _zscore_within_group(values: pd.Series, groups: Optional[pd.Series] = None) -> pd.Series:
    """Z-score normalize values within each detection group, clipping negatives to 0."""
    values = pd.to_numeric(values, errors="coerce").fillna(0.0).astype(float)
    if groups is None:
        std = values.std()
        if pd.isna(std) or std == 0:
            return pd.Series(0.0, index=values.index)
        return ((values - values.mean()) / std).clip(lower=0.0)

    groups = pd.Series(groups, index=values.index)
    normalized = pd.Series(0.0, index=values.index, dtype=float)
    for _, idx in groups.groupby(groups, observed=False).groups.items():
        group_vals = values.loc[idx]
        std = group_vals.std()
        if pd.isna(std) or std == 0:
            normalized.loc[idx] = 0.0
        else:
            normalized.loc[idx] = ((group_vals - group_vals.mean()) / std).clip(lower=0.0)
    return normalized


def _normalize_score(
    values: pd.Series,
    method: str,
    groups: Optional[pd.Series] = None,
) -> pd.Series:
    """Normalize an evidence score using the configured method."""
    values = pd.to_numeric(values, errors="coerce").fillna(0.0).astype(float)
    if method == "raw":
        return values
    if method == "rank":
        return _rank_normalize_within_group(values, groups=groups)
    if method == "zscore":
        return _zscore_within_group(values, groups=groups)

    # method == "minmax" (default fallback)
    if groups is None:
        min_v = values.min()
        max_v = values.max()
        if pd.isna(max_v) or max_v <= min_v:
            return pd.Series(0.0, index=values.index)
        return (values - min_v) / (max_v - min_v)

    groups = pd.Series(groups, index=values.index)
    normalized = pd.Series(0.0, index=values.index, dtype=float)
    for _, idx in groups.groupby(groups, observed=False).groups.items():
        normalized.loc[idx] = _normalize_score(values.loc[idx], method="minmax", groups=None)
    return normalized


def _compute_merged_doublet_score(
    adata: AnnData,
    algorithm_score_col: str,
    heuristic_score_col: str,
    strategy: str = "weighted_average",
    algo_weight: float = 0.6,
    algorithm_score_normalization: str = "rank",
    heuristic_score_normalization: str = "rank",
    merge_rank_normalize: bool = True,
    groups: Optional[pd.Series] = None,
) -> pd.Series:
    """Compute the merged doublet evidence score from normalized inputs."""
    algo_scores = _normalize_score(
        adata.obs[algorithm_score_col], method=algorithm_score_normalization, groups=groups
    )
    heur_scores = _normalize_score(
        adata.obs[heuristic_score_col], method=heuristic_score_normalization, groups=groups
    )

    if strategy == "weighted_average":
        final_score = (algo_weight * algo_scores) + ((1 - algo_weight) * heur_scores)
    elif strategy == "max_score":
        final_score = pd.DataFrame({"algo": algo_scores, "heur": heur_scores}).max(axis=1)
    elif strategy == "heuristic_boost":
        final_score = algo_scores + (heur_scores * 0.5)
    elif strategy == "calibrated_threshold":
        # For calibrated threshold, use the rank-normalized merged score directly.
        final_score = (algo_weight * algo_scores) + ((1 - algo_weight) * heur_scores)
    else:
        log.warning(
            f"Unknown enhanced merge strategy '{strategy}', falling back to 'weighted_average'."
        )
        final_score = (algo_weight * algo_scores) + ((1 - algo_weight) * heur_scores)

    if merge_rank_normalize:
        return _rank_normalize_within_group(final_score, groups=groups)
    return final_score


def _merge_doublet_predictions(
    adata: AnnData,
    algorithm_score_col: str,
    heuristic_score_col: str,
    strategy: str = "weighted_average",
    algo_weight: float = 0.6,
    algorithm_score_normalization: str = "rank",
    heuristic_score_normalization: str = "rank",
    merge_rank_normalize: bool = True,
    expected_rate: Optional[Union[float, Dict[str, float]]] = 0.1,
    score_threshold: Optional[float] = None,
    sample_key: Optional[str] = None,
    final_label_strategy: str = "expected_rate_rank",
    algorithm_pred_col: Optional[str] = None,
) -> pd.Series:
    """
    Merge algorithmic and heuristic doublet scores for a final prediction.

    Both inputs are normalized within detection groups according to the configured
    methods, then combined. The final score is optionally rank-normalized again
    before thresholding. This keeps scores comparable across samples.
    """
    calibration_groups = (
        adata.obs[sample_key]
        if sample_key is not None and sample_key in adata.obs.columns
        else None
    )
    final_score = _compute_merged_doublet_score(
        adata,
        algorithm_score_col=algorithm_score_col,
        heuristic_score_col=heuristic_score_col,
        strategy=strategy,
        algo_weight=algo_weight,
        algorithm_score_normalization=algorithm_score_normalization,
        heuristic_score_normalization=heuristic_score_normalization,
        merge_rank_normalize=merge_rank_normalize,
        groups=calibration_groups,
    )
    adata.obs[COMBINED_SCORE_COL] = final_score

    doublet_params = (
        adata.uns.setdefault("sclucid", {}).setdefault("qc", {}).setdefault("doublet_params", {})
    )
    doublet_params["score_calibration"] = {
        "schema_version": "doublet_score_calibration_v2",
        "calibration_group_key": sample_key if calibration_groups is not None else None,
        "algorithm_score_input": algorithm_score_col,
        "heuristic_score_input": heuristic_score_col,
        "combined_score": COMBINED_SCORE_COL,
        "algorithm_score_normalization": algorithm_score_normalization,
        "heuristic_score_normalization": heuristic_score_normalization,
        "merge_rank_normalize": merge_rank_normalize,
        "merge_strategy": strategy,
        "final_label_strategy": final_label_strategy,
        "score_semantics": "Scores are rank-like evidence indices, not calibrated probabilities.",
        "algorithm_score_distribution": _score_distribution_stats(
            adata.obs[algorithm_score_col], groups=calibration_groups
        ),
        "heuristic_score_distribution": _score_distribution_stats(
            adata.obs[heuristic_score_col], groups=calibration_groups
        ),
        "merged_score_distribution": _score_distribution_stats(
            final_score, groups=calibration_groups
        ),
    }

    if final_label_strategy == "algorithm_label":
        if algorithm_pred_col is None or algorithm_pred_col not in adata.obs:
            raise ValueError(
                "algorithm_pred_col must exist when final_label_strategy='algorithm_label'"
            )
        doublet_params["score_calibration"]["thresholding"] = {
            "method": "algorithm_label",
            "note": (
                "predicted_doublet follows the algorithm binary call directly. "
                "Merged and heuristic scores are retained as review evidence only."
            ),
        }
        return adata.obs[algorithm_pred_col].fillna(False).astype(bool)

    if score_threshold is not None:
        threshold = score_threshold
        log.info(
            f"Using user-provided doublet score threshold of {threshold:.3f} for merged predictions."
        )
        doublet_params["score_calibration"]["thresholding"] = {
            "method": "score_threshold",
            "score_threshold": float(threshold),
        }
        return final_score > threshold

    if expected_rate is None:
        log.warning("expected_doublet_rate is None, using a default of 0.1 for thresholding.")
        expected_rate = 0.1

    grouped_scores = (
        adata.obs[sample_key]
        if isinstance(expected_rate, dict) and sample_key in adata.obs
        else None
    )
    merged_pred, threshold = _expected_rate_grouped_predictions(
        final_score,
        expected_rate=expected_rate,
        groups=grouped_scores,
    )
    if isinstance(threshold, dict):
        fallback = float(np.mean(list(expected_rate.values()))) if expected_rate else 0.1
        for sample, sample_threshold in threshold.items():
            sample_rate = float(expected_rate.get(sample, fallback))
            log.info(
                "Using final score threshold %.3f for sample '%s' based on expected doublet rate %.3f.",
                sample_threshold,
                sample,
                sample_rate,
            )
        doublet_params["score_calibration"]["thresholding"] = {
            "method": "expected_rate_rank",
            "expected_doublet_rate": expected_rate,
            "thresholds": threshold,
            "note": (
                "Final labels are assigned to the top-ranked merged evidence scores "
                "within detection groups according to expected_doublet_rate."
            ),
        }
        return merged_pred

    log.info(
        f"Using a final score threshold of {threshold:.3f} based on expected doublet rate for merged predictions."
    )
    doublet_params["score_calibration"]["thresholding"] = {
        "method": "expected_rate_rank",
        "expected_doublet_rate": expected_rate,
        "threshold": float(threshold),
        "note": (
            "Final labels are assigned to the top-ranked merged evidence scores "
            "according to expected_doublet_rate."
        ),
    }

    return merged_pred


def _collect_external_doublet_evidence(
    adata: AnnData,
    *,
    external_cols: list[str],
    final_col: str,
    policy: str,
) -> dict[str, object]:
    """Record optional hashing/genotype/manual doublet evidence."""
    present_cols = [col for col in external_cols if col in adata.obs]
    missing_cols = [col for col in external_cols if col not in adata.obs]
    if not present_cols and not missing_cols:
        return {
            "available": False,
            "columns_used": [],
            "missing_columns": [],
            "policy": policy,
            "n_external_doublets": 0,
        }

    external_mask = pd.Series(False, index=adata.obs_names, dtype=bool)
    per_column: dict[str, dict[str, object]] = {}
    for col in present_cols:
        values = adata.obs[col]
        if pd.api.types.is_bool_dtype(values):
            mask = values.fillna(False).astype(bool)
        else:
            text = values.astype(str).str.lower()
            mask = text.isin({"true", "1", "doublet", "multiplet", "positive", "yes"})
        external_mask |= mask
        per_column[col] = {
            "n_positive": int(mask.sum()),
            "fraction_positive": float(mask.mean()) if adata.n_obs else 0.0,
        }

    adata.obs["external_doublet_evidence"] = external_mask
    if policy == "include_in_final" and present_cols:
        adata.obs[final_col] = adata.obs[final_col].astype(bool) | external_mask

    return {
        "available": bool(present_cols),
        "columns_used": present_cols,
        "missing_columns": missing_cols,
        "policy": policy,
        "n_external_doublets": int(external_mask.sum()),
        "fraction_external_doublets": float(external_mask.mean()) if adata.n_obs else 0.0,
        "included_in_final": bool(policy == "include_in_final" and present_cols),
        "per_column": per_column,
    }


def _lineage_mixture_fraction(adata: AnnData) -> tuple[pd.Series, Dict[str, float]]:
    """Estimate heterotypic opportunity from dominant lineage composition."""
    if LINEAGE_SCORES_KEY not in adata.obsm:
        return pd.Series(0.5, index=adata.obs_names, dtype=float), {
            "available": False,
            "global_heterotypic_fraction": 0.5,
            "global_homotypic_fraction": 0.5,
        }
    lineage_scores = adata.obsm[LINEAGE_SCORES_KEY]
    if lineage_scores.empty:
        return pd.Series(0.5, index=adata.obs_names, dtype=float), {
            "available": False,
            "global_heterotypic_fraction": 0.5,
            "global_homotypic_fraction": 0.5,
        }
    dominant = lineage_scores.idxmax(axis=1)
    max_scores = lineage_scores.max(axis=1)
    dominant = dominant.where(max_scores > 0, "Unknown")
    freqs = dominant.value_counts(normalize=True)
    homotypic_fraction = float((freqs**2).sum())
    heterotypic_fraction = float(max(0.0, 1.0 - homotypic_fraction))
    per_cell_heterotypic = dominant.map(lambda lin: max(0.0, 1.0 - float(freqs.get(lin, 0.0))))
    return per_cell_heterotypic.reindex(adata.obs_names).fillna(heterotypic_fraction), {
        "available": True,
        "global_heterotypic_fraction": heterotypic_fraction,
        "global_homotypic_fraction": homotypic_fraction,
        "dominant_lineage_frequencies": {str(k): float(v) for k, v in freqs.items()},
    }


def _add_doublet_risk_decomposition(
    adata: AnnData,
    *,
    algorithm_score_col: str,
    algorithm_pred_col: str,
    heuristic_score_col: str,
    expected_rate: Optional[Union[float, Dict[str, float]]],
    sample_key: str,
) -> Dict[str, object]:
    """Add heterotypic/homotypic risk columns and expected-rate decomposition."""
    calibration_groups = adata.obs[sample_key] if sample_key in adata.obs.columns else None
    algo_score = _normalize_score(
        adata.obs[algorithm_score_col], method="rank", groups=calibration_groups
    )
    heur_score = _normalize_score(
        adata.obs[heuristic_score_col], method="rank", groups=calibration_groups
    )
    combined_score = (
        _normalize_score(adata.obs[COMBINED_SCORE_COL], method="rank", groups=calibration_groups)
        if COMBINED_SCORE_COL in adata.obs
        else _normalize_score(algo_score + heur_score, method="rank", groups=calibration_groups)
    )

    heterotypic_opportunity, lineage_meta = _lineage_mixture_fraction(adata)
    total_expected = _expected_rate_series(adata, expected_rate, sample_key=sample_key)
    expected_heterotypic = total_expected * heterotypic_opportunity
    expected_homotypic = (total_expected - expected_heterotypic).clip(lower=0.0)

    complexity_terms = []
    for col in ["n_genes_by_counts", "total_counts"]:
        if col in adata.obs.columns:
            z = pd.to_numeric(adata.obs[col], errors="coerce")
            std = z.std()
            if pd.notna(std) and std > 0:
                complexity_terms.append(((z - z.mean()) / std).clip(lower=0.0))
    if complexity_terms:
        complexity = pd.concat(complexity_terms, axis=1).mean(axis=1)
        complexity = _normalize_score(complexity, method="rank", groups=calibration_groups)
    else:
        complexity = pd.Series(0.0, index=adata.obs_names)

    heterotypic_risk = _normalize_score(
        0.75 * heur_score + 0.25 * combined_score,
        method="rank",
        groups=calibration_groups,
    )
    homotypic_risk = _normalize_score(
        0.55 * algo_score + 0.35 * complexity + 0.10 * combined_score,
        method="rank",
        groups=calibration_groups,
    )
    # Heterotypic co-expression explains the risk better than homotypic complexity.
    homotypic_risk = (homotypic_risk * (1.0 - 0.5 * heterotypic_risk)).clip(0.0, 1.0)

    adata.obs[ALGORITHM_SCORE_COL] = algo_score
    adata.obs[ALGORITHM_PRED_COL] = adata.obs[algorithm_pred_col].astype(bool)
    adata.obs[HETEROTYPIC_RISK_COL] = heterotypic_risk
    adata.obs[HOMOTYPIC_RISK_COL] = homotypic_risk
    adata.obs[EXPECTED_TOTAL_RATE_COL] = total_expected
    adata.obs[EXPECTED_HETEROTYPIC_RATE_COL] = expected_heterotypic
    adata.obs[EXPECTED_HOMOTYPIC_RATE_COL] = expected_homotypic

    return {
        "schema_version": "doublet_risk_decomposition_v2",
        "score_semantics": {
            ALGORITHM_SCORE_COL: "Normalized algorithm score; not calibrated probability.",
            HEURISTIC_SCORE_COL: "Marker co-expression evidence score; not calibrated probability.",
            COMBINED_SCORE_COL: "Merge score used for thresholding; not calibrated probability.",
            HETEROTYPIC_RISK_COL: "Evidence for cross-lineage doublets.",
            HOMOTYPIC_RISK_COL: "Evidence for same-lineage doublets using algorithm score and transcript complexity.",
        },
        "score_calibration": {
            "calibration_group_key": sample_key if calibration_groups is not None else None,
            "normalization": "rank_within_detection_group",
        },
        "expected_rates": {
            "mean_total": float(total_expected.mean()),
            "mean_heterotypic": float(expected_heterotypic.mean()),
            "mean_homotypic": float(expected_homotypic.mean()),
        },
        "lineage_mixture": lineage_meta,
        "evidence_priority": {
            "algorithm_positive_final_priority": (
                "Allowlisted co-expression can downgrade heuristic-only calls, but algorithm-positive "
                "cells remain final positives unless the merge score threshold excludes them."
            ),
            "heterotypic_sources": ["lineage co-expression", "heuristic_confidence_score"],
            "homotypic_sources": [
                "algorithm score",
                "high gene/UMI complexity within similar lineage",
            ],
        },
    }


def _export_doublet_stats(
    adata: AnnData,
    sample_key: str = "sampleID",
    save_dir: Optional[Union[str, Path]] = None,
    export_csv: bool = True,
    export_xlsx: bool = False,
) -> Dict[str, pd.DataFrame]:
    """
    Export comprehensive doublet statistics per sample and globally.

    This function generates detailed statistical summaries of doublet detection
    results, including counts, percentages, and score distributions.

    Args:
        adata: AnnData object with doublet predictions
        sample_key: Key for sample identification
        save_dir: Directory to save statistics files
        export_csv: Whether to export as CSV files
        export_xlsx: Whether to export as Excel file

    Returns:
        Dictionary containing sample-wise and global statistics DataFrames
    """
    # Identify all doublet-related columns
    doublet_cols = [
        col
        for col in adata.obs.columns
        if any(keyword in col.lower() for keyword in ["doublet", "scrublet", "heuristic"])
    ]

    if not doublet_cols:
        log.warning("No doublet-related columns found in adata.obs")
        return {}

    log.info(f"Found doublet columns: {doublet_cols}")

    # Calculate per-sample statistics
    sample_stats = []
    unique_samples = adata.obs[sample_key].unique()
    if not isinstance(adata.obs[sample_key].dtype, pd.CategoricalDtype):
        unique_samples = sorted(unique_samples)

    for sample in unique_samples:
        sample_mask = adata.obs[sample_key] == sample
        sample_data = adata.obs.loc[sample_mask]

        stats = {
            "sample": sample,
            "total_cells": len(sample_data),
        }

        for col in doublet_cols:
            if col in sample_data.columns:
                col_data = sample_data[col].dropna()
                if (
                    pd.api.types.is_numeric_dtype(col_data)
                    and not pd.api.types.is_bool_dtype(col_data)
                    and col_data.nunique() > 2
                ):
                    # Continuous column (scores)
                    stats[f"{col}_mean"] = col_data.mean()
                    stats[f"{col}_median"] = col_data.median()
                    stats[f"{col}_std"] = col_data.std()
                elif pd.api.types.is_bool_dtype(col_data) or col_data.nunique() <= 2:
                    # Boolean/binary column (predictions)
                    positive_count = col_data.astype(bool).sum()
                    stats[f"{col}_count"] = positive_count
                    stats[f"{col}_percentage"] = (
                        (positive_count / len(sample_data) * 100) if len(sample_data) > 0 else 0
                    )

        sample_stats.append(stats)

    sample_df = pd.DataFrame(sample_stats).set_index("sample")

    global_stats = {"metric": "global", "total_cells": adata.n_obs}
    for col in doublet_cols:
        if col in adata.obs.columns:
            col_data = adata.obs[col].dropna()
            if (
                pd.api.types.is_numeric_dtype(col_data)
                and not pd.api.types.is_bool_dtype(col_data)
                and col_data.nunique() > 2
            ):
                global_stats[f"{col}_mean"] = col_data.mean()
                global_stats[f"{col}_median"] = col_data.median()
                global_stats[f"{col}_std"] = col_data.std()
            elif pd.api.types.is_bool_dtype(col_data) or col_data.nunique() <= 2:
                positive_count = col_data.astype(bool).sum()
                global_stats[f"{col}_count"] = positive_count
                global_stats[f"{col}_percentage"] = (
                    (positive_count / adata.n_obs * 100) if adata.n_obs > 0 else 0
                )

    global_df = pd.DataFrame([global_stats])

    if save_dir:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        if export_csv:
            sample_df.to_csv(save_dir / "doublet_stats_per_sample.csv")
            global_df.to_csv(save_dir / "doublet_stats_global.csv", index=False)
            log.info(f"Exported CSV files to {save_dir}")
        if export_xlsx:
            with pd.ExcelWriter(save_dir / "doublet_stats.xlsx") as writer:
                sample_df.to_excel(writer, sheet_name="per_sample")
                global_df.to_excel(writer, sheet_name="global", index=False)
            log.info(f"Exported Excel file to {save_dir / 'doublet_stats.xlsx'}")

    return {"sample": sample_df, "global": global_df}


# --- Main Functions ---
def predict_doublets(
    adata: AnnData,
    config: DoubletConfig,
    sample_key: str = "sampleID",
    cluster_key: Optional[str] = None,
    **kwargs,
) -> AnnData:
    """
    Enhanced doublet prediction with a clear, config-driven workflow.
    This version integrates a quantitative heuristic score with the algorithmic score for improved accuracy.

    Args:
        adata: AnnData object containing single-cell expression data.
        config: A `DoubletConfig` object that controls the entire workflow.
        sample_key: Key for sample identification in adata.obs.
        cluster_key: Optional key for cluster identification. When provided,
            per-cluster doublet fractions are computed and stored in ``.uns``.

    Returns:
        AnnData object with doublet predictions added to .obs and .obsm.
    """
    # === 1. CONFIGURATION SETUP ===
    cfg = DoubletConfig()
    if config is not None:
        cfg = config.model_copy(deep=True)
        if kwargs:
            cfg = cfg.model_copy(update=kwargs, deep=True)
    elif kwargs:
        cfg = cfg.model_copy(update=kwargs, deep=True)

    log.info("--- Running Final Doublet Prediction Workflow ---")

    # Validate input data
    if sample_key not in adata.obs.columns:
        raise ValueError(f"Sample key '{sample_key}' not found in adata.obs")

    samples = adata.obs[sample_key].unique()
    if len(samples) == 0:
        raise ValueError(f"No samples found for key '{sample_key}'")

    log.info(f"Starting doublet prediction for {adata.n_obs} cells across {len(samples)} samples")
    log.info(f"Configuration: method={cfg.method}, merge_strategy={cfg.merge_strategy}, ")
    log.info(
        "Doublet detection is grouped by '%s'. Use a capture/library/GEM-well key "
        "rather than a broad biological condition when possible.",
        sample_key,
    )

    # Initialize result columns
    algo_score_col = f"{cfg.method}_score"
    algo_pred_col = f"{cfg.method}_predicted"
    adata.obs[algo_score_col] = np.nan
    adata.obs[algo_pred_col] = False

    # Use a dispatcher for multi-algorithm support ---
    ALGORITHM_DISPATCHER = {
        "scrublet": _run_scrublet,
        "scanpy_scrublet": _run_scanpy_scrublet,
        "solo": _run_solo,
        "doubletdetection": _run_doubletdetection,
        "scdblfinder": _run_scdblfinder,
    }
    if cfg.method not in ALGORITHM_DISPATCHER:
        raise ValueError(
            f"Method '{cfg.method}' is not supported. Available: {list(ALGORITHM_DISPATCHER.keys())}"
        )

    # === 2. ALGORITHMIC DETECTION (Per-Sample) ===
    if cfg.run_algorithm:
        log.info(f"Running {cfg.method} doublet detection...")

        for sample in samples:
            log.info(f"Processing sample '{sample}' with {cfg.method}...")
            sample_mask = adata.obs[sample_key] == sample
            data_view = adata[sample_mask]

            if data_view.n_obs < 50:
                log.warning(
                    f"Skipping {sample}: fewer than 50 cells (insufficient for reliable doublet detection)."
                )
                continue

            scores, predicted = ALGORITHM_DISPATCHER[cfg.method](data_view, sample, cfg)

            if scores is not None and predicted is not None:
                adata.obs.loc[sample_mask, algo_score_col] = scores
                adata.obs.loc[sample_mask, algo_pred_col] = predicted
    else:
        log.info("Skipping algorithmic detection as per configuration (run_algorithm=False).")

    # === 3. HEURISTIC DETECTION (Global) ===
    adata.obs[HEURISTIC_PRED_COL] = False
    adata.obs[HEURISTIC_SCORE_COL] = 0.0
    if cfg.use_heuristics:
        log.info("Running quantitative heuristic analysis...")
        # Call the new heuristic function and receive its multiple outputs
        heuristic_pred, lineage_scores_df, heuristic_scores = _run_heuristic(
            adata,
            cfg,
            expected_rate=cfg.expected_doublet_rate,
            sample_key=sample_key if sample_key in adata.obs.columns else None,
        )

        # Store all the new results in the AnnData object
        adata.obsm["lineage_module_scores"] = lineage_scores_df  # Store detailed scores in .obsm
        adata.obs[HEURISTIC_PRED_COL] = heuristic_pred  # Store the binary call for simple stats
        adata.obs[HEURISTIC_SCORE_COL] = heuristic_scores  # Store the informative continuous score
        log.info(
            f"Heuristic analysis complete. Found {heuristic_pred.sum()} potential doublets based on score threshold."
        )

    # === 4. MERGE RESULTS ===
    log.info("Merging algorithmic and heuristic scores for final prediction...")
    merged_pred = _merge_doublet_predictions(
        adata,
        algorithm_score_col=algo_score_col,
        heuristic_score_col=HEURISTIC_SCORE_COL,
        strategy=cfg.merge_strategy,
        expected_rate=cfg.expected_doublet_rate,
        algo_weight=cfg.algorithm_weight,
        algorithm_score_normalization=cfg.algorithm_score_normalization,
        heuristic_score_normalization=cfg.heuristic_score_normalization,
        merge_rank_normalize=cfg.merge_rank_normalize,
        score_threshold=cfg.score_threshold,
        sample_key=sample_key,
        final_label_strategy=cfg.final_label_strategy,
        algorithm_pred_col=algo_pred_col,
    )
    adata.obs[FINAL_PRED_COL] = merged_pred

    external_evidence = _collect_external_doublet_evidence(
        adata,
        external_cols=list(cfg.external_doublet_cols),
        final_col=FINAL_PRED_COL,
        policy=cfg.external_doublet_policy,
    )

    # P4: Biology-aware downgrading for allowlisted co-expression pairs
    if cfg.ignore_coexpression_pairs and LINEAGE_SCORES_KEY in adata.obsm:
        lineage_scores_df = adata.obsm[LINEAGE_SCORES_KEY]
        allowlist_priority_records = []
        for lin1, lin2 in cfg.ignore_coexpression_pairs:
            if lin1 in lineage_scores_df.columns and lin2 in lineage_scores_df.columns:
                allowlist_mask = (lineage_scores_df[lin1] > 0.1) & (lineage_scores_df[lin2] > 0.1)
                # Only downgrade if algorithm did not flag the cell.
                # NaN values in algo_pred_col correspond to skipped samples and must
                # be treated as "not flagged" rather than "not doublet".
                algo_not_flagged = ~adata.obs[algo_pred_col].fillna(False).astype(bool)
                downgrade_mask = allowlist_mask & algo_not_flagged
                n_downgraded = int(downgrade_mask.sum())
                if n_downgraded > 0:
                    adata.obs.loc[downgrade_mask, FINAL_PRED_COL] = False
                    log.info(
                        f"Biology-aware downgrade: {n_downgraded} cells allowlisted "
                        f"({lin1}+{lin2}) demoted from final doublet call."
                    )
                allowlist_priority_records.append(
                    {
                        "pair": [lin1, lin2],
                        "n_allowlisted": int(allowlist_mask.sum()),
                        "n_algorithm_positive": int((allowlist_mask & ~algo_not_flagged).sum()),
                        "n_downgraded_heuristic_only": n_downgraded,
                    }
                )
    else:
        allowlist_priority_records = []

    risk_decomposition = _add_doublet_risk_decomposition(
        adata,
        algorithm_score_col=algo_score_col,
        algorithm_pred_col=algo_pred_col,
        heuristic_score_col=HEURISTIC_SCORE_COL,
        expected_rate=cfg.expected_doublet_rate,
        sample_key=sample_key,
    )
    if allowlist_priority_records:
        risk_decomposition["allowlist_priority"] = allowlist_priority_records

    # P5: Record threshold evidence
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {}).setdefault("doublet_params", {})
    doublet_rate = float(adata.obs[FINAL_PRED_COL].mean())
    expected_rate = cfg.expected_doublet_rate
    expected_rate_float = (
        float(expected_rate)
        if isinstance(expected_rate, (int, float))
        else float(np.mean(list(expected_rate.values())))
        if isinstance(expected_rate, dict)
        else 0.10
    )
    adata.uns["sclucid"]["qc"]["doublet_params"].update(
        {
            "merge_strategy": cfg.merge_strategy,
            "algorithm_weight": cfg.algorithm_weight,
            "expected_doublet_rate": cfg.expected_doublet_rate,
            "score_threshold": cfg.score_threshold,
            "method": cfg.method,
            "detection_group_key": sample_key,
            "final_label_strategy": cfg.final_label_strategy,
            "label_strategy_note": (
                "expected_rate_rank thresholds merged evidence scores by expected rate. "
                "algorithm_label uses the algorithm binary call directly and keeps "
                "heuristic evidence for review."
            ),
            "heuristic_evidence": {
                "enabled": bool(cfg.use_heuristics),
                "semantics": (
                    "Marker co-expression is heterotypic doublet evidence, not a "
                    "standalone truth label or calibrated probability."
                ),
                "marker_source": "manual_marker_configs"
                if cfg.marker_configs is not None
                else "marker_manager",
                "marker_species": cfg.marker_species,
                "marker_tissue": cfg.marker_tissue,
                "ignore_coexpression_pairs": [list(pair) for pair in cfg.ignore_coexpression_pairs],
                "score_normalization": cfg.heuristic_score_normalization,
            },
            "external_doublet_evidence": external_evidence,
            "risk_decomposition": risk_decomposition,
            "threshold_evidence": {
                "heuristic_threshold_method": (
                    "adaptive_expected_rate"
                    if cfg.expected_doublet_rate is not None
                    else "legacy_quantile_90"
                ),
                "heuristic_expected_rate": cfg.expected_doublet_rate,
                "actual_doublet_rate": round(doublet_rate, 4),
                "rate_vs_expected_ratio": (
                    round(doublet_rate / expected_rate_float, 2)
                    if expected_rate_float > 0
                    else None
                ),
            },
        }
    )

    # === 5. SUMMARY STATISTICS ===
    log.info("\n" + "=" * 50)
    log.info("DOUBLET DETECTION SUMMARY")
    log.info("=" * 50)

    total_cells = adata.n_obs

    # Algorithm results
    algo_count = adata.obs[algo_pred_col].sum()
    log.info(f"Algorithm ({cfg.method}): {algo_count} doublets ({algo_count / total_cells:.2%})")

    # Heuristic results
    if cfg.use_heuristics:
        heur_count = adata.obs[HEURISTIC_PRED_COL].sum()
        log.info(f"Heuristic: {heur_count} doublets ({heur_count / total_cells:.2%})")

        # Overlap analysis
        overlap_count = (adata.obs[algo_pred_col] & adata.obs[HEURISTIC_PRED_COL]).sum()
        log.info(f"Overlap: {overlap_count} doublets ({overlap_count / total_cells:.2%})")

    # Final merged results
    final_count = adata.obs[FINAL_PRED_COL].sum()
    log.info(f"Final merged: {final_count} doublets ({final_count / total_cells:.2%})")

    # Per-sample breakdown
    log.info("\nPer-sample statistics:")
    for sample in samples:
        sample_mask = adata.obs[sample_key] == sample
        sample_total = sample_mask.sum()
        sample_doublets = adata.obs[FINAL_PRED_COL][sample_mask].sum()
        sample_rate = sample_doublets / sample_total
        log.info(f"  {sample}: {sample_doublets}/{sample_total} doublets ({sample_rate:.2%})")

    log.info("=" * 50)

    # === 6. CLUSTER-LEVEL DOUBLET AUDIT (P2) ===
    if cluster_key and cluster_key in adata.obs.columns:
        cluster_series = adata.obs[cluster_key].astype(str)
        cluster_summary = []
        for cluster in cluster_series.unique():
            mask = cluster_series == cluster
            n_total = int(mask.sum())
            n_doublets = int(adata.obs.loc[mask, FINAL_PRED_COL].sum())
            fraction = n_doublets / n_total if n_total > 0 else 0.0
            cluster_summary.append(
                {
                    "cluster": str(cluster),
                    "n_cells": n_total,
                    "n_doublets": n_doublets,
                    "doublet_fraction": float(fraction),
                    "flagged": bool(fraction > 0.30 and n_total >= 10),
                }
            )
        cluster_df = pd.DataFrame(cluster_summary)
        adata.uns.setdefault("sclucid", {}).setdefault("qc", {}).setdefault("doublet_params", {})
        adata.uns["sclucid"]["qc"]["doublet_params"]["cluster_summary"] = cluster_df
        n_flagged = int(cluster_df["flagged"].sum())
        if n_flagged > 0:
            log.warning(
                f"{n_flagged} cluster(s) have >30% doublet fraction and may be doublet clusters."
            )

    # === 7. Reporting & Visualization ===
    if cfg.plot_summary:
        save_path = Path(cfg.save_dir) if cfg.save_dir else None
        _plot_doublet_summary(
            adata=adata,
            sample_key=sample_key,
            save_dir=save_path,
            show=cfg.show_plots,
            plot_bar=cfg.plot_bar,
            plot_scatter=cfg.plot_scatter,
            plot_upset=cfg.plot_upset,
        )

    if cfg.export_stats and cfg.save_dir:
        _export_doublet_stats(adata, sample_key, Path(cfg.save_dir))

    log.info("Doublet prediction workflow completed.")
    record_qc_artifact_contract(adata)

    return adata
