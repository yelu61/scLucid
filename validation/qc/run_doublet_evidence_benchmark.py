#!/usr/bin/env python3
"""Run a lightweight demuxlet-grounded doublet evidence benchmark."""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scLucid.qc.config import DoubletConfig
from scLucid.qc.doublet import COMBINED_SCORE_COL, FINAL_PRED_COL, predict_doublets
from validation.dataset_registry import DATASETS


def _counts_matrix(adata: ad.AnnData):
    return adata.layers["counts"] if "counts" in adata.layers else adata.X


def _ensure_metrics(adata: ad.AnnData) -> None:
    X = _counts_matrix(adata)
    if "total_counts" not in adata.obs:
        adata.obs["total_counts"] = np.asarray(X.sum(axis=1)).ravel()
    if "n_genes_by_counts" not in adata.obs:
        if sp.issparse(X):
            adata.obs["n_genes_by_counts"] = np.asarray((X > 0).sum(axis=1)).ravel()
        else:
            adata.obs["n_genes_by_counts"] = (np.asarray(X) > 0).sum(axis=1)


def _metrics(truth: pd.Series, predicted: pd.Series) -> dict[str, Any]:
    tp = int((truth & predicted).sum())
    fp = int((~truth & predicted).sum())
    tn = int((~truth & ~predicted).sum())
    fn = int((truth & ~predicted).sum())
    precision = tp / max(tp + fp, 1)
    recall = tp / max(tp + fn, 1)
    f1 = 2 * precision * recall / max(precision + recall, 1e-12)
    return {
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def _auc(truth: pd.Series, score: pd.Series) -> float | None:
    try:
        from sklearn.metrics import roc_auc_score

        valid = score.notna()
        if valid.sum() == 0 or truth.loc[valid].nunique() < 2:
            return None
        return float(roc_auc_score(truth.loc[valid].astype(bool), score.loc[valid].astype(float)))
    except Exception:
        return None


def _stratified_subset(
    adata: ad.AnnData,
    labels: pd.Series,
    max_cells: int | None,
    seed: int,
) -> ad.AnnData:
    if max_cells is None or max_cells <= 0 or adata.n_obs <= max_cells:
        return adata
    rng = np.random.default_rng(seed)
    keep: list[int] = []
    label_counts = labels.astype(str).value_counts()
    for label, count in label_counts.items():
        idx = np.flatnonzero(labels.astype(str).to_numpy() == label)
        n_label = max(1, int(round(max_cells * count / adata.n_obs)))
        n_label = min(n_label, len(idx))
        keep.extend(rng.choice(idx, size=n_label, replace=False).tolist())
    if len(keep) > max_cells:
        keep = rng.choice(np.array(keep), size=max_cells, replace=False).tolist()
    keep = sorted(set(int(i) for i in keep))
    return adata[keep].copy()


def _heuristic_predictions(adata: ad.AnnData) -> dict[str, pd.Series]:
    obs = adata.obs
    n_cells = adata.n_obs
    expected_rate = float(obs["doublet_ground_truth"].mean())
    n_expected = max(1, int(round(expected_rate * n_cells)))
    score = (
        obs["total_counts"].rank(pct=True).astype(float)
        + obs["n_genes_by_counts"].rank(pct=True).astype(float)
    ) / 2.0
    return {
        "lineage_size_heuristic": score.rank(method="first", ascending=False) <= n_expected,
        "top_counts_equal_rate": obs["total_counts"].rank(method="first", ascending=False)
        <= n_expected,
        "top_genes_equal_rate": obs["n_genes_by_counts"].rank(method="first", ascending=False)
        <= n_expected,
    }


def _method_available(method: str) -> tuple[bool, str]:
    package = {
        "scrublet": "scrublet",
        "scdblfinder": "pyscdblfinder",
        "scdblfinder_python": "pyscdblfinder",
        "scdblfinder_python_pyscdblfinder": "pyscdblfinder",
    }.get(method)
    if package is None:
        return True, ""
    return bool(importlib.util.find_spec(package)), package


def _run_sclucid_doublet_method(
    adata: ad.AnnData,
    method: str,
    expected_rate: float,
    sample_key: str,
    seed: int,
    *,
    use_heuristics: bool = False,
    algorithm_weight: float = 0.7,
    merge_strategy: str = "weighted_average",
) -> dict[str, Any]:
    available, package = _method_available(method)
    if not available:
        return {
            "status": "dependency_missing",
            "package": package,
            "predicted": pd.Series(False, index=adata.obs_names),
            "score": pd.Series(np.nan, index=adata.obs_names),
            "runtime_seconds": 0.0,
            "error": f"Missing optional dependency: {package}",
            "base_method": method,
            "sc_method": method,
            "use_heuristics": use_heuristics,
            "algorithm_weight": algorithm_weight if use_heuristics else None,
            "merge_strategy": merge_strategy if use_heuristics else "",
        }
    start = time.perf_counter()
    try:
        sc_method = "scdblfinder" if method.startswith("scdblfinder_python") else method
        cfg = DoubletConfig(
            method=sc_method,
            run_algorithm=True,
            use_heuristics=use_heuristics,
            merge_strategy=merge_strategy,
            algorithm_weight=algorithm_weight,
            expected_doublet_rate=expected_rate,
            random_state=seed,
            plot_summary=False,
            plot_bar=False,
            plot_scatter=False,
            plot_upset=False,
            show_plots=False,
            export_stats=False,
            scr_plot_umap=False,
            scdblfinder_iter=2,
            scdblfinder_dims=15,
            scdblfinder_include_pcs=10,
            scdblfinder_dbr=expected_rate,
        )
        result = predict_doublets(adata.copy(), config=cfg, sample_key=sample_key)
        pred_col = f"{sc_method}_predicted"
        score_col = f"{sc_method}_score"
        if use_heuristics:
            pred = result.obs[FINAL_PRED_COL]
            score = result.obs.get(COMBINED_SCORE_COL, result.obs.get("doublet_score"))
        else:
            pred = result.obs.get(pred_col, result.obs[FINAL_PRED_COL])
            score = result.obs.get(score_col, result.obs.get("doublet_score"))
        return {
            "status": "ok",
            "package": package,
            "predicted": pred.astype(bool),
            "score": pd.to_numeric(score, errors="coerce"),
            "runtime_seconds": time.perf_counter() - start,
            "error": "",
            "base_method": method,
            "sc_method": sc_method,
            "use_heuristics": use_heuristics,
            "algorithm_weight": algorithm_weight if use_heuristics else None,
            "merge_strategy": merge_strategy if use_heuristics else "",
        }
    except Exception as exc:
        return {
            "status": "failed",
            "package": package,
            "predicted": pd.Series(False, index=adata.obs_names),
            "score": pd.Series(np.nan, index=adata.obs_names),
            "runtime_seconds": time.perf_counter() - start,
            "error": f"{type(exc).__name__}: {exc}",
            "base_method": method,
            "sc_method": method,
            "use_heuristics": use_heuristics,
            "algorithm_weight": algorithm_weight if use_heuristics else None,
            "merge_strategy": merge_strategy if use_heuristics else "",
        }


def _load_r_scdblfinder_reference(path: Path) -> pd.DataFrame:
    reference = pd.read_csv(path)
    lower = {str(col).lower(): col for col in reference.columns}
    cell_col = lower.get("cell") or lower.get("barcode") or lower.get("cell_id")
    score_col = (
        lower.get("score") or lower.get("scdblfinder.score") or lower.get("scdblfinder_score")
    )
    class_col = (
        lower.get("predicted")
        or lower.get("class")
        or lower.get("scdblfinder.class")
        or lower.get("scdblfinder_class")
    )
    if cell_col is None or score_col is None or class_col is None:
        raise ValueError(
            "R scDblFinder reference CSV must contain cell/barcode, "
            "score/scDblFinder.score, and predicted/class/scDblFinder.class columns."
        )
    output = pd.DataFrame(index=reference[cell_col].astype(str))
    output["r_score"] = pd.to_numeric(reference[score_col], errors="coerce").to_numpy()
    class_text = reference[class_col].astype(str).str.lower()
    output["r_predicted"] = class_text.isin(["true", "1", "doublet", "yes"]).to_numpy()
    return output


def _scdblfinder_python_vs_r_reference_rows(
    *,
    dataset: str,
    truth: pd.Series,
    predictions: dict[str, pd.Series],
    scores: dict[str, pd.Series],
    r_reference_path: Path | None,
) -> list[dict[str, Any]]:
    method = "scdblfinder_python_pyscdblfinder"
    if method not in predictions:
        return [
            {
                "dataset": dataset,
                "comparison": "pyscdblfinder_vs_bioconductor_scdblfinder",
                "python_method": method,
                "reference_status": "python_method_missing",
                "n_overlap_cells": 0,
                "risk_note": "Run with --methods scdblfinder_python_pyscdblfinder.",
            }
        ]
    if r_reference_path is None or not r_reference_path.exists():
        return [
            {
                "dataset": dataset,
                "comparison": "pyscdblfinder_vs_bioconductor_scdblfinder",
                "python_method": method,
                "r_reference_path": "" if r_reference_path is None else str(r_reference_path),
                "reference_status": "missing_reference",
                "n_overlap_cells": 0,
                "risk_note": (
                    "Provide a CSV exported from Bioconductor scDblFinder on the same cells "
                    "to quantify parity."
                ),
            }
        ]
    try:
        reference = _load_r_scdblfinder_reference(r_reference_path)
    except Exception as exc:
        return [
            {
                "dataset": dataset,
                "comparison": "pyscdblfinder_vs_bioconductor_scdblfinder",
                "python_method": method,
                "r_reference_path": str(r_reference_path),
                "reference_status": "invalid_reference",
                "n_overlap_cells": 0,
                "risk_note": f"{type(exc).__name__}: {exc}",
            }
        ]

    overlap = predictions[method].index.intersection(reference.index).intersection(truth.index)
    if len(overlap) == 0:
        return [
            {
                "dataset": dataset,
                "comparison": "pyscdblfinder_vs_bioconductor_scdblfinder",
                "python_method": method,
                "r_reference_path": str(r_reference_path),
                "reference_status": "no_overlap",
                "n_overlap_cells": 0,
                "risk_note": "R reference cell IDs do not overlap with benchmark usable cells.",
            }
        ]

    py_pred = predictions[method].loc[overlap].astype(bool)
    r_pred = reference.loc[overlap, "r_predicted"].astype(bool)
    py_score = scores[method].loc[overlap].astype(float)
    r_score = reference.loc[overlap, "r_score"].astype(float)
    truth_overlap = truth.loc[overlap].astype(bool)
    score_spearman = float(py_score.corr(r_score, method="spearman"))
    prediction_agreement = float((py_pred == r_pred).mean())
    prediction_jaccard = float((py_pred & r_pred).sum() / max((py_pred | r_pred).sum(), 1))
    py_metrics = _metrics(truth_overlap, py_pred)
    r_metrics = _metrics(truth_overlap, r_pred)
    review_required = bool(
        score_spearman < 0.8
        or prediction_agreement < 0.9
        or abs(py_metrics["f1"] - r_metrics["f1"]) > 0.1
    )
    return [
        {
            "dataset": dataset,
            "comparison": "pyscdblfinder_vs_bioconductor_scdblfinder",
            "python_method": method,
            "r_reference_path": str(r_reference_path),
            "reference_status": "ok",
            "n_overlap_cells": int(len(overlap)),
            "score_spearman": score_spearman,
            "prediction_agreement": prediction_agreement,
            "prediction_jaccard": prediction_jaccard,
            "python_auc": _auc(truth_overlap, py_score),
            "r_auc": _auc(truth_overlap, r_score),
            "python_f1": py_metrics["f1"],
            "r_f1": r_metrics["f1"],
            "f1_delta_python_minus_r": py_metrics["f1"] - r_metrics["f1"],
            "review_required": review_required,
            "risk_note": (
                "Parity with Bioconductor scDblFinder is below target." if review_required else ""
            ),
        }
    ]


def _scdblfinder_python_vs_r_reference_details(
    *,
    dataset: str,
    obs: pd.DataFrame,
    truth: pd.Series,
    predictions: dict[str, pd.Series],
    scores: dict[str, pd.Series],
    r_reference_path: Path | None,
    min_group_cells: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    method = "scdblfinder_python_pyscdblfinder"
    if method not in predictions or r_reference_path is None or not r_reference_path.exists():
        return [], []
    try:
        reference = _load_r_scdblfinder_reference(r_reference_path)
    except Exception:
        return [], []

    overlap = predictions[method].index.intersection(reference.index).intersection(truth.index)
    if len(overlap) == 0:
        return [], []

    detail = pd.DataFrame(index=overlap)
    detail["truth"] = truth.loc[overlap].astype(bool)
    detail["python_predicted"] = predictions[method].loc[overlap].astype(bool)
    detail["r_predicted"] = reference.loc[overlap, "r_predicted"].astype(bool)
    detail["python_score"] = scores[method].loc[overlap].astype(float)
    detail["r_score"] = reference.loc[overlap, "r_score"].astype(float)
    detail["score_delta_python_minus_r"] = detail["python_score"] - detail["r_score"]
    detail["prediction_disagreement"] = detail["python_predicted"] != detail["r_predicted"]
    for key in ("sample", "donor", "condition", "cell_type", "cell_subtype"):
        if key in obs:
            detail[key] = obs.loc[overlap, key].astype(str)

    group_rows: list[dict[str, Any]] = []
    for group_key in ("sample", "donor", "condition", "cell_type", "cell_subtype"):
        if group_key not in detail:
            continue
        for group, frame in detail.groupby(group_key, observed=True):
            if len(frame) < min_group_cells:
                continue
            py_metrics = _metrics(frame["truth"], frame["python_predicted"])
            r_metrics = _metrics(frame["truth"], frame["r_predicted"])
            group_rows.append(
                {
                    "dataset": dataset,
                    "group_key": group_key,
                    "group": str(group),
                    "n_cells": int(len(frame)),
                    "ground_truth_doublets": int(frame["truth"].sum()),
                    "python_predicted_doublets": int(frame["python_predicted"].sum()),
                    "r_predicted_doublets": int(frame["r_predicted"].sum()),
                    "prediction_agreement": float(
                        (frame["python_predicted"] == frame["r_predicted"]).mean()
                    ),
                    "prediction_disagreement_rate": float(frame["prediction_disagreement"].mean()),
                    "prediction_jaccard": float(
                        (frame["python_predicted"] & frame["r_predicted"]).sum()
                        / max((frame["python_predicted"] | frame["r_predicted"]).sum(), 1)
                    ),
                    "score_spearman": float(
                        frame["python_score"].corr(frame["r_score"], method="spearman")
                    ),
                    "python_auc": _auc(frame["truth"], frame["python_score"]),
                    "r_auc": _auc(frame["truth"], frame["r_score"]),
                    "python_f1": py_metrics["f1"],
                    "r_f1": r_metrics["f1"],
                    "f1_delta_python_minus_r": py_metrics["f1"] - r_metrics["f1"],
                    "review_required": bool(
                        frame["prediction_disagreement"].mean() > 0.1
                        or abs(py_metrics["f1"] - r_metrics["f1"]) > 0.1
                    ),
                }
            )

    disagreement = detail[detail["prediction_disagreement"]].copy()
    disagreement["abs_score_delta"] = disagreement["score_delta_python_minus_r"].abs()
    disagreement = disagreement.sort_values("abs_score_delta", ascending=False).head(500)
    disagreement_rows = []
    for cell, row in disagreement.iterrows():
        record = {
            "dataset": dataset,
            "cell": str(cell),
            "truth": bool(row["truth"]),
            "python_predicted": bool(row["python_predicted"]),
            "r_predicted": bool(row["r_predicted"]),
            "python_score": float(row["python_score"]),
            "r_score": float(row["r_score"]),
            "score_delta_python_minus_r": float(row["score_delta_python_minus_r"]),
        }
        for key in ("sample", "donor", "condition", "cell_type", "cell_subtype"):
            if key in row:
                record[key] = row[key]
        disagreement_rows.append(record)
    return group_rows, disagreement_rows


def _threshold_calibration_rows(
    *,
    dataset: str,
    truth: pd.Series,
    predictions: dict[str, pd.Series],
    scores: dict[str, pd.Series],
    target_recalls: tuple[float, ...] = (0.5, 0.7, 0.8),
) -> list[dict[str, Any]]:
    """Scan score thresholds to quantify benchmark-driven calibration behavior."""
    rows: list[dict[str, Any]] = []
    for method, score in scores.items():
        score = pd.to_numeric(score, errors="coerce")
        valid = score.notna() & truth.notna()
        if valid.sum() == 0 or score.loc[valid].nunique() <= 5:
            continue
        truth_valid = truth.loc[valid].astype(bool)
        score_valid = score.loc[valid].astype(float)
        default_pred = predictions[method].loc[valid].astype(bool)
        default_metrics = _metrics(truth_valid, default_pred)
        default_predicted_rate = float(default_pred.mean())

        thresholds = np.unique(
            np.nanquantile(score_valid.to_numpy(), np.linspace(0.01, 0.99, 99))
        )
        scan_rows: list[dict[str, Any]] = []
        for threshold in thresholds:
            calibrated_pred = score_valid >= threshold
            metric = _metrics(truth_valid, calibrated_pred)
            scan_rows.append(
                {
                    "threshold": float(threshold),
                    "predicted_rate": float(calibrated_pred.mean()),
                    **metric,
                }
            )
        if not scan_rows:
            continue
        scan_df = pd.DataFrame(scan_rows)
        best = scan_df.sort_values(["f1", "recall", "precision"], ascending=False).iloc[0]
        for target_recall in target_recalls:
            eligible = scan_df[scan_df["recall"] >= target_recall]
            if eligible.empty:
                chosen = best
                selection_rule = "best_f1_no_threshold_reaches_target_recall"
            else:
                chosen = eligible.sort_values(
                    ["f1", "precision", "threshold"],
                    ascending=[False, False, False],
                ).iloc[0]
                selection_rule = f"best_f1_at_recall_ge_{target_recall:.2f}"
            rows.append(
                {
                    "dataset": dataset,
                    "method": method,
                    "target_recall": target_recall,
                    "selection_rule": selection_rule,
                    "calibrated_threshold": float(chosen["threshold"]),
                    "calibrated_predicted_rate": float(chosen["predicted_rate"]),
                    "calibrated_precision": float(chosen["precision"]),
                    "calibrated_recall": float(chosen["recall"]),
                    "calibrated_f1": float(chosen["f1"]),
                    "best_f1_threshold": float(best["threshold"]),
                    "best_f1": float(best["f1"]),
                    "default_predicted_rate": default_predicted_rate,
                    "default_precision": float(default_metrics["precision"]),
                    "default_recall": float(default_metrics["recall"]),
                    "default_f1": float(default_metrics["f1"]),
                    "score_auc": _auc(truth_valid, score_valid),
                    "recall_gain_vs_default": float(chosen["recall"] - default_metrics["recall"]),
                    "f1_gain_vs_default": float(chosen["f1"] - default_metrics["f1"]),
                    "review_required": bool(
                        default_metrics["recall"] < target_recall
                        and chosen["recall"] > default_metrics["recall"] + 0.05
                    ),
                    "risk_note": (
                        "Default threshold under-recovers demuxlet doublets; consider benchmark-calibrated threshold."
                        if default_metrics["recall"] < target_recall
                        and chosen["recall"] > default_metrics["recall"] + 0.05
                        else ""
                    ),
                }
            )
    return rows


def _algorithm_weight_recommendation_rows(
    evidence_rows: list[dict[str, Any]],
    *,
    min_f1_gain: float = 0.02,
    max_precision_drop: float = 0.10,
    preferred_weight: float = 0.7,
) -> list[dict[str, Any]]:
    """Recommend an algorithm/heuristic merge weight from benchmark evidence."""
    evidence = pd.DataFrame(evidence_rows)
    if evidence.empty or "uses_heuristics" not in evidence.columns:
        return []

    rows: list[dict[str, Any]] = []
    grouped = evidence[evidence["uses_heuristics"] == True].groupby(  # noqa: E712
        ["dataset", "base_method"],
        dropna=False,
    )
    for (dataset, base_method), merged in grouped:
        merged = merged[merged["method_status"] == "ok"].copy()
        if merged.empty:
            continue
        baseline = evidence[
            (evidence["dataset"] == dataset)
            & (evidence["method"] == base_method)
            & (evidence["uses_heuristics"] == False)  # noqa: E712
            & (evidence["method_status"] == "ok")
        ]
        baseline_row = baseline.iloc[0] if not baseline.empty else None
        merged["auc_sort"] = pd.to_numeric(merged["score_auc"], errors="coerce").fillna(-1.0)
        merged["weight_distance_to_default"] = (
            pd.to_numeric(merged["algorithm_weight"], errors="coerce") - preferred_weight
        ).abs()
        max_f1 = float(merged["f1"].max())
        near_best = merged[merged["f1"] >= max_f1 - 0.005].copy()
        best = near_best.sort_values(
            ["weight_distance_to_default", "recall", "precision", "auc_sort"],
            ascending=[True, False, False, False],
        ).iloc[0]

        base_f1 = float(baseline_row["f1"]) if baseline_row is not None else np.nan
        base_precision = (
            float(baseline_row["precision"]) if baseline_row is not None else np.nan
        )
        base_recall = float(baseline_row["recall"]) if baseline_row is not None else np.nan
        base_auc = (
            None
            if baseline_row is None or pd.isna(baseline_row.get("score_auc"))
            else float(baseline_row.get("score_auc"))
        )
        f1_gain = float(best["f1"] - base_f1) if not pd.isna(base_f1) else np.nan
        precision_delta = (
            float(best["precision"] - base_precision) if not pd.isna(base_precision) else np.nan
        )
        recall_delta = (
            float(best["recall"] - base_recall) if not pd.isna(base_recall) else np.nan
        )
        auc_delta = (
            None
            if base_auc is None or pd.isna(best.get("score_auc"))
            else float(best["score_auc"] - base_auc)
        )

        review_required = bool(
            pd.isna(f1_gain)
            or f1_gain < min_f1_gain
            or (not pd.isna(precision_delta) and precision_delta < -max_precision_drop)
        )
        if pd.isna(f1_gain):
            default_mode = "manual_review"
            risk_note = "Algorithm-only baseline unavailable; review before changing defaults."
        elif f1_gain < min_f1_gain:
            default_mode = "algorithm_only_with_heuristic_review_evidence"
            risk_note = (
                "Heuristic fusion does not materially improve F1 on this dataset; "
                "prefer algorithm-only or keep heuristic as review evidence."
            )
        elif not pd.isna(precision_delta) and precision_delta < -max_precision_drop:
            default_mode = "algorithm_plus_heuristic_with_manual_review"
            risk_note = (
                "Heuristic fusion improves recovery but loses precision; use with manual review."
            )
        else:
            default_mode = "algorithm_plus_heuristic"
            risk_note = (
                "Heuristic fusion improves demuxlet-grounded recovery without a large precision loss."
            )

        rows.append(
            {
                "dataset": dataset,
                "base_method": base_method,
                "recommended_default_mode": default_mode,
                "recommended_method": best["method"],
                "recommended_algorithm_weight": float(best["algorithm_weight"]),
                "recommended_merge_strategy": best["merge_strategy"],
                "baseline_f1": base_f1,
                "baseline_precision": base_precision,
                "baseline_recall": base_recall,
                "baseline_score_auc": base_auc,
                "recommended_f1": float(best["f1"]),
                "recommended_precision": float(best["precision"]),
                "recommended_recall": float(best["recall"]),
                "recommended_score_auc": (
                    None if pd.isna(best.get("score_auc")) else float(best.get("score_auc"))
                ),
                "recommended_predicted_rate": float(best["predicted_rate"]),
                "expected_rate_from_demuxlet": float(best["expected_rate_from_demuxlet"]),
                "f1_delta_vs_algorithm_only": f1_gain,
                "precision_delta_vs_algorithm_only": precision_delta,
                "recall_delta_vs_algorithm_only": recall_delta,
                "score_auc_delta_vs_algorithm_only": auc_delta,
                "review_required": review_required,
                "risk_note": risk_note,
            }
        )
    return rows


def _benchmark_summary_payload(
    *,
    evidence_rows: list[dict[str, Any]],
    parity_rows: list[dict[str, Any]],
    group_rows: list[dict[str, Any]],
    calibration_rows: list[dict[str, Any]],
    weight_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    """Compact evidence payload suitable for QC report attachment."""
    evidence = pd.DataFrame(evidence_rows)
    parity = pd.DataFrame(parity_rows)
    groups = pd.DataFrame(group_rows)
    calibration = pd.DataFrame(calibration_rows)
    weights = pd.DataFrame(weight_rows)
    payload: dict[str, Any] = {"schema_version": "doublet_benchmark_evidence_v1"}
    if not evidence.empty:
        best = evidence.sort_values(["f1", "score_auc"], ascending=False).iloc[0]
        payload["best_method"] = str(best["method"])
        payload["best_method_f1"] = float(best["f1"])
        payload["best_method_auc"] = (
            None if pd.isna(best.get("score_auc")) else float(best.get("score_auc"))
        )
    if not parity.empty and parity.iloc[0].get("reference_status") == "ok":
        row = parity.iloc[0]
        payload["python_r_parity"] = {
            "score_spearman": float(row["score_spearman"]),
            "prediction_agreement": float(row["prediction_agreement"]),
            "prediction_jaccard": float(row["prediction_jaccard"]),
            "review_required": bool(row["review_required"]),
        }
    if not groups.empty:
        top_groups = groups.sort_values(
            ["prediction_disagreement_rate", "n_cells"],
            ascending=[False, False],
        ).head(5)
        payload["top_python_r_disagreement_groups"] = [
            {
                "group_key": str(row["group_key"]),
                "group": str(row["group"]),
                "n_cells": int(row["n_cells"]),
                "prediction_disagreement_rate": float(row["prediction_disagreement_rate"]),
                "score_spearman": float(row["score_spearman"]),
            }
            for _, row in top_groups.iterrows()
        ]
    if not calibration.empty:
        flagged = calibration[calibration["review_required"]].copy()
        if not flagged.empty:
            payload["threshold_calibration_review"] = [
                {
                    "method": str(row["method"]),
                    "target_recall": float(row["target_recall"]),
                    "calibrated_threshold": float(row["calibrated_threshold"]),
                    "default_recall": float(row["default_recall"]),
                    "calibrated_recall": float(row["calibrated_recall"]),
                    "recall_gain_vs_default": float(row["recall_gain_vs_default"]),
                }
                for _, row in flagged.sort_values(
                    ["recall_gain_vs_default", "f1_gain_vs_default"],
                    ascending=False,
                ).head(5).iterrows()
            ]
    if not weights.empty:
        payload["algorithm_weight_recommendations"] = [
            {
                "base_method": str(row["base_method"]),
                "recommended_default_mode": str(row["recommended_default_mode"]),
                "recommended_algorithm_weight": float(row["recommended_algorithm_weight"]),
                "recommended_method": str(row["recommended_method"]),
                "f1_delta_vs_algorithm_only": float(row["f1_delta_vs_algorithm_only"]),
                "recall_delta_vs_algorithm_only": float(row["recall_delta_vs_algorithm_only"]),
                "precision_delta_vs_algorithm_only": float(
                    row["precision_delta_vs_algorithm_only"]
                ),
                "review_required": bool(row["review_required"]),
                "risk_note": str(row["risk_note"]),
            }
            for _, row in weights.sort_values(
                ["review_required", "f1_delta_vs_algorithm_only"],
                ascending=[True, False],
            ).iterrows()
        ]
    return payload


def run(
    output_dir: Path,
    max_cells: int | None,
    seed: int,
    methods: list[str],
    r_scdblfinder_reference: Path | None,
    min_parity_group_cells: int,
    algorithm_weights: list[float],
    merge_strategy: str,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    spec = next(dataset for dataset in DATASETS if dataset.key == "kang2018.pbmc")
    adata = ad.read_h5ad(spec.path)
    original_labels = adata.obs["demuxlet_multiplets"].astype(str)
    adata = _stratified_subset(adata, original_labels, max_cells=max_cells, seed=seed)
    _ensure_metrics(adata)

    labels = adata.obs["demuxlet_multiplets"].astype(str)
    usable = labels.isin(["singlet", "doublet"])
    truth = adata.obs.loc[usable, "doublet_ground_truth"].astype(bool)
    work = adata[usable].copy()
    expected_rate = float(truth.mean())
    predictions: dict[str, pd.Series] = {}
    scores: dict[str, pd.Series] = {}
    method_status: dict[str, dict[str, Any]] = {}

    for method, predicted in _heuristic_predictions(work).items():
        predictions[method] = predicted.astype(bool)
        scores[method] = predicted.astype(float)
        method_status[method] = {
            "status": "ok",
            "package": "",
            "runtime_seconds": 0.0,
            "error": "",
            "fallback_mode": True,
            "base_method": method,
            "uses_heuristics": False,
            "algorithm_weight": None,
            "merge_strategy": "",
            "recommendation_role": "fallback_heuristic_baseline",
        }

    full_predictions: dict[str, pd.Series] = {}
    for method in methods:
        run_specs: list[tuple[str, bool, float | None]] = [(method, False, None)]
        run_specs.extend(
            (
                f"{method}_plus_heuristic_w{weight:.2f}",
                True,
                float(weight),
            )
            for weight in algorithm_weights
        )
        for method_label, use_heuristics, weight in run_specs:
            result = _run_sclucid_doublet_method(
                adata,
                method=method,
                expected_rate=expected_rate,
                sample_key="sample",
                seed=seed,
                use_heuristics=use_heuristics,
                algorithm_weight=0.7 if weight is None else weight,
                merge_strategy=merge_strategy,
            )
            full_predictions[method_label] = result["predicted"].astype(bool)
            predictions[method_label] = result["predicted"].loc[usable].astype(bool)
            scores[method_label] = result["score"].loc[usable]
            method_status[method_label] = {
                "status": result["status"],
                "package": result["package"],
                "runtime_seconds": result["runtime_seconds"],
                "error": result["error"],
                "fallback_mode": result["status"] != "ok",
                "base_method": method,
                "uses_heuristics": use_heuristics,
                "algorithm_weight": weight,
                "merge_strategy": merge_strategy if use_heuristics else "",
                "recommendation_role": (
                    "algorithm_plus_heuristic"
                    if use_heuristics
                    else "algorithm_only_baseline"
                ),
            }

    evidence_rows: list[dict[str, Any]] = []
    for method, predicted in predictions.items():
        metric = _metrics(truth, predicted.astype(bool))
        status = method_status[method]
        evidence_rows.append(
            {
                "dataset": spec.key,
                "method": method,
                "base_method": status["base_method"],
                "method_status": status["status"],
                "optional_dependency": status["package"],
                "recommendation_role": status["recommendation_role"],
                "uses_heuristics": bool(status["uses_heuristics"]),
                "merge_strategy": status["merge_strategy"],
                "algorithm_weight": status["algorithm_weight"],
                "external_evidence": "demuxlet",
                "usable_cells": int(usable.sum()),
                "ambiguous_cells_reported_separately": int((labels == "ambs").sum()),
                "ground_truth_doublets": int(truth.sum()),
                "predicted_doublets": int(predicted.sum()),
                "predicted_rate": float(predicted.mean()),
                "expected_rate_from_demuxlet": expected_rate,
                "score_auc": _auc(truth, scores[method]),
                "heterotypic_risk_proxy": "external_demuxlet_positive",
                "homotypic_risk_proxy": "not_measurable_without_cell-type-pair evidence",
                "runtime_seconds": float(status["runtime_seconds"]),
                "fallback_mode": bool(status["fallback_mode"]),
                "review_required": bool(status["status"] != "ok" or metric["f1"] < 0.5),
                "risk_note": (
                    status["error"]
                    if status["status"] != "ok"
                    else (
                        "Transparent heuristic baseline; use as fallback only."
                        if status["recommendation_role"] == "fallback_heuristic_baseline"
                        else (
                            "Algorithm plus lineage/size heuristic fusion compared against demuxlet labels."
                            if status["uses_heuristics"]
                            else "Algorithmic method compared against demuxlet labels."
                        )
                    )
                ),
                **metric,
            }
        )

    overlap_rows: list[dict[str, Any]] = []
    methods = list(predictions)
    for left in methods:
        for right in methods:
            left_mask = predictions[left].astype(bool)
            right_mask = predictions[right].astype(bool)
            overlap_rows.append(
                {
                    "dataset": spec.key,
                    "left_method": left,
                    "right_method": right,
                    "jaccard": float(
                        (left_mask & right_mask).sum() / max((left_mask | right_mask).sum(), 1)
                    ),
                    "overlap_cells": int((left_mask & right_mask).sum()),
                    "left_cells": int(left_mask.sum()),
                    "right_cells": int(right_mask.sum()),
                }
            )

    label_rows = [
        {
            "dataset": spec.key,
            "demuxlet_label": label,
            "n_cells": int(count),
            "fraction": float(count / max(adata.n_obs, 1)),
        }
        for label, count in labels.value_counts(dropna=False).items()
    ]

    paths = {
        "evidence": output_dir / "doublet_evidence.tsv",
        "overlap": output_dir / "doublet_method_overlap.tsv",
        "labels": output_dir / "doublet_ground_truth_label_counts.tsv",
        "scdblfinder_parity": output_dir / "scdblfinder_python_vs_r_reference.tsv",
        "scdblfinder_parity_by_group": output_dir
        / "scdblfinder_python_vs_r_reference_by_group.tsv",
        "scdblfinder_disagreement_cells": output_dir
        / "scdblfinder_python_vs_r_disagreement_cells.tsv",
        "threshold_calibration": output_dir / "doublet_threshold_calibration.tsv",
        "algorithm_weight_recommendation": output_dir
        / "doublet_algorithm_weight_recommendation.tsv",
        "benchmark_summary": output_dir / "doublet_benchmark_report_summary.json",
    }
    pd.DataFrame(evidence_rows).to_csv(paths["evidence"], sep="\t", index=False)
    pd.DataFrame(overlap_rows).to_csv(paths["overlap"], sep="\t", index=False)
    pd.DataFrame(label_rows).to_csv(paths["labels"], sep="\t", index=False)
    parity_rows = _scdblfinder_python_vs_r_reference_rows(
        dataset=spec.key,
        truth=truth,
        predictions=predictions,
        scores=scores,
        r_reference_path=r_scdblfinder_reference,
    )
    pd.DataFrame(parity_rows).to_csv(paths["scdblfinder_parity"], sep="\t", index=False)
    group_rows, disagreement_rows = _scdblfinder_python_vs_r_reference_details(
        dataset=spec.key,
        obs=adata.obs,
        truth=truth,
        predictions=predictions,
        scores=scores,
        r_reference_path=r_scdblfinder_reference,
        min_group_cells=min_parity_group_cells,
    )
    pd.DataFrame(group_rows).to_csv(paths["scdblfinder_parity_by_group"], sep="\t", index=False)
    pd.DataFrame(disagreement_rows).to_csv(
        paths["scdblfinder_disagreement_cells"], sep="\t", index=False
    )
    calibration_rows = _threshold_calibration_rows(
        dataset=spec.key,
        truth=truth,
        predictions=predictions,
        scores=scores,
    )
    pd.DataFrame(calibration_rows).to_csv(paths["threshold_calibration"], sep="\t", index=False)
    weight_rows = _algorithm_weight_recommendation_rows(evidence_rows)
    pd.DataFrame(weight_rows).to_csv(
        paths["algorithm_weight_recommendation"], sep="\t", index=False
    )
    paths["benchmark_summary"].write_text(
        json.dumps(
            _benchmark_summary_payload(
                evidence_rows=evidence_rows,
                parity_rows=parity_rows,
                group_rows=group_rows,
                calibration_rows=calibration_rows,
                weight_rows=weight_rows,
            ),
            indent=2,
        )
    )
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_outputs/qc_doublet_evidence"),
    )
    parser.add_argument(
        "--max-cells", type=int, default=6000, help="Stratified pilot subset. Use 0 for full data."
    )
    parser.add_argument("--seed", type=int, default=23)
    _DOUBLET_METHOD_CHOICES = (
        "scrublet",
        "solo",
        "doubletdetection",
        "scdblfinder",
        "scdblfinder_python",
        "scdblfinder_python_pyscdblfinder",
    )
    parser.add_argument(
        "--methods",
        nargs="*",
        choices=_DOUBLET_METHOD_CHOICES,
        default=["scrublet", "scdblfinder_python_pyscdblfinder"],
        help=(
            "Algorithmic doublet-detection methods to benchmark. "
            "Heuristic predictions are always emitted as a transparent fallback baseline "
            "and cannot be selected as a primary method."
        ),
    )
    parser.add_argument(
        "--r-scdblfinder-reference",
        type=Path,
        default=None,
        help=(
            "CSV exported from Bioconductor scDblFinder on the same benchmark cells. "
            "Expected columns: cell/barcode, score/scDblFinder.score, "
            "predicted/class/scDblFinder.class."
        ),
    )
    parser.add_argument("--min-parity-group-cells", type=int, default=25)
    parser.add_argument(
        "--algorithm-weights",
        nargs="*",
        type=float,
        default=[0.5, 0.7, 0.85],
        help=(
            "Candidate algorithm weights for algorithm+heuristic weighted-average fusion. "
            "0.7 is the current conservative default."
        ),
    )
    parser.add_argument(
        "--merge-strategy",
        choices=["weighted_average", "max_score", "heuristic_boost"],
        default="weighted_average",
    )
    args = parser.parse_args()
    paths = run(
        args.output_dir,
        max_cells=args.max_cells,
        seed=args.seed,
        methods=args.methods,
        r_scdblfinder_reference=args.r_scdblfinder_reference,
        min_parity_group_cells=args.min_parity_group_cells,
        algorithm_weights=args.algorithm_weights,
        merge_strategy=args.merge_strategy,
    )
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
