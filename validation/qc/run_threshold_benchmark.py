#!/usr/bin/env python3
"""Run a lightweight real-data QC threshold benchmark.

This is the first executable Phase 2 runner. It intentionally focuses on
threshold evidence and retention effects; doublet-method comparisons and
ambient-specific diagnostics should live in separate runners.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scLucid.qc.adaptive_threshold import fit_count_mixture_threshold_model
from scLucid.qc.intelligent_qc import (
    IntelligentQCConfig,
    IntelligentQCRecommender,
    StrategyType,
)
from validation.dataset_registry import DATASETS

MARKER_PANELS: dict[str, tuple[str, ...]] = {
    "immune_t": ("CD3D", "CD3E", "TRAC", "CD4", "CD8A", "NKG7"),
    "myeloid": ("LYZ", "S100A8", "S100A9", "FCGR3A", "MS4A7", "LST1"),
    "b_plasma": ("MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"),
    "epithelial": ("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"),
    "stromal": ("COL1A1", "COL1A2", "DCN", "LUM", "ACTA2"),
    "endothelial": ("PECAM1", "VWF", "KDR", "ENG"),
    "proliferation": ("MKI67", "TOP2A", "STMN1", "UBE2C"),
    "hypoxia_stress": ("VEGFA", "CA9", "DDIT3", "HSPA1A", "JUN"),
}


def _subset(adata: ad.AnnData, max_cells: int | None, seed: int) -> ad.AnnData:
    if max_cells is None or max_cells <= 0 or adata.n_obs <= max_cells:
        return adata
    rng = np.random.default_rng(seed)
    keep = np.sort(rng.choice(adata.n_obs, size=max_cells, replace=False))
    return adata[keep].copy()


def _matrix(adata: ad.AnnData):
    return adata.layers["counts"] if "counts" in adata.layers else adata.X


def _ensure_qc_metrics(adata: ad.AnnData) -> None:
    X = _matrix(adata)
    if "total_counts" not in adata.obs:
        adata.obs["total_counts"] = np.asarray(X.sum(axis=1)).ravel()
    if "n_genes_by_counts" not in adata.obs:
        if sp.issparse(X):
            adata.obs["n_genes_by_counts"] = np.asarray((X > 0).sum(axis=1)).ravel()
        else:
            adata.obs["n_genes_by_counts"] = (np.asarray(X) > 0).sum(axis=1)
    if "pct_counts_mt" not in adata.obs:
        names = pd.Index(adata.var_names.astype(str))
        mt_mask = names.str.startswith("MT-") | names.str.startswith("mt-")
        if mt_mask.any():
            mt_counts = np.asarray(X[:, mt_mask].sum(axis=1)).ravel()
            total = np.asarray(adata.obs["total_counts"], dtype=float)
            adata.obs["pct_counts_mt"] = np.divide(
                mt_counts * 100.0,
                total,
                out=np.zeros_like(total, dtype=float),
                where=total > 0,
            )
        else:
            adata.obs["pct_counts_mt"] = np.nan


def _mad(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return float("nan")
    med = np.median(finite)
    return float(np.median(np.abs(finite - med)))


def _adaptive_thresholds(adata: ad.AnnData, tumor_aware: bool) -> dict[str, Any]:
    """Use scLucid's IntelligentQCRecommender to obtain real adaptive thresholds.

    This delegates to the production recommendation code rather than re-implementing
    a simplified MAD heuristic, so the benchmark compares the actual scLucid pipeline.
    """
    strategy = StrategyType.TUMOR_AWARE if tumor_aware else StrategyType.STANDARD
    recommender = IntelligentQCRecommender(strategy=strategy)
    rec = recommender.recommend(adata, tissue_type="tumor" if tumor_aware else "normal", plot=False)
    return {
        "min_genes": int(rec.min_genes.threshold),
        "min_counts": int(rec.n_counts.threshold),
        "max_mt_percent": float(rec.max_mt_percent.threshold)
        if rec.max_mt_percent.threshold is not None
        else None,
        "recommendation_method": rec.max_mt_percent.method,
        "overall_strategy": rec.overall_strategy.value,
        "overall_confidence": rec.overall_confidence,
    }


def _count_mixture_thresholds(adata: ad.AnnData, tumor_aware: bool) -> dict[str, Any]:
    """Adaptive thresholds with explicit count-mixture modelling for n_genes.

    Uses the production :func:`fit_count_mixture_threshold_model` for min_genes
    and the production :class:`IntelligentQCRecommender` for the remaining
    thresholds, mirroring how a user would combine model-driven recommendations.
    """
    n_genes = np.asarray(adata.obs["n_genes_by_counts"], dtype=float)
    count_fit = fit_count_mixture_threshold_model(
        n_genes,
        direction="lower",
        percentile=10.0,
        model="auto",
        random_state=42,
        fallback=True,
    )
    # Cap count-mixture min_genes so it cannot exceed a conservative MAD bound.
    nmads = 4.0 if tumor_aware else 3.0
    med = float(np.nanmedian(n_genes))
    mad = float(_mad(n_genes))
    conservative_max = max(50.0 if tumor_aware else 100.0, med - nmads * mad)
    min_genes = min(conservative_max, float(count_fit["threshold"]))
    min_genes = max(50.0 if tumor_aware else 100.0, min_genes)

    strategy = StrategyType.TUMOR_AWARE if tumor_aware else StrategyType.STANDARD
    cfg = IntelligentQCConfig(min_genes_model="percentile")
    recommender = IntelligentQCRecommender(strategy=strategy, config=cfg)
    rec = recommender.recommend(adata, tissue_type="tumor" if tumor_aware else "normal", plot=False)

    return {
        "min_genes": int(min_genes),
        "min_counts": int(rec.n_counts.threshold),
        "max_mt_percent": float(rec.max_mt_percent.threshold)
        if rec.max_mt_percent.threshold is not None
        else None,
        "count_model": count_fit.get("model"),
        "count_model_aic": count_fit.get("aic"),
        "count_model_fallback": count_fit.get("fallback_used"),
        "recommendation_method": rec.max_mt_percent.method,
        "overall_strategy": rec.overall_strategy.value,
        "overall_confidence": rec.overall_confidence,
    }


def _strategy_thresholds(strategy: str, adata: ad.AnnData, is_tumor: bool) -> dict[str, Any]:
    if strategy == "scanpy_fixed_threshold":
        return {"min_genes": 200, "min_counts": 0, "max_mt_percent": 20.0}
    if strategy == "seurat_fixed_threshold":
        return {"min_genes": 200, "min_counts": 0, "max_mt_percent": 5.0}
    if strategy == "sclucid_adaptive":
        return _adaptive_thresholds(adata, tumor_aware=False)
    if strategy == "sclucid_tumor_aware":
        thresholds = _adaptive_thresholds(adata, tumor_aware=is_tumor)
        if not is_tumor:
            thresholds["max_mt_percent"] = min(
                thresholds["max_mt_percent"] or 20.0,
                20.0,
            )
        return thresholds
    if strategy == "sclucid_count_adaptive":
        thresholds = _count_mixture_thresholds(adata, tumor_aware=is_tumor)
        if not is_tumor:
            thresholds["max_mt_percent"] = min(
                thresholds["max_mt_percent"] or 20.0,
                20.0,
            )
        return thresholds
    raise ValueError(strategy)


def _apply_thresholds(
    adata: ad.AnnData, thresholds: dict[str, Any]
) -> tuple[pd.Series, dict[str, pd.Series]]:
    obs = adata.obs
    fail: dict[str, pd.Series] = {
        "min_genes": obs["n_genes_by_counts"].astype(float) < float(thresholds["min_genes"]),
        "min_counts": obs["total_counts"].astype(float) < float(thresholds["min_counts"]),
    }
    if thresholds.get("max_mt_percent") is not None and obs["pct_counts_mt"].notna().any():
        fail["max_mt_percent"] = obs["pct_counts_mt"].astype(float) > float(
            thresholds["max_mt_percent"]
        )
    else:
        fail["max_mt_percent"] = pd.Series(False, index=obs.index)
    remove = fail["min_genes"] | fail["min_counts"] | fail["max_mt_percent"]
    return ~remove, fail


def _decision_rows(
    dataset: str,
    strategy: str,
    thresholds: dict[str, Any],
    fail: dict[str, pd.Series],
    n_cells: int,
    is_tumor: bool,
) -> list[dict[str, Any]]:
    rows = []
    if strategy in {"scanpy_fixed_threshold", "seurat_fixed_threshold"}:
        source = "fixed_threshold_baseline"
    elif strategy == "sclucid_tumor_aware":
        source = "sclucid_tumor_aware"
    elif strategy == "sclucid_count_adaptive":
        source = "sclucid_count_adaptive"
    else:
        source = "sclucid_adaptive"
    for param, applied in thresholds.items():
        if param in {"count_model", "count_model_aic", "count_model_fallback"}:
            continue
        affected = int(fail[param].sum()) if param in fail else 0
        review_required = bool(is_tumor and param == "max_mt_percent" and affected > 0)
        biological_guardrail = ""
        if is_tumor and param == "max_mt_percent":
            biological_guardrail = "preserve high-mt malignant/stress/program signal until reviewed"
        evidence = {
            "affected_fraction": affected / max(n_cells, 1),
            "n_cells": n_cells,
        }
        if param == "min_genes" and "count_model" in thresholds:
            evidence["count_model"] = thresholds.get("count_model")
            evidence["count_model_aic"] = thresholds.get("count_model_aic")
            evidence["count_model_fallback"] = thresholds.get("count_model_fallback")
        rows.append(
            {
                "dataset": dataset,
                "strategy": strategy,
                "parameter": param,
                "recommended": applied,
                "applied": applied,
                "source": source,
                "confidence": "baseline" if "fixed" in strategy else "medium",
                "evidence": json.dumps(evidence),
                "review_required": review_required,
                "affected_cells": affected,
                "biological_guardrail": biological_guardrail,
                "risk_note": (
                    "Tumor dataset: high mitochondrial cells may include stressed or malignant biology; inspect marker retention."
                    if review_required
                    else ""
                ),
            }
        )
    return rows


def _retention_rows(
    adata: ad.AnnData, dataset: str, strategy: str, keep: pd.Series
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = [
        {
            "dataset": dataset,
            "strategy": strategy,
            "group_type": "all",
            "group": "all",
            "before": int(adata.n_obs),
            "after": int(keep.sum()),
            "removed": int((~keep).sum()),
            "retention_rate": float(keep.mean()),
            "annotation_status": "ok",
        }
    ]
    for col, group_type in (
        ("sample", "sample"),
        ("condition", "condition"),
        ("cell_type", "cell_type"),
    ):
        if col not in adata.obs:
            rows.append(
                {
                    "dataset": dataset,
                    "strategy": strategy,
                    "group_type": group_type,
                    "group": "__missing__",
                    "before": 0,
                    "after": 0,
                    "removed": 0,
                    "retention_rate": float("nan"),
                    "annotation_status": f"missing_{col}",
                }
            )
            continue
        for group, idx in adata.obs.groupby(col, observed=True).indices.items():
            mask = pd.Series(False, index=adata.obs_names)
            mask.iloc[list(idx)] = True
            before = int(mask.sum())
            after = int((mask & keep).sum())
            rows.append(
                {
                    "dataset": dataset,
                    "strategy": strategy,
                    "group_type": group_type,
                    "group": str(group),
                    "before": before,
                    "after": after,
                    "removed": before - after,
                    "retention_rate": after / max(before, 1),
                    "annotation_status": "ok",
                }
            )
    return rows


def _marker_rows(
    adata: ad.AnnData, dataset: str, strategy: str, keep: pd.Series
) -> list[dict[str, Any]]:
    X = _matrix(adata)
    var_names = pd.Index(adata.var_names.astype(str))
    rows: list[dict[str, Any]] = []
    keep_arr = keep.to_numpy()
    for panel, genes in MARKER_PANELS.items():
        present = [g for g in genes if g in var_names]
        if not present:
            continue
        idx = var_names.get_indexer(present)
        before = np.asarray(X[:, idx].mean(axis=0)).ravel()
        after = (
            np.asarray(X[keep_arr, :][:, idx].mean(axis=0)).ravel()
            if keep_arr.any()
            else np.zeros(len(idx))
        )
        before_mean = float(np.mean(before))
        after_mean = float(np.mean(after))
        rows.append(
            {
                "dataset": dataset,
                "strategy": strategy,
                "marker_panel": panel,
                "genes_detected": len(present),
                "genes": ";".join(present),
                "mean_expression_before": before_mean,
                "mean_expression_after": after_mean,
                "relative_change": (after_mean / before_mean - 1.0) if before_mean > 0 else np.nan,
            }
        )
    return rows


def _strategy_scorecard_rows(
    *,
    retention_rows: list[dict[str, Any]],
    marker_rows: list[dict[str, Any]],
    decision_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    retention = pd.DataFrame(retention_rows)
    markers = pd.DataFrame(marker_rows)
    decisions = pd.DataFrame(decision_rows)
    if retention.empty:
        return []

    score_rows: list[dict[str, Any]] = []
    for (dataset, strategy), group in retention.groupby(["dataset", "strategy"], observed=True):
        all_ret = group[group["group_type"] == "all"]["retention_rate"]
        overall_retention = float(all_ret.iloc[0]) if not all_ret.empty else np.nan
        sample_ret = group[(group["group_type"] == "sample") & (group["before"] >= 25)][
            "retention_rate"
        ]
        celltype_ret = group[(group["group_type"] == "cell_type") & (group["before"] >= 25)][
            "retention_rate"
        ]
        min_sample_retention = (
            float(sample_ret.min()) if not sample_ret.empty else overall_retention
        )
        min_celltype_retention = (
            float(celltype_ret.min()) if not celltype_ret.empty else overall_retention
        )

        marker_group = markers[(markers["dataset"] == dataset) & (markers["strategy"] == strategy)]
        marker_drift = pd.to_numeric(
            marker_group.get("relative_change", pd.Series(dtype=float)), errors="coerce"
        ).abs()
        median_marker_abs_change = (
            float(marker_drift.median()) if marker_drift.notna().any() else np.nan
        )
        marker_fidelity_score = (
            1.0 - min(float(median_marker_abs_change), 1.0)
            if pd.notna(median_marker_abs_change)
            else np.nan
        )

        decision_group = decisions[
            (decisions["dataset"] == dataset) & (decisions["strategy"] == strategy)
        ]
        mt_row = decision_group[decision_group["parameter"] == "max_mt_percent"]
        mt_removed_fraction = (
            float(mt_row["affected_cells"].iloc[0])
            / max(int(json.loads(mt_row["evidence"].iloc[0])["n_cells"]), 1)
            if not mt_row.empty and mt_row["evidence"].iloc[0]
            else 0.0
        )
        is_tumor = bool(mt_row["review_required"].any()) if not mt_row.empty else False
        retention_fairness = float(np.nanmin([min_sample_retention, min_celltype_retention]))
        tumor_safety_score = 1.0 - min(mt_removed_fraction * (2.0 if is_tumor else 1.0), 1.0)
        components = {
            "overall_retention": overall_retention,
            "retention_fairness": retention_fairness,
            "marker_fidelity_score": marker_fidelity_score,
            "tumor_safety_score": tumor_safety_score,
        }
        weights = {
            "overall_retention": 0.30,
            "retention_fairness": 0.25,
            "marker_fidelity_score": 0.25,
            "tumor_safety_score": 0.20,
        }
        valid = {key: value for key, value in components.items() if pd.notna(value)}
        denom = sum(weights[key] for key in valid)
        composite = sum(valid[key] * weights[key] for key in valid) / max(denom, 1e-12)
        score_rows.append(
            {
                "dataset": dataset,
                "strategy": strategy,
                "overall_retention": overall_retention,
                "min_sample_retention": min_sample_retention,
                "min_celltype_retention": min_celltype_retention,
                "retention_fairness": retention_fairness,
                "median_marker_abs_change": median_marker_abs_change,
                "marker_fidelity_score": marker_fidelity_score,
                "mt_removed_fraction": mt_removed_fraction,
                "tumor_safety_score": tumor_safety_score,
                "composite_score": composite,
                "review_required": bool(
                    retention_fairness < 0.5
                    or (pd.notna(marker_fidelity_score) and marker_fidelity_score < 0.5)
                    or tumor_safety_score < 0.5
                ),
                "risk_note": (
                    "Low stratified retention or marker fidelity; inspect threshold decision before applying."
                    if retention_fairness < 0.5
                    or (pd.notna(marker_fidelity_score) and marker_fidelity_score < 0.5)
                    else (
                        "High mt-removal pressure in tumor-aware context; inspect biological fidelity."
                        if tumor_safety_score < 0.5
                        else ""
                    )
                ),
            }
        )

    score_df = pd.DataFrame(score_rows)
    if score_df.empty:
        return score_rows
    score_df["rank_within_dataset"] = (
        score_df.groupby("dataset")["composite_score"]
        .rank(method="first", ascending=False)
        .astype(int)
    )
    score_df["recommended_for_review"] = score_df["rank_within_dataset"] == 1
    return score_df.sort_values(["dataset", "rank_within_dataset"]).to_dict("records")


def _annotate_decision_rows_with_strategy_evidence(
    decision_rows: list[dict[str, Any]],
    scorecard_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Attach policy-level scorecard evidence to per-threshold decisions."""
    score_by_key = {
        (row["dataset"], row["strategy"]): row
        for row in scorecard_rows
        if row.get("dataset") is not None and row.get("strategy") is not None
    }
    annotated: list[dict[str, Any]] = []
    for row in decision_rows:
        item = dict(row)
        score = score_by_key.get((item.get("dataset"), item.get("strategy")), {})
        rank = score.get("rank_within_dataset")
        recommended = bool(score.get("recommended_for_review", False))
        item["strategy_rank"] = rank
        item["recommended_policy"] = recommended
        item["strategy_composite_score"] = score.get("composite_score")
        item["strategy_risk_note"] = score.get("risk_note", "")
        item["decision_narrative"] = (
            f"Rank {rank} policy for this dataset by retention fairness, marker fidelity, "
            "and tumor-safety score."
            if rank is not None
            else "Policy scorecard unavailable for this dataset."
        )
        if recommended:
            item["decision_narrative"] += " Recommended as the benchmark-selected QC policy."
        annotated.append(item)
    return annotated


def run(
    output_dir: Path, datasets: set[str] | None, max_cells: int | None, seed: int
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    decision_rows: list[dict[str, Any]] = []
    retention_rows: list[dict[str, Any]] = []
    marker_rows: list[dict[str, Any]] = []

    for spec in DATASETS:
        if datasets and spec.key not in datasets:
            continue
        if not spec.path.exists() or spec.key == "cellbender_tiny":
            continue
        adata = ad.read_h5ad(spec.path)
        adata = _subset(adata, max_cells=max_cells, seed=seed)
        _ensure_qc_metrics(adata)
        is_tumor = "tumor" in spec.modality_role or any("tumor" in role for role in spec.qc_roles)
        for strategy in (
            "scanpy_fixed_threshold",
            "seurat_fixed_threshold",
            "sclucid_adaptive",
            "sclucid_tumor_aware",
            "sclucid_count_adaptive",
        ):
            thresholds = _strategy_thresholds(strategy, adata, is_tumor=is_tumor)
            keep, fail = _apply_thresholds(adata, thresholds)
            decision_rows.extend(
                _decision_rows(spec.key, strategy, thresholds, fail, adata.n_obs, is_tumor)
            )
            retention_rows.extend(_retention_rows(adata, spec.key, strategy, keep))
            marker_rows.extend(_marker_rows(adata, spec.key, strategy, keep))

    paths = {
        "decision_table": output_dir / "qc_threshold_decision_table.tsv",
        "retention": output_dir / "qc_retention_summary.tsv",
        "marker_fidelity": output_dir / "qc_marker_fidelity.tsv",
        "scorecard": output_dir / "qc_strategy_scorecard.tsv",
        "figure2": output_dir / "figure2_threshold_data.tsv",
    }
    scorecard_rows = _strategy_scorecard_rows(
        retention_rows=retention_rows,
        marker_rows=marker_rows,
        decision_rows=decision_rows,
    )
    decision_rows = _annotate_decision_rows_with_strategy_evidence(
        decision_rows,
        scorecard_rows,
    )
    pd.DataFrame(decision_rows).to_csv(paths["decision_table"], sep="\t", index=False)
    pd.DataFrame(retention_rows).to_csv(paths["retention"], sep="\t", index=False)
    pd.DataFrame(marker_rows).to_csv(paths["marker_fidelity"], sep="\t", index=False)
    pd.DataFrame(scorecard_rows).to_csv(paths["scorecard"], sep="\t", index=False)
    figure_rows = []
    for row in scorecard_rows:
        for metric in (
            "overall_retention",
            "retention_fairness",
            "marker_fidelity_score",
            "tumor_safety_score",
            "composite_score",
        ):
            figure_rows.append(
                {
                    "figure_panel": "2B",
                    "dataset": row["dataset"],
                    "strategy": row["strategy"],
                    "metric": metric,
                    "value": row[metric],
                    "context": json.dumps(
                        {
                            "rank_within_dataset": row["rank_within_dataset"],
                            "review_required": row["review_required"],
                        }
                    ),
                }
            )
    pd.DataFrame(figure_rows).to_csv(paths["figure2"], sep="\t", index=False)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("validation_outputs/qc_threshold_benchmark")
    )
    parser.add_argument("--datasets", nargs="*", default=None)
    parser.add_argument(
        "--max-cells",
        type=int,
        default=5000,
        help="Deterministic subset size per dataset for pilot runs. Use 0 for full data.",
    )
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    paths = run(
        args.output_dir,
        datasets=set(args.datasets) if args.datasets else None,
        max_cells=None if args.max_cells == 0 else args.max_cells,
        seed=args.seed,
    )
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
