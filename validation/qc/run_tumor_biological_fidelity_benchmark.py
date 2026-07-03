#!/usr/bin/env python3
"""Validate tumor-aware QC against biological marker/program preservation."""

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

from scLucid.qc.policy.intelligent_qc import (
    IntelligentQCConfig,
    IntelligentQCRecommender,
    StrategyType,
)
from validation.dataset_registry import DATASETS
from validation.gene_panels import MARKER_PANELS, TUMOR_PROGRAM_PANELS, present_genes

STRATEGIES = (
    "scanpy_fixed_threshold",
    "seurat_fixed_threshold",
    "sclucid_adaptive",
    "sclucid_tumor_aware",
)


def _matrix(adata: ad.AnnData):
    return adata.layers["counts"] if "counts" in adata.layers else adata.X


def _subset(adata: ad.AnnData, max_cells: int | None, seed: int) -> ad.AnnData:
    if max_cells is None or max_cells <= 0 or adata.n_obs <= max_cells:
        return adata
    rng = np.random.default_rng(seed)
    keep = np.sort(rng.choice(adata.n_obs, size=max_cells, replace=False))
    return adata[keep].copy()


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


def _adaptive_thresholds(adata: ad.AnnData, *, tumor_aware: bool) -> dict[str, Any]:
    """Use scLucid's production IntelligentQCRecommender for adaptive thresholds."""
    strategy = StrategyType.TUMOR_AWARE if tumor_aware else StrategyType.STANDARD
    recommender = IntelligentQCRecommender(strategy=strategy)
    rec = recommender.recommend(
        adata,
        tissue_type="tumor" if tumor_aware else "normal",
        plot=False,
    )
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


def _thresholds(strategy: str, adata: ad.AnnData) -> dict[str, Any]:
    if strategy == "scanpy_fixed_threshold":
        return {"min_genes": 200, "min_counts": 0, "max_mt_percent": 20.0}
    if strategy == "seurat_fixed_threshold":
        return {"min_genes": 200, "min_counts": 0, "max_mt_percent": 5.0}
    if strategy == "sclucid_adaptive":
        return _adaptive_thresholds(adata, tumor_aware=False)
    if strategy == "sclucid_tumor_aware":
        return _adaptive_thresholds(adata, tumor_aware=True)
    raise ValueError(strategy)


def _keep_mask(
    adata: ad.AnnData, thresholds: dict[str, Any]
) -> tuple[pd.Series, dict[str, pd.Series]]:
    obs = adata.obs
    fail = {
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


def _malignant_like_retention(
    adata: ad.AnnData, keep: pd.Series, genes: tuple[str, ...]
) -> dict[str, Any]:
    """Compute retention of cells with a malignant-like expression signature.

    Cells whose average expression of the malignant-like panel is above the
    dataset median are treated as a putative malignant-enriched population.
    The metric reports what fraction of these cells survive QC, independent of
    any author-provided malignant label.
    """
    present = present_genes(adata.var_names, genes)
    if len(present) < 2:
        return {
            "malignant_like_cells": 0,
            "malignant_like_retention_rate": np.nan,
            "malignant_like_genes_present": 0,
        }
    X = _matrix(adata)
    idx = pd.Index(adata.var_names.astype(str)).get_indexer(present)
    mean_score = np.asarray(X[:, idx].mean(axis=1)).ravel()
    median_score = float(np.median(mean_score))
    malignant_like = mean_score > median_score
    retained = keep.to_numpy()
    n_malignant_like = int(malignant_like.sum())
    n_retained = int((malignant_like & retained).sum())
    return {
        "malignant_like_cells": n_malignant_like,
        "malignant_like_retention_rate": (
            n_retained / max(n_malignant_like, 1) if n_malignant_like > 0 else np.nan
        ),
        "malignant_like_genes_present": len(present),
    }


def _mean_expression(
    adata: ad.AnnData, genes: list[str], mask: pd.Series | np.ndarray | None = None
) -> float:
    if not genes:
        return float("nan")
    X = _matrix(adata)
    idx = pd.Index(adata.var_names.astype(str)).get_indexer(genes)
    if mask is not None:
        arr_mask = np.asarray(mask, dtype=bool)
        if arr_mask.sum() == 0:
            return float("nan")
        X = X[arr_mask, :]
    values = np.asarray(X[:, idx].mean(axis=0)).ravel()
    return float(np.nanmean(values)) if values.size else float("nan")


def _retention_bias_rows(
    adata: ad.AnnData, dataset: str, strategy: str, keep: pd.Series
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for key, group_type in (
        ("sample", "sample"),
        ("patient", "patient"),
        ("condition", "condition"),
        ("cell_type", "cell_type"),
    ):
        if key not in adata.obs:
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
                    "review_required": False,
                    "review_reason": f"missing_{key}",
                }
            )
            continue
        for group, idx in adata.obs.groupby(key, observed=True).indices.items():
            group_mask = pd.Series(False, index=adata.obs_names)
            group_mask.iloc[list(idx)] = True
            before = int(group_mask.sum())
            after = int((group_mask & keep).sum())
            retention_rate = after / max(before, 1)
            # Distinguish "the group itself is tiny" from "the strategy removed most of it".
            if before < 25:
                review_required = False
                review_reason = "small_group_skipped"
            elif retention_rate < 0.5:
                review_required = True
                review_reason = "strategy_bias_low_retention"
            else:
                review_required = False
                review_reason = ""
            rows.append(
                {
                    "dataset": dataset,
                    "strategy": strategy,
                    "group_type": group_type,
                    "group": str(group),
                    "before": before,
                    "after": after,
                    "removed": before - after,
                    "retention_rate": retention_rate,
                    "review_required": review_required,
                    "review_reason": review_reason,
                }
            )
    return rows


def _strategy_scorecard_rows(
    marker_rows: list[dict[str, Any]],
    program_rows: list[dict[str, Any]],
    bias_rows: list[dict[str, Any]],
    malignant_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    markers = pd.DataFrame(marker_rows)
    programs = pd.DataFrame(program_rows)
    bias = pd.DataFrame(bias_rows)
    malignant = pd.DataFrame(malignant_rows)
    if programs.empty:
        return []
    rows: list[dict[str, Any]] = []
    for (dataset, strategy), program_group in programs.groupby(
        ["dataset", "strategy"], observed=True
    ):
        marker_group = markers[(markers["dataset"] == dataset) & (markers["strategy"] == strategy)]
        bias_group = bias[(bias["dataset"] == dataset) & (bias["strategy"] == strategy)]
        malignant_group = malignant[
            (malignant["dataset"] == dataset) & (malignant["strategy"] == strategy)
        ]
        high_mt_cells = int(program_group["high_mt_removed_cells"].max())
        program_retention = pd.to_numeric(program_group["program_retention_ratio"], errors="coerce")
        high_mt_program = pd.to_numeric(
            program_group["high_mt_removed_program_ratio"], errors="coerce"
        )
        epithelial_marker = marker_group[
            marker_group["marker_panel"].isin(["epithelial", "proliferation", "hypoxia_stress"])
        ]
        high_mt_marker = pd.to_numeric(
            epithelial_marker["high_mt_removed_relative_to_all"], errors="coerce"
        )
        eligible_bias = bias_group[bias_group["before"] >= 25]
        min_group_retention = (
            float(eligible_bias["retention_rate"].min()) if not eligible_bias.empty else np.nan
        )
        mean_program_retention = (
            float(program_retention.mean()) if program_retention.notna().any() else np.nan
        )
        high_mt_biology_signal = (
            float(high_mt_program.max()) if high_mt_program.notna().any() else np.nan
        )
        high_mt_marker_signal = (
            float(high_mt_marker.max()) if high_mt_marker.notna().any() else np.nan
        )
        malignant_like_retention = pd.to_numeric(
            malignant_group.get("malignant_like_retention_rate", pd.Series(dtype=float)),
            errors="coerce",
        )
        malignant_retention_value = (
            float(malignant_like_retention.iloc[0])
            if not malignant_like_retention.empty and malignant_like_retention.notna().any()
            else np.nan
        )
        biological_harm_risk = bool(
            high_mt_cells > 0
            and (
                (pd.notna(high_mt_biology_signal) and high_mt_biology_signal >= 1.0)
                or (pd.notna(high_mt_marker_signal) and high_mt_marker_signal >= 1.0)
            )
        )
        biological_fidelity_score = np.nanmean(
            [
                min(mean_program_retention, 1.0) if pd.notna(mean_program_retention) else np.nan,
                min(min_group_retention, 1.0) if pd.notna(min_group_retention) else np.nan,
                min(malignant_retention_value, 1.0) if pd.notna(malignant_retention_value) else np.nan,
                0.0 if biological_harm_risk else 1.0,
            ]
        )
        review_reasons: list[str] = []
        if biological_harm_risk:
            review_reasons.append("biological_harm_risk")
        if pd.notna(min_group_retention) and min_group_retention < 0.5:
            review_reasons.append("stratified_retention_low")
        if pd.notna(mean_program_retention) and mean_program_retention < 0.8:
            review_reasons.append("program_retention_low")
        if pd.notna(malignant_retention_value) and malignant_retention_value < 0.5:
            review_reasons.append("malignant_like_retention_low")
        rows.append(
            {
                "dataset": dataset,
                "strategy": strategy,
                "min_group_retention": min_group_retention,
                "mean_program_retention_ratio": mean_program_retention,
                "malignant_like_retention_rate": malignant_retention_value,
                "high_mt_removed_cells": high_mt_cells,
                "high_mt_max_program_ratio": high_mt_biology_signal,
                "high_mt_max_marker_ratio": high_mt_marker_signal,
                "biological_harm_risk": biological_harm_risk,
                "biological_fidelity_score": float(biological_fidelity_score),
                "review_required": bool(review_reasons),
                "review_reasons": ";".join(review_reasons),
                "risk_note": (
                    "High-mt removed cells carry tumor/stress/proliferation signal; avoid mechanical deletion."
                    if biological_harm_risk
                    else (
                        "Stratified retention, program retention, or malignant-like retention is low; inspect tumor-aware QC."
                        if review_reasons
                        else ""
                    )
                ),
            }
        )
    score_df = pd.DataFrame(rows)
    score_df["rank_within_dataset"] = (
        score_df.groupby("dataset")["biological_fidelity_score"]
        .rank(method="first", ascending=False)
        .astype(int)
    )
    score_df["recommended_for_review"] = score_df["rank_within_dataset"] == 1
    return score_df.sort_values(["dataset", "rank_within_dataset"]).to_dict("records")


def _narrative_rows(
    marker_rows: list[dict[str, Any]],
    program_rows: list[dict[str, Any]],
    bias_rows: list[dict[str, Any]],
    scorecard_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Build reviewer-facing tumor-aware QC evidence narratives."""
    markers = pd.DataFrame(marker_rows)
    programs = pd.DataFrame(program_rows)
    bias = pd.DataFrame(bias_rows)
    if not scorecard_rows:
        return []
    scorecard = pd.DataFrame(scorecard_rows)
    rows: list[dict[str, Any]] = []
    for _, score in scorecard.iterrows():
        dataset = score["dataset"]
        strategy = score["strategy"]
        marker_group = markers[(markers["dataset"] == dataset) & (markers["strategy"] == strategy)]
        program_group = programs[
            (programs["dataset"] == dataset) & (programs["strategy"] == strategy)
        ]
        bias_group = bias[(bias["dataset"] == dataset) & (bias["strategy"] == strategy)]

        protected_marker_group = marker_group[
            marker_group["marker_panel"].isin(["epithelial", "proliferation", "hypoxia_stress"])
        ]
        high_mt_marker = pd.to_numeric(
            protected_marker_group["high_mt_removed_relative_to_all"],
            errors="coerce",
        )
        high_mt_program = pd.to_numeric(
            program_group["high_mt_removed_program_ratio"],
            errors="coerce",
        )
        program_retention = pd.to_numeric(
            program_group["program_retention_ratio"],
            errors="coerce",
        )
        eligible_bias = bias_group[bias_group["before"] >= 25]
        worst_bias = (
            eligible_bias.sort_values("retention_rate").iloc[0].to_dict()
            if not eligible_bias.empty
            else {}
        )
        high_mt_signal_panels = sorted(
            protected_marker_group.loc[
                pd.to_numeric(
                    protected_marker_group["high_mt_removed_relative_to_all"],
                    errors="coerce",
                )
                >= 1.0,
                "marker_panel",
            ]
            .astype(str)
            .unique()
            .tolist()
        )
        high_mt_signal_programs = sorted(
            program_group.loc[
                pd.to_numeric(program_group["high_mt_removed_program_ratio"], errors="coerce")
                >= 1.0,
                "program",
            ]
            .astype(str)
            .unique()
            .tolist()
        )
        recommended = bool(score.get("recommended_for_review", False))
        rank = int(score["rank_within_dataset"])
        harm_risk = bool(score.get("biological_harm_risk", False))
        low_retention = bool(score.get("review_reasons", ""))
        malignant_low = (
            pd.notna(score.get("malignant_like_retention_rate"))
            and score.get("malignant_like_retention_rate") < 0.5
        )
        if recommended:
            narrative = (
                f"Rank {rank} policy by tumor biological fidelity; recommended for review "
                "because it best balances stratified retention, marker/program preservation, "
                "malignant-like cell retention, and high-mt biological-signal risk in this dataset."
            )
        else:
            narrative = (
                f"Rank {rank} policy by tumor biological fidelity; compare against the "
                "recommended policy before using this threshold set."
            )
        if harm_risk:
            narrative += (
                " High-mt removed cells retain tumor/stress/proliferation signal, so mechanical "
                "deletion should be avoided or explicitly justified."
            )
        if low_retention:
            narrative += (
                " One or more retention dimensions are below target: "
                f"{score.get('review_reasons', '')}."
            )
        if malignant_low:
            narrative += (
                " Malignant-like cells show low retention; verify epithelial/malignant marker preservation."
            )

        rows.append(
            {
                "dataset": dataset,
                "strategy": strategy,
                "strategy_rank": rank,
                "recommended_policy": recommended,
                "biological_fidelity_score": score.get("biological_fidelity_score"),
                "mean_program_retention_ratio": (
                    float(program_retention.mean()) if program_retention.notna().any() else np.nan
                ),
                "max_high_mt_program_ratio": (
                    float(high_mt_program.max()) if high_mt_program.notna().any() else np.nan
                ),
                "max_high_mt_marker_ratio": (
                    float(high_mt_marker.max()) if high_mt_marker.notna().any() else np.nan
                ),
                "high_mt_removed_cells": score.get("high_mt_removed_cells"),
                "high_mt_signal_programs": ";".join(high_mt_signal_programs),
                "high_mt_signal_marker_panels": ";".join(high_mt_signal_panels),
                "worst_group_type": worst_bias.get("group_type"),
                "worst_group": worst_bias.get("group"),
                "worst_group_retention": worst_bias.get("retention_rate"),
                "biological_harm_risk": harm_risk,
                "review_required": bool(score.get("review_required", False)),
                "risk_note": score.get("risk_note", ""),
                "decision_narrative": narrative,
            }
        )
    return rows


def run(
    output_dir: Path, datasets: set[str] | None, max_cells: int | None, seed: int
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    marker_rows: list[dict[str, Any]] = []
    program_rows: list[dict[str, Any]] = []
    bias_rows: list[dict[str, Any]] = []
    malignant_rows: list[dict[str, Any]] = []
    figure_rows: list[dict[str, Any]] = []

    tumor_specs = [
        spec
        for spec in DATASETS
        if spec.path.exists()
        and (datasets is None or spec.key in datasets)
        and ("tumor" in spec.modality_role or any("tumor" in role for role in spec.qc_roles))
    ]
    for spec in tumor_specs:
        adata = _subset(ad.read_h5ad(spec.path), max_cells=max_cells, seed=seed)
        _ensure_qc_metrics(adata)
        for strategy in STRATEGIES:
            thresholds = _thresholds(strategy, adata)
            keep, fail = _keep_mask(adata, thresholds)
            mt_removed = fail["max_mt_percent"]
            bias_rows.extend(_retention_bias_rows(adata, spec.key, strategy, keep))
            malignant_genes = TUMOR_PROGRAM_PANELS["epithelial_malignant_like"]
            malignant_ret = _malignant_like_retention(adata, keep, malignant_genes)
            malignant_rows.append(
                {
                    "dataset": spec.key,
                    "strategy": strategy,
                    **malignant_ret,
                }
            )
            for panel, genes in MARKER_PANELS.items():
                present = present_genes(adata.var_names, genes)
                before = _mean_expression(adata, present)
                kept = _mean_expression(adata, present, keep)
                high_mt_removed = _mean_expression(adata, present, mt_removed)
                marker_rows.append(
                    {
                        "dataset": spec.key,
                        "strategy": strategy,
                        "marker_panel": panel,
                        "genes_present": len(present),
                        "genes": ";".join(present),
                        "mean_expression_before": before,
                        "mean_expression_kept": kept,
                        "mean_expression_high_mt_removed": high_mt_removed,
                        "kept_relative_change": (
                            kept / before - 1.0 if before and before > 0 else np.nan
                        ),
                        "high_mt_removed_relative_to_all": (
                            high_mt_removed / before if before and before > 0 else np.nan
                        ),
                        "max_mt_threshold": thresholds.get("max_mt_percent"),
                        "high_mt_removed_cells": int(mt_removed.sum()),
                        "review_required": bool(
                            panel in {"epithelial", "proliferation", "hypoxia_stress"}
                            and int(mt_removed.sum()) > 0
                        ),
                    }
                )
            for program, genes in TUMOR_PROGRAM_PANELS.items():
                present = present_genes(adata.var_names, genes)
                before = _mean_expression(adata, present)
                kept = _mean_expression(adata, present, keep)
                high_mt_removed = _mean_expression(adata, present, mt_removed)
                row = {
                    "dataset": spec.key,
                    "strategy": strategy,
                    "program": program,
                    "genes_present": len(present),
                    "genes": ";".join(present),
                    "program_mean_before": before,
                    "program_mean_kept": kept,
                    "program_mean_high_mt_removed": high_mt_removed,
                    "program_retention_ratio": kept / before if before and before > 0 else np.nan,
                    "high_mt_removed_program_ratio": (
                        high_mt_removed / before if before and before > 0 else np.nan
                    ),
                    "high_mt_removed_cells": int(mt_removed.sum()),
                }
                row["review_required"] = bool(
                    row["high_mt_removed_cells"] > 0
                    and pd.notna(row["high_mt_removed_program_ratio"])
                    and row["high_mt_removed_program_ratio"] >= 1.0
                )
                program_rows.append(row)
                figure_rows.append(
                    {
                        "figure_panel": "2C",
                        "dataset": spec.key,
                        "strategy": strategy,
                        "metric": f"{program}_retention_ratio",
                        "value": row["program_retention_ratio"],
                        "context": json.dumps(
                            {
                                "genes_present": len(present),
                                "high_mt_removed_cells": int(mt_removed.sum()),
                            }
                        ),
                    }
                )

    paths = {
        "marker_retention": output_dir / "tumor_marker_retention.tsv",
        "program_retention": output_dir / "tumor_program_retention.tsv",
        "retention_bias": output_dir / "sample_celltype_retention_bias.tsv",
        "malignant_like_retention": output_dir / "malignant_like_retention.tsv",
        "scorecard": output_dir / "tumor_qc_strategy_scorecard.tsv",
        "narrative": output_dir / "tumor_qc_biological_fidelity_narrative.tsv",
        "figure2": output_dir / "figure2_tumor_fidelity_data.tsv",
    }
    pd.DataFrame(marker_rows).to_csv(paths["marker_retention"], sep="\t", index=False)
    pd.DataFrame(program_rows).to_csv(paths["program_retention"], sep="\t", index=False)
    pd.DataFrame(bias_rows).to_csv(paths["retention_bias"], sep="\t", index=False)
    pd.DataFrame(malignant_rows).to_csv(paths["malignant_like_retention"], sep="\t", index=False)
    scorecard_rows = _strategy_scorecard_rows(marker_rows, program_rows, bias_rows, malignant_rows)
    pd.DataFrame(scorecard_rows).to_csv(paths["scorecard"], sep="\t", index=False)
    narrative_rows = _narrative_rows(marker_rows, program_rows, bias_rows, scorecard_rows)
    pd.DataFrame(narrative_rows).to_csv(paths["narrative"], sep="\t", index=False)
    for row in scorecard_rows:
        figure_rows.append(
            {
                "figure_panel": "2D",
                "dataset": row["dataset"],
                "strategy": row["strategy"],
                "metric": "biological_fidelity_score",
                "value": row["biological_fidelity_score"],
                "context": json.dumps(
                    {
                        "rank_within_dataset": row["rank_within_dataset"],
                        "biological_harm_risk": row["biological_harm_risk"],
                        "high_mt_removed_cells": row["high_mt_removed_cells"],
                        "recommended_policy": row["recommended_for_review"],
                    }
                ),
            }
        )
    for row in narrative_rows:
        figure_rows.append(
            {
                "figure_panel": "2E",
                "dataset": row["dataset"],
                "strategy": row["strategy"],
                "metric": "tumor_qc_review_required",
                "value": int(bool(row["review_required"])),
                "context": json.dumps(
                    {
                        "recommended_policy": row["recommended_policy"],
                        "biological_harm_risk": row["biological_harm_risk"],
                        "worst_group_type": row["worst_group_type"],
                        "worst_group": row["worst_group"],
                        "decision_narrative": row["decision_narrative"],
                    }
                ),
            }
        )
    pd.DataFrame(figure_rows).to_csv(paths["figure2"], sep="\t", index=False)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("validation_outputs/qc_tumor_fidelity")
    )
    parser.add_argument("--datasets", nargs="*", help="Optional dataset keys to include.")
    parser.add_argument(
        "--max-cells",
        type=int,
        default=0,
        help=(
            "Subset size per dataset. Default is 0 (full data) so the benchmark reflects "
            "real retention behavior. Use a positive value only for pilot runs."
        ),
    )
    parser.add_argument("--seed", type=int, default=13)
    args = parser.parse_args()
    paths = run(
        args.output_dir,
        datasets=set(args.datasets) if args.datasets else None,
        max_cells=args.max_cells,
        seed=args.seed,
    )
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
