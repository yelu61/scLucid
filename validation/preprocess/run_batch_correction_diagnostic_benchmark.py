#!/usr/bin/env python3
"""Build batch-correction recommendation evidence on real datasets."""

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
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scLucid.preprocess.config import IntegrationConfig
from scLucid.preprocess.integrate import (
    batch_correction,
    diagnose_integration_risk,
    evaluate_integration,
)
from validation.dataset_registry import DATASETS


def _value_counts(series: pd.Series) -> dict[str, int]:
    counts = series.astype(str).value_counts(dropna=False)
    return {str(key): int(value) for key, value in counts.items()}


def _cramers_v(x: pd.Series, y: pd.Series) -> float:
    table = pd.crosstab(x.astype(str), y.astype(str))
    if table.empty or min(table.shape) < 2:
        return 0.0
    observed = table.to_numpy(dtype=float)
    total = observed.sum()
    row_sum = observed.sum(axis=1, keepdims=True)
    col_sum = observed.sum(axis=0, keepdims=True)
    expected = row_sum @ col_sum / max(total, 1.0)
    chi2 = np.divide(
        (observed - expected) ** 2,
        expected,
        out=np.zeros_like(observed),
        where=expected > 0,
    ).sum()
    phi2 = chi2 / max(total, 1.0)
    r, k = observed.shape
    return float(np.sqrt(phi2 / max(min(k - 1, r - 1), 1)))


def _normalized_entropy(counts: pd.Series) -> float:
    values = counts.to_numpy(dtype=float)
    total = values.sum()
    if total <= 0 or len(values) <= 1:
        return 0.0
    p = values / total
    entropy = -np.sum(p * np.log(p + 1e-12))
    return float(entropy / np.log(len(values)))


def _depth_ratio(adata: ad.AnnData, key: str) -> float | None:
    if "total_counts" not in adata.obs or key not in adata.obs:
        return None
    medians = adata.obs.groupby(key, observed=True)["total_counts"].median().dropna()
    if len(medians) < 2:
        return None
    return float(medians.max() / max(float(medians.min()), 1.0))


def _primary_biology_key(adata: ad.AnnData) -> str:
    for key in ("cell_type", "condition", "cell_subtype"):
        if key in adata.obs and adata.obs[key].astype(str).nunique(dropna=False) > 1:
            return key
    return ""


def _candidate_batch_keys(adata: ad.AnnData) -> list[str]:
    return [
        key
        for key in ("donor", "patient", "sample", "sampleID")
        if key in adata.obs and adata.obs[key].astype(str).nunique(dropna=False) > 1
    ]


def _is_tumor_dataset(spec) -> bool:
    return "tumor" in spec.modality_role or any("tumor" in role for role in spec.preprocess_roles)


def _subset(adata: ad.AnnData, max_cells: int | None, seed: int) -> ad.AnnData:
    if max_cells is None or max_cells <= 0 or adata.n_obs <= max_cells:
        return adata.copy()
    rng = np.random.default_rng(seed)
    keep = np.sort(rng.choice(adata.n_obs, size=max_cells, replace=False))
    return adata[keep].copy()


def _method_dependency(method: str) -> str:
    return {
        "harmony": "harmonypy",
        "bbknn": "bbknn",
        "scvi": "scvi",
    }.get(method, "")


def _method_available(method: str) -> tuple[bool, str]:
    package = _method_dependency(method)
    if not package:
        return True, ""
    return bool(importlib.util.find_spec(package)), package


def _prepare_pca(adata: ad.AnnData, n_top_genes: int = 2000, n_pcs: int = 30) -> ad.AnnData:
    work = adata.copy()
    if "counts" in work.layers:
        work.X = work.layers["counts"].copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    if work.n_vars > n_top_genes:
        variance = (
            np.asarray(work.X.toarray()).var(axis=0)
            if hasattr(work.X, "toarray")
            else np.asarray(work.X).var(axis=0)
        )
        keep = np.argsort(-variance, kind="mergesort")[:n_top_genes]
        work = work[:, np.sort(keep)].copy()
    sc.pp.scale(work, max_value=10)
    sc.tl.pca(work, n_comps=min(n_pcs, work.n_obs - 1, work.n_vars - 1), svd_solver="arpack")
    return work


def _comparison_score(evaluation: dict[str, Any]) -> float | None:
    batch = evaluation.get("batch_silhouette")
    label = evaluation.get("label_silhouette")
    if batch is None and label is None:
        return None
    if batch is None:
        return float(label)
    if label is None:
        return float(batch)
    return float((batch + label) / 2.0)


def _run_method_comparison(
    adata: ad.AnnData,
    *,
    dataset: str,
    batch_key: str,
    label_key: str,
    method: str,
    max_epochs: int,
    model_save_dir: Path | None = None,
) -> dict[str, Any]:
    available, package = _method_available(method)
    start = time.perf_counter()
    if not available:
        return {
            "dataset": dataset,
            "batch_key": batch_key,
            "method": method,
            "method_status": "dependency_missing",
            "optional_dependency": package,
            "embedding_key": "",
            "runtime_seconds": 0.0,
            "batch_silhouette": np.nan,
            "label_silhouette": np.nan,
            "graph_connectivity": np.nan,
            "overall_score": np.nan,
            "production_risk_level": "unknown",
            "production_risk_score": np.nan,
            "production_warnings": "",
            "scvi_model_saved": False,
            "scvi_model_path": "",
            "recommendation": "not_run_dependency_missing",
            "risk_note": f"Missing optional dependency: {package}",
        }

    try:
        work = _prepare_pca(adata)
        scvi_model_path = ""
        if method == "no_correction":
            embedding_key = "X_pca"
        elif method == "bbknn":
            cfg = IntegrationConfig(
                method="bbknn", batch_key=batch_key, use_rep="X_pca", evaluate=False, plot=False
            )
            batch_correction(work, config=cfg, force=True)
            embedding_key = "X_pca"
        else:
            output_key = f"X_{method}"
            scvi_params: dict[str, Any] = {
                "n_latent": 15,
                "max_epochs": max_epochs,
            }
            harmony_params: dict[str, Any] = {
                "max_iter_harmony": 50,
                "theta": 2.0,
                "lambda_val": 1.0,
            }
            if method == "scvi" and model_save_dir is not None:
                scvi_model_path = str(
                    model_save_dir / f"{dataset}_{batch_key}_{method}_scvi_model"
                )
                scvi_params.update({
                    "save_model": True,
                    "model_path": scvi_model_path,
                })
            cfg = IntegrationConfig(
                method=method,
                batch_key=batch_key,
                use_rep="X_pca",
                output_key=output_key,
                evaluate=False,
                plot=False,
                scvi_params=scvi_params,
                harmony_params=harmony_params,
            )
            batch_correction(work, config=cfg, force=True)
            embedding_key = output_key

        evaluation = evaluate_integration(
            work,
            batch_key=batch_key,
            label_key=label_key or None,
            integration_method=method,
            use_rep=embedding_key,
            methods=["silhouette", "graph_connectivity"],
            plot=False,
        )
        score = _comparison_score(evaluation)
        recommendation = "candidate"
        risk_note = ""
        if (
            evaluation.get("label_silhouette") is not None
            and evaluation.get("label_silhouette") < 0
        ):
            recommendation = "review_overcorrection"
            risk_note = "Biology label silhouette is negative after this representation; inspect marker fidelity."

        production_risk: dict[str, Any] = {}
        try:
            production_risk = diagnose_integration_risk(
                work,
                batch_key=batch_key,
                condition_key="condition" if "condition" in work.obs else None,
                biology_columns=[label_key] if label_key else [],
                label_key=label_key or None,
                before_rep="X_pca",
                after_rep=embedding_key,
                tumor=False,
            )
        except Exception as exc:
            production_risk = {
                "risk_level": "unknown",
                "risk_score": np.nan,
                "warnings": [f"diagnose_integration_risk failed: {exc}"],
            }

        return {
            "dataset": dataset,
            "batch_key": batch_key,
            "method": method,
            "method_status": "ok",
            "optional_dependency": package,
            "embedding_key": embedding_key,
            "runtime_seconds": time.perf_counter() - start,
            "batch_silhouette": evaluation.get("batch_silhouette", np.nan),
            "label_silhouette": evaluation.get("label_silhouette", np.nan),
            "graph_connectivity": evaluation.get("graph_connectivity", np.nan),
            "overall_score": score,
            "production_risk_level": production_risk.get("risk_level", "unknown"),
            "production_risk_score": production_risk.get("risk_score", np.nan),
            "production_warnings": "; ".join(production_risk.get("warnings", [])),
            "scvi_model_saved": bool(
                work.uns.get("sclucid", {})
                .get("preprocess", {})
                .get("integration", {})
                .get("scvi", {})
                .get("model_saved", False)
            ),
            "scvi_model_path": scvi_model_path,
            "recommendation": recommendation,
            "risk_note": risk_note,
        }
    except Exception as exc:
        return {
            "dataset": dataset,
            "batch_key": batch_key,
            "method": method,
            "method_status": "failed",
            "optional_dependency": package,
            "embedding_key": "",
            "runtime_seconds": time.perf_counter() - start,
            "batch_silhouette": np.nan,
            "label_silhouette": np.nan,
            "graph_connectivity": np.nan,
            "overall_score": np.nan,
            "production_risk_level": "unknown",
            "production_risk_score": np.nan,
            "production_warnings": "",
            "scvi_model_saved": False,
            "scvi_model_path": "",
            "recommendation": "not_run_failed",
            "risk_note": f"{type(exc).__name__}: {exc}",
        }


def _recommendation(
    *,
    n_groups: int,
    biology_association: float,
    condition_association: float,
    depth_ratio: float | None,
    is_tumor: bool,
) -> tuple[str, str, bool]:
    if n_groups < 2:
        return "no_correction", "Only one batch group is present.", False
    if biology_association >= 0.7:
        return (
            "diagnostic_only",
            "Batch key is strongly associated with biological labels; correction risks removing true structure.",
            True,
        )
    if condition_association >= 0.7:
        return (
            "diagnostic_only",
            "Batch key is strongly associated with biological condition; compare embeddings before any correction.",
            True,
        )
    if depth_ratio is not None and depth_ratio >= 2.0:
        return (
            "correction_candidate_review_depth",
            "Batch groups have large sequencing-depth differences; correction may help but depth/biology should be reviewed.",
            True,
        )
    return (
        "correction_candidate",
        "Multiple batch groups with limited observed biology confounding; opt-in correction can be benchmarked.",
        False,
    )


def run(
    output_dir: Path,
    datasets: set[str] | None,
    methods: list[str],
    max_cells: int | None,
    seed: int,
    scvi_max_epochs: int,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, Any]] = []
    mixing_rows: list[dict[str, Any]] = []
    risk_rows: list[dict[str, Any]] = []
    figure_rows: list[dict[str, Any]] = []
    method_rows: list[dict[str, Any]] = []

    model_save_dir = output_dir / "scvi_models"
    model_save_dir.mkdir(parents=True, exist_ok=True)

    for spec in DATASETS:
        if datasets and spec.key not in datasets:
            continue
        if not spec.preprocess_roles or not spec.path.exists():
            continue
        adata = ad.read_h5ad(spec.path, backed="r")
        biology_key = _primary_biology_key(adata)
        condition_present = "condition" in adata.obs
        is_tumor = _is_tumor_dataset(spec)
        candidate_keys = _candidate_batch_keys(adata)
        for batch_key in candidate_keys:
            n_groups = int(adata.obs[batch_key].astype(str).nunique(dropna=False))
            biology_assoc = (
                _cramers_v(adata.obs[batch_key], adata.obs[biology_key]) if biology_key else 0.0
            )
            condition_assoc = (
                _cramers_v(adata.obs[batch_key], adata.obs["condition"])
                if condition_present
                else 0.0
            )
            depth = _depth_ratio(adata, batch_key)
            recommendation, rationale, review_required = _recommendation(
                n_groups=n_groups,
                biology_association=biology_assoc,
                condition_association=condition_assoc,
                depth_ratio=depth,
                is_tumor=is_tumor,
            )
            group_counts = _value_counts(adata.obs[batch_key])
            group_sizes = pd.Series(group_counts)
            balance = _normalized_entropy(group_sizes)

            # Also run production diagnose_integration_risk on the full dataset for comparison.
            production_risk: dict[str, Any] = {}
            try:
                production_risk = diagnose_integration_risk(
                    adata,
                    batch_key=batch_key,
                    condition_key="condition" if condition_present else None,
                    biology_columns=[biology_key] if biology_key else [],
                    label_key=biology_key or None,
                    before_rep="X_pca" if "X_pca" in adata.obsm else None,
                    after_rep=None,
                    tumor=is_tumor,
                )
            except Exception as exc:
                production_risk = {
                    "risk_level": "unknown",
                    "risk_score": np.nan,
                    "warnings": [f"diagnose_integration_risk failed: {exc}"],
                }

            row = {
                "dataset": spec.key,
                "batch_key": batch_key,
                "biology_key": biology_key,
                "n_batch_groups": n_groups,
                "batch_balance_entropy": balance,
                "batch_biology_cramers_v": biology_assoc,
                "batch_condition_cramers_v": condition_assoc,
                "depth_median_ratio": depth,
                "recommendation": recommendation,
                "production_risk_level": production_risk.get("risk_level", "unknown"),
                "production_risk_score": production_risk.get("risk_score", np.nan),
                "production_warnings": "; ".join(production_risk.get("warnings", [])),
                "review_required": review_required,
                "rationale": rationale,
                "group_sizes": json.dumps(group_counts, ensure_ascii=False),
            }
            summary_rows.append(row)
            mixing_rows.append(
                {
                    "dataset": spec.key,
                    "batch_key": batch_key,
                    "comparison": "no_correction_vs_opt_in_correction",
                    "batch_mixing_need_proxy": 1.0 - balance,
                    "biological_conservation_risk_proxy": max(
                        biology_assoc, condition_assoc if is_tumor else 0.0
                    ),
                    "recommended_next_step": recommendation,
                }
            )
            risk_rows.append(
                {
                    "dataset": spec.key,
                    "batch_key": batch_key,
                    "risk_type": "overcorrection",
                    "risk_score": max(biology_assoc, condition_assoc if is_tumor else 0.0),
                    "risk_level": (
                        "high"
                        if review_required
                        else "moderate" if recommendation != "correction_candidate" else "low"
                    ),
                    "risk_note": rationale,
                }
            )
            figure_rows.append(
                {
                    "figure_panel": "3C",
                    "dataset": spec.key,
                    "batch_key": batch_key,
                    "metric": "batch_biology_cramers_v",
                    "value": biology_assoc,
                    "context": json.dumps(
                        {
                            "recommendation": recommendation,
                            "condition_association": condition_assoc,
                            "depth_ratio": depth,
                        }
                    ),
                }
            )
        if candidate_keys:
            adata_full = ad.read_h5ad(spec.path)
            adata_subset = _subset(adata_full, max_cells=max_cells, seed=seed)
            label_key = _primary_biology_key(adata_subset)
            for batch_key in candidate_keys[:2]:
                for method in methods:
                    method_rows.append(
                        _run_method_comparison(
                            adata_subset,
                            dataset=spec.key,
                            batch_key=batch_key,
                            label_key=label_key,
                            method=method,
                            max_epochs=scvi_max_epochs,
                            model_save_dir=model_save_dir,
                        )
                    )
        adata.file.close()

    paths = {
        "summary": output_dir / "batch_diagnostic_summary.tsv",
        "mixing": output_dir / "batch_mixing_vs_biology.tsv",
        "risk": output_dir / "overcorrection_risk.tsv",
        "method_comparison": output_dir / "batch_method_comparison.tsv",
        "figure3": output_dir / "figure3_batch_data.tsv",
    }
    pd.DataFrame(summary_rows).to_csv(paths["summary"], sep="\t", index=False)
    pd.DataFrame(mixing_rows).to_csv(paths["mixing"], sep="\t", index=False)
    pd.DataFrame(risk_rows).to_csv(paths["risk"], sep="\t", index=False)
    pd.DataFrame(method_rows).to_csv(paths["method_comparison"], sep="\t", index=False)
    pd.DataFrame(figure_rows).to_csv(paths["figure3"], sep="\t", index=False)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("validation_outputs/work/preprocess_batch_diagnostic")
    )
    parser.add_argument("--datasets", nargs="*", help="Optional dataset keys to include.")
    parser.add_argument(
        "--methods", nargs="*", default=["no_correction", "harmony", "bbknn", "scvi"]
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=2000,
        help="Pilot subset size for method comparisons. Use 0 for full data.",
    )
    parser.add_argument("--seed", type=int, default=19)
    parser.add_argument("--scvi-max-epochs", type=int, default=20)
    args = parser.parse_args()
    paths = run(
        args.output_dir,
        datasets=set(args.datasets) if args.datasets else None,
        methods=args.methods,
        max_cells=args.max_cells,
        seed=args.seed,
        scvi_max_epochs=args.scvi_max_epochs,
    )
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
