#!/usr/bin/env python3
"""Benchmark annotation accuracy against author-provided labels.

This runner treats author labels as an external reference (not a gold standard)
and reports major-lineage agreement, confusion matrices, and annotation
confidence calibration. It is intentionally lightweight: it clusters the data,
runs scLucid's marker-based annotation, and compares to the reference.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scLucid.analysis.annotation.workflow import run_annotation
from scLucid.analysis.config import AnnotationConfig
from validation.dataset_registry import DATASETS


def _major_lineage(label: str) -> str:
    """Map a detailed label to a major lineage for lenient agreement."""
    label = str(label).lower()
    if any(token in label for token in ("t cell", "t-cell", "cd4", "cd8", "nk")):
        return "lymphoid"
    if any(token in label for token in ("b cell", "b-cell", "plasma")):
        return "lymphoid"
    if any(token in label for token in ("myeloid", "mono", "macrophage", "dc ", "dendritic")):
        return "myeloid"
    if any(token in label for token in ("epithelial", "malignant", "tumor")):
        return "epithelial"
    if any(token in label for token in ("fibroblast", "stromal", "caf", "endothelial")):
        return "stromal"
    return "other"


def _prepare(adata: ad.AnnData, n_hvg: int = 2000, seed: int = 42) -> ad.AnnData:
    work = adata.copy()
    if "counts" in work.layers:
        work.X = work.layers["counts"].copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    sc.pp.highly_variable_genes(work, n_top_genes=min(n_hvg, work.n_vars - 1), flavor="seurat")
    work.raw = work
    sc.pp.scale(work, max_value=10)
    sc.tl.pca(work, n_comps=min(30, work.n_obs - 1, work.n_vars - 1), svd_solver="arpack")
    sc.pp.neighbors(work, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(work, resolution=1.0, random_state=seed, flavor="leidenalg")
    return work


def _accuracy_rows(reference: pd.Series, predicted: pd.Series) -> list[dict[str, Any]]:
    ref_major = reference.astype(str).apply(_major_lineage)
    pred_major = predicted.astype(str).apply(_major_lineage)
    exact = (reference.astype(str) == predicted.astype(str)).mean()
    major = (ref_major == pred_major).mean()
    return [
        {
            "metric": "exact_label_accuracy",
            "value": float(exact),
        },
        {
            "metric": "major_lineage_accuracy",
            "value": float(major),
        },
        {
            "metric": "n_reference_labels",
            "value": int(reference.nunique()),
        },
        {
            "metric": "n_predicted_labels",
            "value": int(predicted.nunique()),
        },
    ]


def _confusion_rows(reference: pd.Series, predicted: pd.Series) -> list[dict[str, Any]]:
    table = pd.crosstab(reference.astype(str), predicted.astype(str))
    rows = []
    for ref_label in table.index:
        best_pred = table.loc[ref_label].idxmax()
        best_count = int(table.loc[ref_label].max())
        total = int(table.loc[ref_label].sum())
        rows.append(
            {
                "reference_label": ref_label,
                "predicted_label": best_pred,
                "n_cells": total,
                "best_match_count": best_count,
                "best_match_fraction": best_count / max(total, 1),
            }
        )
    return rows


def run(
    output_dir: Path,
    datasets: set[str] | None,
    max_cells: int | None,
    seed: int,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    accuracy_rows: list[dict[str, Any]] = []
    confusion_rows: list[dict[str, Any]] = []
    confidence_rows: list[dict[str, Any]] = []

    for spec in DATASETS:
        if datasets and spec.key not in datasets:
            continue
        if not spec.path.exists():
            continue
        if "cell_type" not in spec.annotation_obs:
            continue
        adata = ad.read_h5ad(spec.path)
        if "cell_type" not in adata.obs:
            continue
        if max_cells is not None and max_cells > 0 and adata.n_obs > max_cells:
            rng = np.random.default_rng(seed)
            keep = np.sort(rng.choice(adata.n_obs, size=max_cells, replace=False))
            adata = adata[keep].copy()

        work = _prepare(adata, seed=seed)
        config = AnnotationConfig(
            cluster_key="leiden",
            key_added="sclucid_annotation",
            final_method="max_score",
            marker_species="human",
            marker_tissue=spec.tissue.lower(),
            plot=False,
        )
        work = run_annotation(work, config=config)

        ref = work.obs["cell_type"]
        pred = work.obs["sclucid_annotation"]
        for row in _accuracy_rows(ref, pred):
            accuracy_rows.append({"dataset": spec.key, **row})
        for row in _confusion_rows(ref, pred):
            confusion_rows.append({"dataset": spec.key, **row})

        evidence = (
            work.uns.get("sclucid", {})
            .get("analysis", {})
            .get("annotation", {})
            .get("sclucid_annotation_evidence", pd.DataFrame())
        )
        if not evidence.empty and "annotation_confidence" in evidence.columns:
            confidence_rows.append(
                {
                    "dataset": spec.key,
                    "mean_annotation_confidence": float(
                        pd.to_numeric(evidence["annotation_confidence"], errors="coerce").mean()
                    ),
                    "n_clusters": int(evidence.shape[0]),
                    "review_required": bool(
                        (evidence["annotation_confidence"] < 0.5).any()
                        if "annotation_confidence" in evidence.columns
                        else False
                    ),
                }
            )

    paths = {
        "accuracy": output_dir / "annotation_accuracy.tsv",
        "confusion": output_dir / "annotation_confusion.tsv",
        "confidence": output_dir / "annotation_confidence.tsv",
    }
    pd.DataFrame(accuracy_rows).to_csv(paths["accuracy"], sep="\t", index=False)
    pd.DataFrame(confusion_rows).to_csv(paths["confusion"], sep="\t", index=False)
    pd.DataFrame(confidence_rows).to_csv(paths["confidence"], sep="\t", index=False)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("validation_outputs/analysis_annotation_accuracy")
    )
    parser.add_argument("--datasets", nargs="*", help="Optional dataset keys to include.")
    parser.add_argument(
        "--max-cells", type=int, default=3000, help="Pilot subset size. Use 0 for full data."
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
