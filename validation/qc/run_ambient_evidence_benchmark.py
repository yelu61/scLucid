#!/usr/bin/env python3
"""Run a lightweight ambient/empty-droplet evidence benchmark.

The benchmark always validates the contract on ``cellbender_tiny``. If the user
provides a raw 10x directory via ``--raw-10x-dir``, it also runs diagnostics and
ambient correction on a real full matrix and attempts CellBender/scAR/SoupX-like
backends when available. Real-data results are still marked as preliminary until
calibrated against validated ground-truth empty droplets.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.dataset_registry import DATASETS


def _counts_matrix(adata: ad.AnnData):
    return adata.layers["counts"] if "counts" in adata.layers else adata.X


def _ensure_metrics(adata: ad.AnnData) -> None:
    X = _counts_matrix(adata)
    if "barcode_total_counts" not in adata.obs:
        adata.obs["barcode_total_counts"] = np.asarray(X.sum(axis=1)).ravel()
    if "n_genes_by_counts" not in adata.obs:
        if sp.issparse(X):
            adata.obs["n_genes_by_counts"] = np.asarray((X > 0).sum(axis=1)).ravel()
        else:
            adata.obs["n_genes_by_counts"] = (np.asarray(X) > 0).sum(axis=1)


def _rank_summary(values: pd.Series, mask: pd.Series) -> dict[str, Any]:
    selected = values[mask]
    if selected.empty:
        return {"median": np.nan, "q05": np.nan, "q95": np.nan}
    return {
        "median": float(selected.median()),
        "q05": float(selected.quantile(0.05)),
        "q95": float(selected.quantile(0.95)),
    }


def _available(backend: str) -> bool:
    if backend == "cellbender":
        return bool(importlib.util.find_spec("cellbender")) or bool(
            importlib.util.find_spec("shutil")
        )
    if backend == "soupx":
        return bool(importlib.util.find_spec("soupx"))
    if backend == "scar":
        return bool(importlib.util.find_spec("scar"))
    return False


def _run_real_raw_evidence(
    raw_10x_dir: Path,
    max_barcodes: int | None,
    seed: int,
) -> dict[str, Any]:
    """Load a raw 10x matrix and run ambient diagnostics + correction."""
    from scLucid.qc.ambient import (
        correct_ambient_rna_linear,
        diagnose_ambient_rna,
        diagnose_empty_droplets,
    )
    from scLucid.qc.ambient_backends import correct_ambient_rna
    from scLucid.utils import read_10x

    adata = read_10x(raw_10x_dir, var_names="gene_symbols", make_unique=True)
    _ensure_metrics(adata)
    if max_barcodes is not None and max_barcodes > 0 and adata.n_obs > max_barcodes:
        rng = np.random.default_rng(seed)
        keep = np.sort(rng.choice(adata.n_obs, size=max_barcodes, replace=False))
        adata = adata[keep].copy()
        _ensure_metrics(adata)

    empty_diag = diagnose_empty_droplets(adata)
    ambient_diag = diagnose_ambient_rna(adata)

    correction: dict[str, Any] = {}
    try:
        correction = correct_ambient_rna(
            adata,
            method="auto",
            backend="auto",
            output_layer="ambient_corrected",
        )
    except Exception as exc:
        correction = {
            "corrected": False,
            "backend": "auto_failed",
            "error": str(exc),
            "risk_note": "Auto backend failed; falling back to linear correction.",
        }
        try:
            linear = correct_ambient_rna_linear(
                adata, output_layer="ambient_corrected"
            )
            correction.update(linear)
            correction["backend_fallback"] = "linear"
        except Exception as linear_exc:
            correction["linear_error"] = str(linear_exc)

    return {
        "dataset": f"raw_10x_{raw_10x_dir.name}",
        "raw_10x_dir": str(raw_10x_dir),
        "n_barcodes": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "empty_droplet_available": empty_diag.get("available", False),
        "empty_droplet_risk_level": empty_diag.get("risk_level", "unknown"),
        "empty_droplet_risk_score": empty_diag.get("risk_score", np.nan),
        "n_putative_empty_droplets": empty_diag.get("n_putative_empty_droplets", 0),
        "ambient_risk_level": ambient_diag.get("risk_level", "unknown"),
        "ambient_risk_score": ambient_diag.get("risk_score", np.nan),
        "corrected": correction.get("corrected", False),
        "correction_backend": correction.get("backend", correction.get("backend_fallback", "")),
        "removed_counts": correction.get("removed_counts", np.nan),
        "removed_fraction": correction.get("removed_fraction", np.nan),
        "mean_rho": correction.get("mean_rho", np.nan),
        "review_required": bool(
            empty_diag.get("risk_level") in {"moderate", "high"}
            or ambient_diag.get("risk_level") in {"moderate", "high"}
            or correction.get("review_required", False)
        ),
        "risk_note": (
            "Real raw 10x ambient evidence is preliminary until compared against validated "
            "empty-droplet calls or CellBender/scAR/SoupX reference output."
        ),
        "cellbender_available": _available("cellbender"),
        "soupx_available": _available("soupx"),
        "scar_available": _available("scar"),
    }


def run(
    output_dir: Path,
    raw_10x_dir: Path | None = None,
    max_barcodes: int | None = None,
    seed: int = 42,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    spec = next(dataset for dataset in DATASETS if dataset.key == "cellbender_tiny")
    adata = ad.read_h5ad(spec.path)
    _ensure_metrics(adata)

    likely_cell = adata.obs["likely_cell"].astype(bool)
    likely_empty = adata.obs["likely_empty_droplet"].astype(bool)
    counts = adata.obs["barcode_total_counts"].astype(float)
    genes = adata.obs["n_genes_by_counts"].astype(float)
    rank = (
        adata.obs["barcode_rank"].astype(float)
        if "barcode_rank" in adata.obs
        else counts.rank(method="first", ascending=False).astype(float)
    )

    evidence_rows = [
        {
            "dataset": spec.key,
            "diagnostic": "empty_droplet_contract",
            "n_barcodes": int(adata.n_obs),
            "likely_cells": int(likely_cell.sum()),
            "likely_empty_droplets": int(likely_empty.sum()),
            "counts_layer_present": bool("counts" in adata.layers),
            "median_counts_likely_cell": float(counts[likely_cell].median()),
            "median_counts_likely_empty": float(counts[likely_empty].median()),
            "cell_to_empty_median_count_ratio": float(
                counts[likely_cell].median() / max(counts[likely_empty].median(), 1.0)
            ),
            "median_genes_likely_cell": float(genes[likely_cell].median()),
            "median_genes_likely_empty": float(genes[likely_empty].median()),
            "review_required": False,
            "risk_note": (
                "Tiny fixture validates diagnostic contract only; use a full raw matrix for performance claims."
            ),
        }
    ]

    if raw_10x_dir is not None and raw_10x_dir.exists():
        evidence_rows.append(_run_real_raw_evidence(raw_10x_dir, max_barcodes, seed))

    rank_rows = []
    for label, mask in (
        ("likely_cell", likely_cell),
        ("likely_empty_droplet", likely_empty),
    ):
        count_summary = _rank_summary(counts, mask)
        gene_summary = _rank_summary(genes, mask)
        rank_summary = _rank_summary(rank, mask)
        rank_rows.append(
            {
                "dataset": spec.key,
                "barcode_class": label,
                "n_barcodes": int(mask.sum()),
                "count_summary": json.dumps(count_summary),
                "gene_summary": json.dumps(gene_summary),
                "rank_summary": json.dumps(rank_summary),
            }
        )

    # Lightweight correction residual assessment
    residual_rows = []
    try:
        from scLucid.qc.ambient import correct_ambient_rna_linear

        correction = correct_ambient_rna_linear(
            adata, output_layer="ambient_corrected", empty_droplet_key="likely_empty_droplet"
        )
        residual_rows.append(
            {
                "dataset": spec.key,
                "method": "linear_background_subtraction",
                "corrected": correction.get("corrected"),
                "removed_counts": correction.get("removed_counts"),
                "removed_fraction": correction.get("removed_fraction"),
                "mean_rho": correction.get("mean_rho"),
                "residual_ambient_score": correction.get("residual_ambient_score"),
                "review_required": correction.get("review_required"),
            }
        )
    except Exception as exc:
        residual_rows.append(
            {
                "dataset": spec.key,
                "method": "linear_background_subtraction",
                "corrected": False,
                "error": str(exc),
            }
        )

    paths = {
        "evidence": output_dir / "ambient_evidence.tsv",
        "barcode_classes": output_dir / "ambient_barcode_class_summary.tsv",
        "correction_residual": output_dir / "ambient_correction_residual.tsv",
    }
    pd.DataFrame(evidence_rows).to_csv(paths["evidence"], sep="\t", index=False)
    pd.DataFrame(rank_rows).to_csv(paths["barcode_classes"], sep="\t", index=False)
    pd.DataFrame(residual_rows).to_csv(paths["correction_residual"], sep="\t", index=False)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_outputs/qc_ambient_evidence"),
    )
    parser.add_argument(
        "--raw-10x-dir",
        type=Path,
        default=None,
        help="Optional path to a Cell Ranger raw_feature_bc_matrix directory for real-data evidence.",
    )
    parser.add_argument(
        "--max-barcodes",
        type=int,
        default=None,
        help="Optional cap on barcodes to load from the raw 10x matrix.",
    )
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    paths = run(
        args.output_dir,
        raw_10x_dir=args.raw_10x_dir,
        max_barcodes=args.max_barcodes,
        seed=args.seed,
    )
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
