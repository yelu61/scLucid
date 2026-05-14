"""Lightweight validation scaffold for QC and preprocessing maturity.

This module intentionally validates engineering and workflow maturity, not
scientific superiority over external pipelines. Comparative benchmarking should
be built on top of this scaffold once the analysis module is equally mature.
"""

from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData

from .contracts import Modules, UnsKeys, validate_review_summary_schema

VALIDATION_SCAFFOLD_SCHEMA_VERSION = "0.1"
VALIDATION_SCOPE = "qc_preprocess_lightweight"
COMPARATIVE_READINESS_LABEL = "ready_for_comparative_validation"


def _now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _as_plain(value: Any) -> Any:
    """Convert numpy/pandas scalar values into JSON-friendly values."""
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    return value


def _get_sclucid_namespace(adata: AnnData) -> dict[str, Any]:
    namespace = adata.uns.get("sclucid", {})
    return namespace if isinstance(namespace, dict) else {}


def _stage_namespace(adata: AnnData, stage: str) -> dict[str, Any]:
    namespace = _get_sclucid_namespace(adata).get(stage, {})
    return namespace if isinstance(namespace, dict) else {}


def _review_summary(adata: AnnData, stage: str) -> dict[str, Any]:
    review = _stage_namespace(adata, stage).get(UnsKeys.REVIEW_SUMMARY, {})
    if isinstance(review, dict) and isinstance(review.get("data"), dict):
        return review["data"]
    return review if isinstance(review, dict) else {}


def _list_like(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, dict) and isinstance(value.get("data"), list):
        return value["data"]
    return [value]


def _stage_warnings(adata: AnnData, stage: str) -> list[Any]:
    stage_ns = _stage_namespace(adata, stage)
    review = _review_summary(adata, stage)
    warnings = []
    warnings.extend(_list_like(review.get("warnings")))
    warnings.extend(_list_like(stage_ns.get("warnings")))
    return warnings


def _fraction_from_bool_column(adata: AnnData, column: str) -> float | None:
    if column not in adata.obs.columns or adata.n_obs == 0:
        return None
    values = adata.obs[column]
    try:
        return float(values.astype(bool).mean())
    except Exception:
        return None


def _hvg_summary(adata: AnnData) -> tuple[int | None, str | None, list[str]]:
    """Return HVG count, canonical/fallback key, and all discovered HVG columns."""
    hvg_columns = [
        col
        for col in adata.var.columns
        if col == "highly_variable" or str(col).startswith("highly_variable_")
    ]
    preferred = [
        "highly_variable",
        "highly_variable_selected",
        "highly_variable_scanpy_seurat_v3",
        "highly_variable_scanpy_seurat",
    ]
    selected_key = next((key for key in preferred if key in adata.var.columns), None)
    if selected_key is None and hvg_columns:
        selected_key = hvg_columns[0]
    if selected_key is None:
        return None, None, hvg_columns
    try:
        return int(adata.var[selected_key].astype(bool).sum()), selected_key, hvg_columns
    except Exception:
        return None, selected_key, hvg_columns


def _review_completeness(adata: AnnData, stage: str) -> dict[str, Any]:
    review = _review_summary(adata, stage)
    if not review:
        return {
            "present": False,
            "valid_schema": False,
            "issues": [f"{stage} review_summary is missing"],
        }
    result = validate_review_summary_schema(review, module=stage, raise_on_error=False)
    issues = list(getattr(result, "errors", []))
    return {
        "present": True,
        "valid_schema": bool(getattr(result, "valid", len(issues) == 0)),
        "issues": issues,
    }


def _status(condition: bool, *, warn_if_false: bool = False) -> str:
    if condition:
        return "pass"
    return "warn" if warn_if_false else "fail"


def _row(metric: str, value: Any, status: str, interpretation: str) -> dict[str, Any]:
    return {
        "metric": metric,
        "value": _as_plain(value),
        "status": status,
        "interpretation": interpretation,
    }


def build_qc_preprocess_validation(
    adata: AnnData,
    *,
    run_manifest: dict[str, Any] | None = None,
    dataset_role: str | None = None,
    workflow_name: str | None = None,
) -> dict[str, Any]:
    """Build a lightweight QC/preprocess validation manifest.

    The manifest answers whether the current run is auditable, reproducible, and
    ready for later comparative validation. It deliberately does not claim that
    scLucid outperforms Scanpy, Seurat, or other standard workflows.
    """
    run_manifest = run_manifest or {}
    input_shape = run_manifest.get("input_shape", {})
    input_cells = input_shape.get("n_cells")
    if input_cells is None:
        input_cells = run_manifest.get("subset_n_cells")
    final_cells = int(adata.n_obs)
    retention_fraction = run_manifest.get("retention_fraction")
    if retention_fraction is None and input_cells:
        retention_fraction = final_cells / max(int(input_cells), 1)

    qc_review = _review_completeness(adata, Modules.QC)
    pp_review = _review_completeness(adata, Modules.PREPROCESS)
    qc_warnings = _stage_warnings(adata, Modules.QC)
    pp_warnings = _stage_warnings(adata, Modules.PREPROCESS)

    hvg_count, hvg_key, hvg_columns = _hvg_summary(adata)
    preprocess_ns = _stage_namespace(adata, Modules.PREPROCESS)
    hvg_stability_available = "hvg_stability" in json.dumps(
        preprocess_ns, default=str
    ).lower()

    layer_contract = {
        "counts_present": "counts" in adata.layers,
        "normalized_present": "normalized" in adata.layers,
        "scaled_present": "scaled" in adata.layers,
        "raw_present": adata.raw is not None,
    }
    representation_contract = {
        "pca_present": "X_pca" in adata.obsm,
        "neighbors_present": "neighbors" in adata.uns,
        "umap_present": "X_umap" in adata.obsm,
    }
    qc_metrics = {
        "retention_fraction": retention_fraction,
        "warning_count": len(qc_warnings),
        "low_quality_fraction": _fraction_from_bool_column(adata, "low_quality"),
        "doublet_fraction": _fraction_from_bool_column(adata, "predicted_doublet"),
    }
    preprocess_metrics = {
        "hvg_count": hvg_count,
        "hvg_key": hvg_key,
        "canonical_hvg_available": "highly_variable" in adata.var.columns,
        "hvg_columns": hvg_columns,
        "hvg_stability_available": hvg_stability_available,
        "layer_contract": layer_contract,
        "representation_contract": representation_contract,
        "warning_count": len(pp_warnings),
    }

    table = [
        _row(
            "cell_retention_fraction",
            retention_fraction,
            _status(retention_fraction is not None and retention_fraction > 0),
            "QC produced a non-empty retained cell set; dataset-specific optimality is not inferred.",
        ),
        _row(
            "qc_warning_count",
            len(qc_warnings),
            "pass",
            "Warnings are counted for audit review, not treated as failure by default.",
        ),
        _row(
            "low_quality_fraction_available",
            qc_metrics["low_quality_fraction"],
            _status(qc_metrics["low_quality_fraction"] is not None, warn_if_false=True),
            "Presence of low-quality calls makes QC filtering auditable.",
        ),
        _row(
            "doublet_fraction_available",
            qc_metrics["doublet_fraction"],
            _status(qc_metrics["doublet_fraction"] is not None, warn_if_false=True),
            "Presence of doublet calls makes doublet handling auditable.",
        ),
        _row(
            "counts_layer_present",
            layer_contract["counts_present"],
            _status(layer_contract["counts_present"]),
            "Raw count preservation is required for reproducible preprocessing and DE choices.",
        ),
        _row(
            "normalized_layer_present",
            layer_contract["normalized_present"],
            _status(layer_contract["normalized_present"]),
            "Normalized layer availability documents the preprocessing handoff.",
        ),
        _row(
            "raw_present",
            layer_contract["raw_present"],
            _status(layer_contract["raw_present"], warn_if_false=True),
            "AnnData.raw is recommended for marker and DE inspection.",
        ),
        _row(
            "hvg_count",
            hvg_count,
            _status(hvg_count is not None and hvg_count > 0),
            "HVG selection produced a discoverable HVG column; canonical highly_variable is preferred.",
        ),
        _row(
            "pca_neighbors_umap_available",
            all(representation_contract.values()),
            _status(all(representation_contract.values())),
            "PCA, graph, and UMAP outputs are present for downstream analysis handoff.",
        ),
        _row(
            "qc_review_summary_valid",
            qc_review["valid_schema"],
            _status(qc_review["valid_schema"]),
            "QC review_summary satisfies the shared scLucid review contract.",
        ),
        _row(
            "preprocess_review_summary_valid",
            pp_review["valid_schema"],
            _status(pp_review["valid_schema"]),
            "Preprocess review_summary satisfies the shared scLucid review contract.",
        ),
    ]
    blocking_failures = [item for item in table if item["status"] == "fail"]
    readiness_status = (
        COMPARATIVE_READINESS_LABEL if not blocking_failures else "needs_scaffold_fix"
    )

    return {
        "schema_version": VALIDATION_SCAFFOLD_SCHEMA_VERSION,
        "scope": VALIDATION_SCOPE,
        "generated_at": _now_iso(),
        "workflow_name": workflow_name or run_manifest.get("workflow") or "unknown",
        "dataset_role": dataset_role or run_manifest.get("dataset_role") or "unknown",
        "claim_boundary": (
            "This scaffold validates auditability, reproducibility, and workflow "
            "maturity only. It does not claim scientific superiority over standard workflows."
        ),
        "readiness_status": readiness_status,
        "ready_for_comparative_validation": readiness_status == COMPARATIVE_READINESS_LABEL,
        "input_shape": input_shape,
        "final_shape": {"n_cells": final_cells, "n_genes": int(adata.n_vars)},
        "qc_metrics": qc_metrics,
        "preprocess_metrics": preprocess_metrics,
        "review_summary_completeness": {
            Modules.QC: qc_review,
            Modules.PREPROCESS: pp_review,
        },
        "compact_validation_table": table,
        "blocking_failures": blocking_failures,
        "next_phase": (
            "After analysis review summaries reach parity with QC/preprocess, "
            "extend this scaffold into qc_preprocess_analysis_validation."
        ),
    }


def write_validation_outputs(
    validation: dict[str, Any],
    output_dir: str | Path,
    *,
    basename: str = "qc_preprocess_validation",
) -> dict[str, str]:
    """Write validation JSON and compact CSV table to ``output_dir``."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    json_path = output_path / f"{basename}.json"
    csv_path = output_path / f"{basename}_table.csv"
    json_path.write_text(json.dumps(validation, indent=2, default=str), encoding="utf-8")

    table = validation.get("compact_validation_table", [])
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["metric", "value", "status", "interpretation"]
        )
        writer.writeheader()
        for row in table:
            writer.writerow(
                {
                    "metric": row.get("metric"),
                    "value": row.get("value"),
                    "status": row.get("status"),
                    "interpretation": row.get("interpretation"),
                }
            )

    return {"json": str(json_path), "table_csv": str(csv_path)}


def validation_table_to_dataframe(validation: dict[str, Any]) -> pd.DataFrame:
    """Return the compact validation table as a DataFrame."""
    return pd.DataFrame(validation.get("compact_validation_table", []))


__all__ = [
    "COMPARATIVE_READINESS_LABEL",
    "VALIDATION_SCAFFOLD_SCHEMA_VERSION",
    "VALIDATION_SCOPE",
    "build_qc_preprocess_validation",
    "validation_table_to_dataframe",
    "write_validation_outputs",
]
