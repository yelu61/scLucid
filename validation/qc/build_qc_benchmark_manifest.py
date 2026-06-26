#!/usr/bin/env python3
"""Build the Phase 2 QC benchmark readiness manifest.

The script is deliberately lightweight: it does not run QC. It verifies that
local datasets expose the metadata needed to make QC benchmark claims credible,
then writes tables that downstream benchmark runners can consume.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.dataset_registry import DATASETS


THRESHOLD_DECISION_COLUMNS = (
    "dataset",
    "strategy",
    "parameter",
    "recommended",
    "applied",
    "source",
    "confidence",
    "evidence",
    "review_required",
    "affected_cells",
    "risk_note",
)

QC_STRATEGIES = (
    "scanpy_fixed_threshold",
    "seurat_fixed_threshold",
    "sclucid_adaptive",
    "sclucid_tumor_aware",
)

DOUBLET_STRATEGIES = (
    "scrublet",
    "pyscdblfinder_or_scdblfinder",
    "sclucid_lineage_heuristic",
    "external_demuxlet",
)


def _value_counts(adata: ad.AnnData, column: str) -> dict[str, int]:
    if column not in adata.obs:
        return {}
    counts = adata.obs[column].astype(str).value_counts(dropna=False)
    return {str(k): int(v) for k, v in counts.items()}


def _contract_row(spec, adata: ad.AnnData) -> dict[str, Any]:
    missing_required = [col for col in spec.required_obs if col not in adata.obs]
    missing_annotation = [col for col in spec.annotation_obs if col not in adata.obs]
    has_counts_layer = "counts" in adata.layers
    has_dataset_uns = "sclucid" in adata.uns and "dataset" in adata.uns["sclucid"]
    return {
        "dataset": spec.key,
        "path": str(spec.path),
        "exists": True,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "has_counts_layer": bool(has_counts_layer),
        "var_names_unique": bool(adata.var_names.is_unique),
        "has_dataset_uns": bool(has_dataset_uns),
        "missing_required_obs": ";".join(missing_required),
        "missing_annotation_obs": ";".join(missing_annotation),
        "samples": json.dumps(_value_counts(adata, "sample"), ensure_ascii=False),
        "conditions": json.dumps(_value_counts(adata, "condition"), ensure_ascii=False),
        "cell_types_top": json.dumps(
            dict(list(_value_counts(adata, "cell_type").items())[:10]),
            ensure_ascii=False,
        ),
        "doublet_labels": json.dumps(
            _value_counts(adata, "demuxlet_multiplets"),
            ensure_ascii=False,
        ),
        "ambient_labels": json.dumps(
            _value_counts(adata, "likely_empty_droplet"),
            ensure_ascii=False,
        ),
        "qc_roles": ";".join(spec.qc_roles),
        "figure2_panels": ";".join(spec.figure2_panels),
        "benchmark_notes": spec.benchmark_notes,
    }


def build_manifest(output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    readiness_rows: list[dict[str, Any]] = []
    strategy_rows: list[dict[str, Any]] = []
    figure_rows: list[dict[str, Any]] = []
    issue_rows: list[dict[str, Any]] = []

    for spec in DATASETS:
        if not spec.path.exists():
            readiness_rows.append(
                {
                    "dataset": spec.key,
                    "path": str(spec.path),
                    "exists": False,
                    "benchmark_notes": "Local h5ad is missing.",
                }
            )
            issue_rows.append(
                {
                    "dataset": spec.key,
                    "severity": "blocker",
                    "issue": "missing_h5ad",
                    "detail": str(spec.path),
                }
            )
            continue

        adata = ad.read_h5ad(spec.path, backed="r")
        row = _contract_row(spec, adata)
        readiness_rows.append(row)

        if row["missing_required_obs"]:
            issue_rows.append(
                {
                    "dataset": spec.key,
                    "severity": "warning",
                    "issue": "missing_required_obs",
                    "detail": row["missing_required_obs"],
                }
            )
        if not row["has_counts_layer"]:
            issue_rows.append(
                {
                    "dataset": spec.key,
                    "severity": "blocker",
                    "issue": "missing_counts_layer",
                    "detail": "layers['counts'] is required for QC benchmarks.",
                }
            )

        for strategy in QC_STRATEGIES:
            strategy_rows.append(
                {
                    "dataset": spec.key,
                    "strategy": strategy,
                    "run_group": "qc_threshold",
                    "required_for_phase2": "tumor_aware_qc" in spec.qc_roles
                    or spec.key == "pbmc3k",
                    "notes": spec.benchmark_notes,
                }
            )
        if "doublet_ground_truth" in spec.qc_roles:
            for strategy in DOUBLET_STRATEGIES:
                strategy_rows.append(
                    {
                        "dataset": spec.key,
                        "strategy": strategy,
                        "run_group": "doublet",
                        "required_for_phase2": True,
                        "notes": "Use demuxlet singlet/doublet as external evidence; report ambs separately.",
                    }
                )
        if "ambient_rna" in spec.qc_roles:
            strategy_rows.append(
                {
                    "dataset": spec.key,
                    "strategy": "sclucid_ambient_empty_droplet_diagnostic",
                    "run_group": "ambient",
                    "required_for_phase2": True,
                    "notes": "Diagnostic contract fixture; not a biological ambient correction benchmark.",
                }
            )

        for panel in spec.figure2_panels:
            figure_rows.append(
                {
                    "figure_panel": panel,
                    "dataset": spec.key,
                    "evidence_role": ";".join(spec.qc_roles),
                    "planned_metric_table": {
                        "2A": "qc_workflow_decision_table.tsv",
                        "2B": "retention_marker_fidelity.tsv",
                        "2C": "tumor_program_retention.tsv",
                        "2D": "doublet_ambient_evidence.tsv",
                    }.get(panel, "qc_summary.tsv"),
                }
            )
        adata.file.close()

    paths = {
        "readiness": output_dir / "qc_dataset_readiness.tsv",
        "strategies": output_dir / "qc_strategy_matrix.tsv",
        "threshold_schema": output_dir / "qc_threshold_decision_schema.tsv",
        "figure2_plan": output_dir / "figure2_panel_plan.tsv",
        "issues": output_dir / "qc_metadata_contract_issues.tsv",
    }
    pd.DataFrame(readiness_rows).to_csv(paths["readiness"], sep="\t", index=False)
    pd.DataFrame(strategy_rows).to_csv(paths["strategies"], sep="\t", index=False)
    pd.DataFrame({"column": THRESHOLD_DECISION_COLUMNS}).to_csv(
        paths["threshold_schema"], sep="\t", index=False
    )
    pd.DataFrame(figure_rows).to_csv(paths["figure2_plan"], sep="\t", index=False)
    pd.DataFrame(issue_rows).to_csv(paths["issues"], sep="\t", index=False)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_outputs/qc_manifest"),
    )
    args = parser.parse_args()
    paths = build_manifest(args.output_dir)
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
