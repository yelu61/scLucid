#!/usr/bin/env python3
"""Build the Phase 3 preprocess/analysis benchmark readiness manifest."""

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


LAYER_TRANSITION_COLUMNS = (
    "dataset",
    "workflow",
    "step",
    "input_layer",
    "output_slot",
    "adata_X_semantics_before",
    "adata_X_semantics_after",
    "raw_semantics",
    "parameters",
    "evidence_key",
    "review_required",
)

PREPROCESS_COMPARISONS = (
    "sclucid_standard_no_integration",
    "sclucid_standard_with_diagnostic_only",
    "sclucid_harmony_opt_in",
    "scanpy_standard",
)


def _value_counts(adata: ad.AnnData, column: str) -> dict[str, int]:
    if column not in adata.obs:
        return {}
    counts = adata.obs[column].astype(str).value_counts(dropna=False)
    return {str(k): int(v) for k, v in counts.items()}


def _row(spec, adata: ad.AnnData) -> dict[str, Any]:
    batch_keys = [col for col in ("sample", "donor", "patient", "condition") if col in adata.obs]
    return {
        "dataset": spec.key,
        "path": str(spec.path),
        "exists": True,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "has_counts_layer": "counts" in adata.layers,
        "var_names_unique": bool(adata.var_names.is_unique),
        "candidate_batch_keys": ";".join(batch_keys),
        "preprocess_roles": ";".join(spec.preprocess_roles),
        "figure3_panels": ";".join(spec.figure3_panels),
        "samples": json.dumps(_value_counts(adata, "sample"), ensure_ascii=False),
        "donors": json.dumps(_value_counts(adata, "donor"), ensure_ascii=False),
        "conditions": json.dumps(_value_counts(adata, "condition"), ensure_ascii=False),
        "cell_types_top": json.dumps(
            dict(list(_value_counts(adata, "cell_type").items())[:10]),
            ensure_ascii=False,
        ),
        "benchmark_notes": spec.benchmark_notes,
    }


def build_manifest(output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    readiness_rows: list[dict[str, Any]] = []
    comparison_rows: list[dict[str, Any]] = []
    figure_rows: list[dict[str, Any]] = []
    marker_rows: list[dict[str, Any]] = []
    issue_rows: list[dict[str, Any]] = []

    for spec in DATASETS:
        if not spec.preprocess_roles:
            continue
        if not spec.path.exists():
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
        row = _row(spec, adata)
        readiness_rows.append(row)
        if not row["has_counts_layer"]:
            issue_rows.append(
                {
                    "dataset": spec.key,
                    "severity": "blocker",
                    "issue": "missing_counts_layer",
                    "detail": "layers['counts'] is required for layer-contract validation.",
                }
            )

        for comparison in PREPROCESS_COMPARISONS:
            run_required = comparison in {
                "sclucid_standard_no_integration",
                "sclucid_standard_with_diagnostic_only",
            }
            integration_candidate = any(
                role in spec.preprocess_roles
                for role in (
                    "batch_diagnostic",
                    "patient_integration_diagnostic",
                    "stimulation_batch_diagnostic",
                )
            )
            comparison_rows.append(
                {
                    "dataset": spec.key,
                    "comparison": comparison,
                    "required_for_phase3": bool(run_required or integration_candidate),
                    "primary_batch_key": (
                        "donor"
                        if "donor" in adata.obs
                        else "patient"
                        if "patient" in adata.obs
                        else "sample"
                        if "sample" in adata.obs
                        else ""
                    ),
                    "biological_key_to_protect": (
                        "condition"
                        if "condition" in adata.obs
                        else "cell_type"
                        if "cell_type" in adata.obs
                        else ""
                    ),
                    "notes": spec.benchmark_notes,
                }
            )

        if spec.annotation_obs:
            marker_rows.append(
                {
                    "dataset": spec.key,
                    "marker_panel_type": "author_cell_type_markers",
                    "source": "obs cell_type/cell_subtype labels",
                    "planned_metrics": "marker_inclusion_rate;cluster_marker_recovery;cell_type_retention",
                    "notes": "Use broad lineage markers first; subtype markers are secondary.",
                }
            )
        if "tumor_marker_preservation" in spec.preprocess_roles or "marker_preservation" in spec.preprocess_roles:
            marker_rows.append(
                {
                    "dataset": spec.key,
                    "marker_panel_type": "tumor_programs",
                    "source": "scLucid marker/gene-set resources",
                    "planned_metrics": "tumor_program_gene_retention;marker_fidelity_after_hvg",
                    "notes": "Do not treat HVG exclusion as failure if retained-gene policy protects the program.",
                }
            )

        for panel in spec.figure3_panels:
            figure_rows.append(
                {
                    "figure_panel": panel,
                    "dataset": spec.key,
                    "evidence_role": ";".join(spec.preprocess_roles),
                    "planned_metric_table": {
                        "3A": "layer_contract_report.tsv",
                        "3B": "hvg_marker_preservation.tsv",
                        "3C": "batch_clustering_stability.tsv",
                        "3D": "inference_semantics_guardrails.tsv",
                    }.get(panel, "preprocess_summary.tsv"),
                }
            )
        adata.file.close()

    paths = {
        "readiness": output_dir / "preprocess_dataset_readiness.tsv",
        "comparisons": output_dir / "preprocess_comparison_matrix.tsv",
        "layer_schema": output_dir / "layer_transition_schema.tsv",
        "marker_plan": output_dir / "marker_preservation_plan.tsv",
        "figure3_plan": output_dir / "figure3_panel_plan.tsv",
        "issues": output_dir / "preprocess_metadata_contract_issues.tsv",
    }
    pd.DataFrame(readiness_rows).to_csv(paths["readiness"], sep="\t", index=False)
    pd.DataFrame(comparison_rows).to_csv(paths["comparisons"], sep="\t", index=False)
    pd.DataFrame({"column": LAYER_TRANSITION_COLUMNS}).to_csv(
        paths["layer_schema"], sep="\t", index=False
    )
    pd.DataFrame(marker_rows).to_csv(paths["marker_plan"], sep="\t", index=False)
    pd.DataFrame(figure_rows).to_csv(paths["figure3_plan"], sep="\t", index=False)
    pd.DataFrame(issue_rows).to_csv(paths["issues"], sep="\t", index=False)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_outputs/work/preprocess_manifest"),
    )
    args = parser.parse_args()
    paths = build_manifest(args.output_dir)
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
