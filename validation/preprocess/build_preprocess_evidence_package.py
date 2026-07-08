#!/usr/bin/env python3
"""Build a unified preprocess evidence package from validation outputs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.dataset_registry import DATASETS


DEFAULT_INPUT_ROOT = Path("validation_outputs")
DEFAULT_OUTPUT_DIR = Path("validation_outputs/preprocess_evidence_package")


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def _json_context(value: Any) -> dict[str, Any]:
    if isinstance(value, dict):
        return value
    if not isinstance(value, str) or not value:
        return {}
    try:
        parsed = json.loads(value)
        return parsed if isinstance(parsed, dict) else {}
    except Exception:
        return {}


def _source_rows_from_module_table(
    table: pd.DataFrame,
    *,
    evidence_domain: str,
    source_file: str,
    default_panel: str,
    metric_columns: tuple[str, ...],
    method_columns: tuple[str, ...] = ("strategy", "method", "batch_key", "stability_type"),
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if table.empty:
        return rows
    for _, row in table.iterrows():
        panel = row.get("module_panel", row.get("figure_panel", default_panel))
        method = next((row.get(col) for col in method_columns if col in row and pd.notna(row.get(col))), "")
        context = _json_context(row.get("context", ""))
        for metric in metric_columns:
            if metric not in row or pd.isna(row.get(metric)):
                continue
            rows.append(
                {
                    "module_panel": panel,
                    "evidence_domain": evidence_domain,
                    "dataset": row.get("dataset"),
                    "strategy_or_method": method,
                    "metric": metric,
                    "value": row.get(metric),
                    "context": json.dumps(context, sort_keys=True),
                    "source_file": source_file,
                }
            )
    return rows


def _dataset_coverage_rows() -> list[dict[str, Any]]:
    return [
        {
            "dataset": spec.key,
            "path": str(spec.path),
            "available": bool(spec.path.exists()),
            "tissue": spec.tissue,
            "modality_role": spec.modality_role,
            "preprocess_roles": ";".join(spec.preprocess_roles),
            "preprocess_panels": ";".join(spec.figure3_panels),
            "benchmark_notes": spec.benchmark_notes,
        }
        for spec in DATASETS
        if spec.preprocess_roles
    ]


def _claim_scorecard(
    *,
    layer_contract: pd.DataFrame,
    hvg_summary: pd.DataFrame,
    batch_summary: pd.DataFrame,
    graph_stability: pd.DataFrame,
    dataset_coverage: pd.DataFrame,
) -> list[dict[str, Any]]:
    available = dataset_coverage[dataset_coverage["available"]].copy()
    layer_ok = not layer_contract.empty
    hvg_datasets = int(hvg_summary["dataset"].nunique()) if "dataset" in hvg_summary else 0
    batch_datasets = int(batch_summary["dataset"].nunique()) if "dataset" in batch_summary else 0
    graph_datasets = int(graph_stability["dataset"].nunique()) if "dataset" in graph_stability else 0

    return [
        {
            "claim_id": "preprocess_layer_contract",
            "claim": "Preprocessing emits auditable counts-to-graph layer contracts.",
            "evidence_status": "supported" if layer_ok else "missing",
            "supporting_datasets": int(layer_contract["dataset"].nunique()) if "dataset" in layer_contract else 0,
            "key_metric": "layer_contract_table_present",
            "key_value": bool(layer_ok),
            "review_required": not layer_ok,
            "limitations": "Layer evidence verifies semantics and handoff, not biological optimality.",
            "next_action": "Keep consuming QC handoff counts-layer recommendations in real-project runs.",
        },
        {
            "claim_id": "hvg_marker_program_preservation",
            "claim": "HVG selection exposes marker/program preservation evidence.",
            "evidence_status": "supported" if hvg_datasets >= 3 else "partial",
            "supporting_datasets": hvg_datasets,
            "key_metric": "hvg_strategy_dataset_count",
            "key_value": hvg_datasets,
            "review_required": hvg_datasets < 3,
            "limitations": "Marker/program panels are proxies and depend on curation coverage.",
            "next_action": "Extend marker/program preservation checks to active project contexts.",
        },
        {
            "claim_id": "batch_correction_diagnostic_guardrail",
            "claim": "Batch correction remains diagnostic/opt-in and reports overcorrection risk.",
            "evidence_status": "supported" if batch_datasets >= 2 else "partial",
            "supporting_datasets": batch_datasets,
            "key_metric": "batch_diagnostic_dataset_count",
            "key_value": batch_datasets,
            "review_required": True,
            "limitations": "Batch metrics require biological-label context before being interpreted as correction quality.",
            "next_action": "Use real projects to decide which warnings are useful versus noisy.",
        },
        {
            "claim_id": "graph_handoff_stability",
            "claim": "PCA/neighbors/clustering handoff has stability evidence.",
            "evidence_status": "supported" if graph_datasets >= 2 else "partial",
            "supporting_datasets": graph_datasets,
            "key_metric": "graph_stability_dataset_count",
            "key_value": graph_datasets,
            "review_required": graph_datasets < 2,
            "limitations": "Stability does not imply annotation correctness.",
            "next_action": "Connect graph stability with downstream annotation acceptance.",
        },
        {
            "claim_id": "preprocess_dataset_coverage",
            "claim": "Preprocess evidence uses multiple real datasets and roles.",
            "evidence_status": "supported" if int(available["dataset"].nunique()) >= 5 else "partial",
            "supporting_datasets": int(available["dataset"].nunique()),
            "key_metric": "available_preprocess_validation_datasets",
            "key_value": int(available["dataset"].nunique()),
            "review_required": False,
            "limitations": "Not every dataset exercises every preprocess risk.",
            "next_action": "Keep dataset roles explicit as validation inventory evolves.",
        },
    ]


def _write_markdown_report(
    *,
    output_path: Path,
    claim_rows: list[dict[str, Any]],
    dataset_rows: list[dict[str, Any]],
) -> None:
    lines = [
        "# Preprocess Evidence Package",
        "",
        "This report summarizes validation evidence for preprocessing layer semantics,",
        "HVG preservation, batch-correction diagnostics, and graph handoff stability.",
        "",
        "## Claim Scorecard",
        "",
        "| Claim | Status | Key metric | Key value | Review required | Next action |",
        "|---|---|---|---:|---|---|",
    ]
    for row in claim_rows:
        lines.append(
            "| {claim_id} | {status} | {metric} | {value} | {review} | {next_action} |".format(
                claim_id=row["claim_id"],
                status=row["evidence_status"],
                metric=row["key_metric"],
                value=row["key_value"],
                review=row["review_required"],
                next_action=str(row["next_action"]).replace("|", "/"),
            )
        )
    lines.extend(
        [
            "",
            "## Dataset Coverage",
            "",
            "| Dataset | Available | Role | Preprocess panels | Notes |",
            "|---|---:|---|---|---|",
        ]
    )
    for row in dataset_rows:
        lines.append(
            "| {dataset} | {available} | {role} | {panels} | {notes} |".format(
                dataset=row["dataset"],
                available=row["available"],
                role=row["modality_role"],
                panels=row["preprocess_panels"],
                notes=str(row["benchmark_notes"]).replace("|", "/"),
            )
        )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_package(input_root: Path, output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    layer_dir = input_root / "preprocess_layer_contract"
    hvg_dir = input_root / "preprocess_hvg_preservation"
    batch_dir = input_root / "preprocess_batch_diagnostic"
    graph_dir = input_root / "preprocess_graph_stability"

    layer_contract = _read_tsv(layer_dir / "layer_contract_report.tsv")
    hvg_summary = _read_tsv(hvg_dir / "hvg_strategy_summary.tsv")
    hvg_source = _read_tsv(hvg_dir / "figure3_hvg_data.tsv")
    batch_summary = _read_tsv(batch_dir / "batch_diagnostic_summary.tsv")
    batch_source = _read_tsv(batch_dir / "figure3_batch_data.tsv")
    graph_stability = _read_tsv(graph_dir / "pca_neighbors_stability.tsv")
    graph_source = _read_tsv(graph_dir / "figure3_graph_data.tsv")

    dataset_rows = _dataset_coverage_rows()
    dataset_coverage = pd.DataFrame(dataset_rows)
    source_rows: list[dict[str, Any]] = []
    source_rows.extend(
        _source_rows_from_module_table(
            hvg_source,
            evidence_domain="hvg_marker_program_preservation",
            source_file=str(hvg_dir / "figure3_hvg_data.tsv"),
            default_panel="preprocess_hvg",
            metric_columns=("inclusion_rate", "genes_present", "genes_retained", "hvg_set_size"),
        )
    )
    source_rows.extend(
        _source_rows_from_module_table(
            batch_source,
            evidence_domain="batch_correction_diagnostic",
            source_file=str(batch_dir / "figure3_batch_data.tsv"),
            default_panel="preprocess_batch",
            metric_columns=("value",),
        )
    )
    source_rows.extend(
        _source_rows_from_module_table(
            graph_source,
            evidence_domain="graph_handoff_stability",
            source_file=str(graph_dir / "figure3_graph_data.tsv"),
            default_panel="preprocess_graph",
            metric_columns=("value",),
        )
    )

    claim_rows = _claim_scorecard(
        layer_contract=layer_contract,
        hvg_summary=hvg_summary,
        batch_summary=batch_summary,
        graph_stability=graph_stability,
        dataset_coverage=dataset_coverage,
    )

    paths = {
        "preprocess_source_data": output_dir / "preprocess_source_data.tsv",
        "claim_scorecard": output_dir / "preprocess_claim_scorecard.tsv",
        "dataset_coverage": output_dir / "preprocess_dataset_coverage.tsv",
        "evidence_report": output_dir / "preprocess_evidence_report.md",
    }
    pd.DataFrame(source_rows).to_csv(paths["preprocess_source_data"], sep="\t", index=False)
    pd.DataFrame(claim_rows).to_csv(paths["claim_scorecard"], sep="\t", index=False)
    dataset_coverage.to_csv(paths["dataset_coverage"], sep="\t", index=False)
    _write_markdown_report(
        output_path=paths["evidence_report"],
        claim_rows=claim_rows,
        dataset_rows=dataset_rows,
    )
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    paths = build_package(args.input_root, args.output_dir)
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
