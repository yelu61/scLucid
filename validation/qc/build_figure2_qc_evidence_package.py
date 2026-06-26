#!/usr/bin/env python3
"""Build a unified Figure 2 QC evidence package from validation outputs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.dataset_registry import DATASETS


DEFAULT_INPUT_ROOT = Path("validation_outputs")
DEFAULT_OUTPUT_DIR = Path("validation_outputs/qc_figure2_package")


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


def _source_rows_from_figure_table(
    table: pd.DataFrame,
    *,
    evidence_domain: str,
    source_file: str,
    panel_map: dict[str, str] | None = None,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if table.empty:
        return rows
    for _, row in table.iterrows():
        context = _json_context(row.get("context", ""))
        source_panel = str(row.get("figure_panel"))
        figure_panel = panel_map.get(source_panel, source_panel) if panel_map else source_panel
        rows.append(
            {
                "figure_panel": figure_panel,
                "evidence_domain": evidence_domain,
                "dataset": row.get("dataset"),
                "strategy_or_method": row.get("strategy"),
                "metric": row.get("metric"),
                "value": row.get("value"),
                "context": json.dumps(context, sort_keys=True),
                "source_file": source_file,
            }
        )
    return rows


def _doublet_source_rows(root: Path) -> list[dict[str, Any]]:
    output_dir = root / "qc_doublet_evidence"
    evidence = _read_tsv(output_dir / "doublet_evidence.tsv")
    weights = _read_tsv(output_dir / "doublet_algorithm_weight_recommendation.tsv")
    parity = _read_tsv(output_dir / "scdblfinder_python_vs_r_reference.tsv")
    calibration = _read_tsv(output_dir / "doublet_threshold_calibration.tsv")
    ambient = _read_tsv(root / "qc_ambient_evidence" / "ambient_evidence.tsv")

    rows: list[dict[str, Any]] = []
    for _, row in evidence.iterrows() if not evidence.empty else []:
        for metric in ("precision", "recall", "f1", "score_auc", "predicted_rate"):
            rows.append(
                {
                    "figure_panel": "2D",
                    "evidence_domain": "doublet_demuxlet_benchmark",
                    "dataset": row.get("dataset"),
                    "strategy_or_method": row.get("method"),
                    "metric": metric,
                    "value": row.get(metric),
                    "context": json.dumps(
                        {
                            "recommendation_role": row.get("recommendation_role"),
                            "uses_heuristics": row.get("uses_heuristics"),
                            "algorithm_weight": row.get("algorithm_weight"),
                            "review_required": row.get("review_required"),
                        },
                        sort_keys=True,
                    ),
                    "source_file": str(output_dir / "doublet_evidence.tsv"),
                }
            )

    for _, row in weights.iterrows() if not weights.empty else []:
        for metric in (
            "f1_delta_vs_algorithm_only",
            "precision_delta_vs_algorithm_only",
            "recall_delta_vs_algorithm_only",
            "recommended_algorithm_weight",
        ):
            rows.append(
                {
                    "figure_panel": "2D",
                    "evidence_domain": "doublet_algorithm_weight_recommendation",
                    "dataset": row.get("dataset"),
                    "strategy_or_method": row.get("base_method"),
                    "metric": metric,
                    "value": row.get(metric),
                    "context": json.dumps(
                        {
                            "recommended_default_mode": row.get("recommended_default_mode"),
                            "recommended_method": row.get("recommended_method"),
                            "review_required": row.get("review_required"),
                            "risk_note": row.get("risk_note"),
                        },
                        sort_keys=True,
                    ),
                    "source_file": str(output_dir / "doublet_algorithm_weight_recommendation.tsv"),
                }
            )

    for _, row in parity.iterrows() if not parity.empty else []:
        for metric in (
            "score_spearman",
            "prediction_agreement",
            "prediction_jaccard",
            "python_f1",
            "r_f1",
        ):
            rows.append(
                {
                    "figure_panel": "2D",
                    "evidence_domain": "doublet_python_r_parity",
                    "dataset": row.get("dataset"),
                    "strategy_or_method": row.get("python_method"),
                    "metric": metric,
                    "value": row.get(metric),
                    "context": json.dumps(
                        {
                            "reference_status": row.get("reference_status"),
                            "review_required": row.get("review_required"),
                            "risk_note": row.get("risk_note"),
                        },
                        sort_keys=True,
                    ),
                    "source_file": str(output_dir / "scdblfinder_python_vs_r_reference.tsv"),
                }
            )

    if not calibration.empty:
        flagged = calibration[calibration.get("review_required", False).astype(bool)]
        for _, row in flagged.iterrows():
            rows.append(
                {
                    "figure_panel": "2D",
                    "evidence_domain": "doublet_threshold_calibration",
                    "dataset": row.get("dataset"),
                    "strategy_or_method": row.get("method"),
                    "metric": f"recall_gain_at_target_{row.get('target_recall')}",
                    "value": row.get("recall_gain_vs_default"),
                    "context": json.dumps(
                        {
                            "calibrated_threshold": row.get("calibrated_threshold"),
                            "default_recall": row.get("default_recall"),
                            "calibrated_recall": row.get("calibrated_recall"),
                            "risk_note": row.get("risk_note"),
                        },
                        sort_keys=True,
                    ),
                    "source_file": str(output_dir / "doublet_threshold_calibration.tsv"),
                }
            )

    for _, row in ambient.iterrows() if not ambient.empty else []:
        rows.append(
            {
                "figure_panel": "2A",
                "evidence_domain": "ambient_contract",
                "dataset": row.get("dataset"),
                "strategy_or_method": row.get("diagnostic"),
                "metric": "cell_to_empty_median_count_ratio",
                "value": row.get("cell_to_empty_median_count_ratio"),
                "context": json.dumps(
                    {
                        "review_required": row.get("review_required"),
                        "risk_note": row.get("risk_note"),
                    },
                    sort_keys=True,
                ),
                "source_file": str(root / "qc_ambient_evidence" / "ambient_evidence.tsv"),
            }
        )
    return rows


def _dataset_coverage_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for spec in DATASETS:
        rows.append(
            {
                "dataset": spec.key,
                "path": str(spec.path),
                "available": bool(spec.path.exists()),
                "tissue": spec.tissue,
                "modality_role": spec.modality_role,
                "qc_roles": ";".join(spec.qc_roles),
                "figure2_panels": ";".join(spec.figure2_panels),
                "benchmark_notes": spec.benchmark_notes,
            }
        )
    return rows


def _claim_scorecard(
    *,
    threshold_scorecard: pd.DataFrame,
    threshold_decisions: pd.DataFrame,
    tumor_scorecard: pd.DataFrame,
    doublet_evidence: pd.DataFrame,
    weight_recommendations: pd.DataFrame,
    parity: pd.DataFrame,
    ambient: pd.DataFrame,
    dataset_coverage: pd.DataFrame,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    available = dataset_coverage[dataset_coverage["available"]].copy()
    threshold_datasets = (
        int(threshold_scorecard["dataset"].nunique()) if not threshold_scorecard.empty else 0
    )
    threshold_strategies = (
        int(threshold_scorecard["strategy"].nunique()) if not threshold_scorecard.empty else 0
    )
    required_decision_cols = {
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
        "biological_guardrail",
        "risk_note",
        "decision_narrative",
    }
    decision_cols_present = required_decision_cols.issubset(set(threshold_decisions.columns))
    rows.append(
        {
            "claim_id": "qc_decision_auditability",
            "claim": "QC emits reviewer-facing adaptive/sample/tumor-aware threshold decisions.",
            "evidence_status": "supported" if threshold_datasets >= 5 and decision_cols_present else "partial",
            "supporting_datasets": threshold_datasets,
            "supporting_methods": threshold_strategies,
            "key_metric": "decision_table_schema_complete",
            "key_value": bool(decision_cols_present),
            "review_required": not decision_cols_present,
            "limitations": "Decision quality is proxy-based unless dataset-specific QC labels exist.",
            "next_action": "Keep adding dataset-specific biological harm proxies and reviewer narratives.",
        }
    )

    tumor_datasets = int(tumor_scorecard["dataset"].nunique()) if not tumor_scorecard.empty else 0
    tumor_aware = (
        tumor_scorecard[tumor_scorecard["strategy"] == "sclucid_tumor_aware"].copy()
        if not tumor_scorecard.empty
        else pd.DataFrame()
    )
    tumor_aware_harm_rate = (
        float(tumor_aware["biological_harm_risk"].astype(bool).mean())
        if not tumor_aware.empty and "biological_harm_risk" in tumor_aware
        else np.nan
    )
    mean_fidelity = (
        float(pd.to_numeric(tumor_aware["biological_fidelity_score"], errors="coerce").mean())
        if not tumor_aware.empty
        else np.nan
    )
    tumor_status = (
        "supported"
        if tumor_datasets >= 3 and pd.notna(mean_fidelity) and mean_fidelity >= 0.8
        else "partial"
    )
    rows.append(
        {
            "claim_id": "tumor_aware_biological_fidelity",
            "claim": "Tumor-aware QC reduces mechanical removal of high-mt biological signal.",
            "evidence_status": tumor_status,
            "supporting_datasets": tumor_datasets,
            "supporting_methods": int(tumor_scorecard["strategy"].nunique())
            if not tumor_scorecard.empty
            else 0,
            "key_metric": "mean_sclucid_tumor_aware_biological_fidelity_score",
            "key_value": mean_fidelity,
            "review_required": bool(tumor_status != "supported" or tumor_aware_harm_rate > 0),
            "limitations": "Marker/program retention is a biological-fidelity proxy, not a cell-level QC gold standard.",
            "next_action": "Convert PDAC/NSCLC/CRC results into Figure 2C/2E reviewer panels.",
        }
    )

    methods = set(doublet_evidence.get("method", pd.Series(dtype=str)).astype(str))
    has_doublet_core = {
        "scrublet",
        "scdblfinder_python_pyscdblfinder",
    }.issubset(methods)
    has_weight_rec = not weight_recommendations.empty
    has_parity = (
        not parity.empty and str(parity.iloc[0].get("reference_status", "")) == "ok"
    )
    best_f1 = (
        float(pd.to_numeric(doublet_evidence["f1"], errors="coerce").max())
        if not doublet_evidence.empty and "f1" in doublet_evidence
        else np.nan
    )
    rows.append(
        {
            "claim_id": "doublet_evidence_calibration",
            "claim": "Doublet reporting separates primary algorithm calls, heuristic review evidence, calibration, and Python/R parity.",
            "evidence_status": "supported" if has_doublet_core and has_weight_rec and has_parity else "partial",
            "supporting_datasets": int(doublet_evidence["dataset"].nunique())
            if not doublet_evidence.empty
            else 0,
            "supporting_methods": len(methods),
            "key_metric": "best_demuxlet_grounded_f1",
            "key_value": best_f1,
            "review_required": True,
            "limitations": "Kang demuxlet labels mainly validate genotype-detectable heterotypic donor doublets.",
            "next_action": "Add homotypic/synthetic or HTO evidence and tumor/solid-tissue doublet stress tests.",
        }
    )

    ambient_contract_ready = not ambient.empty and bool(
        ambient.get("counts_layer_present", pd.Series([False])).astype(bool).any()
    )
    rows.append(
        {
            "claim_id": "ambient_diagnostic_contract",
            "claim": "Ambient RNA diagnostic plumbing and empty-droplet contract are available.",
            "evidence_status": "contract_only" if ambient_contract_ready else "missing",
            "supporting_datasets": int(ambient["dataset"].nunique()) if not ambient.empty else 0,
            "supporting_methods": int(ambient["diagnostic"].nunique()) if not ambient.empty else 0,
            "key_metric": "cell_to_empty_median_count_ratio",
            "key_value": float(ambient["cell_to_empty_median_count_ratio"].iloc[0])
            if ambient_contract_ready and "cell_to_empty_median_count_ratio" in ambient
            else np.nan,
            "review_required": True,
            "limitations": "CellBender tiny validates contract only; it does not support performance claims.",
            "next_action": "Add a full raw 10x matrix with empty droplets and SoupX/CellBender reference.",
        }
    )

    rows.append(
        {
            "claim_id": "phase2_dataset_coverage",
            "claim": "Phase 2 QC evidence uses at least five real datasets with distinct roles.",
            "evidence_status": "supported" if int(available["dataset"].nunique()) >= 5 else "partial",
            "supporting_datasets": int(available["dataset"].nunique()),
            "supporting_methods": int(available["modality_role"].nunique()),
            "key_metric": "available_validation_datasets",
            "key_value": int(available["dataset"].nunique()),
            "review_required": False,
            "limitations": "Coverage is broad, but not every claim has a true external ground truth.",
            "next_action": "Maintain DATASETS.md and this scorecard together as datasets evolve.",
        }
    )
    return rows


def _write_markdown_report(
    *,
    output_path: Path,
    claim_rows: list[dict[str, Any]],
    dataset_rows: list[dict[str, Any]],
) -> None:
    claims = pd.DataFrame(claim_rows)
    datasets = pd.DataFrame(dataset_rows)
    lines = [
        "# Figure 2 QC Evidence Package",
        "",
        "This report summarizes the current evidence surface for scLucid QC. It is a",
        "review-oriented source-data index, not a claim that every QC decision has a",
        "cell-level gold standard.",
        "",
        "## Claim Scorecard",
        "",
        "| Claim | Status | Key metric | Key value | Review required | Next action |",
        "|---|---|---|---:|---|---|",
    ]
    for _, row in claims.iterrows():
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
            "| Dataset | Available | Role | Figure 2 panels | Notes |",
            "|---|---:|---|---|---|",
        ]
    )
    for _, row in datasets.iterrows():
        lines.append(
            "| {dataset} | {available} | {role} | {panels} | {notes} |".format(
                dataset=row["dataset"],
                available=row["available"],
                role=row["modality_role"],
                panels=row["figure2_panels"],
                notes=str(row["benchmark_notes"]).replace("|", "/"),
            )
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- Supported claims are evidence-backed by the current validation outputs.",
            "- Partial claims have useful evidence but should remain conservative in docs.",
            "- Contract-only claims validate API/data plumbing but not performance.",
            "",
        ]
    )
    output_path.write_text("\n".join(lines), encoding="utf-8")


def build_package(input_root: Path, output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    threshold_dir = input_root / "qc_threshold_benchmark"
    tumor_dir = input_root / "qc_tumor_fidelity"
    doublet_dir = input_root / "qc_doublet_evidence"
    ambient_dir = input_root / "qc_ambient_evidence"

    threshold_scorecard = _read_tsv(threshold_dir / "qc_strategy_scorecard.tsv")
    threshold_decisions = _read_tsv(threshold_dir / "qc_threshold_decision_table.tsv")
    threshold_figure = _read_tsv(threshold_dir / "figure2_threshold_data.tsv")
    tumor_scorecard = _read_tsv(tumor_dir / "tumor_qc_strategy_scorecard.tsv")
    tumor_figure = _read_tsv(tumor_dir / "figure2_tumor_fidelity_data.tsv")
    doublet_evidence = _read_tsv(doublet_dir / "doublet_evidence.tsv")
    weight_recommendations = _read_tsv(
        doublet_dir / "doublet_algorithm_weight_recommendation.tsv"
    )
    parity = _read_tsv(doublet_dir / "scdblfinder_python_vs_r_reference.tsv")
    ambient = _read_tsv(ambient_dir / "ambient_evidence.tsv")

    dataset_rows = _dataset_coverage_rows()
    dataset_coverage = pd.DataFrame(dataset_rows)
    source_rows = []
    source_rows.extend(
        _source_rows_from_figure_table(
            threshold_figure,
            evidence_domain="threshold_decision_quality",
            source_file=str(threshold_dir / "figure2_threshold_data.tsv"),
        )
    )
    source_rows.extend(
        _source_rows_from_figure_table(
            tumor_figure,
            evidence_domain="tumor_biological_fidelity",
            source_file=str(tumor_dir / "figure2_tumor_fidelity_data.tsv"),
            panel_map={"2D": "2E"},
        )
    )
    source_rows.extend(_doublet_source_rows(input_root))

    claim_rows = _claim_scorecard(
        threshold_scorecard=threshold_scorecard,
        threshold_decisions=threshold_decisions,
        tumor_scorecard=tumor_scorecard,
        doublet_evidence=doublet_evidence,
        weight_recommendations=weight_recommendations,
        parity=parity,
        ambient=ambient,
        dataset_coverage=dataset_coverage,
    )

    paths = {
        "figure2_source_data": output_dir / "figure2_qc_source_data.tsv",
        "claim_scorecard": output_dir / "qc_claim_scorecard.tsv",
        "dataset_coverage": output_dir / "qc_dataset_coverage.tsv",
        "evidence_report": output_dir / "qc_evidence_report.md",
    }
    pd.DataFrame(source_rows).to_csv(paths["figure2_source_data"], sep="\t", index=False)
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
