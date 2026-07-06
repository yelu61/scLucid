"""Review-summary and sidecar helpers for QC workflows."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
from anndata import AnnData

from ...utils import (
    UnsKeys,
    normalize_review_summary,
    record_artifact,
    save_result,
    validate_review_summary_schema,
)
from ...utils.context import is_tumor_context, resolve_cell_type_key
from ..ambient import (
    AMBIENT_CORRECTED_COUNTS_LAYER,
    correct_ambient_rna_linear,
    diagnose_ambient_rna,
    diagnose_empty_droplets,
    infer_ambient_input_context,
    record_ambient_layer_contract,
)
from ..ambient_backends import correct_ambient_rna as correct_ambient_rna_unified
from ..artifacts import record_benchmark_review
from ..config import QCWorkflowConfig
from ..policy.benchmark import evaluate_qc_benchmark, export_qc_benchmark_report
from ..trace import enrich_qc_review_summary, validate_qc_review_summary

log = logging.getLogger(__name__)


def _is_tumor_aware(tissue_type: Optional[str]) -> bool:
    return is_tumor_context(tissue_type)


def _diff_qc_recommendations(
    recommendation: Any,
    original_config: QCWorkflowConfig,
) -> Dict[str, Any]:
    """Compare recommended values against the original user config.

    This captures genuine user-vs-recommendation divergence.
    """
    if recommendation is None:
        return {}
    diffs: Dict[str, Any] = {}
    rec_dict = recommendation.to_dict() if hasattr(recommendation, "to_dict") else {}
    cfg_dict = original_config.to_dict()

    mapping = {
        "min_genes": ("marking_config", "thresholds", "min_genes"),
        "max_mt_percent": ("marking_config", "thresholds", "pc_mt"),
        "n_counts": ("marking_config", "thresholds", "min_counts"),
        "doublet_threshold": ("doublet_config", "score_threshold"),
    }
    explicit_field_checks = {
        "min_genes": [
            ("marking_config",),
            ("thresholds",),
            ("min_genes",),
        ],
        "max_mt_percent": [
            ("marking_config",),
            ("thresholds",),
            ("pc_mt",),
        ],
        "n_counts": [
            ("marking_config",),
            ("thresholds",),
            ("min_counts",),
        ],
        "doublet_threshold": [
            ("doublet_config",),
            ("score_threshold",),
        ],
    }

    def _is_explicit_user_path(config_obj: Any, fields: list[tuple[str, ...]]) -> bool:
        current = config_obj
        for field_path in fields:
            field_name = field_path[0]
            if current is None or field_name not in getattr(current, "model_fields_set", set()):
                return False
            current = getattr(current, field_name, None)
        return True

    for param_name, path in mapping.items():
        rec_val = None
        param_rec = rec_dict.get(param_name)
        if isinstance(param_rec, dict):
            rec_val = param_rec.get("threshold")

        if not _is_explicit_user_path(original_config, explicit_field_checks[param_name]):
            continue

        actual_val = cfg_dict
        for key in path:
            if isinstance(actual_val, dict):
                actual_val = actual_val.get(key)
            else:
                actual_val = None
                break

        if rec_val is not None and actual_val is not None and rec_val != actual_val:
            diffs[param_name] = {"recommended": rec_val, "actual": actual_val}

    return diffs


def _build_qc_review_summary(
    config: QCWorkflowConfig,
    original_config: QCWorkflowConfig,
    recommendation: Any,
    sample_thresholds: Dict[str, Any],
    filtering_summary: Dict[str, Any],
    warnings: List[str],
    ambient_summary: Optional[Dict[str, Any]] = None,
    empty_droplet_summary: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Build a human-reviewable summary of the QC run.

    This distills the full trace into the artifacts a reviewer needs:
    what was recommended, what was actually applied, what the user
    overrode, per-sample thresholds, and any tumor-aware cautions.
    """
    summary: Dict[str, Any] = {}

    # --- Recommendation summary ---
    rec_summary: Dict[str, Any] = {"available": recommendation is not None}
    if recommendation is not None:
        rec_dict = recommendation.to_dict() if hasattr(recommendation, "to_dict") else {}
        rec_summary["overall_strategy"] = rec_dict.get("overall_strategy", "unknown")
        rec_summary["overall_confidence"] = rec_dict.get("overall_confidence")
        rec_summary["data_quality_score"] = rec_dict.get("data_quality_score")
        rec_summary["concerns"] = rec_dict.get("concerns", [])
        key_thresholds: Dict[str, Any] = {}
        for param, rec_key, path in [
            ("min_genes", "min_genes", ("marking_config", "thresholds", "min_genes")),
            ("max_mt_percent", "max_mt_percent", ("marking_config", "thresholds", "pc_mt")),
            ("n_counts", "n_counts", ("marking_config", "thresholds", "min_counts")),
            ("doublet_threshold", "doublet_threshold", ("doublet_config", "score_threshold")),
        ]:
            rec_val = (
                rec_dict.get(rec_key, {}).get("threshold")
                if isinstance(rec_dict.get(rec_key), dict)
                else None
            )
            cfg_val = original_config.to_dict()
            for key in path:
                cfg_val = cfg_val.get(key) if isinstance(cfg_val, dict) else None
            key_thresholds[param] = {
                "recommended": rec_val,
                "user_provided": cfg_val,
            }
        rec_summary["key_thresholds"] = key_thresholds
    summary["recommendation_summary"] = rec_summary

    # --- Applied threshold summary ---
    th = config.marking_config.thresholds
    summary["applied_threshold_summary"] = {
        "min_genes": th.min_genes,
        "max_genes": th.max_genes,
        "min_counts": th.min_counts,
        "max_counts": th.max_counts,
        "pc_mt": th.pc_mt,
        "pc_hb": th.pc_hb,
        "nmads": th.nmads,
    }

    # --- User override summary ---
    overrides = _diff_qc_recommendations(recommendation, original_config)
    summary["user_override_summary"] = {
        "overrides_detected": bool(overrides),
        "details": overrides,
        "note": (
            "User-specified thresholds take precedence over recommendations. "
            "Empty details means the user accepted all recommendations or no recommendation was generated."
        ),
    }

    # --- Sample-level threshold summary ---
    n_samples = len(sample_thresholds)
    summary["sample_threshold_summary"] = {
        "mode": config.threshold_mode,
        "n_samples_with_thresholds": n_samples,
        "per_sample": (
            {
                sample: {
                    metric: {
                        "lower": (
                            round(vals["lower"], 2)
                            if isinstance(vals.get("lower"), (int, float))
                            else vals.get("lower")
                        ),
                        "upper": (
                            round(vals["upper"], 2)
                            if isinstance(vals.get("upper"), (int, float))
                            else vals.get("upper")
                        ),
                    }
                    for metric, vals in thresholds.items()
                }
                for sample, thresholds in sample_thresholds.items()
            }
            if sample_thresholds
            else {}
        ),
        "note": (
            "Per-sample thresholds are only computed in hierarchical/independent mode with >1 sample. "
            "Pooled mode uses a single global threshold."
        ),
    }

    # --- Tumor-aware summary ---
    is_tumor = _is_tumor_aware(config.tissue_type)
    tumor_notes: List[str] = []
    if is_tumor:
        tumor_notes.append(
            "Tumor-aware QC is active: elevated mitochondrial content is flagged rather than hard-filtered."
        )
        if "outlier_mt" not in config.filter_config.criteria_to_filter:
            tumor_notes.append(
                "Mitochondrial outlier filtering was disabled for this tumor dataset."
            )
        if config.marking_config.thresholds.pc_mt is not None:
            tumor_notes.append(
                "The mitochondrial threshold is retained as a warning signal for review and reporting."
            )
    tumor_warnings = [note for note in tumor_notes if "disabled" in note or "warning" in note]
    summary["tumor_aware_summary"] = {
        "enabled": is_tumor,
        "tissue_type": config.tissue_type,
        "notes": tumor_notes,
        "warnings": tumor_warnings,
        "filtering_criteria": list(config.filter_config.criteria_to_filter),
        "mitochondrial_filtering_enabled": "outlier_mt" in config.filter_config.criteria_to_filter,
    }

    # --- Filtering summary ---
    fs = filtering_summary if isinstance(filtering_summary, dict) else {}
    summary["filtering_summary"] = {
        "initial_cells": fs.get("initial_cells"),
        "final_cells": fs.get("final_cells"),
        "removed_cells": fs.get("removed_cells"),
        "removed_fraction": fs.get("removed_fraction"),
        "criteria_used": fs.get("criteria_used", config.filter_config.criteria_to_filter),
        "criteria_counts": fs.get("criteria_counts", {}),
        "review_criteria_counts": fs.get("review_criteria_counts", {}),
    }

    # --- Warnings ---
    summary["warnings"] = warnings
    if ambient_summary is not None:
        summary["ambient_rna_summary"] = ambient_summary
    if empty_droplet_summary is not None:
        summary["empty_droplet_summary"] = empty_droplet_summary

    return summary


def _sync_ambient_corrected_layer_to_output(
    *,
    source: AnnData,
    target: AnnData,
    output_layer: str,
) -> bool:
    """Copy an ambient-corrected layer from the review source to the final output."""
    if output_layer not in source.layers:
        return False
    if not target.var_names.equals(source.var_names):
        return False
    try:
        corrected = source.layers[output_layer]
        target.layers[output_layer] = corrected[target.obs_names, :].copy()
        return True
    except Exception:
        try:
            row_index = source.obs_names.get_indexer(target.obs_names)
            if np.any(row_index < 0):
                return False
            target.layers[output_layer] = source.layers[output_layer][row_index, :].copy()
            return True
        except Exception:
            return False


def _store_qc_trace(
    adata: AnnData,
    config: QCWorkflowConfig,
    original_config: QCWorkflowConfig,
    recommendation: Any,
    sample_thresholds: Dict[str, Any],
    filtering_summary: Dict[str, Any],
    warnings: List[str],
    steps_executed: Optional[List[str]] = None,
    adata_before_filtering: Optional[AnnData] = None,
) -> None:
    """Store unified QC trace under adata.uns['sclucid']['qc']."""
    filtering_summary = dict(filtering_summary or {})
    review_input = adata_before_filtering if adata_before_filtering is not None else adata
    if (
        _is_tumor_aware(config.tissue_type)
        and config.marking_config.thresholds.pc_mt is not None
        and "pct_counts_mt" in review_input.obs
    ):
        mt_values = np.asarray(review_input.obs["pct_counts_mt"], dtype=float)
        mt_threshold = float(config.marking_config.thresholds.pc_mt)
        review_counts = dict(filtering_summary.get("review_criteria_counts", {}) or {})
        review_counts["outlier_mt"] = int(np.sum(mt_values >= mt_threshold))
        filtering_summary["review_criteria_counts"] = review_counts

    n_samples = int(adata.obs[config.sample_key].nunique()) if config.sample_key in adata.obs else 1
    context = {
        "sample_key": config.sample_key,
        "threshold_mode": config.threshold_mode,
        "n_samples": n_samples,
        "tissue_type": config.tissue_type,
        "use_recommendations": config.use_recommendations,
    }
    save_result(
        adata,
        "qc",
        "context",
        context,
    )
    if recommendation is not None:
        save_result(adata, "qc", "recommendation", recommendation.to_dict())
    save_result(adata, "qc", "original_config", original_config.to_dict())
    save_result(adata, "qc", "applied_config", config.to_dict())
    save_result(
        adata, "qc", "user_overrides", _diff_qc_recommendations(recommendation, original_config)
    )
    save_result(adata, "qc", "sample_thresholds", sample_thresholds)
    save_result(adata, "qc", "filtering_summary", filtering_summary)
    save_result(adata, "qc", "warnings", warnings)
    benchmark_summary = None
    if adata_before_filtering is not None:
        benchmark_summary = evaluate_qc_benchmark(
            adata_before_filtering,
            adata,
            tissue_type=config.tissue_type,
            tissue=config.tissue,
            sample_key=config.sample_key,
            cell_type_key=_detect_cell_type_key(adata_before_filtering),
        )
        save_result(adata, "qc", "benchmark_summary", benchmark_summary)
        record_benchmark_review(adata, benchmark_summary=benchmark_summary)

    ambient_input = adata_before_filtering if adata_before_filtering is not None else adata
    try:
        ambient_input_context = infer_ambient_input_context(ambient_input)
        ambient_summary = diagnose_ambient_rna(ambient_input)
        ambient_summary["input_context"] = ambient_input_context
        save_result(adata, "qc", "ambient_rna_summary", ambient_summary)
        empty_droplet_summary = diagnose_empty_droplets(ambient_input)
        save_result(adata, "qc", "empty_droplet_summary", empty_droplet_summary)
        if ambient_summary.get("risk_level") in {"moderate", "high"}:
            warnings.append(
                "Ambient RNA risk is "
                f"{ambient_summary.get('risk_level')}; inspect ambient_rna_summary "
                "and consider Python backends such as CellBender or scAR."
            )
            if config.ambient_correction in {"linear", "auto"}:
                try:
                    if config.ambient_correction == "auto":
                        correction_summary = correct_ambient_rna_unified(
                            ambient_input,
                            method="auto",
                            output_layer=AMBIENT_CORRECTED_COUNTS_LAYER,
                            empty_droplet_key="likely_empty_droplet"
                            if "likely_empty_droplet" in ambient_input.obs.columns
                            else None,
                        )
                    else:
                        correction_summary = correct_ambient_rna_linear(
                            ambient_input,
                            output_layer=AMBIENT_CORRECTED_COUNTS_LAYER,
                            empty_droplet_key="likely_empty_droplet"
                            if "likely_empty_droplet" in ambient_input.obs.columns
                            else None,
                        )
                    if correction_summary.get("corrected"):
                        correction_summary["output_layer_synced_to_filtered_adata"] = (
                            _sync_ambient_corrected_layer_to_output(
                                source=ambient_input,
                                target=adata,
                                output_layer=str(
                                    correction_summary.get(
                                        "output_layer",
                                        AMBIENT_CORRECTED_COUNTS_LAYER,
                                    )
                                ),
                            )
                        )
                        backend = correction_summary.get("backend", "linear")
                        warnings.append(
                            f"Applied {backend} ambient RNA correction to layer "
                            f"'{correction_summary.get('output_layer')}'. "
                            f"Removed {correction_summary.get('removed_counts', 0):.0f} "
                            f"counts (mean rho = {correction_summary.get('mean_rho', 0):.3f})."
                        )
                    save_result(adata, "qc", "ambient_correction_summary", correction_summary)
                    record_ambient_layer_contract(
                        adata,
                        input_context=ambient_input_context,
                        correction_summary=correction_summary,
                        output_layer=AMBIENT_CORRECTED_COUNTS_LAYER,
                    )
                except Exception as corr_exc:
                    warnings.append(f"Ambient RNA correction failed: {corr_exc}")
                    record_ambient_layer_contract(
                        adata,
                        input_context=ambient_input_context,
                        correction_summary={
                            "corrected": False,
                            "reason": f"ambient correction failed: {corr_exc}",
                            "review_required": True,
                        },
                        output_layer=AMBIENT_CORRECTED_COUNTS_LAYER,
                    )
            else:
                record_ambient_layer_contract(
                    adata,
                    input_context=ambient_input_context,
                    correction_summary={"corrected": False, "reason": "correction_not_requested"},
                    output_layer=AMBIENT_CORRECTED_COUNTS_LAYER,
                )
        if "ambient_layer_contract" not in adata.uns.get("sclucid", {}).get("qc", {}):
            record_ambient_layer_contract(
                adata,
                input_context=ambient_input_context,
                correction_summary={
                    "corrected": False,
                    "reason": "ambient_risk_below_correction_threshold",
                },
                output_layer=AMBIENT_CORRECTED_COUNTS_LAYER,
            )
        save_result(adata, "qc", "warnings", warnings)
    except Exception as exc:
        warnings.append(f"Ambient RNA diagnostics failed: {exc}")
        ambient_summary = {
            "available": False,
            "risk_level": "unknown",
            "reason": f"ambient diagnostic failed: {exc}",
            "review_required": True,
        }
        empty_droplet_summary = {
            "available": False,
            "risk_level": "unknown",
            "reason": f"empty-droplet diagnostic failed: {exc}",
            "review_required": True,
        }
        save_result(adata, "qc", "ambient_rna_summary", ambient_summary)
        save_result(adata, "qc", "empty_droplet_summary", empty_droplet_summary)
        save_result(adata, "qc", "warnings", warnings)

    # Build and store the review-facing summary
    base_review_summary = _build_qc_review_summary(
        config,
        original_config,
        recommendation,
        sample_thresholds,
        filtering_summary,
        warnings,
        ambient_summary=ambient_summary,
        empty_droplet_summary=empty_droplet_summary,
    )
    if benchmark_summary is not None:
        base_review_summary["benchmark_summary"] = benchmark_summary

    review_summary = normalize_review_summary(
        enrich_qc_review_summary(
            base_review_summary,
            adata=adata,
            config=config,
            original_config=original_config,
            recommendation=recommendation,
            sample_thresholds=sample_thresholds,
            filtering_summary=filtering_summary,
            warnings=warnings,
            context=context,
            steps_executed=steps_executed,
            adata_before_filtering=adata_before_filtering,
        ),
        module="qc",
        workflow_name="standard",
        adata=adata,
        steps_executed=steps_executed or [],
        config=config.to_dict(),
        warnings=warnings,
    )
    validate_review_summary_schema(review_summary, module="qc", raise_on_error=True)
    validate_qc_review_summary(review_summary, raise_on_error=True)
    save_result(adata, "qc", UnsKeys.REVIEW_SUMMARY, review_summary)
    return review_summary


def _detect_cell_type_key(adata: AnnData) -> Optional[str]:
    """Detect a likely cell type annotation column for benchmark stratification."""
    return resolve_cell_type_key(adata)


def _md_cell(value: Any) -> str:
    """Render a value safely inside a Markdown table cell."""
    if value is None:
        return ""
    return str(value).replace("\n", " ").replace("|", "/")


def _export_qc_review_summary(
    review_summary: Dict[str, Any],
    save_dir: Path,
    adata: Optional[AnnData] = None,
) -> None:
    """Export review summary as JSON and Markdown sidecars."""
    save_dir.mkdir(parents=True, exist_ok=True)

    # JSON sidecar
    json_path = save_dir / "qc_review_summary.json"
    json_path.write_text(json.dumps(review_summary, indent=2, default=str), encoding="utf-8")

    # Markdown sidecar
    readiness = review_summary.get("qc_readiness", {})
    fs = review_summary.get("filtering_summary", {})
    benchmark = review_summary.get("benchmark_summary", {})
    benchmark_assessment = (
        benchmark.get("assessment", {}) if isinstance(benchmark, dict) else {}
    )
    benchmark_guide = (
        benchmark_assessment.get("interpretation_guide", {})
        if isinstance(benchmark_assessment, dict)
        else {}
    )
    action_items = review_summary.get("review_action_items", [])
    if isinstance(action_items, dict):
        action_items = list(action_items.values())
    top_action = action_items[0] if action_items and isinstance(action_items[0], dict) else {}
    _removed_frac = fs.get("removed_fraction")
    _removed_frac_str = f"{_removed_frac:.1%}" if isinstance(_removed_frac, (int, float)) else "N/A"

    md_lines = [
        "# QC Review Summary",
        "",
        "## Executive Summary",
        "",
        f"- **Readiness**: {readiness.get('status')} ({readiness.get('verdict')})",
        f"- **Cells retained**: {fs.get('final_cells')} of {fs.get('initial_cells')} ({_removed_frac_str} removed)",
        f"- **Benchmark**: {benchmark.get('status', 'not_available')} / {benchmark_assessment.get('risk_level', 'unknown')} risk",
        f"- **Main benchmark risk**: {benchmark_guide.get('main_risk', 'No benchmark guide was generated.')}",
        f"- **Reviewer next step**: {top_action.get('action', benchmark_guide.get('next_step', 'Archive QC outputs with the analysis record.'))}",
        "",
        "## Recommendation Summary",
        "",
    ]
    rec = review_summary.get("recommendation_summary", {})
    if rec.get("available"):
        md_lines.append(f"- **Strategy**: {rec.get('overall_strategy', 'unknown')}")
        md_lines.append(f"- **Confidence**: {rec.get('overall_confidence')}")
        md_lines.append(f"- **Data Quality Score**: {rec.get('data_quality_score')}")
        if rec.get("concerns"):
            md_lines.append("- **Concerns**:")
            for c in rec["concerns"]:
                md_lines.append(f"  - {c}")
        md_lines.append("")
        md_lines.append("| Parameter | Recommended | User Provided |")
        md_lines.append("|-----------|-------------|---------------|")
        for param, vals in rec.get("key_thresholds", {}).items():
            md_lines.append(
                f"| {param} | {vals.get('recommended')} | {vals.get('user_provided')} |"
            )
    else:
        md_lines.append(
            "- No recommendation was generated (recommendations disabled or engine failed)."
        )
    md_lines.append("")

    md_lines.extend(
        [
            "## QC Readiness",
            "",
            f"- **Status**: {readiness.get('status')}",
            f"- **Score**: {readiness.get('score')}",
            f"- **Verdict**: {readiness.get('verdict')}",
            "",
        ]
    )
    if readiness.get("blockers"):
        md_lines.append("- **Blockers**:")
        for blocker in readiness.get("blockers", []):
            md_lines.append(f"  - {blocker}")
    if readiness.get("review_reasons"):
        md_lines.append("- **Review reasons**:")
        for reason in readiness.get("review_reasons", []):
            md_lines.append(f"  - {reason}")
    md_lines.append("")

    if action_items:
        md_lines.extend(
            [
                "## Review Action Items",
                "",
                "| Priority | Action | Rationale |",
                "|----------|--------|-----------|",
            ]
        )
        for item in action_items:
            md_lines.append(
                "| {priority} | {action} | {rationale} |".format(
                    priority=_md_cell(item.get("priority")),
                    action=_md_cell(item.get("action")),
                    rationale=_md_cell(item.get("rationale")),
                )
            )
        md_lines.append("")

    md_lines.extend(
        [
            "## Threshold Reviewer Table",
            "",
            "| Parameter | Recommended | Applied | Source | Confidence | Affected Cells | Risk Note | Review Required |",
            "|-----------|-------------|---------|--------|------------|----------------|-----------|-----------------|",
        ]
    )
    decision_table = review_summary.get(
        "threshold_reviewer_table",
        review_summary.get("decision_table", []),
    )
    if isinstance(decision_table, dict):
        decision_table = list(decision_table.values())
    for row in decision_table:
        md_lines.append(
            "| {parameter} | {recommended} | {applied} | {source} | {confidence} | {affected} | {risk_note} | {review_required} |".format(
                parameter=_md_cell(row.get("parameter")),
                recommended=_md_cell(row.get("recommended")),
                applied=_md_cell(row.get("applied")),
                source=_md_cell(row.get("source")),
                confidence=_md_cell(row.get("confidence")),
                affected=_md_cell(row.get("affected_cells")),
                risk_note=_md_cell(row.get("risk_note")),
                review_required=_md_cell(row.get("review_required")),
            )
        )
    md_lines.append("")

    qc_decisions = review_summary.get("qc_decision_summary", {})
    if qc_decisions:
        decision_counts = qc_decisions.get("decision_counts", {})
        evidence_summary = qc_decisions.get("evidence_summary", {})
        md_lines.extend(
            [
                "## QC Decision Summary",
                "",
                f"- **Policy**: {qc_decisions.get('policy')}",
                f"- **Cells**: {qc_decisions.get('n_cells')}",
                f"- **Review required cells**: {qc_decisions.get('review_required_cells')}",
                f"- **Risk note**: {qc_decisions.get('risk_note')}",
                "",
                "| Decision | Cells |",
                "|----------|-------|",
            ]
        )
        if isinstance(decision_counts, dict):
            for decision, count in decision_counts.items():
                md_lines.append(f"| {decision} | {count} |")
        md_lines.extend(["", "| Evidence | Cells |", "|----------|-------|"])
        if isinstance(evidence_summary, dict):
            for evidence, count in evidence_summary.items():
                md_lines.append(f"| {evidence} | {count} |")
        md_lines.append("")

    health = review_summary.get("output_health", {})
    md_lines.extend(
        [
            "## Output Health",
            "",
            f"- **Status**: {health.get('status')}",
            f"- **Cells**: {health.get('n_cells')}",
            f"- **Genes**: {health.get('n_genes')}",
        ]
    )
    issues = health.get("issues", [])
    if isinstance(issues, dict):
        issues = list(issues.values())
    if issues:
        md_lines.append("- **Issues**:")
        for issue in issues:
            md_lines.append(f"  - {issue}")
    md_lines.append("")

    if benchmark:
        retention = benchmark.get("retention", {})
        marker = benchmark.get("marker_fidelity", {})
        md_lines.extend(
            [
                "## Benchmark Summary",
                "",
                f"- **Profile**: {benchmark.get('profile_label')} ({benchmark.get('profile')})",
                f"- **Status**: {benchmark.get('status')}",
                f"- **Risk level**: {benchmark_assessment.get('risk_level', 'unknown')}",
                f"- **Retention rate**: {retention.get('retention_rate')}",
                f"- **Marker fidelity**: {marker.get('overall_marker_fidelity')}",
                f"- **Interpretation**: {benchmark_assessment.get('summary', 'No benchmark assessment was generated.')}",
                f"- **Next step**: {benchmark_guide.get('next_step', 'Archive QC benchmark outputs with the analysis record.')}",
                "",
            ]
        )

    md_lines.extend(
        [
            "## Applied Thresholds",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
        ]
    )
    for param, val in review_summary.get("applied_threshold_summary", {}).items():
        md_lines.append(f"| {param} | {val} |")
    md_lines.append("")

    ov = review_summary.get("user_override_summary", {})
    md_lines.extend(
        [
            "## User Overrides",
            "",
            f"- **Overrides detected**: {ov.get('overrides_detected', False)}",
        ]
    )
    if ov.get("details"):
        md_lines.append("- **Details**:")
        for param, vals in ov["details"].items():
            md_lines.append(
                f"  - {param}: recommended={vals.get('recommended')}, user={vals.get('actual')}"
            )
    md_lines.append("")

    st = review_summary.get("sample_threshold_summary", {})
    md_lines.extend(
        [
            "## Sample-Level Thresholds",
            "",
            f"- **Mode**: {st.get('mode')}",
            f"- **Samples with thresholds**: {st.get('n_samples_with_thresholds', 0)}",
            "",
        ]
    )
    if st.get("per_sample"):
        md_lines.append("```json")
        md_lines.append(json.dumps(st["per_sample"], indent=2, default=str))
        md_lines.append("```")
    md_lines.append("")

    ta = review_summary.get("tumor_aware_summary", {})
    md_lines.extend(
        [
            "## Tumor-Aware QC",
            "",
            f"- **Enabled**: {ta.get('enabled', False)}",
        ]
    )
    if ta.get("notes"):
        for note in ta["notes"]:
            md_lines.append(f"- {note}")
    md_lines.append("")

    doublet_summary = review_summary.get("doublet_evidence_summary", {})
    if doublet_summary:
        final_doublets = (
            doublet_summary.get("predictions", {}).get("predicted_doublet", {})
            if isinstance(doublet_summary.get("predictions"), dict)
            else {}
        )
        predicted_fraction = final_doublets.get("fraction")
        predicted_fraction_str = (
            f"{predicted_fraction:.1%}" if isinstance(predicted_fraction, (int, float)) else "N/A"
        )
        md_lines.extend(
            [
                "## Doublet Evidence",
                "",
                f"- **Status**: {doublet_summary.get('status')}",
                f"- **Predicted doublets**: {final_doublets.get('count', 'N/A')} ({predicted_fraction_str})",
                f"- **Review required**: {doublet_summary.get('review_required', False)}",
            ]
        )
        notes = doublet_summary.get("notes", [])
        if isinstance(notes, dict):
            notes = list(notes.values())
        for note in notes:
            md_lines.append(f"- {note}")
        method_keys = doublet_summary.get("method_metadata_keys", [])
        if method_keys:
            md_lines.append(f"- **Method metadata**: {method_keys}")
        benchmark_decision = doublet_summary.get("benchmark_decision", {})
        if isinstance(benchmark_decision, dict) and benchmark_decision:
            md_lines.extend(
                [
                    "- **Benchmark decision**:",
                    f"  - recommended_default_mode: {benchmark_decision.get('recommended_default_mode', 'N/A')}",
                    f"  - recommended_primary_method: {benchmark_decision.get('recommended_primary_method', 'N/A')}",
                    f"  - recommended_algorithm_weight: {benchmark_decision.get('recommended_algorithm_weight', 'N/A')}",
                    f"  - review_required: {benchmark_decision.get('review_required', False)}",
                ]
            )
            if benchmark_decision.get("risk_note"):
                md_lines.append(f"  - risk_note: {benchmark_decision['risk_note']}")
        benchmark_evidence = doublet_summary.get("benchmark_evidence", {})
        if benchmark_evidence:
            md_lines.append("- **Benchmark evidence**:")
            if isinstance(benchmark_evidence, dict):
                for key, value in benchmark_evidence.items():
                    md_lines.append(f"  - {key}: {value}")
            else:
                md_lines.append(f"  - {benchmark_evidence}")
        md_lines.append("")

    downstream = review_summary.get("downstream_preprocess_recommendations", {})
    md_lines.extend(
        [
            "## Downstream Preprocess Recommendations",
            "",
            f"- **Status**: {downstream.get('status')}",
            f"- **Ready for preprocess**: {downstream.get('ready_for_preprocess')}",
            "",
        ]
    )
    _recs = downstream.get("recommendations", [])
    if isinstance(_recs, dict):
        _recs = list(_recs.values())
    for item in _recs:
        md_lines.append(
            "- **{target}** ({priority}): {recommendation}".format(
                target=item.get("target"),
                priority=item.get("priority"),
                recommendation=item.get("recommendation"),
            )
        )
    md_lines.append("")

    md_lines.extend(
        [
            "## Filtering Results",
            "",
            f"- **Initial cells**: {fs.get('initial_cells')}",
            f"- **Final cells**: {fs.get('final_cells')}",
            f"- **Removed**: {fs.get('removed_cells')} ({_removed_frac_str})",
            f"- **Criteria used**: {fs.get('criteria_used', [])}",
            "",
        ]
    )

    warnings = review_summary.get("warnings", [])
    if isinstance(warnings, dict):
        warnings = list(warnings.values())
    if warnings:
        md_lines.extend(
            [
                "## Warnings",
                "",
            ]
        )
        for w in warnings:
            md_lines.append(f"- {w}")
        md_lines.append("")

    md_path = save_dir / "qc_review_summary.md"
    md_path.write_text("\n".join(md_lines), encoding="utf-8")
    benchmark_report_paths: Dict[str, Path] = {}
    benchmark_summary = review_summary.get("benchmark_summary")
    if isinstance(benchmark_summary, dict):
        try:
            benchmark_report_paths = export_qc_benchmark_report(benchmark_summary, save_dir)
        except Exception as e:  # pragma: no cover - sidecar export should not fail QC
            log.warning("Failed to export QC benchmark report: %s", e)

    if adata is not None:
        record_artifact(
            adata,
            "qc",
            "qc_review_summary_json",
            str(json_path),
            kind="json",
            description="QC review summary JSON sidecar",
        )
        record_artifact(
            adata,
            "qc",
            "qc_review_summary_md",
            str(md_path),
            kind="md",
            description="QC review summary Markdown sidecar",
        )
        for key, path in benchmark_report_paths.items():
            artifact_path = Path(path)
            record_artifact(
                adata,
                "qc",
                f"qc_benchmark_{key}",
                str(artifact_path),
                kind=artifact_path.suffix.lstrip(".") or "file",
                description=f"QC benchmark {key} sidecar",
            )
    log.info(f"QC review summary exported to {json_path} and {md_path}")
