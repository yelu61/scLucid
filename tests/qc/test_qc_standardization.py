"""Tests for benchmark-grade QC review schema and trace outputs."""

import json

from scLucid import qc
from scLucid.qc.config import (
    DoubletConfig,
    MarkingConfig,
    MetricsReportingConfig,
    QCThresholds,
    QCWorkflowConfig,
)
from scLucid.qc.trace import (
    QC_REQUIRED_REVIEW_SECTIONS,
    QC_TRACE_SCHEMA_VERSION,
    validate_qc_review_summary,
)
from tests.fixtures.synthetic_data import generate_minimal_adata


def _make_adata():
    return generate_minimal_adata(n_cells=240, n_genes=700)


def _make_config(save_dir=None) -> QCWorkflowConfig:
    return QCWorkflowConfig(
        sample_key="sampleID",
        species="human",
        save_dir=str(save_dir) if save_dir is not None else None,
        use_parallel=False,
        use_recommendations=False,
        threshold_mode="pooled",
        metrics_reporting_config=MetricsReportingConfig(
            show_plots=False,
            plot_top_genes=False,
            plot_violin=False,
            plot_scatter=False,
            export_stats=False,
            print_stats=False,
        ),
        marking_config=MarkingConfig(
            show_plots=False,
            plot_outliers=False,
            thresholds=QCThresholds(min_genes=0, min_counts=0, pc_mt=100.0),
        ),
        doublet_config=DoubletConfig(
            run_algorithm=False,
            use_heuristics=False,
            show_plots=False,
            plot_summary=False,
            plot_bar=False,
            plot_scatter=False,
            plot_upset=False,
            export_stats=False,
        ),
        filter_config={
            "criteria_to_filter": ["outlier_min_genes", "outlier_mt"],
            "combination_logic": "any",
        },
    )


def _qc_review(adata):
    return adata.uns["sclucid"]["qc"]["review_summary"]["data"]


def _as_rows(value):
    return list(value.values()) if isinstance(value, dict) else value


def test_qc_review_summary_has_benchmark_schema():
    adata = qc.run_standard_qc(_make_adata(), config=_make_config(), show_progress=False)
    review = _qc_review(adata)

    assert validate_qc_review_summary(review) == []
    assert QC_REQUIRED_REVIEW_SECTIONS.issubset(review.keys())
    assert review["qc_schema_version"] == QC_TRACE_SCHEMA_VERSION
    assert review["execution_trace"]["qc_schema_version"] == QC_TRACE_SCHEMA_VERSION
    assert review["execution_trace"]["steps_executed"] == ["qc_metrics", "filtering"]
    assert [stage["stage"] for stage in _as_rows(review["policy_flow"])] == [
        "profile_dataset",
        "propose_candidate_thresholds",
        "score_biological_risk",
        "choose_recommend_policy",
        "emit_reviewer_table",
        "optionally_apply",
    ]
    assert "recommended_threshold_summary" in review
    assert "downstream_preprocess_recommendations" in review
    assert "qc_readiness" in review
    assert "review_action_items" in review
    assert "reproducibility_manifest" in review


def test_qc_decision_table_is_machine_readable():
    adata = qc.run_standard_qc(_make_adata(), config=_make_config(), show_progress=False)
    review = _qc_review(adata)
    rows = {row["parameter"]: row for row in review["decision_table"]}

    assert {"min_genes", "max_mt_percent", "doublet_threshold"}.issubset(rows)
    assert rows["min_genes"]["applied"] == 0
    assert rows["min_genes"]["source"] == "default_or_config"
    assert rows["max_mt_percent"]["applied"] == 100.0
    assert rows["max_mt_percent"]["is_filtering_enabled"] is True
    reviewer_columns = {
        "metric",
        "recommended",
        "applied",
        "source",
        "confidence",
        "evidence",
        "review_required",
        "affected_cells",
        "biological_guardrail",
        "risk_note",
    }
    for row in rows.values():
        assert reviewer_columns.issubset(row)
        assert isinstance(row["affected_cells"], int)

    reviewer_table = review["qc_reviewer_table"]
    reviewer_rows = {row["item"]: row for row in reviewer_table}
    assert {"min_genes", "ambient_risk", "stress_high", "predicted_doublet"}.issubset(
        reviewer_rows
    )
    unified_columns = {
        "recommended_value",
        "applied_value",
        "source",
        "confidence",
        "affected_cells",
        "biological_risk_note",
        "review_required",
    }
    for row in reviewer_table:
        assert unified_columns.issubset(row)
        assert isinstance(row["affected_cells"], int)


def test_qc_output_health_and_evidence_chain_are_actionable():
    adata = qc.run_standard_qc(_make_adata(), config=_make_config(), show_progress=False)
    review = _qc_review(adata)

    assert review["output_health"]["status"] == "ok"
    assert review["output_health"]["n_cells"] == adata.n_obs
    assert review["output_health"]["missing_required_obs_metrics"] == []
    stages = [item["stage"] for item in review["evidence_chain"]]
    assert stages == [
        "recommendation",
        "threshold_application",
        "sample_thresholds",
        "filtering",
        "filtering_policy",
        "retention_audit",
        "ambient_evidence",
        "post_annotation_qc_review",
        "qc_benchmark_scorecard",
        "output_health",
    ]


def test_qc_review_summary_records_downstream_preprocess_recommendations():
    adata = qc.run_standard_qc(_make_adata(), config=_make_config(), show_progress=False)
    review = _qc_review(adata)
    downstream = review["downstream_preprocess_recommendations"]

    assert downstream["ready_for_preprocess"] is True
    assert downstream["status"] == "ready"
    targets = {item["target"] for item in downstream["recommendations"]}
    assert {"counts_layer", "normalization"}.issubset(targets)
    assert downstream["input_assumptions"]["sample_key"] == "sampleID"


def test_qc_review_summary_records_doublet_evidence():
    adata = _make_adata()
    adata.obs["predicted_doublet"] = False
    adata.obs.iloc[:12, adata.obs.columns.get_loc("predicted_doublet")] = True
    adata.obs["combined_doublet_score"] = 0.05
    adata.obs.iloc[:12, adata.obs.columns.get_loc("combined_doublet_score")] = 0.9

    adata = qc.run_standard_qc(adata, config=_make_config(), show_progress=False)
    review = _qc_review(adata)
    doublet_summary = review["doublet_evidence_summary"]

    assert doublet_summary["status"] == "available"
    assert doublet_summary["predictions"]["predicted_doublet"]["count"] == 12
    assert doublet_summary["scores"]["combined_doublet_score"]["max"] == 0.9

    compact = qc.summarize_qc_review_summary(review)
    assert compact["doublet_status"] == "available"
    assert compact["predicted_doublets"] == 12
    assert compact["benchmark_status"] in {"pass", "review_required", "fail"}
    assert compact["benchmark_next_step"]
    assert compact["top_review_action"]


def test_qc_review_summary_records_ambient_post_annotation_and_scorecard():
    adata = _make_adata()
    adata.obs["cell_type"] = ["T"] * 80 + ["B"] * 80 + ["Myeloid"] * 80
    adata.obs["ambient_fraction"] = 0.01
    adata.obs.iloc[:20, adata.obs.columns.get_loc("ambient_fraction")] = 0.8
    adata.obs["cell_probability"] = 0.98

    adata = qc.run_standard_qc(adata, config=_make_config(), show_progress=False)
    review = _qc_review(adata)

    assert review["ambient_evidence_summary"]["schema_version"] == "ambient_evidence_summary_v1"
    assert review["ambient_evidence_summary"]["cell_probability"]["available"] is True
    assert review["post_annotation_qc_review"]["schema_version"] == "post_annotation_qc_review_v1"
    assert review["post_annotation_qc_review"]["cell_type_key"] == "cell_type"
    assert review["qc_benchmark_scorecard"]["schema_version"] == "qc_benchmark_scorecard_v1"
    assert review["qc_benchmark_scorecard"]["rows"]

    compact = qc.summarize_qc_review_summary(review)
    assert compact["ambient_status"] in {"available", "not_run"}
    assert compact["cell_probability_available"] is True
    assert compact["post_annotation_qc_available"] is True
    assert compact["qc_benchmark_scorecard_status"] in {"pass", "review_required", "partial"}


def test_qc_review_summary_includes_attached_doublet_benchmark_evidence():
    adata = _make_adata()
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {})[
        "doublet_benchmark_evidence"
    ] = {
        "schema_version": "doublet_benchmark_evidence_v1",
        "best_method": "scdblfinder_python_pyscdblfinder",
        "best_method_f1": 0.61,
        "best_method_auc": 0.88,
        "algorithm_weight_recommendations": [
            {
                "base_method": "scdblfinder_python_pyscdblfinder",
                "recommended_default_mode": "algorithm_only_with_heuristic_review_evidence",
                "recommended_algorithm_weight": 0.7,
                "recommended_method": "scdblfinder_python_pyscdblfinder_plus_heuristic_w0.70",
                "f1_delta_vs_algorithm_only": 0.009,
                "precision_delta_vs_algorithm_only": 0.079,
                "recall_delta_vs_algorithm_only": -0.09,
                "review_required": True,
                "risk_note": "Heuristic fusion does not materially improve F1.",
            }
        ],
        "threshold_calibration_review": [
            {"method": "scrublet", "recall_gain_vs_default": 0.5}
        ],
    }

    adata = qc.run_standard_qc(adata, config=_make_config(), show_progress=False)
    review = _qc_review(adata)
    benchmark = review["doublet_evidence_summary"]["benchmark_evidence"]

    assert benchmark["schema_version"] == "doublet_benchmark_evidence_v1"
    assert benchmark["best_method"] == "scdblfinder_python_pyscdblfinder"
    decision = review["doublet_evidence_summary"]["benchmark_decision"]
    assert decision["recommended_default_mode"] == (
        "algorithm_only_with_heuristic_review_evidence"
    )
    assert decision["recommended_primary_method"] == "scdblfinder_python_pyscdblfinder"
    assert decision["recommended_algorithm_weight"] == 0.7
    assert decision["review_required"] is True

    compact = qc.summarize_qc_review_summary(review)
    assert compact["doublet_recommended_default_mode"] == (
        "algorithm_only_with_heuristic_review_evidence"
    )
    assert compact["doublet_recommended_algorithm_weight"] == 0.7


def test_recommend_qc_policy_is_diagnostic_only_and_apply_consumes_policy():
    adata = _make_adata()
    original_obs_cols = set(adata.obs.columns)

    policy = qc.recommend_qc_policy(
        adata,
        config=_make_config(),
        show_progress=False,
    )

    assert set(adata.obs.columns) == original_obs_cols
    assert policy["schema_version"] == "qc_policy_bundle_v1"
    assert policy["mode"] == "recommend_only"
    assert policy["decision_table"]
    assert policy["policy_flow"][-1]["stage"] == "optionally_apply"
    assert policy["filtering_policy_summary"]["final_filter_basis"] == (
        "legacy_threshold_filtering"
    )
    assert policy["recommended_execution"]["entrypoint"] == "scLucid.qc.run_iterative_qc"
    assert policy["recommended_execution"]["final_filter_policy"] == "decision_remove"

    applied = qc.apply_qc_policy(
        adata,
        policy=policy,
        show_progress=False,
    )
    review = _qc_review(applied)
    assert review["decision_table"]
    assert applied.n_obs <= adata.n_obs


def test_qc_filtering_policy_summary_flags_legacy_filtering_for_review():
    adata = qc.run_standard_qc(_make_adata(), config=_make_config(), show_progress=False)
    review = _qc_review(adata)
    policy = review["qc_filtering_policy_summary"]

    assert policy["qc_decision_filter_mode"] == "off"
    assert policy["final_filter_basis"] == "legacy_threshold_filtering"
    assert policy["review_required"] is True
    assert "reviewer-first" in policy["risk_note"]
    assert "qc_filtering_policy_summary" in {
        item["evidence_key"] for item in review["review_action_items"]
    }

    compact = qc.summarize_qc_review_summary(review)
    assert compact["final_filter_basis"] == "legacy_threshold_filtering"
    assert compact["filtering_policy_review_required"] is True


def test_qc_retention_audit_summary_records_stratified_review():
    adata = _make_adata()
    adata.obs["condition"] = ["control"] * 120 + ["treated"] * 120

    adata = qc.run_standard_qc(adata, config=_make_config(), show_progress=False)
    review = _qc_review(adata)
    retention = review["qc_retention_audit_summary"]

    assert retention["schema_version"] == "qc_retention_audit_summary_v1"
    assert retention["available"] is True
    assert retention["retention_review_required"] is True
    assert retention["group_keys_reviewed"]["sample"] == "sampleID"
    assert retention["group_keys_reviewed"]["condition"] == "condition"
    assert "annotation" in retention["group_keys_reviewed"]
    sample_rows = _as_rows(retention["tables"]["sample"])
    condition_rows = _as_rows(retention["tables"]["condition"])
    assert sample_rows
    assert {row["group"] for row in condition_rows} == {"control", "treated"}
    assert all("retention_rate" in row for row in sample_rows + condition_rows)
    assert "qc_retention_audit_summary" in {
        item["evidence_key"] for item in review["review_action_items"]
    }

    compact = qc.summarize_qc_review_summary(review)
    assert compact["retention_audit_available"] is True
    assert compact["retention_audit_review_required"] is True
    assert compact["retention_audit_group_keys"]["condition"] == "condition"


def test_qc_readiness_and_reproducibility_manifest_are_reviewable():
    adata = qc.run_standard_qc(_make_adata(), config=_make_config(), show_progress=False)
    review = _qc_review(adata)

    assert review["qc_readiness"]["status"] in {"ready", "review_required"}
    assert review["qc_readiness"]["score"] > 0
    assert review["review_action_items"]
    priorities = {item["priority"] for item in review["review_action_items"]}
    assert "required" in priorities

    manifest = review["reproducibility_manifest"]
    assert manifest["workflow"] == "run_standard_qc"
    assert manifest["steps_executed"] == ["qc_metrics", "filtering"]
    assert manifest["data_shape"]["n_obs"] == adata.n_obs
    assert manifest["context"]["sample_key"] == "sampleID"
    assert manifest["applied_thresholds"]["min_genes"] == 0


def test_qc_review_summary_contains_shared_evidence_bundle():
    adata = qc.run_standard_qc(_make_adata(), config=_make_config(), show_progress=False)
    review = _qc_review(adata)
    bundle = review["evidence_bundle"]

    assert bundle["module"] == "qc"
    assert bundle["stage"] == "run_standard_qc"
    assert bundle["status"] == review["qc_readiness"]["status"]
    assert bundle["reproducibility"]["workflow"] == "run_standard_qc"
    assert len(bundle["decisions"]) == len(review["decision_table"])
    assert any(decision["parameter"] == "min_genes" for decision in bundle["decisions"])
    assert bundle["related_review_keys"] == [
        "decision_table",
        "qc_reviewer_table",
        "evidence_chain",
        "qc_handoff_readiness",
        "qc_readiness",
        "review_action_items",
        "reproducibility_manifest",
        "qc_filtering_policy_summary",
        "qc_retention_audit_summary",
        "ambient_evidence_summary",
        "post_annotation_qc_review",
        "qc_benchmark_scorecard",
        "benchmark_summary",
    ]
    assert any(item["name"] == "qc_benchmark_assessment" for item in bundle["evidence_chain"])


def test_qc_review_export_includes_benchmark_schema(tmp_path):
    output_dir = tmp_path / "qc"
    _ = qc.run_standard_qc(
        _make_adata(),
        config=_make_config(save_dir=output_dir),
        show_progress=False,
    )

    payload = json.loads((output_dir / "qc_review_summary.json").read_text())
    markdown = (output_dir / "qc_review_summary.md").read_text()
    benchmark_markdown = (output_dir / "qc_benchmark.md").read_text()
    assert payload["qc_schema_version"] == QC_TRACE_SCHEMA_VERSION
    assert isinstance(payload["decision_table"], list)
    assert isinstance(payload["evidence_chain"], list)
    assert "Executive Summary" in markdown
    assert "Reviewer next step" in markdown
    assert "How To Read This Benchmark" in benchmark_markdown


def test_validate_qc_review_summary_reports_missing_sections():
    errors = validate_qc_review_summary({"decision_table": [], "evidence_chain": []})
    assert errors
    assert "missing required sections" in errors[0]
