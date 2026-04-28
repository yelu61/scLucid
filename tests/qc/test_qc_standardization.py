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


def test_qc_review_summary_has_benchmark_schema():
    adata = qc.run_standard_qc(_make_adata(), config=_make_config(), show_progress=False)
    review = _qc_review(adata)

    assert validate_qc_review_summary(review) == []
    assert QC_REQUIRED_REVIEW_SECTIONS.issubset(review.keys())
    assert review["qc_schema_version"] == QC_TRACE_SCHEMA_VERSION
    assert review["execution_trace"]["qc_schema_version"] == QC_TRACE_SCHEMA_VERSION
    assert review["execution_trace"]["steps_executed"] == ["qc_metrics", "filtering"]
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
        "evidence_chain",
        "qc_readiness",
        "review_action_items",
        "reproducibility_manifest",
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
    assert payload["qc_schema_version"] == QC_TRACE_SCHEMA_VERSION
    assert isinstance(payload["decision_table"], list)
    assert isinstance(payload["evidence_chain"], list)


def test_validate_qc_review_summary_reports_missing_sections():
    errors = validate_qc_review_summary({"decision_table": [], "evidence_chain": []})
    assert errors
    assert "missing required sections" in errors[0]
