"""Tests for benchmark-grade preprocessing review summaries."""

from scLucid.preprocess.config import (
    GraphConfig,
    HVGConfig,
    IntegrationConfig,
    ScalingConfig,
    WorkflowConfig,
)
from scLucid.preprocess.trace import (
    PREPROCESS_REQUIRED_REVIEW_SECTIONS,
    validate_preprocessing_review_summary,
)
from scLucid.preprocess.workflow import run_preprocessing
from tests.fixtures.synthetic_data import generate_minimal_adata


def _make_config(*, run_integration: bool = False) -> WorkflowConfig:
    return WorkflowConfig(
        hvg=HVGConfig(n_top_genes=100, flavor="seurat"),
        scaling=ScalingConfig(vars_to_regress=[]),
        graph=GraphConfig(n_pcs=10, n_neighbors=5),
        integration=IntegrationConfig(
            method="harmony" if run_integration else None,
            batch_key="sampleID",
            harmony_params={"max_iter_harmony": 2, "theta": 2.0},
        ),
        run_regression=False,
        run_integration=run_integration,
        run_neighbors=False,
    )


def _review(adata):
    return adata.uns["sclucid"]["preprocess"]["review_summary"]


def test_run_preprocessing_review_summary_has_benchmark_sections():
    adata = generate_minimal_adata(n_cells=160, n_genes=320)
    out = run_preprocessing(
        adata,
        config=_make_config(),
        steps=["normalization", "hvg_selection", "subset_hvg", "scaling", "pca"],
        show_progress=False,
    )

    review = _review(out)
    assert validate_preprocessing_review_summary(review) == []
    assert PREPROCESS_REQUIRED_REVIEW_SECTIONS.issubset(review)
    assert review["module_maturity"]["module"] == "preprocess"
    assert review["module_maturity"]["status"] in {
        "complete",
        "review_required",
        "incomplete",
    }
    assert review["qc_input_context"]["available"] is False
    assert review["applied_parameter_summary"]["normalization"]["output_layer"] == "normalized"
    assert review["applied_parameter_summary"]["hvg_selection"]["requested_n_top_genes"] == 100
    assert review["layer_transition_summary"]["raw_present"] is False
    step_evidence = review["step_evidence_summary"]
    assert step_evidence["status_counts"]["complete"] >= 5
    assert {item["step"] for item in step_evidence["steps"]} >= {
        "normalization",
        "hvg_selection",
        "scaling",
        "pca",
    }
    hvg_step = next(item for item in step_evidence["steps"] if item["step"] == "hvg_selection")
    assert hvg_step["output"]["n_hvg_selected"] > 0
    assert "hvg_selection_evidence_summary" in hvg_step["audit_fields"]
    assert review["hvg_selection_evidence_summary"]["status"] == "ok"
    assert review["hvg_selection_evidence_summary"]["n_hvg_selected"] > 0
    assert review["downstream_analysis_recommendations"]["ready_for_analysis"] is True
    assert review["preprocess_readiness"]["status"] in {"ready", "review_required"}


def test_run_preprocessing_review_summary_contains_evidence_bundle():
    adata = generate_minimal_adata(n_cells=160, n_genes=320)
    out = run_preprocessing(
        adata,
        config=_make_config(),
        steps=["normalization", "hvg_selection", "subset_hvg", "scaling", "pca"],
        show_progress=False,
    )
    bundle = _review(out)["evidence_bundle"]

    assert bundle["module"] == "preprocess"
    assert bundle["stage"] == "run_preprocessing"
    assert bundle["status"] == _review(out)["preprocess_readiness"]["status"]
    assert any(item["name"] == "hvg_selection_evidence_summary" for item in bundle["evidence_chain"])
    assert any(item["name"] == "step_evidence_summary" for item in bundle["evidence_chain"])
    assert "applied_parameter_summary" in bundle["related_review_keys"]
    assert "step_evidence_summary" in bundle["related_review_keys"]


def test_tumor_preprocessing_records_batch_correction_warning():
    adata = generate_minimal_adata(n_cells=180, n_genes=360)
    out = run_preprocessing(
        adata,
        config=_make_config(run_integration=True),
        steps=["normalization", "hvg_selection", "subset_hvg", "scaling", "pca", "batch_correction"],
        tissue_type="lung_tumor",
        show_progress=False,
    )

    review = _review(out)
    warnings = review["tumor_aware_batch_correction_warnings"]
    assert warnings["enabled"] is True
    assert warnings["batch_correction_applied"] is True
    assert warnings["warnings"]
    assert any(
        item["evidence_key"] == "tumor_aware_batch_correction_warnings.warnings"
        for item in review["review_action_items"]
    )


def test_preprocess_module_maturity_and_compact_summary():
    import scLucid as scl

    adata = generate_minimal_adata(n_cells=160, n_genes=320)
    out = run_preprocessing(
        adata,
        config=_make_config(),
        steps=["normalization", "hvg_selection", "subset_hvg", "scaling", "pca"],
        show_progress=False,
    )
    review = _review(out)

    validation = scl.pp.validate_preprocess_module_completeness(out)
    assert validation["valid"] is True
    assert validation["maturity"]["module"] == "preprocess"

    compact = scl.pp.summarize_preprocess_review_summary(review)
    assert compact["module"] == "preprocess"
    assert compact["n_hvg_selected"] == review["hvg_selection_evidence_summary"]["n_hvg_selected"]
    assert compact["actual_n_pcs"] == review["applied_parameter_summary"]["pca"]["actual_n_pcs"]
    assert compact["step_status_counts"]["complete"] >= 5


def test_preprocess_module_completeness_detects_missing_result():
    import scLucid as scl

    adata = generate_minimal_adata(n_cells=80, n_genes=120)
    result = scl.pp.validate_preprocess_module_completeness(adata)

    assert result["valid"] is False
    assert any("review_summary" in issue for issue in result["issues"])


def test_preprocess_module_contract_is_public():
    import scLucid as scl

    contract = scl.pp.get_preprocess_module_contract()

    assert contract["module"] == "preprocess"
    assert "scLucid.preprocess.run_preprocessing" in contract["stable_entrypoints"]
    assert "layer_transition_summary" in contract["required_review_sections"]
    assert "step_evidence_summary" in contract["required_review_sections"]
    assert contract["step_evidence_key"] == "step_evidence_summary"
    assert "adata.layers['normalized']" in contract["expected_outputs"]


def test_preprocess_records_qc_input_context_when_qc_exists():
    import scLucid as scl

    adata = generate_minimal_adata(n_cells=160, n_genes=320)
    qc_config = scl.qc.QCWorkflowConfig(
        save_dir=None,
        use_recommendations=False,
        use_parallel=False,
        metrics_reporting_config=scl.qc.MetricsReportingConfig(show_plots=False),
        marking_config=scl.qc.MarkingConfig(show_plots=False, plot_outliers=False),
        doublet_config=scl.qc.DoubletConfig(
            run_algorithm=False,
            use_heuristics=False,
            show_plots=False,
        ),
        filter_config={"criteria_to_filter": ["predicted_doublet"]},
    )
    adata = scl.qc.run_standard_qc(adata, config=qc_config, show_progress=False)
    out = run_preprocessing(
        adata,
        config=_make_config(),
        steps=["normalization", "hvg_selection", "subset_hvg", "scaling", "pca"],
        show_progress=False,
    )

    qc_context = _review(out)["qc_input_context"]
    assert qc_context["available"] is True
    assert qc_context["qc_readiness_status"] in {"ready", "review_required", "blocked"}
    assert qc_context["counts_layer_present"] is True
