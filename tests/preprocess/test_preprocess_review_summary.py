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
    norm_policy = review["normalization_decision_policy"]
    assert norm_policy["schema_version"] == "normalization_decision_policy_v1"
    assert norm_policy["applied_method"] == "standard"
    assert norm_policy["claim_level"] == "standard_preprocessing"
    assert norm_policy["recommended_input_layer"] == "counts"
    method_semantics = review["preprocess_method_semantics"]
    assert method_semantics["schema_version"] == "preprocess_method_semantics_v1"
    assert method_semantics["claim_level_counts"]["standard_preprocessing"] >= 1
    semantics_row_values = (
        method_semantics["rows"].values()
        if isinstance(method_semantics["rows"], dict)
        else method_semantics["rows"]
    )
    semantics_rows = {row["source_key"]: row for row in semantics_row_values}
    assert semantics_rows["normalization"]["claim_level"] == "standard_preprocessing"
    assert review["applied_parameter_summary"]["hvg_selection"]["requested_n_top_genes"] == 100
    layer_contract = review["preprocess_layer_contract"]
    assert layer_contract["canonical_flow"] == "counts -> normalized -> raw -> HVG -> scaled -> PCA -> graph"
    assert layer_contract["normalized_layer"] == "normalized"
    assert layer_contract["raw_source_layer"] == "normalized"
    assert layer_contract["recommended_counts_layer"] == "counts"
    assert {stage["stage"] for stage in layer_contract["stage_contracts"]} >= {
        "counts",
        "normalized",
        "raw",
        "HVG",
        "scaled",
        "PCA",
        "graph",
    }
    assert review["layer_transition_summary"]["raw_present"] is False
    layer_table = review["layer_transition_table"]
    assert {row["step"] for row in layer_table} >= {
        "input",
        "normalization",
        "hvg_selection",
        "scaling",
        "pca",
    }
    normalization_row = next(row for row in layer_table if row["step"] == "normalization")
    assert normalization_row["output_slot"] == "layers['normalized']; adata.X if update_X"
    assert "adata_X_semantics_before" in normalization_row
    assert "raw_semantics" in normalization_row
    assert normalization_row["review_required"] is False
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
    decision_summary = review["preprocess_decision_summary"]
    assert decision_summary["schema_version"] == "preprocess_decision_summary_v1"
    assert decision_summary["canonical_flow"] == "counts -> normalized -> raw -> HVG -> scaled -> PCA -> graph"
    decision_rows = {row["step"]: row for row in decision_summary["decisions"]}
    assert {"normalization", "hvg_selection", "pca", "batch_correction", "neighbors_umap"}.issubset(
        decision_rows
    )
    assert decision_rows["normalization"]["decision"] == "use"
    assert decision_summary["primary_downstream_representation"] == "X_pca"
    reviewer_table = review["preprocess_reviewer_table"]
    reviewer_rows = {row["item"]: row for row in reviewer_table}
    assert {"normalization", "hvg_selection", "regression", "batch_correction"}.issubset(
        reviewer_rows
    )
    required_reviewer_columns = {
        "recommended_value",
        "applied_value",
        "source",
        "confidence",
        "affected_representation",
        "preprocess_decision",
        "review_required",
        "biological_risk_note",
    }
    for row in reviewer_table:
        assert required_reviewer_columns.issubset(row)
        assert "claim_levels" in row
        assert "scientific_semantics" in row
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
    assert any(item["name"] == "normalization_decision_policy" for item in bundle["evidence_chain"])
    assert any(item["name"] == "preprocess_layer_contract" for item in bundle["evidence_chain"])
    assert any(item["name"] == "step_evidence_summary" for item in bundle["evidence_chain"])
    assert any(item["name"] == "preprocess_method_semantics" for item in bundle["evidence_chain"])
    assert any(item["name"] == "preprocess_decision_summary" for item in bundle["evidence_chain"])
    assert "applied_parameter_summary" in bundle["related_review_keys"]
    assert "normalization_decision_policy" in bundle["related_review_keys"]
    assert "preprocess_layer_contract" in bundle["related_review_keys"]
    assert "preprocess_decision_summary" in bundle["related_review_keys"]
    assert "preprocess_reviewer_table" in bundle["related_review_keys"]
    assert "preprocess_method_semantics" in bundle["related_review_keys"]
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
    assert warnings["auto_decide"] is False
    assert warnings["evaluate"] is False
    assert warnings["warnings"]
    assert any("auto_decide" in item for item in warnings["warnings"])
    assert any("evaluate" in item for item in warnings["warnings"])
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
    assert compact["canonical_layer_flow"] == "counts -> normalized -> raw -> HVG -> scaled -> PCA -> graph"
    assert compact["recommended_counts_layer"] == "counts"
    assert compact["normalization_applied_method"] == "standard"
    assert compact["method_semantics_status"] == "ok"
    assert compact["method_claim_level_counts"]["standard_preprocessing"] >= 1
    assert compact["step_status_counts"]["complete"] >= 5
    assert compact["preprocess_decision_counts"]["use"] >= 4
    assert compact["primary_downstream_representation"] == "X_pca"


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
    assert "preprocess_layer_contract" in contract["required_review_sections"]
    assert "layer_transition_table" in contract["required_review_sections"]
    assert contract["layer_contract_key"] == "preprocess_layer_contract"
    assert contract["layer_transition_table_key"] == "layer_transition_table"
    assert "step_evidence_summary" in contract["required_review_sections"]
    assert contract["step_evidence_key"] == "step_evidence_summary"
    assert "preprocess_decision_summary" in contract["required_review_sections"]
    assert "preprocess_reviewer_table" in contract["required_review_sections"]
    assert "preprocess_method_semantics" in contract["required_review_sections"]
    assert "normalization_decision_policy" in contract["required_review_sections"]
    assert contract["normalization_policy_key"] == "normalization_decision_policy"
    assert contract["decision_summary_key"] == "preprocess_decision_summary"
    assert contract["reviewer_table_key"] == "preprocess_reviewer_table"
    assert contract["method_semantics_key"] == "preprocess_method_semantics"
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


def test_preprocess_deviance_hvg_uses_counts_layer_in_workflow():
    adata = generate_minimal_adata(n_cells=160, n_genes=320)
    config = _make_config()
    config.hvg = HVGConfig(method="deviance", n_top_genes=100, exclude_gene_types=[])

    out = run_preprocessing(
        adata,
        config=config,
        steps=["normalization", "hvg_selection", "subset_hvg", "scaling", "pca"],
        show_progress=False,
    )
    review = _review(out)

    hvg_summary = review["hvg_selection_evidence_summary"]
    assert hvg_summary["method"] == "deviance"
    assert hvg_summary["input_layer"] == "counts"
    assert hvg_summary["method_report"]["backend"] == "deviance_poisson_approx"
