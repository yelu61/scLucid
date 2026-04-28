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
    assert review["applied_parameter_summary"]["normalization"]["output_layer"] == "normalized"
    assert review["applied_parameter_summary"]["hvg_selection"]["requested_n_top_genes"] == 100
    assert review["layer_transition_summary"]["raw_present"] is False
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
    assert "applied_parameter_summary" in bundle["related_review_keys"]


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
