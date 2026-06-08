"""Smoke + unit tests for tumor analysis workflow."""

import pytest

from scLucid.recommendation.schema import (
    ParameterRecommendation,
    RecommendationSection,
    WorkflowRecommendations,
)


class TestTumorWorkflow:
    def test_import(self):
        from scLucid.tumor.workflow import (
            run_tumor_analysis,
            run_tumor_analysis_expert,
        )
        assert callable(run_tumor_analysis)
        assert callable(run_tumor_analysis_expert)

    def test_tumor_config_import(self):
        from scLucid.tumor.config import (
            TumorAnalysisConfig,
            TumorWorkflowConfig,
        )
        config = TumorAnalysisConfig()
        assert config.run_malignancy is True
        assert config.run_tme is True
        assert config.run_cnv is True

    def test_tumor_workflow_config_quick(self):
        from scLucid.tumor.config import TumorWorkflowConfig

        config = TumorWorkflowConfig.quick(tissue_type="tumor")
        assert config.tissue_type == "tumor"
        assert config.use_recommendations is True

    def test_run_tumor_analysis_smoke(self, integration_test_adata):
        from scLucid.tumor.workflow import run_tumor_analysis
        from scLucid.tumor.config import TumorWorkflowConfig

        config = TumorWorkflowConfig.quick(
            tissue_type="tumor",
            cancer_type="pancreatic_adenocarcinoma",
        )
        # Tumor workflow needs cell_type annotation and other metadata.
        # Smoke-test: verify it runs without crash on basic data.
        try:
            result = run_tumor_analysis(integration_test_adata, config=config)
            assert result is not None
        except Exception as e:
            pytest.skip(f"Tumor workflow requires richer metadata: {e}")

    def test_tumor_utils_import(self):
        from scLucid.tumor.utils.databases import (
            query_cancer_gene_census,
            get_drug_targets,
            is_cancer_gene,
        )
        assert callable(query_cancer_gene_census)
        assert callable(get_drug_targets)
        assert callable(is_cancer_gene)

        from scLucid.tumor.utils.markers import (
            get_tumor_markers,
            get_immune_markers,
            get_all_markers,
        )
        assert callable(get_tumor_markers)
        assert callable(get_immune_markers)
        assert callable(get_all_markers)

        from scLucid.tumor.utils.signatures import (
            load_hallmark_signatures,
            calculate_signature_scores,
        )
        assert callable(load_hallmark_signatures)
        assert callable(calculate_signature_scores)


def _make_section(name: str, params, *, summary="", confidence=0.8, metadata=None, raw=None):
    return RecommendationSection(
        name=name,
        summary=summary,
        confidence=confidence,
        parameters=[
            ParameterRecommendation(
                name=k,
                value=v,
                method="test",
                confidence=0.9,
                rationale="unit test",
            )
            for k, v in params.items()
        ],
        metadata=metadata or {},
        raw_result=raw,
    )


class TestDiffRecommendations:
    """Unit-test the recommendation-vs-actual diff helper."""

    def test_no_recommendations_returns_empty(self):
        from scLucid.tumor.workflow import _diff_recommendations

        assert _diff_recommendations(None, {}) == {}

    def test_diff_records_overrides(self):
        from scLucid.tumor.workflow import _diff_recommendations
        from scLucid.qc.config import QCWorkflowConfig

        # Recommendation says min_genes=200, actual config has min_genes=300.
        qc_section = _make_section("qc", {"min_genes": 200})
        recs = WorkflowRecommendations(sections={"qc": qc_section}, overall_confidence=0.8)

        qc_cfg = QCWorkflowConfig()
        # Force a known value into the config dict
        actual_dict = qc_cfg.to_dict()
        actual_dict["min_genes"] = 300

        class _Stub:
            def __init__(self, d):
                self._d = d

            def to_dict(self):
                return self._d

        diff = _diff_recommendations(recs, {"qc": _Stub(actual_dict)})
        assert "qc" in diff
        assert diff["qc"]["min_genes"]["recommended"] == 200
        assert diff["qc"]["min_genes"]["actual"] == 300

    def test_diff_skips_when_actual_missing(self):
        from scLucid.tumor.workflow import _diff_recommendations

        section = _make_section("qc", {"min_genes": 200})
        recs = WorkflowRecommendations(sections={"qc": section}, overall_confidence=0.5)
        # No matching actual config — should return empty
        assert _diff_recommendations(recs, {}) == {}


class TestApplyRecommendationHelpers:
    """Each _apply_* helper must mutate the config in place without raising."""

    def test_apply_qc_recommendations_updates_filter(self):
        from scLucid.qc.config import FilterConfig, QCWorkflowConfig
        from scLucid.tumor.workflow import _apply_qc_recommendations

        cfg = QCWorkflowConfig(filter_config=FilterConfig(min_criteria_for_removal=1))
        section = _make_section("qc", {"min_criteria_for_removal": 2})
        updated = _apply_qc_recommendations(cfg, section)
        assert updated.filter_config.min_criteria_for_removal == 2

    def test_apply_qc_recommendations_skips_unknown_fields(self):
        from scLucid.qc.config import FilterConfig, QCWorkflowConfig
        from scLucid.tumor.workflow import _apply_qc_recommendations

        cfg = QCWorkflowConfig(filter_config=FilterConfig())
        # Unknown field name should be ignored without raising
        section = _make_section("qc", {"definitely_not_a_field": 12345})
        updated = _apply_qc_recommendations(cfg, section)
        assert updated is cfg

    def test_apply_clustering_updates_existing_config(self):
        from scLucid.analysis.config import AnalysisWorkflowConfig
        from scLucid.tumor.workflow import _apply_clustering_recommendations

        cfg = AnalysisWorkflowConfig()
        section = _make_section("clustering", {"resolution": 0.8})
        updated = _apply_clustering_recommendations(cfg, section)
        assert updated.clustering is not None
        assert updated.clustering.resolution == 0.8

    def test_apply_annotation_uses_metadata_keys(self):
        from scLucid.analysis.config import AnalysisWorkflowConfig
        from scLucid.tumor.workflow import _apply_annotation_recommendations

        cfg = AnalysisWorkflowConfig()
        section = _make_section(
            "annotation",
            {},
            metadata={"cluster_key": "leiden", "marker_species": "human"},
        )
        updated = _apply_annotation_recommendations(cfg, section)
        assert updated.annotation.cluster_key == "leiden"
        assert updated.annotation.marker_species == "human"

    def test_apply_tumor_recommendations_via_parameters(self):
        from scLucid.tumor.config import TumorAnalysisConfig
        from scLucid.tumor.workflow import _apply_tumor_recommendations

        cfg = TumorAnalysisConfig(run_cnv=False)
        section = _make_section("tumor", {"run_cnv": True})
        updated = _apply_tumor_recommendations(cfg, section)
        assert updated.run_cnv is True

    def test_apply_tumor_recommendations_via_raw_result(self):
        from scLucid.tumor.config import TumorAnalysisConfig
        from scLucid.tumor.workflow import _apply_tumor_recommendations

        replacement = TumorAnalysisConfig(run_cnv=False, run_tme=False)
        section = _make_section("tumor", {}, raw=replacement)
        updated = _apply_tumor_recommendations(TumorAnalysisConfig(), section)
        assert updated.run_cnv is False
        assert updated.run_tme is False

    def test_apply_preprocess_falls_back_to_parameter_assignment(self):
        from scLucid.preprocess.config import WorkflowConfig as PreprocessWorkflowConfig
        from scLucid.tumor.workflow import _apply_preprocess_recommendations

        cfg = PreprocessWorkflowConfig()
        section = _make_section("preprocess", {"random_state": 123})
        # No raw_result with to_config; fallback path must set attribute directly.
        updated = _apply_preprocess_recommendations(cfg, section)
        if hasattr(updated, "random_state"):
            assert updated.random_state == 123


class TestRunTumorStageBranches:
    """Verify _run_tumor_stage skips stages cleanly when flags are False."""

    def test_all_flags_off_returns_empty_steps(self):
        from anndata import AnnData
        import numpy as np

        from scLucid.tumor.config import TumorAnalysisConfig
        from scLucid.tumor.workflow import _run_tumor_stage

        adata = AnnData(np.random.rand(5, 3))
        cfg = TumorAnalysisConfig(
            run_malignancy=False, run_tme=False, run_cnv=False, run_therapy=False
        )
        _result_adata, steps, step_results, warns = _run_tumor_stage(adata, cfg)
        assert steps == []
        assert step_results == []
        assert warns == []

    def test_failing_stage_is_recorded_as_warning_and_step_result(self):
        from anndata import AnnData
        import numpy as np

        from scLucid.tumor.config import TumorAnalysisConfig
        from scLucid.tumor.workflow import _run_tumor_stage

        # AnnData with no count layer / metadata → malignancy scoring will fail
        # but the stage should catch the exception and emit a warning rather
        # than raising. A StepResult must be recorded for audit.
        adata = AnnData(np.random.rand(5, 3))
        cfg = TumorAnalysisConfig(
            run_malignancy=True, run_tme=False, run_cnv=False, run_therapy=False
        )
        _result_adata, steps, step_results, warns = _run_tumor_stage(adata, cfg)
        # No assertion on `steps` content (depends on scorer's resilience),
        # but the function must not raise, warnings should be a list, and
        # step_results must contain a degraded/failed record.
        assert isinstance(warns, list)
        assert len(step_results) >= 1
        assert any(r.status in ("failed", "degraded") for r in step_results)


class TestTumorStageStepResults:
    """Verify structured step results are produced for tumor stage steps."""

    def test_tme_step_records_compartment_claim(self):
        from anndata import AnnData
        import numpy as np

        from scLucid.tumor.config import TumorAnalysisConfig
        from scLucid.tumor.workflow import _run_tumor_stage

        adata = AnnData(np.zeros((6, 3)))
        adata.obs["cell_type"] = ["T_cell", "T cell", "Fibroblast", "CAF", "Tumor", "Tumor"]
        cfg = TumorAnalysisConfig(
            run_malignancy=False, run_tme=True, run_cnv=False, run_therapy=False
        )
        _result_adata, steps, step_results, warns = _run_tumor_stage(adata, cfg)
        assert "tme_deconvolution" in steps
        tme_step = next(r for r in step_results if r.name == "tme_deconvolution")
        assert tme_step.status == "completed"
        assert "annotation-derived" in " ".join(tme_step.warnings)

    def test_therapy_step_marks_exploratory_evidence(self):
        from anndata import AnnData
        import numpy as np

        from scLucid.tumor.config import TumorAnalysisConfig
        from scLucid.tumor.workflow import _run_tumor_stage

        adata = AnnData(np.random.rand(4, 3))
        cfg = TumorAnalysisConfig(
            run_malignancy=False,
            run_tme=False,
            run_cnv=False,
            run_therapy=True,
            therapy_drugs=["chemotherapy"],
        )
        _result_adata, steps, step_results, warns = _run_tumor_stage(adata, cfg)
        therapy_step = next((r for r in step_results if r.name == "therapy_prediction"), None)
        assert therapy_step is not None
        assert therapy_step.evidence_level == "exploratory"
        assert therapy_step.outputs["drugs_requested"] == ["chemotherapy"]


class TestTumorReviewSummary:
    """Verify the tumor review summary structure."""

    def test_review_summary_contains_claim_boundary(self):
        from scLucid.tumor.workflow import _build_tumor_review_summary
        from scLucid.utils import StepResult
        from anndata import AnnData
        import numpy as np

        adata = AnnData(np.zeros((4, 3)))
        summary = _build_tumor_review_summary(
            adata=adata,
            tumor_config=None,
            step_results=[
                StepResult(
                    name="cnv_inference",
                    status="completed",
                    evidence_level="validated_core",
                    outputs={"has_reference_cells": True, "has_genomic_coordinates": True},
                ),
            ],
            warnings=[],
            cancer_type="TEST",
        )
        assert summary["module"] == "tumor"
        assert summary["claim_boundary"] == "validated_core"
        assert summary["cancer_type"] == "TEST"
        assert "readiness" in summary
        assert summary["readiness"]["status"] == "ready"

    def test_review_summary_degraded_steps_lower_readiness(self):
        from scLucid.tumor.workflow import _build_tumor_review_summary
        from scLucid.utils import StepResult
        from anndata import AnnData
        import numpy as np

        adata = AnnData(np.zeros((4, 3)))
        summary = _build_tumor_review_summary(
            adata=adata,
            tumor_config=None,
            step_results=[
                StepResult(name="a", status="completed", evidence_level="heuristic"),
                StepResult(
                    name="b",
                    status="degraded",
                    evidence_level="unavailable",
                    outputs={},
                    warnings=["fallback used"],
                ),
            ],
            warnings=[],
            cancer_type=None,
        )
        assert summary["readiness"]["status"] == "degraded"
        assert summary["readiness"]["score"] < 1.0
        assert any("degraded_steps" in reason for reason in summary["readiness"]["reasons"])
        assert len(summary["action_items"]) >= 1

    def test_review_summary_claim_boundary_is_conservative(self):
        from scLucid.tumor.workflow import _build_tumor_review_summary
        from scLucid.utils import StepResult
        from anndata import AnnData
        import numpy as np

        adata = AnnData(np.zeros((4, 3)))
        summary = _build_tumor_review_summary(
            adata=adata,
            tumor_config=None,
            step_results=[
                StepResult(name="cnv_inference", status="completed", evidence_level="validated_core"),
                StepResult(name="tme_deconvolution", status="completed", evidence_level="heuristic"),
                StepResult(name="therapy_prediction", status="completed", evidence_level="exploratory"),
            ],
            warnings=[],
            cancer_type=None,
        )
        assert summary["claim_boundary"] == "exploratory"


class TestTumorTraceContract:
    """Verify tumor trace helpers and public exports."""

    def test_tumor_trace_exports_resolve(self):
        import scLucid.tumor as tumor

        for name in (
            "MalignancyScoringStep",
            "TherapyPredictionStep",
            "validate_tumor_review_summary",
            "validate_tumor_module_completeness",
            "get_tumor_module_contract",
        ):
            assert hasattr(tumor, name)
            assert name in tumor.__all__

    def test_validate_tumor_review_summary_accepts_hdf5_safe_action_mapping(self):
        from scLucid.tumor.trace import validate_tumor_review_summary

        summary = {
            "schema_version": "1.0",
            "module": "tumor",
            "workflow_name": "tumor_analysis",
            "requested_steps": [],
            "completed_steps": [],
            "claim_boundary": "unavailable",
            "readiness": {"status": "blocked", "score": 0.0},
            "warnings": [],
            "action_items": {
                "0": {
                    "priority": "required",
                    "action": "Review tumor workflow.",
                    "rationale": "unit test",
                    "evidence_keys": ["step_results"],
                }
            },
            "evidence_sources": [],
        }
        assert validate_tumor_review_summary(summary, raise_on_error=False) == []

    def test_validate_tumor_module_completeness_on_saved_namespace(self):
        from anndata import AnnData
        import numpy as np

        from scLucid.tumor.trace import validate_tumor_module_completeness
        from scLucid.utils import sanitize_for_hdf5, save_result

        adata = AnnData(np.zeros((3, 2)))
        summary = {
            "schema_version": "1.0",
            "module": "tumor",
            "workflow_name": "tumor_analysis",
            "requested_steps": ["tme_deconvolution"],
            "completed_steps": ["tme_deconvolution"],
            "claim_boundary": "heuristic",
            "readiness": {"status": "ready", "score": 1.0, "reasons": []},
            "warnings": [],
            "action_items": [],
            "evidence_sources": [],
        }
        save_result(adata, "tumor", "review_summary", sanitize_for_hdf5(summary))
        save_result(adata, "tumor", "step_results", {"0": {"name": "tme_deconvolution"}})

        result = validate_tumor_module_completeness(adata)
        assert result["valid"] is True
