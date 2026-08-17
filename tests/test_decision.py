"""Tests for project planning and unified run review."""

import numpy as np
from anndata import AnnData


def _adata(n_obs=24, n_vars=12):
    counts = np.random.default_rng(12).poisson(2, size=(n_obs, n_vars)).astype(int)
    adata = AnnData(counts.copy())
    adata.layers["counts"] = counts.copy()
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.var_names = [f"gene_{i}" for i in range(n_vars)]
    return adata


def test_plan_analysis_exposes_tumor_design_and_paired_structure():
    import scLucid as scl

    adata = _adata()
    adata.obs["sample"] = np.repeat(["s1", "s2", "s3", "s4"], 6)
    adata.obs["patient"] = np.repeat(["p1", "p1", "p2", "p2"], 6)
    adata.obs["condition"] = np.repeat(["pre", "post", "pre", "post"], 6)
    adata.obs["batch"] = np.repeat(["b1", "b2", "b1", "b2"], 6)

    plan = scl.plan_analysis(
        adata,
        context=scl.ProjectContext(
            dataset_type="tumor_tissue",
            sample_key="sample",
            condition_key="condition",
            batch_key="batch",
            experimental_unit_key="patient",
            study_objective="paired treatment response",
        ),
    )

    assert plan.profile == "multi_sample_tumor"
    assert plan.context.paired_key == "patient"
    assert plan.status == "REVIEW"
    decisions = {item.decision: item for item in plan.decisions}
    assert decisions["raw_counts_source"].status == "READY"
    assert decisions["sample_level_inference"].applied["replicates_per_condition"] == {
        "post": 2,
        "pre": 2,
    }
    assert decisions["malignancy_boundary"].status == "REVIEW"


def test_plan_analysis_blocks_comparative_project_without_sample_metadata():
    import scLucid as scl

    adata = _adata()
    plan = scl.plan_analysis(
        adata,
        context={
            "dataset_type": "tumor_tissue",
            "study_objective": "treatment response comparison",
        },
    )

    assert plan.status == "BLOCKED"
    assert "sample_key" in plan.required_metadata
    assert plan.ready_to_run is False


def test_sample_identity_is_not_silently_reused_as_technical_batch():
    import scLucid as scl

    adata = _adata()
    adata.obs["sample"] = np.repeat(["s1", "s2", "s3", "s4"], 6)
    plan = scl.plan_analysis(
        adata,
        context=scl.ProjectContext(dataset_type="pbmc_or_blood", sample_key="sample"),
    )

    assert plan.context.batch_key is None
    integration = next(item for item in plan.decisions if item.decision == "integration_policy")
    assert integration.status == "READY"
    assert integration.recommended == "no integration by default"


def test_review_run_normalizes_stage_status_actions_and_decisions():
    import scLucid as scl

    adata = _adata()
    adata.uns["sclucid"] = {
        "analysis_context": {
            "dataset_type": "tumor_tissue",
            "sample_key": "sample",
            "experimental_unit_key": "sample",
        },
        "qc": {
            "review_summary": {
                "qc_readiness": {"status": "ready", "score": 90},
                "evidence_bundle": {
                    "decisions": [
                        {
                            "parameter": "max_mt_percent",
                            "recommended": 20,
                            "applied": 20,
                            "source": "recommendation",
                            "confidence": 0.8,
                        }
                    ],
                    "action_items": [],
                },
            }
        },
        "preprocess": {
            "review_summary": {
                "preprocess_readiness": {
                    "status": "review_required",
                    "score": 75,
                    "review_reasons": ["integration_needs_review"],
                },
                "review_action_items": [
                    {
                        "priority": "review",
                        "action": "Compare integrated and unintegrated embeddings.",
                        "rationale": "Protect condition biology.",
                        "evidence_key": "integration_diagnostics",
                    }
                ],
            }
        },
    }

    review = scl.review_run(adata)

    assert review.overall_status == "REVIEW"
    stage_status = {stage.stage: stage.status for stage in review.stages}
    assert stage_status == {
        "qc": "READY",
        "preprocess": "REVIEW",
        "analysis": "NOT_RUN",
        "tumor": "NOT_RUN",
    }
    assert review.show_next_actions() == ["Compare integrated and unintegrated embeddings."]
    table = review.to_frame()
    assert {"stage", "status", "decision", "next_action", "rerun_scope"}.issubset(table.columns)
    assert adata.uns["sclucid"]["run_review"]["overall_status"] == "REVIEW"


def test_run_pipeline_accepts_and_records_analysis_plan(monkeypatch):
    import scLucid as scl

    adata = _adata()
    adata.obs["sample"] = "s1"
    plan = scl.plan_analysis(
        adata,
        stages=["qc"],
        context=scl.ProjectContext(dataset_type="pbmc_or_blood", sample_key="sample"),
    )

    def fake_qc(input_adata, **kwargs):
        input_adata.uns.setdefault("sclucid", {}).setdefault("qc", {}).update(
            {
                "workflow_config": {"species": "human"},
                "review_summary": {"qc_readiness": {"status": "ready", "score": 100}},
            }
        )
        return input_adata

    monkeypatch.setattr(scl, "run_standard_qc", fake_qc)
    result = scl.run_pipeline(adata, plan=plan, show_progress=False)

    assert result.uns["sclucid"]["analysis_plan"]["profile"] == "baseline"
    assert result.uns["sclucid"]["analysis_context"]["dataset_type"] == "pbmc_or_blood"
    assert result.uns["sclucid"]["run_review"]["overall_status"] == "READY"


def test_run_pipeline_rejects_blocked_plan_by_default():
    import pytest

    import scLucid as scl

    adata = _adata()
    plan = scl.plan_analysis(
        adata,
        context=scl.ProjectContext(
            dataset_type="tumor_tissue",
            study_objective="treatment response comparison",
        ),
    )

    with pytest.raises(ValueError, match="Analysis plan is BLOCKED"):
        scl.run_pipeline(adata, plan=plan, show_progress=False)


def test_decision_api_is_exported_at_package_root():
    import scLucid as scl

    for name in (
        "ProjectContext",
        "AnalysisPlan",
        "RunReview",
        "plan_analysis",
        "review_run",
        "run_qc",
    ):
        assert hasattr(scl, name)


def test_plan_and_run_review_render_in_audit_report(tmp_path):
    import scLucid as scl
    from scLucid.utils.sanitize import sanitize_for_hdf5

    adata = _adata()
    adata.obs["sample"] = "s1"
    plan = scl.plan_analysis(
        adata,
        context=scl.ProjectContext(dataset_type="pbmc_or_blood", sample_key="sample"),
    )
    adata.uns["sclucid"] = {"analysis_plan": sanitize_for_hdf5(plan.to_dict())}
    scl.review_run(adata)

    output = scl.export_audit_report(adata, tmp_path / "decision_audit.html")
    html = output.read_text(encoding="utf-8")

    assert "Analysis Plan" in html
    assert "Run Decision Review" in html
    assert "Unified decision status" in html


def test_run_review_survives_h5ad_roundtrip(tmp_path):
    import scanpy as sc

    import scLucid as scl

    adata = _adata()
    adata.uns["sclucid"] = {
        "qc": {"review_summary": {"qc_readiness": {"status": "ready", "score": 90}}}
    }
    scl.review_run(adata)
    path = tmp_path / "review_roundtrip.h5ad"
    adata.write_h5ad(path)

    restored = sc.read_h5ad(path)
    review = scl.review_run(restored, store=False)

    assert review.overall_status == "READY"
    assert {stage.stage: stage.status for stage in review.stages}["qc"] == "READY"
