"""Tests for dataset context detection and propagation."""

import numpy as np
from anndata import AnnData


def _adata(n_obs=20, n_vars=30):
    counts = np.random.default_rng(0).poisson(1, size=(n_obs, n_vars)).astype(int)
    adata = AnnData(X=counts)
    adata.layers["counts"] = counts.copy()
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.var_names = [f"gene_{i}" for i in range(n_vars)]
    return adata


def test_normalize_dataset_type_aliases():
    from scLucid.utils.context import is_tumor_context, normalize_dataset_type

    assert normalize_dataset_type("PBMC") == "pbmc_or_blood"
    assert normalize_dataset_type("lung_tumor") == "tumor_tissue"
    assert normalize_dataset_type("cell line") == "cell_line"
    assert normalize_dataset_type("Visium") == "spatial"
    assert normalize_dataset_type("multi_sample") == "unknown"
    assert is_tumor_context("lung_tumor")
    assert is_tumor_context("pan-cancer")
    assert is_tumor_context("malignant sample")


def test_infer_analysis_context_detects_obs_and_spatial():
    from scLucid.utils.context import infer_analysis_context

    adata = _adata()
    adata.obs["tissue_type"] = "normal"
    adata.obs["sampleID"] = ["s1"] * 10 + ["s2"] * 10
    adata.obsm["spatial"] = np.random.default_rng(1).normal(size=(adata.n_obs, 2))

    context = infer_analysis_context(adata)

    assert context.dataset_type == "normal_tissue"
    assert context.sample_key == "sampleID"
    assert context.is_multi_sample is True
    assert context.is_spatial is True
    assert context.enables_tumor_module is False


def test_resolve_cell_type_key_uses_canonical_downstream_order():
    from scLucid.utils.context import resolve_cell_lineage_key, resolve_cell_type_key

    adata = _adata()
    adata.obs["cell_type_auto"] = "auto"
    assert resolve_cell_type_key(adata) == "cell_type_auto"

    adata.obs["cell_type"] = "manual"
    assert resolve_cell_type_key(adata) == "cell_type"
    assert resolve_cell_type_key(adata, preferred="cell_type_auto") == "cell_type_auto"

    adata.obs["celltype_lineage_auto"] = "immune"
    assert resolve_cell_lineage_key(adata) == "celltype_lineage_auto"


def test_multi_sample_is_context_axis_not_dataset_type():
    from scLucid.utils.context import infer_analysis_context

    adata = _adata()
    context = infer_analysis_context(adata, dataset_type="multi_sample")

    assert context.dataset_type == "unknown"
    assert context.is_multi_sample is True
    assert any("sample structure" in note for note in context.notes)


def test_run_pipeline_stores_context_and_passes_tumor_aware_qc(monkeypatch):
    import scLucid as scl

    adata = _adata()
    seen = {}

    def fake_qc(input_adata, **kwargs):
        seen.update(kwargs)
        input_adata.uns.setdefault("sclucid", {}).setdefault("qc", {})["workflow_config"] = {
            "species": "human"
        }
        return input_adata

    monkeypatch.setattr(scl, "run_standard_qc", fake_qc)
    out = scl.run_pipeline(
        adata,
        stages=["qc"],
        dataset_type="tumor_tissue",
        tissue_type="lung_tumor",
        show_progress=False,
    )

    assert seen["tissue_type"] == "lung_tumor"
    assert out.uns["sclucid"]["analysis_context"]["dataset_type"] == "tumor_tissue"
    assert out.uns["sclucid"]["pipeline_context"]["dataset_type"] == "tumor_tissue"


def test_run_pipeline_passes_plan_design_context_to_analysis(monkeypatch):
    """The executable analysis stage must receive the design recorded by the plan."""
    import scLucid as scl

    adata = _adata(n_obs=24)
    adata.obs["orig.ident"] = np.repeat(["s1", "s2", "s3", "s4"], 6)
    adata.obs["cohort"] = adata.obs["orig.ident"].map(
        {"s1": "pre", "s2": "pre", "s3": "post", "s4": "post"}
    )
    adata.obs["patient"] = adata.obs["orig.ident"].map(
        {"s1": "p1", "s2": "p2", "s3": "p1", "s4": "p2"}
    )
    adata.obsm["X_pca"] = np.zeros((adata.n_obs, 2))
    adata.uns["sclucid"] = {"preprocess": {}}
    plan = scl.plan_analysis(
        adata,
        stages=["analysis"],
        context=scl.ProjectContext(
            dataset_type="pbmc_or_blood",
            sample_key="orig.ident",
            condition_key="cohort",
            experimental_unit_key="patient",
            study_objective="paired comparison",
        ),
    )
    seen = {}

    def fake_analysis(input_adata, **kwargs):
        seen["context"] = kwargs["context"]
        input_adata.uns["sclucid"]["analysis"] = {
            "review_summary": {"analysis_readiness": {"status": "ready", "score": 100}}
        }
        return input_adata

    monkeypatch.setattr(scl, "run_standard_analysis", fake_analysis)
    result = scl.run_pipeline(adata, plan=plan, show_progress=False)

    assert seen["context"].sample_key == "orig.ident"
    assert seen["context"].condition_key == "cohort"
    assert seen["context"].experimental_unit_key == "patient"
    assert seen["context"].paired_key == "patient"
    assert result.uns["sclucid"]["analysis_context"]["sample_key"] == "orig.ident"


def test_recommendation_tumor_section_disabled_for_pbmc():
    from scLucid.recommendation.config import RecommendationConfig
    from scLucid.recommendation.engine import RecommendationEngine

    adata = _adata(n_obs=600)
    adata.obs["cell_type"] = "T cell"

    engine = RecommendationEngine(config=RecommendationConfig(modules=["tumor"]))
    recs = engine.recommend(adata, dataset_type="pbmc_or_blood")
    tumor = recs.get_section("tumor")

    assert tumor is not None
    assert tumor.raw_result.run_malignancy is False
    assert tumor.raw_result.run_tme is False
    assert tumor.raw_result.run_cnv is False
    assert tumor.metadata["dataset_type"] == "pbmc_or_blood"


def test_recommendation_tumor_section_enabled_for_explicit_cancer_context():
    from scLucid.recommendation.config import RecommendationConfig
    from scLucid.recommendation.engine import RecommendationEngine

    adata = _adata(n_obs=600)
    adata.obs["cell_type"] = "epithelial"

    engine = RecommendationEngine(config=RecommendationConfig(modules=["tumor"]))
    recs = engine.recommend(
        adata,
        context={"is_multi_sample": True},
        cancer_type="melanoma",
    )
    tumor = recs.get_section("tumor")

    assert tumor is not None
    assert tumor.metadata["dataset_type"] == "unknown"
    assert recs.context["is_multi_sample"] is True
    assert tumor.get_parameter("run_malignancy").method != "dataset_type_gate"
