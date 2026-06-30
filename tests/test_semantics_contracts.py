"""Unified schema/contract tests for preprocessing and analysis semantics.

These tests verify that the review summaries produced by the core workflows
contain the explicit scientific-claim metadata required for reproducible,
auditable single-cell analysis.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scLucid.analysis import (
    AnalysisWorkflowConfig,
    AnnotationConfig,
    ClusteringConfig,
    DifferentialConfig,
    ProportionConfig,
    PseudobulkDEConfig,
    analyze_celltype_proportion,
    run_pseudobulk_de,
)
from scLucid.preprocess import PreprocessingWorkflowConfig
from scLucid.qc import (
    DoubletConfig,
    FilterConfig,
    MarkingConfig,
    MetricsReportingConfig,
    QCThresholds,
    QCWorkflowConfig,
)
from scLucid.utils.contracts import validate_review_summary_schema


def _synthetic_adata(n_cells: int = 200, n_genes: int = 500) -> AnnData:
    """Return a reproducible synthetic count matrix for fast pipeline tests."""
    rng = np.random.default_rng(42)
    counts = rng.poisson(2, size=(n_cells, n_genes)).astype(float)
    adata = AnnData(X=counts)
    adata.layers["counts"] = counts.copy()
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]
    adata.obs["sampleID"] = "sample_1"
    return adata


def _fast_qc_config() -> QCWorkflowConfig:
    """Build a minimal QC config that avoids slow algorithmic steps."""
    return QCWorkflowConfig(
        save_dir=None,
        species="human",
        tissue_type="normal_tissue",
        use_parallel=False,
        n_jobs=1,
        metrics_reporting_config=MetricsReportingConfig(
            plot_violin=False,
            plot_scatter=False,
            plot_top_genes=False,
            show_plots=False,
            export_stats=False,
            export_xlsx=False,
        ),
        marking_config=MarkingConfig(
            thresholds=QCThresholds(min_genes=20, pc_mt=30.0),
            plot_outliers=False,
            show_plots=False,
        ),
        doublet_config=DoubletConfig(run_algorithm=False),
        filter_config=FilterConfig(
            criteria_to_filter=["outlier_min_genes", "outlier_mt", "outlier_qc_metrics"],
            combination_logic="threshold",
            min_criteria_for_removal=2,
        ),
    )


def _fast_preprocess_config(n_top_genes: int = 200) -> PreprocessingWorkflowConfig:
    """Build a minimal preprocessing config for small synthetic data."""
    cfg = PreprocessingWorkflowConfig.quick(
        n_top_genes=n_top_genes,
        run_regression=False,
        run_integration=False,
        save_dir=None,
        n_jobs=1,
    )
    cfg.normalization.plot = False
    cfg.normalization.report = False
    cfg.graph.n_pcs = 15
    cfg.graph.n_neighbors = 10
    return cfg


def _fast_analysis_config() -> AnalysisWorkflowConfig:
    """Build a minimal analysis config with fixed resolution."""
    cfg = AnalysisWorkflowConfig(save_dir=None, n_jobs=1)
    cfg.clustering = ClusteringConfig(
        method="leiden",
        resolution=0.8,
        use_rep="X_pca",
        key_added="leiden_clusters",
        plot=False,
    )
    cfg.de = DifferentialConfig(
        groupby="leiden_clusters",
        method="wilcoxon",
        use_raw=True,
        key_added="rank_genes_groups",
    )
    cfg.annotation = AnnotationConfig(
        cluster_key="leiden_clusters",
        marker_species="human",
        run_celltypist=False,
        run_scoring=True,
        final_method="combined",
        key_added="cell_type_auto",
    )
    return cfg


def _run_minimal_pipeline(n_cells: int = 200, n_genes: int = 500) -> AnnData:
    """Run QC, preprocessing, and analysis on synthetic data."""
    import scLucid as scl

    adata = _synthetic_adata(n_cells=n_cells, n_genes=n_genes)
    adata = scl.run_pipeline(
        adata,
        stages=["qc", "preprocess", "analysis"],
        dataset_type="normal_tissue",
        species="human",
        qc_config=_fast_qc_config(),
        preprocess_config=_fast_preprocess_config(n_top_genes=min(200, n_genes)),
        analysis_config=_fast_analysis_config(),
        show_progress=False,
    )
    return adata


def test_preprocess_method_semantics_review_summary():
    """Preprocess review summary must document normalization, HVG, scaling,
    integration, layer transitions, and the interpretability of each layer."""
    adata = _run_minimal_pipeline()
    review = adata.uns["sclucid"]["preprocess"]["review_summary"]

    assert validate_review_summary_schema(review, module="preprocess").valid

    # Core parameter sections must exist.
    applied = review["applied_parameter_summary"]
    assert "normalization" in applied
    assert "hvg_selection" in applied
    assert "scaling" in applied
    assert "batch_correction" in applied

    # Normalization method and log-transform status.
    norm = applied["normalization"]
    assert "method" in norm
    assert norm["method"] is not None
    assert "log_transformed" in norm
    assert isinstance(norm["log_transformed"], bool)

    # HVG method and requested number.
    hvg = applied["hvg_selection"]
    assert "method" in hvg
    assert "requested_n_top_genes" in hvg
    assert hvg["requested_n_top_genes"] is not None

    # Scaling method and zero-center status.
    scaling = applied["scaling"]
    assert "method" in scaling
    assert scaling["method"] is not None
    assert "zero_center" in scaling
    assert isinstance(scaling["zero_center"], bool)

    # Integration method is documented even when disabled.
    batch = applied["batch_correction"]
    assert "method" in batch
    if batch["executed"]:
        assert batch["method"] is not None

    # Layer transition table with canonical counts -> normalized -> scaled flow.
    assert "layer_transition_table" in review
    table = review["layer_transition_table"]
    assert isinstance(table, list) and len(table) > 0
    stages = {row.get("canonical_stage") for row in table}
    assert {"counts", "normalized", "scaled"}.issubset(stages)

    # Explicit interpretability claims.
    decision_summary = review["preprocess_decision_summary"]
    scaling_decision = next(
        (d for d in decision_summary["decisions"] if d.get("step") == "scaling"), None
    )
    assert scaling_decision is not None
    assert "PCA/graph" in scaling_decision["risk_note"] or "expression-level" in scaling_decision["risk_note"]

    norm_policy = review["normalization_decision_policy"]
    downstream_note = norm_policy.get("downstream_note", "")
    assert "counts" in downstream_note and "pseudobulk" in downstream_note


def test_analysis_inference_semantics_review_summary():
    """Analysis review summary must document inference policy and claim levels
    for every DE and proportion result."""
    adata = _run_minimal_pipeline()
    review = adata.uns["sclucid"]["analysis"]["review_summary"]

    assert validate_review_summary_schema(review, module="analysis").valid

    # Required high-level sections.
    assert "analysis_inference_policy" in review
    assert "analysis_claim_level_summary" in review

    policy = review["analysis_inference_policy"]
    assert policy.get("claim_boundary") is not None
    assert "marker_discovery" in policy
    assert "condition_de" in policy

    claim_summary = review["analysis_claim_level_summary"]
    assert "outputs" in claim_summary
    outputs_list = claim_summary["outputs"]
    if isinstance(outputs_list, dict):
        outputs_list = list(outputs_list.values())
    outputs = {out["output"]: out for out in outputs_list}
    assert "cluster_marker_discovery" in outputs
    assert "condition_de" in outputs
    assert "cell_level_compare" in outputs
    assert "celltype_proportion" in outputs

    # Every DE result produced by the workflow must have a valid inference level
    # and cell-level results must carry a pseudoreplication warning.
    de_ns = adata.uns["sclucid"]["analysis"].get("de", {})
    allowed_de_levels = {
        "cell_level_marker_discovery",
        "sample_level",
        "exploratory_cell_level",
    }
    for key, value in de_ns.items():
        if key.endswith("_params"):
            continue
        if not isinstance(value, pd.DataFrame):
            continue
        if value.empty:
            continue
        assert "inference_level" in value.columns, f"DE result {key!r} missing inference_level"
        levels = set(value["inference_level"].dropna().unique())
        assert levels.issubset(allowed_de_levels), f"Unexpected levels in {key}: {levels}"
        if any(level in {"cell_level_marker_discovery", "exploratory_cell_level"} for level in levels):
            assert "pseudoreplication_warning" in value.columns, f"{key} missing pseudoreplication_warning"
            assert value["pseudoreplication_warning"].any()

    # Run a minimal proportion test so we can inspect proportion semantics.
    _add_synthetic_conditions(adata)
    prop_df, stat_df = analyze_celltype_proportion(
        adata,
        method="pseudobulk",
        sample_col="sampleID",
        condition_col="condition",
        celltype_col="cell_type_auto",
        config=ProportionConfig(
            celltype_col="cell_type_auto",
            sample_col="sampleID",
            condition_col="condition",
            test_method="clr-t-test",
            plot_types=[],
            out_dir=None,
        ),
    )
    assert not stat_df.empty
    assert "inference_level" in stat_df.columns
    allowed_prop_levels = {"sample_level", "exploratory_legacy_proportion"}
    assert set(stat_df["inference_level"].dropna().unique()).issubset(allowed_prop_levels)


def _add_synthetic_conditions(adata: AnnData) -> None:
    """Create two conditions with two samples each for closed-loop DE/proportion."""
    n = adata.n_obs
    rng = np.random.default_rng(7)
    conditions = np.where(rng.random(n) < 0.5, "A", "B")
    adata.obs["condition"] = conditions
    # Two samples per condition.
    sample_ids = [f"{cond}_{1 if rng.random() < 0.5 else 2}" for cond in conditions]
    adata.obs["sampleID"] = sample_ids


def test_pseudobulk_de_publication_inference_semantics():
    """A sample-level pseudobulk DE with replicates must be tagged as
    publication-ready at the sample level."""
    adata = _run_minimal_pipeline()
    _add_synthetic_conditions(adata)

    result = run_pseudobulk_de(
        adata,
        config=PseudobulkDEConfig(
            sample_col="sampleID",
            condition_key="condition",
            contrasts=[("A", "B")],
            layer="counts",
            method="auto",
            min_cells_per_sample=5,
            min_samples_per_condition=1,
        ),
    )

    assert not result.empty
    assert "inference_level" in result.columns
    assert set(result["inference_level"].unique()).issubset({"sample_level", "exploratory_cell_level"})
    assert result["valid_for_publication_inference"].any()
    assert any(
        level == "sample_level" and valid
        for level, valid in zip(result["inference_level"], result["valid_for_publication_inference"])
    )
