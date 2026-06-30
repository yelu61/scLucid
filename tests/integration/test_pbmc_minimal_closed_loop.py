"""Minimal PBMC golden-path closed-loop integration test.

This test exercises the full QC -> preprocess -> analysis pipeline on a
subsampled PBMC3k dataset and verifies that downstream publication-aware
DE and proportion analyses can be run.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import scanpy as sc

import scLucid as scl
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


REPO_ROOT = Path(__file__).parents[2]
DATA_PATH = REPO_ROOT / "data" / "pbmc3k.h5ad"


@pytest.mark.integration
def test_pbmc_minimal_closed_loop(tmp_path):
    """Run a fast PBMC pipeline and verify closed-loop comparative analysis."""
    assert DATA_PATH.exists(), f"PBMC test data not found at {DATA_PATH}"

    adata = sc.read_h5ad(str(DATA_PATH))
    if adata.n_obs > 300:
        rng = np.random.default_rng(42)
        idx = rng.choice(adata.n_obs, size=300, replace=False)
        adata = adata[idx].copy()
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    else:
        adata.layers["counts"] = adata.X.copy()

    initial_cells = adata.n_obs

    qc_config = QCWorkflowConfig(
        save_dir=None,
        species="human",
        tissue_type="normal_tissue",
        use_recommendations=False,
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
            thresholds=QCThresholds(min_genes=700, pc_mt=20.0),
            plot_outliers=False,
            show_plots=False,
        ),
        doublet_config=DoubletConfig(run_algorithm=False),
        filter_config=FilterConfig(
            criteria_to_filter=["outlier_min_genes"],
            combination_logic="threshold",
            min_criteria_for_removal=1,
        ),
    )

    preprocess_config = PreprocessingWorkflowConfig.quick(
        n_top_genes=500,
        run_regression=False,
        run_integration=False,
        save_dir=None,
        n_jobs=1,
    )
    preprocess_config.normalization.plot = False
    preprocess_config.normalization.report = False
    preprocess_config.graph.n_pcs = 15
    preprocess_config.graph.n_neighbors = 10

    analysis_config = AnalysisWorkflowConfig(save_dir=None, n_jobs=1)
    analysis_config.clustering = ClusteringConfig(
        method="leiden",
        resolution=0.8,
        use_rep="X_pca",
        key_added="leiden_clusters",
        plot=False,
    )
    analysis_config.de = DifferentialConfig(
        groupby="leiden_clusters",
        method="wilcoxon",
        use_raw=True,
        key_added="rank_genes_groups",
    )
    analysis_config.annotation = AnnotationConfig(
        cluster_key="leiden_clusters",
        marker_species="human",
        run_celltypist=False,
        run_scoring=True,
        final_method="combined",
        key_added="cell_type_auto",
    )

    adata = scl.run_pipeline(
        adata,
        stages=["qc", "preprocess", "analysis"],
        dataset_type="pbmc_or_blood",
        species="human",
        qc_config=qc_config,
        preprocess_config=preprocess_config,
        analysis_config=analysis_config,
        show_progress=False,
    )

    # --- QC retention sanity check ---
    retention = adata.n_obs / initial_cells
    assert 0.50 <= retention <= 0.98, f"Retention {retention:.2f} outside expected range"

    # --- Major PBMC cell types are present ---
    assert "cell_type_auto" in adata.obs.columns
    labels = set(adata.obs["cell_type_auto"].astype(str).str.lower().unique())
    major_types = ["t cell", "b cell", "monocyte", "nk", "cd4", "cd8"]
    found = [ct for ct in major_types if any(ct in label for label in labels)]
    assert len(found) >= 2, f"Expected at least 2 major PBMC types, found {found} in {labels}"

    # --- Sample-level pseudobulk DE between two arbitrary groups ---
    rng = np.random.default_rng(7)
    adata.obs["condition"] = np.where(rng.random(adata.n_obs) < 0.5, "A", "B")
    adata.obs["sampleID"] = [
        f"{cond}_{1 if rng.random() < 0.5 else 2}" for cond in adata.obs["condition"]
    ]

    de_result = run_pseudobulk_de(
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
    assert not de_result.empty
    assert de_result["valid_for_publication_inference"].any()
    assert any(
        level == "sample_level" and valid
        for level, valid in zip(de_result["inference_level"], de_result["valid_for_publication_inference"])
    )

    # --- Sample-level proportion test between the same groups ---
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
    assert "sample_level" in stat_df["inference_level"].values
    assert stat_df["valid_for_publication_inference"].any()

    # --- Required review sections under adata.uns["sclucid"] ---
    for stage in ("qc", "preprocess", "analysis"):
        assert stage in adata.uns["sclucid"], f"Missing sclucid namespace for {stage}"
        assert "review_summary" in adata.uns["sclucid"][stage], f"Missing review_summary for {stage}"

    # Sanity check that the analysis review summary records inference metadata.
    analysis_review = adata.uns["sclucid"]["analysis"]["review_summary"]
    assert "analysis_inference_policy" in analysis_review
    assert "analysis_claim_level_summary" in analysis_review
