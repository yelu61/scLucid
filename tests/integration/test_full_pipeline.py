"""Integration tests for core scLucid pipeline entry points."""

import pytest

from tests.fixtures.synthetic_data import generate_minimal_adata


@pytest.mark.integration
def test_qc_metrics_pipeline_step():
    """QC metrics should run on synthetic data without cell-cycle scoring."""
    from scLucid.qc.metrics import calculate_qc_metric
    from scLucid.qc.config import MetricsReportingConfig

    adata = generate_minimal_adata(n_cells=120, n_genes=300)
    adata.obs["sampleID"] = "sample_1"

    reporting_cfg = MetricsReportingConfig(
        plot_violin=False,
        plot_scatter=False,
        plot_top_genes=False,
        show_plots=False,
        export_stats=False,
        export_xlsx=False,
    )

    out = calculate_qc_metric(
        adata,
        sample_key="sampleID",
        reporting_config=reporting_cfg,
        calculate_cell_cycle=False,
    )

    assert "n_genes_by_counts" in out.obs
    assert "total_counts" in out.obs
    assert "pct_counts_in_top_20_genes" in out.obs


@pytest.mark.integration
def test_preprocessing_normalization_step():
    """Preprocessing normalization workflow step should populate normalized layer."""
    from scLucid.preprocess.workflow import run_preprocessing
    from scLucid.preprocess.config import WorkflowConfig

    adata = generate_minimal_adata(n_cells=120, n_genes=300)

    cfg = WorkflowConfig()
    cfg.normalization.plot = False
    cfg.normalization.report = False

    out = run_preprocessing(
        adata,
        config=cfg,
        steps=["normalization"],
        show_progress=False,
    )

    assert "normalized" in out.layers


@pytest.mark.integration
def test_analysis_clustering_pipeline_step():
    """Analysis clustering should run once PCA/neighbors are available."""
    sc = pytest.importorskip("scanpy")

    from scLucid.analysis.clustering import cluster_cells
    from scLucid.analysis.config import ClusteringConfig

    adata = generate_minimal_adata(n_cells=120, n_genes=300)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)

    cfg = ClusteringConfig(
        method="leiden",
        resolution=0.5,
        key_added="leiden_clusters",
        plot=False,
    )
    out = cluster_cells(adata, config=cfg)

    assert "leiden_clusters" in out.obs
