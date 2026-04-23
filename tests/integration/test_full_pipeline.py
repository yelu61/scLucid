"""Integration tests for core scLucid pipeline entry points."""

import pytest

from tests.fixtures.synthetic_data import generate_minimal_adata


@pytest.mark.integration
def test_qc_metrics_pipeline_step():
    """QC metrics should run on synthetic data without cell-cycle scoring."""
    from scLucid.qc.config import MetricsReportingConfig
    from scLucid.qc.metrics import calculate_qc_metric

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
    from scLucid.preprocess.config import WorkflowConfig
    from scLucid.preprocess.workflow import run_preprocessing

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


@pytest.mark.integration
def test_preprocess_review_summary_created():
    """Preprocessing should create a review summary in adata.uns."""
    from scLucid.preprocess.config import WorkflowConfig
    from scLucid.preprocess.workflow import run_preprocessing

    adata = generate_minimal_adata(n_cells=120, n_genes=300)

    cfg = WorkflowConfig()
    cfg.normalization.plot = False
    cfg.normalization.report = False
    cfg.hvg.n_top_genes = 100

    out = run_preprocessing(
        adata,
        config=cfg,
        steps=["normalization", "hvg_selection", "subset_hvg"],
        show_progress=False,
    )

    assert "sclucid" in out.uns
    assert "preprocess" in out.uns["sclucid"]
    assert "review_summary" in out.uns["sclucid"]["preprocess"]

    summary = out.uns["sclucid"]["preprocess"]["review_summary"]
    assert "steps_executed" in summary
    assert "data_shape" in summary
    assert summary["steps_executed"] == ["normalization", "hvg_selection", "subset_hvg"]


@pytest.mark.integration
def test_run_pipeline_skips_qc_then_preprocess_fails():
    """Running preprocess without QC should raise a clear error."""
    import scLucid as scl

    adata = generate_minimal_adata(n_cells=120, n_genes=300)
    # Remove any pre-existing sclucid results to simulate fresh data
    if "sclucid" in adata.uns:
        del adata.uns["sclucid"]

    with pytest.raises(RuntimeError, match="QC was skipped"):
        scl.run_pipeline(adata, stages=["preprocess"])


@pytest.mark.integration
def test_run_pipeline_skips_preprocess_then_analysis_fails():
    """Running analysis without preprocessing should raise a clear error."""
    import scLucid as scl

    adata = generate_minimal_adata(n_cells=120, n_genes=300)
    if "sclucid" in adata.uns:
        del adata.uns["sclucid"]

    with pytest.raises(RuntimeError, match="preprocessing was skipped"):
        scl.run_pipeline(adata, stages=["analysis"])


@pytest.mark.integration
def test_run_pipeline_strips_stage_kwarg_prefixes(monkeypatch):
    """Stage-prefixed kwargs should be forwarded without their prefixes."""
    import scLucid as scl

    adata = generate_minimal_adata(n_cells=120, n_genes=300)
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
        qc_steps=["qc_metrics"],
        qc_skip_steps=None,
        show_progress=False,
    )

    assert out is adata
    assert seen["steps"] == ["qc_metrics"]
    assert seen["skip_steps"] is None
    assert "qc_steps" not in seen
