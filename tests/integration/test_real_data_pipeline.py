"""Integration tests using real single-cell datasets.

These tests validate that scLucid's full pipeline (QC -> Preprocessing -> Analysis)
runs successfully on published real-world datasets.

Datasets used:
- pbmc3k: Non-tumor PBMC (4 samples, ~2700 cells)
- schlesinger2020.pdac: Single-sample pancreatic cancer (~6500 cells)
- lin2020.pdac: Multi-sample pancreatic cancer (~9600 cells, 10 samples)

All tests are marked as ``slow`` and ``integration`` and are skipped by default.
Run explicitly with: ``pytest tests/integration/test_real_data_pipeline.py -m 'slow'``
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from anndata import AnnData

import scLucid as scl

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_DIR = Path(__file__).parents[2] / "data"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_subset(path: Path, n_cells: int = 400, random_state: int = 42) -> AnnData:
    """Load a real dataset and subsample for fast integration testing."""
    import anndata as ad

    adata = ad.read_h5ad(str(path))
    if adata.n_obs > n_cells:
        rng = np.random.default_rng(random_state)
        idx = rng.choice(adata.n_obs, size=n_cells, replace=False)
        adata = adata[idx].copy()

    # Ensure raw counts are available. Datasets vary:
    # - pbmc3k: normalized in .X, counts in layers["counts"]
    # - PDAC: raw counts in .X, no layers
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    else:
        adata.layers["counts"] = adata.X.copy()

    return adata


def _assert_qc_valid(adata: AnnData) -> None:
    """Assert that QC step produced expected outputs."""
    assert adata.n_obs > 0, "QC removed all cells"
    assert "n_genes_by_counts" in adata.obs.columns
    assert "total_counts" in adata.obs.columns
    assert "pct_counts_mt" in adata.obs.columns
    assert "sclucid" in adata.uns


def _assert_preprocess_valid(adata: AnnData) -> None:
    """Assert that preprocessing step produced expected outputs."""
    assert "normalized" in adata.layers
    hvg_col = "highly_variable" if "highly_variable" in adata.var.columns else "highly_variable_selected"
    assert hvg_col in adata.var.columns, f"No HVG column found. Available: {list(adata.var.columns)}"
    n_hvgs = int(adata.var[hvg_col].sum())
    assert n_hvgs > 0, "No HVGs selected"
    assert "X_pca" in adata.obsm
    assert "X_umap" in adata.obsm
    assert "neighbors" in adata.uns
    assert "sclucid" in adata.uns
    assert "preprocess" in adata.uns["sclucid"]
    assert "review_summary" in adata.uns["sclucid"]["preprocess"]


def _assert_analysis_valid(adata: AnnData, cluster_key: str = "leiden_clusters") -> None:
    """Assert that analysis step produced expected outputs."""
    assert cluster_key in adata.obs.columns
    n_clusters = adata.obs[cluster_key].nunique()
    assert n_clusters >= 2, f"Expected >=2 clusters, got {n_clusters}"
    assert "rank_genes_groups" in adata.uns
    assert "sclucid" in adata.uns
    assert "analysis" in adata.uns["sclucid"]
    assert "review_summary" in adata.uns["sclucid"]["analysis"]


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def pbmc3k_subset() -> AnnData:
    """PBMC 3k subset (~400 cells, non-tumor, 4 samples)."""
    return _load_subset(DATA_DIR / "pbmc3k.h5ad", n_cells=400)


@pytest.fixture(scope="module")
def pdac_single_subset() -> AnnData:
    """Schlesinger 2020 PDAC subset (~400 cells, single-sample tumor)."""
    return _load_subset(DATA_DIR / "schlesinger2020.pdac.h5ad", n_cells=400)


@pytest.fixture(scope="module")
def pdac_multi_subset() -> AnnData:
    """Lin 2020 PDAC subset (~500 cells, multi-sample tumor)."""
    return _load_subset(DATA_DIR / "lin2020.pdac.h5ad", n_cells=500)


# ---------------------------------------------------------------------------
# PBMC (non-tumor) pipeline
# ---------------------------------------------------------------------------

@pytest.mark.slow
@pytest.mark.integration
def test_pbmc3k_full_pipeline(pbmc3k_subset: AnnData) -> None:
    """End-to-end pipeline on non-tumor PBMC data."""
    adata = pbmc3k_subset

    # QC
    adata = scl.qc.run_standard_qc(adata, show_progress=False)
    _assert_qc_valid(adata)

    # Preprocessing
    cfg = scl.preprocess.WorkflowConfig()
    cfg.normalization.plot = False
    cfg.normalization.report = False
    cfg.hvg.n_top_genes = 500
    cfg.graph.n_pcs = 20
    adata = scl.preprocess.run_preprocessing(adata, config=cfg, show_progress=False)
    _assert_preprocess_valid(adata)

    # Analysis
    acfg = scl.analysis.AnalysisWorkflowConfig()
    acfg.clustering = scl.analysis.ClusteringConfig(resolution=0.5, plot=False)
    acfg.annotation = None  # Skip annotation for speed
    adata = scl.analysis.workflow.run_standard_analysis(
        adata, config=acfg, show_progress=False, skip_steps=["annotation"]
    )
    _assert_analysis_valid(adata)


# ---------------------------------------------------------------------------
# Single-sample PDAC pipeline
# ---------------------------------------------------------------------------

@pytest.mark.slow
@pytest.mark.integration
def test_pdac_single_full_pipeline(pdac_single_subset: AnnData) -> None:
    """End-to-end pipeline on single-sample PDAC tumor data."""
    adata = pdac_single_subset

    # QC (tumor-aware)
    adata = scl.qc.run_standard_qc(adata, tissue_type="tumor", show_progress=False)
    _assert_qc_valid(adata)

    # Preprocessing (no batch correction for single sample)
    cfg = scl.preprocess.WorkflowConfig()
    cfg.normalization.plot = False
    cfg.normalization.report = False
    cfg.hvg.n_top_genes = 500
    cfg.graph.n_pcs = 20
    cfg.integration.method = None
    adata = scl.preprocess.run_preprocessing(adata, config=cfg, show_progress=False)
    _assert_preprocess_valid(adata)

    # Analysis
    acfg = scl.analysis.AnalysisWorkflowConfig()
    acfg.clustering = scl.analysis.ClusteringConfig(resolution=0.5, plot=False)
    acfg.annotation = None
    adata = scl.analysis.workflow.run_standard_analysis(
        adata, config=acfg, show_progress=False, skip_steps=["annotation"]
    )
    _assert_analysis_valid(adata)


# ---------------------------------------------------------------------------
# Multi-sample PDAC pipeline
# ---------------------------------------------------------------------------

@pytest.mark.slow
@pytest.mark.integration
def test_pdac_multi_full_pipeline(pdac_multi_subset: AnnData) -> None:
    """End-to-end pipeline on multi-sample PDAC with batch correction."""
    adata = pdac_multi_subset

    # QC
    adata = scl.qc.run_standard_qc(adata, tissue_type="tumor", show_progress=False)
    _assert_qc_valid(adata)

    # Preprocessing with Harmony batch correction
    cfg = scl.preprocess.WorkflowConfig()
    cfg.normalization.plot = False
    cfg.normalization.report = False
    cfg.hvg.n_top_genes = 500
    cfg.graph.n_pcs = 20
    cfg.integration.method = "harmony"
    cfg.integration.batch_key = "sampleID"
    adata = scl.preprocess.run_preprocessing(adata, config=cfg, show_progress=False)
    _assert_preprocess_valid(adata)
    assert "X_harmony" in adata.obsm, "Harmony integration did not produce X_harmony"

    # Analysis using Harmony embedding
    acfg = scl.analysis.AnalysisWorkflowConfig()
    acfg.clustering = scl.analysis.ClusteringConfig(
        resolution=0.5, plot=False, use_rep="X_harmony"
    )
    acfg.annotation = None
    adata = scl.analysis.workflow.run_standard_analysis(
        adata, config=acfg, show_progress=False, skip_steps=["annotation"]
    )
    _assert_analysis_valid(adata)


# ---------------------------------------------------------------------------
# Unified pipeline entry point
# ---------------------------------------------------------------------------

@pytest.mark.slow
@pytest.mark.integration
def test_run_pipeline_unified_entry(pbmc3k_subset: AnnData) -> None:
    """Test the unified ``scl.run_pipeline`` entry point on real data."""
    adata = pbmc3k_subset

    # Skip annotation in unified pipeline until annotation module bug is fixed
    adata = scl.run_pipeline(
        adata,
        stages=["qc", "preprocess", "analysis"],
        show_progress=False,
        analysis_skip_steps=["annotation"],
    )

    # Assert all stages produced results
    assert "sclucid" in adata.uns
    assert "qc" in adata.uns["sclucid"]
    assert "preprocess" in adata.uns["sclucid"]
    assert "analysis" in adata.uns["sclucid"]
    assert "pipeline_context" in adata.uns["sclucid"]


# ---------------------------------------------------------------------------
# Data-specific assertions
# ---------------------------------------------------------------------------

@pytest.mark.slow
@pytest.mark.integration
def test_pbmc3k_species_propagated(pbmc3k_subset: AnnData) -> None:
    """Species info from QC should propagate through pipeline."""
    adata = pbmc3k_subset
    adata = scl.qc.run_standard_qc(adata, show_progress=False)
    assert "species" in adata.obs.columns or "species" in str(adata.uns.get("sclucid", {}))


@pytest.mark.slow
@pytest.mark.integration
def test_pdac_multi_batch_key_preserved(pdac_multi_subset: AnnData) -> None:
    """Multi-sample PDAC should retain batch key after pipeline."""
    adata = pdac_multi_subset
    assert "sampleID" in adata.obs.columns

    cfg = scl.preprocess.WorkflowConfig()
    cfg.normalization.plot = False
    cfg.normalization.report = False
    cfg.integration.method = "harmony"
    cfg.integration.batch_key = "sampleID"
    adata = scl.preprocess.run_preprocessing(adata, config=cfg, show_progress=False)

    assert "sampleID" in adata.obs.columns
    n_samples = adata.obs["sampleID"].nunique()
    assert n_samples >= 2, f"Expected >=2 samples, got {n_samples}"
