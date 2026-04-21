"""Tests for the analysis workflow."""

from unittest.mock import patch

import numpy as np
import scanpy as sc

from scLucid.analysis.workflow import run_standard_analysis
from scLucid.analysis.config import AnalysisWorkflowConfig, AnnotationConfig, ClusteringConfig


def _make_preprocessed_adata(n_obs=200, n_vars=500):
    """Create minimal preprocessed AnnData suitable for analysis."""
    import anndata
    counts = np.random.poisson(5, size=(n_obs, n_vars)).astype(np.float32)
    adata = anndata.AnnData(X=counts)
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.var_names = [f"gene_{i}" for i in range(n_vars)]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=50, flavor="seurat")
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata)
    adata.raw = adata
    return adata


def test_run_standard_analysis_passes_annotation_config():
    """Verify that a Pydantic AnnotationConfig object is passed through to run_annotation."""
    adata = _make_preprocessed_adata()
    custom_annotation = AnnotationConfig(
        cluster_key="leiden_clusters",
        key_added="my_cell_type",
        run_scoring=False,
        final_method="celltypist",
        celltypist_model="Immune_All_Low.pkl",
    )
    workflow_config = AnalysisWorkflowConfig(
        clustering=ClusteringConfig(key_added="leiden_clusters"),
        annotation=custom_annotation,
    )

    with patch("scLucid.analysis.workflow.run_annotation") as mock_run_annotation:
        mock_run_annotation.return_value = adata
        _ = run_standard_analysis(adata, config=workflow_config)

        assert mock_run_annotation.called, "run_annotation was not called"
        _, kwargs = mock_run_annotation.call_args
        passed_config = kwargs.get("config")
        assert isinstance(passed_config, AnnotationConfig), (
            f"Expected AnnotationConfig, got {type(passed_config)}"
        )
        assert passed_config.key_added == "my_cell_type"
        assert passed_config.run_scoring is False
        assert passed_config.final_method == "celltypist"


def test_run_standard_analysis_passes_dict_annotation_config():
    """Verify that a dict-based annotation config is converted and passed through."""
    adata = _make_preprocessed_adata()
    workflow_config = AnalysisWorkflowConfig(
        clustering=ClusteringConfig(key_added="leiden_clusters"),
        annotation={
            "key_added": "dict_cell_type",
            "run_scoring": False,
            "final_method": "celltypist",
        },
    )

    with patch("scLucid.analysis.workflow.run_annotation") as mock_run_annotation:
        mock_run_annotation.return_value = adata
        _ = run_standard_analysis(adata, config=workflow_config)

        assert mock_run_annotation.called
        _, kwargs = mock_run_annotation.call_args
        passed_config = kwargs.get("config")
        assert isinstance(passed_config, AnnotationConfig)
        assert passed_config.key_added == "dict_cell_type"
