"""Tests for the analysis workflow."""

from unittest.mock import patch

import numpy as np
import pytest
import scanpy as sc

from scLucid.analysis.config import (
    AnalysisWorkflowConfig,
    AnnotationConfig,
    ClusteringConfig,
)
from scLucid.analysis.workflow import (
    WorkflowError,
    compare_clustering_resolutions,
    run_custom_analysis,
    run_standard_analysis,
)


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
        assert isinstance(
            passed_config, AnnotationConfig
        ), f"Expected AnnotationConfig, got {type(passed_config)}"
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


def test_run_standard_analysis_marker_step_requires_existing_cluster_key():
    """Marker-only runs should fail clearly when no cluster labels exist."""
    adata = _make_preprocessed_adata()
    workflow_config = AnalysisWorkflowConfig(
        clustering=ClusteringConfig(key_added="missing_clusters"),
        annotation=None,
        characterize=False,
    )

    with pytest.raises(WorkflowError, match="requires clustering results"):
        run_standard_analysis(
            adata,
            config=workflow_config,
            steps=["markers"],
            show_progress=False,
        )


class TestRunStandardAnalysisErrorHandling:
    """Test error recovery and resume behavior."""

    def test_error_recovery_requires_save_dir(self):
        """ValueError when error_recovery=True and on_error='save' without recovery_save_dir."""
        adata = _make_preprocessed_adata()
        with pytest.raises(ValueError, match="recovery_save_dir"):
            run_standard_analysis(
                adata,
                error_recovery=True,
                on_error="save",
            )

    def test_on_error_raise_propagates(self):
        """Default on_error='raise' propagates exceptions."""
        adata = _make_preprocessed_adata()
        with patch("scLucid.analysis.workflow.cluster_cells") as mock_cluster:
            mock_cluster.side_effect = RuntimeError("cluster fail")
            with pytest.raises(WorkflowError, match="cluster fail"):
                run_standard_analysis(adata)

    def test_workflow_stores_config(self):
        """Executed config stored in adata.uns."""
        adata = _make_preprocessed_adata()
        config = AnalysisWorkflowConfig.quick(run_annotation=False)
        result = run_standard_analysis(adata, config=config)
        assert "sclucid" in result.uns
        assert "workflow_config" in result.uns["sclucid"]["analysis"]
        assert "steps_executed" in result.uns["sclucid"]["analysis"]


class TestRunCustomAnalysis:
    """Test run_custom_analysis entry point."""

    def test_runs_clustering_step(self):
        """Custom analysis with clustering step only."""
        adata = _make_preprocessed_adata()
        result = run_custom_analysis(
            adata,
            steps=["clustering"],
            step_configs={"clustering": {"method": "leiden", "resolution": 0.8}},
        )
        assert "leiden_clusters" in result.obs.columns

    def test_runs_resolution_step(self):
        """Custom analysis with resolution search step."""
        adata = _make_preprocessed_adata()
        result = run_custom_analysis(
            adata,
            steps=["resolution"],
        )
        # Resolution search stores results in uns
        assert "sclucid" in result.uns

    def test_unknown_step_warns(self):
        """Unknown steps are warned and skipped."""
        adata = _make_preprocessed_adata()
        result = run_custom_analysis(adata, steps=["unknown_step"])
        # Should not raise, just warn
        assert isinstance(result, type(adata))

    def test_scoring_step_requires_gene_sets(self):
        """Scoring step warns and skips if no gene_sets provided."""
        adata = _make_preprocessed_adata()
        result = run_custom_analysis(adata, steps=["scoring"])
        assert isinstance(result, type(adata))

    def test_markers_step_prefers_leiden_clusters_default(self):
        """Marker workflow defaults to the repository's cluster key."""
        adata = _make_preprocessed_adata()
        adata.obs["leiden_clusters"] = np.where(np.arange(adata.n_obs) % 2 == 0, "0", "1").astype(
            object
        )
        adata.obs["leiden_clusters"] = adata.obs["leiden_clusters"].astype("category")

        with patch("scLucid.analysis.workflow.find_markers") as mock_find_markers:
            mock_find_markers.return_value = {"0": []}
            run_custom_analysis(adata, steps=["markers"])

            passed_config = mock_find_markers.call_args.args[1]
            assert passed_config.groupby == "leiden_clusters"

    def test_characterization_step_passes_save_dir_and_cluster_key(self, tmp_path):
        """Characterization workflow forwards save_dir and uses leiden_clusters by default."""
        adata = _make_preprocessed_adata()
        adata.obs["leiden_clusters"] = np.where(np.arange(adata.n_obs) % 2 == 0, "0", "1").astype(
            object
        )
        adata.obs["leiden_clusters"] = adata.obs["leiden_clusters"].astype("category")

        with patch("scLucid.analysis.workflow.characterize_clusters") as mock_characterize:
            mock_characterize.return_value = adata
            run_custom_analysis(
                adata,
                steps=["characterization"],
                save_dir=tmp_path,
            )

            kwargs = mock_characterize.call_args.kwargs
            assert kwargs["groupby"] == "leiden_clusters"
            assert kwargs["save_path"] == tmp_path


class TestCompareClusteringResolutions:
    """Test compare_clustering_resolutions utility."""

    def test_returns_dataframe(self):
        """Returns a DataFrame with expected columns."""
        adata = _make_preprocessed_adata(n_obs=100)
        df = compare_clustering_resolutions(
            adata,
            resolutions=[0.4, 0.8, 1.2],
            show_progress=False,
        )
        assert len(df) == 3
        assert "resolution" in df.columns
        assert "n_clusters" in df.columns

    def test_custom_metrics(self):
        """Only requested metrics are computed."""
        adata = _make_preprocessed_adata(n_obs=100)
        df = compare_clustering_resolutions(
            adata,
            resolutions=[0.5, 1.0],
            metrics=["n_clusters"],
            show_progress=False,
        )
        assert "n_clusters" in df.columns
        assert "silhouette" not in df.columns

    def test_saves_to_file(self, tmp_path):
        """Results saved when save_path provided."""
        adata = _make_preprocessed_adata(n_obs=100)
        out = tmp_path / "resolutions.csv"
        df = compare_clustering_resolutions(
            adata,
            resolutions=[0.5],
            save_path=out,
            show_progress=False,
        )
        assert out.exists()


class TestWorkflowKwargs:
    """Test kwargs override behavior in workflow functions."""

    def test_kwargs_override_config_values(self):
        """run_standard_analysis kwargs override top-level config fields."""
        adata = _make_preprocessed_adata()
        config = AnalysisWorkflowConfig.quick(run_annotation=False)
        # Override save_dir via kwargs (top-level field)
        with patch("scLucid.analysis.workflow.cluster_cells") as mock_cluster:
            patched = adata.copy()
            patched.obs["leiden_clusters"] = "0"
            mock_cluster.return_value = patched
            result = run_standard_analysis(adata, config=config, save_dir="./override")
            assert result.uns["sclucid"]["analysis"]["workflow_config"]["save_dir"] == "./override"
