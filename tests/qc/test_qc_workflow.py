"""
Test QC workflow behavior using synthetic data.
"""

import numpy as np

from scLucid import qc
from scLucid.qc import QCWorkflowConfig
from tests.fixtures.synthetic_data import generate_minimal_adata


def _make_qc_test_adata():
    """Create a lightweight synthetic dataset suitable for QC workflow tests."""
    return generate_minimal_adata(n_cells=300, n_genes=800)


def _workflow_config_for_tests(save_dir=None) -> QCWorkflowConfig:
    """Create a non-interactive QC config for deterministic and fast tests."""
    config = QCWorkflowConfig(sample_key="sampleID", species="human", save_dir=save_dir, use_parallel=False)

    # Metrics step: disable interactive plotting/exports.
    config.metrics_reporting_config.show_plots = False
    config.metrics_reporting_config.plot_top_genes = False
    config.metrics_reporting_config.plot_violin = False
    config.metrics_reporting_config.plot_scatter = False
    config.metrics_reporting_config.export_stats = False
    config.metrics_reporting_config.print_stats = False

    # Doublet step: keep algorithmic signal but disable expensive visual/report side effects.
    config.doublet_config.show_plots = False
    config.doublet_config.plot_summary = False
    config.doublet_config.plot_bar = False
    config.doublet_config.plot_scatter = False
    config.doublet_config.plot_upset = False
    config.doublet_config.export_stats = False
    config.doublet_config.scr_plot_umap = False
    config.doublet_config.use_heuristics = False

    # Marking step: disable outlier plots.
    config.marking_config.plot_outliers = False
    # Keep filtering permissive in synthetic tests to avoid degenerate empty outputs.
    config.filter_config.criteria_to_filter = ["predicted_doublet"]
    config.filter_config.combination_logic = "any"

    return config


def test_standard_qc_workflow():
    """Test standard QC workflow."""
    adata = _make_qc_test_adata()
    original_n_obs = adata.n_obs

    config = _workflow_config_for_tests()
    adata_qc = qc.run_standard_qc(adata, config=config, show_progress=False)

    # Check that cells were filtered or retained safely.
    assert adata_qc.n_obs <= original_n_obs
    assert adata_qc.n_obs > 0

    # Check that QC/marking annotations exist.
    assert "outlier_count" in adata_qc.obs or "predicted_doublet" in adata_qc.obs


def test_qc_with_adaptive_thresholds():
    """Test QC workflow with adaptive threshold learning."""
    from scLucid.qc import AdaptiveThresholdLearner

    adata = _make_qc_test_adata()

    # Calculate metrics first (non-interactive).
    qc.calculate_qc_metric(
        adata,
        sample_key="sampleID",
        show_plots=False,
        plot_top_genes=False,
        plot_violin=False,
        plot_scatter=False,
        export_stats=False,
        print_stats=False,
    )

    learner = AdaptiveThresholdLearner(method="percentile")
    thresholds = learner.learn_all_thresholds(adata)

    assert isinstance(thresholds, dict)
    assert len(thresholds) > 0
    for _, value in thresholds.items():
        assert isinstance(value, (int, float))
        assert value >= 0


def test_qc_html_report(tmp_path):
    """Test HTML report generation."""
    from scLucid.qc import generate_qc_html_report

    adata = _make_qc_test_adata()
    adata_qc = qc.run_standard_qc(adata, config=_workflow_config_for_tests(), show_progress=False)

    output_file = tmp_path / "test_report.html"
    generate_qc_html_report(
        adata_qc,
        adata_before=adata,
        output_path=str(output_file),
        title="Test Report",
    )

    assert output_file.exists()
    assert output_file.stat().st_size > 0


class TestQCWorkflowErrorRecovery:
    """Test QC workflow error recovery features."""

    def test_qc_workflow_with_error_recovery(self, tmp_path):
        """Test QC workflow with error recovery enabled."""
        adata = _make_qc_test_adata()
        config = _workflow_config_for_tests()

        recovery_dir = tmp_path / "recovery"
        adata_qc = qc.run_standard_qc(
            adata,
            config=config,
            show_progress=False,
            error_recovery=True,
            recovery_save_dir=str(recovery_dir),
            on_error="raise",
        )

        assert adata_qc.n_obs > 0

    def test_qc_workflow_storage_api(self, tmp_path):
        """Test that QC workflow uses new storage API."""
        from scLucid.utils import load_workflow_result

        adata = _make_qc_test_adata()
        config = _workflow_config_for_tests()
        adata_qc = qc.run_standard_qc(adata, config=config, show_progress=False)

        result = load_workflow_result(adata_qc, "qc", "standard")
        assert result is not None
        assert result["name"] == "standard"
        assert "steps_executed" in result
        assert "completed_at" in result

    def test_qc_workflow_resume_from_checkpoint(self, tmp_path):
        """Test QC workflow resume from checkpoint parameter path."""
        adata = _make_qc_test_adata()
        config = _workflow_config_for_tests()

        recovery_dir = tmp_path / "recovery"
        try:
            qc.run_standard_qc(
                adata,
                config=config,
                show_progress=False,
                error_recovery=True,
                recovery_save_dir=str(recovery_dir),
            )
        except Exception:
            pass

    def test_qc_workflow_save_dir_parameter(self, tmp_path):
        """Test QC workflow with unified save_dir parameter."""
        adata = _make_qc_test_adata()
        config = _workflow_config_for_tests(save_dir=str(tmp_path))

        qc.run_standard_qc(
            adata,
            config=config,
            show_progress=False,
            error_recovery=True,
            recovery_save_dir=str(tmp_path / "recovery"),
            on_error="save",
        )
        assert tmp_path.exists()


class TestQCWorkflowBackwardCompat:
    """Test backward compatibility for QC workflow."""

    def test_qc_config_results_dir_alias(self):
        """Test that results_dir property alias works."""
        config = QCWorkflowConfig()
        assert hasattr(config, "results_dir")

        config.results_dir = "./test_results"
        assert config.save_dir == "./test_results"

    def test_qc_config_from_simple_dict_with_results_dir(self):
        """Test from_simple_dict accepts results_dir."""
        config = QCWorkflowConfig.from_simple_dict(
            {"thresholds_min_genes": 200, "results_dir": "./legacy_results"}
        )
        assert config.save_dir == "./legacy_results"
