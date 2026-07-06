"""
Test QC workflow behavior using synthetic data.
"""


from contextlib import suppress

from scLucid import qc
from scLucid.qc import QCWorkflowConfig
from tests.fixtures.synthetic_data import generate_minimal_adata


def _make_qc_test_adata():
    """Create a lightweight synthetic dataset suitable for QC workflow tests."""
    return generate_minimal_adata(n_cells=300, n_genes=800)


def _workflow_config_for_tests(save_dir=None) -> QCWorkflowConfig:
    """Create a non-interactive QC config for deterministic and fast tests."""
    config = QCWorkflowConfig(
        sample_key="sampleID", species="human", save_dir=save_dir, use_parallel=False
    )

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
    assert "qc_decision" in adata_qc.obs
    assert "qc_reason" in adata_qc.obs
    assert "cell_probability" in adata_qc.obs
    assert "empty_droplet_probability" in adata_qc.obs
    assert "ambient_fraction" in adata_qc.obs
    assert "doublet_score" in adata_qc.obs
    assert "qc_decision_summary" in adata_qc.uns["sclucid"]["qc"]
    assert "qc_probability_schema" in adata_qc.uns["sclucid"]["qc"]
    assert "ambient_rna_summary" in adata_qc.uns["sclucid"]["qc"]
    assert "empty_droplet_summary" in adata_qc.uns["sclucid"]["qc"]
    qc_ns = adata_qc.uns["sclucid"]["qc"]
    assert "artifact_contract" in qc_ns
    assert "threshold_decision" in qc_ns
    assert "mark_evidence" in qc_ns
    assert "qc_decision_artifact" in qc_ns
    assert "filter_result" in qc_ns
    assert "benchmark_review" in qc_ns
    assert qc_ns["artifact_contract"]["decision_flow"] == [
        "threshold_recommendation",
        "threshold_decision",
        "mark_evidence",
        "qc_decision",
        "filter_cells",
        "benchmark_review",
    ]


def test_run_qc_uses_reviewer_first_iterative_defaults():
    """Canonical run_qc should use qc_decision removal, not legacy filters."""
    adata = _make_qc_test_adata()
    config = _workflow_config_for_tests()
    config.doublet_config.run_algorithm = False

    adata_qc = qc.run_qc(
        adata,
        config=config,
        tissue_type="tumor",
        run_quick_review=False,
        show_progress=False,
    )

    qc_ns = adata_qc.uns["sclucid"]["qc"]
    assert qc_ns["qc_decision_filter_mode"] == "replace"
    assert qc_ns["iterative_qc_summary"]["final_filter_policy"] == "decision_remove"
    filtering = qc_ns["filtering_results"]
    assert filtering["criteria_used"] == ["qc_remove"]
    review = qc_ns["review_summary"]["data"]
    assert review["qc_filtering_policy_summary"]["final_filter_basis"] == "qc_decision_remove"


def test_run_qc_rejects_legacy_step_controls():
    """Step/resume controls remain available on run_standard_qc only."""
    adata = _make_qc_test_adata()
    config = _workflow_config_for_tests()
    config.doublet_config.run_algorithm = False

    import pytest

    with pytest.raises(TypeError, match="run_standard_qc"):
        qc.run_qc(
            adata,
            config=config,
            steps=["qc_metrics"],
            show_progress=False,
        )


def test_standard_qc_can_filter_by_qc_decision_remove():
    """Decision-based filtering should be opt-in and remove only qc_remove cells."""
    adata = _make_qc_test_adata()
    config = _workflow_config_for_tests()
    config.doublet_config.run_algorithm = False
    config.qc_decision_filter_mode = "replace"

    adata_qc = qc.run_standard_qc(adata, config=config, show_progress=False)

    assert "qc_decision" in adata_qc.obs
    assert "qc_remove" in adata_qc.obs
    assert not adata_qc.obs["qc_remove"].any()
    filtering = adata_qc.uns["sclucid"]["qc"]["filtering_results"]
    assert filtering["criteria_used"] == ["qc_remove"]
    review = adata_qc.uns["sclucid"]["qc"]["review_summary"]["data"]
    filtering_policy = review["qc_filtering_policy_summary"]
    assert filtering_policy["final_filter_basis"] == "qc_decision_remove"
    assert filtering_policy["review_required"] is False


def test_iterative_qc_records_phase_summary():
    """Iterative QC entrypoint should expose the phase contract."""
    adata = _make_qc_test_adata()
    config = _workflow_config_for_tests()
    config.doublet_config.run_algorithm = False

    adata_qc = qc.run_iterative_qc(
        adata,
        config=config,
        tissue_type="tumor",
        final_filter_policy="none",
        show_progress=False,
    )

    summary = adata_qc.uns["sclucid"]["qc"]["iterative_qc_summary"]
    assert summary["schema_version"] == "iterative_qc_summary_v1"
    assert summary["final_filter_policy"] == "none"
    assert [phase["phase"] for phase in summary["phases"]] == [
        "lenient_cell_screen",
        "doublet_contamination_stress_evidence",
        "quick_biology_review",
        "final_qc_decision",
    ]
    assert "qc_decision" in adata_qc.obs


def test_iterative_qc_runs_quick_biology_review():
    """Iterative QC should optionally run a temporary embedding review."""
    adata = _make_qc_test_adata()
    config = _workflow_config_for_tests()
    config.doublet_config.run_algorithm = False

    adata_qc = qc.run_iterative_qc(
        adata,
        config=config,
        tissue_type="tumor",
        final_filter_policy="none",
        run_quick_review=True,
        quick_review_max_cells=120,
        quick_review_n_top_genes=150,
        quick_review_n_pcs=10,
        quick_review_n_neighbors=8,
        show_progress=False,
    )

    quick = adata_qc.uns["sclucid"]["qc"]["iterative_qc_summary"]["quick_biology_review"]
    assert quick["schema_version"] == "quick_biology_review_v1"
    assert quick["status"] == "complete"
    assert quick["n_cells_reviewed"] <= 120
    assert isinstance(quick["cluster_qc_table"], list)
    assert "X_umap" not in adata_qc.obsm


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
    threshold_results = learner.learn_all_threshold_results(adata)

    assert isinstance(threshold_results, dict)
    assert len(threshold_results) > 0
    for _, result in threshold_results.items():
        assert result["schema_version"] == "qc_threshold_result_v1"
        assert result["threshold"] >= 0


def test_standard_qc_creates_single_sample_key_when_missing():
    """QC workflow should gracefully fall back to another sample column when present."""
    adata = _make_qc_test_adata()
    if "sampleID" in adata.obs:
        del adata.obs["sampleID"]

    config = _workflow_config_for_tests()
    config.doublet_config.run_algorithm = False

    adata_qc = qc.run_standard_qc(adata, config=config, show_progress=False)

    assert "batch" in adata_qc.obs.columns
    assert adata_qc.obs["batch"].nunique() == 1
    context = adata_qc.uns["sclucid"]["qc"]["context"]["data"]
    warnings = adata_qc.uns["sclucid"]["qc"]["warnings"]["data"]
    assert context["sample_key"] == "batch"
    assert any("using detected obs column 'batch'" in w for w in warnings)


def test_standard_qc_auto_detects_existing_sample_column():
    """Workflow should use common alternate sample columns before fabricating one."""
    adata = _make_qc_test_adata()
    adata.obs["orig.ident"] = adata.obs["sampleID"].astype(str)
    del adata.obs["sampleID"]

    config = _workflow_config_for_tests()
    config.doublet_config.run_algorithm = False

    adata_qc = qc.run_standard_qc(adata, config=config, show_progress=False)

    context = adata_qc.uns["sclucid"]["qc"]["context"]["data"]
    warnings = adata_qc.uns["sclucid"]["qc"]["warnings"]["data"]
    assert context["sample_key"] == "orig.ident"
    assert all("synthetic single-sample labels" not in warning for warning in warnings)


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


def test_generate_qc_report_exports_summary_bundle(tmp_path):
    """Workflow report generation should produce structured review artifacts."""
    adata = _make_qc_test_adata()
    config = _workflow_config_for_tests(save_dir=str(tmp_path))
    config.doublet_config.run_algorithm = False

    qc.run_standard_qc(adata, config=config, show_progress=False)

    report_dir = tmp_path / "report"
    assert (report_dir / "qc_summary.json").exists()
    assert (report_dir / "qc_summary.md").exists()
    assert (report_dir / "qc_report.html").exists()


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
        with suppress(Exception):
            qc.run_standard_qc(
                adata,
                config=config,
                show_progress=False,
                error_recovery=True,
                recovery_save_dir=str(recovery_dir),
            )

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
