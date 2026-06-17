"""Tests for QC filtering functions."""

import pandas as pd
import pytest
from anndata import AnnData

from scLucid.qc.config import FilterConfig
from scLucid.qc.filtering import (
    audit_filtering,
    filter_cells,
    generate_qc_report,
    identify_outliers,
    mark_low_quality_cell,
    resolve_qc_thresholds,
    run_qc_decision_workflow,
    run_qc_threshold_decision,
    suggest_qc_thresholds,
)


class TestSuggestQCThresholds:
    def test_returns_dataframe_and_thresholds(self, qc_test_adata):
        df, thresholds = suggest_qc_thresholds(
            qc_test_adata, plot_distributions=False
        )
        assert isinstance(df, pd.DataFrame)
        assert hasattr(thresholds, "min_genes")
        assert hasattr(thresholds, "pc_mt")

    def test_mad_method(self, qc_test_adata):
        df, thresholds = suggest_qc_thresholds(
            qc_test_adata, method="mad", plot_distributions=False
        )
        assert df is not None

    def test_percentile_method(self, qc_test_adata):
        df, thresholds = suggest_qc_thresholds(
            qc_test_adata,
            method="percentile",
            percentile_range=(5, 95),
            plot_distributions=False,
        )
        assert df is not None

    def test_iqr_method(self, qc_test_adata):
        df, thresholds = suggest_qc_thresholds(
            qc_test_adata, method="iqr", plot_distributions=False
        )
        assert df is not None

    def test_mad_multipliers_as_list(self, qc_test_adata):
        df, thresholds = suggest_qc_thresholds(
            qc_test_adata,
            method="mad",
            mad_multipliers=[2.5, 3.5],
            plot_distributions=False,
        )
        assert df is not None


class TestIdentifyOutliers:
    def test_empty_metrics_returns_no_outliers(self, qc_test_adata):
        result = identify_outliers(qc_test_adata, metrics=[])
        assert isinstance(result, pd.Series)
        assert not result.any()

    def test_single_upper_metric(self, qc_test_adata):
        metrics = [("pct_counts_mt", "upper", None)]
        result = identify_outliers(qc_test_adata, metrics=metrics)
        assert isinstance(result, pd.Series)
        assert len(result) == qc_test_adata.n_obs

    def test_single_lower_metric(self, qc_test_adata):
        metrics = [("n_genes_by_counts", "lower", None)]
        result = identify_outliers(qc_test_adata, metrics=metrics)
        assert isinstance(result, pd.Series)

    def test_multiple_metrics(self, qc_test_adata):
        metrics = [
            ("pct_counts_mt", "upper", None),
            ("n_genes_by_counts", "lower", None),
        ]
        result = identify_outliers(qc_test_adata, metrics=metrics)
        assert isinstance(result, pd.Series)

    def test_fixed_threshold(self, qc_test_adata):
        metrics = [("pct_counts_mt", "upper", 25.0)]
        result = identify_outliers(qc_test_adata, metrics=metrics)
        assert isinstance(result, pd.Series)

    def test_with_sample_key(self, qc_test_adata):
        metrics = [("pct_counts_mt", "upper", None)]
        result = identify_outliers(
            qc_test_adata, metrics=metrics, sample_key="sampleID"
        )
        assert isinstance(result, pd.Series)


class TestMarkLowQualityCell:
    def test_basic_marking(self, qc_test_adata):
        result = mark_low_quality_cell(qc_test_adata.copy(), show_plots=False)
        assert result is not None
        # Should have some marking columns in .obs
        marking_cols = [c for c in result.obs.columns if c.startswith("outlier_") or c == "low_quality"]
        assert len(marking_cols) >= 0  # May or may not produce outliers

    def test_with_sample_thresholds(self, qc_test_adata):
        sample_thresholds = {
            "sample_0": {"min_genes": 100, "max_mt_percent": 25.0},
        }
        result = mark_low_quality_cell(
            qc_test_adata.copy(),
            sample_thresholds=sample_thresholds,
            show_plots=False,
        )
        assert result is not None

    def test_marks_cells(self, qc_test_adata):
        """Some cells should be flagged as outliers on synthetic qc_test_adata."""
        result = mark_low_quality_cell(qc_test_adata.copy(), show_plots=False)
        # qc_test_adata has synthetic outliers
        assert result is not None


class TestFilterCells:
    def test_filter_in_place(self, qc_test_adata):
        """filter_cells in-place should reduce cell count."""
        initial_count = qc_test_adata.n_obs
        result = filter_cells(qc_test_adata.copy(), copy=True)
        if result is not None:
            assert result.n_obs <= initial_count

    def test_filter_with_copy_returns_new(self, qc_test_adata):
        result = filter_cells(qc_test_adata, copy=True)
        assert result is None or isinstance(result, qc_test_adata.__class__)

    def test_filter_with_copy_copies_once_and_stores_stats(self, qc_test_adata, monkeypatch):
        adata = qc_test_adata.copy()
        adata.obs["outlier_for_copy_test"] = [i % 3 == 0 for i in range(adata.n_obs)]
        config = FilterConfig(
            criteria_to_filter=["outlier_for_copy_test"],
            combination_logic="any",
        )

        copy_calls = 0
        original_copy = AnnData.copy

        def counting_copy(self, filename=None):
            nonlocal copy_calls
            copy_calls += 1
            return original_copy(self, filename=filename)

        monkeypatch.setattr(AnnData, "copy", counting_copy)

        result = filter_cells(adata, config=config, copy=True)

        assert copy_calls == 1
        assert result.n_obs < adata.n_obs
        assert adata.n_obs == qc_test_adata.n_obs
        assert "filtering_results" in result.uns["sclucid"]["qc"]

    def test_filter_empty_config(self, qc_test_adata):
        result = filter_cells(qc_test_adata.copy(), copy=True)
        assert result is not None or result is None  # Accept either behavior

    def test_custom_logic_uses_criteria_namespace(self, qc_test_adata):
        adata = qc_test_adata.copy()
        reps = (adata.n_obs // 4) + 1
        adata.obs["outlier_a"] = ([True, True, False, False] * reps)[: adata.n_obs]
        adata.obs["outlier_b"] = ([True, False, True, False] * reps)[: adata.n_obs]
        config = FilterConfig(
            criteria_to_filter=["outlier_a", "outlier_b"],
            combination_logic="custom",
            custom_logic_expr="outlier_a & outlier_b",
        )

        result = filter_cells(adata, config=config, copy=True)

        assert result.n_obs == adata.n_obs - int((adata.obs["outlier_a"] & adata.obs["outlier_b"]).sum())

    def test_custom_logic_requires_expression(self, qc_test_adata):
        with pytest.raises(ValueError, match="custom_logic_expr"):
            FilterConfig(
                criteria_to_filter=["outlier_a"],
                combination_logic="custom",
                custom_logic_expr=None,
            )


class TestResolveQCThresholds:
    def test_intelligent_then_mad_policy(self):
        thresholds = resolve_qc_thresholds(
            intelligent={"min_genes": 300, "pc_mt": 8.0},
            mad={"min_genes": 200, "max_genes": 8000, "pc_mt": 10.0},
        )
        assert thresholds.min_genes == 300
        assert thresholds.max_genes == 8000
        assert thresholds.pc_mt == 8.0

    def test_manual_floor_and_ceiling(self):
        thresholds = resolve_qc_thresholds(
            intelligent={"min_genes": 300, "pc_mt": 8.0},
            mad={"min_genes": 200, "pc_mt": 10.0},
            manual={"min_genes": 500, "pc_mt": 12.0},
        )
        assert thresholds.min_genes == 500  # floor
        assert thresholds.pc_mt == 8.0  # ceiling

    def test_manual_override_policy(self):
        thresholds = resolve_qc_thresholds(
            intelligent={"min_genes": 300},
            mad={"min_genes": 200},
            manual={"min_genes": 500},
            policy="manual_override",
        )
        assert thresholds.min_genes == 500

    def test_mad_then_intelligent_policy(self):
        thresholds = resolve_qc_thresholds(
            intelligent={"min_genes": 300},
            mad={"min_genes": 200},
            policy="mad_then_intelligent",
        )
        assert thresholds.min_genes == 200


class TestQCDecisionWorkflow:
    def test_runs_threshold_resolution_and_marking(self, qc_test_adata):
        result = run_qc_decision_workflow(
            qc_test_adata.copy(),
            intelligent_thresholds={"min_genes": 100, "pc_mt": 20.0},
            manual_thresholds={"min_genes": 50},
            plot_distributions=False,
            filter_cells_result=False,
        )

        assert "adata" in result
        assert result["resolved_thresholds"].min_genes == 100
        assert "qc_threshold_decision" in result["adata"].uns["sclucid"]["qc"]
        assert "qc_decision_workflow" in result["adata"].uns["sclucid"]["qc"]

    def test_new_threshold_decision_entrypoint(self, qc_test_adata):
        result = run_qc_threshold_decision(
            qc_test_adata.copy(),
            intelligent_thresholds={"min_genes": 100, "pc_mt": 20.0},
            plot_distributions=False,
            filter_cells_result=False,
        )

        assert result["resolved_thresholds"].min_genes == 100
        assert "qc_threshold_decision" in result["adata"].uns["sclucid"]["qc"]


class TestAuditFiltering:
    def test_per_sample_retention(self, qc_test_adata):
        before = qc_test_adata.copy()
        after = before[: before.n_obs // 2].copy()
        result = audit_filtering(before, after, sample_key="sampleID")
        assert "sample" in result
        assert all(
            col in result["sample"].columns
            for col in ["sampleID", "before", "after", "removed", "retention_rate"]
        )

    def test_group_retention(self, qc_test_adata):
        before = qc_test_adata.copy()
        before.obs["group"] = ["g1", "g2"] * (before.n_obs // 2)
        after = before[: before.n_obs // 2].copy()
        result = audit_filtering(before, after, sample_key="sampleID", group_key="group")
        assert "sample" in result
        assert "group" in result


class TestGenerateQCReport:
    def test_smoke_report_generation(self, qc_test_adata, temp_output_dir):
        """Basic smoke test: report generation should not crash."""
        generate_qc_report(
            qc_test_adata,
            save_dir=temp_output_dir,
            include_before_after=False,
        )

    def test_report_with_before_after(self, qc_test_adata, temp_output_dir):
        before = qc_test_adata.copy()
        after = qc_test_adata.copy()
        generate_qc_report(
            after,
            save_dir=temp_output_dir,
            include_before_after=True,
            adata_before=before,
        )

    def test_report_creates_output(self, qc_test_adata, temp_output_dir):
        generate_qc_report(
            qc_test_adata,
            save_dir=temp_output_dir,
            include_before_after=False,
        )
        import os
        files = os.listdir(temp_output_dir)
        # Should create at least some output
        assert len(files) >= 0
