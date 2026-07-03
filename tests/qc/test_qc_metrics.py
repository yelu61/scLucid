"""
Test basic QC metrics calculation.

Converted from examples/quick_qc_test.py
"""


import numpy as np
import pytest

from scLucid.qc import calculate_qc_metric
from tests.fixtures.synthetic_data import generate_minimal_adata


def test_basic_qc_metrics():
    """Test basic QC metrics calculation."""
    adata = generate_minimal_adata(n_cells=1000, n_genes=1000)

    # Calculate QC metrics
    calculate_qc_metric(
        adata,
        sample_key="sampleID",
        show_plots=False,
        plot_top_genes=False,
        plot_violin=False,
        plot_scatter=False,
        export_stats=False,
        print_stats=False,
    )

    # Check required metrics exist
    required_metrics = [
        "n_genes_by_counts",
        "total_counts",
        "log1p_total_counts",
        "pct_counts_mt",
    ]

    for metric in required_metrics:
        assert metric in adata.obs, f"Metric {metric} not found in adata.obs"

    # Check values are reasonable
    assert adata.obs["n_genes_by_counts"].min() >= 0
    assert adata.obs["total_counts"].min() >= 0
    assert adata.obs["pct_counts_mt"].min() >= 0
    assert adata.obs["pct_counts_mt"].max() <= 100
    for n in [20, 50, 100]:
        assert f"pct_counts_in_top_{n}_genes" in adata.obs


def test_qc_metrics_default_percent_top_is_clipped_for_small_panels():
    """Default top-gene fractions should not request more genes than available."""
    adata = generate_minimal_adata(n_cells=100, n_genes=30)
    adata.X = adata.X + 1

    calculate_qc_metric(
        adata,
        sample_key="sampleID",
        show_plots=False,
        plot_top_genes=False,
        plot_violin=False,
        plot_scatter=False,
        export_stats=False,
        print_stats=False,
    )

    assert "pct_counts_in_top_20_genes" in adata.obs
    assert "pct_counts_in_top_50_genes" not in adata.obs
    assert "pct_counts_in_top_100_genes" not in adata.obs


def test_qc_metrics_rejects_unknown_reporting_kwargs():
    """Misspelled reporting options should fail instead of being silently ignored."""
    adata = generate_minimal_adata(n_cells=50, n_genes=100)

    with pytest.raises(TypeError, match="Unknown calculate_qc_metric reporting option"):
        calculate_qc_metric(
            adata,
            sample_key="sampleID",
            show_plots=False,
            plot_top_genes=False,
            plot_violin=False,
            plot_scatter=False,
            export_stats=False,
            print_stats=False,
            plot_top_gene=False,
        )


def test_qc_metrics_uses_existing_var_gene_set_columns():
    """Boolean gene-set columns already present in adata.var should become qc_vars."""
    adata = generate_minimal_adata(n_cells=50, n_genes=100)
    adata.var["cell_cycle"] = False
    adata.var.iloc[:10, adata.var.columns.get_loc("cell_cycle")] = True

    calculate_qc_metric(
        adata,
        sample_key="sampleID",
        show_plots=False,
        plot_top_genes=False,
        plot_violin=False,
        plot_scatter=False,
        export_stats=False,
        print_stats=False,
    )

    assert "pct_counts_cell_cycle" in adata.obs
    assert "cell_cycle" in adata.uns["sclucid"]["qc"]["metrics"]["params"]["qc_vars"]
    provenance = adata.uns["sclucid"]["qc"]["metrics"]["params"]["gene_set_provenance"]
    assert provenance["cell_cycle"]["source"] == "existing_var_column"
    assert provenance["cell_cycle"]["matched_genes"] == 10


def test_qc_metrics_records_count_like_and_plot_metadata():
    """Metrics metadata should make input and plotting provenance explicit."""
    adata = generate_minimal_adata(n_cells=50, n_genes=100)
    matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    adata.X = np.asarray(matrix, dtype=float)
    adata.X = adata.X + 0.25

    calculate_qc_metric(
        adata,
        sample_key="sampleID",
        show_plots=False,
        plot_top_genes=False,
        plot_violin=False,
        plot_scatter=False,
        export_stats=False,
        print_stats=False,
        random_state=123,
    )

    params = adata.uns["sclucid"]["qc"]["metrics"]["params"]
    assert params["matrix_source"] == "X"
    assert params["count_like_input"]["count_like"] is False
    assert "matrix_contains_fractional_positive_values" in params["count_like_input"]["warnings"]
    assert params["plot_sampling_policy"]["random_state"] == 123
    assert "reporting_config" in params


def test_qc_metrics_records_review_only_outlier_and_sample_context_metadata():
    """Metrics warnings should be review-only and sample context should be explicit."""
    adata = generate_minimal_adata(n_cells=80, n_genes=120)
    adata.obs["sampleID"] = ["tumor_a"] * 40 + ["normal_a"] * 40
    adata.obs["sample_role"] = ["tumor"] * 40 + ["normal"] * 40

    calculate_qc_metric(
        adata,
        sample_key="sampleID",
        sample_context_key="sample_role",
        tissue_type="lung_cancer",
        show_plots=False,
        plot_top_genes=False,
        plot_violin=False,
        plot_scatter=False,
        export_stats=False,
        print_stats=False,
    )

    metrics = adata.uns["sclucid"]["qc"]["metrics"]
    params = metrics["params"]
    assert params["sample_context"]["available"] is True
    assert params["sample_context"]["sample_context_key"] == "sample_role"
    assert params["sample_context"]["sample_context_by_sample"]["tumor_a"] == ["tumor"]
    assert params["review_only_warning_thresholds"]["mt_percent"] == 25.0
    assert "Formal filtering labels" in params["review_only_warning_note"]
    assert "review_only_outlier_tips" in metrics
    assert "outlier_qc_metrics" not in adata.obs
    assert "artifact_contract" in adata.uns["sclucid"]["qc"]


def test_qc_metrics_with_samples():
    """Test QC metrics with multiple samples."""
    adata = generate_minimal_adata(n_cells=1000, n_genes=1000)

    # Simulate multiple samples
    import numpy as np

    adata.obs["sampleID"] = np.random.choice(["sample1", "sample2"], adata.n_obs)

    calculate_qc_metric(
        adata,
        sample_key="sampleID",
        show_plots=False,
        plot_top_genes=False,
        plot_violin=False,
        plot_scatter=False,
        export_stats=False,
        print_stats=False,
    )

    # Check sample-specific metrics exist
    assert "sampleID" in adata.obs.columns
