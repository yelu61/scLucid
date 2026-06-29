"""
Test basic QC metrics calculation.

Converted from examples/quick_qc_test.py
"""


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
