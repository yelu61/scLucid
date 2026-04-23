"""
Test basic QC metrics calculation.

Converted from examples/quick_qc_test.py
"""


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
