"""Tests for QC recommendation benchmark metrics and report templates."""

import json

import pytest

from scLucid import qc
from scLucid.qc.benchmark import (
    QC_BENCHMARK_SCHEMA_VERSION,
    compute_marker_fidelity,
    compute_retention_metrics,
    evaluate_qc_benchmark,
    export_qc_benchmark_report,
    infer_qc_benchmark_profile,
)
from scLucid.qc.config import (
    DoubletConfig,
    MarkingConfig,
    MetricsReportingConfig,
    QCThresholds,
    QCWorkflowConfig,
)
from tests.fixtures.data_loader import load_test_data
from tests.fixtures.synthetic_data import SyntheticDataGenerator


def _make_marker_adata(n_cells=240, n_genes=800):
    generator = SyntheticDataGenerator(random_state=61)
    return generator.generate_adata(
        n_cells=n_cells,
        n_genes=n_genes,
        n_cell_types=4,
        n_batches=3,
        with_qc_metrics=True,
        with_cell_types=True,
        with_batches=True,
        sparsity=0.8,
    )


def _make_noninteractive_config(save_dir=None):
    return QCWorkflowConfig(
        sample_key="sampleID",
        tissue_type="pbmc",
        save_dir=str(save_dir) if save_dir is not None else None,
        use_parallel=False,
        use_recommendations=False,
        threshold_mode="pooled",
        metrics_reporting_config=MetricsReportingConfig(
            show_plots=False,
            plot_top_genes=False,
            plot_violin=False,
            plot_scatter=False,
            export_stats=False,
            print_stats=False,
        ),
        marking_config=MarkingConfig(
            show_plots=False,
            plot_outliers=False,
            thresholds=QCThresholds(min_genes=0, min_counts=0, pc_mt=100.0),
        ),
        doublet_config=DoubletConfig(
            run_algorithm=False,
            use_heuristics=False,
            show_plots=False,
            plot_summary=False,
            plot_bar=False,
            plot_scatter=False,
            plot_upset=False,
            export_stats=False,
        ),
        filter_config={
            "criteria_to_filter": ["outlier_min_genes", "outlier_mt"],
            "combination_logic": "any",
        },
    )


def test_infer_qc_benchmark_profile():
    assert infer_qc_benchmark_profile(tissue_type="pbmc") == "pbmc"
    assert infer_qc_benchmark_profile(tissue_type="lung_tumor") == "tumor"
    assert infer_qc_benchmark_profile(dataset_type="cell_line") == "cell_line"
    assert infer_qc_benchmark_profile(tissue_type="kidney") == "tissue"


def test_retention_metrics_are_stratified():
    adata = _make_marker_adata()
    retained = adata[:120].copy()
    metrics = compute_retention_metrics(
        adata,
        retained,
        sample_key="sampleID",
        cell_type_key="cell_type",
    )

    assert metrics["initial_cells"] == 240
    assert metrics["final_cells"] == 120
    assert metrics["retention_rate"] == 0.5
    assert metrics["per_sample"]
    assert metrics["per_cell_type"]


def test_marker_fidelity_uses_synthetic_marker_sets():
    adata = _make_marker_adata()
    retained = adata[:180].copy()
    fidelity = compute_marker_fidelity(adata, retained)

    assert fidelity["available"] is True
    assert fidelity["n_marker_sets_available"] > 0
    assert 0 <= fidelity["overall_marker_fidelity"] <= 1


def test_evaluate_qc_benchmark_has_profile_checks():
    adata = _make_marker_adata()
    retained = adata[:180].copy()
    benchmark = evaluate_qc_benchmark(
        adata,
        retained,
        tissue_type="pbmc",
        sample_key="sampleID",
        cell_type_key="cell_type",
    )

    assert benchmark["schema_version"] == QC_BENCHMARK_SCHEMA_VERSION
    assert benchmark["profile"] == "pbmc"
    assert benchmark["retention"]["retention_rate"] == 0.75
    assert {check["name"] for check in benchmark["checks"]} >= {
        "minimum_retention",
        "maximum_retention",
        "marker_fidelity",
    }
    assert benchmark["assessment"]["risk_level"] in {"low", "moderate", "high", "critical"}
    assert benchmark["assessment"]["summary"]


def test_evaluate_qc_benchmark_flags_stratified_retention_bias():
    adata = _make_marker_adata(n_cells=300)
    retained = adata[adata.obs["sampleID"] != "batch_0"].copy()
    benchmark = evaluate_qc_benchmark(
        adata,
        retained,
        tissue_type="pbmc",
        sample_key="sampleID",
        cell_type_key="cell_type",
    )

    check_names = {check["name"] for check in benchmark["checks"]}
    assert "minimum_sample_retention" in check_names
    assert "sample_retention_spread" in check_names
    assert benchmark["status"] in {"review_required", "fail"}
    assert benchmark["assessment"]["review_required"] is True
    assert benchmark["assessment"]["recommendations"]


def test_evaluate_qc_benchmark_fails_when_retention_collapses():
    adata = _make_marker_adata()
    benchmark = evaluate_qc_benchmark(
        adata,
        adata[:0].copy(),
        tissue_type="pbmc",
        sample_key="sampleID",
    )

    assert benchmark["status"] == "fail"
    assert benchmark["assessment"]["risk_level"] == "critical"
    assert any(
        item["priority"] == "blocking"
        for item in benchmark["assessment"]["recommendations"]
    )


def test_pbmc_real_data_benchmark_smoke():
    try:
        adata = load_test_data("pbmc3k", subsample=300)
    except FileNotFoundError as exc:
        pytest.skip(str(exc))

    retained = adata[:240].copy()
    benchmark = evaluate_qc_benchmark(
        adata,
        retained,
        tissue_type="pbmc",
        sample_key="sampleID",
    )

    assert benchmark["profile"] == "pbmc"
    assert benchmark["retention"]["initial_cells"] == 300
    assert benchmark["retention"]["final_cells"] == 240
    assert benchmark["checks"]


def test_export_qc_benchmark_report(tmp_path):
    adata = _make_marker_adata()
    benchmark = evaluate_qc_benchmark(adata, adata[:180].copy(), tissue_type="pbmc")

    paths = export_qc_benchmark_report(benchmark, tmp_path)
    assert json.loads((tmp_path / "qc_benchmark.json").read_text())["profile"] == "pbmc"
    assert (tmp_path / "qc_benchmark.md").exists()
    assert paths["json"].endswith("qc_benchmark.json")


def test_run_standard_qc_stores_benchmark_summary(tmp_path):
    adata = _make_marker_adata()
    result = qc.run_standard_qc(
        adata,
        config=_make_noninteractive_config(save_dir=tmp_path),
        show_progress=False,
    )
    review = result.uns["sclucid"]["qc"]["review_summary"]["data"]

    assert "benchmark_summary" in review
    assert review["benchmark_summary"]["schema_version"] == QC_BENCHMARK_SCHEMA_VERSION
    assert (tmp_path / "qc_benchmark.json").exists()
    assert (tmp_path / "qc_benchmark.md").exists()
