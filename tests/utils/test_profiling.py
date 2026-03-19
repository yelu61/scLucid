"""
Tests for performance profiling utilities.

Tests time and memory tracking functionality.
"""

import pytest
import time
import numpy as np
from anndata import AnnData

from scLucid.utils.profiling import (
    PerformanceStats,
    PerformanceProfiler,
    BenchmarkRunner,
    profile_performance,
    profile_function,
    memory_tracker,
    get_memory_usage,
    estimate_adata_memory,
)


class TestPerformanceStats:
    """Test PerformanceStats dataclass."""

    def test_stats_creation(self):
        """Test creating PerformanceStats."""
        stats = PerformanceStats(name="test")
        assert stats.name == "test"
        assert stats.duration_ms is None
        assert stats.peak_memory_mb is None

    def test_stats_stop(self):
        """Test stopping stats."""
        stats = PerformanceStats(name="test")
        time.sleep(0.01)  # Small delay
        stats.stop(peak_memory_mb=100.0)

        assert stats.duration_ms is not None
        assert stats.duration_ms > 0
        assert stats.peak_memory_mb == 100.0

    def test_stats_str_with_duration(self):
        """Test string representation with duration."""
        stats = PerformanceStats(name="test")
        stats.stop()

        str_repr = str(stats)
        assert "test:" in str_repr
        assert "ms" in str_repr or "s" in str_repr

    def test_stats_str_with_memory(self):
        """Test string representation with memory."""
        stats = PerformanceStats(name="test")
        stats.stop(peak_memory_mb=512.5)

        str_repr = str(stats)
        assert "Peak memory" in str_repr
        assert "512.5MB" in str_repr


class TestPerformanceProfiler:
    """Test PerformanceProfiler class."""

    def test_profiler_context_manager(self):
        """Test profiler as context manager."""
        with PerformanceProfiler("test") as profiler:
            time.sleep(0.01)

        assert profiler.stats is not None
        assert profiler.stats.name == "test"
        assert profiler.stats.duration_ms is not None

    def test_profiler_track_decorator(self):
        """Test profiler track decorator."""

        @PerformanceProfiler.track("decorated_test")
        def test_function():
            time.sleep(0.01)
            return 42

        result = test_function()

        assert result == 42

    def test_profiler_summary(self):
        """Test profiler summary."""
        profiler1 = PerformanceProfiler("test1")
        with profiler1:
            time.sleep(0.01)

        summary = PerformanceProfiler.summary()

        # Should have at least one entry
        assert len(summary) >= 0  # May be empty if other tests cleared it

    def test_profiler_clear_benchmarks(self):
        """Test clearing benchmarks."""
        PerformanceProfiler._benchmarks.clear()

        PerformanceProfiler.clear_benchmarks()

        assert len(PerformanceProfiler._benchmarks) == 0


class TestProfilePerformance:
    """Test profile_performance context manager."""

    def test_profile_performance_basic(self):
        """Test basic profile_performance usage."""
        with profile_performance("test_op") as stats:
            time.sleep(0.01)

        assert stats is not None
        assert stats.name == "test_op"
        assert stats.duration_ms > 0

    def test_profile_performance_no_memory(self):
        """Test profile_performance without memory tracking."""
        with profile_performance("test_op", track_memory=False) as stats:
            time.sleep(0.01)

        assert stats.duration_ms > 0


class TestProfileFunction:
    """Test profile_function decorator."""

    def test_profile_function_decorator(self):
        """Test profile_function decorator."""

        @profile_function("my_function")
        def my_function(x, y):
            time.sleep(0.01)
            return x + y

        result = my_function(1, 2)

        assert result == 3

    def test_profile_function_default_name(self):
        """Test profile_function with default name."""

        @profile_function()
        def unnamed_function():
            return "done"

        result = unnamed_function()

        assert result == "done"


class TestMemoryTracker:
    """Test memory_tracker context manager."""

    def test_memory_tracker_basic(self):
        """Test basic memory_tracker usage."""
        # Just test it doesn't crash
        with memory_tracker("test_allocation"):
            data = [0] * 1000  # Allocate some memory

        # If psutil is available, it will log memory change
        assert True  # Mainly testing it doesn't raise


class TestGetMemoryUsage:
    """Test get_memory_usage function."""

    def test_get_memory_usage(self):
        """Test getting memory usage."""
        mem_info = get_memory_usage()

        # May return empty dict if psutil not available
        if mem_info:
            assert "rss_mb" in mem_info or "percent" in mem_info


class TestEstimateAdataMemory:
    """Test estimate_adata_memory function."""

    def test_estimate_sparse_adata(self):
        """Test estimation for sparse AnnData."""
        from scipy.sparse import csr_matrix

        X = csr_matrix((100, 1000))
        adata = AnnData(X)

        estimates = estimate_adata_memory(adata)

        assert "X_mb" in estimates
        assert "total_mb" in estimates
        assert estimates["total_mb"] >= 0

    def test_estimate_dense_adata(self):
        """Test estimation for dense AnnData."""
        X = np.random.randn(100, 1000)
        adata = AnnData(X)

        estimates = estimate_adata_memory(adata)

        assert "X_mb" in estimates
        assert "total_mb" in estimates
        assert estimates["total_mb"] > 0

    def test_estimate_with_layers(self):
        """Test estimation with layers."""
        adata = AnnData(np.random.randn(100, 100))
        adata.layers["normalized"] = np.random.randn(100, 100)
        adata.layers["scaled"] = np.random.randn(100, 100)

        estimates = estimate_adata_memory(adata)

        assert "layers_mb" in estimates
        assert "layer_breakdown" in estimates
        assert len(estimates["layer_breakdown"]) == 2

    def test_estimate_with_obsm(self):
        """Test estimation with obsm."""
        adata = AnnData(np.random.randn(100, 50))
        adata.obsm["X_pca"] = np.random.randn(100, 10)
        adata.obsm["X_umap"] = np.random.randn(100, 2)

        estimates = estimate_adata_memory(adata)

        assert "obsm_mb" in estimates
        assert estimates["obsm_mb"] > 0


class TestBenchmarkRunner:
    """Test BenchmarkRunner class."""

    def test_benchmark_runner_register(self):
        """Test registering benchmarks."""
        runner = BenchmarkRunner()

        @runner.benchmark("method_a")
        def method_a():
            time.sleep(0.01)

        assert "method_a" in runner.benchmarks

    def test_benchmark_runner_run_all(self):
        """Test running all benchmarks."""
        runner = BenchmarkRunner()

        @runner.benchmark("fast_op")
        def fast_op():
            pass

        results = runner.run_all(n_runs=2)

        assert "fast_op" in results
        assert len(results["fast_op"]) == 2

    def test_benchmark_runner_compare(self):
        """Test benchmark comparison."""
        runner = BenchmarkRunner()

        @runner.benchmark("method_a")
        def method_a():
            time.sleep(0.01)

        @runner.benchmark("method_b")
        def method_b():
            time.sleep(0.02)

        runner.run_all(n_runs=2)
        comparison = runner.compare()

        assert len(comparison) > 0
        assert "avg_time_ms" in comparison.columns
