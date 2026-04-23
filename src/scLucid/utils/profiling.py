"""
Performance profiling utilities for scLucid.

This module provides tools for monitoring execution time and memory usage:
- Performance profiling decorators
- Context managers for code blocks
- Memory usage tracking
- Benchmark utilities
"""

import functools
import gc
import logging
import time
import tracemalloc
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

import numpy as np
import pandas as pd
from anndata import AnnData

log = logging.getLogger(__name__)


@dataclass
class PerformanceStats:
    """Statistics from a performance profiling run."""

    name: str
    start_time: float = field(default_factory=time.time)
    end_time: Optional[float] = None
    peak_memory_mb: Optional[float] = None
    current_memory_mb: Optional[float] = None
    duration_ms: Optional[float] = None
    annotations: Dict[str, Any] = field(default_factory=dict)

    def stop(self, peak_memory_mb: Optional[float] = None):
        """Stop the timer and record statistics."""
        self.end_time = time.time()
        self.duration_ms = (self.end_time - self.start_time) * 1000
        if peak_memory_mb is not None:
            self.peak_memory_mb = peak_memory_mb

    def __str__(self) -> str:
        parts = [f"{self.name}:"]
        if self.duration_ms is not None:
            if self.duration_ms < 1000:
                parts.append(f"{self.duration_ms:.1f}ms")
            else:
                parts.append(f"{self.duration_ms / 1000:.2f}s")
        if self.peak_memory_mb is not None:
            parts.append(f"Peak memory: {self.peak_memory_mb:.1f}MB")
        return " ".join(parts)


class PerformanceProfiler:
    """
    Performance profiler for tracking execution time and memory usage.

    Can be used as a context manager or decorator.

    Examples:
        >>> # As context manager
        >>> with PerformanceProfiler("HVG selection") as profiler:
        ...     find_hvgs(adata)
        ... print(profiler.stats)

        >>> # As decorator
        >>> @PerformanceProfiler.track("Clustering")
        ... def cluster_cells(adata):
        ...     return sc.tl.leiden(adata)

        >>> # Multiple benchmarks
        >>> profiler = PerformanceProfiler()
        >>> with profiler.benchmark("Step 1"):
        ...     step1()
        >>> with profiler.benchmark("Step 2"):
        ...     step2()
        >>> print(profiler.summary())
    """

    _benchmarks: Dict[str, List[PerformanceStats]] = {}

    def __init__(self, name: Optional[str] = None, track_memory: bool = True):
        """
        Initialize profiler.

        Args:
            name: Name of the profiling run
            track_memory: Whether to track memory usage (requires tracemalloc)
        """
        self.name = name
        self.track_memory = track_memory
        self.stats: Optional[PerformanceStats] = None
        self._tracemalloc_running = False

    def __enter__(self) -> "PerformanceProfiler":
        """Start profiling."""
        self.stats = PerformanceStats(name=self.name or "profile")

        if self.track_memory:
            try:
                gc.collect()  # Clean up before measuring
                tracemalloc.start()
                self._tracemalloc_running = True
            except Exception as e:
                log.debug(f"Could not start memory tracking: {e}")

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Stop profiling and log results."""
        peak_memory = None

        if self._tracemalloc_running:
            try:
                current, peak = tracemalloc.get_traced_memory()
                peak_memory = peak / 1024 / 1024  # Convert to MB
                tracemalloc.stop()
            except Exception as e:
                log.debug(f"Could not get memory stats: {e}")

        self.stats.stop(peak_memory_mb=peak_memory)

        # Log results
        if exc_type is None:
            log.info(f"Performance: {self.stats}")
        else:
            log.warning(f"Performance (failed): {self.stats}")

        # Store benchmark
        if self.name:
            if self.name not in self._benchmarks:
                self._benchmarks[self.name] = []
            self._benchmarks[self.name].append(self.stats)

        return False  # Don't suppress exceptions

    @classmethod
    def track(cls, name: Optional[str] = None, track_memory: bool = True) -> Callable:
        """
        Decorator to track function performance.

        Args:
            name: Name for the profiling run (defaults to function name)
            track_memory: Whether to track memory usage

        Example:
            >>> @PerformanceProfiler.track("Find HVGs")
            ... def find_hvgs(adata):
            ...     return sc.pp.highly_variable_genes(adata)
        """

        def decorator(func: Callable) -> Callable:
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                profiler_name = name or func.__name__
                with cls(profiler_name, track_memory=track_memory):
                    return func(*args, **kwargs)

            return wrapper

        return decorator

    @contextmanager
    def benchmark(self, name: str):
        """Context manager for benchmarking a code block."""
        with self.__class__(name, track_memory=self.track_memory) as p:
            yield p

    @classmethod
    def summary(cls) -> pd.DataFrame:
        """
        Get summary of all benchmarks.

        Returns:
            DataFrame with benchmark statistics
        """
        if not cls._benchmarks:
            return pd.DataFrame()

        rows = []
        for name, stats_list in cls._benchmarks.items():
            durations = [s.duration_ms for s in stats_list if s.duration_ms is not None]
            memories = [s.peak_memory_mb for s in stats_list if s.peak_memory_mb is not None]

            row = {
                "name": name,
                "n_runs": len(stats_list),
                "total_time_ms": sum(durations) if durations else None,
                "avg_time_ms": np.mean(durations) if durations else None,
                "min_time_ms": min(durations) if durations else None,
                "max_time_ms": max(durations) if durations else None,
                "avg_memory_mb": np.mean(memories) if memories else None,
                "peak_memory_mb": max(memories) if memories else None,
            }
            rows.append(row)

        return pd.DataFrame(rows)

    @classmethod
    def clear_benchmarks(cls):
        """Clear all stored benchmarks."""
        cls._benchmarks.clear()


# Convenience functions


@contextmanager
def profile_performance(name: str, track_memory: bool = True):
    """
    Context manager for quick performance profiling.

    Args:
        name: Name for this profiling run
        track_memory: Whether to track memory usage

    Yields:
        PerformanceStats object that gets populated on exit

    Examples:
        >>> with profile_performance("HVG selection") as stats:
        ...     find_hvgs(adata)
        ... print(stats)
        HVG selection: 1250.5ms Peak memory: 512.3MB

        >>> # Without memory tracking
        >>> with profile_performance("Quick step", track_memory=False) as stats:
        ...     do_something()
    """
    profiler = PerformanceProfiler(name, track_memory=track_memory)
    with profiler:
        yield profiler.stats


def profile_function(name: Optional[str] = None, track_memory: bool = True):
    """
    Decorator for profiling function performance.

    Args:
        name: Name for the profiling run (defaults to function name)
        track_memory: Whether to track memory usage

    Example:
        >>> @profile_function("Clustering")
        ... def cluster_cells(adata, resolution=1.0):
        ...     sc.tl.leiden(adata, resolution=resolution)
        ...     return adata
    """
    return PerformanceProfiler.track(name=name, track_memory=track_memory)


# ============================================================================
# Memory profiling utilities
# ============================================================================


def get_memory_usage() -> Dict[str, float]:
    """
    Get current memory usage statistics.

    Returns:
        Dictionary with memory usage in MB
    """
    try:
        import psutil

        process = psutil.Process()
        mem_info = process.memory_info()
        return {
            "rss_mb": mem_info.rss / 1024 / 1024,
            "vms_mb": mem_info.vms / 1024 / 1024,
            "percent": process.memory_percent(),
        }
    except ImportError:
        log.debug("psutil not available, cannot get memory usage")
        return {}


def estimate_adata_memory(adata: AnnData) -> Dict[str, float]:
    """
    Estimate memory usage of an AnnData object.

    Args:
        adata: AnnData object

    Returns:
        Dictionary with memory estimates in MB
    """
    estimates = {}

    # Main matrix
    if hasattr(adata.X, "data"):  # sparse
        estimates["X_mb"] = adata.X.data.nbytes / 1024 / 1024
    else:  # dense
        estimates["X_mb"] = adata.X.nbytes / 1024 / 1024

    # Layers
    layer_memory = {}
    for name, layer in adata.layers.items():
        if hasattr(layer, "data"):
            layer_memory[name] = layer.data.nbytes / 1024 / 1024
        else:
            layer_memory[name] = layer.nbytes / 1024 / 1024
    estimates["layers_mb"] = sum(layer_memory.values())
    estimates["layer_breakdown"] = layer_memory

    # obs and var
    estimates["obs_mb"] = adata.obs.memory_usage(deep=True).sum() / 1024 / 1024
    estimates["var_mb"] = adata.var.memory_usage(deep=True).sum() / 1024 / 1024

    # obsm
    obsm_memory = {}
    for name, matrix in adata.obsm.items():
        obsm_memory[name] = matrix.nbytes / 1024 / 1024
    estimates["obsm_mb"] = sum(obsm_memory.values())

    # Total estimate
    estimates["total_mb"] = (
        estimates["X_mb"]
        + estimates["layers_mb"]
        + estimates["obs_mb"]
        + estimates["var_mb"]
        + estimates["obsm_mb"]
    )

    return estimates


@contextmanager
def memory_tracker(name: str = "memory"):
    """
    Context manager to track memory changes.

    Args:
        name: Name for this tracking run

    Example:
        >>> with memory_tracker("Loading data"):
        ...     adata = sc.read_h5ad("large_file.h5ad")
        Loading data: +512.3MB
    """
    try:
        import psutil

        process = psutil.Process()
        start_mem = process.memory_info().rss / 1024 / 1024

        yield

        end_mem = process.memory_info().rss / 1024 / 1024
        delta = end_mem - start_mem
        sign = "+" if delta >= 0 else ""
        log.info(f"{name}: {sign}{delta:.1f}MB (current: {end_mem:.1f}MB)")

    except ImportError:
        log.debug("psutil not available, skipping memory tracking")
        yield


# ============================================================================
# Benchmark utilities
# ============================================================================


class BenchmarkRunner:
    """
    Run multiple benchmarks and compare results.

    Example:
        >>> runner = BenchmarkRunner()
        >>>
        >>> @runner.benchmark("Method A")
        ... def method_a():
        ...     return process_data(adata, method="a")
        >>>
        >>> @runner.benchmark("Method B")
        ... def method_b():
        ...     return process_data(adata, method="b")
        >>>
        >>> runner.run_all()
        >>> print(runner.compare())
    """

    def __init__(self):
        self.benchmarks: Dict[str, Callable] = {}
        self.results: Dict[str, PerformanceStats] = {}

    def benchmark(self, name: str) -> Callable:
        """Decorator to register a benchmark."""

        def decorator(func: Callable) -> Callable:
            self.benchmarks[name] = func
            return func

        return decorator

    def run_all(self, n_runs: int = 3) -> Dict[str, List[PerformanceStats]]:
        """
        Run all registered benchmarks.

        Args:
            n_runs: Number of times to run each benchmark

        Returns:
            Dictionary mapping benchmark names to lists of stats
        """
        all_results = {}

        for name, func in self.benchmarks.items():
            log.info(f"Running benchmark: {name}")
            results = []

            for i in range(n_runs):
                with PerformanceProfiler(f"{name}_run{i}") as p:
                    func()
                results.append(p.stats)

            all_results[name] = results

        self.results = all_results
        return all_results

    def compare(self) -> pd.DataFrame:
        """Compare benchmark results."""
        if not self.results:
            return pd.DataFrame()

        rows = []
        for name, stats_list in self.results.items():
            durations = [s.duration_ms for s in stats_list]
            memories = [s.peak_memory_mb for s in stats_list if s.peak_memory_mb]

            rows.append(
                {
                    "benchmark": name,
                    "avg_time_ms": np.mean(durations),
                    "std_time_ms": np.std(durations),
                    "avg_memory_mb": np.mean(memories) if memories else None,
                    "speedup": None,  # Will be calculated
                }
            )

        df = pd.DataFrame(rows)

        # Calculate speedup relative to slowest
        if len(df) > 1:
            max_time = df["avg_time_ms"].max()
            df["speedup"] = max_time / df["avg_time_ms"]

        return df


# ============================================================================
# Module exports
# ============================================================================

__all__ = [
    # Core classes
    "PerformanceStats",
    "PerformanceProfiler",
    "BenchmarkRunner",
    # Convenience functions
    "profile_performance",
    "profile_function",
    "memory_tracker",
    # Memory utilities
    "get_memory_usage",
    "estimate_adata_memory",
]
