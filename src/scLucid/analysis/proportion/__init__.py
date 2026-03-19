"""
Cell type proportion analysis submodule.

This module provides comprehensive tools for analyzing cell type proportions,
including multiple statistical methods:
- Pseudo-bulk: Traditional statistical tests (DESeq2, t-test, Wilcoxon)
- scCODA: Bayesian compositional data analysis
- Milo: Neighborhood-based cell-level analysis (TODO)

Usage:
------
>>> from scLucid.analysis import analyze_celltype_proportion, recommend_method
>>>
>>> # Auto-recommend and analyze
>>> method = recommend_method(adata, sample_col="sample", condition_col="condition")
>>> result = analyze_celltype_proportion(adata, method=method)
>>>
>>> # Or specify method directly
>>> prop_df, stat_df = analyze_celltype_proportion(adata, method='pseudobulk')
"""

from .config import ProportionConfig, MethodSelectionConfig
from .pseudobulk import (
    celltype_proportion_analysis,
    compute_celltype_proportion,
    run_statistical_test,
)
# Alias for backward compatibility
pb_analysis = celltype_proportion_analysis
from .stats import export_analysis_data
from .plots import (
    plot_cell_counts,
    plot_proportion_bar,
    plot_box_summary,
    plot_proportion_heatmap,
    plot_celltype_correlation,
    plot_effect_size_volcano,
    plot_proportion_timeseries,
    plot_batch_effect,
)
from .methods import ProportionMethod, recommend_method, compare_methods
from .workflow import (
    analyze_celltype_proportion,
    analyze_all_methods,
)

# Optional scCODA import (may fail if sccoda not installed)
try:
    from .sccoda import run_sccoda, summarize_sccoda
    _sccoda_available = True
except ImportError:
    _sccoda_available = False
    run_sccoda = None
    summarize_sccoda = None

__all__ = [
    # Configuration
    "ProportionConfig",
    "MethodSelectionConfig",
    # Main workflow
    "analyze_celltype_proportion",
    "analyze_all_methods",
    "celltype_proportion_analysis",  # Pseudo-bulk analysis
    "pb_analysis",  # Alias for celltype_proportion_analysis
    # Method selection
    "ProportionMethod",
    "recommend_method",
    "compare_methods",
    # Pseudo-bulk
    "compute_celltype_proportion",
    "run_statistical_test",
    # scCODA (optional)
    "run_sccoda",
    "summarize_sccoda",
    # Utility
    "export_analysis_data",
    # Plotting
    "plot_cell_counts",
    "plot_proportion_bar",
    "plot_box_summary",
    "plot_proportion_heatmap",
    "plot_celltype_correlation",
    "plot_effect_size_volcano",
    "plot_proportion_timeseries",
    "plot_batch_effect",
]
