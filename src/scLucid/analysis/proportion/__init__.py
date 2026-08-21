"""
Cell type proportion analysis submodule.

This module provides comprehensive tools for analyzing cell type proportions,
including multiple statistical methods:
- Pseudo-bulk: sample-level CLR/DESeq2-style tests for compositional proportions
- Legacy raw-proportion tests: retained for exploratory summaries only
- scCODA: Bayesian compositional data analysis

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

from .config import MethodSelectionConfig, ProportionConfig
from .pseudobulk import celltype_proportion_analysis
from .stats import composition_transform, compute_celltype_proportion, run_statistical_test

# Alias for backward compatibility
pb_analysis = celltype_proportion_analysis
from .methods import ProportionMethod, compare_methods, recommend_method
from .plots import (
    plot_batch_effect,
    plot_box_summary,
    plot_cell_counts,
    plot_celltype_alluvial,
    plot_celltype_correlation,
    plot_celltype_variability,
    plot_composition,
    plot_composition_pca,
    plot_composition_shift_bubble,
    plot_composition_shift_effect,
    plot_composition_transform_heatmap,
    plot_diff_stats,
    plot_effect_size_volcano,
    plot_grouped_celltype_counts,
    plot_grouped_proportion_bar,
    plot_individual_boxplots,
    plot_paired_proportion_shifts,
    plot_proportion_bar,
    plot_proportion_heatmap,
    plot_proportion_shifts,
    plot_proportion_timeseries,
    plot_proportion_with_ci,
    summarize_composition_shift,
    transform_composition,
)
from .stats import export_analysis_data
from .workflow import (
    analyze_all_methods,
    analyze_celltype_proportion,
)

# Optional scCODA import (may fail if sccoda not installed)
try:
    from .sccoda import recommend_sccoda_reference, run_sccoda, summarize_sccoda

    _sccoda_available = True
except ImportError:
    _sccoda_available = False
    recommend_sccoda_reference = None
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
    "composition_transform",
    "run_statistical_test",
    # scCODA (optional)
    "recommend_sccoda_reference",
    "run_sccoda",
    "summarize_sccoda",
    # Utility
    "export_analysis_data",
    # Plotting
    "plot_cell_counts",
    "plot_proportion_bar",
    "plot_grouped_celltype_counts",
    "plot_grouped_proportion_bar",
    "plot_celltype_alluvial",
    "plot_box_summary",
    "plot_proportion_heatmap",
    "plot_celltype_correlation",
    "plot_effect_size_volcano",
    "plot_proportion_timeseries",
    "plot_batch_effect",
    "plot_composition",
    "plot_diff_stats",
    "plot_individual_boxplots",
    "plot_proportion_shifts",
    "plot_paired_proportion_shifts",
    "plot_proportion_with_ci",
    "plot_celltype_variability",
    "plot_composition_shift_bubble",
    "plot_composition_shift_effect",
    "plot_composition_transform_heatmap",
    "plot_composition_pca",
    "summarize_composition_shift",
    "transform_composition",
]
