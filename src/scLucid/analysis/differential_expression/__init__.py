"""
Differential Expression and Enrichment Analysis Module.

This module provides comprehensive differential expression analysis and functional
enrichment capabilities for single-cell RNA-seq data.

The module is organized into:
- Core DE functions (find_markers, filter_markers, etc.)
- Enrichment analysis (run_enrichment, batch analysis)
- Visualization (volcano plots, heatmaps)
- High-level workflows (characterize_clusters)
- Result management (ResultManager, save/load)

Usage:
------
>>> from scLucid.analysis import find_markers, run_enrichment
>>> from scLucid.analysis.differential_expression import characterize_clusters
>>>
>>> # Find markers
>>> config = DifferentialConfig(groupby="leiden", method="wilcoxon")
>>> markers = find_markers(adata, config)
>>>
>>> # Run enrichment
>>> enrichment = run_enrichment(markers, EnrichmentConfig())
>>>
>>> # Or use the complete workflow
>>> characterize_clusters(adata, groupby="leiden")
"""

# Core DE functions
from .de_core import (
    find_markers,
    filter_markers,
    compare_groups,
    compare_conditions,
    get_conserved_markers,
)

# Enrichment functions
from .enrichment import (
    run_enrichment,
    export_enrichment_results,
    batch_celltype_deg_enrichment,
)

# High-level workflows
from .de_workflows import (
    characterize_clusters,
    summarize_markers_and_enrichment,
)

# Visualization
from .de_plots import (
    visualize_markers,
    plot_volcano,
    plot_multi_cluster_deg,
)

# Result management
from .de_utils import (
    ResultManager,
    save_results,
    load_results,
)

__all__ = [
    # Core DE
    "find_markers",
    "filter_markers",
    "compare_groups",
    "compare_conditions",
    "get_conserved_markers",
    # Enrichment
    "run_enrichment",
    "export_enrichment_results",
    "batch_celltype_deg_enrichment",
    # Workflows
    "characterize_clusters",
    "summarize_markers_and_enrichment",
    # Visualization
    "visualize_markers",
    "plot_volcano",
    "plot_multi_cluster_deg",
    # Result management
    "ResultManager",
    "save_results",
    "load_results",
]
