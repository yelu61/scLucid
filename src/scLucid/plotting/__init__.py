"""
General-Purpose Visualization Library for scLucid.

This module provides a comprehensive suite of plotting functions to generate
publication-quality visualizations for all stages of single-cell analysis,
from QC and preprocessing to final annotation and downstream analysis.

The module is organized into:
- Embedding plots: UMAP/t-SNE visualizations
- Feature plots: Gene expression, dot plots, violin plots
- Marker plots: Heatmaps, ranked genes, volcano plots
- Advanced plots: Correlation, ridge, co-expression plots
"""

# Embedding plots
# Advanced plots
from .advanced_plots import (
    plot_coexpression,
    plot_differential_abundance,
    plot_feature_correlation,
    plot_ridge,
)

# Annotation evidence plots
from .annotation_plots import (
    export_annotation_report,
    plot_annotation_evidence_panel,
)
from .embedding_plots import (
    plot_embedding,
    plot_faceted_embedding,
)

# Feature/gene expression plots
from .feature_plots import (
    plot_dotplot,
    plot_marker_expression,
    plot_split_violin_with_stats,
    plot_stacked_violin,
)

# Marker plots
from .marker_plots import (
    plot_faceted_feature,
    plot_marker_heatmap,
    plot_ranked_genes,
    plot_volcano,
)

# Define what should be accessible when a user does `from scLucid.plotting import *`
__all__ = [
    "plot_embedding",
    "plot_faceted_embedding",
    "plot_dotplot",
    "plot_stacked_violin",
    "plot_split_violin_with_stats",
    "plot_marker_expression",
    "plot_faceted_feature",
    "plot_marker_heatmap",
    "plot_volcano",
    "plot_ridge",
    "plot_feature_correlation",
    "plot_coexpression",
    "plot_differential_abundance",
    "plot_annotation_evidence_panel",
    "export_annotation_report",
]
