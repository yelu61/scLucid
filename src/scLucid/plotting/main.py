"""
Backward compatibility module.

This file re-exports all plotting functions for backward compatibility.
New code should import directly from scLucid.plotting or specific submodules.

DEPRECATED: Import from scLucid.plotting instead:
    from scLucid.plotting import plot_embedding
"""

# Re-export all functions from submodules for backward compatibility
from .advanced_plots import (
    plot_coexpression,
    plot_differential_abundance,
    plot_feature_correlation,
    plot_ridge,
)
from .annotation_plots import (
    export_annotation_report,
    plot_annotation_evidence_panel,
)
from .embedding_plots import (
    plot_embedding,
    plot_faceted_embedding,
)
from .feature_plots import (
    plot_dotplot,
    plot_marker_expression,
    plot_split_violin_with_stats,
    plot_stacked_violin,
)
from .marker_plots import (
    plot_faceted_feature,
    plot_marker_heatmap,
    plot_ranked_genes,
    plot_volcano,
)

__all__ = [
    "plot_embedding",
    "plot_faceted_embedding",
    "plot_dotplot",
    "plot_stacked_violin",
    "plot_split_violin_with_stats",
    "plot_marker_expression",
    "plot_faceted_feature",
    "plot_marker_heatmap",
    "plot_ranked_genes",
    "plot_volcano",
    "plot_ridge",
    "plot_feature_correlation",
    "plot_coexpression",
    "plot_differential_abundance",
    "plot_annotation_evidence_panel",
    "export_annotation_report",
]
