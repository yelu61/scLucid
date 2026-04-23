"""
CellChat: Python toolkit for inference and analysis of cell-cell communication
from single-cell and spatially resolved transcriptomics

This is a pure Python implementation of CellChat, with no R dependencies.
"""

__version__ = "2.1.0"

# Analysis functions
from .analysis import (
    compute_centrality,
    compute_network_similarity,
    identify_roles,
)

# Comparison functions
from .comparison import (
    compare_cellchat_objects,
    identify_conserved_pathways,
    identify_differential_pathways,
)
from .core import CellChat, CellChatConfig
from .database import CellChatDB, get_default_database
from .utils import create_cellchat_from_scanpy, merge_cellchat_objects

# Visualization functions
from .visualization import (
    plot_bubble,
    plot_chord_diagram,
    plot_circle_network,
    plot_contribution,
    plot_heatmap,
    plot_signaling_gene_expression,
)

__all__ = [
    # Core classes
    "CellChat",
    "CellChatConfig",
    "CellChatDB",
    # Factory functions
    "get_default_database",
    "create_cellchat_from_scanpy",
    "merge_cellchat_objects",
    # Visualization
    "plot_circle_network",
    "plot_chord_diagram",
    "plot_heatmap",
    "plot_bubble",
    "plot_contribution",
    "plot_signaling_gene_expression",
    # Analysis
    "compute_centrality",
    "identify_roles",
    "compute_network_similarity",
    # Comparison
    "compare_cellchat_objects",
    "identify_differential_pathways",
    "identify_conserved_pathways",
]
