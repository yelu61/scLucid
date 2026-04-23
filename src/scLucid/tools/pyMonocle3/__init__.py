"""
pyMonocle3 - Python implementation of Monocle3 for trajectory analysis (R-free)

This package provides a complete, R-free Python implementation of Monocle3
for single-cell trajectory analysis.

Main Components
---------------
- CellDataSet : Core data structure for storing expression and metadata
- preprocess_cds : Normalization and dimensionality reduction
- reduce_dimension : UMAP/tSNE visualization
- cluster_cells : Graph-based clustering (Leiden/Louvain)
- learn_graph : Construct principal graph for trajectory
- order_cells : Calculate pseudotime
- plot_cells : Visualization of cells and trajectories

Example:
-------
>>> from pyMonocle3 import CellDataSet, preprocess_cds, reduce_dimension
>>> from pyMonocle3 import cluster_cells, learn_graph, order_cells, plot_cells
>>>
>>> # Create CellDataSet
>>> cds = CellDataSet(
...     expression_data=expr_matrix,
...     cell_metadata=cell_meta,
...     gene_metadata=gene_meta
... )
>>>
>>> # Run analysis pipeline
>>> cds = preprocess_cds(cds, num_dim=50)
>>> cds = reduce_dimension(cds, reduction_method='UMAP')
>>> cds = cluster_cells(cds)
>>> cds = learn_graph(cds)
>>> cds = order_cells(cds)
>>>
>>> # Visualize
>>> fig, ax = plot_cells(cds, color_cells_by='pseudotime')
"""

__version__ = "0.1.0"
__author__ = "scLucid"

# Core data structure
# Clustering
from .clustering import (
    cluster_cells,
    find_cluster_markers,
    group_cells,
    partition_cells,
)
from .core import (
    CellDataSet,
    create_cds_from_scanpy,
    export_to_scanpy,
    new_cell_data_set,
)

# Differential expression
from .differential import (
    aggregate_gene_expression,
    calculate_gene_modules,
    compare_genes,
    pseudotime_de,
    top_markers,
)

# Dimensionality reduction
from .dimensionality import (
    reduce_dimension,
    run_pca,
    run_umap,
)

# Preprocessing
from .preprocessing import (
    align_cds,
    detect_genes,
    estimate_size_factors,
    preprocess_cds,
)

# Trajectory inference
from .trajectory import (
    choose_graph_segments,
    graph_test,
    learn_graph,
    order_cells,
)

# Utilities
from .utils import (
    convert_to_dense,
    convert_to_sparse,
    detect_sparse_type,
    estimate_memory_usage,
    merge_datasets,
    normalize_expression,
    select_highly_variable_genes,
    subsample_cells,
    validate_cds,
)

# Visualization
from .visualization import (
    plot_cells,
    plot_genes_by_group,
    plot_pseudotime_heatmap,
    plot_trajectory,
)

__all__ = [
    # Core
    "CellDataSet",
    "new_cell_data_set",
    "create_cds_from_scanpy",
    "export_to_scanpy",
    # Preprocessing
    "detect_genes",
    "estimate_size_factors",
    "preprocess_cds",
    "align_cds",
    # Dimensionality
    "reduce_dimension",
    "run_pca",
    "run_umap",
    # Clustering
    "cluster_cells",
    "partition_cells",
    "group_cells",
    "find_cluster_markers",
    # Trajectory
    "learn_graph",
    "order_cells",
    "graph_test",
    "choose_graph_segments",
    # Differential
    "top_markers",
    "aggregate_gene_expression",
    "compare_genes",
    "pseudotime_de",
    "calculate_gene_modules",
    # Visualization
    "plot_cells",
    "plot_genes_by_group",
    "plot_pseudotime_heatmap",
    "plot_trajectory",
    # Utils
    "detect_sparse_type",
    "convert_to_dense",
    "convert_to_sparse",
    "normalize_expression",
    "select_highly_variable_genes",
    "subsample_cells",
    "merge_datasets",
    "validate_cds",
    "estimate_memory_usage",
]
