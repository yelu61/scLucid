"""
pyBayesPrism - Python implementation of BayesPrism for bulk deconvolution (R-free)

BayesPrism is a Bayesian approach for deconvolving bulk RNA-seq data into
cell type-specific expression profiles and proportions.

Main Components
---------------
- PrismConfig : Configuration parameters
- BayesPrismReference : Reference scRNA-seq data handler
- BayesPrism : Main deconvolution class
- BayesPrismEmbedding : NMF-based gene program learning
- GibbsSampler : Optimized Gibbs sampling for posterior inference

Example:
-------
>>> from pyBayesPrism import BayesPrismReference, BayesPrism, PrismConfig
>>>
>>> # Create reference
>>> ref = BayesPrismReference(
...     reference=sc_counts,
...     cell_type_labels=cell_types,
... )
>>>
>>> # Run deconvolution
>>> config = PrismConfig(n_iter=100, burnin=50)
>>> bp = BayesPrism(reference=ref, mixture=bulk_data, config=config)
>>> bp.cleanup_genes(remove_ribo=True, remove_mito=True)
>>> bp.run_deconvolution(n_cores=4)
>>>
>>> # Get results
>>> fractions = bp.get_fraction()
>>> expression = bp.get_expression()
"""

__version__ = "0.1.0"
__author__ = "scLucid"

# Configuration
from .config import (
    DeconvolutionConfig,
    PrismConfig,
    ReferenceConfig,
)
from .core import BayesPrism
from .embedding import BayesPrismEmbedding

# Core classes
from .reference import BayesPrismReference
from .sampling import GibbsSampler

# Utilities
from .utils import (
    batch_correct,
    cleanup_genes,
    compute_correlation,
    compute_rmse,
    find_outlier_genes,
    normalize_expression,
    subsample_cells,
    validate_inputs,
)

# Visualization
from .visualization import (
    plot_correlation,
    plot_cv,
    plot_fraction,
    plot_gene_programs,
    plot_program_usage,
    plot_stacked_bar,
    plot_validation_scatter,
)

__all__ = [
    # Configuration
    "PrismConfig",
    "ReferenceConfig",
    "DeconvolutionConfig",
    # Core classes
    "BayesPrismReference",
    "BayesPrism",
    "BayesPrismEmbedding",
    "GibbsSampler",
    # Visualization
    "plot_fraction",
    "plot_correlation",
    "plot_stacked_bar",
    "plot_gene_programs",
    "plot_program_usage",
    "plot_cv",
    "plot_validation_scatter",
    # Utilities
    "cleanup_genes",
    "find_outlier_genes",
    "compute_correlation",
    "compute_rmse",
    "normalize_expression",
    "batch_correct",
    "subsample_cells",
    "validate_inputs",
]
