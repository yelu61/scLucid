"""
pyDWLS - Python implementation of Dampened Weighted Least Squares (DWLS)

DWLS is a deconvolution method for estimating cell type proportions from
bulk RNA-seq data using single-cell RNA-seq reference data.

Main Components
---------------
- DWLS : Main deconvolution class
- SignatureBuilder : Build cell type signature matrices
- MarkerSelector : Select marker genes for each cell type
- DampenedWLS : Solve weighted least squares with dampening
- CrossValidator : Cross-validation utilities

Example:
-------
>>> from pyDWLS import DWLS
>>>
>>> # Initialize
>>> dwls = DWLS()
>>>
>>> # Build signature matrix
>>> signature = dwls.build_signature_matrix(sc_data, cell_type_labels)
>>>
>>> # Deconvolve
>>> proportions = dwls.deconvolve(bulk_data)
>>>
>>> # Results
>>> print(proportions.head())

References:
----------
- DWLS paper: https://doi.org/10.1016/j.celrep.2019.04.045
"""

__version__ = "0.1.0"
__author__ = "scLucid"

from .core import DWLS
from .markers import MarkerSelector
from .signature import SignatureBuilder
from .solver import DampenedWLS, solve_nnls
from .utils import (
    align_data,
    create_pseudo_bulk,
    filter_genes,
    normalize_data,
)
from .validation import CrossValidator

__all__ = [
    # Main class
    "DWLS",
    # Components
    "SignatureBuilder",
    "DampenedWLS",
    "MarkerSelector",
    "CrossValidator",
    # Utilities
    "solve_nnls",
    "normalize_data",
    "filter_genes",
    "create_pseudo_bulk",
    "align_data",
]
