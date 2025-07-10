"""
Analysis workflows for single-cell RNA-seq data.

This module provides standard and customizable workflows
for end-to-end single-cell RNA-seq data analysis.
"""

# Import and expose key functions from submodules
from .standard import run_standard_workflow
from .custom import create_custom_workflow

# Define what should be accessible when importing from this module
__all__ = [
    "run_standard_workflow",
    "create_custom_workflow"
]