"""
scRNA-toolkit: A comprehensive toolkit for single-cell RNA-seq data analysis
"""

__version__ = "0.1.0"

# Import the main functions from each module
from .qc import calculate_qc_metric, is_low_quality_cell, is_doublet
from .norm import normalize_data
from .hvg import annotate_hvg
from .integrate import integration_scanorama, integration_harmony, integration_scvi
from .infercnv import run_cnv_analysis, find_tumor
from .helper import identify_outliers, merge_data, seed_everything
from .workflow import run_standard_workflow, custom_workflow

# Default Configuration
from .config import DEFAULT_CONFIG, load_config