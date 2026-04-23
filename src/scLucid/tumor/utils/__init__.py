"""
Utilities for tumor analysis.

Provides:
- Tumor marker gene definitions
- Cancer hallmark signatures
- Database interfaces (COSMIC, TCGA, etc.)
"""

from .databases import (
    get_drug_targets,
    query_cancer_gene_census,
    query_tcga_data,
)
from .markers import (
    get_emt_markers,
    get_immune_markers,
    get_proliferation_markers,
    get_stromal_markers,
    get_tumor_markers,
)
from .signatures import (
    HallmarkCalculator,
    calculate_signature_scores,
    load_hallmark_signatures,
)

__all__ = [
    "get_tumor_markers",
    "get_immune_markers",
    "get_stromal_markers",
    "get_proliferation_markers",
    "get_emt_markers",
    "load_hallmark_signatures",
    "calculate_signature_scores",
    "HallmarkCalculator",
    "query_cancer_gene_census",
    "get_drug_targets",
    "query_tcga_data",
]
