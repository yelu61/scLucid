"""
Utilities for tumor analysis.

Provides:
- Tumor marker gene definitions
- Cancer hallmark signatures
- Database interfaces (COSMIC, TCGA, etc.)
"""

from .markers import (
    get_tumor_markers,
    get_immune_markers,
    get_stromal_markers,
    get_proliferation_markers,
    get_emt_markers,
)

from .signatures import (
    load_hallmark_signatures,
    calculate_signature_scores,
    HallmarkCalculator,
)

from .databases import (
    query_cancer_gene_census,
    get_drug_targets,
    query_tcga_data,
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
