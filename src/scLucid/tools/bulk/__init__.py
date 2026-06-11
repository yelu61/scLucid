"""Public API for bulk RNA-seq tools.

Bulk RNA-seq analysis utilities live here as a ``scLucid.tools`` package,
symmetric to ``scLucid.tools.spatial``. The previous home in
``scLucid.analysis.bulk`` is now a backward-compatible shim.
"""

from __future__ import annotations

from ._legacy import differential_abundance, run_deconvolution
from .abundance import correlate_abundance_with_clinical, run_bulk_abundance_test
from .clinical import correlate_abundance_with_clinical as correlate_abundance_with_clinical_alias
from .config import (
    BulkAbundanceConfig,
    BulkClinicalAssociationConfig,
    BulkDEConfig,
    BulkDeconvolutionConfig,
    BulkDiagnosticsConfig,
    BulkNormalizationConfig,
    BulkTraitAssociationConfig,
)
from .deconvolution import deconvolve_bulk
from .deg import run_bulk_de
from .diagnostics import diagnose_bulk_data_quality
from .normalize import estimate_size_factors_median_ratio, normalize_bulk_counts
from .trace import build_bulk_review_summary
from .tumor import (
    associate_tme_with_response,
    bulk_immune_landscape,
    deconvolve_tumor_tme,
    estimate_tumor_purity_from_bulk,
)
from .utils import deduplicate_var_names, filter_bulk_genes

__all__ = [
    "associate_tme_with_response",
    "bulk_immune_landscape",
    "BulkAbundanceConfig",
    "BulkClinicalAssociationConfig",
    "BulkDEConfig",
    "BulkDeconvolutionConfig",
    "BulkDiagnosticsConfig",
    "BulkNormalizationConfig",
    "BulkTraitAssociationConfig",
    "build_bulk_review_summary",
    "correlate_abundance_with_clinical",
    "correlate_abundance_with_clinical_alias",
    "deconvolve_bulk",
    "deconvolve_tumor_tme",
    "deduplicate_var_names",
    "diagnose_bulk_data_quality",
    "differential_abundance",
    "estimate_size_factors_median_ratio",
    "estimate_tumor_purity_from_bulk",
    "filter_bulk_genes",
    "normalize_bulk_counts",
    "run_bulk_abundance_test",
    "run_bulk_de",
    "run_deconvolution",
]
