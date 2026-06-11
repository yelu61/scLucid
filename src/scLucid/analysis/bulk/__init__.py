"""Backward-compatible shim for bulk RNA-seq analysis.

The canonical implementation has moved to ``scLucid.tools.bulk``. This module
re-exports the public API so existing imports continue to work.
"""

from __future__ import annotations

from scLucid.tools.bulk import (  # noqa: F401
    BulkAbundanceConfig,
    BulkClinicalAssociationConfig,
    BulkDEConfig,
    BulkDeconvolutionConfig,
    BulkDiagnosticsConfig,
    BulkNormalizationConfig,
    BulkTraitAssociationConfig,
    build_bulk_review_summary,
    correlate_abundance_with_clinical,
    deconvolve_bulk,
    deduplicate_var_names,
    diagnose_bulk_data_quality,
    differential_abundance,
    estimate_size_factors_median_ratio,
    filter_bulk_genes,
    normalize_bulk_counts,
    run_bulk_abundance_test,
    run_bulk_de,
    run_deconvolution,
)

__all__ = [
    "BulkAbundanceConfig",
    "BulkClinicalAssociationConfig",
    "BulkDEConfig",
    "BulkDeconvolutionConfig",
    "BulkDiagnosticsConfig",
    "BulkNormalizationConfig",
    "BulkTraitAssociationConfig",
    "build_bulk_review_summary",
    "correlate_abundance_with_clinical",
    "deconvolve_bulk",
    "deduplicate_var_names",
    "diagnose_bulk_data_quality",
    "differential_abundance",
    "estimate_size_factors_median_ratio",
    "filter_bulk_genes",
    "normalize_bulk_counts",
    "run_bulk_abundance_test",
    "run_bulk_de",
    "run_deconvolution",
]
