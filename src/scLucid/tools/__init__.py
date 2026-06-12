"""Specialized tools public API for scLucid."""

import warnings
from collections.abc import Iterable
from importlib import import_module
from importlib.util import find_spec

__all__ = []


def _export(
    module: str,
    names: Iterable[str],
    *,
    optional: bool = True,
    requires: Iterable[str] = (),
) -> bool:
    """Import names from a backend submodule without breaking package import.

    Parameters
    ----------
    module : str
        Submodule name under ``scLucid.tools``.
    names : Iterable[str]
        Public symbols to re-export.
    optional : bool, default=True
        Tag for the warning message when an unexpected error occurs.
    requires : Iterable[str], default=()
        External pip packages this backend depends on. If any are not
        installed, the backend is skipped **silently** — declared optional
        deps are not a defect. A warning is still raised if the backend's
        own internals fail to import (i.e., a real bug).
    """
    for dep in requires:
        if find_spec(dep) is None:
            return False

    try:
        loaded = import_module(f"{__name__}.{module}")
    except Exception as exc:
        level = "optional" if optional else "required"
        warnings.warn(
            f"Could not import {level} tools backend '{module}': {exc}",
            ImportWarning,
        )
        return False

    found = False
    for name in names:
        if hasattr(loaded, name):
            globals()[name] = getattr(loaded, name)
            __all__.append(name)
            found = True
    return found


# Python-native tools
_export(
    "bulk",
    [
        "deconvolve_bulk",
        "run_deconvolution",
        "differential_abundance",
        "correlate_abundance_with_clinical",
        "run_bulk_abundance_test",
        "run_bulk_de",
        "diagnose_bulk_data_quality",
        "normalize_bulk_counts",
        "estimate_size_factors_median_ratio",
        "deduplicate_var_names",
        "filter_bulk_genes",
        "build_bulk_review_summary",
        "BulkDiagnosticsConfig",
        "BulkNormalizationConfig",
        "BulkDEConfig",
        "BulkDeconvolutionConfig",
        "BulkAbundanceConfig",
        "BulkClinicalAssociationConfig",
        "BulkTraitAssociationConfig",
    ],
)
_export(
    "cellphonedb",
    [
        "run_cellphonedb",
        "run_cellphonedb_batch",
        "run_cellphonedb_by_group",
        "summarize_cellphonedb",
    ],
)
_export(
    "spatial",
    [
        "diagnose_spatial_data_quality",
        "build_spatial_neighbors",
        "compute_moran_i",
        "compute_spatial_autocorr",
        "find_spatially_variable_genes",
        "find_tissue_zones",
        "subset_spatial_window",
        "crop_visium",
        "rotate_visium",
        "read_visium_10x",
        "build_spatial_review_summary",
        "SpatialDiagnosticsConfig",
        "SpatialNeighborsConfig",
        "SpatialAutocorrConfig",
        "SVGConfig",
        "TissueZonesConfig",
        "SpatialWindowConfig",
        "VisiumIOConfig",
        "run_spatial_analysis",
        "plot_spatial",
        "run_spatial_batch",
        "export_spatial_report",
    ],
)
_export(
    "pySCENIC",
    [
        "analyze_scenic_results",
        "export_scenic_report",
        "run_scenic",
        "run_scenic_batch",
        "run_scenic_by_group",
    ],
)
_export(
    "pyMonocle3",
    [
        "CellDataSet",
        "new_cell_data_set",
        "create_cds_from_scanpy",
        "export_to_scanpy",
        "preprocess_cds",
        "reduce_dimension",
        "cluster_cells",
        "learn_graph",
        "order_cells",
        "graph_test",
        "top_markers",
        "plot_cells",
    ],
)
_export(
    "pyCellChat",
    [
        "CellChat",
        "CellChatDB",
        "get_default_database",
        "create_cellchat_from_scanpy",
        "plot_heatmap",
    ],
)
_export(
    "pyBayesPrism",
    [
        "PrismConfig",
        "BayesPrismReference",
        "BayesPrism",
        "BayesPrismEmbedding",
        "GibbsSampler",
        "plot_fraction",
        "plot_correlation",
        "cleanup_genes",
        "compute_correlation",
        "compute_rmse",
    ],
)
_export(
    "pyDWLS",
    [
        "DWLS",
        "SignatureBuilder",
        "DampenedWLS",
        "MarkerSelector",
        "CrossValidator",
        "solve_nnls",
        "normalize_data",
        "filter_genes",
        "create_pseudo_bulk",
    ],
)
