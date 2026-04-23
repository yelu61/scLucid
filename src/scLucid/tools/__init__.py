"""Specialized tools public API for scLucid."""

import warnings
from collections.abc import Iterable
from importlib import import_module

__all__ = []


def _export(module: str, names: Iterable[str], *, optional: bool = True) -> bool:
    """Import names from a backend submodule without breaking package import."""
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
_export("infercnv", ["find_tumor", "run_cnv_analysis"])
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
