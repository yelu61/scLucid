"""Public API for spatial transcriptomics tools.

Legacy high-level workflows are preserved from ``_legacy.py``. New diagnostic
and statistical utilities are added under this package. Tumor-focused spatial
utilities are available in ``scLucid.tools.spatial.tumor``.
"""

from __future__ import annotations

from .autocorr import compute_moran_i, compute_spatial_autocorr
from .config import (
    SpatialAutocorrConfig,
    SpatialDiagnosticsConfig,
    SpatialNeighborsConfig,
    SpatialWindowConfig,
    SVGConfig,
    TissueZonesConfig,
    VisiumIOConfig,
)
from .diagnostics import diagnose_spatial_data_quality
from .neighbors import build_spatial_neighbors
from .subset import subset_spatial_window
from .svg import find_spatially_variable_genes
from .trace import build_spatial_review_summary
from .tumor import (
    analyze_spatial_niches,
    compute_immune_infiltration_score,
    find_tumor_stroma_boundary,
    spatial_ici_response_signature,
)
from .utils import infer_spatial_platform, validate_spatial_coords
from .visium import crop_visium, read_visium_10x, rotate_visium
from .zones import find_tissue_zones

# Legacy workflow entry points (imports are deferred inside functions)
from ._legacy import (
    export_spatial_report,
    plot_spatial,
    run_spatial_analysis,
    run_spatial_batch,
)

__all__ = [
    "analyze_spatial_niches",
    "build_spatial_neighbors",
    "build_spatial_review_summary",
    "compute_immune_infiltration_score",
    "compute_moran_i",
    "compute_spatial_autocorr",
    "crop_visium",
    "diagnose_spatial_data_quality",
    "export_spatial_report",
    "find_spatially_variable_genes",
    "find_tissue_zones",
    "find_tumor_stroma_boundary",
    "infer_spatial_platform",
    "plot_spatial",
    "read_visium_10x",
    "rotate_visium",
    "run_spatial_analysis",
    "run_spatial_batch",
    "spatial_ici_response_signature",
    "SpatialAutocorrConfig",
    "SpatialDiagnosticsConfig",
    "SpatialNeighborsConfig",
    "SpatialWindowConfig",
    "subset_spatial_window",
    "SVGConfig",
    "TissueZonesConfig",
    "validate_spatial_coords",
    "VisiumIOConfig",
]
