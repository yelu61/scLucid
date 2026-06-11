"""Pydantic configurations for spatial transcriptomics tools."""

from __future__ import annotations

from typing import List, Literal, Optional, Tuple

from pydantic import Field

from ...base_config import SclucidBaseConfig


class SpatialDiagnosticsConfig(SclucidBaseConfig):
    """Configuration for spatial data quality diagnostics."""

    model_config = SclucidBaseConfig.model_config

    spatial_key: str = Field(default="spatial")
    require_image: bool = Field(default=False)
    check_duplicate_coords: bool = Field(default=True)
    min_spots: int = Field(default=10, ge=2)


class SpatialNeighborsConfig(SclucidBaseConfig):
    """Configuration for building spatial neighbor graphs."""

    model_config = SclucidBaseConfig.model_config

    spatial_key: str = Field(default="spatial")
    method: Literal["knn", "radius"] = Field(default="knn")
    n_neigh: int = Field(default=6, ge=1)
    radius: Optional[float] = Field(default=None, ge=0)
    key_added: str = Field(default="spatial_neighbors")


class SpatialAutocorrConfig(SclucidBaseConfig):
    """Configuration for spatial autocorrelation."""

    model_config = SclucidBaseConfig.model_config

    spatial_key: str = Field(default="spatial")
    mode: Literal["moran", "geary"] = Field(default="moran")
    n_permutations: int = Field(default=0, ge=0)
    key_added: str = Field(default="moran_i")


class SVGConfig(SclucidBaseConfig):
    """Configuration for spatially variable gene detection."""

    model_config = SclucidBaseConfig.model_config

    spatial_key: str = Field(default="spatial")
    method: Literal["moran_i", "pearsonr", "prost"] = Field(default="moran_i")
    n_permutations: int = Field(default=100, ge=0)
    alpha: float = Field(default=0.05, gt=0, lt=1)
    layer: Optional[str] = Field(default=None)
    key_added: str = Field(default="spatially_variable")


class TissueZonesConfig(SclucidBaseConfig):
    """Configuration for tissue zone detection."""

    model_config = SclucidBaseConfig.model_config

    n_components: int = Field(default=5, ge=2)
    method: Literal["nmf"] = Field(default="nmf")
    input: Literal["expression", "deconvolution"] = Field(default="deconvolution")
    key_added: str = Field(default="tissue_zones")


class VisiumIOConfig(SclucidBaseConfig):
    """Configuration for Visium IO helpers."""

    model_config = SclucidBaseConfig.model_config

    library_id: Optional[str] = Field(default=None)
    load_images: bool = Field(default=True)


class SpatialWindowConfig(SclucidBaseConfig):
    """Configuration for spatial window subsetting."""

    model_config = SclucidBaseConfig.model_config

    spatial_key: str = Field(default="spatial")
    xlim: Optional[Tuple[float, float]] = Field(default=None)
    ylim: Optional[Tuple[float, float]] = Field(default=None)
