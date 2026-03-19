"""Configuration for proportion analysis submodule."""

from typing import Dict, List, Literal, Optional, Tuple

from pydantic import Field, field_validator

from ...base_config import SclucidBaseConfig


class ProportionConfig(SclucidBaseConfig):
    """Configuration for cell type proportion analysis."""

    # Required fields
    celltype_col: str = Field(description="Column for cell types")
    sample_col: str = Field(description="Column for sample IDs")
    condition_col: str = Field(description="Column for conditions")

    # Optional fields
    pairing_col: Optional[str] = Field(default=None)
    batch_col: Optional[str] = Field(default=None)
    timepoint_col: Optional[str] = Field(default=None)

    auto_configure: bool = Field(default=True)
    test_method: Literal["deseq2", "t-test", "wilcoxon", "anova", "chi-square", "fisher"] = Field(
        default="wilcoxon"
    )

    # Plotting
    plot_types: List[str] = Field(
        default_factory=lambda: ["counts", "bar", "bar_composition", "box", "diff"]
    )

    # Palettes
    ct_palette: Optional[Dict[str, str]] = Field(default=None)
    condition_palette: Optional[Dict[str, str]] = Field(default=None)

    # Output
    out_dir: Optional[str] = Field(default=None)
    export_data: bool = Field(default=True)


class MethodSelectionConfig(SclucidBaseConfig):
    """
    Configuration for automatic method selection.

    Attributes
    ----------
    n_samples_per_group : int
        Number of samples per group
    n_celltypes : int
        Number of cell types
    has_batch_effect : bool
        Whether batch effect exists
    spatial_resolution : bool
        Whether spatial resolution is needed
    min_cells_per_sample : int
        Minimum cells per sample threshold
    celltype_count_threshold : int
        Cell type count threshold for low abundance filtering
    """

    n_samples_per_group: int = Field(
        default=5, ge=1, le=100, description="Number of samples per group"
    )
    n_celltypes: int = Field(default=10, ge=2, description="Number of cell types")
    has_batch_effect: bool = Field(
        default=False, description="Whether batch effect exists"
    )
    spatial_resolution: bool = Field(
        default=False, description="Whether spatial resolution is needed"
    )
    min_cells_per_sample: int = Field(
        default=100, ge=10, description="Minimum cells per sample threshold"
    )
    celltype_count_threshold: int = Field(
        default=10, ge=1, description="Cell type count threshold for low abundance filtering"
    )

    @field_validator("has_batch_effect")
    @classmethod
    def check_batch_effect(cls, v: bool) -> bool:
        """Validate batch effect parameter."""
        if v:
            import logging

            logging.getLogger(__name__).info(
                "Batch effect detection enabled, scCODA method recommended"
            )
        return v
