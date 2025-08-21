"""
Configuration classes for the Quality Control (QC) module.

This file defines the dataclass structures used to configure all aspects
of the QC pipeline, from threshold setting to doublet detection and filtering.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Literal, Optional, Union

__all__ = ["QCThresholds", "FilterConfig", "MarkerConfig", "DoubletConfig"]


@dataclass
class QCThresholds:
    """Configuration for QC filtering thresholds."""

    min_genes: Optional[int] = 200
    max_genes: Optional[int] = None
    min_counts: Optional[int] = None
    max_counts: Optional[int] = None
    pc_mt: Optional[float] = 20.0
    pc_hb: Optional[float] = 20.0
    pc_top_genes: Optional[float] = None
    nmads: float = 5.0
    use_fixed_top_gene_threshold: bool = False

    def to_dict(self) -> Dict[str, Any]:
        """Convert thresholds to dictionary format."""
        return self.__dict__


@dataclass
class FilterConfig:
    """Configuration for cell filtering logic."""

    criteria_to_filter: List[str] = field(
        default_factory=lambda: ["outlier_min_genes", "outlier_mt", "predicted_doublet"]
    )
    combination_logic: Literal["any", "all", "custom"] = "any"
    custom_logic_expr: Optional[str] = None


@dataclass
class MarkerConfig:
    """Configuration for a single marker set in heuristic doublet detection."""

    genes: Union[List[str], str]
    expression_threshold: float = 0.0
    min_genes_required: int = 1
    use_raw: bool = True
    is_regex: bool = field(init=False)

    def __post_init__(self):
        self.is_regex = isinstance(self.genes, str)


@dataclass
class DoubletConfig:
    """Comprehensive configuration for the doublet detection workflow."""

    method: str = "scrublet"
    merge_strategy: Literal[
        "union", "intersection", "algorithm_priority", "heuristic_priority"
    ] = "union"
    n_pcs: int = 30
    expected_doublet_rate: Optional[Union[float, Dict[str, float]]] = None
    use_heuristics: bool = True
    marker_configs: Optional[Dict[str, MarkerConfig]] = None
    min_lineages_for_doublet: int = 2
    marker_species: str = "human"
    marker_tissue: Optional[str] = None
    marker_states: Optional[List[str]] = None
    marker_level: Literal["major", "minor", "all"] = "major"
    plot_umap: bool = True
    export_stats: bool = True
    save_dir: Optional[str] = None
    show_plots: bool = True
