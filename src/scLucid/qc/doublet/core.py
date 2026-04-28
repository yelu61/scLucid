"""Doublet detection core — configuration helpers and simple public API.

Extracted for maintainability.
"""

from __future__ import annotations

import gc
import logging
from pathlib import Path
from typing import Dict, Literal, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData
from upsetplot import plot as upset_plot

from ...utils import get_marker_manager
from ..config import DoubletConfig, MarkerConfig

log = logging.getLogger(__name__)

# Constants for column naming consistency across submodules
LINEAGE_SCORES_KEY = "lineage_module_scores"
HEURISTIC_SCORE_COL = "heuristic_confidence_score"
HEURISTIC_PRED_COL = "heuristic_predicted"
FINAL_PRED_COL = "predicted_doublet"

__all__ = [
    "generate_doublet_rates",
    "create_custom_marker_dict",
]

def _create_doublet_marker_config_from_manager(
    adata: AnnData, cfg: DoubletConfig
) -> Dict[str, MarkerConfig]:
    """
    Create MarkerConfig objects using the marker manager and the main DoubletConfig.

    Args:
        adata: The AnnData object to intersect markers with.
        cfg: The main DoubletConfig object, providing context like species, tissue,
             and default evaluation parameters.

    Returns:
        A dictionary mapping lineage names to MarkerConfig objects.
    """
    try:
        # --- ❗ Select the correct adata object upfront ❗ ---
        # Ensure the intersection is performed on the same data that will be used for scoring.
        adata_for_intersection = (
            adata.raw.to_adata() if cfg.default_use_raw and adata.raw else adata
        )
        log.info(
            f"Performing marker intersection on {'adata.raw' if cfg.default_use_raw and adata.raw else 'adata'}."
        )

        case_sensitive = cfg.marker_species.lower() == "mouse"
        manager = get_marker_manager(
            species=cfg.marker_species,
            tissue=cfg.marker_tissue,
            case_sensitive=case_sensitive,
        )
        # Intersect with the correctly chosen data object
        manager.intersect_with(adata_for_intersection)
        markers_dict = manager.get_doublet_lineage_markers()

        if not markers_dict:
            log.warning(
                "No lineage markers found. Ensure `doublet_lineage = true` is set in your TOML file for desired cell types."
            )
            return {}

    except Exception as e:
        log.error(
            f"Failed to load markers via marker_manager: {e}. Cannot create heuristic configs."
        )
        return {}

    marker_configs = {}
    for lineage, genes in markers_dict.items():
        if genes:
            marker_configs[lineage] = MarkerConfig(
                genes=genes,
                expression_threshold=cfg.default_expression_threshold,
                min_genes_required=cfg.default_min_genes_required,
                use_raw=cfg.default_use_raw,
            )
    log.info(f"Auto-generated {len(marker_configs)} marker configurations for doublet detection.")
    return marker_configs



def generate_doublet_rates(
    adata: AnnData,
    sample_key: str = "sampleID",
    chemistry: str = "v3",
    custom_rate: float = 0.008,
    custom_rate_model: Literal["scale", "fixed"] = "scale",
    max_rate: float = 0.20,
    min_rate: float = 0.001,
) -> Dict[str, float]:
    """
    Automatically generate expected doublet rates based on cell count per sample or platform.

    This function supports three modes:
    1.  Platform Models ('v2', 'v3', 'HT'): Applies 10x-style linear/non-linear scaling
        based on cell count.
    2.  Fixed Models ('BD'): Applies a single, fixed rate to all samples.
    3.  Custom Models ('custom'):
        - If `custom_rate_model='scale'`, applies 10x-style scaling using `custom_rate`
          as the factor per 1000 cells.
        - If `custom_rate_model='fixed'`, applies `custom_rate` as a single fixed rate
          to all samples (ignoring cell counts).

    Args:
        adata: AnnData object containing cell count information.
        sample_key: Column name in adata.obs used to distinguish samples.
        chemistry: The technology platform ('v2', 'v3', 'HT', 'BD', or 'custom').
        custom_rate: The rate to use when `chemistry='custom'`. Interpreted based on
                     `custom_rate_model`. Defaults to 0.008.
        custom_rate_model: Defines how to use `custom_rate` ('scale' or 'fixed').
                           Defaults to 'scale'.
        max_rate: (For 'scale' models) Maximum doublet rate cap.
        min_rate: (For 'scale' models) Minimum doublet rate floor.

    Returns:
        Dictionary mapping sample IDs to calculated doublet rates.
    """
    log.info("Automatically generating doublet rates based on sample chemistry and cell counts...")

    # Define platform-specific models: (model_type, rate_value)
    # 'scale' model_type uses the rate_value as a scaling factor per 1000 cells.
    # 'fixed' model_type uses the rate_value as the final, fixed rate for all samples.
    chemistry_models = {
        "v2": ("scale", 0.007),  # 10x v2 chemistry
        "v3": ("scale", 0.008),  # 10x v3 chemistry
        "HT": ("scale", 0.016),  # 10x High-throughput
        "BD": ("fixed", 0.025),  # BD Rhapsody
    }

    cell_counts = adata.obs[sample_key].value_counts()
    doublet_rates = {}

    model_type = None
    rate_value = None

    # --- 1. Determine Model Type and Rate ---
    if chemistry in chemistry_models:
        model_type, rate_value = chemistry_models[chemistry]
        log.info(
            f"Using known platform model '{chemistry}': type='{model_type}', base_rate={rate_value}"
        )
    elif chemistry == "custom":
        model_type = custom_rate_model
        rate_value = custom_rate
        log.info(f"Using 'custom' model: type='{model_type}', custom_rate={rate_value}")
    else:
        log.warning(f"Unknown chemistry '{chemistry}'. Falling back to default 'v3' model.")
        model_type, rate_value = chemistry_models["v3"]
        chemistry = "v3"  # Set chemistry for scaling logic

    # --- 2. Apply Model to Calculate Rates ---
    if model_type == "fixed":
        log.info(
            f"Applying fixed doublet rate of {rate_value:.4f} to all {len(cell_counts)} samples."
        )
        for sample, n_cells in cell_counts.items():
            doublet_rates[sample] = rate_value
            log.info(f"  - Sample '{sample}': {n_cells} cells -> Doublet rate: {rate_value:.4f}")

    elif model_type == "scale":
        log.info(f"Applying scaling model with base rate of {rate_value:.4f} per 1000 cells.")
        for sample, n_cells in cell_counts.items():
            # Use original logic for 10x non-linear scaling at high cell counts
            if n_cells > 10000 and chemistry in ["v2", "v3"]:
                rate = (0.8 * (n_cells / 1000) * rate_value) + (
                    0.2 * rate_value * (n_cells / 1000) ** 1.5
                )
            else:
                # Standard linear scaling for 'HT', 'custom_scale', and lower counts
                rate = (n_cells / 1000) * rate_value

            # Apply rate constraints
            rate = max(min_rate, min(rate, max_rate))
            doublet_rates[sample] = rate
            log.info(f"  - Sample '{sample}': {n_cells} cells -> Doublet rate: {rate:.4f}")

    else:
        # This case should not be reachable if logic is sound
        raise ValueError(f"Internal error: Unrecognized model_type '{model_type}'")

    return doublet_rates


def create_custom_marker_dict(
    lineage_definitions: Dict[str, Dict], save_path: Optional[Union[str, Path]] = None
) -> Dict[str, MarkerConfig]:
    """
    Create custom marker dictionary from user-defined lineage specifications.

    This function allows users to define their own marker sets with custom
    parameters for specialized doublet detection scenarios.

    Args:
        lineage_definitions: Dictionary defining lineages and their parameters
        save_path: Optional path to save the configuration for future use

    Returns:
        Dictionary mapping lineage names to MarkerConfig objects

    Example:
        lineage_defs = {
            "T_cells": {
                "genes": ["CD3D", "CD3E", "CD8A"],
                "expression_threshold": 0.5,
                "min_genes_required": 1
            },
            "Epithelial": {
                "genes": r"^KRT[0-9]+",  # Regex pattern
                "expression_threshold": 1.0,
                "min_genes_required": 2
            }
        }
        marker_configs = create_custom_marker_dict(lineage_defs)
    """
    config_dict = {}

    for lineage, definition in lineage_definitions.items():
        # Validate required 'genes' field
        if "genes" not in definition:
            raise ValueError(f"Missing 'genes' field for lineage '{lineage}'")

        config_dict[lineage] = MarkerConfig(**definition)

    # Save configuration if requested
    if save_path:
        import json

        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)

        # Convert to serializable format
        serializable_dict = {}
        for lineage, config in config_dict.items():
            serializable_dict[lineage] = {
                "genes": config.genes,
                "expression_threshold": config.expression_threshold,
                "min_genes_required": config.min_genes_required,
                "use_raw": config.use_raw,
            }

        with open(save_path, "w") as f:
            json.dump(serializable_dict, f, indent=2)
        log.info(f"Marker configuration saved to {save_path}")

    return config_dict


