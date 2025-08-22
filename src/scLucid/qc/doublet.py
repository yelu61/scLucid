"""
Enhanced doublet detection for single-cell RNA-seq data.

This module provides comprehensive functions for identifying potential doublet cells
using multiple algorithmic methods and flexible heuristic approaches based on
mutually exclusive lineage marker co-expression. It integrates with a unified
marker management system for consistent and reproducible analysis.
"""

import gc
import json
import logging
from pathlib import Path
from typing import Dict, List, Literal, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
from anndata import AnnData

from ..utils.marker_manager import get_marker_manager
from .config import DoubletConfig, MarkerConfig

log = logging.getLogger(__name__)

# --- Use constants for column names for easier maintenance ---
HEURISTIC_PRED_COL = "heuristic_predicted"
FINAL_PRED_COL = "predicted_doublet"

__all__ = [
    "generate_doublet_rates",
    "predict_doublets",
    "export_doublet_stats",
    "create_custom_marker_dict",
    "create_doublet_marker_config_from_manager",
]


# --- Helper Functions ---
def _get_builtin_markers(species: str, marker_type: str) -> Dict[str, MarkerConfig]:
    """
    Fallback built-in marker definitions when marker manager is unavailable.

    This function provides basic marker sets as a backup when the unified
    marker management system cannot be loaded.

    Args:
        species: Species name ("human" or "mouse")
        marker_type: Type of markers ("major_lineages" or "detailed")

    Returns:
        Dictionary mapping lineage names to MarkerConfig objects
    """
    if species.lower() == "human" and marker_type == "major_lineages":
        return {
            "T_cells": MarkerConfig(
                genes=["CD3D", "CD3E", "CD3G", "CD8A", "CD4"],
                expression_threshold=0.5,
                min_genes_required=1,
            ),
            "B_cells": MarkerConfig(
                genes=["CD19", "MS4A1", "CD79A", "CD79B"],
                expression_threshold=0.5,
                min_genes_required=1,
            ),
            "Myeloid": MarkerConfig(
                genes=["CD14", "CD68", "LYZ", "S100A9", "FCGR3A"],
                expression_threshold=0.5,
                min_genes_required=1,
            ),
            "NK_cells": MarkerConfig(
                genes=["KLRD1", "KLRF1", "NCR1", "GNLY", "NKG7"],
                expression_threshold=0.5,
                min_genes_required=1,
            ),
            "Epithelial": MarkerConfig(
                genes=r"^(KRT|CK)[0-9]+",  # Regex pattern for keratins
                expression_threshold=1.0,
                min_genes_required=2,
            ),
            "Endothelial": MarkerConfig(
                genes=["PECAM1", "VWF", "CDH5", "PLVAP"],
                expression_threshold=0.5,
                min_genes_required=1,
            ),
        }

    elif species.lower() == "mouse" and marker_type == "major_lineages":
        return {
            "T_cells": MarkerConfig(
                genes=["Cd3d", "Cd3e", "Cd3g", "Cd8a", "Cd4"],
                expression_threshold=0.5,
                min_genes_required=1,
            ),
            "B_cells": MarkerConfig(
                genes=["Cd19", "Ms4a1", "Cd79a", "Cd79b"],
                expression_threshold=0.5,
                min_genes_required=1,
            ),
            "Myeloid": MarkerConfig(
                genes=["Cd14", "Cd68", "Lyz2", "S100a9"],
                expression_threshold=0.5,
                min_genes_required=1,
            ),
            "Endothelial": MarkerConfig(
                genes=["Pecam1", "Vwf", "Cdh5"],
                expression_threshold=0.5,
                min_genes_required=1,
            ),
        }

    else:
        log.warning(
            f"No built-in markers for {species} {marker_type}, returning empty dict"
        )
        return {}


def _load_custom_marker_dict(filepath: Union[str, Path]) -> Dict[str, MarkerConfig]:
    """
    Load markers from custom file (JSON/YAML format).

    Args:
        filepath: Path to custom marker configuration file

    Returns:
        Dictionary mapping lineage names to MarkerConfig objects

    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the file format is unsupported
    """

    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Marker file not found: {filepath}")

    with open(filepath, "r") as f:
        if filepath.suffix.lower() == ".json":
            data = json.load(f)
        elif filepath.suffix.lower() in [".yaml", ".yml"]:
            import yaml

            data = yaml.safe_load(f)
        else:
            raise ValueError(f"Unsupported file format: {filepath.suffix}")

    # Convert to MarkerConfig objects
    config_dict = {}
    for lineage, marker_def in data.items():
        if isinstance(marker_def, dict):
            config_dict[lineage] = MarkerConfig(**marker_def)
        else:
            config_dict[lineage] = MarkerConfig(genes=marker_def)

    return config_dict


def _evaluate_lineage_expression(
    adata: AnnData, lineage_name: str, marker_config: MarkerConfig
) -> pd.Series:
    """
    Evaluate expression of markers for a specific lineage.

    This function processes marker genes (either as a list or regex pattern)
    and determines which cells are positive for the given lineage based on
    the specified expression threshold and minimum gene requirements.

    Args:
        adata: AnnData object containing expression data
        lineage_name: Name of the lineage being evaluated
        marker_config: Configuration object defining marker evaluation parameters

    Returns:
        Boolean series indicating cells positive for this lineage
    """
    source_adata = (
        adata.raw.to_adata() if marker_config.use_raw and adata.raw else adata
    )

    if marker_config.is_regex:
        # Handle regex pattern for gene matching
        pattern = marker_config.genes
        matching_genes = source_adata.var_names.str.contains(
            pattern, regex=True, na=False
        )
        valid_genes = source_adata.var_names[matching_genes].tolist()

        if not valid_genes:
            log.warning(
                f"No genes found matching pattern '{pattern}' for lineage '{lineage_name}'"
            )
            return pd.Series(False, index=adata.obs_names)

        log.info(
            f"Found {len(valid_genes)} genes matching pattern '{pattern}' for lineage '{lineage_name}'"
        )

    else:
        # Handle explicit gene list
        valid_genes = [g for g in marker_config.genes if g in source_adata.var_names]
        if not valid_genes:
            log.warning(f"No valid genes found for lineage '{lineage_name}'")
            return pd.Series(False, index=adata.obs_names)

        missing_genes = set(marker_config.genes) - set(valid_genes)
        if missing_genes:
            log.warning(f"Missing genes for lineage '{lineage_name}': {missing_genes}")

    # Extract expression data for valid genes
    expr_data = sc.get.obs_df(source_adata[:, valid_genes], keys=valid_genes)

    # Apply expression threshold to determine gene-level positivity
    expr_binary = expr_data > marker_config.expression_threshold

    # Count number of positive genes per cell
    genes_expressed = expr_binary.sum(axis=1)

    # Determine lineage positivity based on minimum gene requirement
    lineage_positive = genes_expressed >= marker_config.min_genes_required

    positive_count = lineage_positive.sum()
    positive_percentage = positive_count / len(lineage_positive) * 100

    log.info(
        f"Lineage '{lineage_name}': {positive_count} cells positive "
        f"({positive_percentage:.2f}%)"
    )

    return lineage_positive


def _identify_coexpression_doublets(
    adata: AnnData,
    marker_configs: Dict[str, MarkerConfig],
    min_lineages_for_doublet: int = 2,
) -> pd.Series:
    """
    Enhanced co-expression doublet detection with flexible marker configurations.

    This function identifies potential doublets by detecting cells that co-express
    markers from multiple mutually exclusive lineages. This approach is based on
    the biological principle that genuine single cells should predominantly express
    markers from a single lineage.

    Args:
        adata: AnnData object containing expression data
        marker_configs: Dictionary mapping lineage names to MarkerConfig objects
        min_lineages_for_doublet: Minimum number of lineages that must be co-expressed

    Returns:
        Boolean series indicating potential doublet cells
    """
    log.info("Identifying potential doublets via co-expression of exclusive markers...")

    if not marker_configs:
        log.warning(
            "No marker configurations provided. Skipping co-expression heuristic."
        )
        return pd.Series(False, index=adata.obs_names)

    # Evaluate each lineage independently
    lineage_results = {}
    for lineage_name, marker_config in marker_configs.items():
        lineage_results[lineage_name] = _evaluate_lineage_expression(
            adata, lineage_name, marker_config
        )

    # Combine lineage results into a DataFrame
    lineage_df = pd.DataFrame(lineage_results, index=adata.obs_names)

    # Count number of lineages positive per cell
    lineages_per_cell = lineage_df.sum(axis=1)

    # Identify doublets based on co-expression threshold
    coexpression_doublets = lineages_per_cell >= min_lineages_for_doublet

    doublet_count = coexpression_doublets.sum()
    log.info(
        f"Found {doublet_count} potential doublets co-expressing "
        f"markers from >= {min_lineages_for_doublet} lineages."
    )

    # Log detailed statistics for different co-expression levels
    for n_lineages in range(min_lineages_for_doublet, lineage_df.shape[1] + 1):
        count = (lineages_per_cell == n_lineages).sum()
        if count > 0:
            percentage = count / len(lineages_per_cell) * 100
            log.info(
                f"  - Cells expressing {n_lineages} lineages: {count} ({percentage:.2f}%)"
            )

    return coexpression_doublets


def _merge_doublet_predictions(
    adata: AnnData,
    algorithm_col: str,
    heuristic_col: str = HEURISTIC_PRED_COL,
    strategy: str = "union",
) -> pd.Series:
    """
    Merge algorithmic and heuristic doublet predictions using different strategies.

    This function combines predictions from algorithmic methods (e.g., Scrublet)
    with heuristic marker-based predictions using various logical operations.

    Args:
        adata: AnnData object containing prediction results
        algorithm_col: Column name for algorithmic predictions
        heuristic_col: Column name for heuristic predictions
        strategy: Merge strategy ("union", "intersection", "algorithm_priority", "heuristic_priority")

    Returns:
        Boolean series with merged predictions
    """
    algo_pred = adata.obs[algorithm_col].fillna(False).astype(bool)
    heur_pred = adata.obs[heuristic_col].fillna(False).astype(bool)

    if strategy == "union":
        merged = algo_pred | heur_pred
        log.info("Using union strategy: doublets predicted by either method")

    elif strategy == "intersection":
        merged = algo_pred & heur_pred
        log.info("Using intersection strategy: doublets predicted by both methods")

    elif strategy == "algorithm_priority":
        merged = algo_pred.copy()
        log.info("Using algorithm priority: only algorithmic predictions")

    elif strategy == "heuristic_priority":
        merged = heur_pred.copy()
        log.info("Using heuristic priority: only heuristic predictions")

    else:
        raise ValueError(f"Unknown merge strategy: {strategy}")

    # Log detailed merge statistics
    algo_count = algo_pred.sum()
    heur_count = heur_pred.sum()
    merged_count = merged.sum()
    overlap_count = (algo_pred & heur_pred).sum()

    total_cells = len(algo_pred)

    log.info("Merge statistics:")
    log.info(
        f"  - Algorithm predictions: {algo_count} ({algo_count / total_cells:.2%})"
    )
    log.info(
        f"  - Heuristic predictions: {heur_count} ({heur_count / total_cells:.2%})"
    )
    log.info(f"  - Overlap: {overlap_count} ({overlap_count / total_cells:.2%})")
    log.info(f"  - Final merged: {merged_count} ({merged_count / total_cells:.2%})")

    return merged


def _load_marker_dict(
    species: str = "human",
    marker_type: str = "major_lineages",
    custom_path: Optional[Union[str, Path]] = None,
    expression_threshold: float = 0.5,
    min_genes_required: int = 1,
    use_raw: bool = True,
) -> Dict[str, MarkerConfig]:
    """
    Load predefined marker dictionaries using the unified marker manager.

    This function leverages the comprehensive marker management system to load
    appropriate marker sets for doublet detection. It supports both simple
    marker loading and custom configurations.

    Args:
        species: Species name ("human", "mouse")
        marker_type: Type of markers ("major_lineages", "detailed")
        custom_path: Path to custom marker file (overrides species/type)
        expression_threshold: Expression threshold for all markers
        min_genes_required: Minimum genes required for positive lineage
        use_raw: Whether to use raw expression data

    Returns:
        Dictionary mapping lineage names to MarkerConfig objects
    """
    if custom_path:
        return _load_custom_marker_dict(custom_path)

    try:
        # Load the appropriate marker set using unified manager
        if marker_type == "major_lineages":
            manager = get_marker_manager(species=species, case_sensitive=False)
            # Get only major lineages suitable for doublet detection
            markers_dict = manager.get_markers_by_level("major")
        elif marker_type == "detailed":
            manager = get_marker_manager(species=species, case_sensitive=False)
            # Get all markers including subtypes for comprehensive analysis
            markers_dict = manager.get_markers_by_level("all")
        else:
            raise ValueError(f"Unknown marker_type: {marker_type}")

        # Convert to MarkerConfig objects with specified parameters
        marker_configs = {}
        for lineage, genes in markers_dict.items():
            if genes:  # Only include lineages with actual markers
                marker_configs[lineage] = MarkerConfig(
                    genes=genes,
                    expression_threshold=expression_threshold,
                    min_genes_required=min_genes_required,
                    use_raw=use_raw,
                )

        log.info(
            f"Loaded {len(marker_configs)} marker configurations for {species} {marker_type}"
        )
        return marker_configs

    except ImportError:
        log.warning("Marker manager not found, using built-in markers")
        return _get_builtin_markers(species, marker_type)
    except Exception as e:
        log.error(f"Failed to load markers: {e}")
        return _get_builtin_markers(species, marker_type)


def _run_scrublet(
    adata_view: AnnData,
    sample_name: str,
    config: DoubletConfig,
) -> tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Runs the Scrublet algorithm on a single AnnData view."""
    rate = config.expected_doublet_rate
    current_rate = rate.get(sample_name, 0.1) if isinstance(rate, dict) else rate
    actual_n_pcs = min(config.n_pcs, adata_view.n_obs - 1, adata_view.n_vars - 1)

    try:
        scrub = scr.Scrublet(adata_view.X, expected_doublet_rate=current_rate)
        scores, _ = scrub.scrub_doublets(n_prin_comps=actual_n_pcs, verbose=False)
        predicted = scrub.call_doublets(verbose=False)

        doublet_count = sum(predicted)
        doublet_rate = doublet_count / len(predicted)
        log.info(
            f"  Found {doublet_count} potential doublets via Scrublet ({doublet_rate:.2%})"
        )

        if config.plot_umap:
            try:
                scrub.set_embedding(
                    "UMAP", scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3)
                )
                fig = scrub.plot_embedding("UMAP", order_points=True)
                if config.save_dir:
                    save_path = (
                        Path(config.save_dir) / f"{sample_name}_doublets_umap.png"
                    )
                    save_path.parent.mkdir(parents=True, exist_ok=True)
                    fig.savefig(save_path, dpi=300, bbox_inches="tight")
                if config.show_plots:
                    plt.show()
                else:
                    plt.close(fig)
            except Exception as e:
                log.warning(f"Could not generate UMAP for sample {sample_name}: {e}")

        return scores, predicted

    except Exception as e:
        log.error(f"Scrublet failed for sample {sample_name}: {e}")
        return None, None
    finally:
        gc.collect()


# --- Main Functions ---
def generate_doublet_rates(
    adata: AnnData,
    sample_key: str = "sampleID",
    rate_per_1000_cells: float = 0.008,
    max_rate: float = 0.20,
    min_rate: float = 0.001,
) -> Dict[str, float]:
    """
    Automatically generate expected doublet rates based on cell count per sample.

    This function calculates expected doublet rates using the 10x Genomics guideline:
    for every 1000 cells, the multiplet rate increases by approximately 0.8% (0.008).
    This accounts for the fact that higher cell loading increases collision probability.

    Args:
        adata: AnnData object containing cell count information
        sample_key: Column name in adata.obs used to distinguish samples
        rate_per_1000_cells: Expected doublet rate per 1000 cells
                           - 0.008 for standard 3' v3.1 chemistry
                           - 0.016 for high-throughput (HT) kits
        max_rate: Maximum doublet rate cap (prevents unrealistic rates)
        min_rate: Minimum doublet rate floor (ensures some detection sensitivity)

    Returns:
        Dictionary mapping sample IDs to calculated doublet rates

    Example:
        >>> doublet_rates = generate_doublet_rates(adata, sample_key="sampleID")
        >>> print(doublet_rates)
        {'sample_A': 0.04, 'sample_B': 0.08}
    """
    log.info(
        "Automatically generating doublet rates based on cell counts per sample..."
    )

    # Calculate cell counts per sample
    cell_counts = adata.obs[sample_key].value_counts()
    doublet_rates = {}

    for sample, n_cells in cell_counts.items():
        # Apply 10x Genomics linear scaling formula
        rate = (n_cells / 1000) * rate_per_1000_cells

        # Apply rate constraints to prevent unrealistic values
        rate = max(min_rate, min(rate, max_rate))

        doublet_rates[sample] = rate
        log.info(f"  - Sample '{sample}': {n_cells} cells -> Doublet rate: {rate:.4f}")

    return doublet_rates


def create_doublet_marker_config_from_manager(
    species: str = "human",
    tissue: Optional[str] = None,
    states: Optional[List[str]] = None,
    level: Literal["major", "minor", "all"] = "major",
    expression_threshold: float = 0.5,
    min_genes_required: int = 1,
    use_raw: bool = True,
    min_markers_per_type: int = 1,
) -> Dict[str, MarkerConfig]:
    """
    Create MarkerConfig objects using the comprehensive marker manager.

    This function builds marker configurations by combining base, tissue-specific,
    and cell state-specific markers for enhanced doublet detection accuracy.

    Args:
        species: Species name ("human", "mouse")
        tissue: Specific tissue context (e.g., "Lung", "Brain")
        states: Cell states to include (e.g., ["Proliferating", "Hypoxia"])
        level: Cell type level to include ("major", "minor", "all")
        expression_threshold: Expression threshold for markers
        min_genes_required: Minimum genes required for positive lineage
        use_raw: Whether to use raw counts
        min_markers_per_type: Minimum markers required per cell type (filters out sparse types)

    Returns:
        Dictionary mapping lineage names to MarkerConfig objects
    """
    # Build comprehensive marker manager with all requested contexts
    manager = get_marker_manager(
        species=species, tissue=tissue, states=states, case_sensitive=False
    )

    # Filter out cell types with insufficient markers
    if min_markers_per_type > 1:
        removed_types = manager.filter_markers(min_genes_per_type=min_markers_per_type)
        if removed_types:
            log.info(
                f"Filtered out {len(removed_types)} cell types with < {min_markers_per_type} markers"
            )

    # Get markers by specified level
    markers_dict = manager.get_markers_by_level(level)

    # Convert to MarkerConfig objects
    marker_configs = {}
    for lineage, genes in markers_dict.items():
        if len(genes) >= min_markers_per_type:
            marker_configs[lineage] = MarkerConfig(
                genes=genes,
                expression_threshold=expression_threshold,
                min_genes_required=min_genes_required,
                use_raw=use_raw,
            )

    log.info(f"Created {len(marker_configs)} marker configurations from manager")
    return marker_configs


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


def export_doublet_stats(
    adata: AnnData,
    sample_key: str = "sampleID",
    save_dir: Optional[Union[str, Path]] = None,
    export_csv: bool = True,
    export_xlsx: bool = False,
) -> Dict[str, pd.DataFrame]:
    """
    Export comprehensive doublet statistics per sample and globally.

    This function generates detailed statistical summaries of doublet detection
    results, including counts, percentages, and score distributions.

    Args:
        adata: AnnData object with doublet predictions
        sample_key: Key for sample identification
        save_dir: Directory to save statistics files
        export_csv: Whether to export as CSV files
        export_xlsx: Whether to export as Excel file

    Returns:
        Dictionary containing sample-wise and global statistics DataFrames
    """
    # Identify all doublet-related columns
    doublet_cols = [
        col
        for col in adata.obs.columns
        if any(
            keyword in col.lower() for keyword in ["doublet", "scrublet", "heuristic"]
        )
    ]

    if not doublet_cols:
        log.warning("No doublet-related columns found in adata.obs")
        return {}

    log.info(f"Found doublet columns: {doublet_cols}")

    # Calculate per-sample statistics
    sample_stats = []
    for sample in adata.obs[sample_key].unique():
        sample_mask = adata.obs[sample_key] == sample
        sample_data = adata.obs[sample_mask]

        stats = {
            "sample": sample,
            "total_cells": len(sample_data),
        }

        # Process each doublet-related column
        for col in doublet_cols:
            if col in sample_data.columns:
                col_data = sample_data[col].dropna()

                if col_data.dtype == "bool" or set(col_data.unique()).issubset(
                    {0, 1, True, False}
                ):
                    # Boolean/binary column (predictions)
                    positive_count = col_data.sum()
                    stats[f"{col}_count"] = positive_count
                    stats[f"{col}_percentage"] = positive_count / len(sample_data) * 100
                else:
                    # Continuous column (scores)
                    stats[f"{col}_mean"] = col_data.mean()
                    stats[f"{col}_median"] = col_data.median()
                    stats[f"{col}_std"] = col_data.std()
                    stats[f"{col}_q25"] = col_data.quantile(0.25)
                    stats[f"{col}_q75"] = col_data.quantile(0.75)

        sample_stats.append(stats)

    sample_df = pd.DataFrame(sample_stats)

    # Calculate global statistics
    global_stats = {"metric": "global", "total_cells": adata.n_obs}

    for col in doublet_cols:
        if col in adata.obs.columns:
            col_data = adata.obs[col].dropna()

            if col_data.dtype == "bool" or set(col_data.unique()).issubset(
                {0, 1, True, False}
            ):
                positive_count = col_data.sum()
                global_stats[f"{col}_count"] = positive_count
                global_stats[f"{col}_percentage"] = positive_count / adata.n_obs * 100
            else:
                global_stats[f"{col}_mean"] = col_data.mean()
                global_stats[f"{col}_median"] = col_data.median()
                global_stats[f"{col}_std"] = col_data.std()
                global_stats[f"{col}_q25"] = col_data.quantile(0.25)
                global_stats[f"{col}_q75"] = col_data.quantile(0.75)

    global_df = pd.DataFrame([global_stats])

    # Export files if directory specified
    if save_dir:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)

        if export_csv:
            sample_file = Path(save_dir) / "doublet_stats_per_sample.csv"
            global_file = Path(save_dir) / "doublet_stats_global.csv"

            sample_df.to_csv(sample_file, index=False)
            global_df.to_csv(global_file, index=False)

            log.info(f"Exported CSV files to {save_dir}")

        if export_xlsx:
            excel_file = Path(save_dir) / "doublet_stats.xlsx"
            with pd.ExcelWriter(excel_file) as writer:
                sample_df.to_excel(writer, sheet_name="per_sample", index=False)
                global_df.to_excel(writer, sheet_name="global", index=False)

            log.info(f"Exported Excel file to {excel_file}")

    return {"sample": sample_df, "global": global_df}


def predict_doublets(
    adata: AnnData,
    config: Optional[DoubletConfig] = None,
    sample_key: str = "sampleID",
    # --- Allow overriding config with specific params ---
    **kwargs,
) -> AnnData:
    """
    Enhanced doublet prediction with a clear, config-driven workflow.

    This function serves as the main entry point for doublet detection, combining
    algorithmic and heuristic approaches for robust results. For basic use, you
    can rely on default parameters. For advanced control and reproducibility,
    it is highly recommended to create and pass a `DoubletConfig` object.

    Args:
        adata: AnnData object containing single-cell expression data
        config: A `DoubletConfig` object. If provided, it overrides all other parameters
        sample_key: Key for sample identification in adata.obs (used if no config)
        use_heuristics: Whether to use marker-based heuristic detection (used if no config)
        marker_species: Species for default markers (used if no config and no custom markers)
        **kwargs: Additional parameters to override defaults or pass to algorithms
                  (e.g., n_pcs=50, merge_strategy='intersection', rate=0.1)

    Returns:
        AnnData object with doublet predictions added to .obs:
        - {method}_score: Algorithmic doublet scores
        - {method}_predicted: Algorithmic doublet predictions
        - heuristic_predicted: Heuristic doublet predictions (if use_heuristics=True)
        - predicted_doublet: Final merged doublet predictions

    Example:
        # Basic usage
        adata = predict_doublets(adata, marker_species="human")

        # Advanced usage with config
        config = DoubletConfig(
            marker_species="human",
            marker_tissue="Lung",
            merge_strategy="intersection",
            save_dir="./results"
        )
        adata = predict_doublets(adata, config=config)
    """
    # === 1. CONFIGURATION SETUP ===
    # Start with a default config
    base_config = DoubletConfig()

    # If a config object is provided, update the base with its values
    if config is not None:
        log.info("Updating defaults with provided DoubletConfig object.")
        base_config.__dict__.update(config.__dict__)

    # Override with any specific kwargs provided
    if kwargs:
        log.info(f"Overriding config with kwargs: {list(kwargs.keys())}")
        for key, value in kwargs.items():
            if hasattr(base_config, key):
                setattr(base_config, key, value)
            else:
                log.warning(f"Unknown parameter '{key}' ignored.")

    # The final, effective config
    cfg = base_config
    try:
        cfg.validate()  # Assuming this method exists in your DoubletConfig
    except AttributeError:
        log.warning("DoubletConfig has no validate() method. Skipping validation.")

    # Validate input data
    if sample_key not in adata.obs.columns:
        raise ValueError(f"Sample key '{sample_key}' not found in adata.obs")

    samples = adata.obs[sample_key].unique()
    if len(samples) == 0:
        raise ValueError(f"No samples found for key '{sample_key}'")

    log.info(
        f"Starting doublet prediction for {adata.n_obs} cells across {len(samples)} samples"
    )
    log.info(
        f"Configuration: method={cfg.method}, merge_strategy={cfg.merge_strategy}, "  # --- FIX: Use cfg consistently
        f"use_heuristics={cfg.use_heuristics}"
    )

    # Initialize result columns
    algo_score_col = f"{cfg.method}_score"
    algo_pred_col = f"{cfg.method}_predicted"
    adata.obs[algo_score_col] = np.nan
    adata.obs[algo_pred_col] = False
    adata.obs[HEURISTIC_PRED_COL] = False

    # Use a dispatcher for multi-algorithm support ---
    ALGORITHM_DISPATCHER = {
        "scrublet": _run_scrublet
        # "doubletfinder": _run_doubletfinder # Future-ready
    }
    if cfg.method not in ALGORITHM_DISPATCHER:
        raise ValueError(
            f"Method '{cfg.method}' is not supported. Available: {list(ALGORITHM_DISPATCHER.keys())}"
        )

    # === 2. ALGORITHMIC DETECTION (Per-Sample) ===
    log.info(f"Running {cfg.method} doublet detection...")

    for sample in samples:
        log.info(f"Processing sample '{sample}' with {cfg.method}...")
        sample_mask = adata.obs[sample_key] == sample
        data_view = adata[sample_mask]

        if data_view.n_obs < 10:
            log.warning(f"Skipping {sample}: fewer than 10 cells.")
            continue

        scores, predicted = ALGORITHM_DISPATCHER[cfg.method](data_view, sample, cfg)

        if scores is not None and predicted is not None:
            adata.obs.loc[sample_mask, algo_score_col] = scores
            adata.obs.loc[sample_mask, algo_pred_col] = predicted

    # === 3. HEURISTIC DETECTION (Global) ===
    if cfg.use_heuristics:
        log.info("Starting heuristic-based doublet detection...")

        # Load marker configurations if not provided
        if cfg.marker_configs is None:
            log.info("Loading marker configurations...")
            try:
                if cfg.marker_tissue or cfg.marker_states:
                    # Use comprehensive marker manager for enhanced detection
                    cfg.marker_configs = create_doublet_marker_config_from_manager(
                        species=cfg.marker_species,
                        tissue=cfg.marker_tissue,
                        states=cfg.marker_states,
                        level=cfg.marker_level,
                        expression_threshold=0.5,
                        min_genes_required=1,
                        use_raw=True,
                        min_markers_per_type=1,
                    )
                    log.info(
                        f"Loaded enhanced markers for {cfg.marker_species} "
                        f"(tissue: {cfg.marker_tissue}, states: {cfg.marker_states})"
                    )
                else:
                    # Fallback to simple marker loading
                    cfg.marker_configs = _load_marker_dict(
                        species=cfg.marker_species, marker_type="major_lineages"
                    )
                    log.info(f"Loaded basic markers for {cfg.marker_species}")

            except Exception as e:
                log.error(f"Failed to load markers: {e}")
                cfg.marker_configs = {}

        # Run heuristic detection if markers are available
        if cfg.marker_configs:
            adata.obs["heuristic_predicted"] = _identify_coexpression_doublets(
                adata,
                marker_configs=cfg.marker_configs,
                min_lineages_for_doublet=cfg.min_lineages_for_doublet,
            )
        else:
            log.warning(
                "Heuristics enabled, but no valid marker configurations were loaded."
            )
            adata.obs["heuristic_predicted"] = False
    else:
        adata.obs["heuristic_predicted"] = False

    # === 4. MERGE RESULTS ===
    log.info("Merging algorithmic and heuristic predictions...")
    adata.obs[FINAL_PRED_COL] = _merge_doublet_predictions(
        adata,
        algorithm_col=algo_pred_col,
        heuristic_col=HEURISTIC_PRED_COL,
        strategy=cfg.merge_strategy,
    )

    # STORE PARAMS
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {})["doublet_params"] = (
        cfg.__dict__
    )

    # === 5. SUMMARY STATISTICS ===
    log.info("\n" + "=" * 50)
    log.info("DOUBLET DETECTION SUMMARY")
    log.info("=" * 50)

    total_cells = adata.n_obs

    # Algorithm results
    algo_count = adata.obs[algo_pred_col].sum()
    log.info(
        f"Algorithm ({cfg.method}): {algo_count} doublets ({algo_count / total_cells:.2%})"
    )

    # Heuristic results
    if cfg.use_heuristics:
        heur_count = adata.obs[HEURISTIC_PRED_COL].sum()
        log.info(f"Heuristic: {heur_count} doublets ({heur_count / total_cells:.2%})")

        # Overlap analysis
        overlap_count = (adata.obs[algo_pred_col] & adata.obs[HEURISTIC_PRED_COL]).sum()
        log.info(
            f"Overlap: {overlap_count} doublets ({overlap_count / total_cells:.2%})"
        )

    # Final merged results
    final_count = adata.obs[FINAL_PRED_COL].sum()
    log.info(f"Final merged: {final_count} doublets ({final_count / total_cells:.2%})")

    # Per-sample breakdown
    log.info("\nPer-sample statistics:")
    for sample in samples:
        sample_mask = adata.obs[sample_key] == sample
        sample_total = sample_mask.sum()
        sample_doublets = adata.obs[FINAL_PRED_COL][sample_mask].sum()
        sample_rate = sample_doublets / sample_total
        log.info(
            f"  {sample}: {sample_doublets}/{sample_total} doublets ({sample_rate:.2%})"
        )

    log.info("=" * 50)

    # === 6. EXPORT STATISTICS ===
    if cfg.export_stats and cfg.save_dir:
        export_doublet_stats(adata, sample_key, save_dir=cfg.save_dir)

    log.info("Doublet prediction workflow completed successfully.")

    return adata
