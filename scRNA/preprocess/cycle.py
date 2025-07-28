"""
Cell cycle scoring and utility functions for single-cell RNA-seq data.

This module provides functions for scoring cell cycle phases based on
species-specific or user-provided gene lists. It supports human, mouse, and rat
gene symbols and includes automatic species detection capabilities.
"""

import json
import logging
import os
from typing import Dict, List, Literal, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import pkg_resources
import scanpy as sc
from anndata import AnnData

# Configure logging
log = logging.getLogger(__name__)

# Export public functions
__all__ = ["score_cell_cycle", "get_cell_cycle_genes"]


# --- Helper Function to Load Genes ---
def _load_cell_cycle_genes() -> Dict[str, Dict[str, List[str]]]:
    """
    Loads cell cycle genes from the package's resource file.

    Returns:
        Dictionary with species as keys and S/G2M gene lists as values.

    Note:
        The expected format of the JSON file is:
        {
            "human": {
                "s_genes": ["MCM5", "PCNA", ...],
                "g2m_genes": ["AURKA", "BUB1", ...]
            },
            "mouse": { ... },
            "rat": { ... }
        }
    """
    try:
        gene_path = pkg_resources.resource_filename(
            "scRNA", "resources/cell_cycle_genes.json"
        )
        log.debug(f"Loading cell cycle genes from: {gene_path}")

        with open(gene_path, "r") as f:
            genes_dict = json.load(f)

        # Validate the structure of the loaded data
        required_species = ["human", "mouse", "rat"]
        required_keys = ["s_genes", "g2m_genes"]

        for species in required_species:
            if species not in genes_dict:
                log.warning(f"Species '{species}' not found in cell cycle genes file")
                continue

            for key in required_keys:
                if key not in genes_dict[species]:
                    log.warning(f"Key '{key}' not found for species '{species}'")
                elif not isinstance(genes_dict[species][key], list):
                    log.warning(f"'{key}' for '{species}' is not a list")
                elif len(genes_dict[species][key]) == 0:
                    log.warning(f"'{key}' for '{species}' is empty")

        log.info(f"Successfully loaded cell cycle genes for {len(genes_dict)} species")
        return genes_dict

    except FileNotFoundError:
        log.error(
            "cell_cycle_genes.json not found. Please ensure the package is installed correctly."
        )
        # Provide a minimal fallback for human genes to avoid completely failing
        log.warning("Using minimal fallback gene lists for human only")
        return {
            "human": {
                "s_genes": [
                    "MCM5",
                    "PCNA",
                    "TYMS",
                    "FEN1",
                    "MCM2",
                    "MCM4",
                    "RRM1",
                    "UNG",
                    "GINS2",
                    "MCM6",
                    "CDCA7",
                ],
                "g2m_genes": [
                    "HMGB2",
                    "CDK1",
                    "NUSAP1",
                    "UBE2C",
                    "BIRC5",
                    "TPX2",
                    "TOP2A",
                    "NDC80",
                    "CKS2",
                    "NUF2",
                ],
            }
        }
    except json.JSONDecodeError as e:
        log.error(f"Invalid JSON in cell cycle genes file: {e}")
        return {}
    except Exception as e:
        log.error(f"An error occurred while loading cell cycle genes: {e}")
        return {}


# Load genes once when the module is imported
SPECIES_GENES = _load_cell_cycle_genes()


def get_cell_cycle_genes(
    species: Literal["human", "mouse", "rat"] = "human",
) -> Dict[str, List[str]]:
    """
    Returns the list of cell cycle marker genes for a specific species.

    This function provides access to the built-in cell cycle gene lists
    that are used for cell cycle scoring.

    Args:
        species: The species for which to return cell cycle genes.
                 One of "human", "mouse", or "rat".

    Returns:
        Dictionary with keys "s_genes" and "g2m_genes", each containing
        a list of gene symbols.

    Raises:
        ValueError: If the requested species is not available.

    Examples:
        >>> genes = get_cell_cycle_genes("mouse")
        >>> print(f"S-phase genes: {len(genes['s_genes'])}")
        >>> print(f"G2M-phase genes: {len(genes['g2m_genes'])}")
    """
    if not SPECIES_GENES:
        raise RuntimeError("Could not load cell cycle gene lists. Cannot proceed.")

    if species not in SPECIES_GENES:
        available = ", ".join(SPECIES_GENES.keys())
        raise ValueError(
            f"Unknown species: '{species}'. Valid options are: {available}"
        )

    return {
        "s_genes": SPECIES_GENES[species]["s_genes"],
        "g2m_genes": SPECIES_GENES[species]["g2m_genes"],
    }


def _detect_species(adata: AnnData) -> str:
    """
    Attempts to automatically detect the species based on gene presence.

    Args:
        adata: AnnData object containing gene expression data.

    Returns:
        Detected species name ("human", "mouse", "rat") or "unknown".
    """
    # Check for each species
    detection_scores = {}

    for species in SPECIES_GENES:
        # Count how many of the S and G2M genes are present
        s_genes = SPECIES_GENES[species]["s_genes"]
        g2m_genes = SPECIES_GENES[species]["g2m_genes"]

        s_found = sum(1 for gene in s_genes[:20] if gene in adata.var_names)
        g2m_found = sum(1 for gene in g2m_genes[:20] if gene in adata.var_names)

        # Calculate a detection score (percentage of marker genes found)
        total_checked = min(len(s_genes), 20) + min(len(g2m_genes), 20)
        score = (s_found + g2m_found) / total_checked if total_checked > 0 else 0

        detection_scores[species] = score
        log.debug(
            f"Species detection - {species}: score={score:.2f} "
            f"(S: {s_found}/{min(len(s_genes), 20)}, "
            f"G2M: {g2m_found}/{min(len(g2m_genes), 20)})"
        )

    # Find the species with the highest score
    best_species = max(detection_scores, key=detection_scores.get)
    best_score = detection_scores[best_species]

    # Only return a species if the score is reasonably high
    if best_score >= 0.3:  # At least 30% of genes found
        log.info(f"Detected species: {best_species} (confidence: {best_score:.2f})")
        return best_species
    else:
        log.warning(
            f"Could not confidently detect species. Best guess: {best_species} with low confidence ({best_score:.2f})"
        )
        return "unknown"


def _validate_gene_lists(
    adata: AnnData, s_genes: List[str], g2m_genes: List[str]
) -> Tuple[List[str], List[str], bool]:
    """
    Validates and filters gene lists to those present in the data.

    Args:
        adata: AnnData object containing gene expression data.
        s_genes: List of S phase marker genes.
        g2m_genes: List of G2M phase marker genes.

    Returns:
        Tuple containing:
        - Filtered S phase genes
        - Filtered G2M phase genes
        - Boolean indicating whether enough genes were found
    """
    # Filter gene lists to only those present in the data
    s_genes_found = [gene for gene in s_genes if gene in adata.var_names]
    g2m_genes_found = [gene for gene in g2m_genes if gene in adata.var_names]

    # Log how many genes were found
    log.info(f"Found {len(s_genes_found)}/{len(s_genes)} S-phase genes in the dataset")
    log.info(
        f"Found {len(g2m_genes_found)}/{len(g2m_genes)} G2M-phase genes in the dataset"
    )

    # Check if enough genes were found
    has_enough_genes = True
    if len(s_genes_found) < 5 or len(g2m_genes_found) < 5:
        log.warning(
            f"Few cell cycle genes found in data (S: {len(s_genes_found)}, G2M: {len(g2m_genes_found)})."
        )
        log.warning(
            "Scoring may be inaccurate. This could be due to incorrect species, gene ID format, or data filtering."
        )

        if len(s_genes_found) < 3 or len(g2m_genes_found) < 3:
            log.error("Insufficient cell cycle genes found to proceed with scoring.")
            has_enough_genes = False

    return s_genes_found, g2m_genes_found, has_enough_genes


def _plot_cell_cycle(
    adata: AnnData, species: str, save_dir: Optional[str] = None
) -> plt.Figure:
    """
    Generates plots visualizing cell cycle scores and phase distribution.

    Args:
        adata: AnnData object with cell cycle scores.
        species: Species name for the title.
        save_dir: Directory to save the plot.

    Returns:
        The matplotlib Figure object.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), facecolor="white")
    fig.suptitle(
        f"Cell Cycle Analysis ({species.capitalize()})", fontsize=16, fontweight="bold"
    )

    # Scatter plot of scores
    sc.pl.scatter(
        adata,
        x="S_score",
        y="G2M_score",
        color="phase",
        title="Cell Cycle Scores",
        ax=axes[0],
        show=False,
        palette={"G1": "#1f77b4", "S": "#ff7f0e", "G2M": "#2ca02c"},
    )

    # Add decision boundaries as dashed lines
    axes[0].axhline(y=0, color="gray", linestyle="--", alpha=0.5)
    axes[0].axvline(x=0, color="gray", linestyle="--", alpha=0.5)

    # Bar plot of phase distribution
    phase_counts = adata.obs["phase"].value_counts().sort_index()
    colors = {"G1": "#1f77b4", "S": "#ff7f0e", "G2M": "#2ca02c"}
    bar_colors = [colors.get(phase, "#999999") for phase in phase_counts.index]

    axes[1].bar(phase_counts.index, phase_counts.values, color=bar_colors)

    # Add counts as text on bars
    for i, (phase, count) in enumerate(phase_counts.items()):
        pct = count / adata.n_obs * 100
        axes[1].text(
            i,
            count,
            f"{count}\n({pct:.1f}%)",
            ha="center",
            va="bottom",
            fontweight="bold",
        )

    axes[1].set_title("Cell Phase Distribution")
    axes[1].set_ylabel("Number of Cells")
    axes[1].set_ylim(0, phase_counts.max() * 1.15)  # More space for labels
    axes[1].tick_params(axis="x", rotation=0)

    # Add a text box with summary statistics
    phase_pcts = phase_counts / phase_counts.sum() * 100
    stats_text = "\n".join(
        [
            f"Total cells: {adata.n_obs}",
            f"G1: {phase_pcts.get('G1', 0):.1f}%",
            f"S: {phase_pcts.get('S', 0):.1f}%",
            f"G2M: {phase_pcts.get('G2M', 0):.1f}%",
        ]
    )

    # Place text box in the upper right corner
    props = dict(boxstyle="round", facecolor="white", alpha=0.7)
    axes[1].text(
        0.95,
        0.95,
        stats_text,
        transform=axes[1].transAxes,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=props,
        fontsize=10,
    )

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, f"cell_cycle_scores_{species}.png")
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        log.info(f"Saved cell cycle plot to {save_path}")

    return fig


def score_cell_cycle(
    adata: AnnData,
    species: Literal["human", "mouse", "rat"] = "human",
    s_genes: Optional[List[str]] = None,
    g2m_genes: Optional[List[str]] = None,
    layer: Optional[str] = None,
    copy: bool = False,
    plot: bool = True,
    save_dir: Optional[str] = None,
    force: bool = False,
) -> AnnData:
    """
    Scores cell cycle phases (S and G2M) based on canonical marker genes.

    This function uses species-specific gene lists which can be overridden by
    user-provided lists. It also attempts to auto-detect the correct species
    if the default (human) genes are not found.

    Note: This function modifies the AnnData object in place by adding
    'S_score', 'G2M_score', and 'phase' columns to `adata.obs`.

    Args:
        adata: AnnData object.
        species: Species of the dataset. Used to select default gene lists.
        s_genes: Custom list of S phase marker genes. Overrides species default.
        g2m_genes: Custom list of G2M phase marker genes. Overrides species default.
        layer: If provided, use this layer for calculating scores instead of .X
        copy: Whether to return a copy of the AnnData object.
        plot: Whether to generate plots of cell cycle scores and phase distribution.
        save_dir: Directory to save plots. If None, plots are not saved to disk.
        force: Whether to recompute scores if they already exist.

    Returns:
        The modified AnnData object with cell cycle scores.

    Raises:
        ValueError: If the species is unknown or insufficient genes are found.
        RuntimeError: If the cell cycle gene lists couldn't be loaded.

    Examples:
        >>> # Score cell cycle using default human gene lists
        >>> adata = score_cell_cycle(adata)
        >>>
        >>> # Score cell cycle for mouse data
        >>> adata = score_cell_cycle(adata, species="mouse")
        >>>
        >>> # Use custom gene lists
        >>> s_genes = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2"]
        >>> g2m_genes = ["AURKA", "BUB1", "TOP2A", "PLK1"]
        >>> adata = score_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    """
    # Handle copy logic
    if copy:
        adata = adata.copy()

    # Check if scores already exist
    if (
        "S_score" in adata.obs
        and "G2M_score" in adata.obs
        and "phase" in adata.obs
        and not force
    ):
        log.info("Cell cycle scores already exist. Use force=True to recompute.")
        if plot:
            # Still generate plots if requested
            _plot_cell_cycle(adata, species, save_dir)
        return adata

    # --- Gene List Selection and Validation ---
    log.info(f"Scoring cell cycle phases using {species} gene lists")

    if not SPECIES_GENES:
        log.error("Could not load cell cycle gene lists.")
        raise RuntimeError("Failed to load cell cycle gene lists. Cannot proceed.")

    # Determine which species and gene lists to use
    if s_genes is None or g2m_genes is None:
        # If species is unknown or invalid, attempt to auto-detect
        if species not in SPECIES_GENES:
            log.warning(f"Unknown species: '{species}'. Attempting auto-detection.")
            detected_species = _detect_species(adata)
            if detected_species == "unknown":
                available = ", ".join(SPECIES_GENES.keys())
                raise ValueError(
                    f"Unknown species: '{species}' and auto-detection failed. Valid options: {available}"
                )
            current_species = detected_species
        else:
            current_species = species

            # Auto-detect species if default genes aren't found
            if species == "human" and not any(
                g in adata.var_names for g in SPECIES_GENES["human"]["s_genes"][:10]
            ):
                log.info("Testing if mouse or rat genes better match the dataset...")
                for sp in ["mouse", "rat"]:
                    if any(
                        g in adata.var_names for g in SPECIES_GENES[sp]["s_genes"][:10]
                    ):
                        log.warning(
                            f"Default species is 'human', but '{sp}' genes were detected. "
                            f"Switching gene list to '{sp}'."
                        )
                        current_species = sp
                        break

        # Use the appropriate gene lists
        s_genes_list = (
            s_genes
            if s_genes is not None
            else SPECIES_GENES[current_species]["s_genes"]
        )
        g2m_genes_list = (
            g2m_genes
            if g2m_genes is not None
            else SPECIES_GENES[current_species]["g2m_genes"]
        )
    else:
        log.info("Using user-provided S and G2M gene lists.")
        current_species = "custom"
        s_genes_list = s_genes
        g2m_genes_list = g2m_genes

    # Validate and filter gene lists
    s_genes_found, g2m_genes_found, has_enough_genes = _validate_gene_lists(
        adata, s_genes_list, g2m_genes_list
    )

    if not has_enough_genes:
        raise ValueError(
            "Insufficient cell cycle genes found to proceed with scoring. "
            "Provide correct gene lists or check gene naming conventions."
        )

    log.info(
        f"Scoring cell cycle using {len(s_genes_found)} S-phase and {len(g2m_genes_found)} G2M-phase genes."
    )

    # --- Scoring ---
    try:
        sc.tl.score_genes_cell_cycle(
            adata, s_genes=s_genes_found, g2m_genes=g2m_genes_found, layer=layer
        )
        log.info("Cell cycle scoring completed successfully.")
    except Exception as e:
        log.error(f"Cell cycle scoring failed: {str(e)}")
        raise RuntimeError(f"Failed to compute cell cycle scores: {str(e)}")

    # Calculate additional metrics
    adata.obs["cc_diff"] = adata.obs["S_score"] - adata.obs["G2M_score"]

    # Store metadata
    adata.uns["cell_cycle"] = {
        "species_used": current_species,
        "s_genes_used": s_genes_found,
        "g2m_genes_used": g2m_genes_found,
        "date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
    }

    # Calculate phase percentages
    phase_counts = adata.obs["phase"].value_counts()
    phase_pcts = phase_counts / phase_counts.sum() * 100

    log.info("Cell cycle phase distribution:")
    for phase in sorted(phase_counts.index):
        log.info(f"  - {phase}: {phase_counts[phase]} cells ({phase_pcts[phase]:.1f}%)")

    # --- Plotting ---
    if plot:
        try:
            fig = _plot_cell_cycle(adata, current_species, save_dir)
        except Exception as e:
            log.warning(f"Failed to generate cell cycle plots: {str(e)}")

    return adata
