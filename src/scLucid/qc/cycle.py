"""
Cell cycle scoring and utility functions for single-cell RNA-seq data.

This module provides functions for scoring cell cycle phases based on
species-specific or user-provided gene lists. It supports human, mouse, and rat
gene symbols and includes automatic species detection capabilities.
"""

import json
import logging
from importlib import resources
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from anndata import AnnData

# Configure logging
log = logging.getLogger(__name__)

# Export public functions
__all__ = ["score_cell_cycle"]


# --- Helper Functions ---
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
        gene_path = resources.files("scLucid").joinpath("resources/cell_cycle_genes.json")
        log.debug(f"Loading cell cycle genes from: {gene_path}")

        with open(gene_path) as f:
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


def _get_available_species() -> List[str]:
    """Return species keys that have both S and G2M gene lists."""
    available = []
    for species, genes in SPECIES_GENES.items():
        if isinstance(genes, dict) and "s_genes" in genes and "g2m_genes" in genes:
            available.append(species)
    return available


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

    for species in _get_available_species():
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

    # --- Check for ambiguity between top two species ---
    if not detection_scores:
        log.warning("No valid species definitions found in cell cycle gene resource.")
        return "unknown"

    sorted_scores = sorted(detection_scores.items(), key=lambda item: item[1], reverse=True)
    best_species, best_score = sorted_scores[0]

    if len(sorted_scores) > 1:
        second_best_species, second_best_score = sorted_scores[1]
        # If the top score is high but the second is very close, warn the user.
        if best_score > 0.3 and (best_score - second_best_score) < 0.1:
            log.warning(
                f"Ambiguous species detection: '{best_species}' (score: {best_score:.2f}) and "
                f"'{second_best_species}' (score: {second_best_score:.2f}) are very close. "
                "Please verify the species manually."
            )

    if best_score >= 0.3:
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
    log.info(f"Found {len(g2m_genes_found)}/{len(g2m_genes)} G2M-phase genes in the dataset")

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


def _plot_cell_cycle(adata: AnnData, species: str, save_dir: Optional[str] = None) -> plt.Figure:
    """
    Generates plots visualizing cell cycle scores and phase distribution.

    Args:
        adata: AnnData object with cell cycle scores.
        species: Species name for the title.
        save_dir: Directory to save the plot.

    Returns:
        The matplotlib Figure object.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), facecolor="white")
    fig.suptitle(f"Cell Cycle Analysis ({species.capitalize()})", fontsize=16, fontweight="bold")

    # Scatter plot of scores
    sc.pl.scatter(
        adata,
        x="S_score",
        y="G2M_score",
        color="phase",
        title="Cell Cycle Scores",
        ax=axes[0],
        show=False,
    )

    # Add decision boundaries as dashed lines
    axes[0].axhline(y=0, color="gray", linestyle="--", alpha=0.5)
    axes[0].axvline(x=0, color="gray", linestyle="--", alpha=0.5)

    # Bar plot of phase distribution
    phase_counts = adata.obs["phase"].value_counts().sort_index()
    axes[1].bar(phase_counts.index, phase_counts.values, color=["#1f77b4", "#ff7f0e", "#2ca02c"])

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
        Path(save_dir).mkdir(parents=True, exist_ok=True)
        save_path = Path(save_dir) / f"cell_cycle_scores_{species}.png"
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        log.info(f"Saved cell cycle plot to {save_path}")

    return fig


def _get_validated_genes(
    adata: AnnData,
    species: str,
    s_genes: Optional[List[str]],
    g2m_genes: Optional[List[str]],
) -> Tuple[List[str], List[str], str]:
    """
    Selects, validates, and filters cell cycle gene lists.

    This helper centralizes the logic for handling user-provided lists,
    species-based defaults, and auto-detection.

    Returns:
        Tuple of (s_genes_found, g2m_genes_found, species_used)
    """
    if s_genes is not None and g2m_genes is not None:
        log.info("Using user-provided S and G2M gene lists.")
        species_used = "custom"
        s_genes_list, g2m_genes_list = s_genes, g2m_genes
    else:
        # If species is not valid, attempt auto-detection.
        available_species = _get_available_species()
        if species not in available_species:
            log.warning(f"Unknown species: '{species}'. Attempting auto-detection.")
            detected_species = _detect_species(adata)
            if detected_species == "unknown":
                available = ", ".join(available_species)
                raise ValueError(
                    f"Unknown species: '{species}' and auto-detection failed. Valid options: {available}"
                )
            species_used = detected_species
        else:
            species_used = species
            # Even if a valid species is provided, check if another might fit better.
            current_genes_found = any(
                g in adata.var_names for g in SPECIES_GENES[species_used]["s_genes"][:10]
            )
            if not current_genes_found:
                log.warning(
                    f"Very few genes for the specified species '{species_used}' were found. Auto-detecting."
                )
                detected_species = _detect_species(adata)
                if detected_species != "unknown":
                    species_used = detected_species

        log.info(f"Using gene list for species: '{species_used}'")
        s_genes_list = SPECIES_GENES[species_used]["s_genes"]
        g2m_genes_list = SPECIES_GENES[species_used]["g2m_genes"]

    # Final validation against the data
    s_genes_found, g2m_genes_found, has_enough_genes = _validate_gene_lists(
        adata, s_genes_list, g2m_genes_list
    )

    if not has_enough_genes:
        raise ValueError("Insufficient cell cycle genes found to proceed with scoring.")

    return s_genes_found, g2m_genes_found, species_used


# --- Main Functions ---
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
    if "phase" in adata.obs and not force:
        log.info("Cell cycle scores already exist. Use force=True to recompute.")
        if plot:
            species_used = (
                adata.uns.get("sclucid", {})
                .get("qc", {})
                .get("cell_cycle", {})
                .get("species_used", species)
            )
            _plot_cell_cycle(adata, species_used, save_dir)
        return adata

    if not SPECIES_GENES:
        raise RuntimeError("Failed to load cell cycle gene lists. Cannot proceed.")

    # --- All gene selection logic is now in the helper function ---
    s_genes_found, g2m_genes_found, species_used = _get_validated_genes(
        adata, species, s_genes, g2m_genes
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

    # Mark cell cycle genes in var
    adata.var["is_cell_cycle"] = False  # 默认False
    adata.var.loc[s_genes_found + g2m_genes_found, "is_cell_cycle"] = True
    log.info(
        f"Marked {len(s_genes_found + g2m_genes_found)} cell cycle genes in adata.var['is_cell_cycle']."
    )

    # Calculate additional metrics
    adata.obs["cc_diff"] = adata.obs["S_score"] - adata.obs["G2M_score"]

    # Store metadata
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {})["cell_cycle"] = {
        "species_used": species_used,
        "s_genes_used_count": len(s_genes_found),
        "g2m_genes_used_count": len(g2m_genes_found),
        "params": {
            "species_requested": species,
            "custom_s_genes_provided": s_genes is not None,
            "custom_g2m_genes_provided": g2m_genes is not None,
            "layer": layer,
        },
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
            fig = _plot_cell_cycle(adata, species_used, save_dir)
        except Exception as e:
            log.warning(f"Failed to generate cell cycle plots: {str(e)}")

    return adata


def score_cell_cycle_advanced(
    adata: AnnData,
    species: str = "human",
    regress_out: bool = False,  # 新功能
    plot_phase_markers: bool = True,  # 新功能
    **kwargs,
) -> AnnData:
    """
    增强版细胞周期打分，可选回归。
    """
    # 现有打分逻辑
    adata = score_cell_cycle(adata, species, **kwargs)

    # 新功能1: 可选的细胞周期效应回归
    if regress_out:
        import scanpy as sc

        log.info("Regressing out cell cycle effects...")
        sc.pp.regress_out(adata, ["S_score", "G2M_score"])
        adata.uns["cell_cycle_regressed"] = True

    # 新功能2: 每个phase的marker基因表达热图
    if plot_phase_markers:
        _plot_phase_specific_markers(adata, species)

    return adata


def _plot_phase_specific_markers(adata: AnnData, species: str, save_path: Optional[str] = None):
    """
    可视化不同cell cycle phase的marker基因表达。
    """
    import scanpy as sc

    s_genes = SPECIES_GENES[species]["s_genes"][:10]
    g2m_genes = SPECIES_GENES[species]["g2m_genes"][:10]

    marker_genes = s_genes + g2m_genes
    marker_genes = [g for g in marker_genes if g in adata.var_names]

    if len(marker_genes) < 5:
        log.warning("Too few marker genes found for visualization")
        return

    # Create dotplot
    sc.pl.dotplot(
        adata, marker_genes, groupby="phase", dendrogram=True, standard_scale="var", save=save_path
    )
