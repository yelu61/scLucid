"""
Cell cycle scoring and utility functions for single-cell RNA-seq data.

This module provides functions for scoring cell cycle phases based on
species-specific or user-provided gene lists.
"""

import matplotlib.pyplot as plt
import scanpy as sc
import json
import pkg_resources
import os
from typing import Optional, List, Literal

# --- Helper Function to Load Genes ---
def _load_cell_cycle_genes():
    """Loads cell cycle genes from the package's resource file."""
    try:
        gene_path = pkg_resources.resource_filename('scRNA', "resources/cell_cycle_genes.json")
        with open(gene_path, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        print("Error: cell_cycle_genes.json not found. Please ensure the package is installed correctly.")
        return {}
    except Exception as e:
        print(f"An error occurred while loading cell cycle genes: {e}")
        return {}

# Load genes once when the module is imported
SPECIES_GENES = _load_cell_cycle_genes()


def score_cell_cycle(
    adata: sc.AnnData,
    species: Literal["human", "mouse", "rat"] = "human",
    s_genes: Optional[List[str]] = None,
    g2m_genes: Optional[List[str]] = None,
    plot: bool = True,
    save_dir: Optional[str] = None,
) -> sc.AnnData:
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
        plot: Whether to generate plots of cell cycle scores and phase distribution.
        save_dir: Directory to save plots. If None, plots are not saved to disk.

    Returns:
        The modified AnnData object with cell cycle scores.
    """
    # --- Gene List Selection and Validation ---
    if not SPECIES_GENES:
        raise RuntimeError("Could not load cell cycle gene lists. Cannot proceed.")

    if s_genes is None or g2m_genes is None:
        if species not in SPECIES_GENES:
            raise ValueError(f"Unknown species: '{species}'. Valid options are: {', '.join(SPECIES_GENES.keys())}")
        
        # Auto-detect species if default (human) genes are not prominent
        current_species = species
        if species == "human" and not any(g in adata.var_names for g in SPECIES_GENES["human"]["s_genes"][:10]):
            for sp in ["mouse", "rat"]:
                if any(g in adata.var_names for g in SPECIES_GENES[sp]["s_genes"][:10]):
                    print(f"Warning: Default species is 'human', but '{sp}' genes were detected. Switching gene list to '{sp}'.")
                    current_species = sp
                    break
        
        s_genes = s_genes if s_genes is not None else SPECIES_GENES[current_species]["s_genes"]
        g2m_genes = g2m_genes if g2m_genes is not None else SPECIES_GENES[current_species]["g2m_genes"]
    else:
        print("Using user-provided S and G2M gene lists.")
        current_species = "custom"

    # Filter gene lists to only those present in the data
    s_genes_found = [gene for gene in s_genes if gene in adata.var_names]
    g2m_genes_found = [gene for gene in g2m_genes if gene in adata.var_names]

    # Provide detailed warnings if few genes are found
    if len(s_genes_found) < 5 or len(g2m_genes_found) < 5:
        print(f"Warning: Few cell cycle genes found in data (S: {len(s_genes_found)}, G2M: {len(g2m_genes_found)}).")
        print("Scoring may be inaccurate. This could be due to incorrect species, gene ID format, or data filtering.")
        if len(s_genes_found) < 3 or len(g2m_genes_found) < 3:
            raise ValueError("Insufficient cell cycle genes found to proceed with scoring.")

    print(f"Scoring cell cycle using {len(s_genes_found)} S-phase and {len(g2m_genes_found)} G2M-phase genes.")

    # --- Scoring ---
    sc.tl.score_genes_cell_cycle(
        adata,
        s_genes=s_genes_found,
        g2m_genes=g2m_genes_found,
    )

    # Store metadata
    adata.uns['cell_cycle'] = {
        'species_used': current_species,
        's_genes_used': s_genes_found,
        'g2m_genes_used': g2m_genes_found
    }

    # --- Plotting ---
    if plot:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), facecolor='white')
        fig.suptitle(f'Cell Cycle Analysis ({current_species.capitalize()})', fontsize=16)

        # Scatter plot of scores
        sc.pl.scatter(
            adata, x='S_score', y='G2M_score', color='phase',
            title='Cell Cycle Scores', ax=axes[0], show=False
        )
        
        # Bar plot of phase distribution
        phase_counts = adata.obs['phase'].value_counts().sort_index()
        axes[1].bar(phase_counts.index, phase_counts.values, color=['#2ca02c', '#ff7f0e', '#1f77b4'])
        for i, (phase, count) in enumerate(phase_counts.items()):
            axes[1].text(i, count, f' {count}', ha='center', va='bottom')
        axes[1].set_title('Cell Phase Distribution')
        axes[1].set_ylabel('Number of Cells')
        axes[1].set_ylim(0, phase_counts.max() * 1.1)
        axes[1].tick_params(axis='x', rotation=0)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            plt.savefig(os.path.join(save_dir, f"cell_cycle_scores_{current_species}.png"), dpi=300, bbox_inches='tight')
        plt.show()
        plt.close(fig)

    return adata
