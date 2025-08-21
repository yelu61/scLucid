
import logging
import os
from typing import Dict, List, Optional

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

# Configure logging
log = logging.getLogger(__name__)

# Export public functions
__all__ = [
    "run_enrichment",
]

def run_enrichment(
    adata: AnnData,
    groupby: str,
    de_key: str = "rank_genes_groups",
    organism: str = "Human",
    gene_sets: List[str] = ["GO_Biological_Process_2023"],
    n_top_genes: int = 100,
    key_added: str = "enrichment",
    min_genes: int = 10,
    max_genes: int = 500,
    min_enrichment_score: float = 0.0,
    max_padj: float = 0.05,
    background_genes: Optional[List[str]] = None,
    plot: bool = False,
    save_path: Optional[str] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Perform functional enrichment analysis for marker genes in each cluster.

    This function runs pathway/gene set enrichment analysis on marker genes
    to identify biological processes associated with each cell group.

    Args:
        adata: AnnData object with differential expression results
        groupby: Column name in adata.obs that specifies cluster grouping
        rank_genes_key: The key_added value used when running find_markers
        organism: Species, either 'Human' or 'Mouse'
        gene_sets: List of gene sets to use for enrichment analysis (from GSEApy)
        n_top_genes: Number of top marker genes per cluster to use
        key_added: Key to store results in adata.uns
        min_genes: Minimum number of genes required for enrichment analysis
        max_genes: Maximum number of genes to include in enrichment analysis
        min_enrichment_score: Minimum enrichment score to include in results
        max_padj: Maximum adjusted p-value for enrichment results
        background_genes: Optional list of background genes (default: all genes in adata)
        plot: Whether to create enrichment plots
        save_path: Directory to save enrichment plots (if plot=True)

    Returns:
        Dictionary mapping cluster names to DataFrames of enrichment results

    Examples:
        >>> # Basic enrichment analysis with default parameters
        >>> enrichment_results = run_enrichment(adata, groupby='leiden')
        >>>
        >>> # More specific analysis with custom parameters
        >>> enrichment_results = run_enrichment(
        ...     adata,
        ...     groupby='cell_types',
        ...     organism='Mouse',
        ...     gene_sets=['KEGG_2019_Mouse', 'WikiPathways_2019_Mouse'],
        ...     n_top_genes=200,
        ...     plot=True,
        ...     save_path='./enrichment_plots/'
        ... )
    """
    log.info(f"Running functional enrichment analysis for {groupby} groups")

    # Check if gene_sets is a string and convert to list
    if isinstance(gene_sets, str):
        gene_sets = [gene_sets]
        log.info(f"Converted gene_sets to list: {gene_sets}")

    # Check if find_markers results exist
    de_results_key = f"{de_key}_df"
    de_path = ['scrnatk', 'analysis', 'de', de_results_key]
    
    if not (de_path[0] in adata.uns and 
            de_path[1] in adata.uns[de_path[0]] and 
            de_path[2] in adata.uns[de_path[0]][de_path[1]] and
            de_path[3] in adata.uns[de_path[0]][de_path[1]][de_path[2]]):
        raise KeyError(f"DE results not found at adata.uns['scrnatk']['analysis']['de']['{de_results_key}']. Run `find_markers` first.")
    
    marker_df = adata.uns['scrnatk']['analysis']['de'][de_results_key]

    # Get unique groups
    if "group" not in marker_df.columns:
        log.error("Column 'group' not found in marker DataFrame")
        raise ValueError("Column 'group' not found in marker DataFrame")

    clusters = marker_df["group"].unique()
    log.info(f"Analyzing {len(clusters)} groups using {', '.join(gene_sets)} gene sets")

    # Set up background genes if not provided
    if background_genes is None:
        background_genes = list(adata.var_names)
        log.info(f"Using all {len(background_genes)} genes in dataset as background")

    # Create save directory if needed
    if plot and save_path is not None:
        os.makedirs(save_path, exist_ok=True)
        log.info(f"Created directory for saving plots: {save_path}")

    # Run enrichment for each cluster
    enrichment_results = {}
    groups_analyzed = 0
    groups_skipped = 0

    for cluster in clusters:
        try:
            # Extract genes from the marker DataFrame
            cluster_df = marker_df[marker_df["group"] == cluster]

            # Apply sorting and filtering
            gene_list = (
                cluster_df.sort_values("logfoldchanges", ascending=False)["names"]
                .head(n_top_genes)
                .tolist()
            )

            # Check if we have enough genes
            if len(gene_list) < min_genes:
                log.warning(
                    f"Skipping group '{cluster}': only {len(gene_list)} marker genes (min: {min_genes})"
                )
                groups_skipped += 1
                continue

            # Limit to maximum number of genes
            if len(gene_list) > max_genes:
                log.info(
                    f"Limiting group '{cluster}' to {max_genes} genes (from {len(gene_list)})"
                )
                gene_list = gene_list[:max_genes]

            log.info(
                f"Running enrichment for group '{cluster}' with {len(gene_list)} genes"
            )

            # Run enrichment analysis
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=gene_sets,
                organism=organism,
                background=background_genes,
                outdir=None,  # Don't generate output files
                cutoff=max_padj,  # Apply p-value cutoff
            )

            # Extract and filter results
            results = enr.results

            if results.empty:
                log.warning(f"No enrichment results found for group '{cluster}'")
                enrichment_results[cluster] = pd.DataFrame()
                groups_skipped += 1
                continue

            # Filter by minimum enrichment score if requested
            if min_enrichment_score > 0:
                results = results[results["Combined Score"] >= min_enrichment_score]

            # Store filtered results
            enrichment_results[cluster] = results

            # Create plots if requested
            if plot:
                try:
                    # Get top pathways
                    top_pathways = results.head(20)

                    if not top_pathways.empty:
                        # Create plot of top enriched terms
                        plt.figure(
                            figsize=(12, min(10, max(4, len(top_pathways) * 0.3)))
                        )

                        # Plot negative log p-value as bar chart
                        plt.barh(
                            top_pathways["Term"].str.split(" ").str[:5].str.join(" "),
                            -np.log10(top_pathways["Adjusted P-value"]),
                            color="skyblue",
                        )

                        plt.xlabel("-log10(Adjusted P-value)")
                        plt.title(f"Top Enriched Pathways for {cluster}")
                        plt.tight_layout()

                        # Save plot if path provided
                        if save_path is not None:
                            safe_cluster = cluster.replace("/", "_").replace(" ", "_")
                            plt.savefig(
                                f"{save_path}/{safe_cluster}_enrichment.png", dpi=300
                            )
                            plt.close()
                        else:
                            plt.show()

                except Exception as e:
                    log.warning(
                        f"Error creating enrichment plot for '{cluster}': {str(e)}"
                    )

            groups_analyzed += 1

        except Exception as e:
            log.error(f"Error analyzing group '{cluster}': {str(e)}")
            enrichment_results[cluster] = pd.DataFrame()
            groups_skipped += 1

    # Store all results in the AnnData object
    adata.uns.setdefault('scrnatk', {}).setdefault('analysis', {}).setdefault('de', {})
    adata.uns['scrnatk']['analysis']['de'][key_added] = {
        'results': enrichment_results,
        'params': {
            'groupby': groupby,
            'de_key': de_key,
            # ... other params ...
        }
    }

    log.info(
        f"Enrichment analysis complete: {groups_analyzed} groups analyzed, {groups_skipped} groups skipped"
    )
    log.info(f"Results stored in adata.uns['{key_added}']")

    return enrichment_results
