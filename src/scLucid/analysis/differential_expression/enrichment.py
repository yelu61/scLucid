"""
Enrichment analysis functions.

This module provides functional enrichment analysis capabilities:
- run_enrichment: ORA and GSEA enrichment analysis
- export_enrichment_results: Export enrichment results in various formats
- batch_celltype_deg_enrichment: Batch DE + enrichment for multiple cell types
"""

import logging
import pickle
from importlib import resources
from pathlib import Path
from typing import Any, Dict, List, Optional

import gseapy as gp
import pandas as pd
from anndata import AnnData

from ...utils.helpers import sanitize_for_hdf5
from ..config import EnrichmentConfig
from .de_utils import _safe_filename
from .scanpy_compat import standardize_enrichment_cols

log = logging.getLogger(__name__)


def run_enrichment(
    adata: AnnData,
    groupby: Optional[str] = None,
    config: Optional[EnrichmentConfig] = None,
    **kwargs,
) -> Dict[str, Dict[str, pd.DataFrame]]:
    """
    Run functional enrichment analysis (ORA and/or GSEA) using GSEApy.

    Supports:
    - Online mode: Direct API calls to Enrichr/GSEA servers
    - Offline mode: Local GMT files (faster, more stable)
    - Multiple gene set categories per analysis
    - Both ORA (Over-Representation) and GSEA (Prerank) methods

    Args:
        adata: AnnData object with DE results
        groupby: Optional grouping column (usually inferred from DE results)
        config: EnrichmentConfig object
        **kwargs: Override config parameters

    Returns:
        Nested dictionary:
        {
            "cluster1": {
                "ora": DataFrame with ORA results,
                "gsea": DataFrame with GSEA results
            },
            "cluster2": {...},
            ...
        }

    Example:
        >>> config = EnrichmentConfig(
        ...     de_key="rank_genes_groups_filtered_df",
        ...     mode="offline",
        ...     method="both",
        ...     organism="human",
        ...     gene_sets_offline=["hallmark", "go_bp", "reactome"]
        ... )
        >>> enr_results = run_enrichment(adata, groupby="leiden", config=config)
        >>>
        >>> # Access cluster 0 ORA results:
        >>> cluster0_ora = enr_results["0"]["ora"]

    Notes:
        - Offline mode requires GMT files in scLucid/resources/
        - Results stored at: adata.uns['sclucid']['analysis']['de']['{key_added}']
        - Background gene list automatically set to all genes in adata
    """
    # Configuration
    if config is None:
        config = EnrichmentConfig(**kwargs)
    else:
        config = config.model_copy(update=kwargs)

    # Load DE results
    root = adata.uns.get("sclucid", {}).get("analysis", {}).get("de", {})
    if config.de_key not in root:
        raise KeyError(
            f"DE results not found at .uns['sclucid']['analysis']['de']['{config.de_key}']. "
            "Run find_markers() or filter_markers() first."
        )

    marker_df = root[config.de_key]
    if marker_df.empty:
        log.warning("Marker DataFrame is empty. Skipping enrichment analysis.")
        return {}

    # Determine groups
    if groupby and "group" in marker_df.columns:
        group_order = list(pd.unique(marker_df["group"]))
    else:
        group_order = ["all"]
        marker_df = marker_df.copy()
        marker_df["group"] = "all"

    if config.verbose:
        log.info(
            f"Running {config.method.upper()} enrichment for {len(group_order)} groups "
            f"in {config.mode} mode"
        )

    # Background genes
    background_genes = list(adata.var_names)

    # Prepare gene sets
    gmt_files_to_run = {}
    gene_sets_list = (
        config.gene_sets_online if config.mode == "online" else config.gene_sets_offline
    )

    if not isinstance(gene_sets_list, list):
        gene_sets_list = [gene_sets_list]

    if config.mode == "offline":
        # Custom GMT file takes priority
        if config.custom_gene_sets and Path(config.custom_gene_sets).is_file():
            gmt_files_to_run = {"custom": config.custom_gene_sets}
            log.info(f"Using custom gene set file: {config.custom_gene_sets}")
        else:
            # Load built-in GMT files
            for gs_category in gene_sets_list:
                try:
                    filename = f"{config.organism.lower()}_{gs_category}_{config.gmt_version}.gmt"
                    file_path = resources.files("scLucid").joinpath("resources", filename)

                    if file_path.is_file():
                        gmt_files_to_run[gs_category] = str(file_path)
                        log.debug(f"Loaded GMT file: {filename}")
                    else:
                        log.warning(f"GMT file not found: {filename}")

                except Exception as e:
                    log.error(f"Error loading GMT file for '{gs_category}': {e}")

            if not gmt_files_to_run:
                raise FileNotFoundError(
                    f"No valid GMT files found for offline mode. " f"Searched for: {gene_sets_list}"
                )
    else:
        # Online mode: gene set names are used directly
        gmt_files_to_run = {gs: gs for gs in gene_sets_list}

    # Select GSEA ranking column
    rank_col = config.rank_col_gsea
    if config.prefer_score_for_enrichment and "scores" in marker_df.columns:
        rank_col = "scores"
        log.debug("Using 'scores' for GSEA ranking (prefer_score_for_enrichment=True)")
    elif rank_col not in marker_df.columns:
        fallback = "scores" if "scores" in marker_df.columns else "logfoldchanges"
        log.warning(f"GSEA rank column '{rank_col}' not found. Falling back to '{fallback}'")
        rank_col = fallback

    if rank_col not in marker_df.columns:
        raise KeyError(f"GSEA requires a ranking column ('{rank_col}') in the marker DataFrame")

    # Run enrichment for each group
    enrichment_results: Dict[str, Dict[str, pd.DataFrame]] = {}

    for cluster in group_order:
        cluster_results: Dict[str, pd.DataFrame] = {}
        sub = marker_df[marker_df["group"] == cluster]

        if sub.empty:
            log.warning(f"Skipping '{cluster}': no marker genes found in '{config.de_key}'")
            continue

        # === ORA (Over-Representation Analysis) ===
        if config.method in ["ora", "both"]:
            gene_list = (
                sub.sort_values(rank_col, ascending=False)["names"]
                .head(config.n_top_genes_ora)
                .astype(str)
                .tolist()
            )

            if len(gene_list) < config.min_genes_for_ora:
                log.warning(
                    f"Skipping ORA for '{cluster}': only {len(gene_list)} genes "
                    f"(< {config.min_genes_for_ora})"
                )
            else:
                all_ora_results = []

                for category, gmt in gmt_files_to_run.items():
                    try:
                        if config.mode == "online":
                            enr_ora = gp.enrichr(
                                gene_list=gene_list,
                                gene_sets=gmt,
                                organism=config.organism,
                                background=len(background_genes),
                                outdir=None,
                                cutoff=1.0,
                            )
                        else:
                            enr_ora = gp.enrich(
                                gene_list=gene_list,
                                gene_sets=gmt,
                                background=len(background_genes),
                                outdir=None,
                                cutoff=1.0,
                            )

                        res = enr_ora.results.copy()
                        res["gene_set_category"] = category
                        all_ora_results.append(res)

                    except Exception as e:
                        log.error(
                            f"ORA failed for cluster '{cluster}', " f"category '{category}': {e}"
                        )

                if all_ora_results:
                    ora_df = pd.concat(all_ora_results, ignore_index=True)
                    ora_df = standardize_enrichment_cols(ora_df)

                    if "pval_adj" in ora_df.columns:
                        ora_df = ora_df[ora_df["pval_adj"] < config.cutoff_pval]

                    cluster_results["ora"] = ora_df
                else:
                    cluster_results["ora"] = pd.DataFrame()

        # === GSEA (Gene Set Enrichment Analysis) ===
        if config.method in ["gsea", "both"]:
            # Build ranked gene list
            rnk = (
                sub.drop_duplicates(subset="names", keep="first")
                .set_index("names")[rank_col]
                .sort_values(ascending=False)
            )

            if rnk.empty:
                log.warning(f"Skipping GSEA for '{cluster}': no ranked genes available")
            else:
                all_gsea_results = []

                for category, gmt in gmt_files_to_run.items():
                    try:
                        gsea_res = gp.prerank(
                            rnk=rnk,
                            gene_sets=gmt,
                            permutation_num=config.gsea_permutations,
                            min_size=config.gsea_min_size,
                            max_size=config.gsea_max_size,
                            outdir=None,
                            seed=42,
                            processes=4,
                        )

                        res = gsea_res.res2d.copy()
                        res["gene_set_category"] = category
                        all_gsea_results.append(res)

                    except Exception as e:
                        log.error(
                            f"GSEA failed for cluster '{cluster}', " f"category '{category}': {e}"
                        )

                if all_gsea_results:
                    gsea_df = pd.concat(all_gsea_results, ignore_index=True)
                    gsea_df = _standardize_enrichment_cols(gsea_df)

                    if "pval_adj" in gsea_df.columns:
                        gsea_df = gsea_df[gsea_df["pval_adj"] < config.cutoff_pval]

                    cluster_results["gsea"] = gsea_df
                else:
                    cluster_results["gsea"] = pd.DataFrame()

        enrichment_results[str(cluster)] = cluster_results

        # Save individual results
        if config.save_dir:
            Path(config.save_dir).mkdir(parents=True, exist_ok=True)
            for method, res_df in cluster_results.items():
                if not res_df.empty:
                    output_path = (
                        Path(config.save_dir)
                        / f"{_safe_filename(str(cluster))}_{method}_enrichment.csv"
                    )
                    res_df.to_csv(output_path, index=False)

    # Store results
    store_root = adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("de", {})
    store_root[config.key_added] = {
        "results": enrichment_results,
        "params": sanitize_for_hdf5(config.to_dict()),
    }

    if config.verbose:
        log.info(
            f"Enrichment analysis complete: {len(enrichment_results)} groups analyzed. "
            f"Results stored at .uns['...']['{config.key_added}']"
        )

    return enrichment_results


def export_enrichment_results(
    adata: AnnData,
    enrichment_key: str = "enrichment",
    output_path: str = "enrichment_results.xlsx",
) -> None:
    """
    Export enrichment results to Excel with separate sheets per cluster/method.

    Args:
        adata: AnnData with enrichment results
        enrichment_key: Key for enrichment results in adata.uns
        output_path: Output Excel file path

    Example:
        >>> run_enrichment(adata, groupby="leiden", config=enr_config)
        >>> export_enrichment_results(adata, output_path="enrichment.xlsx")
    """
    root = adata.uns.get("sclucid", {}).get("analysis", {}).get("de", {})

    if enrichment_key not in root:
        raise KeyError(
            f"Enrichment results '{enrichment_key}' not found. " "Run run_enrichment() first."
        )

    enr_store = root[enrichment_key]
    if "results" not in enr_store:
        raise ValueError("Enrichment store missing 'results' key")

    enrichment_results = enr_store["results"]

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for cluster, methods_dict in enrichment_results.items():
            for method, df in methods_dict.items():
                if not df.empty:
                    sheet_name = f"{_safe_filename(cluster)}_{method}"[:31]  # Excel limit
                    df.to_excel(writer, sheet_name=sheet_name, index=False)

    log.info(f"Enrichment results exported to {output_path}")


# ==================== Batch Analysis ====================


def batch_celltype_deg_enrichment(
    adata: AnnData,
    celltype_col: str,
    condition_col: str,
    condition1: str,
    condition2: str,
    outdir: str,
    celltypes: Optional[List[str]] = None,
    min_cells: int = 20,
    # DE parameters
    de_method: str = "wilcoxon",
    min_log2fc: float = 0.5,
    max_padj: float = 0.05,
    min_pct: float = 0.1,
    # Enrichment parameters
    run_enrichment_analysis: bool = True,
    gene_sets: Optional[List[str]] = None,
    organism: str = "human",
    enrichment_mode: str = "online",
    # Visualization
    plot_volcano_charts: bool = True,
    plot_enrichment_charts: bool = True,
    # Other
    save_pickle: bool = True,
    verbose: bool = True,
) -> Dict[str, Dict[str, Any]]:
    """
    Batch DEG and enrichment analysis across multiple cell types.

    For each cell type:
    1. Subset cells
    2. Run condition comparison (condition1 vs condition2)
    3. Optional: Run enrichment on up/down-regulated genes
    4. Optional: Generate volcano plots
    5. Save results

    Args:
        adata: AnnData object
        celltype_col: Column in obs with cell type labels
        condition_col: Column in obs with condition labels
        condition1: First condition (numerator in fold change)
        condition2: Second condition (reference/denominator)
        outdir: Output directory for results
        celltypes: Specific cell types to analyze (None = all)
        min_cells: Minimum cells required per cell type

        de_method: DE test method
        min_log2fc: Minimum log2 fold change threshold
        max_padj: Maximum adjusted p-value threshold
        min_pct: Minimum expression percentage

        run_enrichment_analysis: Whether to run enrichment
        gene_sets: Gene sets for enrichment (default: GO_BP)
        organism: Species for enrichment
        enrichment_mode: 'online' or 'offline'

        plot_volcano_charts: Generate volcano plots
        plot_enrichment_charts: Generate enrichment bar plots

        save_pickle: Save complete results as pickle
        verbose: Verbose logging

    Returns:
        Dictionary mapping cell types to results:
        {
            "T_cells": {
                "degs": DataFrame,
                "sig_degs": DataFrame (filtered),
                "enr_up": Enrichr object,
                "enr_down": Enrichr object,
                "n_cells": int
            },
            ...
        }

    Example:
        >>> results = batch_celltype_deg_enrichment(
        ...     adata,
        ...     celltype_col="celltype",
        ...     condition_col="treatment",
        ...     condition1="Treated",
        ...     condition2="Control",
        ...     outdir="./deg_results",
        ...     min_log2fc=1.0,
        ...     max_padj=0.01
        ... )
    """
    # Create output directory
    Path(outdir).mkdir(parents=True, exist_ok=True)

    # Determine cell types
    if celltypes is None:
        celltypes = adata.obs[celltype_col].unique().tolist()

    if verbose:
        log.info(
            f"Batch DEG analysis: {len(celltypes)} cell types, " f"{condition1} vs {condition2}"
        )

    # Default gene sets
    if gene_sets is None:
        gene_sets = ["GO_Biological_Process_2023"]

    results = {}

    for celltype in celltypes:
        safe_celltype = _safe_filename(celltype)

        if verbose:
            log.info(f"\n{'=' * 60}")
            log.info(f"Processing: {celltype}")
            log.info(f"{'=' * 60}")

        # Subset data
        adata_sub = adata[adata.obs[celltype_col] == celltype].copy()
        adata_sub = adata_sub[adata_sub.obs[condition_col].isin([condition1, condition2])]

        n_cells = adata_sub.n_obs

        if n_cells < min_cells:
            log.warning(f"Skipping {celltype}: only {n_cells} cells (< {min_cells})")
            continue

        # Differential expression
        de_config = CompareConditionsConfig(
            groupby=celltype_col,
            group_name=celltype,
            condition_key=condition_col,
            condition1=condition1,
            condition2=condition2,
            method=de_method,
            min_log2fc=min_log2fc,
            max_padj=max_padj,
            min_pct=min_pct,
            save_dir=outdir,
            plot=False,
            verbose=verbose,
        )

        try:
            degs = compare_conditions(adata_sub, config=de_config)
        except Exception as e:
            log.error(f"DE analysis failed for {celltype}: {e}")
            continue

        # Save DEG table
        deg_table_path = Path(outdir) / f"{safe_celltype}_DEG_{condition1}_vs_{condition2}.csv"
        degs.to_csv(deg_table_path, index=False)

        # Volcano plot
        if plot_volcano_charts and len(degs) > 0:
            volcano_path = (
                Path(outdir) / f"{safe_celltype}_volcano_{condition1}_vs_{condition2}.pdf"
            )

            try:
                plot_volcano(
                    degs_df=degs,
                    title=f"DEG: {celltype}",
                    subtitle=f"{condition1} vs {condition2} (n={n_cells} cells)",
                    lfc_threshold=min_log2fc,
                    pval_threshold=max_padj,
                    savepath=str(volcano_path),
                )
            except Exception as e:
                log.warning(f"Volcano plot failed for {celltype}: {e}")

        # Enrichment analysis
        enr_up, enr_down = None, None
        sig_degs = None

        if run_enrichment_analysis and len(degs) > 0:
            sig_degs = degs[
                (degs["pvals_adj"] < max_padj) & (degs["logfoldchanges"].abs() > min_log2fc)
            ]

            up_genes = sig_degs[sig_degs["logfoldchanges"] > min_log2fc]["names"].tolist()

            down_genes = sig_degs[sig_degs["logfoldchanges"] < -min_log2fc]["names"].tolist()

            # Up-regulated enrichment
            if len(up_genes) > 5:
                try:
                    enr_up = gp.enrichr(
                        gene_list=up_genes,
                        gene_sets=gene_sets,
                        organism=organism,
                        outdir=None,
                    )

                    enr_up_df = enr_up.results
                    up_path = (
                        Path(outdir)
                        / f"{safe_celltype}_enrichment_up_{condition1}_vs_{condition2}.csv"
                    )
                    enr_up_df.to_csv(up_path, index=False)

                    # Visualization
                    if plot_enrichment_charts and not enr_up_df.empty:
                        plt.figure(figsize=(10, 8))
                        gp.barplot(
                            enr_up_df,
                            title=f"{celltype} UP ({condition1})",
                            top_term=20,
                            cutoff=1,
                        )
                        plt.tight_layout()
                        plt.savefig(
                            Path(outdir)
                            / f"{safe_celltype}_enrichment_up_{condition1}_vs_{condition2}.pdf",
                            dpi=300,
                        )
                        plt.close()

                except Exception as e:
                    log.warning(f"Enrichment (up) failed for {celltype}: {e}")

            # Down-regulated enrichment
            if len(down_genes) > 5:
                try:
                    enr_down = gp.enrichr(
                        gene_list=down_genes,
                        gene_sets=gene_sets,
                        organism=organism,
                        outdir=None,
                    )

                    enr_down_df = enr_down.results
                    down_path = (
                        Path(outdir)
                        / f"{safe_celltype}_enrichment_down_{condition1}_vs_{condition2}.csv"
                    )
                    enr_down_df.to_csv(down_path, index=False)

                    # Visualization
                    if plot_enrichment_charts and not enr_down_df.empty:
                        plt.figure(figsize=(10, 8))
                        gp.barplot(
                            enr_down_df,
                            title=f"{celltype} DOWN ({condition2})",
                            top_term=20,
                            cutoff=1,
                        )
                        plt.tight_layout()
                        plt.savefig(
                            Path(outdir)
                            / f"{safe_celltype}_enrichment_down_{condition1}_vs_{condition2}.pdf",
                            dpi=300,
                        )
                        plt.close()

                except Exception as e:
                    log.warning(f"Enrichment (down) failed for {celltype}: {e}")

        # Collect results
        results[celltype] = {
            "degs": degs,
            "sig_degs": sig_degs,
            "enr_up": enr_up,
            "enr_down": enr_down,
            "n_cells": n_cells,
        }

    # Save summary
    if save_pickle:
        pickle_path = Path(outdir) / f"all_DEG_enrichment_{condition1}_vs_{condition2}.pkl"
        with open(pickle_path, "wb") as f:
            pickle.dump(results, f, protocol=pickle.HIGHEST_PROTOCOL)
        log.info(f"Complete results saved: {pickle_path}")

    # Summary statistics
    if verbose:
        log.info(f"\n{'=' * 60}")
        log.info("Summary:")
        log.info(f"{'=' * 60}")
        for ct, res in results.items():
            n_degs = len(res["degs"]) if res["degs"] is not None else 0
            log.info(f"{ct}: {n_degs} DEGs, {res['n_cells']} cells")

    return results


# ==================== Advanced Analysis Functions ====================
