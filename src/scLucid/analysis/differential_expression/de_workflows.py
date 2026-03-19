"""
High-level analysis workflows combining DE and enrichment.

This module provides integrated workflows:
- characterize_clusters: Complete cluster characterization with DE + enrichment
- summarize_markers_and_enrichment: Summarize and visualize results
""",

import logging
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple, Union

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from anndata import AnnData

from ...base_config import SclucidBaseConfig
from ..config import DifferentialConfig, EnrichmentConfig, FilterMarkersConfig
from .de_core import find_markers, filter_markers
from .enrichment import run_enrichment

log = logging.getLogger(__name__)


def characterize_clusters(
    adata: AnnData,
    groupby: str,
    de_config: Optional[DifferentialConfig] = None,
    enrichment_config: Optional[EnrichmentConfig] = None,
    key_added: str = "cluster_characterization",
) -> AnnData:
    """
    Comprehensive cluster characterization: DE + enrichment in one step.

    This is a convenience function that:
    1. Runs find_markers() for each cluster
    2. Runs run_enrichment() on the markers
    3. Stores combined results for easy access

    Args:
        adata: AnnData object
        groupby: Column for clustering (e.g., 'leiden')
        de_config: Optional DE configuration
        enrichment_config: Optional enrichment configuration
        key_added: Key for storing combined results

    Returns:
        Modified AnnData with results in adata.uns[key_added]

    Example:
        >>> adata = characterize_clusters(
        ...     adata,
        ...     groupby="leiden",
        ...     de_config=DifferentialConfig(use_raw=True),
        ...     enrichment_config=EnrichmentConfig(mode="offline")
        ... )
        >>> # Access cluster 0 characterization:
        >>> cluster0 = adata.uns["cluster_characterization"]["results"]["0"]
        >>> cluster0_markers = cluster0["top_de_genes"]
        >>> cluster0_pathways = cluster0["enrichment_ora"]
    """
    log.info(f"Characterizing clusters in '{groupby}'...")

    # Set up DE config
    if de_config is None:
        de_config = DifferentialConfig(groupby=groupby, use_raw=True)
    else:
        de_config.groupby = groupby

    de_key = de_config.key_added or "rank_genes_groups"
    de_df_key = f"{de_key}_df"

    # Run DE
    find_markers(adata, config=de_config)

    # Set up enrichment config
    if enrichment_config is None:
        enrichment_config = EnrichmentConfig(
            de_key=de_df_key, mode="offline", method="ora"
        )
    else:
        enrichment_config.de_key = de_df_key

    # Run enrichment
    enrichment_results_dict = run_enrichment(
        adata, groupby=groupby, config=enrichment_config
    )

    # Combine results
    clusters = (
        adata.obs[groupby].cat.categories
        if pd.api.types.is_categorical_dtype(adata.obs[groupby])
        else pd.unique(adata.obs[groupby])
    )

    de_df = adata.uns["sclucid"]["analysis"]["de"][de_df_key]
    characterization_results = {}

    for cluster in clusters:
        enr_res = enrichment_results_dict.get(str(cluster), {})

        characterization_results[str(cluster)] = {
            "top_de_genes": de_df[de_df["group"] == cluster],
            "enrichment_ora": enr_res.get("ora", pd.DataFrame()),
            "enrichment_gsea": enr_res.get("gsea", pd.DataFrame()),
        }

    # Store combined results
    adata.uns[key_added] = {
        "results": characterization_results,
        "params": {
            "groupby": groupby,
            "de_df_key": de_df_key,
            "enrichment_key": enrichment_config.key_added,
            "de_params": de_config.to_dict(),
            "enrichment_params": enrichment_config.to_dict(),
        },
    }

    log.info(f"Cluster characterization complete -> adata.uns['{key_added}']")
    return adata


def summarize_markers_and_enrichment(
    adata: AnnData,
    groupby: str,
    markers_df: Optional[pd.DataFrame] = None,
    enrichment_dict: Optional[Dict[str, pd.DataFrame]] = None,
    markers_key: str = "rank_genes_groups_df",
    enrichment_key: str = "enrichment",
    enrichment_method_to_summarize: Literal["ora", "gsea"] = "ora",
    n_markers: int = 25,
    n_terms: int = 10,
    summary_file: Optional[str] = None,
    sort_markers_by: str = "logfoldchanges",
    enrichment_padj_cutoff: float = 0.05,
) -> Dict[str, str]:
    """
    Generate human-readable Markdown summaries for AI/manual annotation.

    Creates per-cluster summaries with:
    - Top marker genes
    - Top enriched pathways

    Ideal for:
    - Feeding to LLMs for automated annotation
    - Quick manual review
    - Documentation

    Args:
        adata: AnnData with DE and enrichment results
        groupby: Grouping column
        markers_df: Optional marker DataFrame (auto-loaded if None)
        enrichment_dict: Optional enrichment dict (auto-loaded if None)
        markers_key: Key for markers in adata.uns
        enrichment_key: Key for enrichment in adata.uns
        enrichment_method_to_summarize: 'ora' or 'gsea'
        n_markers: Number of top markers to include
        n_terms: Number of top terms to include
        summary_file: Optional output file path
        sort_markers_by: Column for sorting markers
        enrichment_padj_cutoff: P-value cutoff for pathways

    Returns:
        Dictionary mapping cluster names to markdown summaries

    Example:
        >>> summaries = summarize_markers_and_enrichment(
        ...     adata,
        ...     groupby="leiden",
        ...     summary_file="cluster_summaries.md"
        ... )
        >>> print(summaries["0"])
        ### Cluster 0
        **Top Markers**: CD3D, CD3E, CD3G, IL7R, TRAC
        **Top ORA Pathways**: T cell activation, T cell differentiation, ...
    """
    # Load markers
    if markers_df is None:
        try:
            log.info(f"Auto-retrieving markers from '{markers_key}'")
            markers_df = adata.uns["sclucid"]["analysis"]["de"][markers_key]
        except KeyError:
            raise KeyError(
                f"Marker DataFrame not found at "
                f".uns['sclucid']['analysis']['de']['{markers_key}']"
            )

    if markers_df is None or markers_df.empty:
        log.warning("Marker DataFrame is empty. Cannot summarize.")
        markers_df = pd.DataFrame(columns=["group", "names", "logfoldchanges"])

    # Load enrichment
    if enrichment_dict is None:
        enrichment_dict = {}
        try:
            log.info(f"Auto-retrieving enrichment from '{enrichment_key}'")
            enr_store = adata.uns["sclucid"]["analysis"]["de"].get(enrichment_key, {})

            if isinstance(enr_store, dict) and "results" in enr_store:
                enrichment_dict = enr_store["results"]
            else:
                enrichment_dict = enr_store

        except KeyError:
            log.warning("No enrichment data found. Pathways will be 'N/A'")

    # Determine group order
    if groupby in adata.obs:
        group_order = list(pd.unique(adata.obs[groupby].astype(str)))
    else:
        group_order = list(
            pd.unique(markers_df.get("group", pd.Series([], dtype=str)).astype(str))
        )

    # Determine sort column
    sort_col = (
        sort_markers_by
        if sort_markers_by in markers_df.columns
        else "scores"
        if "scores" in markers_df.columns
        else "logfoldchanges"
    )

    if sort_col not in markers_df.columns:
        log.warning(f"Sort column '{sort_col}' not found. Using first column.")
        sort_col = markers_df.columns[0] if not markers_df.empty else "names"

    # Helper: Detect term/pval columns
    def _detect_cols(df: pd.DataFrame) -> Tuple[Optional[str], Optional[str]]:
        term_candidates = [
            "Term",
            "term_name",
            "term",
            "Description",
            "Name",
            "Pathway",
        ]
        pval_candidates = ["pval_adj", "Adjusted P-value", "p.adjust", "padj", "FDR"]

        tcol = next((c for c in term_candidates if c in df.columns), None)
        pcol = next((c for c in pval_candidates if c in df.columns), None)

        return tcol, pcol

    # Build summaries
    summaries: Dict[str, str] = {}
    lines: List[str] = []

    for g in group_order:
        # Top markers
        group_markers = markers_df[markers_df["group"].astype(str) == str(g)]
        top_genes = (
            group_markers.sort_values(sort_col, ascending=False)["names"]
            .head(n_markers)
            .astype(str)
            .tolist()
        )

        # Top pathways
        top_terms: List[str] = []
        cluster_enr_dict = enrichment_dict.get(str(g), {})
        enr_df = cluster_enr_dict.get(enrichment_method_to_summarize, pd.DataFrame())

        if isinstance(enr_df, pd.DataFrame) and not enr_df.empty:
            term_col, pval_col = _detect_cols(enr_df)

            if term_col and pval_col:
                tmp = enr_df.copy()
                tmp[pval_col] = pd.to_numeric(tmp[pval_col], errors="coerce")
                sig = tmp.dropna(subset=[pval_col])
                sig = sig[sig[pval_col] < float(enrichment_padj_cutoff)]

                if not sig.empty:
                    # GSEA: sort by NES (descending)
                    # ORA: sort by p-value (ascending)
                    sort_col_enr = pval_col
                    ascending = True

                    if (
                        enrichment_method_to_summarize == "gsea"
                        and "nes" in sig.columns
                    ):
                        sort_col_enr = "nes"
                        ascending = False

                    top_terms = (
                        sig.sort_values(sort_col_enr, ascending=ascending)[term_col]
                        .head(n_terms)
                        .astype(str)
                        .tolist()
                    )
                else:
                    log.debug(
                        f"Cluster {g}: No pathways below p_adj={enrichment_padj_cutoff}"
                    )
            else:
                log.warning(
                    f"Cluster {g}: Could not detect term/pval columns. "
                    f"Columns: {list(enr_df.columns)}"
                )

        # Format summary
        title = f"### Cluster {g}"
        mk_str = f"**Top Markers**: {', '.join(top_genes) if top_genes else 'N/A'}"
        pt_str = (
            f"**Top {enrichment_method_to_summarize.upper()} Pathways**: "
            f"{', '.join(top_terms) if top_terms else 'N/A'}"
        )

        summary_text = f"{title}\n{mk_str}\n{pt_str}"
        summaries[str(g)] = summary_text
        lines.append(summary_text)

    # Write to file
    if summary_file:
        Path(summary_file).parent.mkdir(parents=True, exist_ok=True)
        content = (
            "\n\n---\n\n".join(lines)
            if lines
            else "# Marker and Enrichment Summary\n\nNo results."
        )

        with open(summary_file, "w", encoding="utf-8") as f:
            f.write(content)

        log.info(f"Summaries exported to {summary_file}")

    return summaries


# ==================== Visualization Functions ====================

