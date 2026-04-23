"""
High-level analysis workflows combining DE and enrichment.

This module provides integrated workflows:
- characterize_clusters: Complete cluster characterization with DE + enrichment
- summarize_markers_and_enrichment: Summarize and visualize results
"""

import logging
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple, Union

import pandas as pd
from anndata import AnnData

from ..config import DifferentialConfig, EnrichmentConfig
from .de_core import find_markers
from .enrichment import run_enrichment

log = logging.getLogger(__name__)


def _resolve_characterization_export_paths(
    save_path: Union[str, Path],
    *,
    key_added: str,
) -> Dict[str, Path]:
    """Resolve stable output paths for characterization review sidecars."""
    target = Path(save_path)
    if target.suffix:
        base_dir = target.parent
        stem = target.stem
        summary_path = target if target.suffix.lower() == ".csv" else target.with_suffix(".csv")
    else:
        base_dir = target
        stem = key_added
        summary_path = base_dir / f"{stem}_summary.csv"

    base_dir.mkdir(parents=True, exist_ok=True)
    return {
        "summary": summary_path,
        "top_markers": base_dir / f"{stem}_top_markers.csv",
        "enrichment": base_dir / f"{stem}_enrichment_summary.csv",
        "markdown": base_dir / f"{stem}_summaries.md",
    }


def _detect_enrichment_columns(df: pd.DataFrame) -> Tuple[Optional[str], Optional[str]]:
    """Detect term and adjusted p-value columns across common enrichment outputs."""
    term_candidates = ["Term", "term_name", "term", "Description", "Name", "Pathway"]
    pval_candidates = ["pval_adj", "Adjusted P-value", "p.adjust", "padj", "FDR"]
    term_col = next((c for c in term_candidates if c in df.columns), None)
    pval_col = next((c for c in pval_candidates if c in df.columns), None)
    return term_col, pval_col


def _build_characterization_tables(
    de_df: pd.DataFrame,
    enrichment_results_dict: Dict[str, Dict[str, pd.DataFrame]],
    *,
    n_top_markers: int,
    n_top_terms: int,
    enrichment_method_to_summarize: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Build notebook-friendly cluster marker and enrichment summary tables."""
    marker_rows: List[Dict[str, object]] = []
    enrichment_rows: List[Dict[str, object]] = []
    marker_columns = ["cluster", "rank", "gene", "logfoldchanges", "scores", "pvals_adj"]
    enrichment_columns = [
        "cluster",
        "method",
        "rank",
        "term",
        "adjusted_p_value",
        "source_columns",
    ]

    if de_df is not None and not de_df.empty:
        if "logfoldchanges" in de_df.columns:
            sort_col = "logfoldchanges"
        elif "scores" in de_df.columns:
            sort_col = "scores"
        else:
            sort_col = "names"

        for group, group_df in de_df.groupby(de_df["group"].astype(str)):
            ranked = group_df.sort_values(sort_col, ascending=False).head(n_top_markers)
            for rank, (_, row) in enumerate(ranked.iterrows(), start=1):
                marker_rows.append(
                    {
                        "cluster": str(group),
                        "rank": rank,
                        "gene": str(row.get("names", "")),
                        "logfoldchanges": row.get("logfoldchanges"),
                        "scores": row.get("scores"),
                        "pvals_adj": row.get("pvals_adj"),
                    }
                )

    for cluster, cluster_results in enrichment_results_dict.items():
        enr_df = cluster_results.get(enrichment_method_to_summarize, pd.DataFrame())
        if not isinstance(enr_df, pd.DataFrame) or enr_df.empty:
            continue
        term_col, pval_col = _detect_enrichment_columns(enr_df)
        if term_col is None:
            continue
        ranked = enr_df.copy()
        if pval_col is not None:
            ranked[pval_col] = pd.to_numeric(ranked[pval_col], errors="coerce")
            ranked = ranked.sort_values(pval_col, ascending=True, na_position="last")
        ranked = ranked.head(n_top_terms)
        for rank, (_, row) in enumerate(ranked.iterrows(), start=1):
            enrichment_rows.append(
                {
                    "cluster": str(cluster),
                    "method": enrichment_method_to_summarize,
                    "rank": rank,
                    "term": str(row.get(term_col, "")),
                    "adjusted_p_value": row.get(pval_col) if pval_col is not None else None,
                    "source_columns": ",".join(map(str, ranked.columns)),
                }
            )

    return (
        pd.DataFrame(marker_rows, columns=marker_columns),
        pd.DataFrame(enrichment_rows, columns=enrichment_columns),
    )


def characterize_clusters(
    adata: AnnData,
    groupby: str,
    de_config: Optional[DifferentialConfig] = None,
    enrichment_config: Optional[EnrichmentConfig] = None,
    key_added: str = "cluster_characterization",
    save_path: Optional[Union[str, Path]] = None,
    n_top_markers: int = 10,
    n_top_terms: int = 5,
    enrichment_method_to_summarize: Literal["ora", "gsea"] = "ora",
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
        save_path: Optional output directory or base CSV path for review sidecars
        n_top_markers: Number of top markers per cluster to export in review tables
        n_top_terms: Number of top enrichment terms per cluster to export
        enrichment_method_to_summarize: Enrichment method to prioritize in summaries

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
        enrichment_config = EnrichmentConfig(de_key=de_df_key, mode="offline", method="ora")
    else:
        enrichment_config.de_key = de_df_key

    # Run enrichment
    enrichment_results_dict = run_enrichment(adata, groupby=groupby, config=enrichment_config)

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

    top_markers_df, enrichment_summary_df = _build_characterization_tables(
        de_df,
        enrichment_results_dict,
        n_top_markers=n_top_markers,
        n_top_terms=n_top_terms,
        enrichment_method_to_summarize=enrichment_method_to_summarize,
    )
    markdown_summaries = summarize_markers_and_enrichment(
        adata,
        groupby=groupby,
        markers_df=de_df,
        enrichment_dict=enrichment_results_dict,
        enrichment_method_to_summarize=enrichment_method_to_summarize,
        n_markers=n_top_markers,
        n_terms=n_top_terms,
    )
    cluster_summary_df = pd.DataFrame(
        {
            "cluster": [str(cluster) for cluster in clusters],
            "n_cells": [
                int((adata.obs[groupby].astype(str) == str(cluster)).sum()) for cluster in clusters
            ],
        }
    )
    if not top_markers_df.empty:
        cluster_summary_df = cluster_summary_df.merge(
            top_markers_df.groupby("cluster")["gene"]
            .agg(lambda genes: ", ".join(pd.Series(genes).dropna().astype(str).tolist()))
            .rename("top_markers"),
            how="left",
            on="cluster",
        )
    else:
        cluster_summary_df["top_markers"] = ""

    if not enrichment_summary_df.empty:
        cluster_summary_df = cluster_summary_df.merge(
            enrichment_summary_df.groupby("cluster")["term"]
            .agg(lambda terms: ", ".join(pd.Series(terms).dropna().astype(str).tolist()))
            .rename("top_pathways"),
            how="left",
            on="cluster",
        )
    else:
        cluster_summary_df["top_pathways"] = ""

    cluster_summary_df = cluster_summary_df.fillna({"top_markers": "", "top_pathways": ""})

    # Store combined results
    adata.uns[key_added] = {
        "results": characterization_results,
        "summary_table": cluster_summary_df,
        "top_markers": top_markers_df,
        "enrichment_summary": enrichment_summary_df,
        "markdown_summary": markdown_summaries,
        "params": {
            "groupby": groupby,
            "de_df_key": de_df_key,
            "enrichment_key": enrichment_config.key_added,
            "n_top_markers": int(n_top_markers),
            "n_top_terms": int(n_top_terms),
            "summary_enrichment_method": enrichment_method_to_summarize,
            "de_params": de_config.to_dict(),
            "enrichment_params": enrichment_config.to_dict(),
        },
    }

    if save_path:
        export_paths = _resolve_characterization_export_paths(save_path, key_added=key_added)
        cluster_summary_df.to_csv(export_paths["summary"], index=False)
        top_markers_df.to_csv(export_paths["top_markers"], index=False)
        enrichment_summary_df.to_csv(export_paths["enrichment"], index=False)
        summarize_markers_and_enrichment(
            adata,
            groupby=groupby,
            markers_df=de_df,
            enrichment_dict=enrichment_results_dict,
            enrichment_method_to_summarize=enrichment_method_to_summarize,
            n_markers=n_top_markers,
            n_terms=n_top_terms,
            summary_file=str(export_paths["markdown"]),
        )
        adata.uns[key_added]["export_paths"] = {
            name: str(path) for name, path in export_paths.items()
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
        group_order = list(pd.unique(markers_df.get("group", pd.Series([], dtype=str)).astype(str)))

    # Determine sort column
    sort_col = (
        sort_markers_by
        if sort_markers_by in markers_df.columns
        else "scores" if "scores" in markers_df.columns else "logfoldchanges"
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

                    if enrichment_method_to_summarize == "gsea" and "nes" in sig.columns:
                        sort_col_enr = "nes"
                        ascending = False

                    top_terms = (
                        sig.sort_values(sort_col_enr, ascending=ascending)[term_col]
                        .head(n_terms)
                        .astype(str)
                        .tolist()
                    )
                else:
                    log.debug(f"Cluster {g}: No pathways below p_adj={enrichment_padj_cutoff}")
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
            "\n\n---\n\n".join(lines) if lines else "# Marker and Enrichment Summary\n\nNo results."
        )

        with open(summary_file, "w", encoding="utf-8") as f:
            f.write(content)

        log.info(f"Summaries exported to {summary_file}")

    return summaries


# ==================== Visualization Functions ====================
