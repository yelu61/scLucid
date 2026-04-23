"""
Plotting functions for single-cell RNA-seq data.
"""

import logging
from typing import Any, Dict, List, Optional, Union

import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.spatial import distance

# Try importing adjustText softly
try:
    from adjustText import adjust_text
except ImportError:
    adjust_text = None

log = logging.getLogger(__name__)


def _subset_adata(adata: sc.AnnData, subset: Optional[pd.Series]) -> sc.AnnData:
    """Helper to safely subset AnnData."""
    if subset is not None:
        if (
            not isinstance(subset, pd.Series)
            or not subset.dtype == bool
            or len(subset) != adata.n_obs
        ):
            raise ValueError("`subset` must be a boolean pandas Series of length adata.n_obs.")
        log.info(f"Subsetting data to {subset.sum()} cells.")
        return adata[subset].copy()
    return adata


def _combine_groupby(adata: sc.AnnData, groupby_main: str, groupby_sub: str) -> str:
    """Helper to create a combined groupby column and return its name."""
    if groupby_main not in adata.obs or groupby_sub not in adata.obs:
        raise ValueError("groupby_main and groupby_sub must be columns in adata.obs.")

    combined_col_name = f"{groupby_main}_{groupby_sub}"
    adata.obs[combined_col_name] = (
        adata.obs[groupby_main].astype(str) + "_" + adata.obs[groupby_sub].astype(str)
    )

    main_cats = adata.obs[groupby_main].cat.categories
    sub_cats = adata.obs[groupby_sub].cat.categories
    combined_order = [f"{m}_{s}" for m in main_cats for s in sub_cats]

    adata.obs[combined_col_name] = pd.Categorical(
        adata.obs[combined_col_name], categories=combined_order, ordered=True
    )
    return combined_col_name


def _get_palette_map(
    adata: sc.AnnData, key: str, palette: Optional[Union[str, Dict]] = None
) -> Dict[str, Any]:
    """
    Robustly resolve color mapping for a categorical column.
    Priority: User Dict > adata.uns > Scanpy default generation.
    """
    if not pd.api.types.is_categorical_dtype(adata.obs[key]):
        return {}

    categories = adata.obs[key].cat.categories

    # 1. User provided dictionary
    if isinstance(palette, dict):
        # Fill missing keys with gray
        return {cat: palette.get(cat, "#cccccc") for cat in categories}

    # 2. Check adata.uns (Scanpy convention: key_colors)
    uns_key = f"{key}_colors"
    if uns_key in adata.uns:
        colors = adata.uns[uns_key]
        if len(colors) >= len(categories):
            return dict(zip(categories, colors[: len(categories)]))

    # 3. Generate new palette (Seaborn/Matplotlib)
    # If user provided a string palette name (e.g., 'tab20'), use it
    palette_name = palette if isinstance(palette, str) else "tab20"
    if len(categories) > 20 and palette_name == "tab20":
        palette_name = "husl"  # Fallback for many categories

    generated_colors = sns.color_palette(palette_name, n_colors=len(categories))
    return dict(zip(categories, generated_colors))


def _sort_genes_within_categories(
    agg_df: pd.DataFrame, marker_dict: Dict[str, List[str]]
) -> List[str]:
    """
    Cluster genes within categories based on their expression across groups.

    Parameters
    ----------
    agg_df : pd.DataFrame
        Aggregated expression data (Groups × Genes)
    marker_dict : Dict[str, List[str]]
        Gene categories

    Returns:
    -------
    List[str]
        Sorted gene names
    """
    sorted_genes = []

    for category, genes in marker_dict.items():
        valid_genes = [g for g in genes if g in agg_df.columns]

        if len(valid_genes) < 3:
            sorted_genes.extend(valid_genes)
            continue

        # Subset dataframe (genes are columns here)
        sub_df = agg_df[valid_genes]
        sub_df_T = sub_df.T  # genes x cells

        sub_df_T = sub_df_T.loc[sub_df_T.var(axis=1) > 0]
        if len(sub_df_T) < 2:
            sorted_genes.extend(valid_genes)
            continue
        try:
            Z = hierarchy.linkage(
                distance.pdist(sub_df_T.values, metric="correlation"), method="average"
            )
            leaves = hierarchy.dendrogram(Z, no_plot=True)["leaves"]
            sorted_genes.extend(sub_df_T.index[leaves].tolist())
        except Exception as e:
            log.warning(f"Clustering failed for '{category}': {e}")
            sorted_genes.extend(valid_genes)

    return sorted_genes


# =============================================================================
# Embedding Visualization Functions
# =============================================================================
