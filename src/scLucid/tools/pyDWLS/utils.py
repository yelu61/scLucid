"""
Utility functions for pyDWLS (R-free)

Pure helpers for gene alignment, expression normalization, gene filtering,
and pseudo-bulk construction used by the DWLS deconvolution pipeline.
"""

import logging
from typing import Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

_NormMethod = Literal["cpm", "log1p", "library_size"]


def align_data(
    sig_df: pd.DataFrame,
    bulk_df: pd.DataFrame,
    min_common_genes: int = 10,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Align signature matrix and bulk data on a shared gene index.

    Performs a case-insensitive intersection of gene names and reindexes
    both inputs to the common gene order. Mirrors the convention used by
    ``scLucid.tools.pyBayesPrism.core.BayesPrism._align_genes``.

    Parameters
    ----------
    sig_df : pd.DataFrame
        Signature matrix (genes x cell types).
    bulk_df : pd.DataFrame
        Bulk expression matrix (genes x samples).
    min_common_genes : int, default=10
        Minimum number of overlapping genes required.

    Returns:
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        ``(sig_aligned, bulk_aligned)`` — both indexed by the same gene list
        (using the casing from ``sig_df``).

    Raises:
    ------
    ValueError
        If fewer than ``min_common_genes`` genes overlap.

    Examples:
    --------
    >>> sig_aligned, bulk_aligned = align_data(signature, bulk)
    >>> assert (sig_aligned.index == bulk_aligned.index).all()
    """
    sig_lookup = {str(g).upper(): g for g in sig_df.index}
    bulk_lookup = {str(g).upper(): g for g in bulk_df.index}
    common_upper = set(sig_lookup) & set(bulk_lookup)

    if len(common_upper) < min_common_genes:
        raise ValueError(
            f"Only {len(common_upper)} common genes between signature ({sig_df.shape[0]}) "
            f"and bulk ({bulk_df.shape[0]}); need at least {min_common_genes}."
        )

    sig_keys = [sig_lookup[g] for g in common_upper]
    bulk_keys = [bulk_lookup[g] for g in common_upper]

    sig_aligned = sig_df.loc[sig_keys].copy()
    bulk_aligned = bulk_df.loc[bulk_keys].copy()
    bulk_aligned.index = sig_aligned.index

    log.info(
        "Aligned signature (%d) and bulk (%d) on %d common genes",
        sig_df.shape[0],
        bulk_df.shape[0],
        sig_aligned.shape[0],
    )

    return sig_aligned, bulk_aligned


def normalize_data(
    data: pd.DataFrame,
    method: _NormMethod = "cpm",
) -> pd.DataFrame:
    """
    Normalize expression data column-wise.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix (genes x samples or cells).
    method : {"cpm", "log1p", "library_size"}, default="cpm"
        Normalization strategy:
        - ``"cpm"``: counts per million (column sums become exactly 1e6).
        - ``"log1p"``: ``log1p`` after CPM normalization.
        - ``"library_size"``: divide by column sums (column sums become 1).

    Returns:
    -------
    pd.DataFrame
        Normalized matrix with the same shape and labels.

    Raises:
    ------
    ValueError
        If method is not recognized or all columns are empty.

    Examples:
    --------
    >>> cpm = normalize_data(counts, method="cpm")
    >>> assert np.allclose(cpm.sum(axis=0), 1e6)
    """
    if method not in ("cpm", "log1p", "library_size"):
        raise ValueError(f"Unknown normalization method: {method!r}")

    col_sums = data.sum(axis=0)
    if (col_sums <= 0).any():
        zero_cols = col_sums[col_sums <= 0].index.tolist()
        raise ValueError(f"Columns with zero total expression: {zero_cols}")

    if method == "cpm":
        return data.div(col_sums, axis=1) * 1e6
    if method == "library_size":
        return data.div(col_sums, axis=1)
    # log1p path
    cpm = data.div(col_sums, axis=1) * 1e6
    return np.log1p(cpm)


def filter_genes(
    data: pd.DataFrame,
    min_cells: int = 5,
    min_expression: float = 1.0,
) -> pd.DataFrame:
    """
    Drop genes that are not expressed above ``min_expression`` in at least
    ``min_cells`` cells/samples.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix (genes x cells/samples).
    min_cells : int, default=5
        Minimum number of columns where the gene must exceed ``min_expression``.
    min_expression : float, default=1.0
        Per-cell expression threshold for counting.

    Returns:
    -------
    pd.DataFrame
        Filtered matrix with rows (genes) failing the threshold removed.

    Examples:
    --------
    >>> filtered = filter_genes(counts, min_cells=10, min_expression=1)
    """
    keep_mask = (data >= min_expression).sum(axis=1) >= min_cells
    kept = data.loc[keep_mask]
    log.info("filter_genes: %d/%d genes retained", kept.shape[0], data.shape[0])
    return kept


def create_pseudo_bulk(
    sc_data: pd.DataFrame,
    cell_type_labels: Union[pd.Series, np.ndarray, list],
    n_cells: int = 50,
    random_state: Optional[int] = None,
) -> Tuple[pd.Series, pd.Series]:
    """
    Build a synthetic bulk sample by summing random single-cell expression.

    Useful for cross-validating signature matrices against known proportions.

    Parameters
    ----------
    sc_data : pd.DataFrame
        Single-cell expression matrix (genes x cells).
    cell_type_labels : array-like
        Cell type label per cell, aligned with ``sc_data.columns``.
    n_cells : int, default=50
        Number of cells to draw per cell type. If a cell type has fewer cells,
        all of them are used.
    random_state : int, optional
        Seed for reproducible cell sampling.

    Returns:
    -------
    Tuple[pd.Series, pd.Series]
        - ``pseudo_bulk`` indexed by gene name, summed expression across the
          sampled cells.
        - ``true_props`` indexed by cell type, summing to 1.0.

    Examples:
    --------
    >>> bulk, props = create_pseudo_bulk(sc_data, labels, n_cells=100, random_state=0)
    >>> assert np.isclose(props.sum(), 1.0)
    """
    labels = pd.Series(cell_type_labels, index=sc_data.columns)
    rng = np.random.default_rng(random_state)

    sampled_columns: list = []
    counts_per_type: dict = {}

    for cell_type in labels.unique():
        type_cells = labels.index[labels == cell_type]
        if len(type_cells) == 0:
            continue
        take = min(n_cells, len(type_cells))
        chosen = rng.choice(type_cells.to_numpy(), size=take, replace=False)
        sampled_columns.extend(chosen.tolist())
        counts_per_type[cell_type] = take

    if not sampled_columns:
        raise ValueError("No cells available to construct pseudo-bulk.")

    pseudo_bulk = sc_data.loc[:, sampled_columns].sum(axis=1)
    pseudo_bulk.name = "pseudo_bulk"

    total = sum(counts_per_type.values())
    true_props = pd.Series(
        {ct: count / total for ct, count in counts_per_type.items()},
        name="true_proportion",
    )

    log.info(
        "create_pseudo_bulk: %d cells sampled across %d cell types",
        total,
        len(counts_per_type),
    )

    return pseudo_bulk, true_props
