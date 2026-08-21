"""
Utility functions for single-cell RNA-seq data analysis.

This module provides common helper functions that can be used
across different parts of the analysis pipeline.
"""

import logging
from contextlib import contextmanager
from typing import Any, Dict, List, Literal, Optional, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from .sanitize import sanitize_for_hdf5

log = logging.getLogger(__name__)


__all__ = [
    "assess_matrix_semantics",
    "use_layer_as_X",
    "sanitize_for_hdf5",
    "subset_adata",
    "subset_from_annotations",
    "merge_obs_metadata",
    "build_metadata_dicts",
    "print_sample_crosstab",
]


def assess_matrix_semantics(
    data: Union[AnnData, np.ndarray, sp.spmatrix],
    *,
    layer: Optional[str] = None,
    semantics: Literal[
        "raw_counts", "normalized", "log_normalized", "scaled", "any"
    ] = "raw_counts",
    require_non_negative: bool = True,
    require_integer: bool = True,
    max_fractional_rate: float = 0.01,
    zero_fraction_threshold: float = 0.20,
    min_max_value: float = 10.0,
    max_cells: int = 2000,
    max_genes: int = 2000,
) -> Dict[str, Any]:
    """Assess whether a matrix conforms to expected single-cell semantics.

    This is the unified entry point for matrix-type checks across scLucid.
    Callers in QC, preprocessing, doublet detection, and ambient correction
    should use this instead of re-implementing sparse/integer/negative checks.

    Parameters
    ----------
    data
        ``AnnData`` object, a dense numpy array, or a scipy sparse matrix.
    layer
        When ``data`` is an ``AnnData``, the layer to inspect. If ``None``,
        ``adata.X`` is used.
    semantics
        Expected semantics of the matrix. ``any`` only performs structural checks.
    require_non_negative
        Reject negative values for ``raw_counts``/``normalized`` semantics.
    require_integer
        Reject non-integer positive values for ``raw_counts`` semantics.
    max_fractional_rate
        Maximum fraction of positive values allowed to be non-integer when
        ``require_integer`` is True.
    zero_fraction_threshold
        Minimum fraction of zeros expected in a sparse raw-count matrix.
    min_max_value
        Minimum maximum value expected for a real count matrix.
    max_cells, max_genes
        Subsample dimensions for large matrices.

    Returns:
    -------
    dict with keys:
        - ``semantics``: requested semantics string.
        - ``is_valid``: bool, whether matrix matches requested semantics.
        - ``is_count_like``: bool, whether matrix looks like raw counts.
        - ``warnings``: list[str], human-readable problems.
        - ``diagnostics``: numeric diagnostics (zero_fraction, max_value, etc.).
        - ``matrix_shape``: tuple.
    """
    # Resolve the matrix from AnnData if needed.
    if isinstance(data, AnnData):
        if layer is None:
            X = data.X
        elif layer in data.layers:
            X = data.layers[layer]
        else:
            return {
                "semantics": semantics,
                "is_valid": False,
                "is_count_like": False,
                "warnings": [f"Layer '{layer}' not found in adata.layers"],
                "diagnostics": {},
                "matrix_shape": data.shape,
            }
    else:
        X = data

    if X is None:
        return {
            "semantics": semantics,
            "is_valid": False,
            "is_count_like": False,
            "warnings": ["Matrix is None"],
            "diagnostics": {},
            "matrix_shape": None,
        }

    n_obs, n_vars = X.shape
    if n_obs == 0 or n_vars == 0:
        return {
            "semantics": semantics,
            "is_valid": False,
            "is_count_like": False,
            "warnings": ["Matrix is empty"],
            "diagnostics": {},
            "matrix_shape": (n_obs, n_vars),
        }

    # Subsample for cheap checks on large matrices.
    if sp.issparse(X):
        sub_cells = min(n_obs, max_cells)
        sub_genes = min(n_vars, max_genes)
        Xs = X[:sub_cells, :sub_genes].copy()
        data_vec = np.asarray(Xs.data)
        zero_fraction = 1.0 - (Xs.nnz / max(1, sub_cells * sub_genes))
    else:
        arr = np.asarray(X)
        sub_cells = min(n_obs, max_cells)
        sub_genes = min(n_vars, max_genes)
        Xs = arr[:sub_cells, :sub_genes]
        data_vec = Xs.ravel()
        zero_fraction = float(np.mean(Xs == 0))

    finite = data_vec[np.isfinite(data_vec)]
    if finite.size == 0:
        return {
            "semantics": semantics,
            "is_valid": False,
            "is_count_like": False,
            "warnings": ["No finite values in matrix"],
            "diagnostics": {"zero_fraction": float(zero_fraction)},
            "matrix_shape": (n_obs, n_vars),
        }

    has_negative = bool(np.any(finite < 0))
    positive = finite[finite > 0]
    fractional_rate = (
        float(np.mean(np.abs(positive - np.round(positive)) > 1e-6)) if positive.size else 0.0
    )
    max_value = float(np.max(finite))
    min_value = float(np.min(finite))

    warnings: List[str] = []
    diagnostics = {
        "matrix_shape": (n_obs, n_vars),
        "subsampled_shape": (sub_cells, sub_genes),
        "has_negative": has_negative,
        "min_value": min_value,
        "max_value": max_value,
        "zero_fraction": float(zero_fraction),
        "fractional_positive_rate": fractional_rate,
    }

    # Count-like check using the historical criteria from the legacy raw-count checker.
    is_count_like = (
        (not has_negative)
        and fractional_rate < max_fractional_rate
        and zero_fraction >= zero_fraction_threshold
        and max_value > min_max_value
    )

    # Validate against requested semantics.
    is_valid = True
    if semantics in ("raw_counts", "normalized", "log_normalized") and require_non_negative:
        if has_negative:
            warnings.append("Matrix contains negative values")
            is_valid = False
    if semantics == "raw_counts" and require_integer:
        if fractional_rate > max_fractional_rate:
            warnings.append(
                f"Matrix contains too many non-integer positive values "
                f"({fractional_rate:.2%} > {max_fractional_rate:.2%})"
            )
            is_valid = False
    if semantics == "raw_counts":
        if zero_fraction < zero_fraction_threshold:
            warnings.append(f"Matrix is too dense for raw counts ({zero_fraction:.2%} zeros)")
            is_valid = False
        if max_value <= min_max_value:
            warnings.append(f"Matrix max value ({max_value}) is too small for raw counts")
            is_valid = False
    if semantics == "scaled":
        # Scaled matrices typically contain negative values and small ranges.
        if not has_negative and max_value > 100:
            warnings.append("Matrix does not look like z-scaled data")
            is_valid = False

    # For "any" semantics we only report diagnostics and count-like status.
    if semantics == "any":
        is_valid = True

    return {
        "semantics": semantics,
        "is_valid": is_valid,
        "is_count_like": is_count_like,
        "warnings": warnings,
        "diagnostics": diagnostics,
        "matrix_shape": (n_obs, n_vars),
    }


@contextmanager
def use_layer_as_X(adata: AnnData, layer: Optional[str]):
    """Context manager to temporarily use a layer as adata.X."""
    if layer is None:
        yield
        return

    if layer not in adata.layers:
        log.warning(f"Layer '{layer}' not found in adata.layers. Using adata.X.")
        yield
        return

    X_backup = adata.X.copy()
    adata.X = adata.layers[layer].copy()
    try:
        yield
    finally:
        # Always restore the original .X
        adata.X = X_backup


def print_sample_crosstab(
    adata: AnnData,
    sample_key: str = "sampleID",
    group_key: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    """Print a crosstab of samples by biological group.

    Useful for quick project overviews before QC or when auditing sample
    balance after filtering.

    Parameters
    ----------
    adata
        AnnData object.
    sample_key
        Column in ``.obs`` identifying samples.
    group_key
        Optional column in ``.obs`` for a biological group. If provided,
        a sample x group crosstab with margins is printed; otherwise the
        per-sample cell counts are printed.

    Returns:
    -------
    pd.DataFrame or None
        The printed crosstab DataFrame.
    """
    if sample_key not in adata.obs.columns:
        raise ValueError(f"sample_key '{sample_key}' not found in adata.obs")

    if group_key and group_key in adata.obs.columns:
        ctab = pd.crosstab(
            adata.obs[sample_key],
            adata.obs[group_key],
            margins=True,
        )
    else:
        ctab = pd.DataFrame(adata.obs[sample_key].value_counts())
        ctab.columns = ["n_cells"]

    print(ctab)
    return ctab


def build_metadata_dicts(
    samples: List[str],
    *,
    group_dict: Optional[Dict[str, Any]] = None,
    batch_dict: Optional[Dict[str, Any]] = None,
    group_key: str = "group",
    batch_key: str = "batch",
    extra_dicts: Optional[Dict[str, Dict[str, Any]]] = None,
    strict: bool = False,
    default_value: Any = None,
) -> Dict[str, Dict[str, Any]]:
    """Build a ``metadata_dicts`` mapping for ``read_10x`` / ``load_10x_data``.

    ``group_dict`` and ``batch_dict`` preserve the historical convenience API.
    ``extra_dicts`` accepts arbitrary sample-level metadata shaped as
    ``{column_name: {sample_id: value}}``. Missing sample values are filled with
    ``default_value`` unless ``strict=True``.
    """
    metadata_dicts: Dict[str, Dict[str, Any]] = {}
    sample_list = list(samples)

    def _build_column(column: str, mapping: Dict[str, Any]) -> None:
        missing = [sample for sample in sample_list if sample not in mapping]
        if strict and missing:
            raise KeyError(f"Metadata column '{column}' is missing values for samples: {missing}")
        metadata_dicts[column] = {
            sample: mapping.get(sample, default_value) for sample in sample_list
        }

    if group_dict is not None and group_key:
        _build_column(group_key, group_dict)
    if batch_dict is not None and batch_key:
        _build_column(batch_key, batch_dict)
    for column, mapping in (extra_dicts or {}).items():
        _build_column(str(column), mapping)

    return metadata_dicts


def subset_adata(
    adata: AnnData,
    filters: Dict[str, Union[Any, List[Any]]],
    keep_raw_genes: bool = True,
    raise_on_empty: bool = True,
) -> AnnData:
    """
    Subset an AnnData object based on metadata criteria, retaining raw gene data.

    This utility function is designed for flexible data exploration, allowing you
    to create a new AnnData object for a specific subset of cells (e.g., one
    sample, or only T cells) that can be re-analyzed from scratch.

    Args:
        adata: The AnnData object to subset.
        filters: Dictionary of metadata filters. Keys are column names in `adata.obs`,
                 and values are the desired value or a list of desired values.
                 Example: `{"sampleID": "sample1", "cell_type": ["T cells", "B cells"]}`
        keep_raw_genes: If True and `adata.raw` exists, the returned object's `.raw`
                        attribute will contain the subset of cells but the full original
                        set of genes, enabling re-running of HVG selection.

    Returns:
        A new, subsetted AnnData object.
    """
    if not isinstance(filters, dict):
        raise TypeError("filters must be a dictionary.")

    initial_cells = adata.n_obs
    combined_mask = pd.Series(True, index=adata.obs_names)

    for key, value in filters.items():
        if key not in adata.obs.columns:
            log.warning(f"Metadata column '{key}' not found in adata.obs. Skipping filter.")
            continue

        if isinstance(value, list):
            mask = adata.obs[key].isin(value)
        else:
            mask = adata.obs[key] == value

        combined_mask &= mask

    final_cells = combined_mask.sum()
    log.info("Subsetting data based on provided filters:")
    log.info(f"  - Initial cells: {initial_cells}")
    log.info(f"  - Final cells after filtering: {final_cells}")

    if final_cells == 0:
        msg = f"No cells remaining after applying filters: {filters}"
        if raise_on_empty:
            raise ValueError(msg)
        else:
            log.warning(msg)
            return AnnData()

    if final_cells < 10:
        log.warning(f"Only {final_cells} cells remaining. Results may be unreliable.")

    # The core slicing operation
    adata_subset = adata[combined_mask, :].copy()

    if keep_raw_genes and adata.raw is not None:
        log.info(f"Subset .raw created, retaining all {adata.raw.n_vars} original genes.")
    elif not keep_raw_genes and adata_subset.raw is not None:
        adata_subset.raw = None

    return adata_subset


def subset_from_annotations(
    adata_target: AnnData,
    adata_source: AnnData,
    filters: Dict[str, Union[Any, List[Any]]],
    columns_to_merge: Union[str, List[str]],
) -> AnnData:
    """
    Subsets a target AnnData object based on annotations from a source object.

    This is a convenience wrapper for the common sub-clustering workflow where
    annotations (e.g., cell types) are generated on a processed object but the
    subsetting needs to be done on an unprocessed object (e.g., containing
    raw counts for all genes).

    Args:
        adata_target: The AnnData object to be subsetted (e.g., the QC'd object).
        adata_source: The AnnData object containing the annotations in its .obs.
                      Must share the same cell indices as adata_target.
        filters: Dictionary of metadata filters to apply. The keys must be present
                 in the `columns_to_merge`.
        columns_to_merge: A column name or list of column names from `adata_source.obs`
                          to merge into `adata_target.obs` before filtering.

    Returns:
        A new, subsetted AnnData object.
    """
    if isinstance(columns_to_merge, str):
        columns_to_merge = [columns_to_merge]

    # --- Step 1: Merge Annotations ---
    log.info(f"Merging annotations for columns: {columns_to_merge} from source object.")

    # Check if columns exist in the source
    missing_cols = [col for col in columns_to_merge if col not in adata_source.obs]
    if missing_cols:
        raise ValueError(f"Columns {missing_cols} not found in the source AnnData object's .obs")

    annotations = adata_source.obs[columns_to_merge]

    # Use a temporary DataFrame to avoid modifying the original adata_target.obs in case of error
    obs_merged = adata_target.obs.join(annotations)

    # Validate that all cells were matched
    if obs_merged[columns_to_merge[0]].isnull().any():
        unmatched_count = obs_merged[columns_to_merge[0]].isnull().sum()
        log.warning(
            f"Found {unmatched_count} cells in the target object that were not present "
            "in the source object's annotations. These will not be selected."
        )

    # Create a temporary AnnData object with the merged obs for filtering
    temp_adata = adata_target.copy()
    temp_adata.obs = obs_merged

    # --- Step 2: Subset ---
    log.info(f"Subsetting target object based on new annotations with filters: {filters}")

    # Now we can call the original, simple subset_adata function
    adata_subset = subset_adata(temp_adata, filters=filters)

    return adata_subset


def merge_obs_metadata(
    adata: AnnData,
    metadata_path: str,
    left_on: Optional[str] = None,  # If None, uses adata.obs.index
    right_on: Optional[str] = None,  # If None, uses metadata_df.index
    how: str = "left",
    handle_duplicates: str = "warn",  # 'warn', 'error', 'overwrite'
) -> AnnData:
    """
    Merges metadata from an external file into the AnnData object's .obs DataFrame.

    Args:
        adata: The AnnData object to modify.
        metadata_path: Path to the metadata file (.csv, .tsv, or .xlsx).
        left_on: Column in adata.obs to join on. If None, uses the index (cell barcodes).
        right_on: Column in the external file to join on. If None, uses the index.
        how: How to perform the merge (e.g., 'left', 'inner'). Defaults to 'left'.

    Returns:
        The AnnData object with merged metadata (modified in place).
    """
    log.info(f"Loading metadata from {metadata_path}")
    if metadata_path.endswith(".csv"):
        meta_df = pd.read_csv(metadata_path)
    elif metadata_path.endswith((".xlsx", ".xls")):
        meta_df = pd.read_excel(metadata_path)
    elif metadata_path.endswith(".tsv"):
        meta_df = pd.read_csv(metadata_path, sep="\t")
    else:
        raise ValueError("Unsupported file format. Please use .csv, .tsv, or .xlsx.")

    initial_cols = set(adata.obs.columns)
    meta_cols = set(meta_df.columns if right_on is None else meta_df.columns.drop(right_on))

    overlapping = initial_cols & meta_cols

    if overlapping:
        if handle_duplicates == "error":
            raise ValueError(f"Columns already exist in adata.obs: {overlapping}")
        elif handle_duplicates == "warn":
            log.warning(
                f"Columns {overlapping} already exist. New columns will be suffixed with '_new'"
            )
            suffixes = ("", "_new")
        elif handle_duplicates == "overwrite":
            log.info(f"Overwriting columns: {overlapping}")
            # Drop existing columns before merge
            adata.obs.drop(columns=overlapping, inplace=True)
            suffixes = ("", "")
    else:
        suffixes = ("", "")

    # Perform merge
    if left_on is None:
        adata.obs = adata.obs.join(
            meta_df.set_index(right_on) if right_on else meta_df,
            how=how,
            rsuffix="_new" if handle_duplicates == "warn" else "",
        )
    else:
        adata.obs = adata.obs.merge(
            meta_df, left_on=left_on, right_on=right_on, how=how, suffixes=suffixes
        )

    return adata
