"""
Utility functions for single-cell RNA-seq data analysis.

This module provides common helper functions that can be used
across different parts of the analysis pipeline.
"""

import gc
import logging
import os
import time
from contextlib import contextmanager
from typing import Any, Dict, List, Optional, Union

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import io

log = logging.getLogger(__name__)

__all__ = [
    "load_10x_data",
    "use_layer_as_X",
    "sanitize_for_hdf5",
    "subset_adata",
    "subset_from_annotations",
    "merge_obs_metadata"
]


def _find_sample_paths(
    base_dir: str, samples: List[str], possible_subpaths: Optional[List[str]] = None
) -> Dict[str, str]:
    """
    Find valid 10x Genomics data paths for a list of samples.

    Args:
        base_dir: Base directory containing the sample folders.
        samples: List of sample IDs to search for.
        possible_subpaths: List of possible subdirectory structures to check.
                          Default paths checked are:
                          - "outs/filtered_feature_bc_matrix"
                          - "filtered_feature_bc_matrix"
                          - "" (sample directory itself)

    Returns:
        Dictionary mapping sample IDs to their valid data paths.
    """
    if possible_subpaths is None:
        possible_subpaths = [
            os.path.join("outs", "filtered_feature_bc_matrix"),
            "filtered_feature_bc_matrix",
            "",  # Use sample directory directly
        ]

    found_paths = {}

    for sample in samples:
        found = False
        for subpath in possible_subpaths:
            full_path = os.path.join(base_dir, sample, subpath)
            if os.path.isdir(full_path):
                # Check if directory contains 10x data files
                mtx_files = ["matrix.mtx", "matrix.mtx.gz"]
                feature_files = [
                    "features.tsv",
                    "features.tsv.gz",
                    "genes.tsv",
                    "genes.tsv.gz",
                ]

                has_mtx = any(
                    os.path.exists(os.path.join(full_path, f)) for f in mtx_files
                )
                has_features = any(
                    os.path.exists(os.path.join(full_path, f)) for f in feature_files
                )

                if has_mtx and has_features:
                    found_paths[sample] = full_path
                    found = True
                    break

        if not found:
            log.warning(f"No valid 10x data path found for sample {sample}")

    return found_paths


def _read_10x_manually(sample_path: str) -> AnnData:
    """
    Manually reads 10x data files as a robust fallback method.

    Args:
        sample_path: Path to the directory containing matrix.mtx.gz, features.tsv.gz, etc.

    Returns:
        An AnnData object.
    """
    log.info(f"Attempting robust manual read from: {sample_path}")

    # --- Find Files (with fallback for different names/compressions) ---
    matrix_file = os.path.join(sample_path, "matrix.mtx.gz")
    if not os.path.exists(matrix_file):
        matrix_file = os.path.join(sample_path, "matrix.mtx")

    features_file = os.path.join(sample_path, "features.tsv.gz")
    if not os.path.exists(features_file):
        features_file = os.path.join(sample_path, "genes.tsv.gz")
    if not os.path.exists(features_file):
        features_file = os.path.join(sample_path, "features.tsv")
    if not os.path.exists(features_file):
        features_file = os.path.join(sample_path, "genes.tsv")

    barcodes_file = os.path.join(sample_path, "barcodes.tsv.gz")
    if not os.path.exists(barcodes_file):
        barcodes_file = os.path.join(sample_path, "barcodes.tsv")

    if not all(os.path.exists(f) for f in [matrix_file, features_file, barcodes_file]):
        raise FileNotFoundError(
            f"Could not find all required 10x files in {sample_path}"
        )

    # --- Read Files with Explicit Type Control ---
    X = io.mmread(matrix_file).T.tocsr()

    features_df = pd.read_csv(
        features_file,
        sep="\t",
        header=None,
        compression="gzip" if features_file.endswith(".gz") else None,
        dtype=str,  # Crucial: ensure all columns are read as strings
    )
    gene_names = features_df[1] if features_df.shape[1] >= 2 else features_df[0]

    barcodes_df = pd.read_csv(
        barcodes_file,
        sep="\t",
        header=None,
        compression="gzip" if barcodes_file.endswith(".gz") else None,
        dtype=str,  # Crucial: ensure barcodes are strings
    )
    barcodes = barcodes_df[0]

    # --- Create and Sanitize AnnData Object ---
    adata = anndata.AnnData(
        X=X,
        obs=pd.DataFrame(index=barcodes.values),
        var=pd.DataFrame(index=gene_names.values),
    )

    adata.var_names_make_unique()  # Ensure gene names are unique
    adata.layers["counts"] = adata.X.copy()

    return adata


def load_10x_data(
    samples: List[str],
    base_dir: Optional[str] = None,
    path_dict: Optional[Dict[str, str]] = None,
    metadata_dicts: Optional[Dict[str, Dict[str, str]]] = None,
    possible_subpaths: Optional[List[str]] = None,
    output_file: Optional[str] = None,
    compression: Optional[str] = "gzip",
    backup_existing: bool = True,
) -> AnnData:
    """
    Load multiple 10x Genomics samples with a robust fallback mechanism.
    First tries the standard scanpy reader, then falls back to a manual method.
    """
    adata_list = []

    if path_dict is None:
        if base_dir is None:
            raise ValueError("Either base_dir or path_dict must be provided")
        log.info(f"Searching for sample paths in {base_dir}")
        path_dict = _find_sample_paths(base_dir, samples, possible_subpaths)

    sample_metadata = {}
    if metadata_dicts:
        for sample in samples:
            sample_metadata[sample] = {}
            for metadata_name, metadata_dict in metadata_dicts.items():
                if sample in metadata_dict:
                    sample_metadata[sample][metadata_name] = metadata_dict[sample]

    valid_samples = [s for s in samples if s in path_dict]
    if len(valid_samples) < len(samples):
        log.warning(
            f"Found valid paths for {len(valid_samples)}/{len(samples)} samples"
        )

    for sample in valid_samples:
        sample_path = path_dict[sample]
        adata = None

        # --- Main method with fallback ---
        try:
            log.info(f"Loading {sample} with standard method from {sample_path}")
            adata = sc.read_10x_mtx(
                sample_path, var_names="gene_symbols", cache=True, make_unique=True
            )
        except Exception as e:
            log.warning(f"Standard method failed for {sample}: {e}")
            log.info(f"Attempting robust fallback method for {sample}...")
            try:
                adata = _read_10x_manually(sample_path)
            except Exception as e2:
                log.error(f"Robust fallback method also failed for {sample}: {e2}")
                continue  # Skip to the next sample

        # --- Post-loading processing (common for both methods) ---
        if adata is not None:
            adata.obs["sampleID"] = sample
            if sample in sample_metadata:
                for meta_key, meta_value in sample_metadata[sample].items():
                    adata.obs[meta_key] = meta_value

            log.info(
                f"Successfully loaded {sample}: {adata.n_obs} cells, {adata.n_vars} genes"
            )
            adata_list.append(adata)
            gc.collect()

    if not adata_list:
        log.error("No samples were loaded successfully.")
        return AnnData()

    log.info(f"Merging {len(adata_list)} samples...")
    combined_adata = anndata.concat(
        adata_list, join="outer", keys=valid_samples, label="batch", index_unique="_"
    )

    log.info(
        f"Combined dataset: {combined_adata.n_obs} cells, {combined_adata.n_vars} genes"
    )

    if output_file:
        if os.path.exists(output_file) and backup_existing:
            backup_file = f"{output_file}.bak.{int(time.time())}"
            log.info(f"File {output_file} exists, creating backup: {backup_file}")
            os.rename(output_file, backup_file)

        log.info(f"Saving combined data to {output_file}")
        combined_adata.write(output_file, compression=compression)
        log.info(f"Data successfully saved to {output_file}")

    combined_adata.layers["counts"] = combined_adata.X.copy()

    return combined_adata


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


def sanitize_for_hdf5(obj):
    """
    Make objects HDF5-compatible by:
    1. Converting tuples to lists
    2. Converting integer keys to strings in dictionaries
    3. Handling other non-HDF5 compatible types
    """
    if isinstance(obj, tuple):
        return [sanitize_for_hdf5(item) for item in obj]
    elif isinstance(obj, list):
        return [sanitize_for_hdf5(item) for item in obj]
    elif isinstance(obj, dict):
        return {str(k): sanitize_for_hdf5(v) for k, v in obj.items()}
    elif isinstance(obj, (int, float, str, bool, np.number, np.bool_)) or obj is None:
        return obj
    else:
        # Try to convert other types to string representation
        try:
            return str(obj)
        except:
            return "Unconvertible object"


def subset_adata(
    adata: AnnData,
    filters: Dict[str, Union[Any, List[Any]]],
    keep_raw_genes: bool = True,
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
            log.warning(
                f"Metadata column '{key}' not found in adata.obs. Skipping filter."
            )
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
        log.warning(
            "No cells remaining after applying filters. Returning an empty AnnData object."
        )
        return AnnData()

    # The core slicing operation
    adata_subset = adata[combined_mask, :].copy()

    # AnnData slicing automatically handles .raw correctly. If we want to ensure
    # the .raw attribute uses the original var, we can explicitly re-assign it.
    if keep_raw_genes and adata.raw is not None:
        # Create a new raw object from the original raw data, but with subsetted cells
        adata_subset.raw = adata.raw[adata_subset.obs_names, :].copy()
        log.info(
            f"Subset .raw created, retaining all {adata.raw.n_vars} original genes."
        )

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
        raise ValueError(
            f"Columns {missing_cols} not found in the source AnnData object's .obs"
        )

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
    log.info(
        f"Subsetting target object based on new annotations with filters: {filters}"
    )

    # Now we can call the original, simple subset_adata function
    adata_subset = subset_adata(temp_adata, filters=filters)

    return adata_subset


def merge_obs_metadata(
    adata: AnnData,
    metadata_path: str,
    left_on: Optional[str] = None, # If None, uses adata.obs.index
    right_on: Optional[str] = None, # If None, uses metadata_df.index
    how: str = "left",
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
    
    # Perform the merge
    if left_on is None: # Join on index
        adata.obs = adata.obs.join(meta_df.set_index(right_on) if right_on else meta_df, how=how)
    else: # Join on a column
        adata.obs = adata.obs.merge(meta_df, left_on=left_on, right_on=right_on, how=how, suffixes=("", "_new"))

    new_cols = set(adata.obs.columns) - initial_cols
    log.info(f"Successfully merged {len(new_cols)} new columns into .obs: {list(new_cols)}")
    
    return adata