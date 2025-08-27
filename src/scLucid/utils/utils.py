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
from typing import Any, Dict, List, Optional, Tuple, Union

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import io
from scipy.stats import median_abs_deviation

log = logging.getLogger(__name__)

__all__ = ["load_10x_data", "use_layer_as_X", "identify_outliers"]


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
    adata = anndata.AnnData(X=X, 
                        obs=pd.DataFrame(index=barcodes.values), 
                        var=pd.DataFrame(index=gene_names.values))

    adata.var_names_make_unique()  # Ensure gene names are unique

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


def _identify_outliers_subset(
    obs_subset: pd.DataFrame,
    metrics: List[Tuple[str, str, Optional[float]]],
    nmads: float = 5.0,
    group_name: str = "global",
) -> pd.Series:
    """
    Internal helper function to identify outliers on a subset of data.
    """
    subset_outliers = pd.Series(False, index=obs_subset.index)

    for metric, direction, threshold in metrics:
        if metric not in obs_subset.columns:
            log.warning(
                f"Metric '{metric}' not found in data for group '{group_name}', skipping."
            )
            continue

        values = obs_subset[metric]
        metric_outliers = pd.Series(False, index=obs_subset.index)

        if threshold is not None:
            # Use fixed threshold
            if direction == "upper":
                metric_outliers = values > threshold
            elif direction == "lower":
                metric_outliers = values < threshold
            elif direction == "both":
                # For fixed threshold, 'both' is not meaningful.
                # A user should provide two separate tuples for upper and lower bounds.
                log.warning(
                    f"Direction 'both' with a fixed threshold is ambiguous for '{metric}'. "
                    "Please provide separate 'upper' and 'lower' tuples if needed. Skipping."
                )
                continue
            else:
                log.warning(
                    f"Invalid direction '{direction}' for '{metric}', skipping."
                )
                continue
        else:
            # Calculate threshold using MAD
            median = np.nanmedian(values)
            mad = median_abs_deviation(values, scale="normal", nan_policy="omit")

            if mad == 0:
                log.warning(
                    f"MAD is zero for '{metric}' in group '{group_name}'. "
                    "Cannot perform outlier detection for this metric."
                )
                continue

            upper_bound = median + nmads * mad
            lower_bound = median - nmads * mad

            if direction == "upper":
                metric_outliers = values > upper_bound
            elif direction == "lower":
                metric_outliers = values < lower_bound
            elif direction == "both":
                metric_outliers = (values > upper_bound) | (values < lower_bound)
            else:
                log.warning(
                    f"Invalid direction '{direction}' for '{metric}', skipping."
                )
                continue

        outlier_count = metric_outliers.sum()
        if outlier_count > 0:
            log.info(
                f"  - Group '{group_name}': Identified {outlier_count} outliers "
                f"({outlier_count / len(values):.2%}) for metric '{metric}' (direction: {direction})"
            )

        subset_outliers |= metric_outliers

    return subset_outliers


def identify_outliers(
    adata: AnnData,
    metrics: List[Tuple[str, str, Optional[float]]],
    sample_key: Optional[str] = None,
    nmads: float = 5.0,
) -> pd.Series:
    """
    Identify outliers based on metrics using median absolute deviation (MAD) or fixed thresholds.

    This function can process multiple metrics and optionally group by sample for per-group
    outlier detection.

    Args:
        adata: AnnData object to check for outliers.
        metrics: List of tuples for outlier detection. Each tuple is (metric, direction, threshold).
                 - metric (str): Column name in `adata.obs`.
                 - direction (str): 'upper', 'lower', or 'both'.
                 - threshold (float, optional): If provided, this fixed value is used as the threshold.
                   If None, the threshold is calculated dynamically using MAD.
        sample_key: If provided, outliers will be identified separately per sample group.
        nmads: Number of median absolute deviations for dynamic outlier detection.

    Returns:
        Boolean pd.Series indicating if a cell is an outlier for any of the specified metrics.
    """
    if not metrics:
        return pd.Series(False, index=adata.obs_names)

    final_outliers = pd.Series(False, index=adata.obs_names)

    if sample_key and sample_key in adata.obs.columns:
        log.info(f"Identifying outliers per group in '{sample_key}'...")
        for sample_id, group_df in adata.obs.groupby(sample_key):
            group_outliers = _identify_outliers_subset(
                group_df, metrics, nmads, group_name=str(sample_id)
            )
            final_outliers[group_outliers.index] = group_outliers
    else:
        log.info("Identifying outliers on the entire dataset...")
        global_outliers = _identify_outliers_subset(
            adata.obs, metrics, nmads, group_name="global"
        )
        final_outliers = global_outliers

    total_count = final_outliers.sum()
    log.info(
        f"Total unique outliers identified: {total_count} ({total_count / len(final_outliers):.2%})"
    )

    return final_outliers


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
