"""
Utility functions for single-cell RNA-seq data analysis.

This module provides common helper functions that can be used
across different parts of the analysis pipeline.
"""

import gc
import logging
import os
from contextlib import contextmanager
from typing import List, Optional, Tuple

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy.stats import median_abs_deviation

log = logging.getLogger(__name__)

__all__ = ["load_10x_data", "use_layer_as_X", "identify_outliers"]


def load_10x_data(
    base_dir: str,
    samples: List[str],
    metadata_file: Optional[str] = None,
    sample_id_col: Optional[str] = None,
) -> AnnData:
    """
    Load multiple 10x Genomics samples and combine them into a single AnnData object.

    Args:
        base_dir: The base directory containing the sample folders.
        samples: A list of sample folder names.
        metadata_file: Path to a CSV/TSV file containing metadata for each sample.
        sample_id_col: The column in the metadata file that matches the sample names.

    Returns:
        A combined AnnData object.
    """
    adata_list = []

    metadata_df = None
    if metadata_file:
        log.info(f"Loading metadata from {metadata_file}")
        metadata_df = pd.read_csv(
            metadata_file, sep=None, engine="python"
        )  # Auto-detect separator
        if sample_id_col is None or sample_id_col not in metadata_df.columns:
            raise ValueError(
                "`sample_id_col` must be a valid column in the metadata file."
            )
        metadata_df.set_index(sample_id_col, inplace=True)

    for sample in samples:
        sample_path = os.path.join(
            base_dir, sample, "outs", "filtered_feature_bc_matrix"
        )

        if not os.path.isdir(sample_path):
            log.warning(f"Directory not found, skipping: {sample_path}")
            continue

        try:
            adata = sc.read_10x_mtx(sample_path, var_names="gene_symbols", cache=True)
            adata.obs["sampleID"] = sample

            # Add metadata if available
            if metadata_df is not None and sample in metadata_df.index:
                for col in metadata_df.columns:
                    adata.obs[col] = metadata_df.loc[sample, col]

            adata.obs_names = [f"{sample}_{bc}" for bc in adata.obs_names]
            log.info(f"Loaded {sample}: {adata.n_obs} cells, {adata.n_vars} genes")
            adata_list.append(adata)
            gc.collect()

        except Exception as e:
            log.error(f"Error processing {sample}: {e}")

    if not adata_list:
        log.error("No samples were loaded successfully.")
        return AnnData()

    log.info(f"Merging {len(adata_list)} samples...")
    combined_adata = anndata.concat(adata_list, join="outer", label="sample_source")
    log.info(
        f"Combined dataset: {combined_adata.n_obs} cells, {combined_adata.n_vars} genes"
    )

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
