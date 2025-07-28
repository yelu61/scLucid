"""
Utility functions for single-cell RNA-seq data analysis.

This module provides common helper functions that can be used
across different parts of the analysis pipeline.
"""

import logging
from contextlib import contextmanager
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import median_abs_deviation

log = logging.getLogger(__name__)

__all__ = ["use_layer_as_X", "identify_outliers"]


@contextmanager
def use_layer_as_X(adata: AnnData, layer: Optional[str]):
    """Context manager to temporarily use a layer as adata.X."""
    if layer is None or layer not in adata.layers:
        # If no layer, do nothing and yield
        yield
        return

    X_backup = adata.X.copy()
    adata.X = adata.layers[layer].copy()
    try:
        yield
    finally:
        # Always restore the original .X
        adata.X = X_backup


def identify_outliers(
    adata: AnnData,
    metrics: List[Tuple[str, str, Optional[float]]],
    sample_key: Optional[str] = None,
    nmads: float = 5.0,
) -> pd.Series:
    """
    Identify outliers based on one or more metrics using the median absolute deviation (MAD).

    This function can process multiple metrics and optionally group by sample.

    Args:
        adata: AnnData object to check for outliers.
        metrics: List of tuples (metric, direction, optional_threshold) for outlier detection.
                 Direction must be 'both', 'upper', or 'lower'.
                 If threshold is None, it will be calculated using MAD.
        sample_key: If provided, outliers will be identified separately per sample.
        nmads: Number of median absolute deviations for outlier detection.

    Returns:
        Boolean pd.Series indicating if a cell is an outlier for any of the metrics.
    """
    if not metrics:
        return pd.Series(False, index=adata.obs_names)

    # Initialize result
    outliers = pd.Series(False, index=adata.obs_names)

    # Process each sample separately if sample_key is provided
    if sample_key is not None and sample_key in adata.obs.columns:
        for sample in adata.obs[sample_key].unique():
            sample_mask = adata.obs[sample_key] == sample
            sample_adata = adata[sample_mask]

            if sample_adata.n_obs == 0:
                continue

            log.info(f"Identifying outliers for sample: {sample}")

            sample_outliers = pd.Series(False, index=sample_adata.obs_names)
            for metric, direction, threshold in metrics:
                if metric not in sample_adata.obs.columns:
                    log.warning(f"Metric '{metric}' not found in adata.obs, skipping")
                    continue

                values = sample_adata.obs[metric]

                if threshold is not None:
                    # Use fixed threshold
                    if direction == "upper":
                        metric_outliers = values > threshold
                    elif direction == "lower":
                        metric_outliers = values < threshold
                    elif direction == "both":
                        # Not implemented for fixed threshold
                        log.warning(
                            f"Direction 'both' not supported with fixed threshold for '{metric}', skipping"
                        )
                        continue
                    else:
                        log.warning(
                            f"Invalid direction '{direction}' for '{metric}', skipping"
                        )
                        continue
                else:
                    # Calculate threshold using MAD
                    median = np.nanmedian(values)
                    mad = median_abs_deviation(values, nan_policy="omit")

                    if mad == 0:
                        log.warning(
                            f"MAD is zero for '{metric}' in sample {sample}, skipping"
                        )
                        continue

                    if direction == "upper":
                        metric_outliers = values > (median + nmads * mad)
                    elif direction == "lower":
                        metric_outliers = values < (median - nmads * mad)
                    elif direction == "both":
                        upper_bound = median + nmads * mad
                        lower_bound = median - nmads * mad
                        metric_outliers = (values > upper_bound) | (
                            values < lower_bound
                        )
                    else:
                        log.warning(
                            f"Invalid direction '{direction}' for '{metric}', skipping"
                        )
                        continue

                outlier_count = metric_outliers.sum()
                log.info(
                    f"  Identified {outlier_count} outliers ({outlier_count / len(values):.2%}) for '{metric}' in sample {sample}"
                )

                sample_outliers |= metric_outliers

            # Update the global result
            outliers.loc[sample_adata.obs_names] = sample_outliers
    else:
        # Process entire dataset together
        for metric, direction, threshold in metrics:
            if metric not in adata.obs.columns:
                log.warning(f"Metric '{metric}' not found in adata.obs, skipping")
                continue

            values = adata.obs[metric]

            if threshold is not None:
                # Use fixed threshold
                if direction == "upper":
                    metric_outliers = values > threshold
                elif direction == "lower":
                    metric_outliers = values < threshold
                elif direction == "both":
                    log.warning(
                        f"Direction 'both' not supported with fixed threshold for '{metric}', skipping"
                    )
                    continue
                else:
                    log.warning(
                        f"Invalid direction '{direction}' for '{metric}', skipping"
                    )
                    continue
            else:
                # Calculate threshold using MAD
                median = np.nanmedian(values)
                mad = median_abs_deviation(values, nan_policy="omit")

                if mad == 0:
                    log.warning(f"MAD is zero for '{metric}', skipping")
                    continue

                if direction == "upper":
                    metric_outliers = values > (median + nmads * mad)
                elif direction == "lower":
                    metric_outliers = values < (median - nmads * mad)
                elif direction == "both":
                    upper_bound = median + nmads * mad
                    lower_bound = median - nmads * mad
                    metric_outliers = (values > upper_bound) | (values < lower_bound)
                else:
                    log.warning(
                        f"Invalid direction '{direction}' for '{metric}', skipping"
                    )
                    continue

            outlier_count = metric_outliers.sum()
            log.info(
                f"Identified {outlier_count} outliers ({outlier_count / len(values):.2%}) for '{metric}'"
            )

            outliers |= metric_outliers

    total_count = outliers.sum()
    log.info(
        f"Total outliers identified: {total_count} ({total_count / len(outliers):.2%})"
    )

    return outliers
