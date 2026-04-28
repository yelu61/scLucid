"""
Adaptive threshold learning for QC metrics.

This module provides machine learning-based approaches to automatically learn
optimal QC thresholds from data distributions.
"""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
from anndata import AnnData
from sklearn.cluster import DBSCAN
from sklearn.mixture import GaussianMixture

from ..runtime import effective_n_jobs

log = logging.getLogger(__name__)

#: Scale factor to convert MAD to approximate standard deviation for a normal distribution.
#: For normally distributed data, std ≈ MAD * 1.4826.
MAD_SCALE_FACTOR: float = 1.4826


def compute_mad_bounds(
    values: np.ndarray,
    nmads: float = 5.0,
    direction: str = "both",
) -> Tuple[float, float]:
    """Compute outlier bounds using Median Absolute Deviation (MAD).

    This is the canonical MAD-based outlier detection implementation used
    throughout scLucid. It replaces duplicated MAD logic in
    ``filtering.py`` and ``adaptive_threshold.py``.

    Parameters
    ----------
    values : np.ndarray
        Input metric values (may contain NaNs).
    nmads : float, default=5.0
        Number of MADs from the median to use as the bound.
    direction : {'upper', 'lower', 'both'}, default='both'
        Which bound(s) to compute.

    Returns:
    -------
    Tuple[float, float]
        ``(lower_bound, upper_bound)``. For ``direction='upper'`` the
        lower bound is ``-inf``; for ``direction='lower'`` the upper
        bound is ``inf``.
    """
    clean = values[~np.isnan(values)]
    if len(clean) == 0:
        return -np.inf, np.inf

    median = float(np.median(clean))
    mad = float(np.median(np.abs(clean - median)))

    if mad == 0:
        log.debug("MAD is zero; bounds collapse to the median.")
        scaled_mad = 0.0
    else:
        scaled_mad = mad * MAD_SCALE_FACTOR

    lower = median - nmads * scaled_mad
    upper = median + nmads * scaled_mad

    if direction == "upper":
        lower = -np.inf
    elif direction == "lower":
        upper = np.inf
    elif direction != "both":
        raise ValueError(f"direction must be 'upper', 'lower', or 'both', got {direction!r}")

    return lower, upper


class AdaptiveThresholdLearner:
    """
    Automatically learn optimal QC thresholds using statistical methods.

    Supports multiple learning strategies:
    - GMM-based: Learn mixture of distributions
    - MAD-based: Median absolute deviation
    - Percentile-based: Statistical percentiles
    - Kernel density: Non-parametric density estimation
    """

    def __init__(
        self,
        method: str = "gmm",
        min_quality_cells: float = 0.5,
        random_state: int = 42,
    ):
        """
        Initialize the adaptive threshold learner.

        Args:
            method: Learning method ('gmm', 'mad', 'percentile', 'kde', 'dbscan')
            min_quality_cells: Minimum fraction of cells to retain
            random_state: Random seed for reproducibility
        """
        self.method = method
        self.min_quality_cells = min_quality_cells
        self.random_state = random_state

        self._learned_thresholds = {}
        self._fitted_models = {}

    def learn_threshold(
        self,
        metric_values: np.ndarray,
        metric_name: str,
        direction: str = "upper",
    ) -> float:
        """
        Learn optimal threshold for a single QC metric.

        Args:
            metric_values: Array of metric values
            metric_name: Name of the metric
            direction: 'upper' (filter high values) or 'lower' (filter low values)

        Returns:
            Learned threshold value
        """
        # Remove NaN and infinite values
        clean_values = metric_values[~np.isnan(metric_values)]
        clean_values = clean_values[~np.isinf(clean_values)]

        if len(clean_values) == 0:
            log.warning(f"No valid values for metric {metric_name}")
            return np.nan

        if self.method == "gmm":
            threshold = self._learn_threshold_gmm(clean_values, direction)
        elif self.method == "mad":
            threshold = self._learn_threshold_mad(clean_values, direction)
        elif self.method == "percentile":
            threshold = self._learn_threshold_percentile(clean_values, direction)
        elif self.method == "kde":
            threshold = self._learn_threshold_kde(clean_values, direction)
        elif self.method == "dbscan":
            threshold = self._learn_threshold_dbscan(clean_values, direction)
        else:
            raise ValueError(f"Unknown method: {self.method}")

        # Apply minimum quality constraint
        threshold = self._apply_min_quality_constraint(clean_values, threshold, direction)

        self._learned_thresholds[metric_name] = threshold

        log.info(
            f"Learned threshold for {metric_name} ({direction}): {threshold:.4f} "
            f"using {self.method}"
        )

        return threshold

    def _learn_threshold_gmm(
        self,
        values: np.ndarray,
        direction: str,
        n_components: int = 2,
    ) -> float:
        """
        Learn threshold using Gaussian Mixture Model.

        Assumes data comes from mixture of quality and low-quality distributions.
        """
        # Reshape for sklearn
        X = values.reshape(-1, 1)

        # Fit GMM
        gmm = GaussianMixture(
            n_components=n_components,
            random_state=self.random_state,
            max_iter=100,
        )

        try:
            gmm.fit(X)
        except Exception as e:
            log.debug(f"GMM fitting failed: {e}, falling back to percentile")
            return self._learn_threshold_percentile(values, direction)

        # Get means and sort
        means = gmm.means_.flatten()
        sorted_idx = np.argsort(means)

        # For upper threshold (filter high values like MT%)
        # Use the boundary between the two components
        if direction == "upper":
            # Threshold between distributions
            if n_components == 2:
                # Use weighted average of means as threshold
                weights = gmm.weights_.flatten()
                threshold = np.sum(means * weights)
            else:
                # Use mean of two largest components
                top_two = sorted_idx[-2:]
                threshold = np.mean(means[top_two])
        else:
            # For lower threshold (filter low values like gene counts)
            if n_components == 2:
                threshold = np.mean(means)
            else:
                # Use mean of two smallest components
                bottom_two = sorted_idx[:2]
                threshold = np.mean(means[bottom_two])

        # Store model for potential use
        self._fitted_models["gmm"] = gmm

        return float(threshold)

    def _learn_threshold_mad(
        self,
        values: np.ndarray,
        direction: str,
        nmads: float = 5.0,
    ) -> float:
        """
        Learn threshold using Median Absolute Deviation.

        Delegates to the canonical ``compute_mad_bounds`` so that the
        same MAD logic is used everywhere in the QC module.
        """
        lower, upper = compute_mad_bounds(values, nmads=nmads, direction=direction)

        if direction == "upper":
            return upper
        else:
            return max(0.0, lower)

    def _learn_threshold_percentile(
        self,
        values: np.ndarray,
        direction: str,
    ) -> float:
        """
        Learn threshold using percentiles.

        Conservative approach based on distribution statistics.
        """
        if direction == "upper":
            # Use 95th percentile for upper threshold
            threshold = np.percentile(values, 95)
        else:
            # Use 5th percentile for lower threshold
            threshold = np.percentile(values, 5)

        return float(threshold)

    def _learn_threshold_kde(
        self,
        values: np.ndarray,
        direction: str,
    ) -> float:
        """
        Learn threshold using Kernel Density Estimation.

        Finds local minima in density as threshold boundaries.
        """
        try:
            from scipy.stats import gaussian_kde

            # Fit KDE
            kde = gaussian_kde(values)
            x_range = np.linspace(values.min(), values.max(), 1000)
            density = kde(x_range)

            # Find local minima
            from scipy.signal import find_peaks

            # Invert density to find minima
            minima_indices, _ = find_peaks(-density, distance=20)

            if len(minima_indices) == 0:
                # No clear minima, fall back to percentile
                return self._learn_threshold_percentile(values, direction)

            if direction == "upper":
                # Use rightmost local minimum
                threshold_idx = minima_indices[-1]
            else:
                # Use leftmost local minimum
                threshold_idx = minima_indices[0]

            threshold = x_range[threshold_idx]

            return float(threshold)

        except Exception as e:
            log.debug(f"KDE failed: {e}, falling back to percentile")
            return self._learn_threshold_percentile(values, direction)

    def _learn_threshold_dbscan(
        self,
        values: np.ndarray,
        direction: str,
    ) -> float:
        """
        Learn threshold using DBSCAN clustering.

        Identifies outliers as low-quality cells.
        """
        try:
            # Reshape for DBSCAN
            X = values.reshape(-1, 1)

            # Run DBSCAN
            dbscan = DBSCAN(eps=values.std() * 0.5, min_samples=max(5, len(values) // 20))
            labels = dbscan.fit_predict(X)

            # Find outliers (label = -1)
            outlier_values = values[labels == -1]

            if len(outlier_values) == 0:
                # No outliers detected, use percentile
                return self._learn_threshold_percentile(values, direction)

            if direction == "upper":
                # Threshold is minimum of upper outliers
                threshold = np.min(outlier_values)
            else:
                # Threshold is maximum of lower outliers
                threshold = np.max(outlier_values)

            return float(threshold)

        except Exception as e:
            log.debug(f"DBSCAN failed: {e}, falling back to percentile")
            return self._learn_threshold_percentile(values, direction)

    def _apply_min_quality_constraint(
        self,
        values: np.ndarray,
        threshold: float,
        direction: str,
    ) -> float:
        """
        Apply minimum quality constraint to threshold.

        Ensures that at least min_quality_cells fraction passes QC.
        """
        n_cells = len(values)
        min_cells_to_keep = int(n_cells * self.min_quality_cells)

        if direction == "upper":
            # At most this fraction can fail
            max_failures = n_cells - min_cells_to_keep

            # Count cells that would fail
            n_failures = np.sum(values > threshold)

            if n_failures > max_failures:
                # Adjust threshold to keep minimum cells
                threshold = np.sort(values)[-max_failures]
        else:
            # At most this fraction can fail
            max_failures = n_cells - min_cells_to_keep

            # Count cells that would fail
            n_failures = np.sum(values < threshold)

            if n_failures > max_failures:
                # Adjust threshold to keep minimum cells
                threshold = np.sort(values)[max_failures]

        return float(threshold)

    def learn_all_thresholds(
        self,
        adata: AnnData,
        metrics: Optional[Dict[str, str]] = None,
    ) -> Dict[str, float]:
        """
        Learn thresholds for all specified QC metrics.

        Args:
            adata: AnnData object with QC metrics
            metrics: Dictionary of {metric_name: direction} pairs
                    If None, uses default metrics

        Returns:
            Dictionary of learned thresholds
        """
        if metrics is None:
            metrics = {
                "log1p_n_genes_by_counts": "lower",
                "log1p_total_counts": "lower",
                "pct_counts_mt": "upper",
                "pct_counts_in_top_20_genes": "upper",
            }

        learned_thresholds = {}

        for metric_name, direction in metrics.items():
            if metric_name not in adata.obs:
                log.warning(f"Metric {metric_name} not found in adata.obs")
                continue

            values = adata.obs[metric_name].values

            try:
                threshold = self.learn_threshold(values, metric_name, direction)
                learned_thresholds[metric_name] = threshold
            except Exception as e:
                log.error(f"Failed to learn threshold for {metric_name}: {e}")
                continue

        return learned_thresholds

    def predict_quality(
        self,
        metric_values: np.ndarray,
        metric_name: str,
        direction: str,
    ) -> np.ndarray:
        """
        Predict quality labels for cells based on learned threshold.

        Args:
            metric_values: Metric values for cells
            metric_name: Name of the metric
            direction: Filter direction

        Returns:
            Boolean array (True = high quality, False = low quality)
        """
        if metric_name not in self._learned_thresholds:
            raise ValueError(f"No learned threshold for {metric_name}")

        threshold = self._learned_thresholds[metric_name]

        if direction == "upper":
            quality = metric_values <= threshold
        else:
            quality = metric_values >= threshold

        return quality


class MultiMetricAdaptiveLearner:
    """
    Adaptive threshold learner that considers multiple metrics jointly.

    Uses multivariate approaches to find optimal threshold combinations.
    """

    def __init__(
        self,
        method: str = "isolation_forest",
        contamination: float = 0.1,
        random_state: int = 42,
    ):
        """
        Initialize multi-metric learner.

        Args:
            method: Method ('isolation_forest', 'local_outlier_factor', 'one_class_svm')
            contamination: Expected fraction of outliers
            random_state: Random seed
        """
        self.method = method
        self.contamination = contamination
        self.random_state = random_state

        self._model = None

    def fit(
        self,
        adata: AnnData,
        metrics: List[str],
    ):
        """
        Fit the multi-metric outlier detection model.

        Args:
            adata: AnnData object with QC metrics
            metrics: List of metric names to use
        """
        # Prepare data matrix
        X = np.column_stack([adata.obs[m].values for m in metrics])

        # Handle missing values
        X = np.nan_to_num(X, nan=0.0)

        if self.method == "isolation_forest":
            from sklearn.ensemble import IsolationForest

            self._model = IsolationForest(
                contamination=self.contamination,
                random_state=self.random_state,
                n_jobs=effective_n_jobs(-1),
            )
            self._model.fit(X)

        elif self.method == "local_outlier_factor":
            from sklearn.neighbors import LocalOutlierFactor

            self._model = LocalOutlierFactor(
                contamination=self.contamination,
                n_neighbors=20,
                n_jobs=effective_n_jobs(-1),
            )
            self._model.fit(X)

        elif self.method == "one_class_svm":
            from sklearn.svm import OneClassSVM

            self._model = OneClassSVM(
                nu=self.contamination,
                kernel="rbf",
            )
            self._model.fit(X)

        else:
            raise ValueError(f"Unknown method: {self.method}")

    def predict(self, adata: AnnData, metrics: List[str]) -> np.ndarray:
        """
        Predict quality labels using fitted model.

        Args:
            adata: AnnData object
            metrics: List of metric names

        Returns:
            Boolean array (True = high quality, False = outlier/low quality)
        """
        if self._model is None:
            raise ValueError("Model not fitted. Call fit() first.")

        # Prepare data
        X = np.column_stack([adata.obs[m].values for m in metrics])
        X = np.nan_to_num(X, nan=0.0)

        # Get predictions
        if self.method == "local_outlier_factor":
            # LOF has fit_predict instead of predict
            predictions = self._model.fit_predict(X)
        else:
            predictions = self._model.predict(X)

        # Convert to boolean (1 = inlier/high quality, -1 = outlier)
        quality = predictions == 1

        return quality

    def fit_predict(
        self,
        adata: AnnData,
        metrics: List[str],
    ) -> np.ndarray:
        """
        Fit model and return predictions.

        Args:
            adata: AnnData object
            metrics: List of metric names

        Returns:
            Boolean array (True = high quality)
        """
        self.fit(adata, metrics)
        return self.predict(adata, metrics)
