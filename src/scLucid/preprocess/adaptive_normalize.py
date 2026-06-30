"""
Adaptive Quality-Aware Normalization for scRNA-seq data.

This module implements novel normalization strategies that adapt to:
- Cell quality (based on QC metrics)
- Cell type-specific RNA content
- Technical batch effects

Key innovations:
1. Quality-stratified normalization
2. Cell-type aware size factor estimation
3. Robust outlier handling
"""

import logging
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Literal, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse
import seaborn as sns
from anndata import AnnData
from sklearn.preprocessing import QuantileTransformer

from .utils import apply_log1p, apply_row_scale, resolve_input_matrix

log = logging.getLogger(__name__)

__all__ = [
    "AdaptiveNormalizationConfig",
    "adaptive_normalize",
    "estimate_cell_size_factors",
    "quality_aware_normalize",
]


# =============================================================================
# Configuration
# =============================================================================


@dataclass
class AdaptiveNormalizationConfig:
    """Configuration for adaptive normalization."""

    # === Basic settings ===
    method: Literal[
        "quality_aware",  # Quality-stratified normalization
        "deconvolution_pool",  # Python-only pooled size factor estimation
        "quantile_regression",  # Quantile regression normalization
    ] = "quality_aware"

    input_layer: str = "counts"
    output_layer: str = "adaptive_normalized"

    # === Quality-aware settings ===
    quality_metrics: List[str] = field(
        default_factory=lambda: ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    )
    n_quality_bins: int = 5  # Stratify cells into N quality bins
    use_quality_weights: bool = True  # Weight cells by quality in normalization

    # === pooled deconvolution settings ===
    pool_size: int = 100
    min_mean: float = 0.1

    # === Quantile regression settings ===
    quantile: float = 0.75  # Which quantile to use as reference
    n_quantile_bins: int = 1000  # Number of bins for quantile transformation

    # === General settings ===
    target_sum: Optional[float] = None  # If None, use median total counts
    log_transform: bool = True
    clip_values: Optional[Tuple[float, float]] = None  # (min, max) for clipping

    # === Diagnostics ===
    plot: bool = True
    save_dir: Optional[str] = None

    def __post_init__(self):
        """Validate configuration."""
        if self.n_quality_bins < 2:
            raise ValueError("n_quality_bins must be >= 2")

        if not 0 < self.quantile < 1:
            raise ValueError("quantile must be in (0, 1)")


# =============================================================================
# Core Functions
# =============================================================================


def estimate_cell_size_factors(
    adata: AnnData,
    method: Literal["median_ratio", "deconvolution"] = "median_ratio",
    layer: Optional[str] = None,
    min_mean: float = 0.1,
    pool_size: int = 100,
) -> np.ndarray:
    """
    Estimate cell-specific size factors for normalization.

    This is more robust than simple total count normalization because:
    1. It uses only stably expressed genes (median-ratio method)
    2. It can pool similar cells with a Python-only deconvolution approximation
    3. It's robust to compositional effects

    Args:
        adata: AnnData object
        method: Size factor estimation method
        layer: Layer to use (if None, use adata.X)
        min_mean: Minimum mean expression for gene inclusion
        pool_size: Pool size for deconvolution method

    Returns:
        Array of size factors (one per cell)
    """
    # Get expression matrix
    if layer is not None:
        X = adata.layers[layer]
    else:
        X = adata.X

    if scipy.sparse.issparse(X):
        X_dense = X.toarray()
    else:
        X_dense = X

    n_cells, n_genes = X_dense.shape

    if method == "median_ratio":
        # DESeq2-style median ratio method
        log.info("Estimating size factors using median-ratio method...")

        # 1. Compute geometric mean per gene (excluding zeros)
        gene_means = np.zeros(n_genes)
        for j in range(n_genes):
            gene_expr = X_dense[:, j]
            nonzero = gene_expr > 0
            if nonzero.sum() > 0:
                gene_means[j] = np.exp(np.mean(np.log(gene_expr[nonzero])))
            else:
                gene_means[j] = 0

        # 2. Select genes with sufficient expression
        valid_genes = gene_means > min_mean

        if valid_genes.sum() < 100:
            log.warning(
                f"Only {valid_genes.sum()} genes with mean > {min_mean}. "
                "Consider lowering min_mean."
            )

        log.info(f"Using {valid_genes.sum()} genes for size factor estimation")

        # 3. Compute size factors
        size_factors = np.zeros(n_cells)
        for i in range(n_cells):
            ratios = X_dense[i, valid_genes] / gene_means[valid_genes]
            ratios = ratios[ratios > 0]  # Remove zeros
            if len(ratios) > 0:
                size_factors[i] = np.median(ratios)
            else:
                size_factors[i] = 1.0

    elif method == "deconvolution":
        # Simplified deconvolution without R dependency
        log.info("Estimating size factors using simplified total-count pooled deconvolution heuristic...")

        # 1. Create cell pools based on similarity
        from sklearn.cluster import KMeans

        # Use total counts for quick clustering
        total_counts = X_dense.sum(axis=1)
        log_counts = np.log1p(total_counts).reshape(-1, 1)

        n_pools = max(5, n_cells // pool_size)
        kmeans = KMeans(n_clusters=n_pools, random_state=42)
        clusters = kmeans.fit_predict(log_counts)

        # 2. Compute pool size factors
        pool_size_factors = np.zeros(n_pools)
        for pool_id in range(n_pools):
            pool_mask = clusters == pool_id
            pool_counts = X_dense[pool_mask, :].sum(axis=0)

            # Normalize by median
            pool_size_factors[pool_id] = np.median(pool_counts[pool_counts > 0])

        # 3. Deconvolve to get cell-level size factors
        size_factors = np.ones(n_cells)
        for pool_id in range(n_pools):
            pool_mask = clusters == pool_id
            pool_cells = np.where(pool_mask)[0]

            # Simple averaging within pool
            for cell_idx in pool_cells:
                cell_counts = X_dense[cell_idx, :]
                nonzero = cell_counts > 0
                if nonzero.sum() > 0:
                    size_factors[cell_idx] = np.median(
                        cell_counts[nonzero] / pool_size_factors[pool_id]
                    )

    else:
        raise ValueError(f"Unknown method: {method}")

    # Normalize size factors to have median = 1
    size_factors = size_factors / np.median(size_factors)

    # Handle extreme values
    size_factors = np.clip(size_factors, 0.1, 10.0)

    log.info(
        f"Size factors - median: {np.median(size_factors):.3f}, "
        f"range: [{size_factors.min():.3f}, {size_factors.max():.3f}]"
    )

    return size_factors


def quality_aware_normalize(
    adata: AnnData,
    quality_metrics: List[str],
    n_bins: int = 5,
    input_layer: str = "counts",
    output_layer: str = "quality_normalized",
    target_sum: Optional[float] = None,
    log_transform: bool = True,
) -> AnnData:
    """
    Quality-stratified normalization.

    Innovation: Different cells are normalized with different strategies
    based on their quality metrics. This prevents low-quality cells from
    distorting the normalization of high-quality cells.

    Algorithm:
    1. Compute composite quality score from multiple metrics
    2. Stratify cells into quality bins
    3. Normalize within each bin using bin-specific size factors
    4. Optionally weight cells by quality in downstream analysis

    Args:
        adata: AnnData object
        quality_metrics: List of QC metric columns in adata.obs
        n_bins: Number of quality bins
        input_layer: Input layer name
        output_layer: Output layer name
        target_sum: Target sum for normalization (if None, use median)
        log_transform: Whether to log-transform

    Returns:
        AnnData with quality-aware normalized data
    """
    warnings.warn(
        "quality_aware_normalize is a low-level algorithm entrypoint; prefer "
        "normalize_data(method='quality_aware', ...) for the public normalization API.",
        FutureWarning,
        stacklevel=2,
    )
    log.info("=" * 60)
    log.info("Quality-Aware Normalization")
    log.info("=" * 60)

    # Validate quality metrics
    missing_metrics = [m for m in quality_metrics if m not in adata.obs.columns]
    if missing_metrics:
        raise ValueError(f"Missing quality metrics in adata.obs: {missing_metrics}")

    # Get expression matrix
    X, _ = resolve_input_matrix(adata, input_layer)
    X = X.copy()

    is_sparse = scipy.sparse.issparse(X)
    if is_sparse:
        X = X.tocsr(copy=True)

    # === 1. Compute composite quality score ===
    log.info(f"Computing composite quality score from: {', '.join(quality_metrics)}")

    quality_scores = np.zeros(adata.n_obs)

    for metric in quality_metrics:
        values = adata.obs[metric].values.astype(np.float64)

        # Use nanmin/nanmax so a single NaN doesn't poison all cells
        vmin = np.nanmin(values)
        vmax = np.nanmax(values)
        denom = (vmax - vmin) + 1e-8

        # Normalize to [0, 1]
        # For metrics like pct_counts_mt: lower is better
        # For metrics like n_genes: higher is better
        if "mt" in metric.lower() or "pct" in metric.lower():
            # Lower is better
            normalized = 1.0 - (values - vmin) / denom
        else:
            # Higher is better
            normalized = (values - vmin) / denom

        # Replace NaN cells (e.g. original NaN) with mid-quality score
        nan_mask = np.isnan(normalized)
        if nan_mask.any():
            normalized[nan_mask] = 0.5

        quality_scores += normalized

    # Average across metrics
    quality_scores /= len(quality_metrics)

    # Store in adata
    adata.obs["quality_score"] = quality_scores

    # === 2. Stratify into quality bins ===
    quality_span = quality_scores.max() - quality_scores.min()
    if quality_span < 1e-10:
        # All quality scores are identical — single bin
        quality_bins = pd.Categorical(
            ["Q1"] * adata.n_obs, categories=[f"Q{i+1}" for i in range(n_bins)], ordered=True
        )
    else:
        # qcut may reduce the effective bin count when duplicates="drop" kicks in,
        # so generate labels after the fact to stay in sync.
        quality_bins = pd.qcut(
            quality_scores, q=n_bins, labels=False, duplicates="drop"
        )
        n_actual = quality_bins.max() + 1  # bins are 0-indexed
        bin_labels = {i: f"Q{i+1}" for i in range(n_bins)}
        quality_bins = pd.Categorical(
            [bin_labels.get(b, f"Q{b+1}") for b in quality_bins],
            categories=[f"Q{i+1}" for i in range(n_actual)],
            ordered=True,
        )
    adata.obs["quality_bin"] = quality_bins

    log.info(f"Stratified cells into {len(quality_bins.categories)} quality bins:")
    for bin_name in quality_bins.categories:
        n_cells = (quality_bins == bin_name).sum()
        log.info(f"  {bin_name}: {n_cells} cells")

    # === 3. Normalize within each bin ===
    row_scale = np.ones(adata.n_obs, dtype=np.float64)

    for bin_name in quality_bins.categories:
        bin_mask = quality_bins == bin_name
        bin_indices = np.where(bin_mask)[0]

        if len(bin_indices) == 0:
            continue

        # Get cells in this bin
        X_bin = X[bin_indices, :]

        # Compute bin-specific size factors
        total_counts_bin = np.asarray(X_bin.sum(axis=1)).ravel()

        if target_sum is None:
            # Use median of this bin
            target_bin = np.median(total_counts_bin)
        else:
            target_bin = target_sum
        if target_bin <= 0:
            target_bin = 1.0

        # Normalize
        size_factors_bin = total_counts_bin / target_bin
        size_factors_bin[size_factors_bin == 0] = 1
        row_scale[bin_indices] = 1.0 / size_factors_bin

        log.info(f"  {bin_name}: target_sum={target_bin:.0f}")

    X_normalized = apply_row_scale(X, row_scale)

    # === 4. Log transform ===
    if log_transform:
        X_normalized = apply_log1p(X_normalized)

    # === 5. Store result ===
    adata.layers[output_layer] = X_normalized

    # === 6. Compute quality weights ===
    # Cells with higher quality get higher weights in downstream analysis
    # This is stored but not automatically applied
    adata.obs["quality_weight"] = quality_scores
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})[
        "quality_aware_normalization_policy"
    ] = {
        "schema_version": "quality_aware_normalization_policy_v1",
        "model_type": "quality_stratified_library_size_heuristic",
        "claim_level": "heuristic_preprocessing",
        "quality_metrics": list(quality_metrics),
        "n_quality_bins": int(n_bins),
        "target_sum": target_sum,
        "log_transform": bool(log_transform),
        "review_note": (
            "This method scales cells within quality-score bins using total counts. "
            "It stores quality_weight for downstream review but does not formally "
            "correct systematic low-quality-cell bias."
        ),
    }

    log.info(f"Quality-aware normalization complete. Stored in layer '{output_layer}'")
    log.info("=" * 60)

    return adata


def adaptive_normalize(
    adata: AnnData,
    config: Optional[AdaptiveNormalizationConfig] = None,
    **kwargs,
) -> AnnData:
    """
    Main entry point for adaptive normalization.

    This function dispatches to specific normalization methods based on config.

    Args:
        adata: AnnData object
        config: Configuration object
        **kwargs: Override config parameters

    Returns:
        AnnData with normalized data
    """
    warnings.warn(
        "adaptive_normalize is a specialized dispatcher; prefer "
        "normalize_data(method=..., ...) for the public normalization API.",
        FutureWarning,
        stacklevel=2,
    )
    # Setup config
    if config is None:
        config = AdaptiveNormalizationConfig()

    # Apply kwargs overrides
    for key, value in kwargs.items():
        if hasattr(config, key):
            setattr(config, key, value)

    log.info(f"Running adaptive normalization with method: {config.method}")

    # === Dispatch to method ===
    if config.method == "quality_aware":
        adata = quality_aware_normalize(
            adata,
            quality_metrics=config.quality_metrics,
            n_bins=config.n_quality_bins,
            input_layer=config.input_layer,
            output_layer=config.output_layer,
            target_sum=config.target_sum,
            log_transform=config.log_transform,
        )
        method_policy = adata.uns.get("sclucid", {}).get("preprocess", {}).get(
            "quality_aware_normalization_policy",
            {},
        )

    elif config.method == "deconvolution_pool":
        # Use Python-only pooled size factors
        size_factors = estimate_cell_size_factors(
            adata,
            method="deconvolution",
            layer=config.input_layer,
            pool_size=config.pool_size,
            min_mean=config.min_mean,
        )

        # Apply size factors
        X, _ = resolve_input_matrix(adata, config.input_layer)
        X = X.copy()
        is_sparse_input = scipy.sparse.issparse(X)
        row_scale = 1.0 / size_factors
        X_normalized = apply_row_scale(X, row_scale)

        if config.log_transform:
            X_normalized = apply_log1p(X_normalized)

        adata.layers[config.output_layer] = X_normalized
        adata.obs["deconvolution_size_factors"] = size_factors
        method_policy = {
            "schema_version": "adaptive_normalization_policy_v1",
            "model_type": "total_count_pooled_deconvolution_heuristic",
            "claim_level": "heuristic_size_factor_estimate",
            "pooling_basis": "log_total_counts_kmeans",
            "pool_size": config.pool_size,
            "min_mean": config.min_mean,
            "review_note": (
                "This Python-only pooled size-factor estimate clusters cells by total counts, "
                "not cell-type composition; it should not be treated as scran-equivalent."
            ),
        }

    elif config.method == "quantile_regression":
        # Quantile normalization
        log.info(f"Applying quantile normalization (quantile={config.quantile})...")

        X, _ = resolve_input_matrix(adata, config.input_layer)
        X = X.copy()

        is_sparse_input = scipy.sparse.issparse(X)
        if is_sparse_input:
            n_obs, n_vars = X.shape
            est_dense_bytes = n_obs * n_vars * 8
            if est_dense_bytes > 8e9:
                log.warning(
                    f"Quantile regression requires dense matrix conversion. "
                    f"Estimated memory: {est_dense_bytes / 1e9:.1f} GB. "
                    f"Consider using method='quality_aware' or subsampling."
                )
            X = X.toarray()

        # Apply quantile transformation
        qt = QuantileTransformer(
            n_quantiles=config.n_quantile_bins, output_distribution="normal", random_state=42
        )

        X_normalized = qt.fit_transform(X)

        if is_sparse_input:
            X_normalized = scipy.sparse.csr_matrix(X_normalized)

        adata.layers[config.output_layer] = X_normalized
        method_policy = {
            "schema_version": "adaptive_normalization_policy_v1",
            "model_type": "quantile_transform_heuristic",
            "claim_level": "heuristic_distribution_transform",
            "review_note": (
                "Quantile transformation can alter biological distributional differences; "
                "use as exploratory preprocessing with downstream sensitivity checks."
            ),
        }

    else:
        raise ValueError(f"Unknown method: {config.method}")

    # === Store metadata ===
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["adaptive_normalization"] = {
        "method": config.method,
        "config": {
            "quality_metrics": config.quality_metrics,
            "n_quality_bins": config.n_quality_bins,
            "target_sum": config.target_sum,
            "log_transform": config.log_transform,
        },
        "output_layer": config.output_layer,
        "model_type": method_policy.get("model_type"),
        "claim_level": method_policy.get("claim_level"),
        "review_note": method_policy.get("review_note"),
        "method_policy": method_policy,
    }

    # === Generate diagnostic plots ===
    if config.plot and config.save_dir:
        _plot_normalization_diagnostics(adata, config)

    return adata


def _plot_normalization_diagnostics(
    adata: AnnData,
    config: AdaptiveNormalizationConfig,
):
    """Generate diagnostic plots for adaptive normalization."""
    save_dir = Path(config.save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # Get data
    X_raw, _ = resolve_input_matrix(adata, config.input_layer)

    X_norm = adata.layers[config.output_layer]

    if scipy.sparse.issparse(X_raw):
        X_raw = X_raw.toarray()
    if scipy.sparse.issparse(X_norm):
        X_norm = X_norm.toarray()

    # 1. Total counts before/after
    ax = axes[0, 0]
    total_before = X_raw.sum(axis=1)
    total_after = X_norm.sum(axis=1)

    ax.hist(total_before, bins=50, alpha=0.5, label="Before", color="blue")
    ax.hist(total_after, bins=50, alpha=0.5, label="After", color="red")
    ax.set_xlabel("Total Counts")
    ax.set_ylabel("Frequency")
    ax.set_title("Total Counts Distribution")
    ax.legend()
    ax.set_yscale("log")

    # 2. Quality score distribution
    if "quality_score" in adata.obs:
        ax = axes[0, 1]
        sns.histplot(adata.obs["quality_score"], bins=50, ax=ax, kde=True)
        ax.set_title("Quality Score Distribution")
        ax.set_xlabel("Quality Score")

    # 3. Quality bins
    if "quality_bin" in adata.obs:
        ax = axes[0, 2]
        adata.obs["quality_bin"].value_counts().plot(kind="bar", ax=ax)
        ax.set_title("Cells per Quality Bin")
        ax.set_ylabel("Number of Cells")
        ax.set_xlabel("Quality Bin")

    # 4. Size factors (if available)
    if "deconvolution_size_factors" in adata.obs:
        ax = axes[1, 0]
        sns.histplot(adata.obs["deconvolution_size_factors"], bins=50, ax=ax, kde=True)
        ax.set_title("Size Factors Distribution")
        ax.axvline(1.0, color="red", linestyle="--", label="Median")
        ax.legend()

    # 5. Gene expression before/after (top genes)
    ax = axes[1, 1]
    gene_means_before = X_raw.mean(axis=0)
    gene_means_after = X_norm.mean(axis=0)

    ax.scatter(gene_means_before, gene_means_after, alpha=0.3, s=1)
    ax.set_xlabel("Mean Expression (Before)")
    ax.set_ylabel("Mean Expression (After)")
    ax.set_title("Gene Expression Shift")
    ax.set_xscale("log")
    ax.set_yscale("log")

    # 6. Variance before/after
    ax = axes[1, 2]
    gene_vars_before = np.var(X_raw, axis=0)
    gene_vars_after = np.var(X_norm, axis=0)

    ax.scatter(gene_vars_before, gene_vars_after, alpha=0.3, s=1)
    ax.set_xlabel("Variance (Before)")
    ax.set_ylabel("Variance (After)")
    ax.set_title("Gene Variance Shift")
    ax.set_xscale("log")
    ax.set_yscale("log")

    plt.tight_layout()
    plt.savefig(save_dir / "adaptive_normalization_diagnostics.png", dpi=300)
    plt.close()

    log.info(f"Saved diagnostic plots to {save_dir}")
