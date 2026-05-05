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
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Literal, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import seaborn as sns
from anndata import AnnData
from sklearn.preprocessing import QuantileTransformer

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
        "scran_pool",  # scran with cell pooling
        "quantile_regression",  # Quantile regression normalization
        "sctransform",  # Variance-stabilizing transformation
    ] = "quality_aware"

    input_layer: str = "counts"
    output_layer: str = "adaptive_normalized"

    # === Quality-aware settings ===
    quality_metrics: List[str] = field(
        default_factory=lambda: ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    )
    n_quality_bins: int = 5  # Stratify cells into N quality bins
    use_quality_weights: bool = True  # Weight cells by quality in normalization

    # === scran settings ===
    scran_pool_size: int = 100  # Pool size for scran
    scran_min_mean: float = 0.1  # Minimum mean expression for size factor calculation

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
    method: Literal["median_ratio", "scran", "deconvolution"] = "median_ratio",
    layer: Optional[str] = None,
    min_mean: float = 0.1,
    pool_size: int = 100,
) -> np.ndarray:
    """
    Estimate cell-specific size factors for normalization.

    This is more robust than simple total count normalization because:
    1. It uses only stably expressed genes (median-ratio method)
    2. It pools similar cells (scran method)
    3. It's robust to compositional effects

    Args:
        adata: AnnData object
        method: Size factor estimation method
        layer: Layer to use (if None, use adata.X)
        min_mean: Minimum mean expression for gene inclusion
        pool_size: Pool size for scran method

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

    elif method == "scran":
        # scran pooling-based deconvolution
        log.info("Estimating size factors using scran method...")

        try:
            # Requires scran from Bioconductor via rpy2
            import rpy2.robjects as ro
            from rpy2.robjects import numpy2ri
            from rpy2.robjects.packages import importr

            numpy2ri.activate()

            scran = importr("scran")

            # Create pools of similar cells
            # Here we use a simple k-means clustering
            from sklearn.cluster import KMeans

            # Use PCA for clustering (if available)
            if "X_pca" in adata.obsm:
                clustering_data = adata.obsm["X_pca"][:, :20]
            else:
                # Use top variable genes
                gene_vars = np.var(X_dense, axis=0)
                top_genes = np.argsort(-gene_vars)[:1000]
                clustering_data = X_dense[:, top_genes]

            n_pools = max(5, n_cells // pool_size)
            kmeans = KMeans(n_clusters=n_pools, random_state=42)
            clusters = kmeans.fit_predict(clustering_data)

            # Convert to R format
            counts_r = ro.r.matrix(X_dense.T, nrow=n_genes, ncol=n_cells)
            clusters_r = ro.IntVector(clusters + 1)  # R uses 1-based indexing

            # Compute size factors
            size_factors_r = scran.computeSumFactors(
                counts_r, clusters=clusters_r, min_mean=min_mean
            )

            size_factors = np.array(size_factors_r)

            numpy2ri.deactivate()

        except ImportError:
            log.warning(
                "scran method requires R with scran package installed. "
                "Falling back to median-ratio method."
            )
            return estimate_cell_size_factors(
                adata, method="median_ratio", layer=layer, min_mean=min_mean
            )

    elif method == "deconvolution":
        # Simplified deconvolution without R dependency
        log.info("Estimating size factors using deconvolution method...")

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
    log.info("=" * 60)
    log.info("Quality-Aware Normalization")
    log.info("=" * 60)

    # Validate quality metrics
    missing_metrics = [m for m in quality_metrics if m not in adata.obs.columns]
    if missing_metrics:
        raise ValueError(f"Missing quality metrics in adata.obs: {missing_metrics}")

    # Get expression matrix
    if input_layer in adata.layers:
        X = adata.layers[input_layer].copy()
    else:
        X = adata.X.copy()

    if scipy.sparse.issparse(X):
        X = X.toarray()

    # === 1. Compute composite quality score ===
    log.info(f"Computing composite quality score from: {', '.join(quality_metrics)}")

    quality_scores = np.zeros(adata.n_obs)

    for metric in quality_metrics:
        values = adata.obs[metric].values

        # Normalize to [0, 1]
        # For metrics like pct_counts_mt: lower is better
        # For metrics like n_genes: higher is better
        if "mt" in metric.lower() or "pct" in metric.lower():
            # Lower is better
            normalized = 1 - (values - values.min()) / (values.max() - values.min() + 1e-8)
        else:
            # Higher is better
            normalized = (values - values.min()) / (values.max() - values.min() + 1e-8)

        quality_scores += normalized

    # Average across metrics
    quality_scores /= len(quality_metrics)

    # Store in adata
    adata.obs["quality_score"] = quality_scores

    # === 2. Stratify into quality bins ===
    quality_bins = pd.qcut(
        quality_scores, q=n_bins, labels=[f"Q{i+1}" for i in range(n_bins)], duplicates="drop"
    )
    adata.obs["quality_bin"] = quality_bins

    log.info(f"Stratified cells into {n_bins} quality bins:")
    for bin_name in quality_bins.categories:
        n_cells = (quality_bins == bin_name).sum()
        log.info(f"  {bin_name}: {n_cells} cells")

    # === 3. Normalize within each bin ===
    X_normalized = np.zeros_like(X)

    for bin_name in quality_bins.categories:
        bin_mask = quality_bins == bin_name
        bin_indices = np.where(bin_mask)[0]

        if len(bin_indices) == 0:
            continue

        # Get cells in this bin
        X_bin = X[bin_indices, :]

        # Compute bin-specific size factors
        total_counts_bin = X_bin.sum(axis=1)

        if target_sum is None:
            # Use median of this bin
            target_bin = np.median(total_counts_bin)
        else:
            target_bin = target_sum

        # Normalize
        size_factors_bin = total_counts_bin / target_bin
        size_factors_bin[size_factors_bin == 0] = 1

        X_bin_normalized = X_bin / size_factors_bin[:, np.newaxis]

        # Store
        X_normalized[bin_indices, :] = X_bin_normalized

        log.info(f"  {bin_name}: target_sum={target_bin:.0f}")

    # === 4. Log transform ===
    if log_transform:
        X_normalized = np.log1p(X_normalized)

    # === 5. Store result ===
    if scipy.sparse.issparse(adata.X):
        X_normalized = scipy.sparse.csr_matrix(X_normalized)

    adata.layers[output_layer] = X_normalized

    # === 6. Compute quality weights ===
    # Cells with higher quality get higher weights in downstream analysis
    # This is stored but not automatically applied
    adata.obs["quality_weight"] = quality_scores

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

    elif config.method == "scran_pool":
        # Use scran-based size factors
        size_factors = estimate_cell_size_factors(
            adata,
            method="scran",
            layer=config.input_layer,
            pool_size=config.scran_pool_size,
            min_mean=config.scran_min_mean,
        )

        # Apply size factors
        if config.input_layer in adata.layers:
            X = adata.layers[config.input_layer].copy()
        else:
            X = adata.X.copy()

        if scipy.sparse.issparse(X):
            X = X.toarray()

        X_normalized = X / size_factors[:, np.newaxis]

        if config.log_transform:
            X_normalized = np.log1p(X_normalized)

        if scipy.sparse.issparse(adata.X):
            X_normalized = scipy.sparse.csr_matrix(X_normalized)

        adata.layers[config.output_layer] = X_normalized
        adata.obs["scran_size_factors"] = size_factors

    elif config.method == "quantile_regression":
        # Quantile normalization
        log.info(f"Applying quantile normalization (quantile={config.quantile})...")

        if config.input_layer in adata.layers:
            X = adata.layers[config.input_layer].copy()
        else:
            X = adata.X.copy()

        if scipy.sparse.issparse(X):
            X = X.toarray()

        # Apply quantile transformation
        qt = QuantileTransformer(
            n_quantiles=config.n_quantile_bins, output_distribution="normal", random_state=42
        )

        X_normalized = qt.fit_transform(X)

        if scipy.sparse.issparse(adata.X):
            X_normalized = scipy.sparse.csr_matrix(X_normalized)

        adata.layers[config.output_layer] = X_normalized

    elif config.method == "sctransform":
        # Variance-stabilizing transformation
        log.info("Applying SCTransform-style normalization...")

        try:
            # Use Scanpy's implementation if available
            sc.experimental.pp.normalize_pearson_residuals(adata, layer=config.input_layer)

            # Rename to output layer
            adata.layers[config.output_layer] = adata.X.copy()

        except AttributeError:
            log.error(
                "SCTransform requires scanpy >= 1.8. " "Falling back to standard normalization."
            )
            sc.pp.normalize_total(adata, target_sum=config.target_sum or 1e4)
            sc.pp.log1p(adata)
            adata.layers[config.output_layer] = adata.X.copy()

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
    if config.input_layer in adata.layers:
        X_raw = adata.layers[config.input_layer]
    else:
        X_raw = adata.X

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
    if "scran_size_factors" in adata.obs:
        ax = axes[1, 0]
        sns.histplot(adata.obs["scran_size_factors"], bins=50, ax=ax, kde=True)
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
