"""
Pure Python implementation of bulk RNA-seq deconvolution tools.

This module provides R-free implementations of:
- BayesPrism: Bayesian cell type proportion inference
- DWLS: Dampened Weighted Least Squares deconvolution
- Bisque: Marker-based deconvolution (simplified Python implementation)

All implementations use only Python/NumPy/SciPy/scikit-learn, no rpy2 required.
"""

import logging
from typing import Literal, Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.optimize import nnls
from scipy.stats import mannwhitneyu, pearsonr, spearmanr, ttest_ind

from .pyBayesPrism import BayesPrismReference, PrismConfig
from .pyDWLS import DWLS

log = logging.getLogger(__name__)


def deconvolve_bulk(
    adata_ref: AnnData,
    bulk_data: pd.DataFrame,
    cell_type_key: str,
    sample_key: str = "sampleID",
    method: Literal["DWLS", "BayesPrism", "Bisque", "NNLS"] = "BayesPrism",
    key_added: str = "bulk_deconvolution",
    **method_kwargs,
) -> AnnData:
    """
    Estimate cell type proportions in bulk RNA-seq data using pure Python.

    This function provides a unified interface to multiple deconvolution methods,
    all implemented in Python without R dependencies.

    Args:
        adata_ref: Single-cell reference data (genes x cells)
        bulk_data: Bulk RNA-seq data (genes x samples)
        cell_type_key: Column in adata_ref.obs with cell type labels
        sample_key: Column in adata_ref.obs with sample IDs (for BayesPrism)
        method: Deconvolution method to use:
            - "BayesPrism": Bayesian inference with Gibbs sampling
            - "DWLS": Dampened Weighted Least Squares
            - "Bisque": Simplified marker-based approach
            - "NNLS": Non-negative Least Squares (baseline)
        key_added: Key for storing results in adata_ref.uns
        **method_kwargs: Method-specific parameters

    Returns:
        adata_ref with deconvolution results in .uns['sclucid']['tools'][key_added]

    Example:
        >>> adata_ref = sc.read_h5ad("sc_reference.h5ad")
        >>> bulk_data = pd.read_csv("bulk_rnaseq.csv", index_col=0)
        >>> adata_ref = deconvolve_bulk(
        ...     adata_ref, bulk_data,
        ...     cell_type_key="cell_type",
        ...     method="BayesPrism",
        ...     n_iter=100
        ... )
    """
    # Find common genes
    common_genes = adata_ref.var_names.intersection(bulk_data.index)
    if len(common_genes) < 100:
        raise ValueError(
            f"Only {len(common_genes)} common genes found. "
            "Need at least 100 common genes for reliable deconvolution."
        )

    log.info(f"Found {len(common_genes)} common genes for deconvolution.")

    # Subset to common genes
    adata_ref_sub = adata_ref[:, common_genes].copy()
    bulk_data_sub = bulk_data.loc[common_genes]

    # Method dispatch
    if method == "BayesPrism":
        proportions_df = _run_bayesprism(
            adata_ref_sub, bulk_data_sub, cell_type_key, sample_key, **method_kwargs
        )
    elif method == "DWLS":
        proportions_df = _run_dwls(adata_ref_sub, bulk_data_sub, cell_type_key, **method_kwargs)
    elif method == "Bisque":
        proportions_df = _run_bisque(adata_ref_sub, bulk_data_sub, cell_type_key, **method_kwargs)
    elif method == "NNLS":
        proportions_df = _run_nnls(adata_ref_sub, bulk_data_sub, cell_type_key)
    else:
        raise ValueError(f"Unknown method: {method}")

    # Store results
    adata_ref.uns.setdefault("sclucid", {}).setdefault("tools", {})
    adata_ref.uns["sclucid"]["tools"][key_added] = {
        "proportions": proportions_df,
        "params": {
            "method": method,
            "n_genes": len(common_genes),
            "n_cell_types": len(adata_ref.obs[cell_type_key].unique()),
            "n_samples": bulk_data.shape[1],
        },
    }

    log.info(f"Deconvolution complete. Results stored in .uns['sclucid']['tools']['{key_added}']")
    return adata_ref


def _run_bayesprism(
    adata_ref: AnnData,
    bulk_data: pd.DataFrame,
    cell_type_key: str,
    sample_key: str,
    n_iter: int = 100,
    n_chains: int = 4,
    burnin: int = 50,
    **kwargs,
) -> pd.DataFrame:
    """Run BayesPrism deconvolution."""
    log.info(f"Running BayesPrism with {n_iter} iterations, {n_chains} chains...")

    # Create reference
    reference = pd.DataFrame(
        adata_ref.X.T if hasattr(adata_ref.X, "toarray") else adata_ref.X.T,
        index=adata_ref.var_names,
        columns=adata_ref.obs_names,
    )

    cell_type_labels = adata_ref.obs[cell_type_key]

    # Initialize BayesPrism
    prism_ref = BayesPrismReference(
        reference=reference, cell_type_labels=cell_type_labels, pseudo_min=1e-8
    )

    # Run deconvolution for each bulk sample
    proportions = {}

    for sample_id in bulk_data.columns:
        bulk_expr = bulk_data[sample_id].values

        # Initialize BayesPrism model
        config = PrismConfig(n_iter=n_iter, n_chains=n_chains, burnin=burnin, **kwargs)

        # Run Gibbs sampling (simplified version)
        # Note: Full implementation would use the complete BayesPrism class
        theta = _bayesprism_gibbs_sample(prism_ref.phi, bulk_expr, config)

        proportions[sample_id] = theta

    proportions_df = pd.DataFrame(proportions, index=prism_ref.cell_types).T
    return proportions_df


def _bayesprism_gibbs_sample(
    phi: np.ndarray, bulk_expr: np.ndarray, config: PrismConfig
) -> np.ndarray:
    """Simplified Gibbs sampling for BayesPrism."""
    n_cell_types = phi.shape[1]

    # Initialize with NNLS
    theta_init, _ = nnls(phi, bulk_expr)
    theta_init = theta_init / (theta_init.sum() + 1e-10)

    # Simple MCMC (simplified for speed)
    theta_samples = []
    theta_curr = theta_init.copy()

    for i in range(config.n_iter + config.burnin):
        # Gibbs update for each cell type proportion
        for k in range(n_cell_types):
            # Sample from Dirichlet-like posterior
            alpha = phi[:, k] @ bulk_expr + 1.0  # Prior + Likelihood
            theta_curr[k] = np.random.gamma(alpha, 1.0)

        # Normalize
        theta_curr = theta_curr / (theta_curr.sum() + 1e-10)

        # Collect samples after burnin
        if i >= config.burnin:
            theta_samples.append(theta_curr.copy())

    # Return mean of samples
    return np.mean(theta_samples, axis=0) if theta_samples else theta_init


def _run_dwls(
    adata_ref: AnnData,
    bulk_data: pd.DataFrame,
    cell_type_key: str,
    dampening_factor: float = 0.1,
    **kwargs,
) -> pd.DataFrame:
    """Run DWLS deconvolution."""
    log.info("Running DWLS deconvolution...")

    # Build signature matrix
    signature_matrix = _build_signature_matrix(adata_ref, cell_type_key)

    # Initialize DWLS
    dwls = DWLS(signature_matrix=signature_matrix, bulk_data=bulk_data)

    # Run deconvolution
    proportions = {}

    for sample_id in bulk_data.columns:
        # Solve with dampening
        theta = _dwls_solve(signature_matrix.values, bulk_data[sample_id].values, dampening_factor)
        proportions[sample_id] = theta

    proportions_df = pd.DataFrame(proportions, index=signature_matrix.columns).T
    return proportions_df


def _dwls_solve(signature: np.ndarray, bulk: np.ndarray, dampening: float) -> np.ndarray:
    """Solve DWLS optimization problem."""
    n_cell_types = signature.shape[1]

    # Initial NNLS solution
    theta, _ = nnls(signature, bulk)

    # Dampening: downweight high-expression genes
    for _ in range(10):  # Iterative refinement
        # Calculate residuals
        predicted = signature @ theta
        residuals = np.abs(bulk - predicted)

        # Dampening weights
        weights = 1.0 / (1.0 + dampening * residuals)

        # Weighted NNLS (approximated by scaling)
        signature_weighted = signature * weights[:, np.newaxis]
        bulk_weighted = bulk * weights

        theta, _ = nnls(signature_weighted, bulk_weighted)

        if theta.sum() == 0:
            break

    # Normalize to sum to 1
    theta = theta / (theta.sum() + 1e-10)
    return theta


def _run_bisque(
    adata_ref: AnnData,
    bulk_data: pd.DataFrame,
    cell_type_key: str,
    marker_genes: Optional[list] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Simplified Bisque-like deconvolution using marker genes.

    This is a simplified implementation that uses marker gene expression
    to estimate cell type proportions.
    """
    log.info("Running Bisque-like deconvolution...")

    cell_types = adata_ref.obs[cell_type_key].unique()

    # Identify marker genes if not provided
    if marker_genes is None:
        marker_genes = _identify_markers(adata_ref, cell_type_key)

    # Subset to marker genes
    available_markers = [g for g in marker_genes if g in adata_ref.var_names]
    if len(available_markers) < 10:
        raise ValueError(f"Only {len(available_markers)} marker genes found")

    adata_marker = adata_ref[:, available_markers]

    # Build signature from markers
    signature = _build_signature_matrix(adata_marker, cell_type_key)

    # Simple linear regression on markers
    proportions = {}
    for sample_id in bulk_data.columns:
        bulk_marker = bulk_data.loc[available_markers, sample_id].values
        theta, _ = nnls(signature.values, bulk_marker)
        theta = theta / (theta.sum() + 1e-10)
        proportions[sample_id] = theta

    proportions_df = pd.DataFrame(proportions, index=cell_types).T
    return proportions_df


def _run_nnls(adata_ref: AnnData, bulk_data: pd.DataFrame, cell_type_key: str) -> pd.DataFrame:
    """Baseline NNLS deconvolution."""
    log.info("Running NNLS baseline deconvolution...")

    signature = _build_signature_matrix(adata_ref, cell_type_key)

    proportions = {}
    for sample_id in bulk_data.columns:
        theta, _ = nnls(signature.values, bulk_data[sample_id].values)
        theta = theta / (theta.sum() + 1e-10)
        proportions[sample_id] = theta

    proportions_df = pd.DataFrame(proportions, index=signature.columns).T
    return proportions_df


def _build_signature_matrix(adata_ref: AnnData, cell_type_key: str) -> pd.DataFrame:
    """Build cell type signature matrix from single-cell reference."""
    cell_types = adata_ref.obs[cell_type_key].unique()

    signatures = {}
    for ct in cell_types:
        mask = adata_ref.obs[cell_type_key] == ct
        # Mean expression per cell type
        if hasattr(adata_ref.X, "toarray"):
            sig = np.array(adata_ref[mask].X.mean(axis=0)).flatten()
        else:
            sig = adata_ref[mask].X.mean(axis=0).flatten()
        signatures[ct] = sig

    signature_df = pd.DataFrame(signatures, index=adata_ref.var_names)
    return signature_df


def _identify_markers(adata_ref: AnnData, cell_type_key: str, n_markers: int = 50) -> list:
    """Identify marker genes for each cell type."""
    from scipy.stats import ttest_ind

    cell_types = adata_ref.obs[cell_type_key].unique()
    markers = []

    for ct in cell_types:
        mask = adata_ref.obs[cell_type_key] == ct

        # t-test for each gene
        pvals = []
        for i in range(adata_ref.n_vars):
            ct_expr = adata_ref[mask, i].X.flatten()
            other_expr = adata_ref[~mask, i].X.flatten()
            _, pval = ttest_ind(ct_expr, other_expr)
            pvals.append(pval)

        # Top markers for this cell type
        top_idx = np.argsort(pvals)[:n_markers]
        markers.extend(adata_ref.var_names[top_idx].tolist())

    return list(set(markers))  # Remove duplicates


def differential_abundance(
    proportions_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    group_col: str,
    group1: str,
    group2: str,
    method: Literal["ttest", "wilcoxon"] = "wilcoxon",
) -> pd.DataFrame:
    """
    Perform differential abundance analysis on deconvolution results.

    Args:
        proportions_df: DataFrame of cell type proportions (samples x cell_types)
        metadata_df: DataFrame with clinical metadata (indexed by sample ID)
        group_col: Column defining groups
        group1: First group
        group2: Second group
        method: Statistical test

    Returns:
        DataFrame with differential abundance results
    """
    # Align data
    data = proportions_df.join(metadata_df, how="inner")

    group1_samples = data[data[group_col] == group1].index
    group2_samples = data[data[group_col] == group2].index

    results = []
    for cell_type in proportions_df.columns:
        scores1 = data.loc[group1_samples, cell_type].dropna()
        scores2 = data.loc[group2_samples, cell_type].dropna()

        if len(scores1) < 2 or len(scores2) < 2:
            continue

        if method == "wilcoxon":
            stat, pval = mannwhitneyu(scores1, scores2, alternative="two-sided")
        else:
            stat, pval = ttest_ind(scores1, scores2)

        results.append(
            {
                "cell_type": cell_type,
                "statistic": stat,
                "pvalue": pval,
                "mean_abundance_group1": scores1.mean(),
                "mean_abundance_group2": scores2.mean(),
                "log2fc_abundance": np.log2(scores1.mean() / (scores2.mean() + 1e-10)),
            }
        )

    results_df = pd.DataFrame(results).sort_values("pvalue")
    return results_df


def correlate_abundance_with_clinical(
    proportions_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    clinical_variable: str,
    method: Literal["pearson", "spearman"] = "spearman",
) -> pd.DataFrame:
    """
    Correlate cell type abundance with continuous clinical variable.

    Args:
        proportions_df: DataFrame of cell type proportions
        metadata_df: DataFrame with clinical metadata
        clinical_variable: Continuous variable to correlate
        method: Correlation method

    Returns:
        DataFrame with correlation results
    """
    data = proportions_df.join(metadata_df, how="inner")

    results = []
    for cell_type in proportions_df.columns:
        subset = data[[cell_type, clinical_variable]].dropna()
        if len(subset) < 5:
            continue

        if method == "pearson":
            corr, pval = pearsonr(subset[cell_type], subset[clinical_variable])
        else:
            corr, pval = spearmanr(subset[cell_type], subset[clinical_variable])

        results.append(
            {
                "cell_type": cell_type,
                "clinical_variable": clinical_variable,
                "correlation_coefficient": corr,
                "pvalue": pval,
            }
        )

    return pd.DataFrame(results).sort_values("pvalue")


# Backward compatibility alias
run_deconvolution = deconvolve_bulk

__all__ = [
    "deconvolve_bulk",
    "run_deconvolution",
    "differential_abundance",
    "correlate_abundance_with_clinical",
]
