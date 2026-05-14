"""
Pure Python implementation of bulk RNA-seq deconvolution tools.

Provides a unified interface to two R-free deconvolution backends:

- BayesPrism: Bayesian cell type proportion inference (Gibbs sampling).
- DWLS: Dampened Weighted Least Squares deconvolution (Tsoucas et al., 2019).

Both backends use only Python/NumPy/SciPy/scikit-learn and require no rpy2
or R installation.
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
    method: Literal["BayesPrism", "DWLS"] = "BayesPrism",
    key_added: str = "bulk_deconvolution",
    **method_kwargs,
) -> AnnData:
    """
    Estimate cell type proportions in bulk RNA-seq data using pure Python.

    Parameters
    ----------
    adata_ref : AnnData
        Single-cell reference data (cells x genes).
    bulk_data : pd.DataFrame
        Bulk RNA-seq data (genes x samples).
    cell_type_key : str
        Column in ``adata_ref.obs`` with cell type labels.
    sample_key : str, default="sampleID"
        Column in ``adata_ref.obs`` with sample IDs (used by BayesPrism).
    method : {"BayesPrism", "DWLS"}, default="BayesPrism"
        Deconvolution backend.

        - ``"BayesPrism"``: Bayesian inference with Gibbs sampling.
        - ``"DWLS"``: Dampened Weighted Least Squares per Tsoucas et al. 2019.
    key_added : str, default="bulk_deconvolution"
        Key under ``adata_ref.uns["sclucid"]["tools"]`` where results are stored.
    **method_kwargs
        Method-specific parameters. For DWLS: ``dampen_factor``, ``n_markers``,
        ``min_cells``, ``method`` (signature aggregation). For BayesPrism:
        ``n_iter``, ``n_chains``, ``burnin``.

    Returns:
    -------
    AnnData
        ``adata_ref`` with deconvolution results stored under
        ``.uns["sclucid"]["tools"][key_added]``.

    Examples:
    --------
    >>> adata_ref = sc.read_h5ad("sc_reference.h5ad")
    >>> bulk_data = pd.read_csv("bulk_rnaseq.csv", index_col=0)
    >>> adata_ref = deconvolve_bulk(
    ...     adata_ref, bulk_data,
    ...     cell_type_key="cell_type",
    ...     method="DWLS",
    ...     dampen_factor=1.0,
    ... )
    """
    common_genes = adata_ref.var_names.intersection(bulk_data.index)
    if len(common_genes) < 100:
        raise ValueError(
            f"Only {len(common_genes)} common genes found. "
            "Need at least 100 common genes for reliable deconvolution."
        )

    log.info(f"Found {len(common_genes)} common genes for deconvolution.")

    adata_ref_sub = adata_ref[:, common_genes].copy()
    bulk_data_sub = bulk_data.loc[common_genes]

    if method == "BayesPrism":
        proportions_df = _run_bayesprism(
            adata_ref_sub, bulk_data_sub, cell_type_key, sample_key, **method_kwargs
        )
    elif method == "DWLS":
        proportions_df = _run_dwls(
            adata_ref_sub, bulk_data_sub, cell_type_key, **method_kwargs
        )
    else:
        raise ValueError(
            f"Unknown method: {method!r}. Supported methods: 'BayesPrism', 'DWLS'."
        )

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

    log.info(
        f"Deconvolution complete. Results stored in .uns['sclucid']['tools']['{key_added}']"
    )
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

    reference = pd.DataFrame(
        adata_ref.X.T if hasattr(adata_ref.X, "toarray") else adata_ref.X.T,
        index=adata_ref.var_names,
        columns=adata_ref.obs_names,
    )

    cell_type_labels = adata_ref.obs[cell_type_key]

    prism_ref = BayesPrismReference(
        reference=reference, cell_type_labels=cell_type_labels, pseudo_min=1e-8
    )

    proportions = {}
    for sample_id in bulk_data.columns:
        bulk_expr = bulk_data[sample_id].values
        config = PrismConfig(n_iter=n_iter, n_chains=n_chains, burnin=burnin, **kwargs)
        theta = _bayesprism_gibbs_sample(prism_ref.phi, bulk_expr, config)
        proportions[sample_id] = theta

    proportions_df = pd.DataFrame(proportions, index=prism_ref.cell_types).T
    return proportions_df


def _bayesprism_gibbs_sample(
    phi: np.ndarray, bulk_expr: np.ndarray, config: PrismConfig
) -> np.ndarray:
    """Simplified Gibbs sampling for BayesPrism."""
    n_cell_types = phi.shape[1]

    theta_init, _ = nnls(phi, bulk_expr)
    theta_init = theta_init / (theta_init.sum() + 1e-10)

    theta_samples = []
    theta_curr = theta_init.copy()

    for i in range(config.n_iter + config.burnin):
        for k in range(n_cell_types):
            alpha = phi[:, k] @ bulk_expr + 1.0
            theta_curr[k] = np.random.gamma(alpha, 1.0)

        theta_curr = theta_curr / (theta_curr.sum() + 1e-10)

        if i >= config.burnin:
            theta_samples.append(theta_curr.copy())

    return np.mean(theta_samples, axis=0) if theta_samples else theta_init


def _run_dwls(
    adata_ref: AnnData,
    bulk_data: pd.DataFrame,
    cell_type_key: str,
    dampen_factor: float = 1.0,
    n_markers: Optional[int] = 50,
    min_cells: int = 10,
    signature_method: Literal["mean", "trimmed_mean"] = "mean",
    **kwargs,
) -> pd.DataFrame:
    """
    Run DWLS deconvolution using the full ``scLucid.tools.pyDWLS.DWLS`` class.

    Parameters
    ----------
    adata_ref : AnnData
        Single-cell reference, already restricted to genes shared with bulk.
    bulk_data : pd.DataFrame
        Bulk expression matrix (genes x samples).
    cell_type_key : str
        Column in ``adata_ref.obs`` with cell type labels.
    dampen_factor : float, default=1.0
        DWLS dampening factor; ``0`` reduces to ordinary NNLS.
    n_markers : int or None, default=50
        Markers per cell type to select. ``None`` skips marker selection and
        uses all common genes for the signature.
    min_cells : int, default=10
        Minimum cells per cell type when building the signature.
    signature_method : {"mean", "trimmed_mean"}, default="mean"
        Aggregation method for the signature matrix.

    Returns:
    -------
    pd.DataFrame
        Cell-type proportions (samples x cell types), each row summing to 1.
    """
    log.info("Running DWLS deconvolution via scLucid.tools.pyDWLS.DWLS")

    sc_expr = adata_ref.X
    if hasattr(sc_expr, "toarray"):
        sc_expr = sc_expr.toarray()
    sc_data = pd.DataFrame(
        sc_expr.T,
        index=adata_ref.var_names,
        columns=adata_ref.obs_names,
    )
    cell_type_labels = pd.Series(
        adata_ref.obs[cell_type_key].values, index=adata_ref.obs_names
    )

    dwls = DWLS(dampen_factor=dampen_factor, use_nonneg=True)

    genes_to_use: Optional[list]
    if n_markers is not None:
        genes_to_use = dwls.select_marker_genes(
            sc_data=sc_data,
            cell_type_labels=cell_type_labels,
            n_markers=n_markers,
        )
        log.info("DWLS selected %d marker genes", len(genes_to_use))
    else:
        genes_to_use = None

    dwls.build_signature_matrix(
        sc_data=sc_data,
        cell_type_labels=cell_type_labels,
        genes_to_use=genes_to_use,
        method=signature_method,
        min_cells=min_cells,
    )

    common_genes = dwls.signature_matrix.index.intersection(bulk_data.index)
    bulk_aligned = bulk_data.loc[common_genes]

    proportions_df = dwls.deconvolve(bulk_aligned, verbose=False)
    return proportions_df


def _build_signature_matrix(adata_ref: AnnData, cell_type_key: str) -> pd.DataFrame:
    """Build cell type signature matrix from single-cell reference."""
    cell_types = adata_ref.obs[cell_type_key].unique()

    signatures = {}
    for ct in cell_types:
        mask = adata_ref.obs[cell_type_key] == ct
        if hasattr(adata_ref.X, "toarray"):
            sig = np.array(adata_ref[mask].X.mean(axis=0)).flatten()
        else:
            sig = adata_ref[mask].X.mean(axis=0).flatten()
        signatures[ct] = sig

    signature_df = pd.DataFrame(signatures, index=adata_ref.var_names)
    return signature_df


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

    Parameters
    ----------
    proportions_df : pd.DataFrame
        Cell type proportions (samples x cell_types).
    metadata_df : pd.DataFrame
        Clinical metadata indexed by sample ID.
    group_col : str
        Column in ``metadata_df`` defining groups.
    group1, group2 : str
        Group values to compare.
    method : {"ttest", "wilcoxon"}, default="wilcoxon"
        Statistical test.

    Returns:
    -------
    pd.DataFrame
        Per-cell-type test statistics, p-values, and abundance summaries.
    """
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
    Correlate cell type abundance with a continuous clinical variable.

    Parameters
    ----------
    proportions_df : pd.DataFrame
        Cell type proportions (samples x cell_types).
    metadata_df : pd.DataFrame
        Clinical metadata indexed by sample ID.
    clinical_variable : str
        Continuous variable in ``metadata_df`` to correlate against.
    method : {"pearson", "spearman"}, default="spearman"
        Correlation method.

    Returns:
    -------
    pd.DataFrame
        Per-cell-type correlation coefficients and p-values.
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
