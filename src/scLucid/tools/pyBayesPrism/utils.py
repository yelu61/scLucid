"""
Utility functions for BayesPrism (R-free)

Helper functions for gene filtering and validation.
"""

import logging
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

log = logging.getLogger(__name__)


def cleanup_genes(
    gene_names: List[str],
    remove_ribo: bool = True,
    remove_mito: bool = True,
    remove_sex: bool = True,
    custom_remove: Optional[List[str]] = None,
) -> List[str]:
    """
    Clean up gene list by removing unwanted gene categories

    Parameters
    ----------
    gene_names : List[str]
        List of gene names
    remove_ribo : bool
        Remove ribosomal genes (RPS*, RPL*)
    remove_mito : bool
        Remove mitochondrial genes (MT-*)
    remove_sex : bool
        Remove sex chromosome genes (XIST, Y*)
    custom_remove : List[str], optional
        Additional gene patterns to remove

    Returns:
    -------
    List[str]
        Cleaned gene list

    Examples:
    --------
    >>> genes = ['RPS14', 'MT-CO1', 'TP53', 'XIST', 'BRCA1']
    >>> clean = cleanup_genes(genes, remove_ribo=True, remove_mito=True)
    >>> print(clean)
    ['TP53', 'BRCA1']
    """
    cleaned_genes = []

    for gene in gene_names:
        gene_upper = gene.upper()
        skip = False

        # Check ribosomal genes
        if remove_ribo:
            if gene_upper.startswith("RPS") or gene_upper.startswith("RPL"):
                skip = True

        # Check mitochondrial genes
        if remove_mito and not skip:
            if gene_upper.startswith("MT-") or gene_upper.startswith("MT"):
                skip = True

        # Check sex chromosome genes
        if remove_sex and not skip:
            if gene_upper == "XIST" or gene_upper.startswith("Y-"):
                skip = True
            if "CHROMOSOME_Y" in gene_upper:
                skip = True

        # Check custom patterns
        if custom_remove and not skip:
            for pattern in custom_remove:
                if pattern.upper() in gene_upper:
                    skip = True
                    break

        if not skip:
            cleaned_genes.append(gene)

    removed = len(gene_names) - len(cleaned_genes)
    log.debug(f"Removed {removed} genes, kept {len(cleaned_genes)}")

    return cleaned_genes


def find_outlier_genes(
    mixture: pd.DataFrame,
    reference_genes: List[str],
    cutoff: float = 0.01,
    min_samples: float = 0.1,
) -> List[str]:
    """
    Find outlier genes not in reference but highly expressed in mixture

    Parameters
    ----------
    mixture : pd.DataFrame
        Bulk expression matrix (genes x samples)
    reference_genes : List[str]
        Reference gene list
    cutoff : float
        Expression threshold relative to max
    min_samples : float
        Minimum fraction of samples exceeding threshold

    Returns:
    -------
    List[str]
        List of outlier genes
    """
    outliers = []
    mix_genes = set(mixture.index)
    ref_genes = set(reference_genes)

    novel_genes = mix_genes - ref_genes

    for gene in novel_genes:
        expr = mixture.loc[gene]
        max_expr = mixture.max().max()

        # Check if highly expressed in enough samples
        high_expr_ratio = (expr > max_expr * cutoff).sum() / len(mixture.columns)

        if high_expr_ratio >= min_samples:
            outliers.append(gene)

    log.info(f"Found {len(outliers)} outlier genes")
    return outliers


def compute_correlation(
    deconv_result: pd.DataFrame,
    true_fraction: pd.DataFrame,
    method: str = "pearson",
) -> pd.DataFrame:
    """
    Compute correlation between deconvolution result and ground truth

    Parameters
    ----------
    deconv_result : pd.DataFrame
        Predicted cell type fractions
    true_fraction : pd.DataFrame
        True cell type fractions
    method : str
        Correlation method ("pearson" or "spearman")

    Returns:
    -------
    pd.DataFrame
        Correlation statistics per cell type
    """
    results = []

    for cell_type in deconv_result.columns:
        if cell_type not in true_fraction.columns:
            continue

        pred = deconv_result[cell_type].values
        true = true_fraction[cell_type].values

        if method == "pearson":
            r, p = pearsonr(pred, true)
        elif method == "spearman":
            r, p = spearmanr(pred, true)
        else:
            raise ValueError(f"Unknown method: {method}")

        results.append(
            {
                "cell_type": cell_type,
                f"{method}_r": r,
                f"{method}_p": p,
            }
        )

    return pd.DataFrame(results)


def compute_rmse(
    predicted: pd.DataFrame,
    actual: pd.DataFrame,
) -> pd.DataFrame:
    """
    Compute root mean squared error

    Parameters
    ----------
    predicted : pd.DataFrame
        Predicted fractions
    actual : pd.DataFrame
        Actual fractions

    Returns:
    -------
    pd.DataFrame
        RMSE per cell type
    """
    results = []

    for cell_type in predicted.columns:
        if cell_type not in actual.columns:
            continue

        pred = predicted[cell_type].values
        true = actual[cell_type].values

        mse = np.mean((pred - true) ** 2)
        rmse = np.sqrt(mse)

        results.append(
            {
                "cell_type": cell_type,
                "rmse": rmse,
            }
        )

    return pd.DataFrame(results)


def normalize_expression(
    expression: pd.DataFrame,
    method: str = "cpm",
) -> pd.DataFrame:
    """
    Normalize expression data

    Parameters
    ----------
    expression : pd.DataFrame
        Raw expression counts
    method : str
        Normalization method ("cpm" or "tpm")

    Returns:
    -------
    pd.DataFrame
        Normalized expression
    """
    if method == "cpm":
        # Counts per million
        lib_sizes = expression.sum(axis=0)
        normalized = expression.div(lib_sizes, axis=1) * 1e6
    elif method == "tpm":
        # Transcripts per million (simplified)
        gene_lengths = pd.Series(1, index=expression.index)  # Placeholder
        rpk = expression.div(gene_lengths, axis=0)
        lib_sizes = rpk.sum(axis=0)
        normalized = rpk.div(lib_sizes, axis=1) * 1e6
    else:
        raise ValueError(f"Unknown method: {method}")

    return normalized


def batch_correct(
    expression: pd.DataFrame,
    batch_labels: pd.Series,
    method: str = "combat",
) -> pd.DataFrame:
    """
    Perform batch correction

    Parameters
    ----------
    expression : pd.DataFrame
        Expression matrix (genes x samples)
    batch_labels : pd.Series
        Batch labels for each sample
    method : str
        Batch correction method ("combat" or "mean_center")

    Returns:
    -------
    pd.DataFrame
        Batch-corrected expression
    """
    if method == "mean_center":
        corrected = expression.copy()

        for batch in batch_labels.unique():
            batch_mask = batch_labels == batch
            batch_mean = expression.loc[:, batch_mask].mean(axis=1)
            overall_mean = expression.mean(axis=1)

            corrected.loc[:, batch_mask] = (
                expression.loc[:, batch_mask].sub(batch_mean, axis=0).add(overall_mean, axis=0)
            )

        return corrected

    elif method == "combat":
        try:
            from combat.pycombat import pycombat

            return pycombat(expression, batch_labels)
        except ImportError:
            raise ImportError("combat package required for ComBat correction")
    else:
        raise ValueError(f"Unknown method: {method}")


def subsample_cells(
    expression: pd.DataFrame,
    cell_labels: pd.Series,
    n_cells_per_type: int = 100,
    random_state: int = 42,
) -> Tuple[pd.DataFrame, pd.Series]:
    """
    Subsample cells to balance cell type representation

    Parameters
    ----------
    expression : pd.DataFrame
        Expression matrix (genes x cells)
    cell_labels : pd.Series
        Cell type labels
    n_cells_per_type : int
        Number of cells per type to keep
    random_state : int
        Random seed

    Returns:
    -------
    Tuple[pd.DataFrame, pd.Series]
        Subsampled expression and labels
    """
    np.random.seed(random_state)

    selected_cells = []

    for cell_type in cell_labels.unique():
        type_cells = cell_labels[cell_labels == cell_type].index

        if len(type_cells) > n_cells_per_type:
            selected = np.random.choice(
                type_cells,
                n_cells_per_type,
                replace=False,
            )
        else:
            selected = type_cells

        selected_cells.extend(selected)

    return expression.loc[:, selected_cells], cell_labels.loc[selected_cells]


def validate_inputs(
    reference_expression,
    mixture_expression: pd.DataFrame,
    cell_type_labels: pd.Series,
) -> bool:
    """
    Validate input data for BayesPrism

    Parameters
    ----------
    reference_expression
        Reference scRNA-seq data
    mixture_expression : pd.DataFrame
        Bulk mixture data
    cell_type_labels : pd.Series
        Cell type labels

    Returns:
    -------
    bool
        True if valid

    Raises:
    ------
    ValueError
        If validation fails
    """
    # Check for empty inputs
    if mixture_expression.empty:
        raise ValueError("Mixture expression is empty")

    if len(cell_type_labels) == 0:
        raise ValueError("Cell type labels are empty")

    # Check for NaN/Inf
    if np.any(np.isnan(mixture_expression.values)):
        raise ValueError("Mixture contains NaN values")

    if np.any(np.isinf(mixture_expression.values)):
        raise ValueError("Mixture contains Inf values")

    # Check for negative values
    if np.any(mixture_expression.values < 0):
        raise ValueError("Mixture contains negative values")

    log.info("Input validation passed")
    return True
