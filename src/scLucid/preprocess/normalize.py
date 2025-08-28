"""
Normalization functions for single-cell RNA-seq data.

This module provides methods for normalizing raw count data and regressing out
unwanted sources of variation, preparing the data for downstream analysis.
"""

import logging
from pathlib import Path
from typing import List, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scipy.sparse
import seaborn as sns
from anndata import AnnData

from .config import NormalizationConfig, ScalingConfig

log = logging.getLogger(__name__)

__all__ = ["normalize_data", "regress_out", "plot_normalization_effect"]


# --- Helper Functions ---
def _diagnose_matrix(data, name="input", max_n=10000):
    """
    Compute and log basic statistics for a matrix or AnnData .X/.layers entry.
    """
    if scipy.sparse.issparse(data):
        arr = data
        total_entries = arr.shape[0] * arr.shape[1]
        nonzeros = arr.nnz
        zeros = total_entries - nonzeros
        mean = arr.mean()
        std = np.sqrt(arr.power(2).mean() - mean**2)
        min_val = arr.min()
        max_val = arr.max()
        zero_frac = zeros / total_entries
    else:
        arr = data
        if arr.shape[0] > max_n or arr.shape[1] > max_n:
            # sample to avoid OOM
            arr = arr[:max_n, :max_n]
        mean = np.mean(arr)
        std = np.std(arr)
        min_val = np.min(arr)
        max_val = np.max(arr)
        zero_frac = np.mean(arr == 0)
    log.info(
        f"[{name}] mean={mean:.2f}, std={std:.2f}, min={min_val:.2f}, max={max_val:.2f}, zero_frac={zero_frac:.2%}"
    )

    return dict(mean=mean, std=std, min=min_val, max=max_val, zero_frac=zero_frac)


def _plot_normalization_global(
    input_data: Union[np.ndarray, scipy.sparse.spmatrix],
    output_data: Union[np.ndarray, scipy.sparse.spmatrix],
    method: str = "standard",
    log_transformed: bool = True,
    save_dir: Optional[str] = None,
    max_cells: int = 40000,
) -> plt.Figure:
    """
    Generate global distribution plots comparing before and after normalization.
    For large data, sample cells to avoid OOM.
    """

    def _sample_rows(matrix, max_n):
        n = matrix.shape[0]
        if n <= max_n:
            return matrix
        idx = np.random.choice(n, max_n, replace=False)
        if scipy.sparse.issparse(matrix):
            return matrix[idx]
        else:
            return matrix[idx, :]

    # Sample before plotting if needed
    input_data = _sample_rows(input_data, max_cells)
    output_data = _sample_rows(output_data, max_cells)

    rc_params = {
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "text.color": "black",
        "axes.labelcolor": "black",
        "axes.edgecolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
    }

    with plt.rc_context(rc_params):
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 5))
        fig.suptitle(
            f"Data Distributions Before and After {method.capitalize()} Normalization",
            fontsize=16,
            color="black",
        )
        # Compute cell sums
        if scipy.sparse.issparse(input_data):
            before_sums = input_data.sum(axis=1).A1
        else:
            before_sums = input_data.sum(axis=1)
        if scipy.sparse.issparse(output_data):
            after_sums = output_data.sum(axis=1).A1
        else:
            after_sums = output_data.sum(axis=1)
        # Plot before normalization
        sns.histplot(
            before_sums,
            bins=100,
            kde=True,
            ax=axes[0],
            color="navy",
        )
        axes[0].set_title("Before Normalization")
        axes[0].set_xlabel("Total Counts per Cell")
        axes[0].set_ylabel("Frequency")
        axes[0].text(
            0.05,
            0.95,
            f"Mean: {before_sums.mean():.1f}\nMedian: {np.median(before_sums):.1f}",
            transform=axes[0].transAxes,
            va="top",
        )

        # Plot after normalization
        sns.histplot(
            after_sums,
            bins=100,
            kde=True,
            ax=axes[1],
            color="crimson",
        )
        title_suffix = " (Log-Transformed)" if log_transformed else ""
        axes[1].set_title(f"After {method.capitalize()} Normalization{title_suffix}")
        axes[1].set_xlabel("Sum of Normalized Values per Cell")
        axes[1].set_ylabel("Frequency")
        axes[1].text(
            0.05,
            0.95,
            f"Mean: {after_sums.mean():.1f}\nMedian: {np.median(after_sums):.1f}",
            transform=axes[1].transAxes,
            va="top",
        )

        plt.tight_layout(rect=[0, 0, 1, 0.96])

        if save_dir:
            save_path = Path(save_dir)
            save_path.mkdir(parents=True, exist_ok=True)
            figure_path = save_path / f"normalization_{method}_global.png"
            plt.savefig(figure_path, dpi=300)

    return fig


def _write_normalization_report(
    report_path: Path,
    stats_before: dict,
    stats_after: dict,
    config: NormalizationConfig,
    n_cells: int,
    n_genes: int,
) -> None:
    """
    Write a simple markdown report comparing normalization before and after.
    """
    with open(report_path, "w") as f:
        f.write("# Normalization Report\n\n")
        f.write(f"**Method:** {config.method}\n\n")
        f.write(f"**Input shape:** {n_cells} cells × {n_genes} genes\n\n")
        f.write("## Input Statistics\n")
        for k, v in stats_before.items():
            f.write(f"- {k}: {v:.3g}\n")
        f.write("\n## Output Statistics\n")
        for k, v in stats_after.items():
            f.write(f"- {k}: {v:.3g}\n")
        f.write("\n## Parameters\n")
        for k, v in config.__dict__.items():
            f.write(f"- {k}: {v}\n")


# --- Main Functions ---
def normalize_data(
    adata: AnnData,
    input_layer: Optional[str] = "counts",
    output_layer: Optional[str] = "normalized",
    config: Optional[NormalizationConfig] = None,
    force: bool = False,
    update_X: bool = True,
    report: bool = False,
    max_cells_plot: int = 40000,
) -> AnnData:
    """
    Normalize raw count data to account for differences in sequencing depth.

    Args:
        adata: The AnnData object to be normalized. Modified in place.
        config: A NormalizationConfig object specifying the method and parameters.
        input_layer: Name of the layer in `adata.layers` containing raw count data.
        output_layer: Name of the layer to store the normalized and log-transformed data.
        force: If True, overwrite the `output_layer` if it already exists.
        update_X: If True, also update adata.X to the newly normalized layer.
        report: Whether to write a markdown report with normalization statistics.
        max_cells_plot: Max cells to use when plotting (for large data).

    Returns:
        The modified AnnData object with the new normalized layer.

    Raises:
        ValueError: If the specified input layer does not exist or method is unknown.
        RuntimeError: If the normalization process fails for any reason.
    """
    if config is None:
        config = NormalizationConfig()

    # Input layer check: allow 'X' for adata.X
    if input_layer == "X":
        source_data = adata.X
        source_name = "adata.X"
    elif input_layer in adata.layers:
        source_data = adata.layers[input_layer]
        source_name = f"adata.layers['{input_layer}']"
    else:
        raise ValueError(
            f"Input layer '{input_layer}' not found in adata.layers nor is it 'X' (adata.X)."
        )

    if output_layer in adata.layers and not force:
        log.info(f"Layer '{output_layer}' already exists. Use force=True to overwrite.")
        return adata

    # --- Step 1: Input Data Diagnosis ---
    log.info(f"Diagnosing input data for normalization from '{source_name}' ...")
    stats_before = _diagnose_matrix(source_data, name="input")
    n_cells, n_genes = source_data.shape

    if np.issubdtype(source_data.dtype, np.floating):
        if scipy.sparse.issparse(source_data):
            arr = source_data.data
        else:
            arr = source_data
        non_int = np.abs(np.modf(arr)[0]).sum()
        if (np.min(arr) < 0) or (non_int > 1e-9):
            log.warning(
                f"Input data in '{source_name}' appears to be already normalized or transformed. "
                "Normalization should be run on raw integer counts."
            )
        if stats_before["zero_frac"] < 0.2:
            log.warning("Input matrix has few zeros; it may not be raw UMI counts.")

    log.info(f"Normalizing data from '{source_name}' using method '{config.method}'.")

    # Step 2: Prepare scanpy normalization arguments from config
    norm_kwargs = {}
    # These parameters exist in config, pass if present
    for param in ["target_sum", "exclude_highly_expressed", "max_fraction"]:
        if hasattr(config, param):
            norm_kwargs[param] = getattr(config, param)

    temp_adata = AnnData(
        X=source_data.copy(), obs=adata.obs.copy(), var=adata.var.copy()
    )
    method_is_log_transformed = False

    try:
        if config.method == "standard":
            log.info(
                f"Applying standard library size normalization (params: {norm_kwargs})"
            )
            sc.pp.normalize_total(temp_adata, inplace=True, **norm_kwargs)
        elif config.method == "scran":
            try:
                import scanpy.external.pp
            except ImportError:
                log.error(
                    "scanpy.external.pp.scran_normalize not found. Install scanpy[external] and rpy2."
                )
                raise RuntimeError(
                    "scran normalization requires scanpy[external] and rpy2. See https://scanpy.readthedocs.io/en/stable/api/scanpy.external.pp.scran_normalize.html"
                )
            log.warning("Method 'scran' requires a correctly configured R environment.")
            scanpy.external.pp.scran_normalize(temp_adata, inplace=True)
            method_is_log_transformed = True  # scran includes log-transform
        elif config.method == "pearson_residuals":
            log.info("Applying Pearson residuals normalization.")
            sc.experimental.pp.normalize_pearson_residuals(temp_adata, inplace=True)
            method_is_log_transformed = True
        elif config.method == "clr":
            log.info("Applying Centered Log-Ratio (CLR) normalization.")
            sc.pp.normalize_total(temp_adata, target_sum=1, inplace=True)
            sc.pp.log1p(temp_adata)
            mean_logs = temp_adata.X.mean(axis=1)
            if scipy.sparse.issparse(temp_adata.X):
                mean_logs = np.array(mean_logs).flatten()
                temp_adata.X = temp_adata.X - mean_logs[:, None]
            else:
                temp_adata.X = temp_adata.X - mean_logs[:, np.newaxis]
            method_is_log_transformed = True
        else:
            valid_methods = ["standard", "scran", "pearson_residuals", "clr"]
            raise ValueError(
                f"Unknown normalization method: '{config.method}'. Choose from {valid_methods}."
            )
    except Exception as e:
        log.error(f"Normalization failed for method '{config.method}': {e}")
        raise RuntimeError(
            "Failed to normalize data. Check dependencies and data format."
        )

    # Step 3: Final log1p transformation if needed
    final_log_transformed = method_is_log_transformed
    if not method_is_log_transformed:
        log.info("Applying log1p transformation.")
        sc.pp.log1p(temp_adata)
        final_log_transformed = True
    else:
        log.info(
            f"Method '{config.method}' includes a log-like transformation; skipping explicit log1p."
        )

    # Step 4: Store results
    adata.layers[output_layer] = temp_adata.X.copy()
    if update_X:
        adata.X = adata.layers[output_layer].copy()

    stats_after = _diagnose_matrix(adata.layers[output_layer], name="normalized")
    # Step 5: Store metadata in .uns for reproducibility
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})[
        "normalization"
    ] = {
        "method": config.method,
        "input_layer": input_layer,
        "output_layer": output_layer,
        "log_transformed": final_log_transformed,
        "params": config.__dict__,
        "input_stats": stats_before,
        "output_stats": stats_after,
        "n_cells": n_cells,
        "n_genes": n_genes,
    }

    # Step 6: Diagnostic plotting
    plot_flag = getattr(config, "plot_global_distribution", False)
    save_dir = getattr(config, "save_dir", None)
    if plot_flag:
        log.info("Generating diagnostic plots for normalization...")
        try:
            _plot_normalization_global(
                source_data,
                adata.layers[output_layer],
                method=config.method,
                log_transformed=final_log_transformed,
                save_dir=save_dir,
                max_cells=max_cells_plot,
            )
            plt.show()
        except Exception as e:
            log.warning(f"Failed to generate normalization plots: {e}")

    # Step 7: Optionally write a markdown report
    if report and save_dir:
        try:
            report_path = Path(save_dir) / "normalization_report.md"
            _write_normalization_report(
                report_path=report_path,
                stats_before=stats_before,
                stats_after=stats_after,
                config=config,
                n_cells=n_cells,
                n_genes=n_genes,
            )
            log.info(f"Normalization report written to {report_path}")
        except Exception as e:
            log.warning(f"Failed to write normalization report: {e}")

    return adata


def regress_out(
    adata: AnnData,
    config: ScalingConfig,
    input_layer: Optional[str] = "log1p_norm",
    output_layer: Optional[str] = "regressed_out",
    force: bool = False,
    update_X: bool = True,
) -> Optional[AnnData]:
    """
    Regress out unwanted sources of variation from gene expression data.

    Args:
        adata: The AnnData object to process.
        config: A ScalingConfig object containing `vars_to_regress`.
        input_layer: Layer containing the data to be corrected (e.g., 'log1p_norm').
        output_layer: Layer to store the regressed-out data.
        force: If True, overwrite the `output_layer` if it already exists.
        update_X: If True, also update adata.X to the regressed layer.

    Returns:
        The modified AnnData object.

    Raises:
        ValueError: If keys are not found in adata.obs or layers are invalid.
    """
    # Parameter Validation
    if not config.vars_to_regress:
        log.info(
            "No variables specified in `config.vars_to_regress`. Skipping regression."
        )
        if input_layer != output_layer:
            adata.layers[output_layer] = adata.layers[input_layer].copy()
        if update_X:
            adata.X = adata.layers[output_layer].copy()
        return adata

    if output_layer in adata.layers and not force:
        log.info(f"Layer '{output_layer}' already exists. Use force=True to overwrite.")
        return adata

    missing_keys = [key for key in config.vars_to_regress if key not in adata.obs]
    if missing_keys:
        raise ValueError(f"Keys not found in adata.obs: {', '.join(missing_keys)}")

    if input_layer not in adata.layers:
        raise ValueError(f"Input layer '{input_layer}' not found in adata.layers.")

    log.info(
        f"Regressing out: {', '.join(config.vars_to_regress)} from layer '{input_layer}'"
    )

    # Use a temporary object for the regression
    temp_adata = AnnData(X=adata.layers[input_layer].copy(), obs=adata.obs.copy())

    try:
        sc.pp.regress_out(temp_adata, keys=config.vars_to_regress)
    except Exception as e:
        log.error(f"Regression failed: {e}")
        raise RuntimeError("Failed to regress out variables.")

    adata.layers[output_layer] = temp_adata.X.copy()
    log.info(f"Regression complete. Results stored in 'adata.layers[{output_layer}]'.")

    if update_X:
        adata.X = adata.layers[output_layer].copy()

    # Store info to .uns for traceability
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["regress_out"] = {
        "vars_to_regress": config.vars_to_regress,
        "input_layer": input_layer,
        "output_layer": output_layer,
        "params": config.__dict__,
    }

    return adata


def plot_normalization_effect(
    adata: AnnData,
    original_layer: str,
    normalized_layer: str,
    log_transformed: bool = True,
    n_top_genes: int = 10,
    gene_subset: Optional[List[str]] = None,
    save_dir: Optional[str] = None,
    max_cells_plot: int = 40000,
) -> plt.Figure:
    """
    Generate comparison plots showing gene expression distributions before and after normalization.

    Args:
        adata: AnnData object with raw and normalized data
        original_layer: Layer containing original data
        normalized_layer: Layer containing normalized data
        log_transformed: Whether the normalized data is log-transformed
        n_top_genes: Number of top genes to show (by mean expression)
        gene_subset: Specific genes to show instead of top expressing genes
        save_dir: Directory to save the generated figure
        max_cells_plot: Max cells to use when plotting (for large data)

    Returns:
        matplotlib Figure object with the comparison plots

    Raises:
        ValueError: If specified layers don't exist
    """
    if original_layer not in adata.layers:
        raise ValueError(f"Original layer '{original_layer}' not found in adata.layers")
    if normalized_layer not in adata.layers:
        raise ValueError(
            f"Normalized layer '{normalized_layer}' not found in adata.layers"
        )

    import scipy.sparse

    def _sample_rows(matrix, max_n):
        n = matrix.shape[0]
        if n <= max_n:
            return matrix
        idx = np.random.choice(n, max_n, replace=False)
        if scipy.sparse.issparse(matrix):
            return matrix[idx]
        else:
            return matrix[idx, :]

    # Sample for plotting if needed
    orig_data = _sample_rows(adata.layers[original_layer], max_cells_plot)
    norm_data = _sample_rows(adata.layers[normalized_layer], max_cells_plot)

    # Select genes to plot
    if gene_subset is not None:
        genes_to_plot = [g for g in gene_subset if g in adata.var_names]
        if not genes_to_plot:
            raise ValueError("None of the specified genes were found in the data")
        if len(genes_to_plot) < len(gene_subset):
            log.warning(
                f"Only {len(genes_to_plot)}/{len(gene_subset)} specified genes were found"
            )
    else:
        # Select top expressed genes
        mean_expr = np.array(orig_data.mean(axis=0)).flatten()
        top_genes_idx = np.argsort(-mean_expr)[:n_top_genes]
        genes_to_plot = adata.var_names[top_genes_idx].tolist()

    n_to_plot = len(genes_to_plot)
    fig, axes = plt.subplots(n_to_plot, 2, figsize=(14, 3 * n_to_plot))

    if n_to_plot == 1:
        axes = np.array([axes])

    def get_data(matrix, gene):
        gene_idx = adata.var.index.get_loc(gene)
        data = matrix[:, gene_idx]
        if scipy.sparse.issparse(data):
            return data.toarray().flatten()
        return data.flatten()

    for i, gene in enumerate(genes_to_plot):
        try:
            before_data = get_data(orig_data, gene)
            after_data = get_data(norm_data, gene)

            sns.histplot(before_data, bins=50, kde=True, ax=axes[i, 0])
            axes[i, 0].set_title(f"{gene} - Before")
            axes[i, 0].set_ylabel("Frequency")
            axes[i, 0].text(
                0.05,
                0.95,
                f"Mean: {before_data.mean():.3f}\nStd: {before_data.std():.3f}\n"
                f"% zeros: {(before_data == 0).mean():.1%}",
                transform=axes[i, 0].transAxes,
                va="top",
            )

            sns.histplot(after_data, bins=50, kde=True, ax=axes[i, 1])
            suffix = " (Log-transformed)" if log_transformed else ""
            axes[i, 1].set_title(f"{gene} - After{suffix}")
            axes[i, 1].text(
                0.05,
                0.95,
                f"Mean: {after_data.mean():.3f}\nStd: {after_data.std():.3f}\n"
                f"% zeros: {(after_data == 0).mean():.1%}",
                transform=axes[i, 1].transAxes,
                va="top",
            )
        except Exception as e:
            log.warning(f"Failed to plot gene {gene}: {str(e)}")
            axes[i, 0].text(0.5, 0.5, f"Error plotting {gene}", ha="center")
            axes[i, 1].text(0.5, 0.5, f"Error plotting {gene}", ha="center")

    plt.tight_layout()

    if save_dir:
        save_path = Path(save_dir)
        save_path.mkdir(parents=True, exist_ok=True)
        figure_path = save_path / "normalization_effect_summary.png"
        plt.savefig(figure_path, dpi=300, bbox_inches="tight")
        log.info(f"Saved normalization effect plot to {figure_path}")

    return fig
