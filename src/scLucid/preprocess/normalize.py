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

from .config import NormalizationConfig

log = logging.getLogger(__name__)

__all__ = ["normalize_data", "plot_normalization_effect"]


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
    config: Optional[NormalizationConfig] = None,
    **kwargs,
) -> AnnData:
    """
    Normalize raw count data using a flexible, config-driven workflow.

    This function supports both a full config object for reproducibility and
    keyword arguments (`**kwargs`) for interactive overrides.

    Args:
        adata: The AnnData object to be normalized. Modified in place.
        config: A NormalizationConfig object. If None, a default config is created.
        **kwargs: Keyword arguments to override parameters in the config object
                  (e.g., `method='scran'`, `target_sum=1e5`, `input_layer='raw'`).

    Returns:
        The modified AnnData object with the new normalized layer.
    """
    # --- 1. Establish the final configuration ---
    if config is None:
        active_config = NormalizationConfig()
    else:
        # Create a copy of the config to avoid modifying the original
        # Filter kwargs to only include valid fields
        valid_fields = set(config.model_fields.keys())
        update_dict = {k: v for k, v in kwargs.items() if k in valid_fields}
        active_config = config.model_copy(update=update_dict)

    # Apply overrides from kwargs, making the function highly interactive
    for key, value in kwargs.items():
        if key in active_config.model_fields:
            setattr(active_config, key, value)
        elif key != "force":  # 'force' is handled separately
            log.warning(f"Ignoring unknown normalization parameter: '{key}'")

    # --- 2. Extract parameters from the final config ---
    input_layer = active_config.input_layer
    output_layer = active_config.output_layer
    force = kwargs.get("force", False)
    report = active_config.report
    plot = active_config.plot
    save_dir = Path(active_config.save_dir) if active_config.save_dir else None

    # --- 3. Input validation and data diagnosis ---
    if input_layer == "X":
        source_data = adata.X
        source_name = "adata.X"
    elif input_layer in adata.layers:
        source_data = adata.layers[input_layer]
        source_name = f"adata.layers['{input_layer}']"
    else:
        raise ValueError(f"Input layer '{input_layer}' not found.")

    if output_layer in adata.layers and not force:
        log.info(f"Layer '{output_layer}' already exists. Use force=True to overwrite.")
        return adata

    log.info(f"Diagnosing input data for normalization from '{source_name}' ...")
    stats_before = _diagnose_matrix(source_data, name="input")
    n_cells, n_genes = source_data.shape

    # Heuristics: warn if looks already transformed
    try:
        if scipy.sparse.issparse(source_data):
            arr = source_data.data
            min_val = source_data.min()
        else:
            arr = source_data
            min_val = np.min(arr)
        if np.issubdtype(arr.dtype, np.floating):
            if min_val < 0:
                log.warning(
                    f"'{source_name}' has negative values. Likely already transformed/residualized."
                )
            # If too few zeros, also warn
            if stats_before["zero_frac"] < 0.2:
                log.warning("Input matrix has few zeros; may not be raw UMI counts.")
    except Exception:
        pass

    log.info(
        f"Normalizing data from '{source_name}' using method '{active_config.method}'."
    )

    # --- 4. Core Normalization Logic ---
    norm_kwargs = {}
    for param in ["target_sum", "exclude_highly_expressed", "max_fraction"]:
        if hasattr(active_config, param):
            norm_kwargs[param] = getattr(active_config, param)

    temp_adata = AnnData(
        X=source_data.copy(), obs=adata.obs.copy(), var=adata.var.copy()
    )
    method_is_log_transformed = False

    try:
        if active_config.method == "standard":
            log.info(
                f"Applying standard library size normalization (params: {norm_kwargs})"
            )
            sc.pp.normalize_total(temp_adata, inplace=True, **norm_kwargs)
        elif active_config.method == "scran":
            try:
                import scanpy.external.pp as scepp
            except ImportError:
                log.error(
                    "scanpy.external.pp.scran_normalize not found. Install scanpy[external] and rpy2."
                )
                raise RuntimeError(
                    "scran normalization requires scanpy[external] and rpy2. See Scanpy docs."
                )
            log.warning("Method 'scran' requires a correctly configured R environment.")
            scepp.scran_normalize(temp_adata, inplace=True)
            method_is_log_transformed = True  # scran yields log-normalized-like output
        elif active_config.method == "pearson_residuals":
            log.info("Applying Pearson residuals normalization (experimental).")
            sc.experimental.pp.normalize_pearson_residuals(temp_adata, inplace=True)
            method_is_log_transformed = True
        elif active_config.method == "clr":
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
                f"Unknown normalization method: '{active_config.method}'. Choose from {valid_methods}."
            )
    except Exception as e:
        log.error(f"Normalization failed for method '{active_config.method}': {e}")
        raise RuntimeError(
            "Failed to normalize data. Check dependencies and data format."
        )

    # --- 5. Log transform and store results ---
    final_log_transformed = method_is_log_transformed
    if not method_is_log_transformed:
        log.info("Applying log1p transformation.")
        sc.pp.log1p(temp_adata)
        final_log_transformed = True

    adata.layers[output_layer] = temp_adata.X.copy()
    if active_config.update_X:
        adata.X = adata.layers[output_layer].copy()

    stats_after = _diagnose_matrix(adata.layers[output_layer], name="normalized")

    # --- 6. Store metadata in .uns ---
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})[
        "normalization"
    ] = {
        "params": active_config.to_dict(),  # Pydantic's built-in serialization
        "input_stats": stats_before,
        "output_stats": stats_after,
        "scanpy_version": getattr(sc, "__version__", "unknown"),
        "log_transformed": final_log_transformed,
        "input_layer": input_layer,
        "output_layer": output_layer,
    }

    # --- 7. Reporting and Plotting ---
    if plot:
        log.info("Generating diagnostic plots for normalization...")
        try:
            fig = _plot_normalization_global(
                source_data,
                adata.layers[output_layer],
                method=active_config.method,
                log_transformed=final_log_transformed,
                save_dir=save_dir,
            )
            plt.show()
        except Exception as e:
            log.warning(f"Failed to generate normalization plots: {e}")

    if report and save_dir:
        try:
            report_path = save_dir / "normalization_report.md"
            _write_normalization_report(
                report_path=report_path,
                stats_before=stats_before,
                stats_after=stats_after,
                config=active_config,
                n_cells=n_cells,
                n_genes=n_genes,
            )
            log.info(f"Normalization report written to {report_path}")
        except Exception as e:
            log.warning(f"Failed to write normalization report: {e}")

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
