"""
Data scaling and regression functions for single-cell RNA-seq data.

This module provides flexible, config-driven functions for scaling and
regressing out covariates, ensuring consistency with the scLucid workflow.
"""

import logging
from importlib.metadata import version
from pathlib import Path
from typing import Dict, List, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import seaborn as sns
from anndata import AnnData

from .config import ScalingConfig, apply_config_overrides
from .utils import validate_matrix_input

log = logging.getLogger(__name__)
MAD_NORMAL_CONSISTENCY_FACTOR = 1.4826

__all__ = [
    "diagnose_cell_cycle_regression",
    "scale_data",
    "regress_out",
    "plot_scaling_effect",
]


def _safe_corr(a: pd.Series, b: pd.Series) -> float:
    values = pd.concat(
        [pd.to_numeric(a, errors="coerce"), pd.to_numeric(b, errors="coerce")],
        axis=1,
    ).dropna()
    if len(values) < 3:
        return 0.0
    if values.iloc[:, 0].std() == 0 or values.iloc[:, 1].std() == 0:
        return 0.0
    return float(values.iloc[:, 0].corr(values.iloc[:, 1]))


def _eta_squared_numeric_by_group(values: pd.Series, groups: pd.Series) -> float:
    df = pd.DataFrame({"value": pd.to_numeric(values, errors="coerce"), "group": groups}).dropna()
    if df.empty or df["group"].nunique() < 2:
        return 0.0
    grand = df["value"].mean()
    ss_total = float(((df["value"] - grand) ** 2).sum())
    if ss_total <= 0:
        return 0.0
    ss_between = 0.0
    for _, sub in df.groupby("group", observed=False):
        ss_between += len(sub) * float((sub["value"].mean() - grand) ** 2)
    return float(max(0.0, min(1.0, ss_between / ss_total)))


def diagnose_cell_cycle_regression(
    adata: AnnData,
    *,
    condition_key: Optional[str] = None,
    batch_key: Optional[str] = None,
    cell_type_key: Optional[str] = None,
    tumor: bool = False,
    proliferation_markers: Optional[List[str]] = None,
    confounding_threshold: float = 0.15,
    technical_threshold: float = 0.20,
    key_added: str = "cell_cycle_regression_diagnostic",
    record: bool = False,
) -> Dict[str, object]:
    """Diagnose whether cell-cycle regression is biologically appropriate.

    The function deliberately does not run regression. It records whether cell
    cycle scores are associated with condition, batch, or cell type so users can
    decide whether regression would remove biology or mitigate technical bias.
    """
    has_scores = {"S_score", "G2M_score"}.issubset(adata.obs.columns)
    result: Dict[str, object] = {
        "schema_version": "cell_cycle_regression_diagnostic_v1",
        "scores_present": bool(has_scores),
        "status": "not_available",
        "metrics": {},
        "warnings": [],
        "recommendation": "Cell-cycle scores were not detected.",
    }
    if not has_scores:
        if record:
            adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})[key_added] = result
        return result

    cc_score = pd.to_numeric(adata.obs["S_score"], errors="coerce") - pd.to_numeric(
        adata.obs["G2M_score"], errors="coerce"
    )
    cycling_fraction = (
        adata.obs["phase"].astype(str).isin(["S", "G2M"]).mean()
        if "phase" in adata.obs.columns
        else float(((adata.obs["S_score"] > 0) | (adata.obs["G2M_score"] > 0)).mean())
    )
    metrics: Dict[str, object] = {"cycling_fraction": float(cycling_fraction)}
    warnings: List[str] = []

    for key_name, key in [
        ("condition", condition_key),
        ("batch", batch_key),
        ("cell_type", cell_type_key),
    ]:
        if key and key in adata.obs.columns:
            groups = adata.obs[key].astype(str)
            eta = _eta_squared_numeric_by_group(cc_score, groups)
            metrics[f"{key_name}_eta2_cc_score"] = eta
            group_sizes = groups.value_counts(dropna=True).astype(int).to_dict()
            if group_sizes:
                min_group_size = min(group_sizes.values())
                max_group_size = max(group_sizes.values())
                imbalance_ratio = (
                    float(max_group_size / min_group_size) if min_group_size > 0 else np.inf
                )
                metrics[f"{key_name}_group_sizes"] = {
                    str(group): int(size) for group, size in group_sizes.items()
                }
                metrics[f"{key_name}_group_size_imbalance_ratio"] = imbalance_ratio
                if imbalance_ratio >= 5.0:
                    warnings.append(
                        f"{key_name} groups are imbalanced; eta-squared should be treated as diagnostic evidence only."
                    )
            if key_name in {"condition", "cell_type"} and eta >= confounding_threshold:
                warnings.append(
                    f"Cell-cycle scores are associated with {key_name}; regression may remove biology."
                )
            if key_name == "batch" and eta >= technical_threshold:
                warnings.append(
                    "Cell-cycle scores are associated with batch; regression may be useful if biology is not confounded."
                )

    if "n_genes_by_counts" in adata.obs.columns:
        metrics["corr_cc_n_genes"] = abs(_safe_corr(cc_score, adata.obs["n_genes_by_counts"]))
    if "total_counts" in adata.obs.columns:
        metrics["corr_cc_total_counts"] = abs(_safe_corr(cc_score, adata.obs["total_counts"]))

    present_markers = []
    markers = proliferation_markers or ["MKI67", "TOP2A", "PCNA", "MCM6", "TYMS"]
    var_upper = {str(g).upper(): str(g) for g in adata.var_names}
    for marker in markers:
        if marker.upper() in var_upper:
            present_markers.append(var_upper[marker.upper()])
    if tumor or present_markers:
        warnings.append(
            "Proliferation may be a tumor or cell-state signal; avoid direct cell-cycle regression unless explicitly justified."
        )

    if any("batch" in warning for warning in warnings) and not any(
        "condition" in warning or "cell_type" in warning or "tumor" in warning.lower()
        for warning in warnings
    ):
        status = "technical_regression_candidate"
        recommendation = (
            "Cell cycle appears batch-associated without strong biology-confounding evidence; "
            "regression can be considered with pre/post preservation checks."
        )
    elif warnings:
        status = "review_required"
        recommendation = (
            "Do not regress cell cycle by default. Review condition, cell type, and tumor-state "
            "confounding before enabling regression."
        )
    else:
        status = "low_risk"
        recommendation = (
            "No strong cell-cycle confounding detected; keep regression disabled unless a specific "
            "technical rationale exists."
        )

    result.update(
        {
            "status": status,
            "metrics": metrics,
            "warnings": warnings,
            "present_proliferation_markers": present_markers,
            "recommendation": recommendation,
        }
    )
    if record:
        adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})[key_added] = result
    return result


# --- Helper functions for different scaling methods ---
def _robust_scale(X: np.ndarray, max_value: Optional[float]) -> np.ndarray:
    """Robustly scales a dense matrix X."""
    gene_medians = np.median(X, axis=0)
    # MAD calculation with normal-consistency factor so robust z-scores share the z-score scale.
    gene_mads = MAD_NORMAL_CONSISTENCY_FACTOR * np.median(np.abs(X - gene_medians), axis=0)
    gene_mads[gene_mads == 0] = 1e-8  # Avoid division by zero

    X_scaled = (X - gene_medians) / gene_mads

    if max_value is not None:
        X_scaled = np.clip(X_scaled, -max_value, max_value)
    return X_scaled


def _minmax_scale(X: np.ndarray) -> np.ndarray:
    """Scales a dense matrix X to a [0, 1] range."""
    gene_mins = np.min(X, axis=0)
    gene_ranges = np.max(X, axis=0) - gene_mins
    gene_ranges[gene_ranges == 0] = 1e-8  # Avoid division by zero

    return (X - gene_mins) / gene_ranges


SPARSE_ROBUST_SCALE_DENSE_THRESHOLD_ELEMENTS = 50_000_000


def _robust_scale_sparse(
    X: scipy.sparse.spmatrix, max_value: Optional[float]
) -> Union[np.ndarray, scipy.sparse.spmatrix]:
    """
    Robust scaling for sparse matrices that includes zeros in median/MAD.

    Zeros are biologically informative in scRNA-seq data, so they must participate
    in the median and MAD calculations. Because subtracting a non-zero median turns
    implicit zeros into non-zero values, the output is dense whenever any feature
    has a non-zero median.

    For matrices below ``SPARSE_ROBUST_SCALE_DENSE_THRESHOLD_ELEMENTS`` elements,
    the input is densified and the dense implementation is used (with a warning).
    For larger matrices, medians and MADs are computed from the sparse structure
    without materializing the full input, then the output is produced in column
    chunks to limit peak memory.
    """
    import warnings

    import scipy.sparse as sp

    if not sp.issparse(X):
        return _robust_scale(X, max_value)  # Fallback to dense version

    X_csc = X.tocsc()
    n_obs, n_genes = X_csc.shape
    n_elements = n_obs * n_genes

    # Compute sparse-aware medians and MADs column by column. This is done first
    # so we can preserve sparsity when every feature has median 0.
    medians = np.empty(n_genes, dtype=np.float64)
    mads = np.empty(n_genes, dtype=np.float64)

    for j in range(n_genes):
        start, end = X_csc.indptr[j], X_csc.indptr[j + 1]
        col_data = X_csc.data[start:end]
        n_nnz = len(col_data)
        n_zeros = n_obs - n_nnz

        if n_nnz == 0:
            medians[j] = 0.0
            mads[j] = 1e-8
            continue

        sorted_data = np.sort(col_data)
        neg_count = int(np.searchsorted(sorted_data, 0.0, side="left"))

        def _value_at(k: int) -> float:
            if k < neg_count:
                return float(sorted_data[k])
            if k < neg_count + n_zeros:
                return 0.0
            return float(sorted_data[k - n_zeros])

        if n_obs % 2 == 1:
            median = _value_at(n_obs // 2)
        else:
            median = (_value_at(n_obs // 2 - 1) + _value_at(n_obs // 2)) / 2.0

        medians[j] = median
        abs_median = abs(median)

        # MAD: median(|x - median|) over the full column, including zeros.
        dev_data = np.abs(sorted_data - median)
        dev_sorted = np.sort(dev_data)
        n_dev_less = int(np.searchsorted(dev_sorted, abs_median, side="left"))
        n_dev_equal = int(np.searchsorted(dev_sorted, abs_median, side="right")) - n_dev_less

        def _dev_value_at(k: int) -> float:
            if k < n_dev_less:
                return float(dev_sorted[k])
            if k < n_dev_less + n_dev_equal + n_zeros:
                return abs_median
            return float(dev_sorted[k - n_zeros])

        if n_obs % 2 == 1:
            mad = _dev_value_at(n_obs // 2)
        else:
            mad = (_dev_value_at(n_obs // 2 - 1) + _dev_value_at(n_obs // 2)) / 2.0

        mads[j] = max(mad * MAD_NORMAL_CONSISTENCY_FACTOR, 1e-8)

    mads[mads == 0] = 1e-8

    # If every median is zero, implicit zeros stay zero and we can keep a sparse
    # output. This is rare for positive count data but handles the edge case.
    if np.allclose(medians, 0.0):
        X_scaled = X_csc.copy()
        if X_scaled.nnz:
            counts_per_gene = np.diff(X_scaled.indptr)
            gene_indices = np.repeat(np.arange(n_genes), counts_per_gene)
            X_scaled.data = (X_scaled.data - medians[gene_indices]) / mads[gene_indices]
            if max_value is not None:
                X_scaled.data = np.clip(X_scaled.data, -max_value, max_value)
        return X_scaled.tocsr()

    # Small matrices: densify and use the dense implementation. This is the
    # simplest correct path and matches the reference implementation exactly.
    if n_elements <= SPARSE_ROBUST_SCALE_DENSE_THRESHOLD_ELEMENTS:
        warnings.warn(
            f"Sparse robust scaling is densifying the {n_obs}x{n_genes} matrix because "
            "zeros must be included in median/MAD. For very large datasets consider "
            "z-score scaling or Harmony/scVI integration.",
            UserWarning,
            stacklevel=2,
        )
        return _robust_scale(X_csc.toarray(), max_value)

    warnings.warn(
        "Robust scaling produced a dense output because at least one feature has a "
        "non-zero median; scaled zero entries are no longer zero. For very large "
        "datasets consider z-score scaling or Harmony/scVI integration.",
        UserWarning,
        stacklevel=2,
    )

    # Large matrices: materialize the dense scaled matrix in column chunks.
    X_scaled = np.empty((n_obs, n_genes), dtype=np.float64)
    chunk_size = max(1, SPARSE_ROBUST_SCALE_DENSE_THRESHOLD_ELEMENTS // n_obs)
    for j_start in range(0, n_genes, chunk_size):
        j_end = min(j_start + chunk_size, n_genes)
        chunk = X_csc[:, j_start:j_end].toarray()
        X_scaled[:, j_start:j_end] = (chunk - medians[j_start:j_end]) / mads[j_start:j_end]

    if max_value is not None:
        X_scaled = np.clip(X_scaled, -max_value, max_value)
    return X_scaled


def _minmax_scale_sparse(X: scipy.sparse.spmatrix) -> scipy.sparse.spmatrix:
    """
    MinMax scaling for sparse matrices.
    """
    import scipy.sparse as sp

    if not sp.issparse(X):
        return _minmax_scale(X)

    X_csc = X.tocsc()
    n_genes = X_csc.shape[1]

    mins = np.zeros(n_genes)
    maxs = np.zeros(n_genes)

    for i in range(n_genes):
        col_data = X_csc.getcol(i).data
        if len(col_data) > 0:
            mins[i] = col_data.min()
            maxs[i] = col_data.max()

    ranges = maxs - mins
    ranges[ranges == 0] = 1e-8

    X_scaled = X_csc.copy()
    for i in range(n_genes):
        col = X_scaled.getcol(i)
        col.data = (col.data - mins[i]) / ranges[i]
        X_scaled[:, i] = col

    return X_scaled.tocsr()


# --- Main Functions ---
def scale_data(
    adata: AnnData,
    config: Optional[ScalingConfig] = None,
    output_layer: Optional[str] = "scaled",
    **kwargs,
) -> AnnData:
    """
    Scales gene expression data using the specified method.

    Operates on adata.X and modifies it in place. Assumes adata has been
    subsetted to the desired features (e.g., HVGs).

    Args:
        adata: AnnData object (will be modified in place).
        config: A ScalingConfig object. If None, a default config is used.
        **kwargs: Keyword arguments to override parameters in the config object
                  (e.g., `max_value=15`, `scale_method='robust'`).

    Returns:
        The modified AnnData object with scaled adata.X.
        If `output_layer` is not None, scaled values are also stored in
        `adata.layers[output_layer]` for downstream compatibility.
    """
    # --- 1. Establish the final configuration ---
    if config is None:
        active_config = ScalingConfig()
    else:
        active_config = apply_config_overrides(config, **kwargs)

    log.info(
        f"Scaling data in .X (shape: {adata.shape}) using '{active_config.scale_method}' method."
    )

    # --- 2. Validate input matrix ---
    validate_matrix_input(adata.X, name="adata.X", allow_negative=True)

    # --- 3. Apply the scaling method ---
    if active_config.regress_in_scale:
        vars_reg = active_config.vars_to_regress_in_scale or active_config.vars_to_regress or []
        vars_reg = [v for v in vars_reg if v]

        if vars_reg:
            missing = [k for k in vars_reg if k not in adata.obs.columns]
            if missing:
                log.warning(f"Vars to regress not found: {missing}")
                vars_reg = [k for k in vars_reg if k in adata.obs.columns]

            if vars_reg:
                # === IMPROVED: Use config setting ===
                input_layer = active_config.input_layer_for_regress

                if input_layer not in adata.layers:
                    raise ValueError(
                        f"Layer '{input_layer}' specified in config.input_layer_for_regress "
                        f"not found in adata.layers. Available: {list(adata.layers.keys())}"
                    )

                X_in = adata.layers[input_layer].copy()
                temp = AnnData(X=X_in, obs=adata.obs.copy(), var=adata.var.copy())

                log.info(f"Regressing {vars_reg} from layer '{input_layer}' before scaling")
                sc.pp.regress_out(temp, keys=vars_reg)

                adata.X = temp.X.copy()

                # Store metadata
                adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})[
                    "regress_inline"
                ] = {
                    "vars_to_regress": vars_reg,
                    "input_layer": input_layer,
                    "timestamp": pd.Timestamp.now().isoformat(),
                }
            else:
                log.info("No valid variables to regress in scale step; skipping inline regression.")
        else:
            log.info("regress_in_scale=True but no variables provided; skipping inline regression.")

    # --- 3. Scale the data ---
    try:
        if active_config.scale_method == "zscore":
            sc.pp.scale(adata, max_value=active_config.max_value, zero_center=True)

        elif active_config.scale_method == "robust":
            if scipy.sparse.issparse(adata.X):
                log.info("Using sparse-aware robust scaling")
                adata.X = _robust_scale_sparse(adata.X, max_value=active_config.max_value)
            else:
                adata.X = _robust_scale(adata.X, max_value=active_config.max_value)

        elif active_config.scale_method == "minmax":
            if scipy.sparse.issparse(adata.X):
                log.info("Using sparse-aware minmax scaling")
                adata.X = _minmax_scale_sparse(adata.X)
            else:
                adata.X = _minmax_scale(adata.X)

        else:
            raise ValueError(
                f"Unknown scale_method '{active_config.scale_method}'. "
                "Expected one of: zscore, robust, minmax."
            )

    except Exception as e:
        raise RuntimeError(f"[preprocess] Scaling failed: {e}. Check input data format.") from e

    # Backward compatibility: persist scaled matrix into a named layer.
    if output_layer:
        adata.layers[output_layer] = adata.X.copy()

    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["scaling"] = {
        "params": active_config.to_dict(),  # Pydantic's built-in serialization
        "output_layer": output_layer,
        "zero_center": active_config.scale_method == "zscore",
        "model_type": (
            "gaussian_consistent_mad_robust_zscore"
            if active_config.scale_method == "robust"
            else active_config.scale_method
        ),
        "claim_level": "standard_preprocessing",
        "mad_consistency_factor": (
            MAD_NORMAL_CONSISTENCY_FACTOR if active_config.scale_method == "robust" else None
        ),
        "review_note": (
            "Robust scaling uses median and MAD with the 1.4826 normal-consistency factor."
            if active_config.scale_method == "robust"
            else "Scaling is intended for PCA/graph construction, not expression-level interpretation."
        ),
    }

    log.info("Scaling complete. adata.X has been updated.")
    return adata


def regress_out(
    adata: AnnData,
    config: Optional[ScalingConfig] = None,
    input_layer: str = "normalized",
    output_layer: str = "regressed",
    **kwargs,
) -> AnnData:
    """
    Regress out unwanted sources of variation from gene expression data.

    Notes:
    - Expects input_layer to be log-normalized-like.
    - Variables must exist in adata.obs; typical covariates:
      ['total_counts', 'pct_counts_mt', 'S_score', 'G2M_score', 'cc_diff'].
    """
    if config is None:
        active_config = ScalingConfig()
    else:
        active_config = apply_config_overrides(config, **kwargs)

    vars_to_regress = list(active_config.vars_to_regress or [])
    if not vars_to_regress:
        log.info("No variables specified for regression. Skipping.")
        if input_layer != output_layer and input_layer in adata.layers:
            adata.layers[output_layer] = adata.layers[input_layer].copy()
        return adata

    missing_keys = [k for k in vars_to_regress if k not in adata.obs.columns]
    if missing_keys:
        log.warning(
            f"Variables to regress not found in adata.obs: {missing_keys}. "
            "Proceeding with available variables."
        )
        vars_to_regress = [k for k in vars_to_regress if k in adata.obs.columns]
        if not vars_to_regress:
            log.info("No valid variables left to regress. Skipping.")
            if input_layer != output_layer and input_layer in adata.layers:
                adata.layers[output_layer] = adata.layers[input_layer].copy()
            return adata

    if input_layer not in adata.layers:
        raise ValueError(f"Input layer '{input_layer}' not found in adata.layers.")

    cc_diagnostic = None
    if {"S_score", "G2M_score", "phase"} & set(vars_to_regress):
        cc_diagnostic = diagnose_cell_cycle_regression(adata)

    log.info(f"Regressing out: {', '.join(vars_to_regress)} from layer '{input_layer}'")
    temp_adata = AnnData(
        X=adata.layers[input_layer].copy(), obs=adata.obs.copy(), var=adata.var.copy()
    )
    try:
        sc.pp.regress_out(temp_adata, keys=vars_to_regress)
    except Exception as e:
        log.error(f"regress_out failed: {e}")
        raise

    adata.layers[output_layer] = temp_adata.X.copy()
    log.info(f"Regression complete. Results stored in adata.layers['{output_layer}'].")

    # Metadata
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["regress"] = {
        "input_layer": input_layer,
        "output_layer": output_layer,
        "vars_to_regress": vars_to_regress,
        "scanpy_version": version("scanpy"),
        "cell_cycle_regression_diagnostic": cc_diagnostic,
    }
    return adata


def plot_scaling_effect(
    adata: AnnData,
    original_data: Union[np.ndarray, scipy.sparse.spmatrix],
    scaled_layer: str = "scaled",
    n_genes: int = 5,
    gene_subset: Optional[List[str]] = None,
    save_dir: Optional[str] = None,
) -> plt.Figure:
    """
    Plot the effect of scaling on gene expression distributions.

    This function creates before/after distribution plots to visualize
    how scaling affects gene expression values.

    Args:
        adata: AnnData object with scaled data
        original_data: Original data before scaling
        scaled_layer: Layer containing scaled data
        n_genes: Number of top variable genes to plot
        gene_subset: Specific genes to plot instead of top variable
        save_dir: Directory to save the plot

    Returns:
        matplotlib Figure object
    """
    if scaled_layer not in adata.layers:
        raise ValueError(f"Scaled layer '{scaled_layer}' not found in adata.layers")

    if gene_subset is not None:
        genes_to_plot = [g for g in gene_subset if g in adata.var_names]
        if not genes_to_plot:
            raise ValueError("None of the specified genes were found in the data")
    else:
        # Find top variable genes in original data
        data = original_data.toarray() if scipy.sparse.issparse(original_data) else original_data
        gene_vars = np.var(data, axis=0)
        top_idx = np.argsort(-gene_vars)[:n_genes]
        genes_to_plot = adata.var_names[top_idx].tolist()

    fig, axes = plt.subplots(len(genes_to_plot), 2, figsize=(12, 3 * len(genes_to_plot)))
    fig.suptitle("Effect of Scaling on Gene Expression Distributions", fontsize=16)
    if len(genes_to_plot) == 1:
        axes = np.array([axes])

    def _get_gene_data(matrix, gene_name):
        idx = adata.var.index.get_loc(gene_name)
        data = matrix[:, idx]
        if scipy.sparse.issparse(data):
            return data.toarray().flatten()
        return data.flatten()

    for i, gene in enumerate(genes_to_plot):
        try:
            orig = _get_gene_data(original_data, gene)
            scaled = _get_gene_data(adata.layers[scaled_layer], gene)

            sns.histplot(orig, bins=30, kde=True, ax=axes[i, 0])
            axes[i, 0].set_title(f"{gene} - Before Scaling")
            axes[i, 0].text(
                0.05,
                0.95,
                f"Mean: {np.mean(orig):.2f}\nStd: {np.std(orig):.2f}",
                transform=axes[i, 0].transAxes,
                va="top",
            )
            sns.histplot(scaled, bins=30, kde=True, ax=axes[i, 1])
            axes[i, 1].set_title(f"{gene} - After Scaling")
            axes[i, 1].text(
                0.05,
                0.95,
                f"Mean: {np.mean(scaled):.2f}\nStd: {np.std(scaled):.2f}",
                transform=axes[i, 1].transAxes,
                va="top",
            )
        except Exception as e:
            log.warning(f"Failed to plot gene {gene}: {str(e)}")
            axes[i, 0].text(0.5, 0.5, f"Error plotting {gene}", ha="center", va="center")
            axes[i, 1].text(0.5, 0.5, f"Error plotting {gene}", ha="center", va="center")

    plt.tight_layout(rect=[0, 0.03, 1, 0.97])

    if save_dir:
        save_path = Path(save_dir)
        save_path.mkdir(parents=True, exist_ok=True)
        figure_path = save_path / "scaling_effect.png"
        plt.savefig(figure_path, dpi=300)
        log.info(f"Saved scaling effect plot to {figure_path}")

    return fig
