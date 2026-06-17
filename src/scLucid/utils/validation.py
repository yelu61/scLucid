"""
Validation utilities for scLucid.

This module provides comprehensive validation functions for:
- AnnData object structure and contents
- Configuration objects
- Analysis results
- Layer consistency
"""

import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp
from anndata import AnnData

from .contracts import ContractError, StageName, validate_stage_contract

log = logging.getLogger(__name__)


class ValidationError(Exception):
    """Exception raised for validation errors."""

    def __init__(self, message: str, field: Optional[str] = None):
        super().__init__(message)
        self.field = field


def validate_adata(
    adata: AnnData,
    required_layers: Optional[List[str]] = None,
    required_obs: Optional[List[str]] = None,
    required_var: Optional[List[str]] = None,
    required_obsm: Optional[List[str]] = None,
    required_uns: Optional[List[str]] = None,
    check_counts: bool = False,
    check_normalized: bool = False,
    raise_on_error: bool = True,
) -> Dict[str, Any]:
    """
    Validate AnnData object structure and contents.

    Args:
        adata: AnnData object to validate
        required_layers: List of required layer names
        required_obs: List of required obs column names
        required_var: List of required var column names
        required_obsm: List of required obsm key names
        required_uns: List of required uns key names
        check_counts: Whether to verify raw counts exist
        check_normalized: Whether to verify normalized data exists
        raise_on_error: If True, raise ValidationError on failure. If False, return dict with results.

    Returns:
        Dictionary with validation results:
        {
            "valid": bool,
            "errors": List[str],
            "warnings": List[str],
            "layer_info": Dict,
            "shape": Tuple[int, int]
        }

    Examples:
        >>> # Validate QC input
        >>> validate_adata(adata, required_layers=["counts"])

        >>> # Validate preprocessing input
        >>> validate_adata(
        ...     adata,
        ...     required_layers=["counts"],
        ...     check_counts=True
        ... )

        >>> # Validate analysis input
        >>> validate_adata(
        ...     adata,
        ...     required_obsm=["X_pca"],
        ...     required_obs=["leiden"]
        ... )
    """
    errors = []
    warnings = []
    layer_info = {}

    # Check basic structure
    if adata is None:
        errors.append("AnnData object is None")
        return _make_result(False, errors, warnings, layer_info, (0, 0))

    # Check dimensions
    if adata.n_obs == 0:
        errors.append("AnnData has 0 cells (n_obs=0)")
    if adata.n_vars == 0:
        errors.append("AnnData has 0 genes (n_vars=0)")

    # Check required layers
    if required_layers:
        for layer in required_layers:
            if layer not in adata.layers:
                errors.append(f"Required layer '{layer}' not found")
            else:
                layer_data = adata.layers[layer]
                layer_info[layer] = {
                    "shape": layer_data.shape,
                    "dtype": str(layer_data.dtype),
                    "sparse": hasattr(layer_data, "nnz"),
                }

    # Check required obs columns
    if required_obs:
        for col in required_obs:
            if col not in adata.obs.columns:
                errors.append(f"Required obs column '{col}' not found")

    # Check required var columns
    if required_var:
        for col in required_var:
            if col not in adata.var.columns:
                errors.append(f"Required var column '{col}' not found")

    # Check required obsm keys
    if required_obsm:
        for key in required_obsm:
            if key not in adata.obsm:
                errors.append(f"Required obsm key '{key}' not found")
            else:
                obsm_data = adata.obsm[key]
                if obsm_data.shape[0] != adata.n_obs:
                    errors.append(
                        f"obsm['{key}'] has wrong shape: {obsm_data.shape[0]} rows, "
                        f"expected {adata.n_obs}"
                    )

    # Check required uns keys
    if required_uns:
        for key in required_uns:
            if key not in adata.uns:
                errors.append(f"Required uns key '{key}' not found")

    # Check counts layer
    if check_counts:
        if "counts" not in adata.layers and not _is_counts_matrix(adata.X):
            errors.append("No counts data found. Expected 'counts' layer or integer X matrix")
        elif "counts" in adata.layers:
            if not _is_counts_matrix(adata.layers["counts"]):
                warnings.append(
                    "'counts' layer may not contain raw counts (contains non-integer values)"
                )

    # Check normalized layer
    if check_normalized:
        # AnnData always exposes a `.raw` attribute (it may be None), so a
        # plain hasattr check is not informative. Detect either a populated
        # raw slot or a normalized layer.
        has_normalized_layer = "normalized" in adata.layers
        has_raw = getattr(adata, "raw", None) is not None
        if not has_normalized_layer and not has_raw:
            warnings.append(
                "No normalized data found. Expected 'normalized' layer or .raw attribute"
            )

    # Additional sanity checks
    if adata.obs_names.duplicated().any():
        n_dup = adata.obs_names.duplicated().sum()
        errors.append(f"Found {n_dup} duplicate cell names in obs_names")

    if adata.var_names.duplicated().any():
        n_dup = adata.var_names.duplicated().sum()
        errors.append(f"Found {n_dup} duplicate gene names in var_names")

    result = _make_result(
        valid=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        layer_info=layer_info,
        shape=(adata.n_obs, adata.n_vars),
    )

    if raise_on_error and errors:
        raise ValidationError(f"Validation failed: {'; '.join(errors)}")

    return result


def validate_config(
    config: Any,
    required_fields: Optional[List[str]] = None,
    field_types: Optional[Dict[str, type]] = None,
    raise_on_error: bool = True,
) -> Dict[str, Any]:
    """
    Validate configuration object.

    Args:
        config: Configuration object to validate
        required_fields: List of required field names
        field_types: Dict mapping field names to expected types
        raise_on_error: If True, raise ValidationError on failure

    Returns:
        Dictionary with validation results

    Examples:
        >>> validate_config(
        ...     config,
        ...     required_fields=["method", "n_top_genes"],
        ...     field_types={"n_top_genes": int}
        ... )
    """
    errors = []

    if config is None:
        errors.append("Config is None")
        result = {"valid": False, "errors": errors, "warnings": []}
        if raise_on_error:
            raise ValidationError("Config is None")
        return result

    # Check required fields
    if required_fields:
        for field in required_fields:
            if not hasattr(config, field):
                errors.append(f"Required field '{field}' not found in config")
            elif getattr(config, field) is None:
                errors.append(f"Required field '{field}' is None")

    # Check field types
    if field_types:
        for field, expected_type in field_types.items():
            if hasattr(config, field):
                value = getattr(config, field)
                if value is not None and not isinstance(value, expected_type):
                    errors.append(
                        f"Field '{field}' has wrong type: {type(value).__name__}, "
                        f"expected {expected_type.__name__}"
                    )

    result = {"valid": len(errors) == 0, "errors": errors, "warnings": []}

    if raise_on_error and errors:
        raise ValidationError(f"Config validation failed: {'; '.join(errors)}")

    return result


def validate_analysis_results(
    adata: AnnData,
    analysis_type: str,
    raise_on_error: bool = True,
) -> Dict[str, Any]:
    """
    Validate that required analysis results exist in AnnData.

    Args:
        adata: AnnData object
        analysis_type: Type of analysis ("qc", "preprocess", "clustering", "annotation")
        raise_on_error: If True, raise ValidationError on failure

    Returns:
        Dictionary with validation results

    Examples:
        >>> # Validate QC was run
        >>> validate_analysis_results(adata, "qc")

        >>> # Validate preprocessing was run
        >>> validate_analysis_results(adata, "preprocess")

        >>> # Validate clustering was run
        >>> validate_analysis_results(adata, "clustering")
    """
    errors = []
    warnings = []

    if "sclucid" not in adata.uns:
        errors.append("No sclucid results found. Run analysis first.")
        result = {"valid": False, "errors": errors, "warnings": warnings}
        if raise_on_error:
            raise ValidationError("No sclucid results found")
        return result

    sclucid_data = adata.uns["sclucid"]

    if analysis_type == "qc":
        if "qc" not in sclucid_data:
            errors.append("QC results not found")
        else:
            required_qc_cols = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
            missing_cols = [c for c in required_qc_cols if c not in adata.obs.columns]
            if missing_cols:
                errors.append(f"Missing QC metric columns: {missing_cols}")

    elif analysis_type == "preprocess":
        if "preprocess" not in sclucid_data:
            errors.append("Preprocessing results not found")
        else:
            if "normalized" not in adata.layers:
                errors.append("'normalized' layer not found")
            if "X_pca" not in adata.obsm:
                errors.append("PCA results (obsm['X_pca']) not found")

    elif analysis_type == "clustering":
        if "leiden" not in adata.obs.columns and "louvain" not in adata.obs.columns:
            errors.append("Clustering results not found (no 'leiden' or 'louvain' in obs)")

    elif analysis_type == "annotation":
        if "cell_type" not in adata.obs.columns:
            errors.append("Cell type annotations not found (no 'cell_type' in obs)")

    elif analysis_type == "markers":
        if "rank_genes_groups" not in adata.uns:
            errors.append("Marker gene results not found (uns['rank_genes_groups'])")

    else:
        warnings.append(f"Unknown analysis type: {analysis_type}")

    result = {"valid": len(errors) == 0, "errors": errors, "warnings": warnings}

    if raise_on_error and errors:
        raise ValidationError(f"Analysis validation failed: {'; '.join(errors)}")

    return result


def check_layer_consistency(adata: AnnData, layers: Optional[List[str]] = None) -> Dict[str, Any]:
    """
    Check consistency between different layers.

    Args:
        adata: AnnData object
        layers: List of layers to check. If None, checks all layers.

    Returns:
        Dictionary with consistency check results

    Examples:
        >>> check_layer_consistency(adata, ["counts", "normalized", "scaled"])
    """
    if layers is None:
        layers = list(adata.layers.keys())

    errors = []
    warnings = []
    shapes = {}

    base_shape = (adata.n_obs, adata.n_vars)

    for layer in layers:
        if layer not in adata.layers:
            errors.append(f"Layer '{layer}' not found")
            continue

        layer_data = adata.layers[layer]
        layer_shape = layer_data.shape
        shapes[layer] = layer_shape

        if layer_shape != base_shape:
            errors.append(
                f"Layer '{layer}' has inconsistent shape: {layer_shape}, expected {base_shape}"
            )

    # Check for all-NaN or all-zero layers
    for layer in layers:
        if layer not in adata.layers:
            continue
        layer_data = adata.layers[layer]

        if hasattr(layer_data, "nnz"):  # sparse
            if layer_data.nnz == 0:
                warnings.append(f"Layer '{layer}' is all zeros (sparse)")
        else:  # dense
            if np.all(layer_data == 0):
                warnings.append(f"Layer '{layer}' is all zeros")
            elif np.all(np.isnan(layer_data)):
                errors.append(f"Layer '{layer}' is all NaN")

    return {
        "consistent": len(errors) == 0,
        "errors": errors,
        "warnings": warnings,
        "shapes": shapes,
        "expected_shape": base_shape,
    }


# ============================================================================
# Internal helper functions
# ============================================================================


def is_raw_count_matrix(
    X,
    *,
    max_cells: int = 2000,
    max_genes: int = 2000,
    zero_fraction_threshold: float = 0.20,
    max_fractional_rate: float = 0.01,
    min_max_value: float = 10.0,
) -> Tuple[bool, Dict[str, Any]]:
    """Lightweight check for raw-count-like matrix semantics.

    Subsamples large matrices to avoid OOM. Returns a diagnostic dict so callers
    can explain why a matrix was rejected.

    Parameters
    ----------
    X
        Array-like or sparse matrix. Typically ``adata.X`` or
        ``adata.layers['counts']``.
    max_cells, max_genes
        Subsample dimensions for the check. Larger values are more accurate but
        slower and more memory-hungry.
    zero_fraction_threshold
        Minimum fraction of zeros expected in a sparse count matrix.
    max_fractional_rate
        Maximum fraction of positive values that are non-integer. Count matrices
        may contain a small amount of rounding noise from upstream tools; values
        above this threshold are treated as suspicious.
    min_max_value
        Minimum maximum value expected for a real count matrix. Very small max
        values suggest the matrix is already normalized or scaled.

    Returns:
    -------
    (is_raw_like, diagnostics)
        ``is_raw_like`` is True when the matrix looks like raw counts.
        ``diagnostics`` contains the supporting numbers.

    Examples:
    --------
    >>> is_raw, info = is_raw_count_matrix(adata.layers["counts"])
    >>> if not is_raw:
    ...     print(f"Warning: counts layer looks suspicious: {info}")
    """
    if sp.issparse(X):
        n_cells, n_genes = X.shape
        sub_cells = min(n_cells, max_cells)
        sub_genes = min(n_genes, max_genes)
        Xs = X[:sub_cells, :sub_genes].copy()
        data = Xs.data
        zero_fraction = 1.0 - (Xs.nnz / max(1, sub_cells * sub_genes))
    else:
        arr = np.asarray(X)
        n_cells, n_genes = arr.shape
        sub_cells = min(n_cells, max_cells)
        sub_genes = min(n_genes, max_genes)
        Xs = arr[:sub_cells, :sub_genes]
        data = Xs.ravel()
        zero_fraction = float(np.mean(Xs == 0))

    data = np.asarray(data)
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        return False, {"reason": "no_finite_values"}

    has_negative = bool(np.any(finite < 0))
    positive = finite[finite > 0]
    fractional_rate = (
        float(np.mean(np.abs(positive - np.round(positive)) > 1e-6))
        if positive.size
        else 0.0
    )
    max_value = float(np.max(finite))

    raw_like = (
        (not has_negative)
        and fractional_rate < max_fractional_rate
        and zero_fraction >= zero_fraction_threshold
        and max_value > min_max_value
    )

    diagnostics = {
        "has_negative": has_negative,
        "fractional_positive_rate": fractional_rate,
        "zero_fraction": float(zero_fraction),
        "max_value": max_value,
        "subsampled_shape": (sub_cells, sub_genes),
    }
    return raw_like, diagnostics


def _is_counts_matrix(data) -> bool:
    """Internal bool compatibility shim around ``is_raw_count_matrix``."""
    try:
        raw_like, _ = is_raw_count_matrix(
            data,
            max_cells=64,
            max_genes=512,
            zero_fraction_threshold=0.0,
            min_max_value=0.0,
        )
        return bool(raw_like)
    except Exception:
        return False


def _make_result(
    valid: bool, errors: List[str], warnings: List[str], layer_info: Dict, shape: tuple
) -> Dict[str, Any]:
    """Create standardized validation result dict."""
    return {
        "valid": valid,
        "errors": errors,
        "warnings": warnings,
        "layer_info": layer_info,
        "shape": shape,
    }


# ============================================================================
# Convenience validation functions
# ============================================================================


def assert_qc_ready(adata: AnnData) -> None:
    """Assert that adata is ready for QC (has raw counts)."""
    validate_adata(adata, check_counts=True, raise_on_error=True)
    try:
        validate_stage_contract(adata, "qc", when="input", raise_on_error=True)
    except ContractError as exc:
        raise ValidationError(str(exc)) from exc


def assert_preprocessing_ready(adata: AnnData) -> None:
    """Assert that adata has completed QC and is ready for preprocessing."""
    try:
        validate_stage_contract(adata, "preprocess", when="input", raise_on_error=True)
    except ContractError as exc:
        raise ValidationError(str(exc)) from exc


def assert_analysis_ready(adata: AnnData) -> None:
    """Assert that adata is ready for analysis (clustering, annotation)."""
    try:
        validate_stage_contract(adata, "analysis", when="input", raise_on_error=True)
    except ContractError as exc:
        raise ValidationError(str(exc)) from exc


def validate_workflow_contract(
    adata: AnnData,
    stage: StageName,
    *,
    when: str,
    raise_on_error: bool = True,
) -> Dict[str, Any]:
    """Validate one workflow stage against the canonical scLucid contract."""
    result = validate_stage_contract(
        adata,
        stage,
        when=when,  # type: ignore[arg-type]
        raise_on_error=False,
    )
    if raise_on_error and not result.valid:
        raise ValidationError(str(ContractError(result)))
    return result.to_dict()


__all__ = [
    "ValidationError",
    "validate_adata",
    "validate_config",
    "validate_analysis_results",
    "check_layer_consistency",
    "is_raw_count_matrix",
    "assert_qc_ready",
    "assert_preprocessing_ready",
    "assert_analysis_ready",
    "validate_workflow_contract",
]
