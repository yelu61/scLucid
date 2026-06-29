"""Shared preprocessing utilities."""

import logging
from typing import Tuple, Union

from anndata import AnnData
import numpy as np
import scipy.sparse

log = logging.getLogger(__name__)

__all__ = [
    "validate_matrix_input",
    "resolve_input_matrix",
    "apply_row_scale",
    "apply_log1p",
]


def validate_matrix_input(
    data: Union[np.ndarray, scipy.sparse.spmatrix],
    name: str = "input",
    *,
    allow_negative: bool = True,
) -> None:
    """
    Validate a matrix before preprocessing.

    Checks:
    - Non-empty shape
    - Finite values (no NaN or Inf)
    - Non-negative values (optional, default allows negatives)

    Args:
        data: Input matrix
        name: Descriptive name for error messages
        allow_negative: If False, raises on negative values

    Raises:
        ValueError: On validation failure
    """
    if data.shape[0] == 0 or data.shape[1] == 0:
        raise ValueError(f"{name} is empty with shape {data.shape}.")

    if scipy.sparse.issparse(data):
        values = data.data
        min_val = data.min() if data.nnz > 0 else 0.0
    else:
        values = np.asarray(data)
        min_val = np.min(values)

    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} contains NaN or Inf values.")

    if not allow_negative and min_val < 0:
        raise ValueError(f"{name} contains negative values. Use raw non-negative counts as input.")


def resolve_input_matrix(
    adata: AnnData,
    input_layer: str,
) -> Tuple[Union[np.ndarray, scipy.sparse.spmatrix], str]:
    """Resolve the configured input matrix with graceful counts fallback."""
    if input_layer == "X":
        return adata.X, "adata.X"
    if input_layer in adata.layers:
        return adata.layers[input_layer], f"adata.layers['{input_layer}']"

    if input_layer == "counts":
        log.warning(
            "Layer 'counts' not found. Falling back to adata.X for normalization. "
            "Consider creating adata.layers['counts'] = adata.X.copy() for reproducibility."
        )
        return adata.X, "adata.X (fallback from missing counts layer)"

    available = list(adata.layers.keys())
    raise ValueError(
        f"Input layer '{input_layer}' not found. Available layers: {available or '[]'}."
    )


def apply_row_scale(
    matrix: Union[np.ndarray, scipy.sparse.spmatrix],
    row_scale: np.ndarray,
) -> Union[np.ndarray, scipy.sparse.spmatrix]:
    """Scale each row of a matrix by the provided factors."""
    row_scale = np.asarray(row_scale, dtype=float).ravel()
    if scipy.sparse.issparse(matrix):
        scaled = matrix.tocsr(copy=True)
        scaled.data *= np.repeat(row_scale, np.diff(scaled.indptr))
        return scaled
    return np.asarray(matrix, dtype=float) * row_scale[:, np.newaxis]


def apply_log1p(
    matrix: Union[np.ndarray, scipy.sparse.spmatrix],
) -> Union[np.ndarray, scipy.sparse.spmatrix]:
    """Apply log1p while preserving sparse structure when possible."""
    if scipy.sparse.issparse(matrix):
        transformed = matrix.copy()
        transformed.data = np.log1p(transformed.data)
        transformed.eliminate_zeros()
        return transformed
    return np.log1p(np.asarray(matrix, dtype=float))
