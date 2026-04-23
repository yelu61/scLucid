"""
Shared preprocessing utilities.
"""

import logging
from typing import Union

import numpy as np
import scipy.sparse

log = logging.getLogger(__name__)


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
