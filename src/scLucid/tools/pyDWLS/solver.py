"""
DWLS solver for bulk RNA-seq deconvolution.

Implements the dampened weighted least squares optimization of Tsoucas et al.
(2019, Cell Reports) for inferring cell-type proportions from a bulk
expression vector and a single-cell-derived signature matrix.

The classical WLS step minimizes ``|| W (b - S theta) ||^2`` with
``w_i = 1 / mu_i^2`` (``mu = S theta``); the dampening cap controls how
aggressively highly expressed genes are downweighted, preventing a handful of
saturating transcripts from dominating the fit.
"""

import logging
from typing import Tuple

import numpy as np
from scipy.optimize import nnls

log = logging.getLogger(__name__)


def solve_nnls(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Solve a non-negative least-squares problem ``min ||A x - b||`` subject to
    ``x >= 0`` and return only the solution vector.

    Thin convenience wrapper around :func:`scipy.optimize.nnls`.

    Parameters
    ----------
    A : np.ndarray
        Design matrix of shape ``(m, n)``.
    b : np.ndarray
        Target vector of length ``m``.

    Returns:
    -------
    np.ndarray
        Non-negative solution vector of length ``n``.

    Examples:
    --------
    >>> x = solve_nnls(signature, bulk_sample)
    """
    x, _residual = nnls(A, b)
    return x


def _normalize_to_simplex(theta: np.ndarray) -> np.ndarray:
    """Project a non-negative vector onto the unit simplex."""
    total = float(theta.sum())
    if total <= 0.0:
        # All-zero degenerate case: spread mass uniformly so downstream
        # consumers still receive a valid proportion vector.
        n = theta.shape[0]
        return np.full(n, 1.0 / n) if n > 0 else theta
    return theta / total


class DampenedWLS:
    """
    Dampened Weighted Least Squares solver for cell-type proportion inference.

    The solver alternates between (a) computing predicted expression
    ``mu = S theta`` for the current proportions and (b) re-solving NNLS with
    weights ``w = 1 / mu^2`` capped at a quantile controlled by
    ``dampen_factor``. The cap is what gives DWLS its robustness against
    highly expressed marker genes.

    Parameters
    ----------
    dampen_factor : float, default=1.0
        Controls how aggressively raw WLS weights are capped before reweighting.
        ``0.0`` disables iterative reweighting entirely and returns the
        ordinary NNLS estimate. Larger values flatten the weight distribution
        (cap at a lower quantile), making the fit closer to ordinary least
        squares; typical R-package defaults sit near ``1.0``.
    use_nonneg : bool, default=True
        Enforce non-negative proportions. The current implementation always
        uses NNLS; this flag is reserved for a future unconstrained fallback.
    max_iter : int, default=10
        Maximum number of reweighting iterations.
    tol : float, default=1e-5
        L2 convergence tolerance on the proportion vector between iterations.

    Examples:
    --------
    >>> solver = DampenedWLS(dampen_factor=1.0)
    >>> proportions = solver.solve(signature_matrix, bulk_sample)
    >>> assert np.isclose(proportions.sum(), 1.0)
    """

    def __init__(
        self,
        dampen_factor: float = 1.0,
        use_nonneg: bool = True,
        max_iter: int = 10,
        tol: float = 1e-5,
    ):
        if dampen_factor < 0:
            raise ValueError(f"dampen_factor must be >= 0, got {dampen_factor}")
        if max_iter < 1:
            raise ValueError(f"max_iter must be >= 1, got {max_iter}")

        self.dampen_factor = float(dampen_factor)
        self.use_nonneg = bool(use_nonneg)
        self.max_iter = int(max_iter)
        self.tol = float(tol)

    def solve(self, signature: np.ndarray, bulk: np.ndarray) -> np.ndarray:
        """
        Estimate cell-type proportions for a single bulk sample.

        Parameters
        ----------
        signature : np.ndarray
            Signature matrix of shape ``(n_genes, n_cell_types)``.
        bulk : np.ndarray
            Bulk expression vector of length ``n_genes``.

        Returns:
        -------
        np.ndarray
            Non-negative proportion vector of length ``n_cell_types`` summing
            to 1.

        Raises:
        ------
        ValueError
            If the input shapes are incompatible.
        """
        signature = np.asarray(signature, dtype=float)
        bulk = np.asarray(bulk, dtype=float).ravel()
        if signature.ndim != 2:
            raise ValueError(
                f"signature must be 2-D, got shape {signature.shape}"
            )
        if bulk.shape[0] != signature.shape[0]:
            raise ValueError(
                f"signature has {signature.shape[0]} genes but bulk has "
                f"{bulk.shape[0]}"
            )

        theta = solve_nnls(signature, bulk)

        if self.dampen_factor == 0.0:
            return _normalize_to_simplex(theta)

        for it in range(self.max_iter):
            predicted = signature @ theta
            predicted = np.maximum(predicted, 1e-8)

            raw_weights = 1.0 / np.square(predicted)
            cap = self._weight_cap(raw_weights)
            weights = np.minimum(raw_weights, cap)

            sqrt_w = np.sqrt(weights)
            A_w = signature * sqrt_w[:, np.newaxis]
            b_w = bulk * sqrt_w

            theta_new = solve_nnls(A_w, b_w)

            shift = float(np.linalg.norm(theta_new - theta))
            theta = theta_new
            log.debug(
                "DampenedWLS iter %d: shift=%.3e, theta_sum=%.3f",
                it,
                shift,
                float(theta.sum()),
            )
            if shift < self.tol:
                break

            if theta.sum() == 0.0:
                break

        return _normalize_to_simplex(theta)

    def _weight_cap(self, raw_weights: np.ndarray) -> float:
        """
        Compute the upper bound on per-gene weights.

        Larger ``dampen_factor`` → smaller quantile → tighter cap → flatter
        weight distribution. ``dampen_factor=1`` caps at the 90th percentile;
        ``dampen_factor=2`` caps at the 95th percentile of inverse-mu^2 values.
        """
        quantile = 1.0 - 0.1 / max(self.dampen_factor, 1e-3)
        quantile = float(np.clip(quantile, 0.5, 0.999))
        return float(np.quantile(raw_weights, quantile))


__all__ = [
    "DampenedWLS",
    "solve_nnls",
]
