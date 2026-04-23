"""
Gibbs sampling for BayesPrism (R-free)

Optimized Gibbs sampler for Bayesian deconvolution.
"""

import logging
from typing import Optional, Tuple

import numpy as np
from numba import njit, prange

log = logging.getLogger(__name__)


@njit(parallel=True, cache=True)
def _sample_Z_numba(
    mixture: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
    Z: np.ndarray,
) -> np.ndarray:
    """
    Numba-optimized Z sampling

    Samples cell type-specific expression Z from multinomial distribution.
    Parallelized over genes for performance.

    Parameters
    ----------
    mixture : np.ndarray
        Bulk expression vector (n_genes,)
    theta : np.ndarray
        Cell type proportions (n_cell_types,)
    phi : np.ndarray
        Reference profile (n_genes, n_cell_types)
    Z : np.ndarray
        Output array (n_genes, n_cell_types)

    Returns:
    -------
    np.ndarray
        Updated Z matrix
    """
    n_genes = mixture.shape[0]
    n_cell_types = theta.shape[0]

    for g in prange(n_genes):
        if mixture[g] <= 0:
            Z[g, :] = 0
            continue

        # Compute probabilities for each cell type
        probs = np.empty(n_cell_types)
        prob_sum = 0.0

        for k in range(n_cell_types):
            p = phi[g, k] * theta[k]
            probs[k] = p
            prob_sum += p

        # Normalize
        if prob_sum > 0:
            for k in range(n_cell_types):
                probs[k] /= prob_sum
        else:
            probs[:] = 1.0 / n_cell_types

        # Multinomial sampling
        # For efficiency, use cumulative probabilities
        cumsum = np.empty(n_cell_types)
        cumsum[0] = probs[0]
        for k in range(1, n_cell_types):
            cumsum[k] = cumsum[k - 1] + probs[k]

        # Draw samples
        n_reads = int(mixture[g])
        Z[g, :] = 0

        for _ in range(n_reads):
            u = np.random.random()
            for k in range(n_cell_types):
                if u <= cumsum[k]:
                    Z[g, k] += 1
                    break

    return Z


@njit(cache=True)
def _sample_theta_numba(
    Z: np.ndarray,
    alpha_prior: np.ndarray,
) -> np.ndarray:
    """
    Numba-optimized theta sampling from Dirichlet

    Parameters
    ----------
    Z : np.ndarray
        Cell type expression (n_genes, n_cell_types)
    alpha_prior : np.ndarray
        Dirichlet prior parameters (n_cell_types,)

    Returns:
    -------
    np.ndarray
        Sampled theta proportions
    """
    n_cell_types = Z.shape[1]

    # Sum expression per cell type
    cell_type_counts = np.zeros(n_cell_types)
    for k in range(n_cell_types):
        s = 0.0
        for g in range(Z.shape[0]):
            s += Z[g, k]
        cell_type_counts[k] = s

    # Posterior parameters
    alpha_posterior = alpha_prior + cell_type_counts

    # Sample from Dirichlet using Gamma samples
    theta = np.empty(n_cell_types)
    theta_sum = 0.0

    for k in range(n_cell_types):
        if alpha_posterior[k] > 0:
            # Gamma sampling
            theta[k] = np.random.gamma(alpha_posterior[k], 1.0)
        else:
            theta[k] = 0.0
        theta_sum += theta[k]

    # Normalize
    if theta_sum > 0:
        for k in range(n_cell_types):
            theta[k] /= theta_sum
    else:
        # Uniform if all zero
        for k in range(n_cell_types):
            theta[k] = 1.0 / n_cell_types

    return theta


class GibbsSampler:
    """
    Gibbs sampler for BayesPrism deconvolution

    Implements blocked Gibbs sampling to estimate cell type proportions
    and cell type-specific expression.

    Parameters
    ----------
    n_iter : int
        Number of sampling iterations
    burnin : int
        Number of burn-in iterations
    thinning : int
        Thinning interval for samples
    use_numba : bool
        Whether to use Numba JIT compilation

    Attributes:
    ----------
    theta_samples_ : np.ndarray
        Posterior samples of theta (n_iter - burnin, n_cell_types)
    Z_samples_ : np.ndarray
        Posterior samples of Z

    Examples:
    --------
    >>> sampler = GibbsSampler(n_iter=100, burnin=50)
    >>> theta_samples, Z_samples = sampler.sample(
    ...     mixture=bulk_expr,
    ...     theta_init=theta0,
    ...     phi=ref_phi,
    ... )
    """

    def __init__(
        self,
        n_iter: int = 100,
        burnin: int = 50,
        thinning: int = 1,
        use_numba: bool = True,
    ):
        self.n_iter = n_iter
        self.burnin = burnin
        self.thinning = thinning
        self.use_numba = use_numba

        # Storage for samples
        self.theta_samples_ = None
        self.Z_samples_ = None

    def sample(
        self,
        mixture: np.ndarray,
        theta_init: np.ndarray,
        phi: np.ndarray,
        alpha_prior: Optional[np.ndarray] = None,
        verbose: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Run Gibbs sampling

        Parameters
        ----------
        mixture : np.ndarray
            Bulk expression (n_genes,)
        theta_init : np.ndarray
            Initial cell type proportions
        phi : np.ndarray
            Reference profile (n_genes, n_cell_types)
        alpha_prior : np.ndarray, optional
            Dirichlet prior parameters (default: uniform)
        verbose : bool
            Whether to print progress

        Returns:
        -------
        theta_samples : np.ndarray
            Posterior samples of theta
        Z_samples : np.ndarray
            Posterior samples of Z
        """
        n_genes, n_cell_types = phi.shape

        # Initialize
        theta = theta_init.copy()
        Z = np.zeros((n_genes, n_cell_types))

        if alpha_prior is None:
            alpha_prior = np.ones(n_cell_types)

        # Calculate number of samples to store
        n_store = (self.n_iter - self.burnin) // self.thinning
        if n_store <= 0:
            n_store = 1

        theta_samples = np.zeros((n_store, n_cell_types))
        Z_samples = np.zeros((n_store, n_genes, n_cell_types))

        sample_idx = 0

        # Run Gibbs sampling
        for iteration in range(self.n_iter + self.burnin):
            # Sample Z | theta, X
            if self.use_numba:
                Z = _sample_Z_numba(mixture, theta, phi, Z)
            else:
                Z = self._sample_Z_python(mixture, theta, phi)

            # Sample theta | Z
            if self.use_numba:
                theta = _sample_theta_numba(Z, alpha_prior)
            else:
                theta = self._sample_theta_python(Z, alpha_prior)

            # Store samples after burnin
            if iteration >= self.burnin and (iteration - self.burnin) % self.thinning == 0:
                if sample_idx < n_store:
                    theta_samples[sample_idx] = theta
                    Z_samples[sample_idx] = Z
                    sample_idx += 1

            if verbose and iteration % 50 == 0:
                log.info(f"Gibbs iteration {iteration}/{self.n_iter + self.burnin}")

        self.theta_samples_ = theta_samples
        self.Z_samples_ = Z_samples

        return theta_samples, Z_samples

    def _sample_Z_python(
        self,
        mixture: np.ndarray,
        theta: np.ndarray,
        phi: np.ndarray,
    ) -> np.ndarray:
        """Python fallback for Z sampling"""
        n_genes, n_cell_types = phi.shape
        Z = np.zeros((n_genes, n_cell_types))

        for g in range(n_genes):
            if mixture[g] == 0:
                continue

            probs = phi[g, :] * theta
            probs = probs / (probs.sum() + 1e-10)

            Z[g, :] = np.random.multinomial(int(mixture[g]), probs)

        return Z

    def _sample_theta_python(
        self,
        Z: np.ndarray,
        alpha_prior: np.ndarray,
    ) -> np.ndarray:
        """Python fallback for theta sampling"""
        cell_type_counts = Z.sum(axis=0)
        alpha_posterior = alpha_prior + cell_type_counts

        return np.random.dirichlet(alpha_posterior)

    def get_posterior_mean(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get posterior mean estimates

        Returns:
        -------
        theta_mean : np.ndarray
            Mean cell type proportions
        Z_mean : np.ndarray
            Mean cell type expression
        """
        if self.theta_samples_ is None:
            raise ValueError("Must run sample() first")

        return (
            self.theta_samples_.mean(axis=0),
            self.Z_samples_.mean(axis=0),
        )

    def get_credible_intervals(
        self,
        level: float = 0.95,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get credible intervals for theta

        Parameters
        ----------
        level : float
            Confidence level

        Returns:
        -------
        lower : np.ndarray
            Lower bounds
        upper : np.ndarray
            Upper bounds
        """
        if self.theta_samples_ is None:
            raise ValueError("Must run sample() first")

        alpha = (1 - level) / 2
        lower = np.percentile(self.theta_samples_, alpha * 100, axis=0)
        upper = np.percentile(self.theta_samples_, (1 - alpha) * 100, axis=0)

        return lower, upper
