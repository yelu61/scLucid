"""
Adaptive threshold learning for QC metrics.

This module provides machine learning-based approaches to automatically learn
optimal QC thresholds from data distributions.
"""

import logging
import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from anndata import AnnData
from scipy import optimize, stats
from scipy.special import gammaln
from sklearn.cluster import DBSCAN
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors

from ..runtime import effective_n_jobs

log = logging.getLogger(__name__)

#: Scale factor to convert MAD to approximate standard deviation for a normal distribution.
#: For normally distributed data, std ≈ MAD * 1.4826.
MAD_SCALE_FACTOR: float = 1.4826

#: Maximum number of observations used for k-distance estimation in DBSCAN eps selection.
#: This avoids O(N^2) k-NN computations on very large datasets.
KDISTANCE_MAX_SAMPLES: int = 20000


def _nb_cdf(y: np.ndarray, mu: float, alpha: float) -> np.ndarray:
    """CDF of the negative binomial distribution at integer points.

    Parameterized via mean ``mu`` and overdispersion ``alpha`` where
    variance = mu + alpha * mu^2.
    """
    if mu <= 0 or alpha <= 0:
        return np.zeros_like(y, dtype=float)
    if alpha < 1e-8:
        return stats.poisson.cdf(y, mu)
    p = 1.0 / (1.0 + alpha * mu)
    r = 1.0 / alpha
    # scipy nbinom uses n=r, p=p
    return stats.nbinom.cdf(y, r, p)


def _zinb_cdf(y: np.ndarray, pi: float, mu: float, alpha: float) -> np.ndarray:
    """CDF of the zero-inflated negative binomial distribution.

    ``F(y) = pi + (1 - pi) * F_nb(y)`` for y >= 0.  At y < 0 the CDF is 0.
    """
    pi = float(np.clip(pi, 1e-12, 1.0 - 1e-12))
    y = np.asarray(y, dtype=float)
    result = np.empty_like(y, dtype=float)
    result[y < 0] = 0.0
    nonneg = y >= 0
    result[nonneg] = pi + (1.0 - pi) * _nb_cdf(y[nonneg], mu, alpha)
    return result


def _zinb_threshold(
    percentile: float,
    pi: float,
    mu: float,
    alpha: float,
    y_max: Optional[float] = None,
) -> float:
    """Invert the ZINB CDF to obtain a percentile threshold.

    Uses ``scipy.optimize.brentq`` on the integer support ``[0, y_max]``.
    The target percentile is clipped to the achievable range of the ZINB CDF
    (``[pi, 1]``) because values below ``pi`` correspond to the structural-zero
    atom at y=0.
    """
    target = float(np.clip(percentile / 100.0, pi + 1e-12, 1.0 - 1e-12))
    if y_max is None:
        y_max = float(_nb_threshold(mu, alpha, 99.9))
    y_max = max(1.0, y_max)

    lo, hi = 0.0, y_max
    if _zinb_cdf(np.array([lo]), pi, mu, alpha)[0] >= target:
        return 0.0
    if _zinb_cdf(np.array([hi]), pi, mu, alpha)[0] < target:
        return float(hi)

    try:
        root = optimize.brentq(
            lambda yy: float(_zinb_cdf(np.array([yy]), pi, mu, alpha)[0] - target),
            lo,
            hi,
            xtol=1e-6,
        )
    except Exception as exc:
        log.debug(f"ZINB threshold inversion failed: {exc}; falling back to NB threshold")
        return _nb_threshold(mu, alpha, percentile)
    return float(root)


def _check_statsmodels() -> bool:
    """Return True if statsmodels is available for ZINB fitting."""
    try:
        import statsmodels.api as sm  # noqa: F401
        from statsmodels.discrete.count_model import ZeroInflatedNegativeBinomialP  # noqa: F401

        return True
    except Exception:
        return False


def _nb_loglik_logparams(log_params: np.ndarray, y: np.ndarray) -> float:
    """Negative log-likelihood for negative binomial in log parameter space.

    ``log_params = [log(mu), log(alpha)]`` keeps the optimizer in an
    unconstrained domain and improves numerical stability when ``alpha`` is
    close to zero.
    """
    mu = float(np.exp(log_params[0]))
    alpha = float(np.exp(log_params[1]))
    if mu <= 0 or alpha <= 0:
        return 1e12
    p = 1.0 / (1.0 + alpha * mu)
    r = 1.0 / alpha
    log_nb = (
        gammaln(y + r)
        - gammaln(r)
        - gammaln(y + 1)
        + r * np.log(p)
        + y * np.log(1.0 - p)
    )
    return -float(np.sum(log_nb))


def _nb_loglik(params: np.ndarray, y: np.ndarray) -> float:
    """Negative log-likelihood for negative binomial.

    Kept for backward compatibility; new optimization uses
    :func:`_nb_loglik_logparams`.
    """
    mu, alpha = params
    if mu <= 0 or alpha <= 0:
        return 1e12
    return _nb_loglik_logparams(np.log([mu, alpha]), y)


def _zinb_loglik(params: np.ndarray, y: np.ndarray) -> float:
    """Negative log-likelihood for zero-inflated negative binomial."""
    pi, mu, alpha = params
    pi = float(np.clip(pi, 1e-6, 1 - 1e-6))
    if mu <= 0 or alpha <= 0:
        return 1e12
    p = 1.0 / (1.0 + alpha * mu)
    r = 1.0 / alpha
    # NB pmf in log form
    log_nb = (
        gammaln(y + r)
        - gammaln(r)
        - gammaln(y + 1)
        + r * np.log(p)
        + y * np.log(1.0 - p)
    )
    loglik_zero = np.log(pi + (1.0 - pi) * np.exp(log_nb[y == 0]))
    loglik_pos = np.log(1.0 - pi) + log_nb[y > 0]
    ll = np.empty_like(y, dtype=float)
    ll[y == 0] = loglik_zero
    ll[y > 0] = loglik_pos
    return -float(np.sum(ll))


def _zinb_loglik_logparams(log_params: np.ndarray, y: np.ndarray) -> float:
    """ZINB negative log-likelihood with log-transformed ``mu`` and ``alpha``."""
    pi = float(1.0 / (1.0 + np.exp(-log_params[0])))
    pi = float(np.clip(pi, 1e-12, 1.0 - 1e-12))
    mu = float(np.exp(log_params[1]))
    alpha = float(np.exp(log_params[2]))
    return _zinb_loglik(np.array([pi, mu, alpha]), y)


def _fit_count_distribution(
    values: np.ndarray,
    *,
    model: str = "auto",
    random_state: int = 42,
) -> Dict[str, Any]:
    """Fit a count distribution to discrete QC metrics such as n_genes.

    Parameters
    ----------
    values
        1-D array of non-negative counts (e.g. n_genes per cell).
    model
        One of ``"auto"``, ``"zinb"``, ``"nb"`` (Poisson-Gamma / negative
        binomial), ``"poisson"``.  ``"auto"`` tries ZINB, NB and Poisson and
        picks the best by AIC.
    random_state
        Random seed used for reproducibility of any stochastic fitting steps.

    Returns
    -------
    dict with keys ``model``, ``params``, ``aic``, ``bic``, ``threshold_percentile_func``,
    ``is_success``.
    """
    y = np.asarray(values, dtype=float)
    y = y[np.isfinite(y) & (y >= 0)]
    if len(y) == 0:
        return {"model": "none", "params": None, "is_success": False}

    rng = np.random.default_rng(random_state)
    y_mean = float(np.mean(y))
    y_var = float(np.var(y))
    y_alpha = max(1e-6, (y_var - y_mean) / max(y_mean**2, 1e-12))

    candidates: List[Tuple[str, Any]] = []

    def _nb_threshold(mu: float, alpha: float, percentile: float) -> float:
        p = 1.0 / (1.0 + alpha * mu)
        r = 1.0 / alpha
        return float(stats.nbinom.ppf(percentile / 100.0, r, p))

    def _zinb_threshold_wrapper(pi: float, mu: float, alpha: float, percentile: float) -> float:
        return float(_zinb_threshold(percentile, pi, mu, alpha, y_max=None))

    def _poisson_threshold(mu: float, percentile: float) -> float:
        return float(stats.poisson.ppf(percentile / 100.0, mu))

    # Pre-screen: equidispersion / underdispersion favors Poisson.
    # NB becomes degenerate when alpha -> 0; fitting it wastes optimizer time
    # and can produce pathological thresholds.
    y_mean = float(np.mean(y))
    y_var = float(np.var(y))
    dispersion_ratio = y_var / max(y_mean, 1e-12)
    force_poisson = dispersion_ratio <= 1.05

    # 1. Negative Binomial (Poisson-Gamma)
    if not force_poisson:
        try:
            res_nb = _minimize_nb(y, y_mean, y_alpha, rng)
            if res_nb is not None:
                mu_nb, alpha_nb = res_nb
                ll_nb = -_nb_loglik(res_nb, y)
                k_nb = 2  # mu, alpha
                aic_nb = 2 * k_nb - 2 * ll_nb
                candidates.append(
                    (
                        "nb",
                        {
                            "params": {"mu": mu_nb, "alpha": alpha_nb},
                            "aic": aic_nb,
                            "threshold_func": lambda pct, mu=mu_nb, alpha=alpha_nb: _nb_threshold(
                                mu, alpha, pct
                            ),
                        },
                    )
                )
        except Exception as exc:
            log.debug(f"NB fitting failed: {exc}")

    # 2. Zero-inflated NB
    if model in ("auto", "zinb") and _check_statsmodels():
        try:
            import statsmodels.api as sm
            from statsmodels.discrete.count_model import ZeroInflatedNegativeBinomialP

            # statsmodels needs an exog; use intercept only
            exog = np.ones((len(y), 1))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                zinb = ZeroInflatedNegativeBinomialP(
                    y, exog, exog_infl=exog, missing="drop"
                )
                zinb_res = zinb.fit(disp=0, maxiter=100)
            pi_zinb = float(1.0 / (1.0 + np.exp(-zinb_res.params["inflate_const"])))
            mu_zinb = float(np.exp(zinb_res.params["const"]))
            alpha_zinb = float(zinb_res.params["alpha"])
            ll_zinb = float(zinb_res.llf)
            k_zinb = 3
            aic_zinb = 2 * k_zinb - 2 * ll_zinb
            candidates.append(
                (
                    "zinb",
                    {
                        "params": {"pi": pi_zinb, "mu": mu_zinb, "alpha": alpha_zinb},
                        "aic": aic_zinb,
                        "threshold_func": lambda pct, pi=pi_zinb, mu=mu_zinb, alpha=alpha_zinb: _zinb_threshold_wrapper(
                            pi, mu, alpha, pct
                        ),
                    },
                )
            )
        except Exception as exc:
            log.debug(f"statsmodels ZINB fitting failed: {exc}")
    elif model in ("auto", "zinb") and not force_poisson:
        try:
            res_zinb = _minimize_zinb(y, y_mean, y_alpha, rng)
            if res_zinb is not None:
                pi_z, mu_z, alpha_z = res_zinb
                ll_zinb = -_zinb_loglik(res_zinb, y)
                k_zinb = 3
                aic_zinb = 2 * k_zinb - 2 * ll_zinb
                candidates.append(
                    (
                        "zinb",
                        {
                            "params": {"pi": pi_z, "mu": mu_z, "alpha": alpha_z},
                            "aic": aic_zinb,
                            "threshold_func": lambda pct, pi=pi_z, mu=mu_z, alpha=alpha_z: _zinb_threshold_wrapper(
                                pi, mu, alpha, pct
                            ),
                        },
                    )
                )
        except Exception as exc:
            log.debug(f"Fallback ZINB fitting failed: {exc}")

    # 3. Poisson
    if model in ("auto", "poisson"):
        mu_pois = float(np.mean(y))
        ll_pois = float(np.sum(stats.poisson.logpmf(y, mu_pois)))
        k_pois = 1
        aic_pois = 2 * k_pois - 2 * ll_pois
        candidates.append(
            (
                "poisson",
                {
                    "params": {"mu": mu_pois},
                    "aic": aic_pois,
                    "threshold_func": lambda pct, mu=mu_pois: _poisson_threshold(mu, pct),
                },
            )
        )

    # Equidispersion tie-break: if Poisson AIC is within 10 of the best,
    # prefer the simpler model to avoid degenerate NB thresholds.
    if candidates:
        best_aic = min(info["aic"] for _, info in candidates)
        candidates = [
            (name, info)
            for name, info in candidates
            if not (name == "poisson" and info["aic"] > best_aic + 10)
        ]

    if not candidates:
        return {"model": "none", "params": None, "is_success": False}

    best = min(candidates, key=lambda x: x[1]["aic"])
    return {
        "model": best[0],
        "params": best[1]["params"],
        "aic": best[1]["aic"],
        "threshold_func": best[1]["threshold_func"],
        "is_success": True,
        "all_aic": {name: info["aic"] for name, info in candidates},
    }


def _minimize_nb(
    y: np.ndarray, y_mean: float, y_alpha: float, rng: np.random.Generator
) -> Optional[np.ndarray]:
    """Fit NB via differential evolution in log parameter space.

    Optimizes ``[log(mu), log(alpha)]`` to avoid boundary issues when
    ``alpha`` is close to zero.  Uses the data-driven mean and a small set
    of dispersed initial points.
    """
    best_res = None
    best_ll = np.inf

    # Bound log(mu) within roughly [mean/10, mean*10] and log(alpha) within
    # [-12, 2] which covers the vast majority of scRNA-seq count regimes.
    log_mu_lo = np.log(max(1.0, y_mean * 0.1))
    log_mu_hi = np.log(max(10.0, y_mean * 5.0))
    bounds = [(log_mu_lo, log_mu_hi), (-12.0, 2.0)]

    inits = []
    for init_mu in [y_mean, y_mean * 0.5, y_mean * 1.5]:
        for init_alpha in [max(1e-6, y_alpha), 0.1, 1.0]:
            inits.append(np.log([max(init_mu, 1e-6), max(init_alpha, 1e-12)]))

    for init in inits:
        try:
            res = optimize.minimize(
                lambda p: _nb_loglik_logparams(p, y),
                x0=init,
                bounds=bounds,
                method="L-BFGS-B",
                options={"maxiter": 1000},
            )
            if res.success and res.fun < best_ll:
                best_ll = res.fun
                best_res = np.exp(res.x)
        except Exception:
            continue
    return best_res


def _minimize_zinb(
    y: np.ndarray, y_mean: float, y_alpha: float, rng: np.random.Generator
) -> Optional[np.ndarray]:
    """Fit ZINB via differential evolution from several random starts.

    Optimizes ``[logit(pi), log(mu), log(alpha)]`` for numerical stability.
    """
    best_res = None
    best_ll = np.inf

    log_mu_lo = np.log(max(1.0, y_mean * 0.1))
    log_mu_hi = np.log(max(10.0, y_mean * 5.0))
    bounds = [(-6.0, 0.0), (log_mu_lo, log_mu_hi), (-12.0, 2.0)]

    inits = []
    for init_pi in [0.05, 0.15, 0.3]:
        for init_mu in [y_mean, y_mean * 0.7, y_mean * 1.3]:
            for init_alpha in [max(1e-6, y_alpha), 0.1, 1.0]:
                inits.append(
                    np.array([
                        np.log(init_pi / (1.0 - init_pi)),
                        np.log(max(init_mu, 1e-6)),
                        np.log(max(init_alpha, 1e-12)),
                    ])
                )

    for init in inits:
        try:
            res = optimize.minimize(
                lambda p: _zinb_loglik_logparams(p, y),
                x0=init,
                bounds=bounds,
                method="L-BFGS-B",
                options={"maxiter": 1000},
            )
            if res.success and res.fun < best_ll:
                best_ll = res.fun
                logit_pi, log_mu, log_alpha = res.x
                pi = 1.0 / (1.0 + np.exp(-logit_pi))
                best_res = np.array([pi, np.exp(log_mu), np.exp(log_alpha)])
        except Exception:
            continue
    return best_res


def fit_count_mixture_threshold_model(
    values: np.ndarray,
    *,
    direction: str = "lower",
    percentile: float = 10.0,
    model: str = "auto",
    random_state: int = 42,
    fallback: bool = True,
) -> Dict[str, Any]:
    """Fit a count distribution and derive a direction-aware threshold.

    This is preferred over ``fit_gmm_threshold_model`` for discrete QC metrics
    such as ``n_genes`` because count distributions (ZINB/NB/Poisson) respect
    the non-negative integer support of the data.

    Parameters
    ----------
    values
        1-D array of counts.
    direction
        ``"lower"`` for lower-bound thresholds (e.g. min_genes) or
        ``"upper"`` for upper-bound thresholds.
    percentile
        Percentile of the fitted distribution used for the threshold.  For
        ``direction="lower"`` a low percentile (5-10) gives a conservative
        min_genes; for ``direction="upper"`` a high percentile (90-95) gives a
        max_genes/counts threshold.
    model
        Count model to fit; see :func:`_fit_count_distribution`.
    random_state
        Seed for reproducible optimization.
    fallback
        If True and count fitting fails, fall back to
        :func:`fit_gmm_threshold_model`.

    Returns
    -------
    dict with keys ``threshold``, ``model``, ``params``, ``aic``, ``method``,
    ``fallback_used``, ``is_success``.
    """
    if direction == "upper":
        percentile = 100.0 - percentile

    fit = _fit_count_distribution(values, model=model, random_state=random_state)
    if fit.get("is_success"):
        threshold = float(fit["threshold_func"](percentile))
        threshold = max(0.0, threshold)
        return {
            "threshold": threshold,
            "model": fit["model"],
            "params": fit["params"],
            "aic": fit.get("aic"),
            "all_aic": fit.get("all_aic"),
            "method": f"count_mixture ({fit['model']})",
            "fallback_used": False,
            "is_success": True,
        }

    if fallback:
        log.debug("Count mixture fitting failed; falling back to GMM threshold model.")
        values_arr = np.asarray(values)
        if values_arr.size == 0:
            return {
                "threshold": np.nan,
                "model": "none",
                "params": None,
                "aic": None,
                "method": "count mixture (empty input fallback)",
                "fallback_used": True,
                "is_success": True,
            }
        gmm_fit = fit_gmm_threshold_model(
            values_arr, direction=direction, random_state=random_state
        )
        return {
            "threshold": gmm_fit["threshold"],
            "model": "gmm",
            "params": None,
            "aic": None,
            "method": "gmm (count mixture fallback)",
            "fallback_used": True,
            "is_success": True,
            "gmm_fit": gmm_fit,
        }

    return {
        "threshold": np.nan,
        "model": "none",
        "params": None,
        "aic": None,
        "method": "count mixture (failed)",
        "fallback_used": False,
        "is_success": False,
    }


def fit_bimodal_gmm_threshold_model(
    values: np.ndarray,
    *,
    direction: str = "upper",
    random_state: int = 42,
    min_separation: float = 2.0,
) -> Dict[str, Any]:
    """Fit a 2-component GMM and return the Bayes-optimal crossing point.

    Useful for metrics with a clear bimodal structure, e.g. mitochondrial
    percentage in tumor samples (low-MT viable cells vs. high-MT damaged or
    metabolically stressed cells).

    Parameters
    ----------
    values
        1-D metric values.
    direction
        ``"upper"`` returns the crossing point as an upper-bound threshold;
        ``"lower"`` returns it as a lower-bound threshold.
    random_state
        Seed for the GMM fit.
    min_separation
        Minimum separation (in units of the smaller component std) required to
        accept the bimodal model.  If separation is smaller, the function falls
        back to a single Gaussian at the mean.

    Returns
    -------
    dict with ``threshold``, ``n_components``, ``separation``, ``method``,
    ``is_bimodal``, ``model``.
    """
    X = np.asarray(values).reshape(-1, 1)
    X = X[np.isfinite(X).ravel()]
    if len(X) < 10:
        return {
            "threshold": np.nan,
            "n_components": 0,
            "separation": 0.0,
            "method": "bimodal_gmm",
            "is_bimodal": False,
            "model": None,
        }

    gmm = GaussianMixture(
        n_components=2,
        random_state=random_state,
        covariance_type="full",
        reg_covar=1e-6,
        max_iter=200,
        n_init=5,
    )
    try:
        gmm.fit(X)
    except Exception as exc:
        log.debug(f"Bimodal GMM fitting failed: {exc}")
        return {
            "threshold": np.nan,
            "n_components": 0,
            "separation": 0.0,
            "method": "bimodal_gmm",
            "is_bimodal": False,
            "model": None,
        }

    means = gmm.means_.ravel()
    covs = gmm.covariances_.ravel()
    stds = np.sqrt(np.maximum(covs, 1e-12))
    sorted_idx = np.argsort(means)
    low_idx, high_idx = sorted_idx[0], sorted_idx[1]
    separation = float((means[high_idx] - means[low_idx]) / max(stds[low_idx], 1e-12))

    x_grid = np.linspace(float(X.min()), float(X.max()), 2000).reshape(-1, 1)
    proba = gmm.predict_proba(x_grid)
    diff = proba[:, low_idx] - proba[:, high_idx]
    sign_changes = np.where(np.diff(np.sign(diff)))[0]
    if len(sign_changes) > 0:
        sc = sign_changes[0]
        x0, x1 = x_grid[sc, 0], x_grid[sc + 1, 0]
        d0, d1 = diff[sc], diff[sc + 1]
        crossing = float(x0 - d0 * (x1 - x0) / (d1 - d0)) if d1 != d0 else float(x0)
    else:
        crossing = float(np.mean(means))

    is_bimodal = separation >= min_separation
    if not is_bimodal:
        # Fallback: single Gaussian at upper tail
        threshold = float(means[high_idx] + 2.0 * stds[high_idx])
    else:
        threshold = crossing

    return {
        "threshold": threshold,
        "n_components": 2,
        "separation": separation,
        "method": "bimodal_gmm",
        "is_bimodal": is_bimodal,
        "model": gmm,
        "crossing": crossing,
        "component_means": means.tolist(),
        "component_stds": stds.tolist(),
    }


def compute_kdistance_eps(
    values: np.ndarray,
    k: int = 5,
    random_state: Optional[Union[int, np.random.Generator]] = None,
    max_samples: int = KDISTANCE_MAX_SAMPLES,
) -> float:
    """Estimate DBSCAN ``eps`` from the k-distance graph elbow.

    Parameters
    ----------
    values
        1-D metric values.
    k
        Number of nearest neighbours to consider.
    random_state
        Random seed or generator used when subsampling large inputs.
    max_samples
        Maximum number of observations to use for the k-NN computation.  For
        datasets larger than this, a stratified random subsample is drawn to
        keep the operation fast and memory-efficient.

    Returns
    -------
    float
        Estimated eps.  Falls back to half the std if k-NN fails.
    """
    clean = np.asarray(values)[np.isfinite(values)].reshape(-1, 1)
    if len(clean) <= k + 1:
        return float(np.std(clean) * 0.5) if len(clean) > 1 else 1.0

    rng = np.random.default_rng(random_state)

    # Subsample very large datasets to avoid O(N^2) k-NN.
    if len(clean) > max_samples:
        idx = rng.choice(len(clean), size=max_samples, replace=False)
        clean_sub = clean[idx]
    else:
        clean_sub = clean

    try:
        neigh = NearestNeighbors(n_neighbors=min(k + 1, len(clean_sub)))
        neigh.fit(clean_sub)
        distances, _ = neigh.kneighbors(clean_sub)
        k_dist = np.sort(distances[:, k])
        # Use a robust low percentile to avoid including the same cluster
        eps = float(np.percentile(k_dist, 4.0))
        # Guard against pathologically small eps on dense or subsampled data.
        std_eps = float(np.std(clean)) * 0.5
        if eps < 1e-12 * std_eps or not np.isfinite(eps):
            eps = max(std_eps * 0.05, 1e-12)
        return eps
    except Exception as exc:
        log.debug(f"k-distance eps estimation failed: {exc}")
        return float(np.std(clean) * 0.5)


def fit_gmm_threshold_model(
    values: np.ndarray,
    *,
    direction: str,
    random_state: int = 42,
    n_components: Optional[int] = None,
    max_components: Optional[int] = None,
    covariance_type: str = "full",
    reg_covar: float = 1e-6,
    max_iter: int = 200,
    n_init: int = 3,
) -> Dict[str, Any]:
    """Fit a 1D Gaussian mixture and derive a direction-aware threshold.

    Returns a dictionary containing the selected ``GaussianMixture`` model,
    threshold, selected component count, and BIC diagnostics.
    """
    X = values.reshape(-1, 1)
    active_n_components = n_components
    best_bic = np.inf
    null_bic = None

    if active_n_components is None:
        candidate_max = max_components or 2
        candidate_max = max(1, candidate_max)
        best_n = 1
        for k in range(1, candidate_max + 1):
            gmm_k = GaussianMixture(
                n_components=k,
                random_state=random_state,
                covariance_type=covariance_type,
                reg_covar=reg_covar,
                max_iter=max_iter,
                n_init=n_init,
            )
            try:
                gmm_k.fit(X)
                bic_k = gmm_k.bic(X)
                if bic_k < best_bic:
                    best_bic = bic_k
                    best_n = k
            except Exception:
                continue
        active_n_components = best_n

    gmm = GaussianMixture(
        n_components=active_n_components,
        random_state=random_state,
        covariance_type=covariance_type,
        reg_covar=reg_covar,
        max_iter=max_iter,
        n_init=n_init,
    )
    gmm.fit(X)
    bic = gmm.bic(X)

    if active_n_components > 1:
        try:
            null_gmm = GaussianMixture(
                n_components=1,
                random_state=random_state,
                covariance_type=covariance_type,
                reg_covar=reg_covar,
                max_iter=max_iter,
                n_init=n_init,
            )
            null_gmm.fit(X)
            null_bic = null_gmm.bic(X)
        except Exception:
            null_bic = None

    means = gmm.means_.flatten()
    sorted_idx = np.argsort(means)
    x_min, x_max = float(values.min()), float(values.max())
    pad = 0.1 * (x_max - x_min) if x_max > x_min else 1.0
    x_grid = np.linspace(x_min - pad, x_max + pad, 2000).reshape(-1, 1)
    proba = gmm.predict_proba(x_grid)

    crossings = []
    for i in range(active_n_components - 1):
        j = sorted_idx[i]
        k = sorted_idx[i + 1]
        diff = proba[:, j] - proba[:, k]
        sign_changes = np.where(np.diff(np.sign(diff)))[0]
        for sc in sign_changes:
            x0, x1 = x_grid[sc, 0], x_grid[sc + 1, 0]
            d0, d1 = diff[sc], diff[sc + 1]
            if d1 != d0:
                crossing_x = x0 - d0 * (x1 - x0) / (d1 - d0)
                crossings.append((crossing_x, j, k))

    if crossings:
        if direction == "upper":
            target_comp = sorted_idx[-1]
            relevant = [c for c in crossings if target_comp in (c[1], c[2])] or crossings
            threshold = float(max(c[0] for c in relevant))
        else:
            target_comp = sorted_idx[0]
            relevant = [c for c in crossings if target_comp in (c[1], c[2])] or crossings
            threshold = float(min(c[0] for c in relevant))
    else:
        if direction == "upper":
            rel = sorted_idx[-2:] if active_n_components > 1 else sorted_idx[-1:]
        else:
            rel = sorted_idx[:2] if active_n_components > 1 else sorted_idx[:1]
        threshold = float(np.mean(means[rel]))

    return {
        "model": gmm,
        "threshold": threshold,
        "n_components": int(active_n_components),
        "bic": float(bic),
        "null_bic": float(null_bic) if null_bic is not None else None,
    }


def compute_mad_bounds(
    values: np.ndarray,
    nmads: float = 5.0,
    direction: str = "both",
) -> Tuple[float, float]:
    """Compute outlier bounds using Median Absolute Deviation (MAD).

    This is the canonical MAD-based outlier detection implementation used
    throughout scLucid. It replaces duplicated MAD logic in
    ``filtering.py`` and ``adaptive_threshold.py``.

    Parameters
    ----------
    values : np.ndarray
        Input metric values (may contain NaNs).
    nmads : float, default=5.0
        Number of MADs from the median to use as the bound.
    direction : {'upper', 'lower', 'both'}, default='both'
        Which bound(s) to compute.

    Returns:
    -------
    Tuple[float, float]
        ``(lower_bound, upper_bound)``. For ``direction='upper'`` the
        lower bound is ``-inf``; for ``direction='lower'`` the upper
        bound is ``inf``.
    """
    clean = values[~np.isnan(values)]
    if len(clean) == 0:
        return -np.inf, np.inf

    median = float(np.median(clean))
    mad = float(np.median(np.abs(clean - median)))

    if mad == 0:
        log.debug("MAD is zero; bounds collapse to the median.")
        scaled_mad = 0.0
    else:
        scaled_mad = mad * MAD_SCALE_FACTOR

    lower = median - nmads * scaled_mad
    upper = median + nmads * scaled_mad

    if direction == "upper":
        lower = -np.inf
    elif direction == "lower":
        upper = np.inf
    elif direction != "both":
        raise ValueError(f"direction must be 'upper', 'lower', or 'both', got {direction!r}")

    return lower, upper


class AdaptiveThresholdLearner:
    """
    Automatically learn optimal QC thresholds using statistical methods.

    Supports multiple learning strategies:
    - GMM-based: Learn mixture of distributions
    - MAD-based: Median absolute deviation
    - Percentile-based: Statistical percentiles
    - Kernel density: Non-parametric density estimation
    """

    def __init__(
        self,
        method: str = "gmm",
        min_quality_cells: float = 0.5,
        random_state: int = 42,
        dbscan_min_samples_max: int = 50,
        dbscan_subsample_size: int = KDISTANCE_MAX_SAMPLES,
    ):
        """
        Initialize the adaptive threshold learner.

        Args:
            method: Learning method ('gmm', 'mad', 'percentile', 'kde', 'dbscan')
            min_quality_cells: Minimum fraction of cells to retain
            random_state: Random seed for reproducibility
            dbscan_min_samples_max: Upper cap for DBSCAN ``min_samples``.
            dbscan_subsample_size: Maximum observations used for k-distance eps
                estimation on large datasets.
        """
        self.method = method
        self.min_quality_cells = min_quality_cells
        self.random_state = random_state
        self.dbscan_min_samples_max = max(5, int(dbscan_min_samples_max))
        self.dbscan_subsample_size = max(1000, int(dbscan_subsample_size))

        self._learned_thresholds = {}
        self._fitted_models = {}

    def learn_threshold(
        self,
        metric_values: np.ndarray,
        metric_name: str,
        direction: str = "upper",
    ) -> float:
        """
        Learn optimal threshold for a single QC metric.

        Args:
            metric_values: Array of metric values
            metric_name: Name of the metric
            direction: 'upper' (filter high values) or 'lower' (filter low values)

        Returns:
            Learned threshold value
        """
        # Remove NaN and infinite values
        clean_values = metric_values[~np.isnan(metric_values)]
        clean_values = clean_values[~np.isinf(clean_values)]

        if len(clean_values) == 0:
            log.warning(f"No valid values for metric {metric_name}")
            return np.nan

        if self.method == "gmm":
            threshold = self._learn_threshold_gmm(clean_values, direction)
        elif self.method == "mad":
            threshold = self._learn_threshold_mad(clean_values, direction)
        elif self.method == "percentile":
            threshold = self._learn_threshold_percentile(clean_values, direction)
        elif self.method == "kde":
            threshold = self._learn_threshold_kde(clean_values, direction)
        elif self.method == "dbscan":
            threshold = self._learn_threshold_dbscan(clean_values, direction)
        else:
            raise ValueError(f"Unknown method: {self.method}")

        # Apply minimum quality constraint
        threshold = self._apply_min_quality_constraint(clean_values, threshold, direction)

        self._learned_thresholds[metric_name] = threshold

        log.info(
            f"Learned threshold for {metric_name} ({direction}): {threshold:.4f} "
            f"using {self.method}"
        )

        return threshold

    def _learn_threshold_gmm(
        self,
        values: np.ndarray,
        direction: str,
        n_components: int = 2,
    ) -> float:
        """
        Learn threshold using Gaussian Mixture Model.

        Uses the Bayes-optimal decision boundary (posterior-probability crossing
        point) between the relevant components instead of an unweighted mean of
        means.  This is stable even when the component variances differ.
        """
        try:
            fit = fit_gmm_threshold_model(
                values,
                direction=direction,
                random_state=self.random_state,
                n_components=n_components,
                covariance_type="full",
                reg_covar=1e-6,
                max_iter=200,
                n_init=3,
            )
        except Exception as e:
            log.debug(f"GMM fitting failed: {e}, falling back to percentile")
            return self._learn_threshold_percentile(values, direction)
        self._fitted_models["gmm"] = fit["model"]
        return float(fit["threshold"])

    def _learn_threshold_mad(
        self,
        values: np.ndarray,
        direction: str,
        nmads: float = 5.0,
    ) -> float:
        """
        Learn threshold using Median Absolute Deviation.

        Delegates to the canonical ``compute_mad_bounds`` so that the
        same MAD logic is used everywhere in the QC module.
        """
        lower, upper = compute_mad_bounds(values, nmads=nmads, direction=direction)

        if direction == "upper":
            return upper
        else:
            return max(0.0, lower)

    def _learn_threshold_percentile(
        self,
        values: np.ndarray,
        direction: str,
    ) -> float:
        """
        Learn threshold using percentiles.

        Conservative approach based on distribution statistics.
        """
        if direction == "upper":
            # Use 95th percentile for upper threshold
            threshold = np.percentile(values, 95)
        else:
            # Use 5th percentile for lower threshold
            threshold = np.percentile(values, 5)

        return float(threshold)

    def _learn_threshold_kde(
        self,
        values: np.ndarray,
        direction: str,
    ) -> float:
        """
        Learn threshold using Kernel Density Estimation.

        Finds local minima in density as threshold boundaries.
        """
        try:
            from scipy.stats import gaussian_kde

            # Fit KDE
            kde = gaussian_kde(values)
            x_range = np.linspace(values.min(), values.max(), 1000)
            density = kde(x_range)

            # Find local minima
            from scipy.signal import find_peaks

            # Invert density to find minima
            minima_indices, _ = find_peaks(-density, distance=20)

            if len(minima_indices) == 0:
                # No clear minima, fall back to percentile
                return self._learn_threshold_percentile(values, direction)

            if direction == "upper":
                # Use rightmost local minimum
                threshold_idx = minima_indices[-1]
            else:
                # Use leftmost local minimum
                threshold_idx = minima_indices[0]

            threshold = x_range[threshold_idx]

            return float(threshold)

        except Exception as e:
            log.debug(f"KDE failed: {e}, falling back to percentile")
            return self._learn_threshold_percentile(values, direction)

    def _learn_threshold_dbscan(
        self,
        values: np.ndarray,
        direction: str,
    ) -> float:
        """
        Learn threshold using DBSCAN clustering.

        Identifies outliers as low-quality cells.  ``eps`` is estimated from
        the k-distance graph and ``min_samples`` is capped so the algorithm
        remains stable on large datasets.
        """
        try:
            clean = values[np.isfinite(values)]
            if len(clean) == 0:
                return self._learn_threshold_percentile(values, direction)

            X = clean.reshape(-1, 1)
            eps = compute_kdistance_eps(
                clean,
                k=5,
                random_state=self.random_state,
                max_samples=self.dbscan_subsample_size,
            )
            # Use sub-linear scaling for min_samples so DBSCAN does not require
            # an impossibly dense neighbourhood on very large datasets.  The
            # square-root scaling keeps local-density semantics while capping
            # the absolute value.
            min_samples = min(
                self.dbscan_min_samples_max,
                max(5, int(np.sqrt(len(clean)))),
            )

            dbscan = DBSCAN(eps=eps, min_samples=min_samples)
            labels = dbscan.fit_predict(X)

            # If every point is noise the eps estimate was too small; fall back
            # to a percentile-based threshold rather than declaring all cells
            # outliers.
            if np.all(labels == -1):
                return self._learn_threshold_percentile(values, direction)

            outlier_values = clean[labels == -1]

            if len(outlier_values) == 0:
                # No outliers detected, use percentile
                return self._learn_threshold_percentile(values, direction)

            if direction == "upper":
                # Threshold is minimum of upper outliers
                threshold = np.min(outlier_values)
            else:
                # Threshold is maximum of lower outliers
                threshold = np.max(outlier_values)

            return float(threshold)

        except Exception as e:
            log.debug(f"DBSCAN failed: {e}, falling back to percentile")
            return self._learn_threshold_percentile(values, direction)

    def _apply_min_quality_constraint(
        self,
        values: np.ndarray,
        threshold: float,
        direction: str,
    ) -> float:
        """
        Apply minimum quality constraint to threshold.

        Ensures that at least min_quality_cells fraction passes QC.
        """
        n_cells = len(values)
        min_cells_to_keep = int(n_cells * self.min_quality_cells)

        if direction == "upper":
            # At most this fraction can fail
            max_failures = n_cells - min_cells_to_keep

            # Count cells that would fail
            n_failures = np.sum(values > threshold)

            if n_failures > max_failures:
                # Adjust threshold to keep minimum cells
                threshold = np.sort(values)[-max_failures]
        else:
            # At most this fraction can fail
            max_failures = n_cells - min_cells_to_keep

            # Count cells that would fail
            n_failures = np.sum(values < threshold)

            if n_failures > max_failures:
                # Adjust threshold to keep minimum cells
                threshold = np.sort(values)[max_failures]

        return float(threshold)

    def learn_all_thresholds(
        self,
        adata: AnnData,
        metrics: Optional[Dict[str, str]] = None,
    ) -> Dict[str, float]:
        """
        Learn thresholds for all specified QC metrics.

        Args:
            adata: AnnData object with QC metrics
            metrics: Dictionary of {metric_name: direction} pairs
                    If None, uses default metrics

        Returns:
            Dictionary of learned thresholds
        """
        if metrics is None:
            metrics = {
                "log1p_n_genes_by_counts": "lower",
                "log1p_total_counts": "lower",
                "pct_counts_mt": "upper",
                "pct_counts_in_top_20_genes": "upper",
            }

        learned_thresholds = {}

        for metric_name, direction in metrics.items():
            if metric_name not in adata.obs:
                log.warning(f"Metric {metric_name} not found in adata.obs")
                continue

            values = adata.obs[metric_name].values

            try:
                threshold = self.learn_threshold(values, metric_name, direction)
                learned_thresholds[metric_name] = threshold
            except Exception as e:
                log.error(f"Failed to learn threshold for {metric_name}: {e}")
                continue

        return learned_thresholds

    def predict_quality(
        self,
        metric_values: np.ndarray,
        metric_name: str,
        direction: str,
    ) -> np.ndarray:
        """
        Predict quality labels for cells based on learned threshold.

        Args:
            metric_values: Metric values for cells
            metric_name: Name of the metric
            direction: Filter direction

        Returns:
            Boolean array (True = high quality, False = low quality)
        """
        if metric_name not in self._learned_thresholds:
            raise ValueError(f"No learned threshold for {metric_name}")

        threshold = self._learned_thresholds[metric_name]

        if direction == "upper":
            quality = metric_values <= threshold
        else:
            quality = metric_values >= threshold

        return quality


class MultiMetricAdaptiveLearner:
    """
    Adaptive threshold learner that considers multiple metrics jointly.

    Uses multivariate approaches to find optimal threshold combinations.
    """

    def __init__(
        self,
        method: str = "isolation_forest",
        contamination: float = 0.1,
        random_state: int = 42,
    ):
        """
        Initialize multi-metric learner.

        Args:
            method: Method ('isolation_forest', 'local_outlier_factor', 'one_class_svm')
            contamination: Expected fraction of outliers
            random_state: Random seed
        """
        self.method = method
        self.contamination = contamination
        self.random_state = random_state

        self._model = None

    def fit(
        self,
        adata: AnnData,
        metrics: List[str],
    ):
        """
        Fit the multi-metric outlier detection model.

        Args:
            adata: AnnData object with QC metrics
            metrics: List of metric names to use
        """
        # Prepare data matrix
        X = np.column_stack([adata.obs[m].values for m in metrics])

        # Handle missing values
        X = np.nan_to_num(X, nan=0.0)

        if self.method == "isolation_forest":
            from sklearn.ensemble import IsolationForest

            self._model = IsolationForest(
                contamination=self.contamination,
                random_state=self.random_state,
                n_jobs=effective_n_jobs(-1),
            )
            self._model.fit(X)

        elif self.method == "local_outlier_factor":
            from sklearn.neighbors import LocalOutlierFactor

            self._model = LocalOutlierFactor(
                contamination=self.contamination,
                n_neighbors=20,
                n_jobs=effective_n_jobs(-1),
            )
            self._model.fit(X)

        elif self.method == "one_class_svm":
            from sklearn.svm import OneClassSVM

            self._model = OneClassSVM(
                nu=self.contamination,
                kernel="rbf",
            )
            self._model.fit(X)

        else:
            raise ValueError(f"Unknown method: {self.method}")

    def predict(self, adata: AnnData, metrics: List[str]) -> np.ndarray:
        """
        Predict quality labels using fitted model.

        Args:
            adata: AnnData object
            metrics: List of metric names

        Returns:
            Boolean array (True = high quality, False = outlier/low quality)
        """
        if self._model is None:
            raise ValueError("Model not fitted. Call fit() first.")

        # Prepare data
        X = np.column_stack([adata.obs[m].values for m in metrics])
        X = np.nan_to_num(X, nan=0.0)

        # Get predictions
        if self.method == "local_outlier_factor":
            # LOF has fit_predict instead of predict
            predictions = self._model.fit_predict(X)
        else:
            predictions = self._model.predict(X)

        # Convert to boolean (1 = inlier/high quality, -1 = outlier)
        quality = predictions == 1

        return quality

    def fit_predict(
        self,
        adata: AnnData,
        metrics: List[str],
    ) -> np.ndarray:
        """
        Fit model and return predictions.

        Args:
            adata: AnnData object
            metrics: List of metric names

        Returns:
            Boolean array (True = high quality)
        """
        self.fit(adata, metrics)
        return self.predict(adata, metrics)
