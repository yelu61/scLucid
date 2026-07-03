"""Runtime compatibility shims for the unmaintained scrublet package.

scrublet (https://github.com/swolock/scrublet) has not been updated for modern
NumPy, SciPy, and Matplotlib. This module applies minimal, idempotent patches
at import time so that scLucid can keep using scrublet as an optional
algorithm without maintaining a fork.

Shims applied:
1. ``scipy.sparse.spmatrix.A`` is aliased to ``.toarray()`` if it is missing or
   deprecated.
2. ``scrublet.helper_functions.tot_counts_norm`` is replaced with a version that
   does not rely on ``np.matrix``/``np.asarray(...)[0, :]`` indexing.
3. ``scrublet.helper_functions.darken_cmap`` and ``custom_cmap`` are replaced
   with versions that use the modern Matplotlib
   ``LinearSegmentedColormap.from_list`` signature.
4. ``scrublet.Scrublet.call_doublets`` is wrapped so that bare ``except:`` does
   not swallow ``KeyboardInterrupt``/``SystemExit``.

These shims are applied lazily inside ``_run_scrublet`` and only once per
process. They are guarded by ``try/except ImportError`` so the module is safe to
import even when scrublet is not installed.
"""

from __future__ import annotations

import logging
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse
from matplotlib.colors import ListedColormap

log = logging.getLogger(__name__)

# Sentinel to ensure shims are applied only once per process.
_SCRUBLET_SHIMS_APPLIED = False


def _ensure_scrublet_compatibility() -> bool:
    """Return whether Scrublet may need NumPy compatibility handling."""
    return not hasattr(np.ndarray, "ptp")


def _coerce_scrublet_array(values, *, expected_len: int, name: str, sample_name: str):
    """Return a 1-D numpy array when a Scrublet output has the expected length."""
    if values is None:
        return None

    arr = np.asarray(values).ravel()
    if arr.shape[0] != expected_len:
        log.warning(
            "Ignoring Scrublet %s for sample '%s': expected %d values, got %d.",
            name,
            sample_name,
            expected_len,
            arr.shape[0],
        )
        return None
    return arr


def _scrublet_scores_degenerate(scores) -> bool:
    """Return True when Scrublet scores carry no usable ranking signal."""
    if scores is None:
        return False
    arr = np.asarray(scores, dtype=float).ravel()
    finite = arr[np.isfinite(arr)]
    if finite.size < 2:
        return True
    return bool(np.nanmax(finite) - np.nanmin(finite) <= 1e-12)


def _scrub_doublets(scrub, *, n_prin_comps: int, use_approx_neighbors: bool):
    """Call scrub_doublets with exact/approx neighbor control when supported."""
    try:
        return scrub.scrub_doublets(
            n_prin_comps=n_prin_comps,
            verbose=False,
            use_approx_neighbors=use_approx_neighbors,
        )
    except TypeError:
        return scrub.scrub_doublets(n_prin_comps=n_prin_comps, verbose=False)


def _plot_scrublet_embedding_fallback(
    embedding: np.ndarray,
    scores: np.ndarray,
    predicted: np.ndarray,
    *,
    title: str,
):
    """Plot Scrublet embedding without relying on scrublet's NumPy-1.x plotting code."""
    embedding = np.asarray(embedding)
    scores = np.asarray(scores, dtype=float).ravel()
    predicted = np.asarray(predicted, dtype=bool).ravel()

    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    x = embedding[:, 0]
    y = embedding[:, 1]
    x_range = float(np.max(x) - np.min(x)) if x.size else 0.0
    y_range = float(np.max(y) - np.min(y)) if y.size else 0.0
    x_pad = x_range * 0.05
    y_pad = y_range * 0.05
    xlim = (float(np.min(x)) - x_pad, float(np.max(x)) + x_pad)
    ylim = (float(np.min(y)) - y_pad, float(np.max(y)) + y_pad)
    order = np.argsort(scores)

    axes[0].scatter(
        x[order],
        y[order],
        s=5,
        c=predicted[order],
        cmap=ListedColormap(["#BDBDBD", "#000000"]),
        edgecolors="none",
    )
    axes[0].set_title("Predicted doublets")

    scatter = axes[1].scatter(
        x[order],
        y[order],
        s=5,
        c=scores[order],
        cmap="Reds",
        edgecolors="none",
    )
    axes[1].set_title("Doublet score")
    fig.colorbar(scatter, ax=axes[1])

    for ax in axes:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def _compute_scrublet_umap_embedding(manifold_obs: np.ndarray, random_state: int = 61):
    """Compute a UMAP embedding from scrublet's manifold observations.

    Falls back to a direct ``umap.UMAP`` call when scrublet's ``get_umap`` is
    unavailable or incompatible with the installed NumPy version.
    """
    try:
        import scrublet as scr

        embedding = scr.get_umap(manifold_obs, 10, min_dist=0.3)
        if embedding is not None and embedding.shape[0] == manifold_obs.shape[0]:
            return np.asarray(embedding)
    except Exception as exc:
        log.debug("scrublet.get_umap failed (%s); falling back to umap-learn.", exc)

    try:
        import umap

        reducer = umap.UMAP(
            n_neighbors=10,
            min_dist=0.3,
            n_components=2,
            random_state=random_state,
        )
        return reducer.fit_transform(np.asarray(manifold_obs))
    except Exception as exc:
        log.warning("Could not compute UMAP embedding via umap-learn: %s", exc)
        raise


def _patch_sparse_A() -> None:
    """Alias ``scipy.sparse.spmatrix.A`` to ``.toarray()`` if absent.

    Modern SciPy deprecates and, in 2.x, removes the ``.A`` attribute on sparse
    matrices. scrublet uses ``.A`` in several places (e.g.
    ``E.sum(1).A.squeeze()``). Restoring the alias avoids breaking scrublet
    without modifying its source.
    """
    if not hasattr(scipy.sparse.spmatrix, "A"):
        try:
            scipy.sparse.spmatrix.A = property(lambda self: self.toarray())
            log.debug("Patched scipy.sparse.spmatrix.A -> .toarray()")
        except Exception as exc:  # pragma: no cover
            log.warning("Could not patch scipy.sparse.spmatrix.A: %s", exc)


def _patch_tot_counts_norm() -> None:
    """Replace ``scrublet.helper_functions.tot_counts_norm`` with a safe version.

    The original implementation relies on ``np.asarray(...)[0, :]`` behaving
    like a 2-D ``np.matrix`` slice. With modern SciPy/NumPy, ``.sum(axis=0)``
    may return a 1-D array, making ``[0, :]`` raise an ``IndexError``.
    """
    try:
        import scrublet.helper_functions as hf
    except ImportError:
        return

    if getattr(hf.tot_counts_norm, "_sclucid_patched", False):
        return

    original = hf.tot_counts_norm

    def tot_counts_norm(
        E: Any,
        total_counts: Any = None,
        exclude_dominant_frac: float = 1.0,
        included: Any = None,
        target_total: Any = None,
    ) -> Any:
        """Cell-level total counts normalization, excluding overly abundant genes."""
        if included is None:
            included = []

        E = E.tocsc()
        ncell = E.shape[0]
        if total_counts is None:
            if len(included) == 0:
                if exclude_dominant_frac == 1:
                    tots_use = E.sum(axis=1)
                else:
                    tots = np.asarray(E.sum(axis=1)).ravel()
                    wtmp = scipy.sparse.lil_matrix((ncell, ncell))
                    with np.errstate(divide="ignore", invalid="ignore"):
                        wtmp.setdiag(np.where(tots > 0, 1.0 / tots, 0.0))
                    frac = (wtmp * E).tocsr()
                    included = np.asarray(~((frac > exclude_dominant_frac).sum(axis=0) > 0)).ravel()
                    tots_use = E[:, included].sum(axis=1)
                    log.debug(
                        "scrublet normalization excluded %i dominant genes",
                        int(np.sum(~included)),
                    )
            else:
                tots_use = E[:, included].sum(axis=1)
        else:
            tots_use = total_counts.copy()

        tots_use = np.asarray(tots_use).ravel()

        if target_total is None:
            target_total = float(np.mean(tots_use))

        with np.errstate(divide="ignore", invalid="ignore"):
            scale = np.where(tots_use > 0, float(target_total) / tots_use, 0.0)

        w = scipy.sparse.lil_matrix((ncell, ncell))
        w.setdiag(scale)
        Enorm = w * E

        return Enorm.tocsc()

    tot_counts_norm._sclucid_patched = True
    tot_counts_norm._sclucid_original = original
    hf.tot_counts_norm = tot_counts_norm
    log.debug("Patched scrublet.helper_functions.tot_counts_norm")


def _patch_cmaps() -> None:
    """Patch scrublet colormap helpers for modern Matplotlib.

    ``LinearSegmentedColormap.from_list`` now requires a string ``name`` as the
    first argument and accepts ``N`` as a keyword. The old scrublet calls pass
    an integer name, which fails or warns on recent Matplotlib.
    """
    try:
        import matplotlib.colors as mcolors
        import scrublet.helper_functions as hf
    except ImportError:
        return

    if getattr(hf.darken_cmap, "_sclucid_patched", False):
        return

    original_darken = hf.darken_cmap
    original_custom = hf.custom_cmap

    def darken_cmap(cmap: Any, scale_factor: float) -> Any:
        cdat = np.zeros((cmap.N, 4))
        for ii in range(cdat.shape[0]):
            curcol = cmap(ii)
            cdat[ii, 0] = curcol[0] * scale_factor
            cdat[ii, 1] = curcol[1] * scale_factor
            cdat[ii, 2] = curcol[2] * scale_factor
            cdat[ii, 3] = 1
        return mcolors.LinearSegmentedColormap.from_list(f"{cmap.name}_darkened", cdat, N=cmap.N)

    def custom_cmap(rgb_list: Any) -> Any:
        rgb_list = np.array(rgb_list)
        return mcolors.LinearSegmentedColormap.from_list(
            "custom_scrublet_cmap", rgb_list, N=rgb_list.shape[0]
        )

    darken_cmap._sclucid_patched = True
    darken_cmap._sclucid_original = original_darken
    custom_cmap._sclucid_patched = True
    custom_cmap._sclucid_original = original_custom

    hf.darken_cmap = darken_cmap
    hf.custom_cmap = custom_cmap
    log.debug("Patched scrublet.helper_functions.darken_cmap and custom_cmap")


def _patch_call_doublets() -> None:
    """Wrap ``Scrublet.call_doublets`` so bare ``except:`` does not catch interrupts.

    The original method uses a bare ``except:`` when automatic threshold
    detection fails, which also catches ``KeyboardInterrupt`` and
    ``SystemExit``. We wrap it to reraise base exceptions and only convert
    ordinary exceptions into a ``None`` prediction (which scLucid already
    handles with fallback thresholding).
    """
    try:
        from scrublet import Scrublet
    except ImportError:
        return

    if getattr(Scrublet.call_doublets, "_sclucid_patched", False):
        return

    original_call_doublets = Scrublet.call_doublets

    def call_doublets(self: Any, threshold: Any = None, verbose: bool = True) -> Any:
        try:
            if threshold is None:
                return original_call_doublets(self, verbose=verbose)
            return original_call_doublets(self, threshold=threshold, verbose=verbose)
        except TypeError as exc:
            try:
                return original_call_doublets(self, verbose=verbose)
            except (KeyboardInterrupt, SystemExit):
                raise
            except Exception:
                log.debug(
                    "scrublet.call_doublets signature fallback failed after %s; returning None",
                    type(exc).__name__,
                )
                self.predicted_doublets_ = None
                return None
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as exc:
            log.debug("scrublet.call_doublets raised %s; returning None", type(exc).__name__)
            self.predicted_doublets_ = None
            return None

    call_doublets._sclucid_patched = True
    call_doublets._sclucid_original = original_call_doublets
    Scrublet.call_doublets = call_doublets
    log.debug("Patched scrublet.Scrublet.call_doublets")


def apply_scrublet_compatibility_shims() -> None:
    """Apply all scrublet compatibility shims idempotently."""
    global _SCRUBLET_SHIMS_APPLIED
    if _SCRUBLET_SHIMS_APPLIED:
        return

    _patch_sparse_A()
    _patch_tot_counts_norm()
    _patch_cmaps()
    _patch_call_doublets()

    _SCRUBLET_SHIMS_APPLIED = True
    log.debug("Scrublet compatibility shims applied")
