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

import numpy as np
import scipy.sparse

log = logging.getLogger(__name__)

# Sentinel to ensure shims are applied only once per process.
_SCRUBLET_SHIMS_APPLIED = False


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
                    print("Excluded %i genes from normalization" % (np.sum(~included)))
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
        return mcolors.LinearSegmentedColormap.from_list(
            f"{cmap.name}_darkened", cdat, N=cmap.N
        )

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
            return original_call_doublets(self, threshold=threshold, verbose=verbose)
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
