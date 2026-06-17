"""External doublet detection algorithm wrappers.

Extracted from core.py for maintainability.
"""

from __future__ import annotations

import gc
import logging
from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
import scipy.sparse as sparse
from anndata import AnnData

from ..config import DoubletConfig
from ._scrublet_compat import apply_scrublet_compatibility_shims

log = logging.getLogger(__name__)


def _raw_count_guard(adata: AnnData, *, sample_name: str, method: str) -> bool:
    """Return True when the active matrix looks compatible with raw UMI counts."""
    X = adata.X
    if X is None:
        log.warning("Skipping %s for sample '%s': expression matrix is missing.", method, sample_name)
        return False

    values = X.data if sparse.issparse(X) else np.asarray(X).ravel()
    if values.size == 0:
        log.warning("Skipping %s for sample '%s': expression matrix is empty.", method, sample_name)
        return False

    finite = np.asarray(values[np.isfinite(values)], dtype=float)
    if finite.size == 0:
        log.warning("Skipping %s for sample '%s': matrix has no finite values.", method, sample_name)
        return False
    if np.min(finite) < 0:
        log.warning(
            "Skipping %s for sample '%s': matrix contains negative values and is not raw counts.",
            method,
            sample_name,
        )
        return False

    positive = finite[finite > 0]
    if positive.size:
        fractional_fraction = float(np.mean(np.abs(positive - np.rint(positive)) > 1e-6))
        if fractional_fraction > 0.01:
            log.warning(
                "Skipping %s for sample '%s': %.1f%% of positive values are fractional; "
                "doublet algorithms require raw UMI-like counts.",
                method,
                sample_name,
                fractional_fraction * 100.0,
            )
            return False

    total_entries = adata.n_obs * adata.n_vars
    if total_entries > 0:
        nonzero = X.nnz if sparse.issparse(X) else int(np.count_nonzero(X))
        zero_fraction = 1.0 - (float(nonzero) / float(total_entries))
        if zero_fraction < 0.05:
            log.warning(
                "Skipping %s for sample '%s': matrix has very few zeros (%.1f%%), "
                "suggesting normalized/transformed input rather than raw counts.",
                method,
                sample_name,
                zero_fraction * 100.0,
            )
            return False

    return True


def _ensure_scrublet_compatibility():
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


def _expected_rate_topk_predictions(scores, expected_rate: float) -> Tuple[np.ndarray, float]:
    """Flag the top expected-rate cells by score, robust to tied quantiles."""
    scores = np.asarray(scores, dtype=float).ravel()
    finite_mask = np.isfinite(scores)
    predicted = np.zeros(scores.shape[0], dtype=bool)
    if scores.size == 0 or not np.any(finite_mask):
        return predicted, float("nan")

    rate = max(0.0, min(1.0, float(expected_rate)))
    n_expected = int(np.ceil(rate * scores.size))
    if rate > 0 and n_expected == 0:
        n_expected = 1
    n_expected = min(n_expected, int(np.sum(finite_mask)))
    if n_expected <= 0:
        return predicted, float(np.nanmax(scores[finite_mask]))

    finite_indices = np.flatnonzero(finite_mask)
    finite_scores = scores[finite_mask]
    order = np.argsort(finite_scores, kind="mergesort")
    selected = finite_indices[order[-n_expected:]]
    predicted[selected] = True
    threshold = float(np.min(scores[selected]))
    return predicted, threshold


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


def _run_scrublet(
    adata_view: AnnData,
    sample_name: str,
    config: DoubletConfig,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Run Scrublet algorithm for doublet detection on a single AnnData view.
    Returns (scores, predicted) arrays.
    """
    _patched = _ensure_scrublet_compatibility()

    rate = config.expected_doublet_rate
    current_rate = rate.get(sample_name, 0.1) if isinstance(rate, dict) else rate
    if current_rate is None:  # Handle case where rate is not provided for a sample
        log.warning(f"No doublet rate provided for sample '{sample_name}', using default of 0.1.")
        current_rate = 0.1

    # Data-quality guard: skip samples with too few features or too-low counts
    if not _raw_count_guard(adata_view, sample_name=sample_name, method="Scrublet"):
        return None, None

    if adata_view.n_vars < 100:
        log.warning(
            f"Skipping Scrublet for sample '{sample_name}': only {adata_view.n_vars} genes "
            f"(minimum 100 required for reliable doublet detection)."
        )
        return None, None
    _cell_sums = np.array(adata_view.X.sum(axis=1)).ravel()
    median_counts = float(np.median(_cell_sums))
    if median_counts < 200:
        log.warning(
            f"Skipping Scrublet for sample '{sample_name}': median UMI count {median_counts:.0f} "
            f"is too low (minimum 200 required)."
        )
        return None, None

    actual_n_pcs = min(config.scr_n_pcs, adata_view.n_obs - 1, adata_view.n_vars - 1)

    # Apply compatibility shims for modern NumPy/SciPy/Matplotlib before importing scrublet.
    apply_scrublet_compatibility_shims()

    try:
        import scrublet as scr

        scrub = scr.Scrublet(adata_view.X, expected_doublet_rate=current_rate)
        scores, predicted = _scrub_doublets(
            scrub,
            n_prin_comps=actual_n_pcs,
            use_approx_neighbors=True,
        )
        scores = _coerce_scrublet_array(
            scores,
            expected_len=adata_view.n_obs,
            name="scores",
            sample_name=sample_name,
        )
        predicted = _coerce_scrublet_array(
            predicted,
            expected_len=adata_view.n_obs,
            name="predictions",
            sample_name=sample_name,
        )

        if _scrublet_scores_degenerate(scores):
            log.warning(
                "Scrublet approximate-neighbor scores are degenerate for sample '%s'; "
                "rerunning with exact nearest neighbors.",
                sample_name,
            )
            scrub = scr.Scrublet(adata_view.X, expected_doublet_rate=current_rate)
            scores, predicted = _scrub_doublets(
                scrub,
                n_prin_comps=actual_n_pcs,
                use_approx_neighbors=False,
            )
            scores = _coerce_scrublet_array(
                scores,
                expected_len=adata_view.n_obs,
                name="exact-neighbor scores",
                sample_name=sample_name,
            )
            predicted = _coerce_scrublet_array(
                predicted,
                expected_len=adata_view.n_obs,
                name="exact-neighbor predictions",
                sample_name=sample_name,
            )

        if predicted is None:
            try:
                predicted = scrub.call_doublets(verbose=False)
            except Exception as e:
                log.warning(f"Scrublet call_doublets failed for sample '{sample_name}': {e}")
            predicted = _coerce_scrublet_array(
                predicted,
                expected_len=adata_view.n_obs,
                name="call_doublets predictions",
                sample_name=sample_name,
            )

        if scores is None:
            scores = _coerce_scrublet_array(
                getattr(scrub, "doublet_scores_obs_", None),
                expected_len=adata_view.n_obs,
                name="doublet_scores_obs_",
                sample_name=sample_name,
            )

        if predicted is None:
            predicted = _coerce_scrublet_array(
                getattr(scrub, "predicted_doublets_", None),
                expected_len=adata_view.n_obs,
                name="predicted_doublets_",
                sample_name=sample_name,
            )

        if predicted is None and scores is not None:
            predicted, threshold = _expected_rate_topk_predictions(scores, float(current_rate))
            log.warning(
                "Scrublet did not return binary predictions for sample '%s'; "
                "falling back to expected-rate top-score threshold %.4f.",
                sample_name,
                threshold,
            )

        if predicted is None:
            log.warning("Scrublet produced no usable scores or predictions for sample '%s'.", sample_name)
            return None, None
        predicted = np.asarray(predicted, dtype=bool).ravel()

        doublet_count = sum(predicted)
        doublet_rate = doublet_count / len(predicted)
        log.info(f"  Found {doublet_count} potential doublets via Scrublet ({doublet_rate:.2%})")

        if config.scr_plot_umap:
            try:
                embedding = _compute_scrublet_umap_embedding(
                    np.asarray(scrub.manifold_obs_),
                    random_state=config.random_state,
                )
                scrub.set_embedding("UMAP", embedding)
                try:
                    before_figs = set(plt.get_fignums())
                    fig, ax = scrub.plot_embedding("UMAP", order_points=True)
                except Exception as plot_exc:
                    for fig_num in set(plt.get_fignums()) - before_figs:
                        plt.close(fig_num)
                    log.warning(
                        "Scrublet native embedding plot failed for sample %s (%s); "
                        "using scLucid fallback plotter.",
                        sample_name,
                        plot_exc,
                    )
                    fig = _plot_scrublet_embedding_fallback(
                        scrub._embeddings["UMAP"],
                        scores,
                        predicted,
                        title=f"{sample_name} Scrublet doublets",
                    )
                if config.save_dir:
                    umap_dir = Path(config.save_dir) / "doublet_umaps"
                    umap_dir.mkdir(parents=True, exist_ok=True)
                    save_path = umap_dir / f"{sample_name}_doublets_umap.png"
                    fig.savefig(save_path, dpi=300, bbox_inches="tight")
                    log.info("Saved Scrublet UMAP for sample '%s' to %s", sample_name, save_path)
                if config.show_plots:
                    plt.show()
                else:
                    plt.close(fig)
            except Exception as e:
                log.warning("Could not generate UMAP for sample %s: %s", sample_name, e)

        return scores, predicted

    except Exception as e:
        log.error(f"Scrublet failed for sample {sample_name}: {e}")
        return None, None
    finally:
        if _patched:
            log.debug("Scrublet NumPy compatibility flag set; no global ndarray patch was applied.")
        gc.collect()


def _run_solo(
    adata_view: AnnData,
    sample_name: str,
    config: DoubletConfig,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Run Solo (from scvi-tools) for doublet detection on a single AnnData view.

    Args:
        adata_view: AnnData object for a single sample.
        sample_name: Name of the sample being processed.
        config: The main DoubletConfig object.

    Returns:
        A tuple containing (scores, predicted) numpy arrays.
    """
    try:
        import scvi
        import torch
    except ImportError:
        log.error(
            "scvi-tools is not installed. Please install it to use the 'solo' method: pip install scvi-tools"
        )
        return None, None

    # --- 1. Data Preparation ---
    adata_solo = adata_view.copy()
    if config.solo_use_raw and adata_solo.raw:
        log.info("  Using 'adata.raw' for Solo as configured.")
        adata_solo = adata_solo.raw.to_adata()
    else:
        log.info("  Using 'adata.X' for Solo.")

    if not _raw_count_guard(adata_solo, sample_name=sample_name, method="Solo"):
        return None, None

    log.info("  Setting up AnnData for scvi-tools model...")
    scvi.model.SCVI.setup_anndata(adata_solo)

    # --- 2. VAE Model Training ---
    log.info("  Training the underlying VAE model...")
    vae_model = scvi.model.SCVI(adata_solo)

    use_gpu_flag = torch.cuda.is_available() and config.solo_use_gpu
    accelerator = "gpu" if use_gpu_flag else "cpu"
    devices = 1 if use_gpu_flag else "auto"

    vae_model.train(
        max_epochs=config.solo_n_epochs,
        accelerator=accelerator,
        devices=devices,
        plan_kwargs={"lr": config.solo_learning_rate},
        early_stopping=True,
        early_stopping_patience=20,
        enable_progress_bar=False,
        logger=False,
    )

    # --- 3. Solo Model Training and Prediction ---
    log.info("  Training the Solo model for doublet detection...")
    solo_model = scvi.external.SOLO.from_scvi_model(vae_model)
    solo_model.train(
        accelerator=accelerator,
        devices=devices,
        enable_progress_bar=False,
        logger=False,
    )

    log.info("  Predicting doublets with Solo...")

    # In newer scvi-tools, .predict() returns a Series of labels
    predictions_series = solo_model.predict(soft=False)
    predicted = (predictions_series == "doublet").values

    # The scores retrieval has changed across scvi-tools versions.
    # Try the newer API first, then fall back to legacy APIs.
    scores = None
    try:
        scores_df = solo_model.get_scores()
        if isinstance(scores_df, pd.DataFrame) and "doublet_scores" in scores_df.columns:
            scores = scores_df["doublet_scores"].values
        elif isinstance(scores_df, pd.Series):
            scores = scores_df.values
    except AttributeError:
        log.debug("solo_model.get_scores() not available, trying fallback APIs.")

    if scores is None:
        try:
            # Fallback: predict(soft=True) returns probability-like scores
            soft_predictions = solo_model.predict(soft=True)
            if isinstance(soft_predictions, pd.Series):
                scores = soft_predictions.values
            else:
                scores = np.asarray(soft_predictions)
        except Exception as e:
            log.warning(f"Could not retrieve Solo scores via fallback: {e}. Using binary predictions as scores.")
            scores = predicted.astype(float)

    # #############################################

    doublet_count = sum(predicted)
    doublet_rate = doublet_count / len(predicted) if len(predicted) > 0 else 0
    log.info(f"  Found {doublet_count} potential doublets via Solo ({doublet_rate:.2%})")

    if use_gpu_flag and config.solo_clear_cache:
        torch.cuda.empty_cache()

    return scores, predicted


def _run_doubletdetection(
    adata_view: AnnData,
    sample_name: str,
    config: DoubletConfig,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Run DoubletDetection algorithm on a single AnnData view.

    Args:
        adata_view: AnnData object for a single sample.
        sample_name: Name of the sample being processed.
        config: The main DoubletConfig object.

    Returns:
        A tuple containing (scores, predicted) numpy arrays.
    """
    try:
        import doubletdetection as dd
    except ImportError:
        log.error(
            "DoubletDetection is not installed. Please install it to use the 'doubletdetection' method: pip install doubletdetection"
        )
        return None, None

    # --- 1. Data Preparation ---
    # DoubletDetection works best with raw counts.
    adata_dd = adata_view.copy()
    if config.dd_use_raw and adata_dd.raw:
        log.info("  Using 'adata.raw' for DoubletDetection as configured.")
        raw_counts = adata_dd.raw.X
    else:
        log.info("  Using 'adata.X' for DoubletDetection.")
        raw_counts = adata_dd.X

    # Ensure matrix is not sparse, or convert if necessary (dd prefers numpy array)
    if hasattr(raw_counts, "toarray"):
        raw_counts = raw_counts.toarray()

    log.info("  Running DoubletDetection classifier...")

    # --- 2. Run DoubletDetection ---
    try:
        clf = dd.BoostClassifier(
            n_components=config.dd_n_components,
            n_top_var_genes=config.dd_n_top_var_genes,
        )

        # Fit and predict
        results = clf.fit(raw_counts).predict(
            p_thresh=config.dd_p_thresh, voter_thresh=config.dd_voter_thresh
        )

        # Extract scores and predictions
        # 'results' is the binary prediction array (1=doublet, 0=singlet)
        predicted = results == 1
        # The 'scores' are stored in the classifier object
        scores = clf.doublet_score()

        doublet_count = sum(predicted)
        doublet_rate = doublet_count / len(predicted) if len(predicted) > 0 else 0
        log.info(
            f"  Found {doublet_count} potential doublets via DoubletDetection ({doublet_rate:.2%})"
        )

        return scores, predicted

    except Exception as e:
        log.error(f"DoubletDetection failed for sample {sample_name}: {e}")
        return None, None
    finally:
        gc.collect()


def _run_scdblfinder(
    adata_view: AnnData,
    sample_name: str,
    config: DoubletConfig,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Run scDblFinder (via the pure-Python port pyscdblfinder) on a single sample.

    Parameters
    ----------
    adata_view : AnnData
        Per-sample expression data. Raw counts are expected in ``.X`` or in the
        layer configured by ``config.scdblfinder_use_raw``.
    sample_name : str
        Name of the sample being processed.
    config : DoubletConfig
        Configuration object; relevant fields are ``scdblfinder_*`` and
        ``expected_doublet_rate``.

    Returns:
    -------
    Tuple of (scores, predicted) numpy arrays, or (None, None) on failure.
    """
    try:
        from pyscdblfinder import ScDblFinder
    except ImportError:
        log.error(
            "pyscdblfinder is not installed. Please install it to use the 'scdblfinder' method: "
            "pip install pyscdblfinder"
        )
        return None, None

    # Data-quality guard
    adata_sdf = adata_view.copy()
    if config.scdblfinder_use_raw and adata_sdf.raw is not None:
        log.info("  Using 'adata.raw' for scDblFinder as configured.")
        adata_sdf = adata_sdf.raw.to_adata()
    else:
        log.info("  Using 'adata.X' for scDblFinder.")

    if not _raw_count_guard(adata_sdf, sample_name=sample_name, method="scDblFinder"):
        return None, None

    if adata_sdf.n_vars < 100:
        log.warning(
            f"Skipping scDblFinder for sample '{sample_name}': only {adata_sdf.n_vars} genes "
            f"(minimum 100 required for reliable doublet detection)."
        )
        return None, None

    # Determine expected doublet rate. scDblFinder-specific dbr takes priority
    # over the workflow-level expected rate when provided.
    expected_rate = (
        config.scdblfinder_dbr
        if config.scdblfinder_dbr is not None
        else config.expected_doublet_rate
    )
    if isinstance(expected_rate, dict):
        expected_rate = expected_rate.get(sample_name)

    try:
        log.info(f"  Running scDblFinder for sample '{sample_name}'...")
        finder = ScDblFinder(adata_sdf, random_state=config.random_state)
        finder.run(
            dbr=expected_rate,
            n_features=config.scdblfinder_nfeatures,
            dims=config.scdblfinder_dims,
            k=config.scdblfinder_k,
            include_pcs=config.scdblfinder_include_pcs,
            iter=config.scdblfinder_iter,
            verbose=False,
        )

        score_col = next(
            (
                col
                for col in ("scDblFinder_score", "scDblFinder.score", "scdblfinder_score")
                if col in adata_sdf.obs
            ),
            None,
        )
        class_col = next(
            (
                col
                for col in ("scDblFinder_class", "scDblFinder.class", "scdblfinder_class")
                if col in adata_sdf.obs
            ),
            None,
        )
        if score_col is None or class_col is None:
            raise KeyError(
                "pyscdblfinder did not add expected score/class columns to adata.obs"
            )

        scores = adata_sdf.obs[score_col].to_numpy(dtype=float)
        predicted = adata_sdf.obs[class_col].astype(str).str.lower().eq("doublet").to_numpy()

        doublet_count = int(predicted.sum())
        doublet_rate = doublet_count / len(predicted) if len(predicted) > 0 else 0.0
        log.info(f"  Found {doublet_count} potential doublets via scDblFinder ({doublet_rate:.2%})")

        return scores, predicted

    except Exception as e:
        log.error(f"scDblFinder failed for sample '{sample_name}': {e}")
        return None, None
    finally:
        gc.collect()
