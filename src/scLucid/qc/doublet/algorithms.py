"""External doublet detection algorithm wrappers.

Extracted from core.py for maintainability.
"""

from __future__ import annotations

import gc
import logging
from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sparse
from anndata import AnnData

from ..config import DoubletConfig

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

    try:
        import scrublet as scr

        scrub = scr.Scrublet(adata_view.X, expected_doublet_rate=current_rate)
        scores, predicted = scrub.scrub_doublets(n_prin_comps=actual_n_pcs, verbose=False)

        if predicted is None:
            try:
                predicted = scrub.call_doublets(verbose=False)
            except Exception as e:
                log.warning(f"Scrublet call_doublets failed for sample '{sample_name}': {e}")

        if predicted is None and scores is not None:
            threshold = np.quantile(scores, max(0.0, min(1.0, 1.0 - float(current_rate))))
            log.warning(
                "Scrublet did not return binary predictions for sample '%s'; "
                "falling back to expected-rate score quantile threshold %.4f.",
                sample_name,
                threshold,
            )
            predicted = np.asarray(scores) > threshold

        if predicted is None:
            log.warning("Scrublet produced no usable scores or predictions for sample '%s'.", sample_name)
            return None, None

        doublet_count = sum(predicted)
        doublet_rate = doublet_count / len(predicted)
        log.info(f"  Found {doublet_count} potential doublets via Scrublet ({doublet_rate:.2%})")

        if config.scr_plot_umap:
            try:
                scrub.set_embedding("UMAP", scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
                fig, ax = scrub.plot_embedding("UMAP", order_points=True)
                if config.save_dir:
                    save_path = Path(config.save_dir) / f"{sample_name}_doublets_umap.png"
                    save_path.parent.mkdir(parents=True, exist_ok=True)
                    fig.savefig(save_path, dpi=300, bbox_inches="tight")
                if config.show_plots:
                    plt.show()
                else:
                    plt.close(fig)
            except Exception as e:
                log.warning(f"Could not generate UMAP for sample {sample_name}: {e}")

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
