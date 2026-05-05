"""External doublet detection algorithm wrappers.

Extracted from core.py for maintainability.
"""

from __future__ import annotations

import gc
import importlib
import logging
from typing import Optional, Tuple

import numpy as np
import pandas as pd
from anndata import AnnData

from ..config import DoubletConfig

log = logging.getLogger(__name__)

def _run_scrublet(
    adata_view: AnnData,
    sample_name: str,
    config: DoubletConfig,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Run Scrublet algorithm for doublet detection on a single AnnData view.
    Returns (scores, predicted) arrays.
    """
    if not hasattr(np.ndarray, "ptp"):
        np.ndarray.ptp = np.ptp

    rate = config.expected_doublet_rate
    current_rate = rate.get(sample_name, 0.1) if isinstance(rate, dict) else rate
    if current_rate is None:  # Handle case where rate is not provided for a sample
        log.warning(f"No doublet rate provided for sample '{sample_name}', using default of 0.1.")
        current_rate = 0.1

    # Data-quality guard: skip samples with too few features or too-low counts
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
        scores, _ = scrub.scrub_doublets(n_prin_comps=actual_n_pcs, verbose=False)
        predicted = scrub.call_doublets(verbose=False)

        if predicted is None:
            log.warning(
                f"Scrublet call_doublets returned None for sample '{sample_name}' "
                f"(simulated doublets may be too similar to real cells). "
                f"Falling back to heuristic-only for this sample."
            )
            return scores, np.zeros(scores.shape, dtype=bool) if scores is not None else None

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
        # Add a check to prevent excessive console output from the trainer
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

    # The scores are now retrieved using .get_scores()
    scores_df = solo_model.get_scores()
    scores = scores_df["doublet_scores"].values

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


