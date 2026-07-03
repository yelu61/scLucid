"""
Batch effect correction and data integration for single-cell RNA-seq.

This module provides high-level and low-level wrappers for common integration
methods: Harmony, Scanorama, scVI, BBKNN, ComBat. It ensures consistent API,
robust logging, and complete traceability for reproducible single-cell workflows.
"""

import logging

logging.getLogger("harmonypy").setLevel(logging.ERROR)
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

import harmonypy as hm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from scLucid.plotting.plotting_utils import _show_or_close

from .config import IntegrationConfig, apply_config_overrides

# Logging config
log = logging.getLogger(__name__)

__all__ = [
    "batch_correction",
    "decide_integration",
    "diagnose_integration_risk",
    "detect_integration_confounding",
    "evaluate_integration",
]


def _cramers_v(x: pd.Series, y: pd.Series) -> float:
    table = pd.crosstab(x.astype(str), y.astype(str))
    if table.empty or min(table.shape) < 2:
        return 0.0
    from scipy.stats import chi2_contingency

    chi2 = chi2_contingency(table)[0]
    n = table.to_numpy().sum()
    phi2 = chi2 / max(1, n)
    r, k = table.shape
    return float(np.sqrt(phi2 / max(1, min(k - 1, r - 1))))


def _silhouette_score(
    adata: AnnData,
    *,
    rep: Optional[str],
    labels: pd.Series,
    metric: str = "euclidean",
    sample_size: Optional[int] = None,
    random_state: int = 42,
) -> float:
    """Shared sampled silhouette helper for integration diagnostics."""
    from sklearn.metrics import silhouette_score

    if rep is None:
        X = adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
    else:
        if rep not in adata.obsm:
            raise KeyError(f"Representation '{rep}' not found in adata.obsm")
        X = adata.obsm[rep]
    labels = labels.astype(str)
    if labels.nunique() < 2:
        raise ValueError("Silhouette requires at least two label groups")
    active_sample_size = sample_size or min(adata.n_obs, max(100, min(5000, adata.n_obs)))
    return float(
        silhouette_score(
            X,
            labels,
            metric=metric,
            sample_size=min(active_sample_size, adata.n_obs),
            random_state=random_state,
        )
    )


def _interpret_integration_evaluation(results: Dict[str, Any]) -> Dict[str, Any]:
    """Translate integration metrics into a compact reviewer-facing summary."""
    warnings: List[str] = []
    strengths: List[str] = []

    batch_sil = results.get("batch_silhouette")
    if batch_sil is not None:
        if batch_sil >= -0.05:
            strengths.append("Batch silhouette suggests acceptable batch mixing.")
        elif batch_sil < -0.25:
            warnings.append("Batch silhouette remains strongly separated; inspect residual batch effects.")
        else:
            warnings.append("Batch silhouette is intermediate; compare integrated and unintegrated embeddings.")

    local_batch_acceptance = results.get("local_batch_chi2_acceptance")
    if local_batch_acceptance is not None:
        if local_batch_acceptance >= 0.7:
            strengths.append(
                "Local batch chi-square acceptance suggests local batch mixing is acceptable."
            )
        elif local_batch_acceptance < 0.3:
            warnings.append(
                "Local batch chi-square acceptance is low; local neighborhoods may still be batch-biased."
            )
        else:
            warnings.append(
                "Local batch chi-square acceptance is intermediate; review local batch mixing visually."
            )

    label_sil = results.get("label_silhouette")
    if label_sil is not None:
        if label_sil >= 0.15:
            strengths.append("Label silhouette suggests biological labels remain separated.")
        elif label_sil < 0.05:
            warnings.append("Label silhouette is low; check for overcorrection or loss of biological structure.")

    graph_connectivity = results.get("graph_connectivity")
    if graph_connectivity is not None:
        if graph_connectivity >= 0.8:
            strengths.append("Graph connectivity suggests labels remain connected after integration.")
        elif graph_connectivity < 0.5:
            warnings.append("Graph connectivity is low; label neighborhoods may be fragmented.")

    if warnings:
        status = "review"
        recommendation = "Review integration with unintegrated controls and biological marker preservation."
    else:
        status = "acceptable"
        recommendation = "Integration metrics do not show obvious problems; still validate marker preservation."

    return {
        "status": status,
        "warnings": warnings,
        "strengths": strengths,
        "recommendation": recommendation,
        "metric_notes": {
            "batch_silhouette": "Stored as negative silhouette; values closer to 0 are generally better for batch mixing.",
            "local_batch_chi2_acceptance": "Higher values indicate better local batch mixing.",
            "label_silhouette": "Higher values indicate better preservation of biological labels.",
            "graph_connectivity": "Higher values indicate better within-label graph connectivity.",
        },
    }


def detect_integration_confounding(
    adata: AnnData,
    batch_key: str,
    biology_columns: List[str],
    *,
    association_threshold: float = 0.8,
) -> List[str]:
    """Check whether ``batch_key`` is confounded with biology columns.

    Strong batch-biology association can cause integration to remove biological
    signal. One-to-one mappings are reported explicitly because they are the
    clearest high-risk case.

    Parameters
    ----------
    adata
        AnnData object.
    batch_key
        Column in ``.obs`` used for batch correction.
    biology_columns
        Columns representing biological signal that should not be removed.

    Returns:
    -------
    List[str]
        Problem descriptions. Empty list if no confounding is detected.
    """
    problems = []
    if batch_key not in adata.obs.columns:
        return [f"batch_key '{batch_key}' not in adata.obs"]

    for col in biology_columns:
        if col not in adata.obs.columns:
            continue
        association = _cramers_v(adata.obs[batch_key], adata.obs[col])
        if association >= association_threshold:
            problems.append(
                f"{batch_key} is strongly associated with {col} "
                f"(Cramer's V={association:.2f}); integration may remove biological signal"
            )
        mapping = adata.obs[[batch_key, col]].dropna().drop_duplicates()
        per_batch = mapping.groupby(batch_key)[col].nunique()
        if len(per_batch) > 0 and (per_batch <= 1).all():
            n_batch = adata.obs[batch_key].nunique()
            n_col = adata.obs[col].nunique()
            if n_batch == n_col:
                problems.append(
                    f"{batch_key} is one-to-one with {col} — integration would remove biological signal"
                )
    return problems


def diagnose_integration_risk(
    adata: AnnData,
    *,
    batch_key: str,
    condition_key: Optional[str] = None,
    biology_columns: Optional[List[str]] = None,
    label_key: Optional[str] = None,
    tumor: bool = False,
    clone_key: Optional[str] = None,
    state_key: Optional[str] = None,
    before_rep: str = "X_pca",
    after_rep: Optional[str] = None,
    confounding_threshold: float = 0.8,
    preservation_drop_threshold: float = 0.15,
    key_added: str = "integration_risk",
    record: bool = False,
) -> Dict[str, Union[str, float, bool, List[str], Dict[str, float]]]:
    """Diagnose when batch integration may be scientifically unsafe."""
    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not in adata.obs")

    warnings: List[str] = []
    metrics: Dict[str, float] = {}
    risk_score = 0.0

    if condition_key and condition_key in adata.obs.columns:
        confounding = _cramers_v(adata.obs[batch_key], adata.obs[condition_key])
        metrics["batch_condition_cramers_v"] = confounding
        if confounding >= confounding_threshold:
            warnings.append(
                "Batch and biological condition are strongly confounded; integration may remove biology."
            )
            risk_score += 0.4

    for problem in detect_integration_confounding(
        adata, batch_key, biology_columns or []
    ):
        warnings.append(problem)
        risk_score += 0.4

    if label_key and label_key in adata.obs.columns and before_rep in adata.obsm and after_rep in adata.obsm:
        try:
            labels = adata.obs[label_key].astype(str)
            sample_size = min(adata.n_obs, max(100, min(5000, adata.n_obs)))
            before_sil = _silhouette_score(
                adata, rep=before_rep, labels=labels, sample_size=sample_size
            )
            after_sil = _silhouette_score(
                adata, rep=after_rep, labels=labels, sample_size=sample_size
            )
            drop = float(before_sil - after_sil)
            metrics["label_silhouette_before"] = float(before_sil)
            metrics["label_silhouette_after"] = float(after_sil)
            metrics["label_silhouette_drop"] = drop
            if drop >= preservation_drop_threshold:
                warnings.append(
                    "Cell type/state separation decreased after integration; inspect overcorrection."
                )
                risk_score += 0.3
        except Exception as exc:
            warnings.append(f"Label preservation diagnostic failed: {exc}")

    tumor_keys = [key for key in [clone_key, state_key] if key and key in adata.obs.columns]
    if tumor or tumor_keys:
        warnings.append(
            "Tumor data can contain clone/state structure that Harmony/scVI-like integration may smooth away."
        )
        risk_score += 0.2
        for key in tumor_keys:
            metrics[f"n_{key}"] = float(adata.obs[key].nunique())

    if not warnings:
        recommendation = "Integration risk appears low; still compare pre/post cell type preservation."
        risk_level = "low"
    else:
        risk_level = "high" if risk_score >= 0.6 else "moderate"
        recommendation = "Run integration only with explicit biological-preservation checks and unintegrated controls."

    result = {
        "available": True,
        "risk_level": risk_level,
        "risk_score": float(min(1.0, risk_score)),
        "batch_key": batch_key,
        "condition_key": condition_key,
        "biology_columns": biology_columns or [],
        "label_key": label_key,
        "before_rep": before_rep,
        "after_rep": after_rep,
        "metrics": metrics,
        "warnings": warnings,
        "recommendation": recommendation,
    }
    if record:
        adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {}).setdefault(
            "integration", {}
        )[key_added] = result
    return result


def decide_integration(
    adata: AnnData,
    *,
    batch_key: str,
    run_integration: Union[bool, Literal["auto"]] = "auto",
    risk: Optional[Dict[str, Any]] = None,
    biology_columns: Optional[List[str]] = None,
    condition_key: Optional[str] = None,
    tumor: bool = False,
    before_rep: str = "X_pca",
    low_risk_only: bool = True,
) -> Tuple[bool, List[str], Optional[Dict[str, Any]]]:
    """Decide whether to run batch integration.

    Parameters
    ----------
    adata
        AnnData object with PCA already computed.
    batch_key
        Column in ``.obs`` identifying batches.
    run_integration
        ``True``/``False`` forces the decision. ``"auto"`` uses
        ``risk`` when provided, otherwise runs ``diagnose_integration_risk``.
    risk
        Optional precomputed result from ``diagnose_integration_risk``. Passing
        this keeps diagnosis and decision explicitly separated.
    biology_columns
        Columns representing biological signal that should not be removed by
        integration.
    condition_key
        Optional condition column checked for confounding.
    tumor
        Whether the data is tumor tissue (raises integration caution).
    before_rep
        Embedding key used for risk diagnosis.
    low_risk_only
        In ``"auto"`` mode, only run integration when risk is ``"low"``.

    Returns:
    -------
    (run, warnings, risk_dict)
        ``run`` is the boolean decision, ``warnings`` explains the decision, and
        ``risk_dict`` is the full risk assessment (None when ``run_integration``
        is a plain bool).

    Examples:
    --------
    >>> run, warnings, risk = decide_integration(
    ...     adata, batch_key="batch", biology_columns=["group"], run_integration="auto"
    ... )
    """
    if run_integration is False:
        return False, ["disabled by user"], None

    if run_integration is True:
        return True, [], None

    # auto mode
    warnings: List[str] = []

    if batch_key not in adata.obs.columns:
        return False, [f"batch_key '{batch_key}' not found in adata.obs"], None

    if adata.obs[batch_key].nunique() < 2:
        return False, [f"single technical batch in '{batch_key}' — no integration needed"], None

    biology_columns = biology_columns or []
    if risk is None:
        risk = diagnose_integration_risk(
            adata,
            batch_key=batch_key,
            condition_key=condition_key,
            biology_columns=biology_columns,
            tumor=tumor,
            before_rep=before_rep,
        )

    risk_level = risk.get("risk_level", "high")
    warnings.extend(risk.get("warnings", []))
    if low_risk_only and risk_level != "low":
        warnings.append(
            f"Integration risk level is '{risk_level}'; skipping integration in low_risk_only mode."
        )
        return False, warnings, risk

    if risk_level == "high":
        warnings.append(
            "Integration risk is high; run only with explicit biological-preservation checks."
        )
        return False, warnings, risk

    return True, warnings, risk


# ==============================================================================
# Low-Level Integration Wrappers (private, not in __all__)
# ==============================================================================
def _integrate_harmony(
    adata: AnnData,
    covariate_keys: Union[str, List[str]],
    basis: str = "X_pca",
    embedding_key: str = "X_harmony",
    max_iter_harmony: int = 50,
    theta: float = 2.0,
    lambda_val: float = 1.0,
    sigma: float = 0.1,
    random_state: int = 42,
    plot_convergence: bool = False,
    copy: bool = False,
    check_convergence: bool = True,
    convergence_threshold: float = 1e-4,
    **kwargs,
) -> AnnData:
    """
    Wrapper for Harmony batch correction using the core `harmonypy` library.
    Adds result to adata.obsm[embedding_key].
    This version supports single or multiple covariates.
    """
    if isinstance(covariate_keys, str):
        covariate_keys = [covariate_keys]

    missing_keys = [key for key in covariate_keys if key not in adata.obs]
    if missing_keys:
        raise ValueError(f"Covariate keys {missing_keys} not found in adata.obs")

    if basis not in adata.obsm:
        raise ValueError(f"Basis '{basis}' not found in adata.obsm (run PCA first)")

    if copy:
        adata = adata.copy()
    n_covariates = len(covariate_keys)
    log.info(
        f"Running Harmony integration on {n_covariates} covariate(s): {', '.join(covariate_keys)}"
    )
    log.info(
        f"Harmony params: basis={basis}, theta={theta}, lambda={lambda_val}, sigma={sigma}, max_iter={max_iter_harmony}"
    )
    # Initialize convergence tracking variables
    final_change: Optional[float] = None
    converged: Optional[bool] = None

    # Run Harmony using the harmonypy core function
    harmony_out = hm.run_harmony(
        data_mat=adata.obsm[basis],
        meta_data=adata.obs,
        vars_use=covariate_keys,
        theta=theta,
        lamb=lambda_val,
        sigma=sigma,
        max_iter_harmony=max_iter_harmony,
        random_state=random_state,
        plot_convergence=plot_convergence,
        **kwargs,
    )

    # === Validation ===
    Z_corr = harmony_out.Z_corr.T

    # 1. Check for NaN/Inf
    if np.any(~np.isfinite(Z_corr)):
        n_invalid = np.sum(~np.isfinite(Z_corr))
        raise RuntimeError(
            f"[preprocess] Harmony integration failed: output contains {n_invalid} NaN/Inf values. "
            "This may indicate convergence failure. Try adjusting theta or lambda."
        )

    # 2. Check convergence
    if check_convergence and hasattr(harmony_out, "objective_history"):
        obj_history = harmony_out.objective_history

        if len(obj_history) >= 2:
            final_change = abs(obj_history[-1] - obj_history[-2])
            converged = final_change <= convergence_threshold

            if not converged:
                log.warning(
                    f"⚠️  Harmony may not have fully converged. "
                    f"Final objective change: {final_change:.2e} > threshold {convergence_threshold:.2e}. "
                    f"Consider increasing max_iter_harmony (current: {max_iter_harmony})."
                )
            else:
                log.info(f"✓ Harmony converged. Final change: {final_change:.2e}")
        elif len(obj_history) < 2:
            log.info("Harmony converged in < 2 iterations. Setting convergence as True.")
            converged = True  # if less than 2 iterations, consider converged

    # 3. Check variance preservation
    input_var = np.var(adata.obsm[basis], axis=0).sum()
    output_var = np.var(Z_corr, axis=0).sum()
    var_ratio = output_var / input_var

    if var_ratio < 0.5:
        log.warning(
            f"⚠️  Harmony reduced total variance by {(1-var_ratio)*100:.1f}%. "
            "This may indicate over-correction. Consider reducing theta."
        )

    log.info(f"Variance retention: {var_ratio:.1%}")

    # Store result
    adata.obsm[embedding_key] = Z_corr

    # Store metadata for reproducibility
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {}).setdefault("integration", {})[
        "harmony"
    ] = {
        "covariate_keys": covariate_keys,
        "params": {
            "theta": theta,
            "lamb": lambda_val,
            "sigma": sigma,
            "max_iter_harmony": max_iter_harmony,
            "random_state": random_state,
        },
        "input_dims": adata.obsm[basis].shape[1],
        "output_dims": Z_corr.shape[1],
        "variance_retention": float(var_ratio),
        "converged": converged,
        "date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
    }

    log.info(
        f"Harmony integration finished: result stored in .obsm['{embedding_key}'] with shape {adata.obsm[embedding_key].shape}"
    )

    return adata


def _integrate_scanorama(
    adata: AnnData,
    batch_key: str,
    hvg: Optional[List[str]] = None,
    dims: int = 50,
    embedding_key: str = "X_scanorama",
    layer: Optional[str] = None,
    correct_expression: bool = False,
    return_corrected_expression: bool = False,
    knn: int = 20,
    sigma: float = 15,
    approx: bool = True,
    alpha: float = 0.1,
    batch_size: int = 5000,
    copy: bool = False,
    **kwargs,
) -> AnnData:
    """
    Wrapper for Scanorama integration. Use HVGs if provided, else all genes.
    Adds embedding to adata.obsm[embedding_key].
    """
    try:
        import scanorama
    except ImportError:
        log.error("Scanorama requires: pip install scanorama")
        raise ImportError("Please install scanorama")

    if batch_key not in adata.obs:
        raise ValueError(f"batch_key '{batch_key}' not in adata.obs")
    batches = adata.obs[batch_key].unique()
    n_batches = len(batches)
    if n_batches < 2:
        raise ValueError("Scanorama requires >=2 batches.")
    if copy:
        adata = adata.copy()

    # Split by batch
    adatas_list = []
    for b in batches:
        batch_data = adata[adata.obs[batch_key] == b].copy()
        if layer is not None:
            if layer not in batch_data.layers:
                raise ValueError(f"Layer '{layer}' not found in adata")
            batch_data.X = batch_data.layers[layer].copy()
        adatas_list.append(batch_data)

    # Subset to HVGs if provided
    if hvg is not None:
        genes_to_use = list(set(hvg).intersection(*[set(ad.var_names) for ad in adatas_list]))
        if len(genes_to_use) == 0:
            raise ValueError("No HVGs found in all batches.")
    else:
        genes_to_use = list(set.intersection(*[set(ad.var_names) for ad in adatas_list]))
    for i, ad in enumerate(adatas_list):
        adatas_list[i] = ad[:, genes_to_use].copy()
    log.info(f"Scanorama: {n_batches} batches, {len(genes_to_use)} genes.")

    # Integrate
    scanorama.integrate_scanpy(
        adatas_list,
        dimred=dims,
        knn=knn,
        sigma=sigma,
        approx=approx,
        alpha=alpha,
        batch_size=batch_size,
        **kwargs,
    )
    # Assemble embedding
    order = np.concatenate([np.where(adata.obs[batch_key] == b)[0] for b in batches])
    integrated = np.zeros((adata.shape[0], dims))
    idx = 0
    for i, b in enumerate(batches):
        n = adatas_list[i].n_obs
        integrated[order[idx : idx + n]] = adatas_list[i].obsm["X_scanorama"]
        idx += n
    adata.obsm[embedding_key] = integrated

    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {}).setdefault("integration", {})[
        "scanorama"
    ] = {
        "batch_key": batch_key,
        "params": {
            "knn": knn,
            "sigma": sigma,
            "approx": approx,
            "alpha": alpha,
            "batch_size": batch_size,
            "n_hvgs": len(genes_to_use),
            "layer": layer,
        },
        "output_dims": dims,
        "date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
    }
    log.info(
        f"Scanorama integration finished: {embedding_key} shape {adata.obsm[embedding_key].shape}"
    )
    return adata


def _integrate_scvi(
    adata: AnnData,
    batch_key: str,
    layer: Optional[str] = "counts",
    n_layers: int = 2,
    n_latent: int = 30,
    batch_size: int = 2560,
    max_epochs: int = 1000,
    embedding_key: str = "X_scVI",
    gene_likelihood: str = "nb",
    save_model: bool = False,
    model_path: Optional[str] = None,
    plan_kwargs: Optional[dict] = None,
    copy: bool = False,
    **kwargs,
) -> AnnData:
    """
    Wrapper for scVI integration. Adds latent to adata.obsm[embedding_key].
    """
    try:
        import scvi
    except ImportError:
        log.error("scvi-tools required: pip install scvi-tools")
        raise ImportError("Please install scvi-tools")

    if batch_key not in adata.obs:
        raise ValueError(f"batch_key '{batch_key}' not in adata.obs")
    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        log.warning("scVI is designed for >1 batch.")
    if copy:
        adata = adata.copy()
    if layer is not None and layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' not in adata.layers")

    scvi.model.SCVI.setup_anndata(adata, layer=layer, batch_key=batch_key)
    plan_kwargs = plan_kwargs or {}
    model = scvi.model.SCVI(
        adata, n_layers=n_layers, n_latent=n_latent, gene_likelihood=gene_likelihood
    )
    model.train(
        batch_size=batch_size,
        max_epochs=max_epochs,
        early_stopping=True,
        plan_kwargs=plan_kwargs,
        **kwargs,
    )
    adata.obsm[embedding_key] = model.get_latent_representation()
    if save_model:
        if not model_path:
            raise ValueError("model_path must be provided when save_model=True")
        # Convert string path to a Path object
        save_path = Path(model_path)
        # Create the parent directory if it doesn't exist
        save_path.parent.mkdir(parents=True, exist_ok=True)
        # Save the model
        model.save(str(save_path), overwrite=True)
        log.info(f"scVI model saved to: {save_path}")
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {}).setdefault("integration", {})[
        "scvi"
    ] = {
        "batch_key": batch_key,
        "params": {
            "n_layers": n_layers,
            "n_latent": n_latent,
            "batch_size": batch_size,
            "max_epochs": max_epochs,
            "gene_likelihood": gene_likelihood,
            "layer": layer,
        },
        "output_dims": n_latent,
        "model_saved": save_model,
        "model_path": model_path if save_model else None,
        "date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
    }
    log.info(f"scVI integration finished: {embedding_key} shape {adata.obsm[embedding_key].shape}")
    return adata


def _integrate_scanvi(
    adata: AnnData,
    batch_key: str,
    labels_key: str,
    layer: Optional[str] = "counts",
    n_layers: int = 2,
    n_latent: int = 30,
    batch_size: int = 2560,
    max_epochs: int = 500,
    early_stopping_patience: int = 20,
    embedding_key: str = "X_scANVI",
    gene_likelihood: str = "nb",
    save_model: bool = False,
    model_path: Optional[str] = None,
    plan_kwargs: Optional[dict] = None,
    copy: bool = False,
    **kwargs,
) -> AnnData:
    """
    Wrapper for scANVI (scVI with label supervision) integration.

    scANVI leverages cell type labels to improve batch correction,
    particularly useful when biological populations are confounded with batch.

    Args:
        adata: AnnData object.
        batch_key: Column in adata.obs for batch labels.
        labels_key: Column in adata.obs for cell type labels (required).
        layer: Layer with raw counts.
        n_layers: Number of hidden layers in VAE.
        n_latent: Latent dimension.
        batch_size: Training batch size.
        max_epochs: Maximum training epochs.
        early_stopping_patience: Early stopping patience.
        embedding_key: Key for storing the embedding in adata.obsm.
        gene_likelihood: Gene likelihood model.
        save_model: Whether to save the trained model.
        model_path: Path to save the model.
        plan_kwargs: Additional training plan kwargs.
        copy: Whether to copy adata.
        **kwargs: Additional kwargs passed to scANVI.

    Returns:
        AnnData with scANVI embedding in adata.obsm[embedding_key].
    """
    try:
        import scvi
    except ImportError:
        log.error("scvi-tools required: pip install scvi-tools")
        raise ImportError("Please install scvi-tools")

    if batch_key not in adata.obs:
        raise ValueError(f"batch_key '{batch_key}' not in adata.obs")
    if labels_key not in adata.obs:
        raise ValueError(
            f"labels_key '{labels_key}' not in adata.obs. "
            f"scANVI requires cell type labels. Available columns: {list(adata.obs.columns)}"
        )

    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        log.warning("scANVI is designed for >1 batch.")

    if copy:
        adata = adata.copy()

    if layer is not None and layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' not in adata.layers")

    # Ensure labels are categorical
    if not hasattr(adata.obs[labels_key], "cat"):
        adata.obs[labels_key] = adata.obs[labels_key].astype("category")

    scvi.model.SCANVI.setup_anndata(
        adata,
        layer=layer,
        batch_key=batch_key,
        labels_key=labels_key,
        unlabeled_category="Unknown",
    )

    plan_kwargs = plan_kwargs or {}
    model = scvi.model.SCANVI(
        adata,
        n_layers=n_layers,
        n_latent=n_latent,
        gene_likelihood=gene_likelihood,
    )

    # Train with early stopping
    model.train(
        batch_size=batch_size,
        max_epochs=max_epochs,
        early_stopping=True,
        early_stopping_patience=early_stopping_patience,
        plan_kwargs=plan_kwargs,
        **kwargs,
    )

    adata.obsm[embedding_key] = model.get_latent_representation()

    if save_model:
        if not model_path:
            raise ValueError("model_path must be provided when save_model=True")
        save_path = Path(model_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        model.save(str(save_path), overwrite=True)
        log.info(f"scANVI model saved to: {save_path}")

    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {}).setdefault("integration", {})[
        "scanvi"
    ] = {
        "batch_key": batch_key,
        "labels_key": labels_key,
        "params": {
            "n_layers": n_layers,
            "n_latent": n_latent,
            "batch_size": batch_size,
            "max_epochs": max_epochs,
            "early_stopping_patience": early_stopping_patience,
            "gene_likelihood": gene_likelihood,
            "layer": layer,
        },
        "output_dims": n_latent,
        "model_saved": save_model,
        "model_path": model_path if save_model else None,
        "date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
    }
    log.info(f"scANVI integration finished: {embedding_key} shape {adata.obsm[embedding_key].shape}")
    return adata


def _integrate_bbknn(
    adata: AnnData,
    batch_key: str,
    use_rep: str = "X_pca",
    neighbors_within_batch: int = 3,
    n_pcs: Optional[int] = None,
    trim: Optional[int] = None,
    copy: bool = False,
    **kwargs,
) -> AnnData:
    """
    Wrapper for BBKNN batch-corrected graph.
    """
    try:
        from scanpy.external.pp import bbknn
    except ImportError:
        log.error("BBKNN requires: pip install bbknn")
        raise ImportError("Please install bbknn")
    if batch_key not in adata.obs:
        raise ValueError(f"batch_key '{batch_key}' not in adata.obs")
    if use_rep not in adata.obsm:
        raise ValueError(f"use_rep '{use_rep}' not in adata.obsm")
    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        log.warning("BBKNN is designed for >1 batch.")
    if copy:
        adata = adata.copy()
    if n_pcs is None:
        n_pcs = adata.obsm[use_rep].shape[1]
    bbknn(
        adata,
        batch_key=batch_key,
        use_rep=use_rep,
        neighbors_within_batch=neighbors_within_batch,
        n_pcs=n_pcs,
        trim=trim,
        **kwargs,
    )
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {}).setdefault("integration", {})[
        "bbknn"
    ] = {
        "batch_key": batch_key,
        "params": {
            "neighbors_within_batch": neighbors_within_batch,
            "n_pcs": n_pcs,
            "trim": trim,
            "use_rep": use_rep,
        },
        "date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
    }
    log.info("BBKNN integration finished; neighborhood graph updated.")
    return adata


def _integrate_combat(
    adata: AnnData,
    batch_key: str,
    layer: Optional[str] = None,
    covariates: Optional[List[str]] = None,
    inplace: bool = True,
    output_layer: str = "combat_corrected",
    copy: bool = False,
    force_dense: bool = False,
    **kwargs,
) -> AnnData:
    """
    Wrapper for ComBat batch correction.

    Args:
        force_dense: If True, allow ComBat to densify a sparse input matrix.
            If False and the input is sparse, a warning is issued for small
            matrices and a ``MemoryError`` is raised for large matrices
            (>50 million elements or non-zeros), because ComBat requires a
            dense array. Consider Harmony or scVI integration for large datasets.
    """
    if batch_key not in adata.obs:
        raise ValueError(f"batch_key '{batch_key}' not in adata.obs")
    if covariates is not None:
        missing_covariates = [c for c in covariates if c not in adata.obs.columns]
        if missing_covariates:
            raise ValueError(f"Covariates not found in adata.obs: {missing_covariates}")
    if layer is not None and layer not in adata.layers:
        raise ValueError(f"layer '{layer}' not found in adata.layers")
    if copy:
        adata = adata.copy()
    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        log.warning("ComBat is designed for >1 batch.")
        return adata
    if layer is not None:
        X_combat = adata.layers[layer].copy()
    else:
        X_combat = adata.X.copy()
    import scipy.sparse

    if scipy.sparse.issparse(X_combat):
        n_elements = X_combat.shape[0] * X_combat.shape[1]
        n_nonzero = X_combat.nnz
        large_threshold = 50_000_000
        is_large = n_elements > large_threshold or n_nonzero > large_threshold

        if is_large and not force_dense:
            raise MemoryError(
                f"ComBat would densify a {X_combat.shape[0]}x{X_combat.shape[1]} sparse matrix "
                f"({n_elements} elements, {n_nonzero} non-zeros). Set force_dense=True to "
                f"proceed, or use Harmony/scVI integration for large datasets."
            )

        if is_large:
            log.warning(
                "ComBat is densifying a large sparse matrix (%d elements, %d non-zeros). "
                "Consider Harmony or scVI integration for large datasets.",
                n_elements,
                n_nonzero,
            )
        else:
            log.warning(
                "ComBat is densifying a sparse input matrix. "
                "Consider Harmony or scVI integration for large datasets."
            )

        X_combat = X_combat.toarray()
    import anndata

    temp_adata = anndata.AnnData(X=X_combat, obs=adata.obs.copy())
    sc.pp.combat(temp_adata, key=batch_key, covariates=covariates, inplace=True, **kwargs)
    X_corrected = temp_adata.X
    if inplace:
        adata.X = X_corrected
    if output_layer:
        adata.layers[output_layer] = X_corrected
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {}).setdefault("integration", {})[
        "combat"
    ] = {
        "batch_key": batch_key,
        "params": {
            "covariates": covariates,
            "layer": layer,
            "output_layer": output_layer,
        },
        "date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
    }
    log.info("ComBat integration finished.")
    return adata


def _compute_local_batch_chi2_acceptance(
    X: np.ndarray,
    batch_labels: pd.Series,
    n_neighbors: int = 25,
    alpha: float = 0.05,
    n_sample_cells: int = 1000,
) -> Dict[str, Any]:
    """
    Compute an approximate local batch-mixing score using chi-square tests.

    This is a k-BET-like heuristic: for each sampled cell, the batch composition
    of its local neighborhood is compared to the global batch distribution using
    a chi-square goodness-of-fit test. It is useful as an audit diagnostic, but
    should not be reported as an exact reproduction of the original kBET package.

    Returns:
        Dict with 'rejection_rate' and 'acceptance_rate'.
    """
    from scipy.stats import chi2
    from sklearn.neighbors import NearestNeighbors

    n_cells = X.shape[0]

    # Sample cells if dataset is large
    rng = np.random.default_rng(42)
    if n_cells > n_sample_cells:
        sample_idx = rng.choice(n_cells, n_sample_cells, replace=False)
        X_sample = X[sample_idx]
    else:
        X_sample = X

    # Build k-NN graph
    nn = NearestNeighbors(n_neighbors=n_neighbors + 1)  # +1 to exclude self
    nn.fit(X)
    distances, indices = nn.kneighbors(X_sample)
    indices = indices[:, 1:]  # Remove self

    # Get global batch distribution
    batch_categories = sorted(batch_labels.unique())
    global_freq = batch_labels.value_counts(normalize=True).reindex(batch_categories).values

    # Test each neighborhood
    rejections = 0
    valid_tests = 0

    for i in range(len(X_sample)):
        # Get batch composition of neighborhood
        neighbor_batches = batch_labels.iloc[indices[i]].values

        # Count batches in neighborhood
        observed = np.array([(neighbor_batches == batch).sum() for batch in batch_categories])

        # Expected counts under null hypothesis
        expected = global_freq * n_neighbors

        # Skip if expected counts are too small
        if (expected < 5).any():
            continue

        # Chi-square test
        chi2_stat = np.sum((observed - expected) ** 2 / expected)

        # Degrees of freedom
        df = len(batch_categories) - 1

        # p-value
        p_value = 1 - chi2.cdf(chi2_stat, df)

        valid_tests += 1

        if p_value < alpha:
            rejections += 1

    if valid_tests == 0:
        log.warning("Local batch chi-square: No valid tests could be performed")
        return {
            "rejection_rate": np.nan,
            "acceptance_rate": np.nan,
            "n_tests": 0,
            "n_neighbors": int(n_neighbors),
            "alpha": float(alpha),
            "model_type": "chi_square_local_batch_approximation",
            "claim_level": "batch_mixing_diagnostic_heuristic",
            "review_note": (
                "No valid chi-square neighborhood tests were available; local batch "
                "mixing should be treated as unavailable for this run."
            ),
        }

    rejection_rate = rejections / valid_tests
    acceptance_rate = 1 - rejection_rate

    return {
        "rejection_rate": rejection_rate,
        "acceptance_rate": acceptance_rate,
        "n_tests": valid_tests,
        "n_neighbors": int(n_neighbors),
        "alpha": float(alpha),
        "model_type": "chi_square_local_batch_approximation",
        "claim_level": "batch_mixing_diagnostic_heuristic",
        "review_note": (
            "Local batch mixing is approximated with chi-square goodness-of-fit tests "
            "over sampled neighborhoods; inspect alongside silhouette, graph connectivity, "
            "and biology preservation evidence."
        ),
    }


def _compute_kbet_score(*args: Any, **kwargs: Any) -> Dict[str, Any]:
    """Deprecated alias for :func:`_compute_local_batch_chi2_acceptance`."""
    import warnings

    warnings.warn(
        "_compute_kbet_score is deprecated; use _compute_local_batch_chi2_acceptance.",
        FutureWarning,
        stacklevel=2,
    )
    return _compute_local_batch_chi2_acceptance(*args, **kwargs)


def _compute_kbet_reference(
    adata: AnnData,
    batch_key: str,
    use_rep: Optional[str] = None,
    n_neighbors: int = 25,
    **kwargs: Any,
) -> Dict[str, Any]:
    """
    Compute kBET using ``scib.metrics.kBET`` when available.

    Falls back to the local chi-square approximation with a clear warning if
    scib is not installed or if the reference computation fails.

    Returns:
        Dict matching the contract of :func:`_compute_local_batch_chi2_acceptance`,
        with an additional ``reference_backend`` key indicating whether the scib
        reference or the approximation was used.
    """
    try:
        import scib
    except ImportError:
        log.warning(
            "scib is not installed; falling back to the local chi-square approximation "
            "for kBET. Install scib for the reference implementation."
        )
        X = adata.obsm[use_rep] if use_rep and use_rep in adata.obsm else adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        result = _compute_local_batch_chi2_acceptance(
            X, adata.obs[batch_key], n_neighbors=n_neighbors, **kwargs
        )
        result["reference_backend"] = "chi_square_approximation"
        return result

    try:
        # scib.metrics.kBET returns a float acceptance rate by default.
        acceptance_rate = float(
            scib.metrics.kBET(
                adata,
                batch_key=batch_key,
                label_key=batch_key,
                type_=None if use_rep is None else "embed",
                embed=use_rep,
                k0=n_neighbors,
                **kwargs,
            )
        )
        return {
            "rejection_rate": 1.0 - acceptance_rate,
            "acceptance_rate": acceptance_rate,
            "n_tests": None,
            "n_neighbors": int(n_neighbors),
            "alpha": float(kwargs.get("alpha", 0.05)),
            "model_type": "scib_kbet_reference",
            "claim_level": "batch_mixing_diagnostic_reference",
            "reference_backend": "scib",
            "review_note": (
                "kBET computed via scib.metrics.kBET reference implementation."
            ),
        }
    except Exception as exc:
        log.warning(
            "scib kBET reference computation failed (%s); falling back to the local "
            "chi-square approximation.",
            exc,
        )
        X = adata.obsm[use_rep] if use_rep and use_rep in adata.obsm else adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        result = _compute_local_batch_chi2_acceptance(
            X, adata.obs[batch_key], n_neighbors=n_neighbors, **kwargs
        )
        result["reference_backend"] = "chi_square_approximation"
        return result


# ==============================================================================
# High-level entry: batch_correction
# ==============================================================================


def batch_correction(
    adata: AnnData,
    config: Optional[IntegrationConfig] = None,
    **kwargs,
) -> AnnData:
    """
    High-level, config-driven wrapper for batch correction/integration.

    Args:
        adata: AnnData object.
        config: An IntegrationConfig object. If None, a default config is used.
        **kwargs: Keyword arguments to override parameters in the config object
                  (e.g., `method='scanorama'`, `use_rep='X_pca'`).

    Returns:
        AnnData with integration result in .obsm.
    """
    # --- 1. Establish the final configuration ---
    if config is None:
        active_config = IntegrationConfig()
    else:
        active_config = apply_config_overrides(config, ignored_keys={"force"}, **kwargs)

    # --- 2. Extract parameters from the final config ---
    method = active_config.method.lower() if isinstance(active_config.method, str) else active_config.method
    batch_key = active_config.batch_key
    use_rep = active_config.use_rep
    output_key = active_config.output_key or f"X_{method}"
    force = kwargs.get("force", False)
    plot = active_config.plot
    save_dir = Path(active_config.save_dir) if active_config.save_dir else None

    # --- 3. ❗ ENHANCED Input validation ❗ ---
    if not method or not batch_key:
        log.info("`method` or `batch_key` not specified in config. Skipping batch correction.")
        return adata

    # Handle both string and list for batch_key
    if isinstance(batch_key, str):
        keys_to_check = [batch_key]
    elif isinstance(batch_key, list):
        keys_to_check = batch_key
    else:
        raise TypeError(f"batch_key must be a string or a list of strings, not {type(batch_key)}")

    missing_keys = [key for key in keys_to_check if key not in adata.obs]
    if missing_keys:
        # If the default batch_key is missing, try common alternatives
        if len(keys_to_check) == 1 and keys_to_check[0] == "sampleID":
            for alt in ["sample", "Sample", "batch", "Batch"]:
                if alt in adata.obs:
                    log.warning(
                        f"batch_key 'sampleID' not found. Using detected alternative: '{alt}'"
                    )
                    batch_key = alt
                    keys_to_check = [alt]
                    missing_keys = []
                    break
        if missing_keys:
            raise ValueError(
                f"batch_key(s) {missing_keys} not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )
    # --- END OF ENHANCEMENT ---

    if output_key in adata.obsm and not force:
        log.info(f"Integration result '{output_key}' already exists. Use force=True to rerun.")
        return adata

    integration_meta = (
        adata.uns.setdefault("sclucid", {})
        .setdefault("preprocess", {})
        .setdefault("integration", {})
    )

    decision_record = None
    if active_config.auto_decide:
        decision_batch_key = batch_key[0] if isinstance(batch_key, list) else batch_key
        risk = diagnose_integration_risk(
            adata,
            batch_key=decision_batch_key,
            condition_key=active_config.condition_key,
            biology_columns=active_config.biology_columns,
            label_key=active_config.label_key,
            tumor=active_config.tumor,
            before_rep=use_rep,
            record=True,
        )
        run, decision_warnings, risk = decide_integration(
            adata,
            batch_key=decision_batch_key,
            run_integration="auto",
            risk=risk,
            biology_columns=active_config.biology_columns,
            condition_key=active_config.condition_key,
            tumor=active_config.tumor,
            before_rep=use_rep,
        )
        decision_record = {
            "run": bool(run),
            "warnings": decision_warnings,
            "risk": risk,
        }
        integration_meta["workflow"] = {
            "params": active_config.to_dict(),
            "method": method,
            "batch_key": batch_key,
            "use_rep": use_rep,
            "output_key": output_key,
            "auto_decide": True,
            "decision": decision_record,
            "decision_warnings": decision_warnings,
            "risk": risk,
        }
        if not run:
            integration_meta["workflow"]["skipped"] = True
            integration_meta["workflow"]["skip_reason"] = "; ".join(decision_warnings)
            log.warning(
                "Skipping integration after auto decision: %s",
                integration_meta["workflow"]["skip_reason"],
            )
            return adata

    # --- 4. Plot before state (if requested) ---
    if plot:
        adata_before = adata.copy()
        if "neighbors" not in adata_before.uns and use_rep in adata_before.obsm:
            sc.pp.neighbors(adata_before, use_rep=use_rep)
        if "X_umap" not in adata_before.obsm:
            sc.tl.umap(adata_before)

    # --- 5. Main integration logic ---
    method_kwargs = {}
    if method == "harmony":
        method_kwargs = dict(active_config.harmony_params)
    elif method == "scvi":
        method_kwargs = dict(active_config.scvi_params)
    elif method == "scanvi":
        method_kwargs = dict(active_config.scanvi_params)
    elif method == "scanorama" and active_config.hvg_key:
        if active_config.hvg_key in adata.var:
            method_kwargs["hvg"] = adata.var_names[adata.var[active_config.hvg_key]].tolist()
        else:
            log.warning(
                f"hvg_key '{active_config.hvg_key}' not found in .var. Running on all genes."
            )

    # Merge user-supplied method_kwargs (takes precedence over method-specific defaults)
    if active_config.method_kwargs:
        overlap = set(method_kwargs.keys()) & set(active_config.method_kwargs.keys())
        if overlap:
            log.info(f"method_kwargs override default params: {overlap}")
        method_kwargs.update(active_config.method_kwargs)

    log.info(f"Running batch correction with method: '{method}'")

    if method == "harmony":
        adata = _integrate_harmony(
            adata,
            covariate_keys=batch_key,
            basis=use_rep,
            embedding_key=output_key,
            **method_kwargs,
        )
    elif method == "scanorama":
        adata = _integrate_scanorama(adata, batch_key, embedding_key=output_key, **method_kwargs)
    elif method == "scvi":
        adata = _integrate_scvi(adata, batch_key, embedding_key=output_key, **method_kwargs)
    elif method == "scanvi":
        if not active_config.scanvi_labels_key:
            raise ValueError(
                "scANVI requires 'scanvi_labels_key' in IntegrationConfig. "
                "Provide a key from adata.obs containing cell type labels."
            )
        adata = _integrate_scanvi(
            adata,
            batch_key=batch_key,
            labels_key=active_config.scanvi_labels_key,
            embedding_key=output_key,
            **method_kwargs,
        )
    elif method == "bbknn":
        adata = _integrate_bbknn(adata, batch_key, use_rep=use_rep, **method_kwargs)
    elif method == "combat":
        adata = _integrate_combat(adata, batch_key, **method_kwargs)
    else:
        raise ValueError(
            f"Unknown integration method '{method}'. "
            "Expected one of: harmony, scanorama, scvi, scanvi, bbknn, combat."
        )

    # --- 6. Store metadata and plot after state ---
    integration_meta["workflow"] = {
        "params": active_config.to_dict(),
        "method": method,
        "batch_key": batch_key,
        "use_rep": use_rep,
        "output_key": output_key,
        "auto_decide": bool(active_config.auto_decide),
        "skipped": False,
    }
    if decision_record is not None:
        integration_meta["workflow"]["decision"] = decision_record
        integration_meta["workflow"]["decision_warnings"] = decision_record["warnings"]
        integration_meta["workflow"]["risk"] = decision_record["risk"]

    log.info(f"Integration complete. Result stored in: adata.obsm['{output_key}']")

    if active_config.evaluate:
        eval_batch_key = batch_key[0] if isinstance(batch_key, list) else batch_key
        eval_rep = output_key if output_key in adata.obsm else use_rep
        integration_meta["evaluation"] = evaluate_integration(
            adata,
            batch_key=eval_batch_key,
            label_key=active_config.label_key,
            integration_method=method,
            use_rep=eval_rep,
            plot=False,
        )

    if plot:
        if method != "bbknn":
            sc.pp.neighbors(adata, use_rep=output_key)
        sc.tl.umap(adata)

        # If batch_key is a list, use the first element for coloring the UMAP
        color_key = batch_key[0] if isinstance(batch_key, list) else batch_key
        log.info(f"Using '{color_key}' for UMAP color annotation.")

        fig, axes = plt.subplots(1, 2, figsize=(16, 7))
        fig.suptitle("Batch Correction Comparison", fontsize=16)
        sc.pl.umap(
            adata_before,
            color=color_key,
            ax=axes[0],
            show=False,
            title=f"Before ({use_rep})",
        )
        sc.pl.umap(
            adata,
            color=color_key,
            ax=axes[1],
            show=False,
            title=f"After {method.capitalize()}",
        )

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

        if save_dir:
            save_dir.mkdir(parents=True, exist_ok=True)
            figure_path = save_dir / f"batch_correction_{method}.png"
            plt.savefig(figure_path, dpi=300)
            log.info(f"Saved correction plot to: {figure_path}")

        _show_or_close(fig)

    return adata


# ==============================================================================
# Integration Evaluation
# ==============================================================================


def evaluate_integration(
    adata: AnnData,
    batch_key: str,
    label_key: Optional[str] = None,
    integration_method: Optional[str] = None,
    use_rep: Optional[str] = None,
    n_neighbors: int = 30,
    metric: str = "euclidean",
    methods: List[str] = ["silhouette", "kbet", "graph_connectivity"],
    plot: bool = True,
    save_path: Optional[str] = None,
    use_reference_kbet: bool = True,
) -> Dict[str, float]:
    """
    Evaluate integration quality using batch/label silhouette, local batch
    chi-square, and graph connectivity. Automatically finds the best embedding
    if not given.

    Args:
        use_reference_kbet: If True, attempt to use ``scib.metrics.kBET`` as the
            reference implementation. Falls back to the local chi-square
            approximation if scib is unavailable or fails.

    Returns:
        Dict of scores (higher is better).
    """
    if batch_key not in adata.obs:
        raise ValueError(f"batch_key '{batch_key}' not in adata.obs")
    if label_key and label_key not in adata.obs:
        log.warning(f"label_key '{label_key}' not in adata.obs, ignoring.")
        label_key = None

    # Guess use_rep if not provided
    if use_rep is None:
        for key in ["X_harmony", "X_scanorama", "X_scVI"]:
            if key in adata.obsm:
                use_rep = key
                break
        if use_rep is None and "X_pca" in adata.obsm:
            use_rep = "X_pca"
    if use_rep is None:
        log.warning("No embedding found, using adata.X")
        X = adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
    else:
        X = adata.obsm[use_rep]
    results = {
        "method": integration_method or use_rep or "unknown",
        "n_batches": adata.obs[batch_key].nunique(),
        "n_cells": adata.n_obs,
    }
    if label_key:
        results["n_labels"] = adata.obs[label_key].nunique()

    # 1. Silhouette score (measures batch mixing)
    if "silhouette" in methods:
        try:
            # Data-adaptive sample size: at least 5k, up to 20k for large datasets,
            # but never more than the full dataset.
            sil_sample_size = min(adata.n_obs, max(5000, min(20000, adata.n_obs // 5)))
            batch_sil = _silhouette_score(
                adata,
                rep=use_rep,
                labels=adata.obs[batch_key],
                metric=metric,
                sample_size=sil_sample_size,
            )
            results["batch_silhouette"] = -batch_sil
            if label_key:
                label_sil = _silhouette_score(
                    adata,
                    rep=use_rep,
                    labels=adata.obs[label_key],
                    metric=metric,
                    sample_size=sil_sample_size,
                )
                results["label_silhouette"] = label_sil
                results["overall_silhouette"] = (results["batch_silhouette"] + label_sil) / 2
            log.info(f"Batch silhouette: {results['batch_silhouette']:.4f}")
        except Exception as e:
            log.warning(f"Silhouette failed: {e}")

    # 2. Local batch chi-square (k-BET-like batch mixing in local neighborhoods)
    if "kbet" in methods:
        try:
            log.info("Computing local batch mixing score...")
            if use_reference_kbet:
                kbet_result = _compute_kbet_reference(
                    adata,
                    batch_key=batch_key,
                    use_rep=use_rep,
                    n_neighbors=n_neighbors,
                    n_sample_cells=min(2000, adata.n_obs),
                )
            else:
                kbet_result = _compute_local_batch_chi2_acceptance(
                    X,
                    adata.obs[batch_key],
                    n_neighbors=n_neighbors,
                    n_sample_cells=min(2000, adata.n_obs),
                )

            results["local_batch_chi2_acceptance"] = kbet_result["acceptance_rate"]
            results["local_batch_chi2_rejection_rate"] = kbet_result["rejection_rate"]
            results["local_batch_chi2_model_type"] = kbet_result.get("model_type")
            results["local_batch_chi2_claim_level"] = kbet_result.get("claim_level")
            results["local_batch_chi2_review_note"] = kbet_result.get("review_note")
            results["local_batch_chi2_n_neighbors"] = kbet_result.get("n_neighbors")
            if "reference_backend" in kbet_result:
                results["local_batch_chi2_reference_backend"] = kbet_result[
                    "reference_backend"
                ]

            # Backward-compatible aliases (deprecated)
            import warnings

            warnings.warn(
                "The 'kbet_acceptance' and 'kbet_rejection_rate' keys are deprecated; "
                "use 'local_batch_chi2_acceptance' and "
                "'local_batch_chi2_rejection_rate' instead.",
                FutureWarning,
                stacklevel=2,
            )
            results["kbet_acceptance"] = kbet_result["acceptance_rate"]
            results["kbet_rejection_rate"] = kbet_result["rejection_rate"]
            results["kbet_model_type"] = kbet_result.get("model_type")
            results["kbet_claim_level"] = kbet_result.get("claim_level")
            results["kbet_review_note"] = kbet_result.get("review_note")
            results["kbet_n_neighbors"] = kbet_result.get("n_neighbors")

            log.info(f"Local batch chi-square acceptance rate: {kbet_result['acceptance_rate']:.4f}")
            log.info(f"Local batch chi-square tests performed: {kbet_result['n_tests']}")

        except Exception as e:
            log.warning(f"Failed to compute local batch mixing: {str(e)}")

    # 3. Graph connectivity (measures label preservation)
    if "graph_connectivity" in methods and label_key is not None:
        try:
            # Build a nearest-neighbor graph
            from sklearn.neighbors import kneighbors_graph

            # Create adjacency matrix
            adjacency = kneighbors_graph(
                X,
                n_neighbors=n_neighbors,
                mode="connectivity",
                metric=metric,
                include_self=False,
            )

            # For each label, compute the largest connected component size
            from scipy.sparse.csgraph import connected_components

            label_categories = (
                adata.obs[label_key].cat.categories
                if hasattr(adata.obs[label_key], "cat")
                else sorted(adata.obs[label_key].unique())
            )
            connectivities = []

            for label in label_categories:
                # Get indices for this label
                label_mask = adata.obs[label_key] == label
                label_indices = np.where(label_mask)[0]

                if len(label_indices) <= 1:
                    continue

                # Extract subgraph for this label
                subgraph = adjacency[label_indices][:, label_indices]

                # Find connected components
                n_components, component_labels = connected_components(subgraph, directed=False)

                # Compute largest component size
                component_sizes = np.bincount(component_labels)
                largest_component = component_sizes.max()

                # Compute connectivity score (ratio of largest component to total cells in this label)
                connectivity = largest_component / len(label_indices)
                connectivities.append(connectivity)

            if connectivities:
                # Average connectivity across all labels
                results["graph_connectivity"] = np.mean(connectivities)
                log.info(f"Graph connectivity: {results['graph_connectivity']:.4f}")
            else:
                log.warning("Could not compute graph connectivity: no valid labels")

        except Exception as e:
            log.warning(f"Failed to compute graph connectivity: {str(e)}")

    # 4. Batch ASW approximation based on nearest within/between-batch distances.
    if "batch_asw" in methods:
        try:
            # Compute average silhouette width per batch
            batch_categories = (
                adata.obs[batch_key].cat.categories
                if hasattr(adata.obs[batch_key], "cat")
                else sorted(adata.obs[batch_key].unique())
            )
            batch_asw_scores = []

            for batch in batch_categories:
                # Get indices for this batch
                batch_mask = adata.obs[batch_key] == batch
                batch_indices = np.where(batch_mask)[0]
                other_indices = np.where(~batch_mask)[0]

                if len(batch_indices) <= 1 or len(other_indices) <= 1:
                    continue

                # Sample for large datasets
                if len(batch_indices) > 1000:
                    batch_indices = np.random.choice(batch_indices, 1000, replace=False)
                if len(other_indices) > 1000:
                    other_indices = np.random.choice(other_indices, 1000, replace=False)

                # Compute distances from batch cells to other batch cells
                from sklearn.metrics import pairwise_distances

                # Within-batch distances
                batch_X = X[batch_indices]
                within_dists = pairwise_distances(batch_X, metric=metric)
                np.fill_diagonal(within_dists, np.inf)  # Exclude self-distances
                avg_within = np.mean(np.min(within_dists, axis=1))

                # Between-batch distances
                other_X = X[other_indices]
                between_dists = pairwise_distances(batch_X, other_X, metric=metric)
                avg_between = np.mean(np.min(between_dists, axis=1))

                # Compute ASW for this batch
                batch_size = len(batch_indices)
                total_cells = adata.n_obs
                weight = batch_size / total_cells

                asw = (avg_between - avg_within) / max(avg_between, avg_within)
                weighted_asw = weight * asw

                batch_asw_scores.append(weighted_asw)

            if batch_asw_scores:
                # Sum weighted ASW scores
                results["batch_asw"] = sum(batch_asw_scores)
                results["batch_asw_model_type"] = "nearest_distance_asw_approximation"
                results["batch_asw_claim_level"] = "batch_separation_diagnostic_heuristic"
                results["batch_asw_review_note"] = (
                    "batch_asw is computed from nearest within/between-batch distances, "
                    "not the standard per-cell silhouette implementation."
                )
                log.info(f"Batch ASW: {results['batch_asw']:.4f}")
            else:
                log.warning("Could not compute batch ASW: insufficient data")

        except Exception as e:
            log.warning(f"Failed to compute batch ASW: {str(e)}")

    results["interpretation"] = _interpret_integration_evaluation(results)
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {}).setdefault("integration", {})[
        "evaluation"
    ] = results

    # Plot
    if plot and results:
        import matplotlib.pyplot as plt

        plot_metrics = [
            k
            for k, value in results.items()
            if k not in ["method", "n_batches", "n_cells", "n_labels", "interpretation"]
            and isinstance(value, (int, float, np.integer, np.floating))
        ]
        if plot_metrics:
            fig, axes = plt.subplots(1, len(plot_metrics), figsize=(4 * len(plot_metrics), 5))
            if len(plot_metrics) == 1:
                axes = [axes]
            for i, m in enumerate(plot_metrics):
                value = results[m]
                axes[i].bar([0], [value], color="skyblue", width=0.5)
                axes[i].set_title(m.replace("_", " ").title())
                axes[i].set_ylim(0, 1)
                axes[i].set_xticks([])
                axes[i].text(
                    0,
                    value / 2,
                    f"{value:.3f}",
                    ha="center",
                    va="center",
                    fontweight="bold",
                )
                axes[i].grid(axis="y", linestyle="--", alpha=0.7)
            fig.suptitle(
                f"Integration Quality: {results.get('method', 'unknown')}",
                fontsize=16,
                fontweight="bold",
            )
            plt.tight_layout(rect=[0, 0.05, 1, 0.95])
            if save_path:
                figure_path = Path(save_path)
                figure_path.parent.mkdir(parents=True, exist_ok=True)
                plt.savefig(figure_path, dpi=300, bbox_inches="tight")
                log.info(f"Saved evaluation plot to {figure_path}")

            _show_or_close(fig)
    return results
