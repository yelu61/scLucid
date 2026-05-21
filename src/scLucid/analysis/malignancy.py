"""Malignancy interpretation bridge for the analysis workflow.

This module keeps malignancy calls evidence-first: epithelial/tumor annotation,
tumor-context marker evidence, optional CNV burden, and optional malignancy
signature scores are combined into reviewable calls rather than treated as a
single ground truth label.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from anndata import AnnData

from ..utils import Manager, get_marker_manager, sanitize_for_hdf5

__all__ = ["run_malignancy_interpretation"]

NORMAL_COMPARTMENT_TERMS = (
    "immune",
    "t cell",
    "b cell",
    "nk",
    "myeloid",
    "monocyte",
    "macrophage",
    "dendritic",
    "neutrophil",
    "mast",
    "plasma",
    "endothelial",
    "fibroblast",
    "stromal",
    "pericyte",
    "smooth muscle",
)

TUMOR_COMPARTMENT_TERMS = (
    "tumor",
    "malignant",
    "cancer",
    "carcinoma",
    "epithelial",
    "ductal",
    "squamous",
    "adenocarcinoma",
)


def run_malignancy_interpretation(
    adata: AnnData,
    *,
    annotation_key: str = "cell_type_auto",
    cluster_key: str | None = None,
    species: str = "human",
    cancer_type: str | None = None,
    marker_manager: Manager | None = None,
    run_cnv: bool = False,
    cnv_score_key: str | None = None,
    cnv_key_added: str = "analysis_cnv",
    reference_labels: Sequence[str] | None = None,
    run_malignancy_score: bool = True,
    malignancy_score_key: str = "analysis_malignancy",
    key_added: str = "malignancy_call",
    score_key: str = "malignancy_interpretation_score",
    threshold: float = 0.55,
    suspect_threshold: float = 0.35,
) -> pd.DataFrame:
    """Interpret malignant-cell evidence and store reviewable calls.

    Parameters are intentionally lightweight so the function can be used inside
    ``run_standard_analysis`` without requiring heavy CNV backends. If
    ``run_cnv=True`` and no CNV score is present, scLucid's expression-based CNV
    estimator is executed and consumed as one evidence source.
    """
    if annotation_key not in adata.obs.columns:
        raise KeyError(f"'{annotation_key}' not found in adata.obs.")
    if cluster_key is not None and cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs.")

    annotation_labels = adata.obs[annotation_key].astype(str)
    cnv_key = _resolve_or_run_cnv(
        adata,
        run_cnv=run_cnv,
        cnv_score_key=cnv_score_key,
        cnv_key_added=cnv_key_added,
        annotation_key=annotation_key,
        reference_labels=reference_labels,
    )
    mal_key = _resolve_or_run_malignancy_score(
        adata,
        run_malignancy_score=run_malignancy_score,
        malignancy_score_key=malignancy_score_key,
        cancer_type=cancer_type,
    )
    tumor_marker_score = _score_tumor_marker_context(
        adata,
        marker_manager=marker_manager,
        species=species,
        cancer_type=cancer_type,
    )
    annotation_prior = annotation_labels.map(_annotation_malignancy_prior).astype(float)
    normal_prior = annotation_labels.map(_annotation_normal_prior).astype(float)

    evidence_parts: list[tuple[str, pd.Series, float]] = [
        ("annotation_prior", annotation_prior, 0.15),
        ("tumor_marker_score", tumor_marker_score, 0.20),
    ]
    if cnv_key is not None:
        evidence_parts.append(("cnv_score", _normalize_series(adata.obs[cnv_key]), 0.40))
    if mal_key is not None:
        evidence_parts.append(
            ("malignancy_signature_score", _normalize_series(adata.obs[mal_key]), 0.25)
        )

    total_weight = sum(weight for _, _, weight in evidence_parts) or 1.0
    combined = sum(values * weight for _, values, weight in evidence_parts) / total_weight
    combined = combined.clip(lower=0.0, upper=1.0)

    calls = pd.Series("unresolved", index=adata.obs_names, dtype=object)
    calls.loc[combined >= threshold] = "malignant"
    calls.loc[(combined >= suspect_threshold) & (combined < threshold)] = "suspect_malignant"
    calls.loc[(combined < suspect_threshold) & (normal_prior >= 0.5)] = "non_malignant"
    calls.loc[(combined < suspect_threshold) & (annotation_prior < 0.2)] = "non_malignant"

    confidence = _call_confidence(combined, calls, threshold=threshold, suspect_threshold=suspect_threshold)

    adata.obs[score_key] = combined.astype(float)
    adata.obs[key_added] = pd.Categorical(
        calls,
        categories=["malignant", "suspect_malignant", "non_malignant", "unresolved"],
    )
    adata.obs[f"{key_added}_confidence"] = confidence.astype(float)
    adata.obs[f"{key_added}_basis"] = pd.Categorical(
        _build_basis_labels(
            annotation_prior=annotation_prior,
            normal_prior=normal_prior,
            cnv_score=_normalize_series(adata.obs[cnv_key]) if cnv_key else None,
            tumor_marker_score=tumor_marker_score,
            malignancy_score=_normalize_series(adata.obs[mal_key]) if mal_key else None,
        )
    )

    summary_df = _build_malignancy_review_table(
        adata,
        annotation_key=annotation_key,
        cluster_key=cluster_key,
        key_added=key_added,
        score_key=score_key,
    )
    interpretation_summary = {
        "available": True,
        "annotation_key": annotation_key,
        "cluster_key": cluster_key,
        "call_key": key_added,
        "score_key": score_key,
        "cnv_score_key": cnv_key,
        "malignancy_score_key": mal_key,
        "threshold": float(threshold),
        "suspect_threshold": float(suspect_threshold),
        "n_malignant": int((calls == "malignant").sum()),
        "n_suspect_malignant": int((calls == "suspect_malignant").sum()),
        "n_non_malignant": int((calls == "non_malignant").sum()),
        "n_unresolved": int((calls == "unresolved").sum()),
        "mean_score": float(combined.mean()) if combined.notna().any() else None,
        "review_required": bool((calls == "suspect_malignant").any() or (calls == "unresolved").any()),
        "evidence_sources": [name for name, _, _ in evidence_parts],
    }
    analysis_ns = adata.uns.setdefault("sclucid", {}).setdefault("analysis", {})
    malignancy_ns = analysis_ns.setdefault("malignancy", {})
    malignancy_ns["malignancy_interpretation_table"] = summary_df
    malignancy_ns["malignancy_interpretation_summary"] = sanitize_for_hdf5(
        interpretation_summary
    )
    malignancy_ns["malignancy_interpretation_params"] = sanitize_for_hdf5(
        {
            "species": species,
            "cancer_type": cancer_type,
            "run_cnv": bool(run_cnv),
            "reference_labels": list(reference_labels or []),
            "run_malignancy_score": bool(run_malignancy_score),
        }
    )
    return summary_df


def _resolve_or_run_cnv(
    adata: AnnData,
    *,
    run_cnv: bool,
    cnv_score_key: str | None,
    cnv_key_added: str,
    annotation_key: str,
    reference_labels: Sequence[str] | None,
) -> str | None:
    if cnv_score_key and cnv_score_key in adata.obs.columns:
        return cnv_score_key
    for candidate in ("cnv_score", f"{cnv_key_added}_score"):
        if candidate in adata.obs.columns:
            return candidate
    if not run_cnv:
        return None

    from ..tumor.cnv.infercnv import infer_cnv

    ref_labels = list(reference_labels or _infer_reference_labels(adata.obs[annotation_key]))
    available_labels = set(adata.obs[annotation_key].astype(str))
    ref_labels = [label for label in ref_labels if str(label) in available_labels]
    infer_cnv(
        adata,
        reference_cells=ref_labels or None,
        reference_key=annotation_key,
        key_added=cnv_key_added,
    )
    resolved = f"{cnv_key_added}_score"
    return resolved if resolved in adata.obs.columns else None


def _resolve_or_run_malignancy_score(
    adata: AnnData,
    *,
    run_malignancy_score: bool,
    malignancy_score_key: str,
    cancer_type: str | None,
) -> str | None:
    if malignancy_score_key in adata.obs.columns:
        return malignancy_score_key
    if "malignancy" in adata.obs.columns:
        return "malignancy"
    if not run_malignancy_score:
        return None
    try:
        from ..tumor.malignancy.scoring import score_malignancy

        score_malignancy(adata, key_added=malignancy_score_key, cancer_type=cancer_type)
    except Exception:
        return None
    return malignancy_score_key if malignancy_score_key in adata.obs.columns else None


def _score_tumor_marker_context(
    adata: AnnData,
    *,
    marker_manager: Manager | None,
    species: str,
    cancer_type: str | None,
) -> pd.Series:
    try:
        mgr = marker_manager or get_marker_manager(
            species=species,
            cancer_type=cancer_type,
            view="tumor_interpretation",
        )
    except Exception:
        return pd.Series(0.0, index=adata.obs_names)

    markers: list[str] = []
    for cell in getattr(mgr, "CELLS", {}).values():
        metadata = getattr(cell, "metadata", {}) or {}
        if metadata.get("use_for_malignancy_interpretation") is False:
            continue
        markers.extend(str(g) for g in getattr(cell, "markers", []) or [])
    markers = list(dict.fromkeys(g for g in markers if g in adata.var_names))
    if not markers:
        return pd.Series(0.0, index=adata.obs_names)

    X = adata[:, markers].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    values = np.asarray(X).mean(axis=1)
    return _normalize_series(pd.Series(values, index=adata.obs_names))


def _annotation_malignancy_prior(label: str) -> float:
    text = str(label).lower()
    if any(term in text for term in ("immune", "endothelial", "fibroblast", "stromal")):
        return 0.0
    if any(term in text for term in TUMOR_COMPARTMENT_TERMS):
        return 0.65
    return 0.1


def _annotation_normal_prior(label: str) -> float:
    text = str(label).lower()
    return 1.0 if any(term in text for term in NORMAL_COMPARTMENT_TERMS) else 0.0


def _infer_reference_labels(labels: pd.Series) -> list[str]:
    unique = labels.astype(str).dropna().unique().tolist()
    return [
        label
        for label in unique
        if _annotation_normal_prior(label) >= 0.5 and _annotation_malignancy_prior(label) == 0.0
    ]


def _normalize_series(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    if numeric.notna().sum() == 0:
        return pd.Series(0.0, index=values.index)
    min_value = float(numeric.min())
    max_value = float(numeric.max())
    if np.isclose(max_value, min_value):
        return pd.Series(0.0, index=values.index)
    return ((numeric - min_value) / (max_value - min_value)).fillna(0.0)


def _call_confidence(
    scores: pd.Series,
    calls: pd.Series,
    *,
    threshold: float,
    suspect_threshold: float,
) -> pd.Series:
    confidence = pd.Series(0.5, index=scores.index, dtype=float)
    malignant = calls == "malignant"
    suspect = calls == "suspect_malignant"
    non_malignant = calls == "non_malignant"
    confidence.loc[malignant] = (scores.loc[malignant] - threshold) / max(1e-6, 1 - threshold)
    confidence.loc[suspect] = 1.0 - (
        (scores.loc[suspect] - suspect_threshold).abs()
        / max(1e-6, threshold - suspect_threshold)
    )
    confidence.loc[non_malignant] = 1.0 - (
        scores.loc[non_malignant] / max(1e-6, suspect_threshold)
    )
    return confidence.clip(lower=0.0, upper=1.0)


def _build_basis_labels(
    *,
    annotation_prior: pd.Series,
    normal_prior: pd.Series,
    cnv_score: pd.Series | None,
    tumor_marker_score: pd.Series,
    malignancy_score: pd.Series | None,
) -> pd.Series:
    labels = []
    for cell in annotation_prior.index:
        evidence = []
        if normal_prior.loc[cell] >= 0.5:
            evidence.append("normal_annotation")
        if annotation_prior.loc[cell] >= 0.5:
            evidence.append("tumor_annotation")
        if cnv_score is not None and cnv_score.loc[cell] >= 0.5:
            evidence.append("cnv_high")
        if tumor_marker_score.loc[cell] >= 0.5:
            evidence.append("tumor_marker_high")
        if malignancy_score is not None and malignancy_score.loc[cell] >= 0.5:
            evidence.append("malignancy_signature_high")
        labels.append("+".join(evidence) if evidence else "weak_evidence")
    return pd.Series(labels, index=annotation_prior.index)


def _build_malignancy_review_table(
    adata: AnnData,
    *,
    annotation_key: str,
    cluster_key: str | None,
    key_added: str,
    score_key: str,
) -> pd.DataFrame:
    group_key = cluster_key if cluster_key and cluster_key in adata.obs.columns else annotation_key
    rows = []
    for group, obs in adata.obs.groupby(group_key, observed=False):
        calls = obs[key_added].astype(str)
        labels = obs[annotation_key].astype(str)
        rows.append(
            {
                "group": str(group),
                "n_cells": int(obs.shape[0]),
                "dominant_annotation": labels.value_counts().index[0] if not labels.empty else "",
                "mean_score": float(pd.to_numeric(obs[score_key], errors="coerce").mean()),
                "malignant_fraction": float((calls == "malignant").mean()),
                "suspect_fraction": float((calls == "suspect_malignant").mean()),
                "non_malignant_fraction": float((calls == "non_malignant").mean()),
                "unresolved_fraction": float((calls == "unresolved").mean()),
                "dominant_call": calls.value_counts().index[0] if not calls.empty else "unresolved",
                "needs_review": bool(
                    ((calls == "suspect_malignant") | (calls == "unresolved")).any()
                ),
            }
        )
    return pd.DataFrame(rows)
