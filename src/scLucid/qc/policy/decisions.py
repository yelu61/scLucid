"""Evidence-based QC decision schema and scoring helpers.

This module turns heterogeneous QC evidence into stable reviewer-facing
columns.  It intentionally separates "mark evidence" from "remove cells":
high mitochondrial fraction, stress, ambient signal, or doublet score can be
biologically meaningful in tumor, CSF, low-RNA immune, and fragile-cell
datasets, so the default output is a decision with confidence and reasons.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Literal

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from ...utils.context import is_tumor_context
from ...utils.helpers import sanitize_for_hdf5
from ..artifacts import record_qc_decision_artifact

QC_DECISION_SCHEMA_VERSION = "qc_decision_schema_v1"
QC_DECISION_VALUES = ("keep", "remove", "review", "sensitivity_only")


GENE_PANELS: dict[str, tuple[str, ...]] = {
    "hemoglobin": ("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2"),
    "platelet": ("PPBP", "PF4", "GP9", "GP1BA", "GP1BB", "ITGA2B", "ITGB3", "TUBB1"),
    "dissociation": (
        "FOS",
        "JUN",
        "JUNB",
        "JUND",
        "EGR1",
        "IER2",
        "IER3",
        "ATF3",
        "HSPA1A",
        "HSPA1B",
        "DNAJB1",
    ),
    "heatshock": ("HSPA1A", "HSPA1B", "HSPA6", "HSPB1", "HSP90AA1", "DNAJB1"),
    "nfkb": ("NFKB1", "NFKBIA", "TNF", "IL1B", "CXCL8", "PTGS2", "ICAM1"),
    "hypoxia": ("HIF1A", "VEGFA", "CA9", "SLC2A1", "PGK1", "LDHA", "BNIP3"),
    "apoptosis": ("BAX", "CASP3", "CASP7", "CASP8", "FAS", "BCL2L11", "BBC3", "CYCS"),
    "inflammatory": ("IL1B", "TNF", "CXCL8", "CXCL10", "NFKBIA", "PTGS2", "ICAM1"),
}


def _matrix_for_scoring(adata: AnnData, layer: str | None = None):
    if layer:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in adata.layers")
        return adata.layers[layer]
    return adata.X


def _present_genes(adata: AnnData, genes: Sequence[str]) -> list[str]:
    var_upper = pd.Index([str(g).upper() for g in adata.var_names])
    present: list[str] = []
    for gene in genes:
        matches = np.flatnonzero(var_upper == gene.upper())
        if matches.size:
            present.append(str(adata.var_names[int(matches[0])]))
    return present


def _row_mean(X) -> np.ndarray:
    if X.shape[1] == 0:
        return np.zeros(X.shape[0], dtype=float)
    if sp.issparse(X):
        return np.asarray(X.mean(axis=1)).ravel().astype(float)
    return np.asarray(X, dtype=float).mean(axis=1)


def _robust_unit_interval(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    out = np.zeros(values.shape[0], dtype=float)
    if not finite.any():
        return out
    lo = float(np.nanpercentile(values[finite], 5))
    hi = float(np.nanpercentile(values[finite], 95))
    if hi <= lo:
        out[finite] = values[finite] > hi
        return out
    out[finite] = np.clip((values[finite] - lo) / (hi - lo), 0.0, 1.0)
    return out


def score_qc_gene_panels(
    adata: AnnData,
    *,
    layer: str | None = None,
    panels: Mapping[str, Sequence[str]] | None = None,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Score contamination/stress panels into stable QC columns.

    Scores are lightweight reviewer signals, not calibrated biological
    probabilities. They are intentionally bounded to 0-1 so they can feed a
    conservative multi-evidence decision engine.
    """
    panel_map = dict(panels or GENE_PANELS)
    X = _matrix_for_scoring(adata, layer)
    var_index = pd.Index(adata.var_names.astype(str))
    rows: list[dict[str, Any]] = []

    for name, genes in panel_map.items():
        score_col = f"{name}_score" if name not in {"dissociation"} else "dissociation_score"
        if score_col in adata.obs and not overwrite:
            matched = [g for g in genes if g in var_index]
            rows.append(
                {
                    "panel": name,
                    "score_column": score_col,
                    "n_genes": len(genes),
                    "n_matched": len(matched),
                    "status": "existing",
                }
            )
            continue
        matched = _present_genes(adata, genes)
        if matched:
            gene_idx = var_index.get_indexer(matched)
            raw = _row_mean(X[:, gene_idx])
            score = _robust_unit_interval(np.log1p(np.maximum(raw, 0.0)))
            status = "scored"
        else:
            score = np.zeros(adata.n_obs, dtype=float)
            status = "no_genes_matched"
        adata.obs[score_col] = score
        rows.append(
            {
                "panel": name,
                "score_column": score_col,
                "n_genes": len(genes),
                "n_matched": len(matched),
                "matched_genes": matched,
                "status": status,
            }
        )

    if "dissociation_score" in adata.obs:
        adata.obs["stress_score"] = np.maximum(
            np.asarray(adata.obs.get("dissociation_score", 0), dtype=float),
            np.asarray(adata.obs.get("heatshock_score", 0), dtype=float),
        )
    elif "heatshock_score" in adata.obs:
        adata.obs["stress_score"] = np.asarray(adata.obs["heatshock_score"], dtype=float)

    qc_ns = adata.uns.setdefault("sclucid", {}).setdefault("qc", {})
    summary = {
        "schema_version": "qc_gene_panel_scores_v1",
        "layer": layer or "X",
        "panels": rows,
        "score_note": (
            "Panel scores are QC review evidence. High stress or hypoxia can be "
            "technical or biological and is not removed automatically."
        ),
    }
    qc_ns["qc_gene_panel_scores"] = sanitize_for_hdf5(summary)
    return summary


def _bool_series(adata: AnnData, col: str) -> pd.Series:
    if col not in adata.obs:
        return pd.Series(False, index=adata.obs_names)
    return adata.obs[col].fillna(False).astype(bool)


def _numeric_series(adata: AnnData, col: str) -> pd.Series:
    if col not in adata.obs:
        return pd.Series(np.nan, index=adata.obs_names, dtype=float)
    return pd.to_numeric(adata.obs[col], errors="coerce")


def _first_numeric_obs(adata: AnnData, candidates: Sequence[str]) -> pd.Series | None:
    """Return the first available numeric obs column from ``candidates``."""
    for col in candidates:
        if col in adata.obs:
            return _numeric_series(adata, col)
    return None


def _ensure_canonical_probability_columns(adata: AnnData) -> dict[str, Any]:
    """Populate stable probability/score columns without inventing evidence.

    Optional tools such as CellBender, EmptyDrops, SoupX/DecontX, or hashing
    pipelines use different column names.  The QC decision schema exposes
    canonical placeholders so reports and downstream code can rely on stable
    columns while still distinguishing unavailable evidence from measured
    evidence via NaN values and provenance metadata.
    """
    provenance: dict[str, Any] = {
        "schema_version": "qc_probability_schema_v1",
        "columns": {},
    }
    mappings = {
        "cell_probability": (
            "cell_probability",
            "cellbender_cell_probability",
            "prob_cell",
            "cell_call_probability",
            "emptydrops_fdr_cell_probability",
        ),
        "empty_droplet_probability": (
            "empty_droplet_probability",
            "empty_probability",
            "background_probability",
            "prob_empty",
        ),
        "doublet_score": (
            "doublet_score",
            "combined_doublet_score",
            "scrublet_score",
            "scdblfinder_score",
            "doubletdetection_score",
        ),
        "ambient_score": (
            "ambient_score",
            "ambient_fraction",
            "contamination_fraction",
            "decontx_contamination",
            "soupx_contamination",
            "cellbender_ambient_fraction",
        ),
        "ambient_fraction": (
            "ambient_fraction",
            "contamination_fraction",
            "decontx_contamination",
            "soupx_contamination",
            "cellbender_ambient_fraction",
        ),
    }
    for canonical, candidates in mappings.items():
        source = next((col for col in candidates if col in adata.obs), None)
        if canonical not in adata.obs:
            values = _first_numeric_obs(adata, candidates)
            if values is None:
                adata.obs[canonical] = np.nan
                source = None
                status = "unavailable"
            else:
                adata.obs[canonical] = values.to_numpy(dtype=float)
                status = "aliased"
        else:
            status = "existing"
            source = canonical
        provenance["columns"][canonical] = {
            "status": status,
            "source": source,
            "available": bool(adata.obs[canonical].notna().any()),
        }
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {})[
        "qc_probability_schema"
    ] = sanitize_for_hdf5(provenance)
    return provenance


def _high_by_quantile(values: pd.Series, *, minimum: float, quantile: float = 0.9) -> pd.Series:
    valid = values.dropna()
    if valid.empty:
        return pd.Series(False, index=values.index)
    threshold = max(float(minimum), float(valid.quantile(quantile)))
    return values.fillna(0.0) >= threshold


def diagnose_stress_sources(
    adata: AnnData,
    *,
    sample_key: str | None = "sampleID",
    cell_type_key: str | None = None,
    stress_score_col: str = "stress_score",
    minimum_cells: int = 10,
    sample_bias_threshold: float = 0.60,
    cell_type_concordance_threshold: float = 0.70,
) -> dict[str, Any]:
    """Determine whether stress signal is technical (sample-driven) or biological.

    Returns a diagnostic summary with flags for:

    - ``sample_driven``: one or a few samples contribute a disproportionate
      fraction of stress-high cells.
    - ``cell_type_concordant``: stress is elevated consistently across many cell
      types, suggesting a shared biological response (e.g., hypoxia, inflammation).
    - ``sample_dominant``: list of samples flagged as drivers.
    - ``dominant_cell_types``: cell types with the highest stress-high fraction.

    The function does not modify ``adata``; it only reads ``.obs`` columns.
    """
    if stress_score_col not in adata.obs.columns:
        return {
            "schema_version": "stress_source_diagnosis_v1",
            "available": False,
            "reason": f"{stress_score_col} not found in adata.obs",
        }

    stress = adata.obs[stress_score_col].fillna(0.0)
    high_mask = stress >= stress.quantile(0.9)
    n_high = int(high_mask.sum())

    result: dict[str, Any] = {
        "schema_version": "stress_source_diagnosis_v1",
        "available": n_high >= minimum_cells,
        "n_high_stress_cells": n_high,
        "sample_driven": False,
        "cell_type_concordant": False,
        "sample_dominant": [],
        "dominant_cell_types": [],
    }

    if n_high < minimum_cells:
        result["reason"] = "too_few_stress_high_cells"
        return result

    # Sample-level bias
    if sample_key is not None and sample_key in adata.obs.columns:
        sample_counts = (
            adata.obs.loc[high_mask, sample_key].astype(str).value_counts(normalize=True)
        )
        if not sample_counts.empty and float(sample_counts.iloc[0]) >= sample_bias_threshold:
            result["sample_driven"] = True
            result["sample_dominant"] = [str(sample_counts.index[0])]
            # Check for additional samples crossing half the bias threshold
            for sample, frac in sample_counts.iloc[1:].items():
                if float(frac) >= sample_bias_threshold * 0.5:
                    result["sample_dominant"].append(str(sample))

    # Cell-type-level concordance
    if cell_type_key is not None and cell_type_key in adata.obs.columns:
        ct_table = adata.obs.groupby(cell_type_key, observed=True).apply(
            lambda g: pd.Series(
                {
                    "n_cells": len(g),
                    "stress_high": int(g[stress_score_col].fillna(0.0) >= g[stress_score_col].quantile(0.9)),
                }
            )
        )
        ct_table["fraction"] = ct_table["stress_high"] / ct_table["n_cells"].clip(lower=1)
        ct_table = ct_table[ct_table["n_cells"] >= minimum_cells].sort_values(
            "fraction", ascending=False
        )
        result["dominant_cell_types"] = [
            {"cell_type": str(idx), "fraction": float(row["fraction"]), "n_cells": int(row["n_cells"])}
            for idx, row in ct_table.head(5).iterrows()
        ]
        # Concordant if >= threshold fraction of cell types (with enough cells) show elevated stress.
        elevated = ct_table[ct_table["fraction"] >= 0.2]
        if len(ct_table) > 0 and (len(elevated) / len(ct_table)) >= cell_type_concordance_threshold:
            result["cell_type_concordant"] = True

    # Interpretation note
    if result["sample_driven"] and not result["cell_type_concordant"]:
        result["interpretation"] = "technical_stress_sample_bias"
        result["recommendation"] = (
            "Stress signal is concentrated in one or a few samples; treat as technical "
            "dissociation/artefact and consider sample-level sensitivity analysis."
        )
    elif result["cell_type_concordant"] and not result["sample_driven"]:
        result["interpretation"] = "biological_stress_response"
        result["recommendation"] = (
            "Stress signal is shared across many cell types; likely reflects a true "
            "biological condition (e.g., hypoxia, inflammation, treatment response)."
        )
    elif result["sample_driven"] and result["cell_type_concordant"]:
        result["interpretation"] = "mixed_stress_pattern"
        result["recommendation"] = (
            "Stress signal has both sample-driven and cross-cell-type components; "
            "review sample quality and biological context before removal."
        )
    else:
        result["interpretation"] = "no_clear_stress_pattern"
        result["recommendation"] = "Stress signal is diffuse; retain cells for review."

    return result


def build_qc_decisions(
    adata: AnnData,
    *,
    tissue_type: str | None = None,
    policy: Literal["conservative", "screening", "strict"] = "conservative",
    score_layer: str | None = None,
    score_panels: bool = True,
    overwrite_scores: bool = False,
    sample_key: str | None = "sampleID",
    cell_type_key: str | None = None,
) -> dict[str, Any]:
    """Build a unified cell-level QC decision table in ``adata.obs``.

    The default conservative policy removes only multi-evidence low-quality
    cells. Single high-MT, stress, ambient, or doublet-like signals are marked
    for review unless they are already part of a strong combined failure.
    """
    if score_panels:
        score_qc_gene_panels(adata, layer=score_layer, overwrite=overwrite_scores)
    _ensure_canonical_probability_columns(adata)

    tumor_context = is_tumor_context(tissue_type)
    obs = adata.obs

    low_counts = _bool_series(adata, "outlier_count") | _bool_series(adata, "outlier_min_counts")
    low_genes = _bool_series(adata, "outlier_min_genes")
    high_mt = _bool_series(adata, "outlier_mt") | _bool_series(adata, "mt_hard_fail")
    low_complexity = _bool_series(adata, "outlier_qc_metrics")
    high_hb = _bool_series(adata, "outlier_hb")
    doublet = _bool_series(adata, "predicted_doublet")

    mt_pct = _numeric_series(adata, "pct_counts_mt")
    if not high_mt.any() and mt_pct.notna().any():
        high_mt = mt_pct >= (35.0 if tumor_context else 25.0)

    stress_high = _high_by_quantile(_numeric_series(adata, "stress_score"), minimum=0.75)
    apoptosis_high = _high_by_quantile(_numeric_series(adata, "apoptosis_score"), minimum=0.75)
    platelet_high = _high_by_quantile(_numeric_series(adata, "platelet_score"), minimum=0.80)
    hb_score_high = _high_by_quantile(_numeric_series(adata, "hemoglobin_score"), minimum=0.80)

    ambient_score = _numeric_series(adata, "ambient_score")
    if "ambient_fraction" in obs:
        ambient_score = _numeric_series(adata, "ambient_fraction")
    ambient_high = _high_by_quantile(ambient_score, minimum=0.30)

    adata.obs["qc_low_counts"] = low_counts.to_numpy(dtype=bool)
    adata.obs["qc_low_genes"] = low_genes.to_numpy(dtype=bool)
    adata.obs["qc_high_mt"] = high_mt.to_numpy(dtype=bool)
    adata.obs["qc_low_complexity"] = low_complexity.to_numpy(dtype=bool)
    adata.obs["qc_high_hb"] = (high_hb | hb_score_high).to_numpy(dtype=bool)
    adata.obs["platelet_contamination"] = platelet_high.to_numpy(dtype=bool)
    adata.obs["hemoglobin_contamination"] = (high_hb | hb_score_high).to_numpy(dtype=bool)
    adata.obs["stress_high"] = stress_high.to_numpy(dtype=bool)
    adata.obs["apoptosis_high"] = apoptosis_high.to_numpy(dtype=bool)
    adata.obs["ambient_risk"] = ambient_high.to_numpy(dtype=bool)

    # Record stress source diagnosis when stress_high is present, so reviewers can
    # distinguish technical dissociation artefacts from biological stress responses.
    stress_diagnosis = diagnose_stress_sources(
        adata,
        sample_key=sample_key if sample_key in adata.obs.columns else None,
        cell_type_key=None,
    )
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {})["stress_source_diagnosis"] = (
        sanitize_for_hdf5(stress_diagnosis)
    )

    evidence_count = (
        low_counts.astype(int)
        + low_genes.astype(int)
        + high_mt.astype(int)
        + low_complexity.astype(int)
        + (high_hb | hb_score_high).astype(int)
        + platelet_high.astype(int)
        + ambient_high.astype(int)
        + apoptosis_high.astype(int)
        + doublet.astype(int)
    )
    # Policy-based minimum evidence for removal. In tumor/stress/CSF/low-RNA
    # contexts we never remove cells based on a single piece of evidence,
    # because high MT, stress scores, or ambient signal can be biological.
    remove_min = 3 if policy == "conservative" else 2 if policy == "screening" else 1
    if tumor_context or tissue_type in ("csf", "snrna", "low_rna", "fragile"):
        remove_min = max(remove_min, 2)
    remove = evidence_count >= remove_min

    biologically_ambiguous = high_mt | stress_high | ambient_high | doublet
    review = biologically_ambiguous & ~remove
    sensitivity = stress_high & ~remove
    if tumor_context:
        # High-MT tumor cells are preserved for review unless they also carry
        # strong multi-evidence failure (>= max(remove_min + 1, 3) issues).
        remove = remove & ~(high_mt & (evidence_count < max(remove_min + 1, 3)))
        review = review | high_mt
        sensitivity = sensitivity | (stress_high & ~remove)

    decision = pd.Series("keep", index=adata.obs_names, dtype=object)
    decision.loc[sensitivity] = "sensitivity_only"
    decision.loc[review] = "review"
    decision.loc[remove] = "remove"

    reasons: list[str] = []
    risk_notes: list[str] = []
    for idx in adata.obs_names:
        cell_reasons: list[str] = []
        if bool(low_counts.loc[idx]):
            cell_reasons.append("low_counts")
        if bool(low_genes.loc[idx]):
            cell_reasons.append("low_genes")
        if bool(high_mt.loc[idx]):
            cell_reasons.append("high_mt")
        if bool(low_complexity.loc[idx]):
            cell_reasons.append("low_complexity")
        if bool((high_hb | hb_score_high).loc[idx]):
            cell_reasons.append("hemoglobin_contamination")
        if bool(platelet_high.loc[idx]):
            cell_reasons.append("platelet_contamination")
        if bool(ambient_high.loc[idx]):
            cell_reasons.append("ambient_risk")
        if bool(stress_high.loc[idx]):
            cell_reasons.append("stress_high")
        if bool(apoptosis_high.loc[idx]):
            cell_reasons.append("apoptosis_high")
        if bool(doublet.loc[idx]):
            cell_reasons.append("doublet")
        reasons.append(";".join(cell_reasons) if cell_reasons else "no_qc_evidence")
        if tumor_context and "high_mt" in cell_reasons:
            risk_notes.append("Tumor context: high MT can reflect malignant/stress biology; review before deletion.")
        elif "stress_high" in cell_reasons and decision.loc[idx] != "remove":
            risk_notes.append("Stress-high cell retained for review/sensitivity analysis.")
        elif "doublet" in cell_reasons and decision.loc[idx] != "remove":
            risk_notes.append("Doublet-like cell marked for review before irreversible removal.")
        else:
            risk_notes.append("")

    confidence = np.clip(evidence_count.astype(float) / max(remove_min, 1), 0.0, 1.0)
    adata.obs["qc_evidence_count"] = evidence_count.astype(int).to_numpy()
    adata.obs["qc_decision"] = pd.Categorical(decision, categories=QC_DECISION_VALUES)
    adata.obs["qc_remove"] = decision.eq("remove").to_numpy(dtype=bool)
    adata.obs["qc_reason"] = reasons
    adata.obs["qc_confidence"] = confidence.to_numpy(dtype=float)
    adata.obs["qc_review_required"] = decision.isin(["review", "sensitivity_only"]).to_numpy()
    adata.obs["qc_biological_risk_note"] = risk_notes
    adata.obs["qc_phase"] = "final_decision"

    summary = summarize_qc_decisions(adata, tissue_type=tissue_type, policy=policy)
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {})["qc_decision_summary"] = (
        sanitize_for_hdf5(summary)
    )
    record_qc_decision_artifact(
        adata,
        summary=summary,
        evidence_columns=[
            "qc_low_counts",
            "qc_low_genes",
            "qc_high_mt",
            "qc_low_complexity",
            "qc_high_hb",
            "platelet_contamination",
            "hemoglobin_contamination",
            "ambient_risk",
            "stress_high",
            "apoptosis_high",
            "predicted_doublet",
        ],
    )
    return summary


def summarize_qc_decisions(
    adata: AnnData,
    *,
    tissue_type: str | None = None,
    policy: str = "conservative",
) -> dict[str, Any]:
    """Summarize unified QC decisions for reports and benchmarks."""
    if "qc_decision" not in adata.obs:
        raise KeyError("adata.obs['qc_decision'] is missing; run build_qc_decisions first.")
    decision_counts = adata.obs["qc_decision"].astype(str).value_counts().to_dict()
    evidence_cols = [
        "qc_low_counts",
        "qc_low_genes",
        "qc_high_mt",
        "qc_low_complexity",
        "qc_high_hb",
        "platelet_contamination",
        "hemoglobin_contamination",
        "ambient_risk",
        "stress_high",
        "apoptosis_high",
        "predicted_doublet",
    ]
    evidence_summary = {
        col: int(adata.obs[col].fillna(False).astype(bool).sum())
        for col in evidence_cols
        if col in adata.obs
    }
    review_required = int(adata.obs.get("qc_review_required", pd.Series(False)).sum())
    return {
        "schema_version": QC_DECISION_SCHEMA_VERSION,
        "policy": policy,
        "tissue_type": tissue_type,
        "tumor_context": bool(is_tumor_context(tissue_type)),
        "n_cells": int(adata.n_obs),
        "decision_counts": {str(k): int(v) for k, v in decision_counts.items()},
        "evidence_summary": evidence_summary,
        "review_required_cells": review_required,
        "remove_cells": int(decision_counts.get("remove", 0)),
        "risk_note": (
            "QC decisions are evidence labels. In tumor/fragile-cell contexts, "
            "review and sensitivity_only cells should be inspected before removal."
        ),
    }


__all__ = [
    "GENE_PANELS",
    "QC_DECISION_SCHEMA_VERSION",
    "QC_DECISION_VALUES",
    "build_qc_decisions",
    "score_qc_gene_panels",
    "summarize_qc_decisions",
]
