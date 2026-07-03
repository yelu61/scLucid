"""Iterative QC quick-review helpers."""

from __future__ import annotations

import logging
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from anndata import AnnData

log = logging.getLogger(__name__)


def _as_bool_obs(adata: AnnData, column: str) -> pd.Series:
    """Return a boolean obs column or an all-False series."""
    if column not in adata.obs:
        return pd.Series(False, index=adata.obs_names)
    return adata.obs[column].fillna(False).astype(bool)


def _as_float_obs(adata: AnnData, column: str) -> pd.Series:
    """Return a numeric obs column or an all-NaN series."""
    if column not in adata.obs:
        return pd.Series(np.nan, index=adata.obs_names, dtype=float)
    return pd.to_numeric(adata.obs[column], errors="coerce")


def _sample_for_quick_review(
    adata: AnnData,
    *,
    max_cells: int,
    random_state: int,
) -> AnnData:
    """Return a deterministic subset for quick biology review."""
    if max_cells <= 0 or adata.n_obs <= max_cells:
        return adata.copy()
    rng = np.random.default_rng(random_state)
    keep = np.sort(rng.choice(adata.n_obs, size=max_cells, replace=False))
    return adata[keep].copy()


def _run_quick_biology_review(
    adata: AnnData,
    *,
    sample_key: str,
    max_cells: int = 2000,
    n_top_genes: int = 1000,
    n_pcs: int = 20,
    n_neighbors: int = 10,
    resolution: float = 0.5,
    random_state: int = 0,
) -> Dict[str, Any]:
    """Run a temporary quick embedding/cluster review for iterative QC."""
    if adata.n_obs < 20 or adata.n_vars < 20:
        return {
            "schema_version": "quick_biology_review_v1",
            "status": "skipped_too_few_cells",
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "review_required": False,
            "review_findings": [],
            "cluster_qc_table": [],
        }

    try:
        tmp = _sample_for_quick_review(adata, max_cells=max_cells, random_state=random_state)
        tmp.X = tmp.X.copy()

        # Avoid double-normalising if the input matrix is already log1p-normalised.
        # A raw count matrix typically has a max value well above 100; log1p(CPM)
        # or log1p-normalised data is usually below 10-20.
        x_max = float(np.asarray(tmp.X).max()) if not scipy.sparse.issparse(tmp.X) else float(tmp.X.max())
        if x_max > 100:
            sc.pp.normalize_total(tmp, target_sum=1e4)
            sc.pp.log1p(tmp)
        else:
            log.debug(
                "Quick biology review: input max value %.2f suggests already-normalised data; "
                "skipping normalize_total/log1p.",
                x_max,
            )
        if tmp.n_vars > 50:
            sc.pp.highly_variable_genes(
                tmp,
                n_top_genes=min(n_top_genes, tmp.n_vars),
                flavor="seurat",
            )
            if "highly_variable" in tmp.var and bool(tmp.var["highly_variable"].any()):
                tmp = tmp[:, tmp.var["highly_variable"].to_numpy()].copy()
        sc.pp.scale(tmp, max_value=10)
        effective_pcs = min(n_pcs, tmp.n_obs - 1, tmp.n_vars - 1)
        if effective_pcs < 2:
            raise ValueError("Too few cells/genes for PCA in quick biology review.")
        sc.tl.pca(tmp, n_comps=effective_pcs, random_state=random_state)
        sc.pp.neighbors(
            tmp,
            n_neighbors=min(n_neighbors, max(2, tmp.n_obs - 1)),
            n_pcs=effective_pcs,
            random_state=random_state,
        )
        sc.tl.umap(tmp, random_state=random_state)
        sc.tl.leiden(
            tmp,
            resolution=resolution,
            key_added="quick_qc_leiden",
            random_state=random_state,
        )

        obs = tmp.obs
        cluster_rows: list[dict[str, Any]] = []
        review_findings: list[dict[str, Any]] = []
        for cluster, group in obs.groupby("quick_qc_leiden", observed=True):
            n_cells = int(group.shape[0])
            qc_remove_frac = float(_as_bool_obs(tmp, "qc_remove").loc[group.index].mean())
            review_frac = float(_as_bool_obs(tmp, "qc_review_required").loc[group.index].mean())
            doublet_frac = float(_as_bool_obs(tmp, "predicted_doublet").loc[group.index].mean())
            stress_frac = float(_as_bool_obs(tmp, "stress_high").loc[group.index].mean())
            ambient_frac = float(_as_bool_obs(tmp, "ambient_risk").loc[group.index].mean())
            mean_mt = _as_float_obs(tmp, "pct_counts_mt").loc[group.index].mean()
            dominant_sample_fraction = None
            dominant_sample = None
            if sample_key in group:
                sample_counts = group[sample_key].astype(str).value_counts(normalize=True)
                if not sample_counts.empty:
                    dominant_sample = str(sample_counts.index[0])
                    dominant_sample_fraction = float(sample_counts.iloc[0])
            flags: list[str] = []
            if qc_remove_frac >= 0.30 or review_frac >= 0.50:
                flags.append("qc_enriched_cluster")
            if stress_frac >= 0.50:
                flags.append("stress_enriched_cluster")
            if doublet_frac >= 0.20:
                flags.append("doublet_enriched_cluster")
            if ambient_frac >= 0.30:
                flags.append("ambient_enriched_cluster")
            if dominant_sample_fraction is not None and dominant_sample_fraction >= 0.80 and n_cells >= 10:
                flags.append("sample_dominated_cluster")

            row = {
                "cluster": str(cluster),
                "n_cells": n_cells,
                "qc_remove_fraction": qc_remove_frac,
                "qc_review_required_fraction": review_frac,
                "doublet_fraction": doublet_frac,
                "stress_high_fraction": stress_frac,
                "ambient_risk_fraction": ambient_frac,
                "mean_pct_counts_mt": None if pd.isna(mean_mt) else float(mean_mt),
                "dominant_sample": dominant_sample,
                "dominant_sample_fraction": dominant_sample_fraction,
                "review_flags": flags,
            }
            cluster_rows.append(row)
            for flag in flags:
                review_findings.append(
                    {
                        "scope": "cluster",
                        "cluster": str(cluster),
                        "finding": flag,
                        "review_required": True,
                    }
                )

        stress_mask = _as_bool_obs(tmp, "stress_high")
        if sample_key in obs and int(stress_mask.sum()) >= 10:
            stress_samples = obs.loc[stress_mask, sample_key].astype(str).value_counts(normalize=True)
            if not stress_samples.empty and float(stress_samples.iloc[0]) >= 0.60:
                review_findings.append(
                    {
                        "scope": "sample",
                        "finding": "stress_high_sample_bias",
                        "sample": str(stress_samples.index[0]),
                        "fraction_of_stress_high_cells": float(stress_samples.iloc[0]),
                        "review_required": True,
                    }
                )

        doublet_mask = _as_bool_obs(tmp, "predicted_doublet")
        if int(doublet_mask.sum()) >= 5:
            doublet_clusters = obs.loc[doublet_mask, "quick_qc_leiden"].astype(str).nunique()
            review_findings.append(
                {
                    "scope": "embedding",
                    "finding": "doublet_boundary_review",
                    "doublet_positive_clusters": int(doublet_clusters),
                    "review_required": bool(doublet_clusters > 1),
                }
            )

        ambient_mask = _as_bool_obs(tmp, "ambient_risk")
        if int(ambient_mask.sum()) >= 5:
            ambient_clusters = obs.loc[ambient_mask, "quick_qc_leiden"].astype(str).nunique()
            if ambient_clusters >= max(2, obs["quick_qc_leiden"].nunique() // 2):
                review_findings.append(
                    {
                        "scope": "embedding",
                        "finding": "ambient_marker_widespread_leakage",
                        "ambient_positive_clusters": int(ambient_clusters),
                        "review_required": True,
                    }
                )

        return {
            "schema_version": "quick_biology_review_v1",
            "status": "complete",
            "n_cells_reviewed": int(tmp.n_obs),
            "n_genes_reviewed": int(tmp.n_vars),
            "cluster_key": "quick_qc_leiden",
            "n_clusters": int(obs["quick_qc_leiden"].nunique()),
            "umap_computed": bool("X_umap" in tmp.obsm),
            "review_required": any(item.get("review_required") for item in review_findings),
            "review_findings": review_findings,
            "cluster_qc_table": cluster_rows,
            "note": (
                "Quick review was computed on a temporary normalized/log1p subset and "
                "does not modify formal preprocessing layers."
            ),
        }
    except Exception as exc:
        return {
            "schema_version": "quick_biology_review_v1",
            "status": "failed",
            "error": str(exc),
            "review_required": True,
            "review_findings": [
                {
                    "scope": "workflow",
                    "finding": "quick_biology_review_failed",
                    "review_required": True,
                }
            ],
            "cluster_qc_table": [],
        }


def _refine_qc_decisions_from_review(
    adata: AnnData,
    quick_biology_review: Dict[str, Any],
    sample_key: str,
) -> Dict[str, Any]:
    """Use quick-biology-review findings to refine qc_decision labels.

    The quick review runs on a temporary embedding and may flag:

    - ``qc_enriched_cluster``: clusters with >30% qc_remove or >50% review_required
    - ``stress_enriched_cluster`` / ``stress_high_sample_bias``: stress signals
    - ``doublet_enriched_cluster`` / ``doublet_boundary_review``: doublet signals
    - ``ambient_enriched_cluster`` / ``ambient_marker_widespread_leakage``: ambient leakage
    - ``sample_dominated_cluster``: technical batch/sample segregation

    This function does **not** delete cells. It upgrades the decision label of
    affected cells to ``review`` or ``sensitivity_only`` when the quick review
    suggests the original decision may be missing biology or over-filtering.
    """
    if quick_biology_review.get("status") != "complete":
        return {"refined": False, "reason": "review_not_complete", "changes": {}}

    findings = quick_biology_review.get("review_findings", [])
    if not findings:
        return {"refined": False, "reason": "no_findings", "changes": {}}

    original_decisions = adata.obs["qc_decision"].copy()
    changes: Dict[str, int] = {}

    # Identify affected samples flagged by review findings.
    affected_samples: set[str] = set()
    for finding in findings:
        if finding.get("finding") == "stress_high_sample_bias":
            affected_samples.add(str(finding.get("sample", "")))

    if affected_samples and sample_key in adata.obs:
        sample_mask = adata.obs[sample_key].astype(str).isin(affected_samples)
        # Stress-high cells in affected samples: keep under sensitivity review.
        stress_mask = _as_bool_obs(adata, "stress_high") & sample_mask
        upgrade_mask = stress_mask & adata.obs["qc_decision"].isin(["keep"])
        if upgrade_mask.any():
            adata.obs.loc[upgrade_mask, "qc_decision"] = "sensitivity_only"
            adata.obs.loc[upgrade_mask, "qc_review_required"] = True
            changes["stress_sample_to_sensitivity"] = int(upgrade_mask.sum())

    # Clusters enriched for QC issues: ensure cells are at least under review.
    cluster_table = quick_biology_review.get("cluster_qc_table", [])
    review_clusters: set[str] = set()
    for row in cluster_table:
        flags = row.get("review_flags", [])
        if any(f in flags for f in ("qc_enriched_cluster", "ambient_enriched_cluster")):
            review_clusters.add(str(row.get("cluster", "")))

    if review_clusters:
        # Map quick-review clusters back to full adata via the sampled subset.
        # The quick review stores its cluster assignments only in the temporary
        # object, so we approximate by upgrading cells that share the flagged
        # QC evidence across the whole dataset.
        for flag in ("qc_high_mt", "ambient_risk", "qc_low_complexity"):
            if flag in adata.obs.columns:
                upgrade_mask = (
                    adata.obs[flag].astype(bool)
                    & adata.obs["qc_decision"].isin(["keep"])
                )
                if upgrade_mask.any():
                    adata.obs.loc[upgrade_mask, "qc_decision"] = "review"
                    adata.obs.loc[upgrade_mask, "qc_review_required"] = True
                    changes.setdefault(f"{flag}_to_review", 0)
                    changes[f"{flag}_to_review"] += int(upgrade_mask.sum())

    # Doublet-enriched clusters: doublet-like cells already flagged by the
    # algorithm should remain under review if they were about to be kept.
    doublet_upgrade = (
        _as_bool_obs(adata, "predicted_doublet")
        & adata.obs["qc_decision"].isin(["keep"])
    )
    if doublet_upgrade.any():
        adata.obs.loc[doublet_upgrade, "qc_decision"] = "review"
        adata.obs.loc[doublet_upgrade, "qc_review_required"] = True
        changes["doublet_to_review"] = int(doublet_upgrade.sum())

    # Recompute qc_remove to respect refined decisions.
    adata.obs["qc_remove"] = adata.obs["qc_decision"].eq("remove").to_numpy(bool)

    n_changed = int((adata.obs["qc_decision"] != original_decisions).sum())
    return {
        "refined": n_changed > 0,
        "n_changed": n_changed,
        "changes": changes,
    }


def _build_iterative_qc_summary(
    adata: AnnData,
    *,
    final_filter_policy: str,
    run_quick_review: bool,
    quick_biology_review: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Build the phase-level audit summary for iterative QC."""
    qc_ns = adata.uns.get("sclucid", {}).get("qc", {})
    decision_summary = qc_ns.get("qc_decision_summary", {})
    filtering_summary = qc_ns.get("filtering_results", {})
    phases = [
        {
            "phase": "lenient_cell_screen",
            "status": "complete",
            "outputs": [
                "total_counts",
                "n_genes_by_counts",
                "pct_counts_mt",
                "outlier_*",
            ],
            "note": "Extreme barcode/cell quality evidence was marked before final filtering.",
        },
        {
            "phase": "doublet_contamination_stress_evidence",
            "status": "complete",
            "outputs": [
                "predicted_doublet",
                "qc_decision",
                "qc_reason",
                "qc_confidence",
                "stress_score",
                "hemoglobin_score",
                "platelet_score",
            ],
            "note": "Evidence columns are review signals; ambiguous biology is not automatically deleted.",
        },
        {
            "phase": "quick_biology_review",
            "status": "not_run"
            if not run_quick_review
            else (quick_biology_review or {}).get("status", "not_run"),
            "outputs": ["quick_biology_review"] if run_quick_review else [],
            "note": (
                "Quick embedding/cluster review was not requested."
                if not run_quick_review
                else "Temporary embedding review summarized QC/stress/doublet/ambient patterns."
            ),
        },
        {
            "phase": "final_qc_decision",
            "status": "complete",
            "outputs": ["qc_decision_summary", "filtering_results"],
            "note": f"Final filter policy: {final_filter_policy}.",
        },
    ]
    return {
        "schema_version": "iterative_qc_summary_v1",
        "final_filter_policy": final_filter_policy,
        "run_quick_review": bool(run_quick_review),
        "phases": phases,
        "quick_biology_review": quick_biology_review or {},
        "qc_decision_summary": decision_summary,
        "filtering_summary": filtering_summary,
        "recommended_next_step": (
            "Run formal preprocessing on retained cells; inspect review/sensitivity_only cells "
            "before irreversible exclusion in tumor or fragile-cell contexts."
        ),
    }
