"""Evidence-calibrated, read-only QC review and explicit policy execution."""

from __future__ import annotations

import hashlib
import uuid
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sparse
from anndata import AnnData
from sklearn.mixture import GaussianMixture

from ...decision import DecisionCard, QCPolicy, RunEvidence
from ...utils.context import AnalysisContext, infer_analysis_context
from ...utils.sanitize import sanitize_for_hdf5

_PROFILE_RULES = {
    "scrna": {
        "catastrophic_min_median_genes": 200.0,
        "catastrophic_max_median_mt": 50.0,
        "catastrophic_max_median_top20": 90.0,
        "global_min_genes": 200.0,
        "global_max_mt": 30.0,
        "global_max_top20": 80.0,
    },
    "snrna": {
        "catastrophic_min_median_genes": 100.0,
        "catastrophic_max_median_mt": 20.0,
        "catastrophic_max_median_top20": 95.0,
        "global_min_genes": 100.0,
        "global_max_mt": 10.0,
        "global_max_top20": 90.0,
    },
}


def _profile_name(context: AnalysisContext) -> str:
    assay = str(context.assay or "scrna").lower().replace("-", "")
    return "snrna" if "nuc" in assay or assay.startswith("sn") else "scrna"


def _fingerprint(adata: AnnData) -> dict[str, Any]:
    digest = hashlib.sha256()
    for name in adata.obs_names.astype(str):
        digest.update(name.encode("utf-8"))
        digest.update(b"\0")
    matrix = _matrix(adata)
    matrix_digest = hashlib.sha256()
    if sparse.issparse(matrix):
        csr = matrix.tocsr()
        matrix_digest.update(np.ascontiguousarray(csr.data).tobytes())
        matrix_digest.update(np.ascontiguousarray(csr.indices).tobytes())
        matrix_digest.update(np.ascontiguousarray(csr.indptr).tobytes())
    else:
        matrix_digest.update(np.ascontiguousarray(matrix).tobytes())
    return {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "obs_names_sha256": digest.hexdigest(),
        "matrix_source": "layers[counts]" if "counts" in adata.layers else "X",
        "matrix_sha256": matrix_digest.hexdigest(),
    }


def _matrix(adata: AnnData):
    return adata.layers.get("counts", adata.X)


def _row_sum(matrix) -> np.ndarray:
    return np.asarray(matrix.sum(axis=1)).ravel().astype(float)


def _detected_genes(matrix) -> np.ndarray:
    return np.asarray((matrix > 0).sum(axis=1)).ravel().astype(float)


def _top_fraction(matrix, n_top: int = 20) -> np.ndarray:
    n_top = min(int(n_top), int(matrix.shape[1]))
    totals = _row_sum(matrix)
    if sparse.issparse(matrix):
        csr = matrix.tocsr()
        top = np.zeros(csr.shape[0], dtype=float)
        for idx in range(csr.shape[0]):
            values = csr.data[csr.indptr[idx] : csr.indptr[idx + 1]]
            if values.size <= n_top:
                top[idx] = values.sum()
            elif values.size:
                top[idx] = np.partition(values, -n_top)[-n_top:].sum()
    else:
        values = np.asarray(matrix)
        if n_top == values.shape[1]:
            top = values.sum(axis=1).astype(float)
        else:
            top = np.partition(values, -n_top, axis=1)[:, -n_top:].sum(axis=1).astype(float)
    return np.divide(100.0 * top, totals, out=np.zeros_like(top), where=totals > 0)


def _mt_mask(adata: AnnData, species: str) -> np.ndarray:
    if "mt" in adata.var:
        return np.asarray(adata.var["mt"].fillna(False), dtype=bool)
    names = pd.Index(adata.var_names.astype(str))
    if str(species).lower().startswith("mouse"):
        return np.asarray(names.str.match(r"^mt-"), dtype=bool)
    return np.asarray(names.str.match(r"^MT-"), dtype=bool)


def _metric_frame(adata: AnnData, context: AnalysisContext) -> tuple[pd.DataFrame, dict[str, Any]]:
    matrix = _matrix(adata)
    total = (
        adata.obs["total_counts"].to_numpy(float)
        if "total_counts" in adata.obs
        else _row_sum(matrix)
    )
    genes = (
        adata.obs["n_genes_by_counts"].to_numpy(float)
        if "n_genes_by_counts" in adata.obs
        else _detected_genes(matrix)
    )
    mt = None
    mt_source = "unavailable"
    if "pct_counts_mt" in adata.obs:
        mt = adata.obs["pct_counts_mt"].to_numpy(float)
        mt_source = "obs[pct_counts_mt]"
    else:
        mt_mask = _mt_mask(adata, context.species)
        if bool(mt_mask.any()):
            mt_counts = _row_sum(matrix[:, mt_mask])
            mt = np.divide(100.0 * mt_counts, total, out=np.zeros_like(total), where=total > 0)
            mt_source = "species_aware_gene_mask"
    top_col = "pct_counts_in_top_20_genes"
    top20 = (
        adata.obs[top_col].to_numpy(float) if top_col in adata.obs else _top_fraction(matrix, 20)
    )
    frame = pd.DataFrame(
        {
            "total_counts": total,
            "n_genes_by_counts": genes,
            "pct_counts_mt": np.nan if mt is None else mt,
            "pct_counts_in_top_20_genes": top20,
        },
        index=adata.obs_names.astype(str),
    )
    return frame, {
        "matrix_source": "layers[counts]" if "counts" in adata.layers else "X",
        "mt_source": mt_source,
        "mt_available": mt is not None,
    }


def _mad(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return 0.0
    return float(np.median(np.abs(finite - np.median(finite))))


def _robust_z(values: np.ndarray, *, direction: str) -> np.ndarray:
    finite = np.isfinite(values)
    result = np.zeros(values.shape, dtype=float)
    if not finite.any():
        return result
    center = float(np.median(values[finite]))
    scale = max(1.4826 * _mad(values[finite]), 1e-8)
    raw = (values - center) / scale
    result[finite] = -raw[finite] if direction == "lower" else raw[finite]
    return result


def _sample_review(
    metrics: pd.DataFrame,
    sample_values: pd.Series,
    rules: dict[str, float],
) -> list[dict[str, Any]]:
    table = (
        metrics.assign(_sample=sample_values.astype(str).to_numpy())
        .groupby("_sample", observed=True)
        .median(numeric_only=True)
    )
    relative = {
        "n_genes_by_counts": _robust_z(table["n_genes_by_counts"].to_numpy(), direction="lower"),
        "pct_counts_mt": _robust_z(table["pct_counts_mt"].to_numpy(), direction="upper"),
        "pct_counts_in_top_20_genes": _robust_z(
            table["pct_counts_in_top_20_genes"].to_numpy(), direction="upper"
        ),
    }
    rows: list[dict[str, Any]] = []
    for idx, sample in enumerate(table.index.astype(str)):
        values = table.iloc[idx]
        absolute_axes = {
            "very_low_median_genes": bool(
                values["n_genes_by_counts"] < rules["catastrophic_min_median_genes"]
            ),
            "very_high_median_mt": bool(
                np.isfinite(values["pct_counts_mt"])
                and values["pct_counts_mt"] > rules["catastrophic_max_median_mt"]
            ),
            "extreme_top20_dominance": bool(
                values["pct_counts_in_top_20_genes"] > rules["catastrophic_max_median_top20"]
            ),
        }
        relative_axes = {
            key: bool(len(table) >= 3 and scores[idx] >= 4.0) for key, scores in relative.items()
        }
        n_absolute = sum(absolute_axes.values())
        n_relative = sum(relative_axes.values())
        if n_absolute >= 2:
            status = "BLOCKED"
        elif n_absolute >= 1 or n_relative >= 1:
            status = "REVIEW"
        else:
            status = "READY"
        rows.append(
            {
                "sample": sample,
                "status": status,
                "n_cells": int((sample_values.astype(str) == sample).sum()),
                "median_n_genes": float(values["n_genes_by_counts"]),
                "median_pct_mt": (
                    float(values["pct_counts_mt"]) if np.isfinite(values["pct_counts_mt"]) else None
                ),
                "median_top20": float(values["pct_counts_in_top_20_genes"]),
                "absolute_axes": absolute_axes,
                "relative_axes": relative_axes,
                "reason": (
                    "At least two independent sample-quality axes indicate catastrophic failure."
                    if status == "BLOCKED"
                    else "A relative or single-axis sample anomaly requires review but is not catastrophic."
                    if status == "REVIEW"
                    else "No joint catastrophic sample-quality pattern was detected."
                ),
            }
        )
    return rows


def _expert_global(metrics: pd.DataFrame, rules: dict[str, float], *, tumor: bool) -> np.ndarray:
    mt_limit = max(rules["global_max_mt"], 40.0) if tumor else rules["global_max_mt"]
    axes = np.column_stack(
        [
            metrics["n_genes_by_counts"].to_numpy() < rules["global_min_genes"],
            metrics["pct_counts_mt"].to_numpy() > mt_limit,
            metrics["pct_counts_in_top_20_genes"].to_numpy() > rules["global_max_top20"],
        ]
    )
    return np.nansum(axes, axis=1) >= 2


def _per_sample_mad(metrics: pd.DataFrame, sample_values: pd.Series) -> np.ndarray:
    flagged = np.zeros(len(metrics), dtype=bool)
    samples = sample_values.astype(str).to_numpy()
    for sample in pd.unique(samples):
        mask = samples == sample
        sub = metrics.iloc[np.where(mask)[0]]
        axes = []
        for column, direction in (
            ("n_genes_by_counts", "lower"),
            ("pct_counts_mt", "upper"),
            ("pct_counts_in_top_20_genes", "upper"),
        ):
            axes.append(_robust_z(sub[column].to_numpy(float), direction=direction) >= 3.0)
        flagged[mask] = np.sum(np.column_stack(axes), axis=1) >= 2
    return flagged


def _joint_mixture(metrics: pd.DataFrame) -> tuple[np.ndarray, str]:
    mt = metrics["pct_counts_mt"].to_numpy(float)
    genes = metrics["n_genes_by_counts"].to_numpy(float)
    valid = np.isfinite(mt) & np.isfinite(genes)
    flagged = np.zeros(len(metrics), dtype=bool)
    if valid.sum() < 100 or np.unique(mt[valid]).size < 5 or np.unique(genes[valid]).size < 5:
        return flagged, "NOT_EVALUABLE"
    transformed = np.column_stack(
        [
            np.log1p(genes[valid]),
            np.log(
                (np.clip(mt[valid], 0.01, 99.99) / 100.0)
                / (1 - np.clip(mt[valid], 0.01, 99.99) / 100.0)
            ),
        ]
    )
    transformed = (transformed - transformed.mean(axis=0)) / np.maximum(
        transformed.std(axis=0), 1e-8
    )
    try:
        model = GaussianMixture(
            n_components=2, covariance_type="full", random_state=0, reg_covar=1e-5
        )
        model.fit(transformed)
        quality_axis = -model.means_[:, 0] + model.means_[:, 1]
        compromised = int(np.argmax(quality_axis))
        reference = 1 - compromised
        genes_shift = float(model.means_[compromised, 0] - model.means_[reference, 0])
        mt_shift = float(model.means_[compromised, 1] - model.means_[reference, 1])
        probabilities = model.predict_proba(transformed)[:, compromised]
        valid_index = np.where(valid)[0]
        flagged[valid_index[probabilities >= 0.75]] = True
        flagged_fraction = float(flagged.mean()) if len(flagged) else 0.0
        if genes_shift > -0.5 or mt_shift < 0.5 or flagged_fraction > 0.25:
            return flagged, "UNRELIABLE_FIT_REVIEW_ONLY"
    except (ValueError, np.linalg.LinAlgError):
        return flagged, "NOT_EVALUABLE"
    return flagged, "EVALUATED_REVIEW_ONLY"


def _multisample_robust(metrics: pd.DataFrame, sample_values: pd.Series) -> np.ndarray:
    samples = sample_values.astype(str).to_numpy()
    axes = np.zeros((len(metrics), 3), dtype=float)
    specs = (
        ("n_genes_by_counts", "lower"),
        ("pct_counts_mt", "upper"),
        ("pct_counts_in_top_20_genes", "upper"),
    )
    for sample in pd.unique(samples):
        mask = samples == sample
        positions = np.where(mask)[0]
        for axis, (column, direction) in enumerate(specs):
            axes[positions, axis] = _robust_z(
                metrics.iloc[positions][column].to_numpy(float), direction=direction
            )
    axes = np.maximum(axes, 0.0)
    distance = np.sqrt(np.sum(axes**2, axis=1))
    return (distance >= 4.0) & ((axes >= 2.0).sum(axis=1) >= 2)


def _program_fraction(adata: AnnData, positions: np.ndarray, genes: tuple[str, ...]) -> np.ndarray:
    lookup = {str(name).upper(): idx for idx, name in enumerate(adata.var_names)}
    gene_positions = [lookup[gene] for gene in genes if gene in lookup]
    if not gene_positions:
        return np.full(len(positions), np.nan)
    matrix = _matrix(adata)[positions, :]
    totals = _row_sum(matrix)
    program = _row_sum(matrix[:, gene_positions])
    return np.divide(100.0 * program, totals, out=np.zeros_like(totals), where=totals > 0)


def _quick_map_review(
    adata: AnnData,
    metrics: pd.DataFrame,
    decisions: np.ndarray,
    rules: dict[str, float],
    sample_values: pd.Series,
    *,
    max_cells: int = 5000,
) -> dict[str, Any]:
    """Build a bounded temporary map for cluster-level sensitivity evidence."""
    if adata.n_obs < 500 or adata.n_vars < 50:
        return {
            "status": "NOT_EVALUABLE",
            "reason": "At least 500 cells and 50 genes are required for the temporary map.",
            "action": "Retain candidate conflicts for review; do not infer cluster evidence.",
        }
    matrix = _matrix(adata)
    if not sparse.issparse(matrix) and min(adata.n_obs, max_cells) * adata.n_vars > 20_000_000:
        return {
            "status": "NOT_EVALUABLE",
            "reason": "Dense temporary map would exceed the bounded review memory budget.",
            "action": "Create a project-level quick map after a safe sparse conversion.",
        }

    samples = sample_values.astype(str).to_numpy()
    rng = np.random.default_rng(0)
    if adata.n_obs <= max_cells:
        positions = np.arange(adata.n_obs)
    else:
        selected: list[int] = []
        unique_samples = pd.unique(samples)
        quota = max(1, max_cells // len(unique_samples))
        for sample in unique_samples:
            available = np.flatnonzero(samples == sample)
            take = min(quota, len(available))
            selected.extend(rng.choice(available, size=take, replace=False).tolist())
        remaining = max_cells - len(selected)
        if remaining > 0:
            pool = np.setdiff1d(np.arange(adata.n_obs), np.asarray(selected), assume_unique=False)
            selected.extend(
                rng.choice(pool, size=min(remaining, len(pool)), replace=False).tolist()
            )
        positions = np.sort(np.asarray(selected, dtype=int))

    values = matrix[positions, :].astype(float)
    totals = _row_sum(values)
    scale = np.divide(1e4, totals, out=np.zeros_like(totals), where=totals > 0)
    if sparse.issparse(values):
        normalized = values.multiply(scale[:, None]).tocsr()
        normalized.data = np.log1p(normalized.data)
        mean = np.asarray(normalized.mean(axis=0)).ravel()
        mean_sq = np.asarray(normalized.multiply(normalized).mean(axis=0)).ravel()
        variance = np.maximum(mean_sq - mean**2, 0.0)
    else:
        normalized = np.log1p(np.asarray(values) * scale[:, None])
        variance = normalized.var(axis=0)
    n_features = min(1000, adata.n_vars)
    feature_positions = np.argsort(variance, kind="mergesort")[-n_features:]
    from sklearn.cluster import MiniBatchKMeans
    from sklearn.decomposition import TruncatedSVD

    n_components = min(20, len(positions) - 1, n_features - 1)
    embedding = TruncatedSVD(n_components=n_components, random_state=0).fit_transform(
        normalized[:, feature_positions]
    )
    n_clusters = min(15, max(3, int(round(np.sqrt(len(positions) / 200)))))
    clusters = MiniBatchKMeans(
        n_clusters=n_clusters,
        random_state=0,
        n_init=5,
        batch_size=min(1024, len(positions)),
    ).fit_predict(embedding)

    sampled_metrics = metrics.iloc[positions].copy()
    sampled_metrics["qc_decision"] = decisions[positions]
    sampled_metrics["cluster"] = clusters.astype(str)
    sampled_metrics["stress_fraction"] = _program_fraction(
        adata,
        positions,
        ("FOS", "JUN", "JUNB", "HSPA1A", "HSPA1B", "DNAJB1", "ATF3"),
    )
    sampled_metrics["apoptosis_fraction"] = _program_fraction(
        adata,
        positions,
        ("BAX", "BAK1", "CASP3", "CASP7", "PMAIP1", "BBC3"),
    )
    if "predicted_doublet" in adata.obs:
        sampled_metrics["predicted_doublet"] = (
            adata.obs["predicted_doublet"].iloc[positions].fillna(False).to_numpy(bool)
        )

    rows: list[dict[str, Any]] = []
    for cluster, frame in sampled_metrics.groupby("cluster", observed=True):
        median_genes = float(frame["n_genes_by_counts"].median())
        median_mt = float(frame["pct_counts_mt"].median())
        median_top20 = float(frame["pct_counts_in_top_20_genes"].median())
        low_quality_axes = sum(
            [
                median_genes < rules["global_min_genes"],
                np.isfinite(median_mt) and median_mt > rules["global_max_mt"],
                median_top20 > rules["global_max_top20"],
            ]
        )
        rows.append(
            {
                "cluster": str(cluster),
                "n_sampled_cells": int(len(frame)),
                "median_n_genes": median_genes,
                "median_pct_mt": median_mt if np.isfinite(median_mt) else None,
                "median_top20": median_top20,
                "remove_fraction": float((frame["qc_decision"] == "REMOVE").mean()),
                "review_fraction": float((frame["qc_decision"] == "REVIEW").mean()),
                "median_stress_fraction": (
                    float(frame["stress_fraction"].median())
                    if frame["stress_fraction"].notna().any()
                    else None
                ),
                "median_apoptosis_fraction": (
                    float(frame["apoptosis_fraction"].median())
                    if frame["apoptosis_fraction"].notna().any()
                    else None
                ),
                "doublet_fraction": (
                    float(frame["predicted_doublet"].mean())
                    if "predicted_doublet" in frame
                    else None
                ),
                "status": "REVIEW" if low_quality_axes >= 2 else "KEEP",
            }
        )
    suspicious = [row["cluster"] for row in rows if row["status"] == "REVIEW"]
    return {
        "status": "REVIEW" if suspicious else "READY",
        "method": "temporary_stratified_log1p_svd_minibatch_kmeans",
        "n_sampled_cells": int(len(positions)),
        "n_features": int(n_features),
        "n_clusters": int(n_clusters),
        "suspicious_clusters": suspicious,
        "rows": rows,
        "action": (
            "Inspect suspicious temporary clusters; these proxy labels never delete cells."
            if suspicious
            else "No jointly low-quality temporary cluster was detected; continue lineage review."
        ),
    }


def recommend_evidence_calibrated_qc(
    adata: AnnData,
    context: AnalysisContext | dict[str, Any],
) -> DecisionCard:
    """Return a read-only, evidence-calibrated QC decision card."""
    resolved = infer_analysis_context(adata, context=context)
    profile = _profile_name(resolved)
    rules = dict(_PROFILE_RULES[profile])
    metrics, provenance = _metric_frame(adata, resolved)

    blockers: list[str] = []
    if resolved.is_multi_sample and not resolved.sample_key:
        blockers.append("Multi-sample QC requires an explicit, valid sample_key.")
    if resolved.sample_key and resolved.sample_key not in adata.obs:
        blockers.append(f"sample_key {resolved.sample_key!r} is absent from adata.obs.")

    if resolved.sample_key and resolved.sample_key in adata.obs:
        sample_values = adata.obs[resolved.sample_key].astype(str)
    else:
        sample_values = pd.Series("__single_sample__", index=adata.obs_names)
    sample_decisions = _sample_review(metrics, sample_values, rules)
    failed_samples = [row["sample"] for row in sample_decisions if row["status"] == "BLOCKED"]
    review_samples = [row["sample"] for row in sample_decisions if row["status"] == "REVIEW"]
    if failed_samples:
        blockers.append("Catastrophic sample-quality pattern: " + ", ".join(failed_samples))

    candidates: list[tuple[str, np.ndarray, str, str]] = []
    candidates.append(
        (
            "expert_global",
            _expert_global(metrics, rules, tumor=resolved.enables_tumor_module),
            "EVALUATED",
            "Protocol-profiled global baseline; not asserted to be universally optimal.",
        )
    )
    candidates.append(
        (
            "per_sample_mad",
            _per_sample_mad(metrics, sample_values),
            "EVALUATED",
            "Robust per-sample baseline requiring two abnormal QC axes.",
        )
    )
    mixture_mask, mixture_status = _joint_mixture(metrics)
    candidates.append(
        (
            "miqc_family",
            mixture_mask,
            mixture_status,
            "Internal joint-mixture sensitivity proxy only; not a miQC implementation and never an automatic REMOVE vote.",
        )
    )
    candidates.append(
        (
            "sampleqc_family",
            _multisample_robust(metrics, sample_values),
            "EXPERIMENTAL_PROXY_REVIEW_ONLY",
            "Sample-shifted robust sensitivity proxy; not a SampleQC implementation and never an automatic REMOVE vote.",
        )
    )

    baseline_masks = [candidates[0][1], candidates[1][1]]
    votes = np.sum(np.column_stack(baseline_masks), axis=1)
    decisions = np.full(len(metrics), "KEEP", dtype=object)
    decisions[votes == 1] = "REVIEW"
    decisions[votes >= 2] = "REMOVE"
    for _, mask, candidate_status, _ in candidates[2:]:
        if candidate_status in {"EVALUATED_REVIEW_ONLY", "EXPERIMENTAL_PROXY_REVIEW_ONLY"}:
            decisions[mask & (decisions == "KEEP")] = "REVIEW"

    doublet_available = "predicted_doublet" in adata.obs
    if doublet_available:
        doublets = adata.obs["predicted_doublet"].fillna(False).to_numpy(bool)
        decisions[doublets & (decisions == "KEEP")] = "REVIEW"
    stress_available = "stress_high" in adata.obs
    if stress_available:
        stress = adata.obs["stress_high"].fillna(False).to_numpy(bool)
        decisions[stress & (decisions == "KEEP")] = "REVIEW"

    quick_map = _quick_map_review(
        adata,
        metrics,
        decisions,
        rules,
        sample_values,
    )

    lineage_review: list[dict[str, Any]] = []
    lineage_key = resolved.cell_type_key
    if lineage_key and lineage_key in adata.obs:
        labels = adata.obs[lineage_key].astype(str).to_numpy()
        rare_limit = max(50, int(np.ceil(0.01 * adata.n_obs)))
        for label in pd.unique(labels):
            mask = labels == label
            n_remove = int(np.sum(decisions[mask] == "REMOVE"))
            if mask.sum() <= rare_limit and n_remove / max(1, mask.sum()) > 0.2:
                positions = np.where(mask & (decisions == "REMOVE"))[0]
                decisions[positions] = "REVIEW"
                lineage_review.append(
                    {
                        "lineage": str(label),
                        "n_cells": int(mask.sum()),
                        "downgraded_remove_to_review": int(len(positions)),
                    }
                )

    missing_evidence: list[str] = []
    if not provenance["mt_available"]:
        missing_evidence.append("Mitochondrial QC metric could not be established.")
    ambient_evaluable = resolved.input_provenance == "unfiltered_droplets"
    if not ambient_evaluable:
        missing_evidence.append(
            "Ambient RNA and cell-calling evidence are NOT_EVALUABLE without unfiltered droplets."
        )
    if not doublet_available:
        missing_evidence.append(
            "Doublet evidence was not supplied; no automatic doublet removal was made."
        )
    if quick_map["status"] == "NOT_EVALUABLE":
        missing_evidence.append("Quick-map cluster sensitivity was NOT_EVALUABLE.")

    candidate_records = [
        {
            "name": name,
            "status": status,
            "flagged_cells": int(mask.sum()),
            "flagged_fraction": float(mask.mean()) if len(mask) else 0.0,
            "note": note,
        }
        for name, mask, status, note in candidates
    ]
    policy_candidate_records = [
        {
            **record,
            "flagged_obs_names": metrics.index[mask].tolist(),
        }
        for record, (_, mask, _, _) in zip(candidate_records, candidates)
    ]
    selector_remove = decisions == "REMOVE"
    counterfactual_comparison = []
    for record, (_, mask, _, _) in zip(candidate_records, candidates):
        counterfactual_comparison.append(
            {
                "name": record["name"],
                "status": record["status"],
                "candidate_flagged_cells": int(mask.sum()),
                "additional_vs_selector": int(np.sum(mask & ~selector_remove)),
                "selector_remove_not_flagged": int(np.sum(selector_remove & ~mask)),
            }
        )
    policy_status = (
        "BLOCKED"
        if blockers
        else "REVIEW"
        if (decisions == "REVIEW").any() or missing_evidence
        else "READY"
    )
    policy_id = f"QCP-{uuid.uuid4().hex[:12]}"
    remove_names = metrics.index[decisions == "REMOVE"].tolist()
    review_names = metrics.index[decisions == "REVIEW"].tolist()
    policy = QCPolicy(
        policy_id=policy_id,
        status=policy_status,
        context=resolved,
        profile=profile,
        sample_key=resolved.sample_key,
        input_fingerprint=_fingerprint(adata),
        sample_decisions=sample_decisions,
        candidate_policies=policy_candidate_records,
        remove_obs_names=remove_names,
        review_obs_names=review_names,
        blockers=blockers,
        missing_evidence=missing_evidence,
        evidence_heads={
            "cell_quality": {"status": "EVALUATED", "candidates": candidate_records},
            "sample_quality": {
                "status": "BLOCKED" if failed_samples else "REVIEW" if review_samples else "READY"
            },
            "doublet": {
                "status": "AVAILABLE" if doublet_available else "NOT_EVALUABLE",
                "action": "review_only",
            },
            "ambient": {
                "status": "EVALUABLE" if ambient_evaluable else "NOT_EVALUABLE",
                "action": "no_cell_removal",
            },
            "stress": {
                "status": "AVAILABLE" if stress_available else "NOT_EVALUABLE",
                "action": "review_only",
            },
            "quick_map": quick_map,
            "lineage_sensitivity": {
                "status": "AVAILABLE" if lineage_key else "NOT_EVALUABLE",
                "rows": lineage_review,
            },
        },
        execution={"filter_remove": True, "decision_key": "qc_decision"},
        claim_boundary={
            "supported": [
                "Candidate disagreement and sample-level catastrophic patterns were audited."
            ],
            "exploratory": [
                "Consensus REMOVE calls are a review policy pending external validation."
            ],
            "unsupported": ["Universal superiority over traditional QC."],
        },
    )
    if policy_status == "BLOCKED":
        next_action = "Resolve or exclude the blocked sample/library before cell-level filtering."
    elif policy_status == "REVIEW":
        next_action = "Inspect REVIEW cells and lineage impact, then explicitly apply the policy."
    else:
        next_action = "Explicitly apply the reviewed QC policy."
    return DecisionCard(
        stage="qc",
        status=policy_status,
        decision="qc_policy",
        recommended="consensus_remove_disagreement_review",
        reason=(
            blockers[0]
            if blockers
            else "Independent QC candidates were compared without an opaque aggregate score."
        ),
        evidence=[
            "sample-level joint catastrophic gate",
            "expert global baseline",
            "per-sample robust MAD baseline",
            "joint genes-mitochondrial mixture candidate",
            "sample-shifted multivariate candidate",
            "bounded temporary cluster sensitivity map",
        ],
        next_action=next_action,
        rerun_scope="QC -> preprocessing -> all downstream stages",
        priority="blocker" if policy_status == "BLOCKED" else "review",
        source="evidence_calibrated_qc_v0.1-draft",
        candidates=candidate_records,
        affected={
            "remove_cells": len(remove_names),
            "review_cells": len(review_names),
            "blocked_samples": failed_samples,
            "review_samples": review_samples,
            "quick_map_review_clusters": quick_map.get("suspicious_clusters", []),
        },
        comparison=counterfactual_comparison,
        uncertainty=missing_evidence,
        missing_evidence=missing_evidence,
        sensitivity=["Reassess REMOVE/REVIEW by provisional lineage before formal analysis."],
        policy=policy,
        details={
            "metric_provenance": provenance,
            "sample_decisions": sample_decisions,
            "profile_source": (
                "scLucid pre-registered protocol safety guard v0.1; these are catastrophic "
                "warning regions, not universal cell-filtering thresholds."
            ),
            "quick_map": quick_map,
        },
    )


def apply_evidence_calibrated_qc(adata: AnnData, policy: QCPolicy) -> RunEvidence:
    """Apply a reviewed QC policy and return execution evidence."""
    if policy.status == "BLOCKED":
        raise RuntimeError("Blocked QC policies cannot be applied: " + "; ".join(policy.blockers))
    if _fingerprint(adata) != policy.input_fingerprint:
        raise ValueError("QC policy input fingerprint does not match this AnnData object.")

    result = adata.copy()
    remove = result.obs_names.astype(str).isin(policy.remove_obs_names)
    review = result.obs_names.astype(str).isin(policy.review_obs_names)
    decisions = np.full(result.n_obs, "KEEP", dtype=object)
    decisions[review] = "REVIEW"
    decisions[remove] = "REMOVE"
    result.obs["qc_decision"] = pd.Categorical(
        decisions, categories=["KEEP", "REVIEW", "REMOVE"], ordered=True
    )
    result.obs["qc_remove"] = remove
    n_before = int(result.n_obs)
    if policy.execution.get("filter_remove", True):
        result = result[~remove, :].copy()
    n_after = int(result.n_obs)

    run_id = f"QCR-{uuid.uuid4().hex[:12]}"
    serializable = {
        "schema_version": "RunEvidence-0.1-draft",
        "run_id": run_id,
        "policy_id": policy.policy_id,
        "status": policy.status,
        "n_cells_before": n_before,
        "n_cells_after": n_after,
        "n_removed": n_before - n_after,
        "n_review": len(policy.review_obs_names),
        "claim_boundary": policy.claim_boundary,
    }
    result.uns.setdefault("sclucid", {}).setdefault("qc", {})["policy_run_evidence"] = (
        sanitize_for_hdf5(serializable)
    )
    return RunEvidence(
        evidence_id=f"E-{policy.policy_id}",
        run_id=run_id,
        stage="qc",
        status=policy.status,
        policy=policy,
        adata=result,
        artifact={
            "path_or_key": 'adata.uns["sclucid"]["qc"]["policy_run_evidence"]',
            "representation": "counts",
        },
        result=serializable,
        supports=["The reviewed policy was applied to the exact fingerprinted input."],
        challenges=[],
        limitations=list(policy.missing_evidence),
    )


__all__ = ["recommend_evidence_calibrated_qc", "apply_evidence_calibrated_qc"]
