"""Evidence-calibrated preprocessing review and explicit policy execution.

The selector is intentionally conservative: it records which methods are
applicable, but it recommends the simple unintegrated baseline until an
empirical comparison demonstrates a Pareto improvement.
"""

from __future__ import annotations

import hashlib
import uuid
from importlib.util import find_spec
from typing import Any, Literal

import numpy as np
import scanpy as sc
import scipy.sparse as sparse
from anndata import AnnData

from ..decision import DecisionCard, PreprocessPolicy, RunEvidence
from ..utils.context import AnalysisContext, infer_analysis_context
from ..utils.sanitize import sanitize_for_hdf5
from .config import HVGConfig, IntegrationConfig, NormalizationConfig
from .hvg import find_hvgs
from .integrate import batch_correction, detect_integration_confounding
from .normalize import normalize_data


def _fingerprint(adata: AnnData) -> dict[str, Any]:
    digest = hashlib.sha256()
    for name in adata.obs_names.astype(str):
        digest.update(name.encode("utf-8"))
        digest.update(b"\0")
    matrix = adata.layers.get("counts", adata.X)
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


def _looks_like_counts(matrix) -> bool:
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    if values.size > 100_000:
        values = values[:100_000]
    return bool(
        values.size
        and np.all(np.isfinite(values))
        and np.min(values) >= 0
        and np.allclose(values, np.rint(values), atol=1e-6)
    )


def _is_umi_assay(context: AnalysisContext) -> bool:
    token = str(context.assay or "").lower().replace("-", "").replace("_", "")
    return not any(label in token for label in ("smartseq", "fullength", "platebased"))


def _integration_review(adata: AnnData, context: AnalysisContext) -> dict[str, Any]:
    batch_key = context.batch_key
    if not batch_key or batch_key not in adata.obs:
        return {
            "status": "NOT_EVALUABLE",
            "decision": "do_not_integrate",
            "reason": "No explicit technical batch key is available; unintegrated is the baseline.",
            "batch_key": batch_key,
            "confounding": [],
        }
    if adata.obs[batch_key].nunique(dropna=True) < 2:
        return {
            "status": "READY",
            "decision": "do_not_integrate",
            "reason": "Only one technical batch is present.",
            "batch_key": batch_key,
            "confounding": [],
        }

    biology_columns = [
        key
        for key in (context.condition_key, context.cell_type_key)
        if key and key in adata.obs and key != batch_key
    ]
    confounding = detect_integration_confounding(
        adata,
        batch_key,
        biology_columns,
    )
    if confounding:
        return {
            "status": "BLOCKED",
            "decision": "do_not_integrate",
            "reason": "Batch is strongly confounded with protected biology.",
            "batch_key": batch_key,
            "confounding": confounding,
        }
    return {
        "status": "REVIEW",
        "decision": "compare_before_integrating",
        "reason": (
            "Multiple batches exist, but no quantitative Pareto comparison has shown "
            "that integration improves batch mixing without unacceptable biology loss."
        ),
        "batch_key": batch_key,
        "confounding": [],
    }


def recommend_preprocess_policy(
    adata: AnnData,
    context: AnalysisContext | dict[str, Any],
    *,
    consumer: Literal[
        "exploration", "annotation", "integration", "expression_inference"
    ] = "exploration",
) -> DecisionCard:
    """Review preprocessing candidates without modifying ``adata``."""
    resolved = infer_analysis_context(adata, context=context)
    counts_source = "layers[counts]" if "counts" in adata.layers else "X"
    counts_matrix = adata.layers.get("counts", adata.X)
    counts_valid = _looks_like_counts(counts_matrix)
    umi_assay = _is_umi_assay(resolved)

    blockers: list[str] = []
    if not counts_valid:
        blockers.append("A non-negative integer count matrix is required in layers['counts'] or X.")
    if resolved.is_multi_sample and not resolved.sample_key:
        blockers.append("Multi-sample preprocessing requires an explicit sample_key.")
    if resolved.sample_key and resolved.sample_key not in adata.obs:
        blockers.append(f"sample_key {resolved.sample_key!r} is absent from adata.obs.")

    normalization_candidates = [
        {
            "name": "library_size_log1p",
            "status": "AVAILABLE" if counts_valid else "BLOCKED",
            "selected": True,
            "space": "normalized_full and discovery baseline",
            "reason": "Simple reference baseline; complexity must demonstrate empirical gain.",
        },
        {
            "name": "analytic_pearson_residuals",
            "status": "AVAILABLE" if counts_valid and umi_assay else "NOT_APPLICABLE",
            "selected": False,
            "space": "discovery only",
            "reason": "Count-model residual representation restricted to UMI discovery tasks.",
        },
        {
            "name": "scran",
            "status": (
                "OPTIONAL_NOT_VERIFIED"
                if find_spec("rpy2") is not None
                else "OPTIONAL_DEPENDENCY_MISSING"
            ),
            "selected": False,
            "space": "normalized_full",
            "reason": "Never silently substituted when the R/scran dependency is unavailable.",
        },
    ]
    feature_candidates = [
        {
            "name": "batch_aware_hvg",
            "status": "AVAILABLE" if counts_valid else "BLOCKED",
            "selected": True,
            "space": "unsupervised discovery",
            "reason": "Uses sample-aware recurrence when a sample key is available.",
        },
        {
            "name": "multinomial_deviance",
            "status": "AVAILABLE" if counts_valid and umi_assay else "NOT_APPLICABLE",
            "selected": False,
            "space": "unsupervised sensitivity",
            "reason": "Requires held-out task evidence before replacing the simple baseline.",
        },
        {
            "name": "protected_marker_union",
            "status": "DISALLOWED_DEFAULT",
            "selected": False,
            "space": "hypothesis sensitivity only",
            "reason": "Prior markers must not shape the default unsupervised clustering space.",
        },
    ]
    integration = _integration_review(adata, resolved)
    if consumer == "integration" and integration["status"] == "BLOCKED":
        blockers.extend(integration["confounding"])

    missing_evidence = [
        "Candidate transformations have not been compared on this dataset.",
        "Feature-selection regret on held-out data has not been measured.",
    ]
    if integration["status"] == "REVIEW":
        missing_evidence.append(
            "Batch removal, biology conservation, rare-population retention, and stability lack a Pareto comparison."
        )

    status = "BLOCKED" if blockers else "REVIEW"
    policy = PreprocessPolicy(
        policy_id=f"PPP-{uuid.uuid4().hex[:12]}",
        status=status,
        context=resolved,
        consumer=consumer,
        input_fingerprint=_fingerprint(adata),
        layer_contract={
            "counts": "layers[counts]",
            "normalized_full": "layers[normalized_full] and raw",
            "discovery_rep": "obsm[X_pca]",
            "integrated_rep": "not_selected",
            "marker_program_source": "layers[normalized_full]",
            "formal_count_model_source": "layers[counts]",
        },
        normalization_method="standard",
        feature_selection_method="scanpy",
        run_integration=False,
        integration_method=None,
        candidates={
            "normalization": normalization_candidates,
            "feature_selection": feature_candidates,
            "integration": [integration],
        },
        blockers=blockers,
        missing_evidence=missing_evidence,
        execution={
            "counts_source": counts_source,
            "n_top_genes": min(2000, max(100, int(adata.n_vars))),
            "sample_aware_hvg_key": resolved.sample_key,
            "integration_allowed": integration["status"] != "BLOCKED",
            "integration_review": integration,
        },
        claim_boundary={
            "supported": ["The representation contract and method applicability were audited."],
            "exploratory": [
                "The simple unintegrated baseline is selected pending empirical comparison."
            ],
            "unsupported": ["Universal superiority of the selected preprocessing strategy."],
        },
    )
    next_action = (
        "Restore a trustworthy count matrix or required sample metadata, then rerun review."
        if status == "BLOCKED"
        else "Apply the simple unintegrated baseline and inspect the resulting RunEvidence."
    )
    return DecisionCard(
        stage="preprocess",
        status=status,
        decision="preprocess_policy",
        recommended="standard_log1p + batch-aware_HVG + unintegrated_PCA",
        reason=(
            blockers[0]
            if blockers
            else "No complex candidate has yet demonstrated a dataset-specific Pareto improvement."
        ),
        evidence=[
            "permanent count-space contract",
            "full-gene interpretation space",
            "unsupervised discovery space",
            "unintegrated baseline",
        ],
        next_action=next_action,
        rerun_scope="preprocessing -> all downstream stages",
        priority="blocker" if status == "BLOCKED" else "review",
        source="evidence_calibrated_preprocess_v0.1-draft",
        candidates=normalization_candidates + feature_candidates + [integration],
        affected={
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "consumer": consumer,
        },
        comparison=[
            {
                "objective": "Pareto comparison",
                "status": "PENDING",
                "metrics": [
                    "batch_removal",
                    "biology_conservation",
                    "rare_population_retention",
                    "stability",
                ],
            }
        ],
        uncertainty=missing_evidence,
        missing_evidence=missing_evidence,
        sensitivity=[
            "Compare Pearson residuals and deviance features without replacing the full-gene interpretation space."
        ],
        policy=policy,
        details={"integration_review": integration, "counts_source": counts_source},
    )


def apply_preprocess_policy(
    adata: AnnData,
    policy: PreprocessPolicy | DecisionCard,
) -> RunEvidence:
    """Execute an explicit preprocessing policy on the fingerprinted input."""
    if isinstance(policy, DecisionCard):
        policy = policy.policy
    if not isinstance(policy, PreprocessPolicy):
        raise TypeError("policy must be a PreprocessPolicy or a DecisionCard containing one.")
    if policy.status == "BLOCKED":
        raise RuntimeError(
            "Blocked preprocessing policies cannot be applied: " + "; ".join(policy.blockers)
        )
    if _fingerprint(adata) != policy.input_fingerprint:
        raise ValueError("Preprocess policy input fingerprint does not match this AnnData object.")
    if policy.run_integration and not policy.execution.get("integration_allowed", False):
        raise RuntimeError(
            "Integration is blocked by the recorded batch-biology confounding review."
        )

    result = adata.copy()
    if "counts" not in result.layers:
        if not _looks_like_counts(result.X):
            raise ValueError("Cannot establish layers['counts'] from a transformed X matrix.")
        result.layers["counts"] = result.X.copy()

    full_norm_method = "scran" if policy.normalization_method == "scran" else "standard"
    normalize_data(
        result,
        config=NormalizationConfig(
            method=full_norm_method,
            input_layer="counts",
            output_layer="normalized_full",
            update_X=False,
            set_raw=False,
            plot=False,
            report=False,
        ),
        force=True,
    )
    raw_source = AnnData(
        X=result.layers["normalized_full"].copy(),
        obs=result.obs.copy(),
        var=result.var.copy(),
    )
    result.raw = raw_source

    discovery_layer = "normalized_full"
    if policy.normalization_method == "pearson_residuals":
        normalize_data(
            result,
            config=NormalizationConfig(
                method="pearson_residuals",
                input_layer="counts",
                output_layer="discovery_residuals",
                update_X=False,
                set_raw=False,
                plot=False,
                report=False,
            ),
            force=True,
        )
        discovery_layer = "discovery_residuals"

    hvg_input = "counts" if policy.feature_selection_method == "deviance" else discovery_layer
    hvg_config = HVGConfig(
        method=policy.feature_selection_method,
        n_top_genes=int(policy.execution.get("n_top_genes", 2000)),
        flavor="auto",
        batch_key=(
            policy.context.sample_key if policy.feature_selection_method == "scanpy" else None
        ),
        sample_key=policy.context.sample_key or "sampleID",
        protect_genes=False,
        protected_gene_presets=[],
        protected_gene_sets={},
        plot=False,
        report=False,
    )
    find_hvgs(
        result,
        config=hvg_config,
        input_layer=hvg_input,
        species=policy.context.species,
        force=True,
        plot=False,
    )
    hvg_key = result.uns["sclucid"]["preprocess"]["hvg"]["output_key"]
    discovery_mask = result.var[hvg_key].fillna(False).to_numpy(bool)
    if int(discovery_mask.sum()) < 2:
        raise RuntimeError("Feature selection produced fewer than two discovery features.")
    result.var["discovery_feature"] = discovery_mask

    discovery = result[:, discovery_mask].copy()
    discovery_source = result.layers[discovery_layer][:, discovery_mask]
    source_was_sparse = sparse.issparse(discovery_source)
    if source_was_sparse:
        # PCA is intentionally fitted on a bounded discovery-feature matrix.
        # Make the required centering allocation explicit so Scanpy cannot
        # silently densify a whole-expression representation.
        discovery.X = discovery_source.toarray()
    else:
        discovery.X = np.asarray(discovery_source).copy()
    sc.pp.scale(discovery, max_value=10.0, zero_center=True)
    n_comps = min(50, discovery.n_obs - 1, discovery.n_vars - 1)
    if n_comps < 1:
        raise RuntimeError("At least two cells and two discovery features are required for PCA.")
    dense_matrix_bytes = int(np.asarray(discovery.X).nbytes)
    pca_output_bytes = int(
        (discovery.n_obs + discovery.n_vars) * n_comps * np.dtype(np.float32).itemsize
    )
    discovery_temporary_contract = {
        "densification_occurred": bool(source_was_sparse),
        "scope": "temporary_discovery_feature_matrix",
        "source_storage": "sparse" if source_was_sparse else "dense",
        "temporary_shape": [int(discovery.n_obs), int(discovery.n_vars)],
        "temporary_dtype": str(np.asarray(discovery.X).dtype),
        "dense_matrix_bytes": dense_matrix_bytes,
        "estimated_peak_bytes": int(2 * dense_matrix_bytes + pca_output_bytes),
        "estimate_method": (
            "two dense discovery buffers plus float32 PCA scores and loadings; "
            "planning estimate, not measured RSS"
        ),
        "bounded_by": "selected discovery features",
        "n_top_genes_requested": int(policy.execution.get("n_top_genes", 2000)),
        "persistent": False,
        "consumer": "PCA_and_neighbor_graph_only",
        "marker_program_eligible": False,
        "expression_inference_eligible": False,
        "persistent_storage": {
            "counts_sparse": bool(sparse.issparse(result.layers["counts"])),
            "normalized_full_sparse": bool(sparse.issparse(result.layers["normalized_full"])),
        },
    }
    sc.tl.pca(discovery, n_comps=n_comps, svd_solver="arpack")
    result.obsm["X_pca"] = discovery.obsm["X_pca"].copy()
    result.uns["pca"] = discovery.uns["pca"].copy()
    full_loadings = np.zeros((result.n_vars, n_comps), dtype=np.float32)
    full_loadings[discovery_mask, :] = np.asarray(discovery.varm["PCs"], dtype=np.float32)
    result.varm["PCs"] = full_loadings

    graph_rep = "X_pca"
    integrated_rep = "not_selected"
    if policy.run_integration:
        if not policy.integration_method or not policy.context.batch_key:
            raise ValueError("Integration requires integration_method and context.batch_key.")
        integrated_rep = f"X_{policy.integration_method.lower()}"
        batch_correction(
            result,
            config=IntegrationConfig(
                method=policy.integration_method,
                batch_key=policy.context.batch_key,
                use_rep="X_pca",
                output_key=integrated_rep,
                auto_decide=False,
                evaluate=False,
                plot=False,
            ),
            force=True,
        )
        graph_rep = integrated_rep

    result.X = result.layers["normalized_full"].copy()
    if result.n_obs >= 3:
        sc.pp.neighbors(
            result,
            use_rep=graph_rep,
            n_neighbors=min(15, result.n_obs - 1),
            n_pcs=min(n_comps, 50),
        )
        sc.tl.umap(result, random_state=0)

    representation_contract = dict(policy.layer_contract)
    representation_contract["integrated_rep"] = integrated_rep
    representation_contract["graph_rep"] = f"obsm[{graph_rep}]"
    representation_contract["discovery_temporary_contract"] = discovery_temporary_contract
    run_id = f"PPR-{uuid.uuid4().hex[:12]}"
    serializable = {
        "schema_version": "RunEvidence-0.1-draft",
        "run_id": run_id,
        "policy_id": policy.policy_id,
        "status": policy.status,
        "consumer": policy.consumer,
        "normalization_method": policy.normalization_method,
        "feature_selection_method": policy.feature_selection_method,
        "hvg_key": hvg_key,
        "n_discovery_features": int(discovery_mask.sum()),
        "representation_contract": representation_contract,
        "discovery_temporary_contract": discovery_temporary_contract,
        "claim_boundary": policy.claim_boundary,
    }
    preprocess_ns = result.uns.setdefault("sclucid", {}).setdefault("preprocess", {})
    preprocess_ns["representation_contract"] = sanitize_for_hdf5(representation_contract)
    preprocess_ns["policy_run_evidence"] = sanitize_for_hdf5(serializable)
    return RunEvidence(
        evidence_id=f"E-{policy.policy_id}",
        run_id=run_id,
        stage="preprocess",
        status=policy.status,
        policy=policy,
        adata=result,
        artifact={
            "path_or_key": 'adata.uns["sclucid"]["preprocess"]["policy_run_evidence"]',
            "representation": graph_rep,
        },
        result=serializable,
        supports=[
            "Counts, full-gene normalized expression, and discovery representation were kept semantically distinct."
        ],
        challenges=[],
        limitations=list(policy.missing_evidence),
    )


__all__ = ["recommend_preprocess_policy", "apply_preprocess_policy"]
