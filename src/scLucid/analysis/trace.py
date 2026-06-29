"""Review-summary enrichment for the analysis workflow."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import pandas as pd
from anndata import AnnData

from scLucid.utils.contracts import _review_payload

from ..utils import sanitize_for_hdf5
from ..utils.evidence import EvidenceBundle, EvidenceItem, ReviewAction, model_to_dict

ANALYSIS_TRACE_SCHEMA_VERSION = "1.0"
ANALYSIS_MODULE_MATURITY_SCHEMA_VERSION = "1.0"

ANALYSIS_REQUIRED_REVIEW_SECTIONS = {
    "analysis_schema_version",
    "analysis_inference_policy",
    "analysis_claim_level_summary",
    "analysis_output_contract",
    "analysis_decision_summary",
    "analysis_reviewer_table",
    "preprocess_input_context",
    "clustering_evidence_summary",
    "annotation_evidence_summary",
    "annotation_consensus_summary",
    "posthoc_qc_review_summary",
    "malignancy_interpretation_summary",
    "analysis_readiness",
    "review_action_items",
    "evidence_bundle",
    "module_maturity",
}

ANALYSIS_STABLE_ENTRYPOINTS = (
    "scLucid.analysis.run_standard_analysis",
    "scLucid.analysis.run_clustering_review",
    "scLucid.analysis.run_annotation_evidence",
    "scLucid.analysis.build_annotation_consensus",
    "scLucid.tumor.malignancy.run_malignancy_interpretation",
    "scLucid.analysis.cluster_cells",
    "scLucid.analysis.find_markers",
    "scLucid.analysis.run_annotation",
)

ANALYSIS_EXPECTED_OUTPUTS = (
    "adata.obs['leiden_clusters']",
    "adata.uns['rank_genes_groups']",
    "adata.uns['sclucid']['analysis']['annotation']['annotation_review_table']",
    "adata.uns['sclucid']['analysis']['malignancy']['malignancy_interpretation_summary']",
    "adata.uns['sclucid']['analysis']['review_summary']",
)


def enrich_analysis_review_summary(
    summary: dict[str, Any],
    *,
    adata: AnnData,
    config: Any,
    successful_steps: list[str],
    cluster_key: str,
) -> dict[str, Any]:
    """Add benchmark-grade review fields to an analysis summary."""
    summary = dict(summary)
    preprocess_context = build_preprocess_input_context(adata)
    clustering = build_clustering_evidence_summary(adata, cluster_key)
    annotation = build_annotation_evidence_summary(adata, config=config)
    consensus = build_annotation_consensus_summary(adata, config=config)
    posthoc_qc = build_posthoc_qc_review_summary(adata, cluster_key=cluster_key)
    malignancy = build_malignancy_interpretation_summary(adata, config=config)
    inference_policy = build_analysis_inference_policy(
        config=config,
        successful_steps=successful_steps,
    )
    output_contract = build_analysis_output_contract(
        adata=adata,
        config=config,
        cluster_key=cluster_key,
        successful_steps=successful_steps,
        preprocess_context=preprocess_context,
    )
    decision_summary = build_analysis_decision_summary(
        adata=adata,
        config=config,
        successful_steps=successful_steps,
        cluster_key=cluster_key,
        inference_policy=inference_policy,
        output_contract=output_contract,
        clustering_summary=clustering,
        annotation_summary=annotation,
        consensus_summary=consensus,
        posthoc_qc_summary=posthoc_qc,
        malignancy_summary=malignancy,
    )
    claim_level_summary = build_analysis_claim_level_summary(
        adata=adata,
        successful_steps=successful_steps,
        clustering_summary=clustering,
        annotation_summary=annotation,
        consensus_summary=consensus,
        decision_summary=decision_summary,
    )
    readiness = build_analysis_readiness_assessment(
        adata=adata,
        successful_steps=successful_steps,
        cluster_key=cluster_key,
        preprocess_context=preprocess_context,
        clustering_summary=clustering,
        annotation_summary=annotation,
        consensus_summary=consensus,
        posthoc_qc_summary=posthoc_qc,
        malignancy_summary=malignancy,
        claim_level_summary=claim_level_summary,
    )
    actions = build_analysis_review_action_items(
        readiness=readiness,
        clustering_summary=clustering,
        annotation_summary=annotation,
        consensus_summary=consensus,
        posthoc_qc_summary=posthoc_qc,
        malignancy_summary=malignancy,
        claim_level_summary=claim_level_summary,
    )

    summary["analysis_schema_version"] = ANALYSIS_TRACE_SCHEMA_VERSION
    summary["analysis_inference_policy"] = inference_policy
    summary["analysis_claim_level_summary"] = claim_level_summary
    summary["analysis_output_contract"] = output_contract
    summary["analysis_decision_summary"] = decision_summary
    summary["analysis_reviewer_table"] = build_analysis_reviewer_table(decision_summary)
    summary["preprocess_input_context"] = preprocess_context
    summary["clustering_evidence_summary"] = clustering
    summary["annotation_evidence_summary"] = annotation
    summary["annotation_consensus_summary"] = consensus
    summary["posthoc_qc_review_summary"] = posthoc_qc
    summary["malignancy_interpretation_summary"] = malignancy
    summary["analysis_readiness"] = readiness
    summary["review_action_items"] = actions
    summary["evidence_bundle"] = build_analysis_evidence_bundle(summary)
    summary["module_maturity"] = build_analysis_module_maturity_assessment(summary)
    return _json_safe(summary)


def get_analysis_module_contract() -> dict[str, Any]:
    """Return the frozen analysis module maturity contract."""
    return {
        "schema_version": ANALYSIS_MODULE_MATURITY_SCHEMA_VERSION,
        "module": "analysis",
        "stable_entrypoints": list(ANALYSIS_STABLE_ENTRYPOINTS),
        "required_review_sections": sorted(ANALYSIS_REQUIRED_REVIEW_SECTIONS),
        "expected_outputs": list(ANALYSIS_EXPECTED_OUTPUTS),
        "canonical_namespace": 'adata.uns["sclucid"]["analysis"]',
        "readiness_key": "analysis_readiness",
        "inference_policy_key": "analysis_inference_policy",
        "claim_level_key": "analysis_claim_level_summary",
        "output_contract_key": "analysis_output_contract",
        "decision_summary_key": "analysis_decision_summary",
        "reviewer_table_key": "analysis_reviewer_table",
        "clustering_evidence_key": "clustering_evidence_summary",
        "annotation_evidence_key": "annotation_evidence_summary",
        "annotation_consensus_key": "annotation_consensus_summary",
        "posthoc_qc_review_key": "posthoc_qc_review_summary",
        "malignancy_interpretation_key": "malignancy_interpretation_summary",
        "preprocess_input_key": "preprocess_input_context",
    }


def build_preprocess_input_context(adata: AnnData) -> dict[str, Any]:
    """Summarize the preprocessing state consumed by analysis."""
    pp_ns = adata.uns.get("sclucid", {}).get("preprocess", {})
    review = _review_payload(pp_ns.get("review_summary", {})) if isinstance(pp_ns, Mapping) else {}
    readiness = review.get("preprocess_readiness", {}) if isinstance(review, Mapping) else {}
    maturity = review.get("module_maturity", {}) if isinstance(review, Mapping) else {}

    return _json_safe(
        {
            "available": bool(review),
            "preprocess_readiness_status": (
                readiness.get("status") if isinstance(readiness, Mapping) else None
            ),
            "preprocess_readiness_score": (
                readiness.get("score") if isinstance(readiness, Mapping) else None
            ),
            "preprocess_maturity_status": (
                maturity.get("status") if isinstance(maturity, Mapping) else None
            ),
            "normalized_layer_present": "normalized" in adata.layers,
            "pca_present": "X_pca" in adata.obsm,
            "neighbors_present": "neighbors" in adata.uns,
            "umap_present": "X_umap" in adata.obsm,
            "hvg_present": "highly_variable" in adata.var,
        }
    )


def build_analysis_inference_policy(
    *,
    config: Any,
    successful_steps: list[str],
) -> dict[str, Any]:
    """Describe conservative inference boundaries for standard analysis."""
    run_clustering_review = bool(getattr(config, "run_clustering_review", True))
    clustering_review_ran = "clustering_review" in successful_steps
    return _json_safe(
        {
            "schema_version": "analysis_inference_policy_v1",
            "claim_boundary": "exploratory_until_reviewed",
            "clustering_review": {
                "recommended": True,
                "enabled": run_clustering_review,
                "executed": clustering_review_ran,
                "can_disable": True,
                "risk_note": (
                    "Clustering resolution review is recommended before final annotation; "
                    "skipping is allowed for speed but should be documented."
                ),
            },
            "marker_discovery": {
                "default_inference_level": "cell_level_marker_discovery",
                "risk_note": (
                    "Cluster marker discovery is useful for annotation and review; it is not "
                    "a publication-grade condition effect test."
                ),
            },
            "condition_de": {
                "recommended_primary_method": "sample_level_pseudobulk",
                "recommended_entrypoint": "scLucid.analysis.run_pseudobulk_de",
                "cell_level_compare_policy": "exploratory_only",
                "fallback_to_cell_level_default": False,
                "risk_note": (
                    "Condition DE should default to sample-level pseudobulk. Cell-level "
                    "comparisons treat cells as independent observations and are exploratory."
                ),
            },
            "cell_level_compare": {
                "inference_level": "exploratory_cell_level",
                "allowed_use": "screening_visualization_and_hypothesis_generation",
                "publication_note": (
                    "Use biological-sample replicates and pseudobulk or mixed models for "
                    "formal condition-level inference."
                ),
            },
        }
    )


def build_analysis_output_contract(
    *,
    adata: AnnData,
    config: Any,
    cluster_key: str,
    successful_steps: list[str],
    preprocess_context: Mapping[str, Any],
) -> dict[str, Any]:
    """Describe the stable analysis outputs and their inference semantics."""
    annotation_key = (
        getattr(config.annotation, "key_added", "cell_type_auto")
        if getattr(config, "annotation", None) is not None
        else "cell_type_auto"
    )
    lineage_key = (
        getattr(config.annotation, "lineage_key", "celltype_lineage_auto")
        if getattr(config, "annotation", None) is not None
        else "celltype_lineage_auto"
    )
    de_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("de", {})
    de_keys = sorted(de_ns.keys()) if isinstance(de_ns, Mapping) else []
    graph_source = (
        "preprocess_neighbors"
        if preprocess_context.get("neighbors_present")
        else "analysis_or_external"
        if "neighbors" in adata.uns
        else "missing"
    )
    stage_contracts = [
        {
            "stage": "preprocess_handoff",
            "required": True,
            "present": bool(preprocess_context.get("pca_present")),
            "primary_slot": "adata.obsm['X_pca']",
            "inference_level": "input_contract",
            "risk_note": "Analysis assumes preprocessing produced PCA and, ideally, a graph.",
        },
        {
            "stage": "clustering",
            "required": True,
            "present": cluster_key in adata.obs,
            "primary_slot": f"adata.obs['{cluster_key}']",
            "inference_level": "unsupervised_structure",
            "risk_note": "Cluster labels are hypotheses until resolution and marker evidence are reviewed.",
        },
        {
            "stage": "marker_discovery",
            "required": False,
            "present": "rank_genes_groups" in adata.uns or any(key.endswith("_df") for key in de_keys),
            "primary_slot": "adata.uns['sclucid']['analysis']['de']",
            "inference_level": "cell_level_marker_discovery",
            "risk_note": "Cell-level markers support annotation and characterization, not formal condition inference.",
        },
        {
            "stage": "annotation",
            "required": False,
            "present": annotation_key in adata.obs,
            "primary_slot": f"adata.obs['{annotation_key}']",
            "canonical_alias": "adata.obs['cell_type']",
            "lineage_slot": f"adata.obs['{lineage_key}']",
            "inference_level": "evidence_consensus",
            "risk_note": "Final labels should be reviewed when marker/reference evidence conflicts.",
        },
        {
            "stage": "condition_de",
            "required": False,
            "present": any("pseudobulk" in key for key in de_keys),
            "primary_slot": "adata.uns['sclucid']['analysis']['de']['pseudobulk_de']",
            "inference_level": "sample_level_when_replicated",
            "risk_note": "Publication-grade condition DE requires biological-sample replication.",
        },
        {
            "stage": "posthoc_qc",
            "required": False,
            "present": True,
            "primary_slot": "review_summary['posthoc_qc_review_summary']",
            "inference_level": "review_gate",
            "risk_note": "Post-hoc QC flags clusters for review without automatically deleting cells.",
        },
    ]
    review_required = [
        row["stage"]
        for row in stage_contracts
        if row.get("required") and not bool(row.get("present"))
    ]
    return _json_safe(
        {
            "schema_version": "analysis_output_contract_v1",
            "canonical_flow": (
                "preprocess handoff -> clustering -> marker discovery -> "
                "annotation consensus -> optional DE/proportion/tumor interpretation"
            ),
            "cluster_key": cluster_key,
            "annotation_key": annotation_key,
            "lineage_key": lineage_key,
            "canonical_annotation_aliases": [
                "cell_type_auto",
                "cell_type",
                "celltype_lineage_auto",
                "celltype_lineage",
            ],
            "graph_source": graph_source,
            "successful_steps": list(successful_steps),
            "de_result_keys": de_keys,
            "stage_contracts": stage_contracts,
            "review_required": bool(review_required),
            "review_required_stages": review_required,
        }
    )


def _decision_row(
    *,
    step: str,
    decision: str,
    recommended: Any,
    applied: Any,
    source: str,
    confidence: float,
    affected_output: str,
    inference_level: str,
    biological_risk_note: str,
    review_required: bool = False,
    reason: str = "",
) -> dict[str, Any]:
    return {
        "step": step,
        "decision": decision,
        "recommended_value": recommended,
        "applied_value": applied,
        "source": source,
        "confidence": float(max(0.0, min(1.0, confidence))),
        "affected_output": affected_output,
        "inference_level": inference_level,
        "biological_risk_note": biological_risk_note,
        "review_required": bool(review_required),
        "reason": reason,
    }


def build_analysis_decision_summary(
    *,
    adata: AnnData,
    config: Any,
    successful_steps: list[str],
    cluster_key: str,
    inference_policy: Mapping[str, Any],
    output_contract: Mapping[str, Any],
    clustering_summary: Mapping[str, Any],
    annotation_summary: Mapping[str, Any],
    consensus_summary: Mapping[str, Any],
    posthoc_qc_summary: Mapping[str, Any],
    malignancy_summary: Mapping[str, Any],
    claim_level_summary: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build one reviewer-facing decision summary for analysis outputs."""
    annotation_key = output_contract.get("annotation_key", "cell_type_auto")
    condition_policy = inference_policy.get("condition_de", {})
    clustering_policy = inference_policy.get("clustering_review", {})
    clustering_review_ran = "clustering_review" in successful_steps
    annotation_present = bool(consensus_summary.get("final_obs_present"))
    annotation_review_required = bool(
        annotation_summary.get("needs_review_clusters", 0)
        or annotation_summary.get("low_confidence_clusters", 0)
        or consensus_summary.get("needs_review_cells", 0)
    )
    posthoc_review_required = bool(posthoc_qc_summary.get("review_required"))
    malignancy_enabled = bool(malignancy_summary.get("enabled"))
    malignancy_review_required = bool(malignancy_summary.get("review_required"))
    de_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("de", {})
    de_keys = sorted(de_ns.keys()) if isinstance(de_ns, Mapping) else []
    has_pseudobulk = any("pseudobulk" in key for key in de_keys)

    decisions = [
        _decision_row(
            step="clustering_review",
            decision="use" if clustering_review_ran else "recommended_not_run",
            recommended=True,
            applied=clustering_review_ran,
            source="analysis_inference_policy",
            confidence=0.9 if clustering_review_ran else 0.65,
            affected_output=cluster_key,
            inference_level="review_gate",
            biological_risk_note=str(clustering_policy.get("risk_note", "")),
            review_required=not clustering_review_ran,
            reason=(
                "Clustering review executed."
                if clustering_review_ran
                else "Clustering review is recommended before final annotation."
            ),
        ),
        _decision_row(
            step="clustering",
            decision="use" if cluster_key in adata.obs else "blocked",
            recommended=cluster_key,
            applied=cluster_key if cluster_key in adata.obs else None,
            source="workflow_config",
            confidence=0.85 if cluster_key in adata.obs else 0.1,
            affected_output=f"adata.obs['{cluster_key}']",
            inference_level="unsupervised_structure",
            biological_risk_note=(
                "Clusters are analysis hypotheses; small or QC-heavy clusters need review."
            ),
            review_required=bool(clustering_summary.get("review_required")),
            reason="; ".join(map(str, clustering_summary.get("review_required", []))),
        ),
        _decision_row(
            step="marker_discovery",
            decision="use" if "markers" in successful_steps else "not_run",
            recommended="cell_level_marker_discovery",
            applied="rank_genes_groups" if "markers" in successful_steps else None,
            source="analysis_inference_policy",
            confidence=0.8 if "markers" in successful_steps else 0.5,
            affected_output="adata.uns['sclucid']['analysis']['de']",
            inference_level="cell_level_marker_discovery",
            biological_risk_note=(
                inference_policy.get("marker_discovery", {}).get("risk_note", "")
                if isinstance(inference_policy.get("marker_discovery"), Mapping)
                else ""
            ),
            review_required=False,
            reason="Marker discovery is for annotation support, not formal condition inference.",
        ),
        _decision_row(
            step="annotation_consensus",
            decision="use" if annotation_present else "review_or_not_run",
            recommended="consensus_label_with_evidence",
            applied=annotation_key if annotation_present else None,
            source="annotation_evidence_summary",
            confidence=float(consensus_summary.get("mean_confidence") or 0.7 if annotation_present else 0.3),
            affected_output=f"adata.obs['{annotation_key}']",
            inference_level="evidence_consensus",
            biological_risk_note=(
                "Consensus labels combine available marker/reference evidence and should be "
                "manually reviewed when confidence is low or evidence conflicts."
            ),
            review_required=annotation_review_required or not annotation_present,
            reason=(
                f"{annotation_summary.get('needs_review_clusters', 0)} cluster(s) require review."
                if annotation_summary.get("needs_review_clusters", 0)
                else f"{annotation_summary.get('low_confidence_clusters', 0)} low-confidence cluster(s) require review."
                if annotation_summary.get("low_confidence_clusters", 0)
                else f"{consensus_summary.get('needs_review_cells', 0)} cell(s) inherit annotation review status."
                if consensus_summary.get("needs_review_cells", 0)
                else "Consensus annotation present."
                if annotation_present
                else "Final annotation column is missing."
            ),
        ),
        _decision_row(
            step="condition_de",
            decision="prefer_pseudobulk" if not has_pseudobulk else "use_pseudobulk_results",
            recommended=condition_policy.get("recommended_primary_method", "sample_level_pseudobulk")
            if isinstance(condition_policy, Mapping)
            else "sample_level_pseudobulk",
            applied="pseudobulk_de" if has_pseudobulk else None,
            source="analysis_inference_policy",
            confidence=0.9 if has_pseudobulk else 0.75,
            affected_output="condition differential expression",
            inference_level="sample_level_condition_inference",
            biological_risk_note=condition_policy.get("risk_note", "")
            if isinstance(condition_policy, Mapping)
            else "",
            review_required=not has_pseudobulk,
            reason=(
                "Sample-level pseudobulk results detected."
                if has_pseudobulk
                else "No pseudobulk result detected; cell-level condition comparisons remain exploratory."
            ),
        ),
        _decision_row(
            step="cell_level_compare",
            decision="exploratory_only",
            recommended="do_not_use_for_publication_condition_inference",
            applied=condition_policy.get("cell_level_compare_policy", "exploratory_only")
            if isinstance(condition_policy, Mapping)
            else "exploratory_only",
            source="analysis_inference_policy",
            confidence=0.95,
            affected_output="compare_groups/compare_conditions",
            inference_level="exploratory_cell_level",
            biological_risk_note=(
                "Cell-level tests can rank hypotheses but cells are not biological replicates."
            ),
            review_required=False,
            reason="Cell-level comparison outputs are explicitly marked exploratory.",
        ),
        _decision_row(
            step="posthoc_qc_review",
            decision="review" if posthoc_review_required else "pass",
            recommended="inspect_flagged_clusters",
            applied=posthoc_review_required,
            source="posthoc_qc_review_summary",
            confidence=0.85,
            affected_output="cluster interpretation",
            inference_level="review_gate",
            biological_risk_note=str(posthoc_qc_summary.get("message", "")),
            review_required=posthoc_review_required,
            reason="; ".join(
                [
                    f"doublet_heavy={posthoc_qc_summary.get('n_doublet_heavy_clusters', 0)}",
                    f"high_mt={posthoc_qc_summary.get('n_high_mitochondrial_clusters', 0)}",
                    f"stress_high={posthoc_qc_summary.get('n_stress_high_clusters', 0)}",
                ]
            ),
        ),
        _decision_row(
            step="malignancy_interpretation",
            decision="review" if malignancy_review_required else "optional",
            recommended="run_in_tumor_workflow_when_relevant",
            applied=malignancy_enabled,
            source="malignancy_interpretation_summary",
            confidence=0.75 if malignancy_enabled else 0.5,
            affected_output="tumor interpretation",
            inference_level="tumor_context_evidence",
            biological_risk_note=(
                "Malignancy calls are integrated evidence for tumor workflows and need "
                "human review before biological claims."
            ),
            review_required=malignancy_review_required,
            reason=(
                "Malignancy interpretation enabled."
                if malignancy_enabled
                else "Not enabled in this analysis run."
            ),
        ),
    ]
    counts = pd.Series([row["decision"] for row in decisions]).value_counts().to_dict()
    return _json_safe(
        {
            "schema_version": "analysis_decision_summary_v1",
            "claim_boundary": inference_policy.get("claim_boundary"),
            "canonical_flow": output_contract.get("canonical_flow"),
            "primary_cluster_key": cluster_key,
            "primary_annotation_key": annotation_key,
            "primary_condition_de_method": (
                condition_policy.get("recommended_primary_method")
                if isinstance(condition_policy, Mapping)
                else "sample_level_pseudobulk"
            ),
            "decisions": decisions,
            "decision_counts": counts,
            "review_required_steps": [
                row["step"] for row in decisions if bool(row.get("review_required"))
            ],
        }
    )


def _split_marker_preview(value: Any, *, limit: int = 8) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [str(item) for item in value[:limit] if str(item)]
    if pd.isna(value):
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()][:limit]


def build_annotation_cluster_evidence_table(adata: AnnData, *, config: Any) -> list[dict[str, Any]]:
    """Build per-cluster annotation evidence rows from stored review outputs.

    The table is an audit artifact: it describes marker/reference/LLM support
    and conflicts, but it does not upgrade labels to formal truth.
    """
    annotation_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("annotation", {})
    review_table = (
        annotation_ns.get("annotation_review_table") if isinstance(annotation_ns, Mapping) else None
    )
    if not isinstance(review_table, pd.DataFrame) or review_table.empty:
        return []

    rows: list[dict[str, Any]] = []
    for _, row in review_table.iterrows():
        confidence = pd.to_numeric(
            pd.Series([row.get("annotation_confidence")]), errors="coerce"
        ).iloc[0]
        confidence_value = float(confidence) if pd.notna(confidence) else None
        needs_review = bool(row.get("needs_review", True))
        warnings = [item.strip() for item in str(row.get("warnings", "")).split(",") if item.strip()]
        conflicts = [item.strip() for item in str(row.get("conflicts", "")).split(",") if item.strip()]
        source_labels = {
            "marker": row.get("marker_label"),
            "reference": row.get("reference_label"),
            "llm": row.get("llm_label"),
        }
        source_confidences = {
            "marker": row.get("marker_confidence"),
            "reference": row.get("reference_confidence"),
            "llm": row.get("llm_confidence"),
        }
        usable_sources = [
            source
            for source, label in source_labels.items()
            if str(label) not in {"", "Unknown", "nan", "None"}
        ]
        contradiction_markers = sorted(
            {
                str(label)
                for label in source_labels.values()
                if str(label) not in {"", "Unknown", "nan", "None"}
                and str(label) != str(row.get("final_label"))
            }
        )
        evidence_status = (
            "conflicting"
            if conflicts or contradiction_markers
            else "supported"
            if usable_sources and not needs_review
            else "limited"
            if usable_sources
            else "missing"
        )
        if confidence_value is not None and confidence_value < 0.5:
            evidence_status = "low_confidence"

        rows.append(
            _json_safe(
                {
                    "cluster": str(row.get("cluster")),
                    "n_cells": int(row.get("n_cells")) if pd.notna(row.get("n_cells")) else 0,
                    "pct_cells": row.get("pct_cells"),
                    "predicted_label": row.get("final_label"),
                    "annotation_confidence": confidence_value,
                    "confidence_level": (
                        "high"
                        if confidence_value is not None and confidence_value >= 0.75
                        else "medium"
                        if confidence_value is not None and confidence_value >= 0.5
                        else "low"
                    ),
                    "claim_level": "evidence_consensus_not_formal_truth",
                    "evidence_status": evidence_status,
                    "positive_marker_support": _split_marker_preview(row.get("top_markers")),
                    "contradictory_labels": contradiction_markers,
                    "reference_model_label": row.get("reference_label"),
                    "reference_model_confidence": source_confidences.get("reference"),
                    "marker_label": row.get("marker_label"),
                    "marker_confidence": source_confidences.get("marker"),
                    "llm_label": row.get("llm_label"),
                    "llm_confidence": source_confidences.get("llm"),
                    "agreement_sources": usable_sources,
                    "decision": row.get("decision"),
                    "warnings": warnings,
                    "conflicts": conflicts,
                    "requires_manual_review": bool(
                        needs_review
                        or evidence_status in {"missing", "limited", "low_confidence", "conflicting"}
                    ),
                    "manual_review_recommendation": (
                        "Accept as provisional label after marker/reference spot-check."
                        if evidence_status == "supported"
                        else "Review marker support, contradictory labels, and cluster QC context before treating this label as final."
                    ),
                }
            )
        )
    return rows


def build_analysis_claim_level_summary(
    *,
    adata: AnnData,
    successful_steps: list[str],
    clustering_summary: Mapping[str, Any],
    annotation_summary: Mapping[str, Any],
    consensus_summary: Mapping[str, Any],
    decision_summary: Mapping[str, Any],
) -> dict[str, Any]:
    """Summarize what each analysis output is scientifically allowed to claim."""
    de_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("de", {})
    de_keys = sorted(de_ns.keys()) if isinstance(de_ns, Mapping) else []
    has_pseudobulk = any("pseudobulk" in key for key in de_keys)
    sccoda_ns = adata.uns.get("sclucid", {}).get("sccoda", {})
    has_sccoda = bool(sccoda_ns) if isinstance(sccoda_ns, Mapping) else False
    annotation_evidence_rows = annotation_summary.get("cluster_evidence_table", [])
    n_annotation_review_rows = (
        len(annotation_evidence_rows) if isinstance(annotation_evidence_rows, list) else 0
    )
    n_annotation_manual_review = (
        sum(bool(row.get("requires_manual_review")) for row in annotation_evidence_rows)
        if isinstance(annotation_evidence_rows, list)
        else int(annotation_summary.get("needs_review_clusters", 0) or 0)
    )

    outputs = [
        {
            "output": "clustering_resolution_recommendation",
            "claim_level": "heuristic_review_recommendation",
            "evidence_status": "available"
            if clustering_summary.get("resolution_review_available")
            else "not_run",
            "recommended_use": "Use to pick a practical clustering resolution for annotation review.",
            "not_allowed_claim": "Do not describe the selected resolution as mathematically optimal.",
            "requires_manual_review": bool(clustering_summary.get("review_required"))
            or not bool(clustering_summary.get("resolution_review_available")),
        },
        {
            "output": "cluster_marker_discovery",
            "claim_level": "exploratory_marker_screen",
            "evidence_status": "available" if "markers" in successful_steps else "not_run",
            "recommended_use": "Use for cluster identity hypotheses and marker visualization.",
            "not_allowed_claim": "Do not use cell-level marker tests as formal condition DE.",
            "requires_manual_review": False,
        },
        {
            "output": "annotation_consensus",
            "claim_level": "evidence_consensus_not_formal_truth",
            "evidence_status": "available"
            if consensus_summary.get("final_obs_present")
            else "missing",
            "recommended_use": "Use as provisional cell labels after marker/reference review.",
            "not_allowed_claim": "Do not treat automated labels as ground truth without review.",
            "requires_manual_review": bool(
                n_annotation_manual_review
                or annotation_summary.get("low_confidence_clusters", 0)
                or not consensus_summary.get("final_obs_present")
            ),
            "cluster_evidence_rows": n_annotation_review_rows,
            "clusters_requiring_manual_review": int(n_annotation_manual_review),
        },
        {
            "output": "condition_de",
            "claim_level": "formal_inference_when_sample_level_replicated"
            if has_pseudobulk
            else "not_formal_until_pseudobulk",
            "evidence_status": "pseudobulk_available" if has_pseudobulk else "missing_pseudobulk",
            "recommended_use": "Use sample-level pseudobulk with biological replicates for condition claims.",
            "not_allowed_claim": "Cell-level comparisons alone cannot support condition-level inference.",
            "requires_manual_review": not has_pseudobulk,
        },
        {
            "output": "cell_level_compare",
            "claim_level": "exploratory_hypothesis_generation",
            "evidence_status": "policy_defined",
            "recommended_use": "Use for marker ranking, visualization, and hypothesis generation.",
            "not_allowed_claim": "Do not report as publication-grade condition effect test.",
            "requires_manual_review": False,
        },
        {
            "output": "celltype_proportion",
            "claim_level": "sample_level_compositional_inference"
            if has_sccoda
            else "exploratory_or_sample_level_summary",
            "evidence_status": "sccoda_available" if has_sccoda else "no_compositional_model_detected",
            "recommended_use": "Prefer sample-level CLR/scCODA-style compositional analysis for condition comparisons.",
            "not_allowed_claim": "Raw cell fractions without sample-level compositional modeling are trend summaries.",
            "requires_manual_review": not has_sccoda,
        },
    ]
    review_required_outputs = [
        row["output"] for row in outputs if bool(row.get("requires_manual_review"))
    ]
    return _json_safe(
        {
            "schema_version": "analysis_claim_level_summary_v1",
            "global_claim_boundary": "heuristic_and_exploratory_until_evidence_review",
            "annotation_cluster_evidence_table": annotation_evidence_rows,
            "outputs": outputs,
            "review_required_outputs": review_required_outputs,
            "recommended_reporting_language": {
                "annotation": "automated evidence-supported labels requiring review",
                "cell_level_de": "exploratory marker screen",
                "condition_de": "sample-level pseudobulk inference when replicated",
                "proportion": "sample-level compositional inference when modeled",
            },
            "decision_summary_counts": decision_summary.get("decision_counts", {}),
        }
    )


def build_analysis_reviewer_table(decision_summary: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Return the single reviewer table for analysis decisions."""
    rows = []
    decision_rows = decision_summary.get("decisions", [])
    if isinstance(decision_rows, Mapping):
        decision_rows = decision_rows.values()
    for row in decision_rows:
        if not isinstance(row, Mapping):
            continue
        rows.append(
            {
                "item": row.get("step"),
                "recommended_value": row.get("recommended_value"),
                "applied_value": row.get("applied_value"),
                "source": row.get("source"),
                "confidence": row.get("confidence"),
                "affected_output": row.get("affected_output"),
                "analysis_decision": row.get("decision"),
                "inference_level": row.get("inference_level"),
                "biological_risk_note": row.get("biological_risk_note"),
                "review_required": row.get("review_required"),
                "review_reason": row.get("reason"),
            }
        )
    return _json_safe(rows)


def build_clustering_evidence_summary(adata: AnnData, cluster_key: str) -> dict[str, Any]:
    """Summarize clustering output and optional resolution-review evidence."""
    clustering_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("clustering", {})
    review_summary = (
        clustering_ns.get("clustering_review_summary", {})
        if isinstance(clustering_ns, Mapping)
        else {}
    )
    cluster_series = adata.obs[cluster_key].astype(str) if cluster_key in adata.obs else None
    counts = cluster_series.value_counts() if cluster_series is not None else pd.Series(dtype=int)
    review_required = []
    if cluster_key not in adata.obs:
        review_required.append(f"cluster_key_missing:{cluster_key}")
    if isinstance(review_summary, Mapping) and review_summary.get("review_required_clusters"):
        review_required.append("resolution_review_flagged_clusters")

    return _json_safe(
        {
            "cluster_key": cluster_key,
            "available": cluster_key in adata.obs,
            "n_clusters": int(counts.shape[0]) if not counts.empty else 0,
            "min_cluster_size": int(counts.min()) if not counts.empty else None,
            "median_cluster_size": float(counts.median()) if not counts.empty else None,
            "resolution_review_available": bool(review_summary),
            "recommended_resolution": review_summary.get("recommended_resolution")
            if isinstance(review_summary, Mapping)
            else None,
            "recommended_cluster_key": review_summary.get("recommended_cluster_key")
            if isinstance(review_summary, Mapping)
            else None,
            "recommendation_rationale": review_summary.get("rationale")
            if isinstance(review_summary, Mapping)
            else "",
            "recommendation_claim_level": "heuristic_review_recommendation",
            "recommendation_not_allowed_claim": (
                "Recommended clustering resolution is not an automatic optimum; it is a "
                "practical review heuristic based on marker and composition evidence."
            ),
            "review_required": review_required,
        }
    )


def build_annotation_evidence_summary(adata: AnnData, *, config: Any) -> dict[str, Any]:
    """Summarize available annotation evidence tables."""
    annotation_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("annotation", {})
    if not isinstance(annotation_ns, Mapping):
        annotation_ns = {}
    review_table = annotation_ns.get("annotation_review_table")
    marker_evidence = annotation_ns.get("marker_annotation_evidence")
    llm_bundle = annotation_ns.get("llm_annotation_bundle")
    annotation_key = (
        getattr(config.annotation, "key_added", "cell_type_auto")
        if getattr(config, "annotation", None) is not None
        else "cell_type_auto"
    )

    n_review_rows = int(review_table.shape[0]) if hasattr(review_table, "shape") else 0
    n_needs_review = 0
    n_low_confidence = 0
    min_confidence = None
    mean_confidence = None
    if isinstance(review_table, pd.DataFrame) and "needs_review" in review_table.columns:
        n_needs_review = int(review_table["needs_review"].fillna(True).astype(bool).sum())
    if isinstance(review_table, pd.DataFrame) and "annotation_confidence" in review_table.columns:
        confidence_values = pd.to_numeric(review_table["annotation_confidence"], errors="coerce")
        if confidence_values.notna().any():
            min_confidence = float(confidence_values.min())
            mean_confidence = float(confidence_values.mean())
            n_low_confidence = int((confidence_values < 0.5).fillna(True).sum())
    cluster_evidence_table = build_annotation_cluster_evidence_table(adata, config=config)
    n_cluster_evidence_review = sum(
        bool(row.get("requires_manual_review")) for row in cluster_evidence_table
    )

    return _json_safe(
        {
            "annotation_key": annotation_key,
            "annotation_obs_present": annotation_key in adata.obs,
            "review_table_available": isinstance(review_table, pd.DataFrame),
            "review_table_rows": n_review_rows,
            "needs_review_clusters": n_needs_review,
            "low_confidence_clusters": n_low_confidence,
            "min_annotation_confidence": min_confidence,
            "mean_annotation_confidence": mean_confidence,
            "confidence_review_threshold": 0.5,
            "marker_evidence_available": isinstance(marker_evidence, pd.DataFrame),
            "llm_bundle_available": isinstance(llm_bundle, Mapping),
            "evidence_methods": list(getattr(config, "annotation_methods", ()) or ()),
            "annotation_level": getattr(config, "annotation_level", None),
            "cluster_evidence_table": cluster_evidence_table,
            "cluster_evidence_rows": len(cluster_evidence_table),
            "cluster_evidence_review_required": int(n_cluster_evidence_review),
            "claim_level": "evidence_consensus_not_formal_truth",
            "recommended_use": "Use automated labels as provisional labels after marker/reference review.",
        }
    )


def build_annotation_consensus_summary(adata: AnnData, *, config: Any) -> dict[str, Any]:
    """Summarize final consensus labels applied to cells."""
    annotation_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("annotation", {})
    review_table = (
        annotation_ns.get("annotation_review_table") if isinstance(annotation_ns, Mapping) else None
    )
    final_key = (
        getattr(config.annotation, "key_added", "cell_type_auto")
        if getattr(config, "annotation", None) is not None
        else "cell_type_auto"
    )
    confidence_key = f"{final_key}_confidence"
    status_key = f"{final_key}_status"

    confidence = None
    if confidence_key in adata.obs:
        values = pd.to_numeric(adata.obs[confidence_key], errors="coerce")
        confidence = float(values.mean()) if values.notna().any() else None

    return _json_safe(
        {
            "final_key": final_key,
            "final_obs_present": final_key in adata.obs,
            "n_final_labels": int(adata.obs[final_key].nunique()) if final_key in adata.obs else 0,
            "mean_confidence": confidence,
            "status_key": status_key if status_key in adata.obs else None,
            "needs_review_cells": (
                int((adata.obs[status_key].astype(str) == "needs_review").sum())
                if status_key in adata.obs
                else 0
            ),
            "review_table_available": isinstance(review_table, pd.DataFrame),
        }
    )


def build_posthoc_qc_review_summary(
    adata: AnnData,
    *,
    cluster_key: str,
    doublet_flag_cols: tuple[str, ...] = (
        "predicted_doublet",
        "scrublet_predicted",
        "doubletdetection_predicted",
        "heuristic_predicted",
    ),
    doublet_fraction_threshold: float = 0.50,
    mt_col: str = "pct_counts_mt",
    mt_mean_threshold: float = 20.0,
    stress_score_cols: tuple[str, ...] = ("stress_score", "dissociation_stress_score"),
    stress_score_threshold: float = 0.50,
) -> dict[str, Any]:
    """Summarize analysis-time QC risks without automatically filtering cells."""
    if cluster_key not in adata.obs:
        return _json_safe(
            {
                "available": False,
                "cluster_key": cluster_key,
                "review_required": False,
                "message": f"Cluster key '{cluster_key}' is missing; post-hoc QC review skipped.",
            }
        )

    present_doublet_cols = tuple(col for col in doublet_flag_cols if col in adata.obs)
    present_stress_cols = tuple(col for col in stress_score_cols if col in adata.obs)
    has_mt = mt_col in adata.obs
    cluster_series = adata.obs[cluster_key].astype(str)
    rows: list[dict[str, Any]] = []

    for cluster, obs in adata.obs.groupby(cluster_series, observed=False):
        row: dict[str, Any] = {
            "cluster": str(cluster),
            "n_cells": int(obs.shape[0]),
            "reasons": [],
        }
        if present_doublet_cols:
            doublet_mask = obs.loc[:, list(present_doublet_cols)].fillna(False).astype(bool).any(axis=1)
            row["doublet_fraction"] = float(doublet_mask.mean())
            if row["doublet_fraction"] >= doublet_fraction_threshold:
                row["reasons"].append("doublet_heavy_cluster")
        else:
            row["doublet_fraction"] = None

        if has_mt:
            mt_values = pd.to_numeric(obs[mt_col], errors="coerce")
            row["mean_pct_counts_mt"] = float(mt_values.mean()) if mt_values.notna().any() else None
            if row["mean_pct_counts_mt"] is not None and row["mean_pct_counts_mt"] >= mt_mean_threshold:
                row["reasons"].append("high_mitochondrial_cluster")
        else:
            row["mean_pct_counts_mt"] = None

        if present_stress_cols:
            stress_values = obs.loc[:, list(present_stress_cols)].apply(
                pd.to_numeric, errors="coerce"
            )
            row["mean_stress_score"] = (
                float(stress_values.max(axis=1).mean())
                if stress_values.notna().any(axis=None)
                else None
            )
            if row["mean_stress_score"] is not None and row["mean_stress_score"] >= stress_score_threshold:
                row["reasons"].append("stress_high_cluster")
        else:
            row["mean_stress_score"] = None

        row["review_required"] = bool(row["reasons"])
        rows.append(row)

    doublet_heavy = [row["cluster"] for row in rows if "doublet_heavy_cluster" in row["reasons"]]
    high_mt = [row["cluster"] for row in rows if "high_mitochondrial_cluster" in row["reasons"]]
    stress_high = [row["cluster"] for row in rows if "stress_high_cluster" in row["reasons"]]
    review_required = bool(doublet_heavy or high_mt or stress_high)

    return _json_safe(
        {
            "available": True,
            "cluster_key": cluster_key,
            "n_clusters_reviewed": len(rows),
            "doublet_columns_used": list(present_doublet_cols),
            "stress_score_columns_used": list(present_stress_cols),
            "mt_column_used": mt_col if has_mt else None,
            "doublet_fraction_threshold": float(doublet_fraction_threshold),
            "mt_mean_threshold": float(mt_mean_threshold),
            "stress_score_threshold": float(stress_score_threshold),
            "n_doublet_heavy_clusters": len(doublet_heavy),
            "n_high_mitochondrial_clusters": len(high_mt),
            "n_stress_high_clusters": len(stress_high),
            "doublet_heavy_clusters": doublet_heavy,
            "high_mitochondrial_clusters": high_mt,
            "stress_high_clusters": stress_high,
            "review_required": review_required,
            "table": rows,
            "message": (
                "Post-hoc QC review found clusters that should be inspected before final "
                "annotation or filtering."
                if review_required
                else "No cluster-level doublet, mitochondrial, or stress flags exceeded review thresholds."
            ),
        }
    )


def build_malignancy_interpretation_summary(adata: AnnData, *, config: Any) -> dict[str, Any]:
    """Summarize optional malignancy interpretation evidence."""
    malignancy_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("malignancy", {})
    stored = (
        malignancy_ns.get("malignancy_interpretation_summary", {})
        if isinstance(malignancy_ns, Mapping)
        else {}
    )
    call_key = getattr(config, "malignancy_key_added", "malignancy_call")
    score_key = getattr(config, "malignancy_score_key", "malignancy_interpretation_score")
    if isinstance(stored, Mapping) and stored:
        summary = dict(stored)
    else:
        summary = {
            "available": call_key in adata.obs,
            "call_key": call_key,
            "score_key": score_key,
            "review_required": False,
            "evidence_sources": [],
        }
    if call_key in adata.obs:
        calls = adata.obs[call_key].astype(str)
        n_cells = int(calls.shape[0])
        suspect_or_malignant = int(
            ((calls == "malignant") | (calls == "suspect_malignant")).sum()
        )
        summary.update(
            {
                "available": True,
                "n_malignant": int((calls == "malignant").sum()),
                "n_suspect_malignant": int((calls == "suspect_malignant").sum()),
                "n_non_malignant": int((calls == "non_malignant").sum()),
                "n_unresolved": int((calls == "unresolved").sum()),
                "n_cells_evaluated": n_cells,
                "malignant_fraction": (
                    float((calls == "malignant").sum() / n_cells) if n_cells else 0.0
                ),
                "suspect_or_malignant_fraction": (
                    float(suspect_or_malignant / n_cells) if n_cells else 0.0
                ),
                "tumor_purity_estimate": (
                    float(suspect_or_malignant / n_cells) if n_cells else 0.0
                ),
                "low_tumor_purity_threshold": float(
                    summary.get("low_tumor_purity_threshold", 0.10)
                ),
            }
        )
        summary["low_tumor_purity_warning"] = bool(
            summary["suspect_or_malignant_fraction"] < summary["low_tumor_purity_threshold"]
        )
    if score_key in adata.obs:
        scores = pd.to_numeric(adata.obs[score_key], errors="coerce")
        summary["mean_score"] = float(scores.mean()) if scores.notna().any() else None
    summary["enabled"] = bool(getattr(config, "run_malignancy_interpretation", False))
    return _json_safe(summary)


def build_analysis_readiness_assessment(
    *,
    adata: AnnData,
    successful_steps: list[str],
    cluster_key: str,
    preprocess_context: Mapping[str, Any],
    clustering_summary: Mapping[str, Any],
    annotation_summary: Mapping[str, Any],
    consensus_summary: Mapping[str, Any],
    posthoc_qc_summary: Mapping[str, Any],
    malignancy_summary: Mapping[str, Any],
    claim_level_summary: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Assess whether analysis outputs are ready for downstream interpretation."""
    score = 1.0
    reasons: list[str] = []
    blockers: list[str] = []

    if not preprocess_context.get("pca_present"):
        blockers.append("preprocess_input_context.pca_present=False")
        score -= 0.35
    if cluster_key not in adata.obs:
        blockers.append(f"cluster_key_missing:{cluster_key}")
        score -= 0.35
    if "markers" in successful_steps and "rank_genes_groups" not in adata.uns:
        reasons.append("markers_step_ran_but_rank_genes_groups_missing")
        score -= 0.15
    if annotation_summary.get("needs_review_clusters", 0):
        reasons.append("annotation_review_required_clusters_present")
        score -= 0.15
    if annotation_summary.get("low_confidence_clusters", 0):
        reasons.append("annotation_low_confidence_clusters_present")
        score -= 0.10
    if consensus_summary.get("needs_review_cells", 0):
        reasons.append("annotation_consensus_needs_review_cells_present")
        score -= 0.10
    if "annotation" in successful_steps and not consensus_summary.get("final_obs_present"):
        reasons.append("annotation_consensus_not_applied")
        score -= 0.15
    if posthoc_qc_summary.get("review_required"):
        reasons.append("posthoc_qc_review_required")
        score -= 0.05
    if "malignancy_interpretation" in successful_steps:
        if not malignancy_summary.get("available"):
            reasons.append("malignancy_interpretation_missing")
            score -= 0.10
        elif malignancy_summary.get("review_required"):
            reasons.append("malignancy_interpretation_review_required")
            score -= 0.05
    if claim_level_summary and claim_level_summary.get("review_required_outputs"):
        reasons.append("claim_level_outputs_require_review")
        score -= 0.05

    score = float(max(0.0, min(1.0, score)))
    if blockers:
        status = "blocked"
    elif reasons:
        status = "review_required"
    else:
        status = "ready"

    return _json_safe(
        {
            "status": status,
            "score": score,
            "blockers": blockers,
            "review_reasons": reasons,
            "ready_for_downstream": status in {"ready", "review_required"},
        }
    )


def build_analysis_review_action_items(
    *,
    readiness: Mapping[str, Any],
    clustering_summary: Mapping[str, Any],
    annotation_summary: Mapping[str, Any],
    consensus_summary: Mapping[str, Any],
    posthoc_qc_summary: Mapping[str, Any],
    malignancy_summary: Mapping[str, Any],
    claim_level_summary: Mapping[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Create human-review action items from analysis evidence."""
    actions: list[ReviewAction] = []
    if readiness.get("blockers"):
        actions.append(
            ReviewAction(
                priority="blocking",
                action=(
                    "Resolve missing preprocessing or clustering outputs before downstream "
                    "analysis."
                ),
                rationale="Analysis readiness has blocking failures.",
                evidence_keys=["analysis_readiness", "preprocess_input_context"],
            )
        )
    if clustering_summary.get("review_required"):
        actions.append(
            ReviewAction(
                priority="review",
                action="Review clustering resolution and flagged clusters before final annotation.",
                rationale="Resolution evidence marked one or more clustering concerns.",
                evidence_keys=["clustering_evidence_summary"],
            )
        )
    if annotation_summary.get("needs_review_clusters", 0):
        actions.append(
            ReviewAction(
                priority="review",
                action="Manually inspect annotation clusters marked as needing review.",
                rationale="Consensus annotation found weak or conflicting evidence.",
                evidence_keys=["annotation_evidence_summary", "annotation_consensus_summary"],
            )
        )
    elif annotation_summary.get("low_confidence_clusters", 0):
        actions.append(
            ReviewAction(
                priority="review",
                action="Review low-confidence annotation clusters before treating labels as final.",
                rationale=(
                    "One or more clusters have annotation_confidence below the review threshold "
                    "even if a final label was assigned."
                ),
                evidence_keys=["annotation_evidence_summary", "annotation_consensus_summary"],
            )
        )
    if not consensus_summary.get("final_obs_present"):
        actions.append(
            ReviewAction(
                priority="optional",
                action="Apply consensus labels when a final cell-type column is required.",
                rationale=(
                    "Annotation evidence exists but no final consensus obs column was detected."
                ),
                evidence_keys=["annotation_consensus_summary"],
            )
        )
    if posthoc_qc_summary.get("review_required"):
        actions.append(
            ReviewAction(
                priority="review",
                action=(
                    "Inspect doublet-heavy, high-mitochondrial, or stress-high clusters before "
                    "final annotation or downstream tumor interpretation."
                ),
                rationale=(
                    "Analysis-time QC review found cluster-level technical-risk patterns; "
                    "these should usually trigger manual review, not automatic deletion."
                ),
                evidence_keys=["posthoc_qc_review_summary"],
            )
        )
    if malignancy_summary.get("enabled") and malignancy_summary.get("review_required"):
        actions.append(
            ReviewAction(
                priority="review",
                action="Review malignant/suspect/unresolved calls before tumor downstream analysis.",
                rationale=(
                    "Malignancy interpretation combines annotation, marker, CNV, and signature "
                    "evidence and found calls that require human confirmation."
                ),
                evidence_keys=["malignancy_interpretation_summary"],
            )
        )
    if claim_level_summary and claim_level_summary.get("review_required_outputs"):
        actions.append(
            ReviewAction(
                priority="review",
                action="Review analysis claim levels before reporting automated results.",
                rationale=(
                    "Some outputs are heuristic, exploratory, or missing formal sample-level "
                    "evidence; report them with the claim level recorded in the audit summary."
                ),
                evidence_keys=["analysis_claim_level_summary"],
            )
        )
    return [model_to_dict(action) for action in actions]


def build_analysis_evidence_bundle(summary: Mapping[str, Any]) -> dict[str, Any]:
    """Build the common evidence bundle for analysis."""
    readiness = summary.get("analysis_readiness", {}) if isinstance(summary, Mapping) else {}
    clustering = (
        summary.get("clustering_evidence_summary", {}) if isinstance(summary, Mapping) else {}
    )
    annotation = (
        summary.get("annotation_evidence_summary", {}) if isinstance(summary, Mapping) else {}
    )
    malignancy = (
        summary.get("malignancy_interpretation_summary", {})
        if isinstance(summary, Mapping)
        else {}
    )
    posthoc_qc = (
        summary.get("posthoc_qc_review_summary", {}) if isinstance(summary, Mapping) else {}
    )
    inference_policy = (
        summary.get("analysis_inference_policy", {}) if isinstance(summary, Mapping) else {}
    )
    claim_level_summary = (
        summary.get("analysis_claim_level_summary", {}) if isinstance(summary, Mapping) else {}
    )
    output_contract = (
        summary.get("analysis_output_contract", {}) if isinstance(summary, Mapping) else {}
    )
    decision_summary = (
        summary.get("analysis_decision_summary", {}) if isinstance(summary, Mapping) else {}
    )

    evidence_chain = [
        EvidenceItem(
            source="context",
            name="analysis_inference_policy",
            value=inference_policy,
            confidence=None,
            rationale=(
                "Defines conservative claim boundaries for clustering, marker discovery, "
                "condition DE, and cell-level comparisons."
            ),
            related_keys=["analysis_inference_policy"],
        ),
        EvidenceItem(
            source="contract",
            name="analysis_claim_level_summary",
            value=claim_level_summary,
            confidence=None,
            rationale=(
                "Separates heuristic, exploratory, evidence-consensus, and formal-inference "
                "claims for analysis outputs."
            ),
            related_keys=["analysis_claim_level_summary"],
        ),
        EvidenceItem(
            source="contract",
            name="analysis_output_contract",
            value=output_contract,
            confidence=None,
            rationale="Defines stable analysis output slots and inference semantics.",
            related_keys=["analysis_output_contract"],
        ),
        EvidenceItem(
            source="recommendation",
            name="analysis_decision_summary",
            value=decision_summary,
            confidence=None,
            rationale="Condenses analysis choices into reviewer-facing decisions.",
            related_keys=["analysis_decision_summary", "analysis_reviewer_table"],
        ),
        EvidenceItem(
            source="context",
            name="preprocess_input_context",
            value=summary.get("preprocess_input_context", {}),
            confidence=None,
            rationale="Analysis depends on preprocessing outputs such as PCA and neighbors.",
            related_keys=["preprocess_input_context"],
        ),
        EvidenceItem(
            source="metric",
            name="clustering_evidence",
            value=clustering,
            confidence=None,
            rationale="Cluster interpretability gates marker discovery and annotation.",
            related_keys=["clustering_evidence_summary"],
        ),
        EvidenceItem(
            source="output_health",
            name="annotation_evidence",
            value=annotation,
            confidence=None,
            rationale="Annotation evidence summarizes marker/reference/data-driven agreement.",
            related_keys=["annotation_evidence_summary", "annotation_consensus_summary"],
        ),
        EvidenceItem(
            source="output_health",
            name="malignancy_interpretation",
            value=malignancy,
            confidence=None,
            rationale=(
                "Malignancy interpretation summarizes optional tumor-context evidence for "
                "downstream tumor analysis."
            ),
            related_keys=["malignancy_interpretation_summary"],
        ),
        EvidenceItem(
            source="output_health",
            name="posthoc_qc_review",
            value=posthoc_qc,
            confidence=None,
            rationale=(
                "Cluster-level doublet, mitochondrial, and stress patterns are reviewed after "
                "analysis to avoid over-filtering biologically plausible tumor states."
            ),
            related_keys=["posthoc_qc_review_summary"],
        ),
    ]
    bundle = EvidenceBundle(
        module="analysis",
        stage="run_standard_analysis",
        status=str(readiness.get("status", "unknown")),
        confidence=readiness.get("score") if isinstance(readiness, Mapping) else None,
        context={"workflow_name": summary.get("workflow_name", "standard")},
        evidence_chain=evidence_chain,
        action_items=[
            ReviewAction(**item)
            for item in summary.get("review_action_items", [])
            if isinstance(item, Mapping)
        ],
        related_review_keys=sorted(ANALYSIS_REQUIRED_REVIEW_SECTIONS),
    )
    return _json_safe(model_to_dict(bundle))


def build_analysis_module_maturity_assessment(summary: Mapping[str, Any]) -> dict[str, Any]:
    """Assess whether an analysis review summary satisfies the module contract."""
    payload = _review_payload(summary)
    required_sections = set(ANALYSIS_REQUIRED_REVIEW_SECTIONS)
    required_sections.discard("module_maturity")
    missing = sorted(required_sections - set(payload.keys()))
    issues = [f"missing_required_section:{key}" for key in missing]
    review_required = []

    readiness = payload.get("analysis_readiness", {}) if isinstance(payload, Mapping) else {}
    if readiness.get("status") == "blocked":
        issues.extend(readiness.get("blockers", []))
    elif readiness.get("status") == "review_required":
        review_required.extend(readiness.get("review_reasons", []))
    output_contract = payload.get("analysis_output_contract", {})
    if isinstance(output_contract, Mapping) and output_contract.get("review_required"):
        review_required.append("analysis_output_contract.review_required=True")
    decision_summary = payload.get("analysis_decision_summary", {})
    if isinstance(decision_summary, Mapping) and decision_summary.get("review_required_steps"):
        review_required.extend(
            f"analysis_decision_review:{step}"
            for step in decision_summary.get("review_required_steps", [])
        )
    claim_level_summary = payload.get("analysis_claim_level_summary", {})
    if isinstance(claim_level_summary, Mapping) and claim_level_summary.get(
        "review_required_outputs"
    ):
        review_required.extend(
            f"analysis_claim_level_review:{output}"
            for output in claim_level_summary.get("review_required_outputs", [])
        )

    if issues:
        status = "incomplete"
    elif review_required:
        status = "review_required"
    else:
        status = "complete"

    return _json_safe(
        {
            "schema_version": ANALYSIS_MODULE_MATURITY_SCHEMA_VERSION,
            "module": "analysis",
            "status": status,
            "issues": issues,
            "review_required": review_required,
            "contract": get_analysis_module_contract(),
            "summary": (
                "Analysis review summary satisfies the benchmark module contract."
                if status == "complete"
                else "Analysis review summary is present but requires review."
                if status == "review_required"
                else "Analysis review summary does not satisfy the benchmark module contract."
            ),
        }
    )


def summarize_analysis_review_summary(summary: Mapping[str, Any]) -> dict[str, Any]:
    """Return a compact user-facing summary of analysis output."""
    payload = _review_payload(summary)
    readiness = payload.get("analysis_readiness", {}) if isinstance(payload, Mapping) else {}
    maturity = payload.get("module_maturity", {}) if isinstance(payload, Mapping) else {}
    inference_policy = (
        payload.get("analysis_inference_policy", {}) if isinstance(payload, Mapping) else {}
    )
    claim_level_summary = (
        payload.get("analysis_claim_level_summary", {}) if isinstance(payload, Mapping) else {}
    )
    output_contract = (
        payload.get("analysis_output_contract", {}) if isinstance(payload, Mapping) else {}
    )
    decision_summary = (
        payload.get("analysis_decision_summary", {}) if isinstance(payload, Mapping) else {}
    )
    clustering = (
        payload.get("clustering_evidence_summary", {}) if isinstance(payload, Mapping) else {}
    )
    annotation = (
        payload.get("annotation_evidence_summary", {}) if isinstance(payload, Mapping) else {}
    )
    consensus = (
        payload.get("annotation_consensus_summary", {}) if isinstance(payload, Mapping) else {}
    )
    malignancy = (
        payload.get("malignancy_interpretation_summary", {})
        if isinstance(payload, Mapping)
        else {}
    )
    posthoc_qc = (
        payload.get("posthoc_qc_review_summary", {}) if isinstance(payload, Mapping) else {}
    )
    return _json_safe(
        {
            "module": "analysis",
            "maturity_status": maturity.get("status"),
            "readiness_status": readiness.get("status"),
            "readiness_score": readiness.get("score"),
            "claim_boundary": inference_policy.get("claim_boundary"),
            "global_claim_boundary": claim_level_summary.get("global_claim_boundary")
            if isinstance(claim_level_summary, Mapping)
            else None,
            "claim_level_review_required_outputs": (
                claim_level_summary.get("review_required_outputs", [])
                if isinstance(claim_level_summary, Mapping)
                else []
            ),
            "condition_de_primary_method": (
                inference_policy.get("condition_de", {}).get("recommended_primary_method")
                if isinstance(inference_policy.get("condition_de"), Mapping)
                else None
            ),
            "cell_level_compare_policy": (
                inference_policy.get("condition_de", {}).get("cell_level_compare_policy")
                if isinstance(inference_policy.get("condition_de"), Mapping)
                else None
            ),
            "canonical_analysis_flow": output_contract.get("canonical_flow"),
            "analysis_decision_counts": decision_summary.get("decision_counts", {}),
            "analysis_review_required_steps": decision_summary.get("review_required_steps", []),
            "primary_annotation_key": decision_summary.get("primary_annotation_key"),
            "cluster_key": clustering.get("cluster_key"),
            "n_clusters": clustering.get("n_clusters"),
            "recommended_resolution": clustering.get("recommended_resolution"),
            "recommended_resolution_claim_level": clustering.get(
                "recommendation_claim_level"
            ),
            "annotation_key": annotation.get("annotation_key"),
            "review_table_rows": annotation.get("review_table_rows"),
            "needs_review_clusters": annotation.get("needs_review_clusters"),
            "annotation_cluster_evidence_rows": annotation.get("cluster_evidence_rows"),
            "annotation_cluster_evidence_review_required": annotation.get(
                "cluster_evidence_review_required"
            ),
            "final_key": consensus.get("final_key"),
            "n_final_labels": consensus.get("n_final_labels"),
            "mean_confidence": consensus.get("mean_confidence"),
            "posthoc_qc_review_required": posthoc_qc.get("review_required"),
            "n_doublet_heavy_clusters": posthoc_qc.get("n_doublet_heavy_clusters"),
            "n_high_mitochondrial_clusters": posthoc_qc.get("n_high_mitochondrial_clusters"),
            "n_stress_high_clusters": posthoc_qc.get("n_stress_high_clusters"),
            "malignancy_enabled": malignancy.get("enabled"),
            "n_malignant": malignancy.get("n_malignant"),
            "n_suspect_malignant": malignancy.get("n_suspect_malignant"),
            "suspect_or_malignant_fraction": malignancy.get("suspect_or_malignant_fraction"),
        }
    )


def validate_analysis_review_summary(
    summary: Mapping[str, Any],
    *,
    raise_on_error: bool = False,
) -> list[str]:
    """Validate analysis-specific review-summary sections."""
    payload = _review_payload(summary)
    errors: list[str] = []
    missing = sorted(ANALYSIS_REQUIRED_REVIEW_SECTIONS - set(payload.keys()))
    if missing:
        errors.append(f"Analysis review summary missing required sections: {missing}")
    bundle = payload.get("evidence_bundle")
    if not isinstance(bundle, Mapping):
        errors.append("Analysis review summary field 'evidence_bundle' must be a mapping.")
    elif bundle.get("module") != "analysis":
        errors.append("Analysis evidence_bundle.module must be 'analysis'.")
    maturity = payload.get("module_maturity")
    if not isinstance(maturity, Mapping):
        errors.append("Analysis review summary field 'module_maturity' must be a mapping.")
    elif maturity.get("module") != "analysis":
        errors.append("Analysis module_maturity.module must be 'analysis'.")
    output_contract = payload.get("analysis_output_contract")
    if not isinstance(output_contract, Mapping):
        errors.append("Analysis review summary field 'analysis_output_contract' must be a mapping.")
    elif output_contract.get("schema_version") != "analysis_output_contract_v1":
        errors.append("Analysis output contract has an unsupported schema version.")
    claim_level_summary = payload.get("analysis_claim_level_summary")
    if not isinstance(claim_level_summary, Mapping):
        errors.append(
            "Analysis review summary field 'analysis_claim_level_summary' must be a mapping."
        )
    elif claim_level_summary.get("schema_version") != "analysis_claim_level_summary_v1":
        errors.append("Analysis claim level summary has an unsupported schema version.")
    decision_summary = payload.get("analysis_decision_summary")
    if not isinstance(decision_summary, Mapping):
        errors.append("Analysis review summary field 'analysis_decision_summary' must be a mapping.")
    elif decision_summary.get("schema_version") != "analysis_decision_summary_v1":
        errors.append("Analysis decision summary has an unsupported schema version.")
    reviewer_table = payload.get("analysis_reviewer_table")
    if isinstance(reviewer_table, Mapping):
        reviewer_rows = list(reviewer_table.values())
    elif isinstance(reviewer_table, list):
        reviewer_rows = reviewer_table
    else:
        reviewer_rows = []
    if not reviewer_rows:
        errors.append(
            "Analysis review summary field 'analysis_reviewer_table' must be a non-empty list or mapping."
        )
    else:
        required_columns = {
            "item",
            "recommended_value",
            "applied_value",
            "source",
            "confidence",
            "affected_output",
            "analysis_decision",
            "inference_level",
            "biological_risk_note",
            "review_required",
        }
        for idx, row in enumerate(reviewer_rows):
            if not isinstance(row, Mapping):
                errors.append(f"Analysis reviewer table row {idx} must be a mapping.")
                continue
            missing_columns = sorted(required_columns - set(row.keys()))
            if missing_columns:
                errors.append(f"Analysis reviewer table row {idx} missing columns: {missing_columns}")
    if errors and raise_on_error:
        raise ValueError("; ".join(errors))
    return errors


def validate_analysis_module_completeness(
    adata: AnnData,
    *,
    require_ready: bool = False,
    raise_on_error: bool = False,
) -> dict[str, Any]:
    """Validate that an AnnData object contains a benchmark-grade analysis result."""
    issues: list[str] = []
    warnings: list[str] = []
    analysis_ns = adata.uns.get("sclucid", {}).get("analysis", {})
    if not isinstance(analysis_ns, Mapping):
        issues.append('Missing or invalid adata.uns["sclucid"]["analysis"] namespace.')
        analysis_ns = {}

    review_summary = analysis_ns.get("review_summary")
    payload = _review_payload(review_summary) if isinstance(review_summary, Mapping) else {}
    if not payload:
        issues.append('Missing adata.uns["sclucid"]["analysis"]["review_summary"].')
        maturity = build_analysis_module_maturity_assessment({})
    else:
        issues.extend(validate_analysis_review_summary(payload))
        maturity = build_analysis_module_maturity_assessment(payload)
        if maturity.get("status") == "incomplete":
            issues.extend(maturity.get("issues", []))
        elif maturity.get("status") == "review_required":
            warnings.extend(maturity.get("review_required", []))

    if not any(key in adata.obs for key in ("leiden_clusters", "leiden")):
        issues.append("Missing canonical analysis cluster column: 'leiden_clusters' or 'leiden'.")

    readiness = payload.get("analysis_readiness", {}) if isinstance(payload, Mapping) else {}
    if require_ready and readiness.get("status") != "ready":
        issues.append(f"Analysis readiness is {readiness.get('status')!r}, expected 'ready'.")

    result = {
        "schema_version": ANALYSIS_MODULE_MATURITY_SCHEMA_VERSION,
        "module": "analysis",
        "valid": len(issues) == 0,
        "status": "valid" if not issues else "invalid",
        "issues": list(dict.fromkeys(str(item) for item in issues)),
        "warnings": list(dict.fromkeys(str(item) for item in warnings)),
        "maturity": maturity,
        "summary": summarize_analysis_review_summary(payload) if payload else {},
    }
    if result["issues"] and raise_on_error:
        raise ValueError("; ".join(result["issues"]))
    return _json_safe(result)


def _json_safe(value: Any) -> Any:
    return sanitize_for_hdf5(value)


__all__ = [
    "ANALYSIS_MODULE_MATURITY_SCHEMA_VERSION",
    "ANALYSIS_REQUIRED_REVIEW_SECTIONS",
    "ANALYSIS_TRACE_SCHEMA_VERSION",
    "build_analysis_decision_summary",
    "build_analysis_inference_policy",
    "build_analysis_module_maturity_assessment",
    "build_analysis_output_contract",
    "build_analysis_reviewer_table",
    "build_malignancy_interpretation_summary",
    "build_posthoc_qc_review_summary",
    "enrich_analysis_review_summary",
    "get_analysis_module_contract",
    "summarize_analysis_review_summary",
    "validate_analysis_module_completeness",
    "validate_analysis_review_summary",
]
