"""Project planning and cross-stage decision review for scLucid.

This module is deliberately conservative.  It does not claim to select a
scientifically optimal workflow from metadata alone.  Instead it exposes the
minimum study-design assumptions required for a defensible first pass and
normalizes the evidence already produced by QC, preprocessing, analysis, and
tumor workflows into one action-oriented review surface.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Literal, Optional

from anndata import AnnData
from pydantic import BaseModel, ConfigDict, Field

from .utils.context import AnalysisContext, ProjectContext, infer_analysis_context
from .utils.contracts import SCLUCID_ROOT, UnsKeys, ensure_sclucid_namespace
from .utils.sanitize import sanitize_for_hdf5

DecisionStatus = Literal["BLOCKED", "REVIEW", "READY"]
StageStatus = Literal["BLOCKED", "REVIEW", "READY", "NOT_RUN"]

RUN_REVIEW_SCHEMA_VERSION = "1.0"
ANALYSIS_PLAN_SCHEMA_VERSION = "1.0"

_STAGE_ORDER = ("qc", "preprocess", "analysis", "tumor")
_RERUN_SCOPE = {
    "qc": "rerun QC, then preprocessing and downstream stages",
    "preprocess": "rerun preprocessing, then analysis and tumor stages",
    "analysis": "rerun analysis and any dependent tumor interpretation",
    "tumor": "rerun only the affected tumor interpretation stage",
}


class _DecisionModel(BaseModel):
    """Small JSON-friendly model base used by the product decision layer."""

    model_config = ConfigDict(extra="ignore", arbitrary_types_allowed=True)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible dictionary."""
        return self.model_dump(
            mode="json",
            fallback=lambda value: value.item() if hasattr(value, "item") else str(value),
        )


class DecisionCard(_DecisionModel):
    """One explicit scientific or workflow decision requiring user visibility."""

    stage: str
    status: DecisionStatus
    decision: str
    recommended: Any = None
    applied: Any = None
    reason: str = ""
    evidence: list[str] = Field(default_factory=list)
    next_action: str = ""
    rerun_scope: str = ""
    priority: str = "review"
    source: str = "project_context"


class AnalysisPlan(_DecisionModel):
    """Conservative pre-run plan derived from project context and AnnData structure."""

    schema_version: str = ANALYSIS_PLAN_SCHEMA_VERSION
    profile: str
    status: DecisionStatus
    stages: list[str]
    context: ProjectContext
    decisions: list[DecisionCard] = Field(default_factory=list)
    required_metadata: list[str] = Field(default_factory=list)
    blockers: list[str] = Field(default_factory=list)
    assumptions: list[str] = Field(default_factory=list)
    highest_risk_failure_modes: list[str] = Field(default_factory=list)
    smallest_next_step: str = "Review the plan before running the workflow."

    @property
    def ready_to_run(self) -> bool:
        """Whether the plan has no structural blocker."""
        return self.status != "BLOCKED"

    def to_frame(self):
        """Return decision cards as a compact pandas DataFrame."""
        import pandas as pd

        return pd.DataFrame([item.to_dict() for item in self.decisions])


class StageReview(_DecisionModel):
    """Normalized readiness for one workflow stage."""

    stage: str
    status: StageStatus
    raw_status: str = "unknown"
    confidence: Optional[float] = None
    reasons: list[str] = Field(default_factory=list)
    claim_boundary: Optional[str] = None
    next_action: str = ""
    rerun_scope: str = ""


class RunReview(_DecisionModel):
    """Unified cross-stage run review returned by :func:`review_run`."""

    schema_version: str = RUN_REVIEW_SCHEMA_VERSION
    overall_status: StageStatus
    context: dict[str, Any] = Field(default_factory=dict)
    stages: list[StageReview] = Field(default_factory=list)
    items: list[DecisionCard] = Field(default_factory=list)
    next_actions: list[str] = Field(default_factory=list)
    assumptions: list[str] = Field(default_factory=list)
    highest_risk_failure_modes: list[str] = Field(default_factory=list)

    @property
    def ready_for_final_interpretation(self) -> bool:
        """Whether every recorded stage is ready without unresolved review."""
        analysis_ready = any(
            stage.stage == "analysis" and stage.status == "READY" for stage in self.stages
        )
        return self.overall_status == "READY" and analysis_ready

    def to_frame(self):
        """Return all stage, action, and parameter decisions as a DataFrame."""
        import pandas as pd

        columns = [
            "stage",
            "status",
            "decision",
            "recommended",
            "applied",
            "reason",
            "evidence",
            "next_action",
            "rerun_scope",
            "priority",
            "source",
        ]
        records = [item.to_dict() for item in self.items]
        return pd.DataFrame(records, columns=columns)

    def show_next_actions(self) -> list[str]:
        """Return the prioritized, deduplicated actions the user should take next."""
        return list(self.next_actions)


def _records(value: Any) -> list[Mapping[str, Any]]:
    """Restore list-like records from native or HDF5-sanitized mappings."""
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [item for item in value if isinstance(item, Mapping)]
    if isinstance(value, Mapping) and value and all(str(key).isdigit() for key in value):
        keys = sorted(value, key=lambda key: int(str(key)))
        return [value[key] for key in keys if isinstance(value[key], Mapping)]
    return []


def _strings(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, Mapping) and value and all(str(key).isdigit() for key in value):
        value = [value[key] for key in sorted(value, key=lambda key: int(str(key)))]
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes, bytearray)):
        value = value.tolist()
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [str(item) for item in value if item not in (None, "")]
    return [str(value)]


def _review_payload(summary: Any) -> Mapping[str, Any]:
    if not isinstance(summary, Mapping):
        return {}
    data = summary.get("data")
    return data if isinstance(data, Mapping) else summary


def _worst_status(statuses: Sequence[str], *, empty: StageStatus = "NOT_RUN") -> StageStatus:
    if "BLOCKED" in statuses:
        return "BLOCKED"
    if "REVIEW" in statuses:
        return "REVIEW"
    if "READY" in statuses:
        return "READY"
    return empty


def _normalize_status(raw_status: Any, *, present: bool) -> StageStatus:
    if not present:
        return "NOT_RUN"
    token = str(raw_status or "unknown").strip().lower().replace("-", "_")
    if token in {
        "blocked",
        "failed",
        "error",
        "invalid",
        "incomplete",
        "not_ready",
        "unavailable",
    }:
        return "BLOCKED"
    if token in {
        "review",
        "review_required",
        "needs_review",
        "degraded",
        "partial",
        "warning",
        "heuristic",
        "exploratory",
        "unknown",
    }:
        return "REVIEW"
    if token in {"ready", "complete", "completed", "valid", "ok", "pass", "passed"}:
        return "READY"
    return "REVIEW"


def _is_comparative_objective(objective: Optional[str]) -> bool:
    token = str(objective or "").lower()
    return any(
        marker in token
        for marker in (
            "compare",
            "comparison",
            "differential",
            "treatment",
            "response",
            "paired",
            "condition",
            "比较",
            "差异",
            "治疗",
            "疗效",
        )
    )


def _resolve_profile(context: AnalysisContext, requested: str) -> str:
    if requested and requested != "auto":
        return requested
    if context.dataset_type == "tumor_tissue" and context.is_multi_sample:
        return "multi_sample_tumor"
    if context.dataset_type == "tumor_tissue":
        return "tumor_conservative"
    if context.dataset_type == "cell_line":
        return "cell_line"
    if _is_comparative_objective(context.study_objective):
        return "treatment_response"
    return "baseline"


def _condition_replicates(adata: AnnData, context: AnalysisContext) -> dict[str, int]:
    condition_key = context.condition_key
    unit_key = context.experimental_unit_key or context.sample_key
    if not condition_key or not unit_key:
        return {}
    if condition_key not in adata.obs.columns or unit_key not in adata.obs.columns:
        return {}
    table = adata.obs[[condition_key, unit_key]].dropna().drop_duplicates()
    return {
        str(condition): int(frame[unit_key].nunique())
        for condition, frame in table.groupby(condition_key, observed=True)
    }


def _batch_condition_confounded(adata: AnnData, context: AnalysisContext) -> bool:
    if not context.batch_key or not context.condition_key:
        return False
    if context.batch_key not in adata.obs.columns or context.condition_key not in adata.obs.columns:
        return False
    table = adata.obs[[context.batch_key, context.condition_key]].dropna().drop_duplicates()
    if table.empty:
        return False
    if table[context.batch_key].nunique() < 2 or table[context.condition_key].nunique() < 2:
        return False
    conditions_per_batch = table.groupby(context.batch_key, observed=True)[
        context.condition_key
    ].nunique()
    batches_per_condition = table.groupby(context.condition_key, observed=True)[
        context.batch_key
    ].nunique()
    return bool((conditions_per_batch <= 1).all() or (batches_per_condition <= 1).all())


def plan_analysis(
    adata: AnnData,
    *,
    context: Optional[AnalysisContext | Mapping[str, Any]] = None,
    profile: str = "auto",
    stages: Optional[Sequence[str]] = None,
    **context_hints: Any,
) -> AnalysisPlan:
    """Build a conservative first-pass plan from project metadata.

    The plan checks structural prerequisites and highlights decisions that
    require biological review.  It intentionally does not auto-apply batch
    correction, final annotation, formal differential inference, or malignant
    labels.
    """
    selected_stages = list(stages or ("qc", "preprocess", "analysis"))
    invalid = sorted(set(selected_stages) - {"qc", "preprocess", "analysis"})
    if invalid:
        raise ValueError(f"Invalid plan stages: {invalid}.")

    inferred = infer_analysis_context(adata, context=context, **context_hints)
    resolved_profile = _resolve_profile(inferred, profile)
    decisions: list[DecisionCard] = []
    required_metadata: list[str] = []

    counts_present = "counts" in adata.layers
    decisions.append(
        DecisionCard(
            stage="qc",
            status="READY" if counts_present else "REVIEW",
            decision="raw_counts_source",
            recommended="adata.layers['counts']",
            applied=(
                "adata.layers['counts']" if counts_present else "adata.X (semantics unverified)"
            ),
            reason=(
                "A dedicated raw-count layer is available."
                if counts_present
                else "No dedicated counts layer is present, so raw-count semantics must be verified."
            ),
            evidence=["adata.layers", "adata.X"],
            next_action=(
                "No action required."
                if counts_present
                else "Verify that adata.X contains raw counts and copy it to layers['counts'] before normalization."
            ),
            rerun_scope=_RERUN_SCOPE["qc"],
            priority="required" if not counts_present else "optional",
        )
    )

    comparative = _is_comparative_objective(inferred.study_objective) or bool(
        inferred.condition_key
    )
    sample_key_present = bool(inferred.sample_key and inferred.sample_key in adata.obs.columns)
    if sample_key_present:
        sample_status: DecisionStatus = "READY"
        sample_reason = f"Biological samples are represented by obs[{inferred.sample_key!r}]."
    elif comparative:
        sample_status = "BLOCKED"
        sample_reason = (
            f"Configured sample_key {inferred.sample_key!r} is absent from adata.obs."
            if inferred.sample_key
            else "A comparative objective requires an explicit biological sample identifier."
        )
        required_metadata.append("sample_key")
    else:
        sample_status = "REVIEW"
        sample_reason = "No biological sample identifier was found; formal sample-level inference is unavailable."
        required_metadata.append("sample_key (required for sample-level inference)")
    decisions.append(
        DecisionCard(
            stage="analysis",
            status=sample_status,
            decision="biological_sample",
            recommended="explicit obs column identifying biological samples",
            applied=inferred.sample_key,
            reason=sample_reason,
            evidence=["adata.obs", "ProjectContext.sample_key"],
            next_action=(
                "Confirm that this column identifies independent biological samples."
                if sample_key_present
                else "Add sample metadata before condition-level DE or composition testing."
            ),
            rerun_scope=_RERUN_SCOPE["analysis"],
            priority="blocking" if sample_status == "BLOCKED" else "review",
        )
    )

    condition_key_present = bool(
        inferred.condition_key and inferred.condition_key in adata.obs.columns
    )
    if comparative and not condition_key_present:
        required_metadata.append("condition_key")
        decisions.append(
            DecisionCard(
                stage="analysis",
                status="BLOCKED",
                decision="comparison_condition",
                recommended="explicit obs column encoding the intended comparison",
                applied=inferred.condition_key,
                reason=(
                    f"Configured condition_key {inferred.condition_key!r} is absent from adata.obs."
                    if inferred.condition_key
                    else "The comparative study objective has no condition_key."
                ),
                evidence=["adata.obs", "ProjectContext.condition_key"],
                next_action="Add or correct the condition metadata before formal comparisons.",
                rerun_scope=_RERUN_SCOPE["analysis"],
                priority="blocking",
            )
        )

    replicates = _condition_replicates(adata, inferred)
    if condition_key_present:
        experimental_unit_present = bool(
            inferred.experimental_unit_key and inferred.experimental_unit_key in adata.obs.columns
        )
        if not experimental_unit_present:
            inference_status: DecisionStatus = "BLOCKED"
            inference_reason = (
                f"Configured experimental_unit_key {inferred.experimental_unit_key!r} is absent from adata.obs."
                if inferred.experimental_unit_key
                else "No experimental unit is available for the requested condition comparison."
            )
            required_metadata.append("experimental_unit_key")
        elif not replicates:
            inference_status = "REVIEW"
            inference_reason = "The condition-by-experimental-unit design could not be summarized."
        elif len(replicates) < 2:
            inference_status = "REVIEW"
            inference_reason = (
                "Fewer than two conditions are represented in the comparison metadata."
            )
        elif min(replicates.values()) < 2:
            inference_status = "REVIEW"
            inference_reason = (
                "At least one condition has fewer than two independent experimental units; "
                "formal inference is not supported by replication."
            )
        else:
            inference_status = "REVIEW"
            inference_reason = (
                "Replicates are present, but the contrast, covariates, and paired structure still "
                "require confirmation before pseudobulk inference."
            )
        decisions.append(
            DecisionCard(
                stage="analysis",
                status=inference_status,
                decision="sample_level_inference",
                recommended="pseudobulk or an explicit sample-aware model",
                applied={
                    "condition_key": inferred.condition_key,
                    "experimental_unit_key": inferred.experimental_unit_key,
                    "paired_key": inferred.paired_key,
                    "replicates_per_condition": replicates,
                },
                reason=inference_reason,
                evidence=["condition-by-experimental-unit table"],
                next_action="Confirm the contrast, replicate counts, pairing, and covariates before formal DE/proportion claims.",
                rerun_scope=_RERUN_SCOPE["analysis"],
                priority="blocking" if inference_status == "BLOCKED" else "review",
            )
        )

    batch_values = (
        int(adata.obs[inferred.batch_key].nunique())
        if inferred.batch_key and inferred.batch_key in adata.obs.columns
        else 0
    )
    confounded = _batch_condition_confounded(adata, inferred)
    decisions.append(
        DecisionCard(
            stage="preprocess",
            status="REVIEW" if batch_values > 1 else "READY",
            decision="integration_policy",
            recommended=(
                "compare unintegrated and candidate integrated representations"
                if batch_values > 1
                else "no integration by default"
            ),
            applied="not selected by project planning",
            reason=(
                "Batch and condition appear confounded; integration could remove biological signal."
                if confounded
                else (
                    "Multiple batches require diagnostic review before correction."
                    if batch_values > 1
                    else "No multi-batch structure requiring correction was detected."
                )
            ),
            evidence=["batch-by-condition table", "unintegrated PCA/UMAP"],
            next_action=(
                "Do not finalize integration until biological preservation and batch mixing are compared."
                if batch_values > 1
                else "Keep integration disabled unless downstream diagnostics show a batch-driven structure."
            ),
            rerun_scope=_RERUN_SCOPE["preprocess"],
            priority="review" if batch_values > 1 else "optional",
        )
    )

    decisions.append(
        DecisionCard(
            stage="analysis",
            status="REVIEW",
            decision="annotation_boundary",
            recommended="broad lineage first, then fine labels with multi-marker evidence",
            applied="first-pass automated annotation",
            reason="Automated labels are evidence for review, not final biological truth.",
            evidence=["positive markers", "negative markers", "reference agreement", "cluster QC"],
            next_action="Review low-confidence and conflicting clusters before promoting labels to final annotations.",
            rerun_scope=_RERUN_SCOPE["analysis"],
            priority="review",
        )
    )

    if inferred.enables_tumor_module:
        decisions.append(
            DecisionCard(
                stage="tumor",
                status="REVIEW",
                decision="malignancy_boundary",
                recommended="CNV plus transcriptional and sample-aware evidence",
                applied="not inferred by project planning",
                reason="Epithelial markers alone cannot establish malignancy.",
                evidence=[
                    "CNV burden",
                    "tumor programs",
                    "normal reference",
                    "inter-patient heterogeneity",
                ],
                next_action="Run and review a separate malignancy interpretation after broad annotation.",
                rerun_scope=_RERUN_SCOPE["tumor"],
                priority="review",
            )
        )

    statuses = [decision.status for decision in decisions]
    status = _worst_status(statuses, empty="READY")
    blockers = list(
        dict.fromkeys(decision.reason for decision in decisions if decision.status == "BLOCKED")
    )
    assumptions = [
        f"Dataset profile: {resolved_profile}.",
        f"Raw-count source: {'layers[counts]' if counts_present else 'adata.X, pending verification'}.",
        (
            f"Experimental unit: obs[{inferred.experimental_unit_key!r}]."
            if inferred.experimental_unit_key
            else "Experimental unit is not yet explicit."
        ),
        "Integrated embeddings, if used, are for representation/clustering; unintegrated expression remains the interpretation source.",
        "Automated annotation and tumor calls remain first-pass until human review is recorded.",
    ]
    risks = [
        "Pseudo-replication if cells are treated as independent samples.",
        "Over-integration if batch is confounded with condition or patient biology.",
        "Overconfident annotation if a single marker or one reference method is treated as truth.",
    ]
    if inferred.enables_tumor_module:
        risks.append(
            "Malignant-cell misclassification if epithelial identity is used without CNV/multi-evidence support."
        )
    if not counts_present:
        risks.append(
            "Irreversible layer ambiguity if raw counts are not preserved before preprocessing."
        )

    if blockers:
        smallest_next_step = (
            "Resolve the blocking project metadata before running formal downstream analysis."
        )
    elif any(item.status == "REVIEW" for item in decisions):
        smallest_next_step = (
            "Confirm the REVIEW decisions, then run the conservative first-pass workflow."
        )
    else:
        smallest_next_step = (
            "Run the first-pass workflow and inspect review_run() before interpretation."
        )

    return AnalysisPlan(
        profile=resolved_profile,
        status=status,  # type: ignore[arg-type]
        stages=selected_stages,
        context=inferred,
        decisions=decisions,
        required_metadata=list(dict.fromkeys(required_metadata)),
        blockers=blockers,
        assumptions=assumptions,
        highest_risk_failure_modes=risks,
        smallest_next_step=smallest_next_step,
    )


def _readiness_for_stage(stage: str, payload: Mapping[str, Any]) -> Mapping[str, Any]:
    keys = {
        "qc": ("qc_readiness", "qc_handoff_readiness"),
        "preprocess": ("preprocess_readiness", "analysis_handoff_readiness"),
        "analysis": ("analysis_readiness",),
        "tumor": ("readiness",),
    }[stage]
    for key in keys:
        value = payload.get(key)
        if isinstance(value, Mapping):
            return value
    bundle = payload.get("evidence_bundle")
    if isinstance(bundle, Mapping):
        return {
            "status": bundle.get("status", "unknown"),
            "score": bundle.get("confidence"),
        }
    return {}


def _readiness_reasons(readiness: Mapping[str, Any], payload: Mapping[str, Any]) -> list[str]:
    reasons: list[str] = []
    for key in ("blockers", "review_reasons", "reasons", "issues", "warnings"):
        reasons.extend(_strings(readiness.get(key)))
    if not reasons:
        reasons.extend(_strings(payload.get("warnings")))
    return list(dict.fromkeys(reasons))


def _confidence(readiness: Mapping[str, Any]) -> Optional[float]:
    value = readiness.get("score", readiness.get("confidence"))
    if not isinstance(value, (int, float)):
        return None
    score = float(value)
    if score > 1:
        score /= 100.0
    return max(0.0, min(1.0, score))


def _action_records(payload: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    bundle = payload.get("evidence_bundle")
    candidates: list[Mapping[str, Any]] = []
    if isinstance(bundle, Mapping):
        candidates.extend(_records(bundle.get("action_items")))
    candidates.extend(_records(payload.get("review_action_items")))
    candidates.extend(_records(payload.get("action_items")))
    unique: list[Mapping[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for item in candidates:
        key = (str(item.get("priority", "review")), str(item.get("action", "")))
        if key in seen or not key[1]:
            continue
        seen.add(key)
        unique.append(item)
    return unique


def _decision_records(payload: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    bundle = payload.get("evidence_bundle")
    candidates: list[Mapping[str, Any]] = []
    if isinstance(bundle, Mapping):
        candidates.extend(_records(bundle.get("decisions")))
    if not candidates:
        candidates.extend(_records(payload.get("decision_table")))
    return candidates


def _decision_status(record: Mapping[str, Any]) -> DecisionStatus:
    recommended = record.get("recommended")
    applied = record.get("applied")
    confidence = record.get("confidence")
    if recommended is not None and applied is None:
        return "REVIEW"
    if isinstance(confidence, (int, float)) and float(confidence) < 0.5:
        return "REVIEW"
    return "READY"


def _stage_review(stage: str, namespace: Any) -> tuple[StageReview, list[DecisionCard]]:
    present = isinstance(namespace, Mapping) and bool(namespace)
    summary = namespace.get(UnsKeys.REVIEW_SUMMARY) if present else None
    payload = _review_payload(summary)
    readiness = _readiness_for_stage(stage, payload) if payload else {}
    raw_status = str(readiness.get("status", "unknown"))
    status = _normalize_status(raw_status, present=bool(payload))
    reasons = _readiness_reasons(readiness, payload)
    if present and not payload:
        reasons.append("The stage namespace exists but no review_summary was recorded.")
    claim_boundary = payload.get("claim_boundary")
    stage_next_action = {
        "BLOCKED": f"Resolve the {stage} blockers before relying on downstream results.",
        "REVIEW": f"Review the flagged {stage} evidence before finalizing this stage.",
        "READY": "No blocking action is recorded for this stage.",
        "NOT_RUN": "Run this stage only if it is required by the project objective.",
    }[status]
    stage_review = StageReview(
        stage=stage,
        status=status,
        raw_status=raw_status,
        confidence=_confidence(readiness),
        reasons=reasons,
        claim_boundary=str(claim_boundary) if claim_boundary is not None else None,
        next_action=stage_next_action,
        rerun_scope=_RERUN_SCOPE[stage],
    )

    items: list[DecisionCard] = [
        DecisionCard(
            stage=stage,
            status="REVIEW" if status == "NOT_RUN" else status,  # type: ignore[arg-type]
            decision="stage_readiness",
            recommended="READY before final interpretation",
            applied=raw_status if status != "NOT_RUN" else "not run",
            reason="; ".join(reasons) if reasons else stage_next_action,
            evidence=[f'adata.uns["{SCLUCID_ROOT}"]["{stage}"]["review_summary"]'],
            next_action=stage_next_action,
            rerun_scope=_RERUN_SCOPE[stage],
            priority=(
                "blocking"
                if status == "BLOCKED"
                else "review" if status == "REVIEW" else "optional"
            ),
            source="stage_readiness",
        )
    ]

    for action in _action_records(payload):
        priority = str(action.get("priority", "review"))
        action_status: DecisionStatus = (
            "BLOCKED" if priority == "blocking" else "READY" if priority == "optional" else "REVIEW"
        )
        items.append(
            DecisionCard(
                stage=stage,
                status=action_status,
                decision="review_action",
                recommended=action.get("action"),
                reason=str(action.get("rationale", "")),
                evidence=_strings(action.get("evidence_keys") or action.get("evidence_key")),
                next_action=str(action.get("action", "")),
                rerun_scope=_RERUN_SCOPE[stage],
                priority=priority,
                source="evidence_bundle.action_items",
            )
        )

    for decision in _decision_records(payload):
        evidence_names = []
        for item in _records(decision.get("evidence")):
            evidence_names.append(str(item.get("name") or item.get("source") or "evidence"))
        decision_name = str(decision.get("parameter") or decision.get("decision") or "parameter")
        decision_status = _decision_status(decision)
        items.append(
            DecisionCard(
                stage=stage,
                status=decision_status,
                decision=decision_name,
                recommended=decision.get("recommended"),
                applied=decision.get("applied"),
                reason=str(
                    decision.get("downstream_impact")
                    or decision.get("rationale")
                    or f"Decision source: {decision.get('source', 'unknown')}."
                ),
                evidence=evidence_names,
                next_action=(
                    "Confirm or override this recommendation and record the reason."
                    if decision_status == "REVIEW"
                    else "No action required unless project context conflicts with the applied value."
                ),
                rerun_scope=_RERUN_SCOPE[stage],
                priority="review" if decision_status == "REVIEW" else "optional",
                source=str(decision.get("source", "evidence_bundle.decisions")),
            )
        )
    return stage_review, items


def review_run(adata: AnnData, *, store: bool = True) -> RunReview:
    """Normalize all recorded module evidence into one decision-oriented review.

    ``review_run`` is safe on partially processed objects.  Missing stages are
    labelled ``NOT_RUN`` rather than treated as automatic failures.
    """
    root = adata.uns.get(SCLUCID_ROOT, {})
    if not isinstance(root, Mapping):
        root = {}

    context_payload = root.get(UnsKeys.ANALYSIS_CONTEXT)
    if not isinstance(context_payload, Mapping):
        context_payload = infer_analysis_context(adata).to_dict()

    stages: list[StageReview] = []
    items: list[DecisionCard] = []
    for stage in _STAGE_ORDER:
        stage_review, stage_items = _stage_review(stage, root.get(stage))
        stages.append(stage_review)
        items.extend(stage_items)

    recorded_statuses = [stage.status for stage in stages if stage.status != "NOT_RUN"]
    overall = _worst_status(recorded_statuses)
    action_items = [
        item
        for item in items
        if item.status in {"BLOCKED", "REVIEW"}
        and item.decision != "stage_readiness"
        and item.next_action
    ]
    if not action_items:
        action_items = [
            item for item in items if item.status in {"BLOCKED", "REVIEW"} and item.next_action
        ]
    action_items.sort(
        key=lambda item: (
            0 if item.status == "BLOCKED" else 1,
            _STAGE_ORDER.index(item.stage),
        )
    )
    next_actions = list(dict.fromkeys(item.next_action for item in action_items))

    context = AnalysisContext.model_validate(context_payload)
    assumptions = [
        (
            f"Formal inference unit is obs[{context.experimental_unit_key!r}]."
            if context.experimental_unit_key
            else "Formal inference unit is not explicit."
        ),
        "Cell-level marker tests are exploratory unless a sample-aware method is recorded.",
        "Automated labels remain first-pass until annotation review is completed.",
    ]
    risks = [
        "Pseudo-replication if cell-level tests are interpreted as sample-level evidence.",
        "Over-integration if batch correction removes condition or patient biology.",
        "Review fatigue if warnings are read without resolving their recorded next actions.",
    ]
    if context.enables_tumor_module:
        risks.append(
            "Malignancy claims require CNV and transcriptional evidence beyond epithelial markers."
        )

    review = RunReview(
        overall_status=overall,
        context=context.to_dict(),
        stages=stages,
        items=items,
        next_actions=next_actions,
        assumptions=assumptions,
        highest_risk_failure_modes=risks,
    )
    if store:
        namespace = ensure_sclucid_namespace(adata)
        namespace[UnsKeys.RUN_REVIEW] = sanitize_for_hdf5(review.to_dict())
    return review


__all__ = [
    "ANALYSIS_PLAN_SCHEMA_VERSION",
    "AnalysisPlan",
    "DecisionCard",
    "ProjectContext",
    "RUN_REVIEW_SCHEMA_VERSION",
    "RunReview",
    "StageReview",
    "plan_analysis",
    "review_run",
]
