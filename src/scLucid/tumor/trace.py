"""Review-summary schema and validation for the tumor workflow.

Mirrors the analysis/qc/preprocess trace pattern: a small frozen contract
that declares what a benchmark-grade tumor review summary must contain,
plus helpers to validate and enrich it.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Dict, List

from anndata import AnnData

from ..utils import sanitize_for_hdf5
from ..utils.evidence import EvidenceBundle, EvidenceItem, ReviewAction, model_to_dict
from scLucid.utils.contracts import _review_payload

TUMOR_TRACE_SCHEMA_VERSION = "1.0"
TUMOR_MODULE_MATURITY_SCHEMA_VERSION = "1.0"

TUMOR_REVIEW_SUMMARY_REQUIRED_KEYS = {
    "schema_version",
    "module",
    "workflow_name",
    "requested_steps",
    "completed_steps",
    "claim_boundary",
    "readiness",
    "warnings",
    "action_items",
    "evidence_sources",
}

TUMOR_STABLE_ENTRYPOINTS = (
    "scLucid.tumor.run_tumor_analysis",
    "scLucid.tumor.run_tumor_analysis_expert",
    "scLucid.tumor.cnv.infer_cnv",
    "scLucid.tumor.malignancy.run_malignancy_interpretation",
    "scLucid.tumor.microenvironment.deconvolve_tme",
    "scLucid.tumor.therapy.predict_therapy_response",
)

TUMOR_EXPECTED_OUTPUTS = (
    'adata.uns["sclucid"]["tumor"]["review_summary"]',
    'adata.uns["sclucid"]["tumor"]["step_results"]',
    'adata.uns["sclucid"]["tumor"]["execution_trace"]',
)


def validate_tumor_review_summary(
    summary: Mapping[str, Any],
    *,
    raise_on_error: bool = False,
) -> List[str]:
    """Validate that a tumor review summary satisfies the schema contract."""
    errors: List[str] = []
    if not isinstance(summary, Mapping):
        errors.append("Tumor review summary must be a mapping.")
        if raise_on_error:
            raise ValueError("; ".join(errors))
        return errors

    missing = sorted(TUMOR_REVIEW_SUMMARY_REQUIRED_KEYS - set(summary.keys()))
    if missing:
        errors.append(f"Tumor review summary missing required keys: {missing}")

    if summary.get("module") != "tumor":
        errors.append("Tumor review summary 'module' must be 'tumor'.")

    if summary.get("schema_version") != TUMOR_TRACE_SCHEMA_VERSION:
        errors.append(
            f"Tumor review summary schema_version must be '{TUMOR_TRACE_SCHEMA_VERSION}'."
        )

    readiness = summary.get("readiness")
    if not isinstance(readiness, Mapping):
        errors.append("Tumor review summary 'readiness' must be a mapping.")
    else:
        if "status" not in readiness:
            errors.append("Tumor readiness must contain 'status'.")
        if "score" not in readiness:
            errors.append("Tumor readiness must contain 'score'.")

    claim = summary.get("claim_boundary")
    if claim not in ("validated_core", "heuristic", "exploratory", "unavailable"):
        errors.append(
            "Tumor claim_boundary must be one of "
            "validated_core/heuristic/exploratory/unavailable."
        )

    action_items = summary.get("action_items")
    if not isinstance(action_items, (list, Mapping)):
        errors.append("Tumor review summary 'action_items' must be a list or mapping.")

    if raise_on_error and errors:
        raise ValueError("; ".join(errors))
    return errors


def enrich_tumor_review_summary(
    summary: Dict[str, Any],
    *,
    adata: AnnData,
    step_results: List[Any],
) -> Dict[str, Any]:
    """Add benchmark-grade fields to a tumor review summary.

    Ensures the summary satisfies ``TUMOR_REVIEW_SUMMARY_REQUIRED_KEYS`` and
    includes an evidence bundle suitable for cross-module audit.
    """
    summary = dict(summary)
    summary.setdefault("schema_version", TUMOR_TRACE_SCHEMA_VERSION)
    summary.setdefault("module", "tumor")
    summary.setdefault("workflow_name", "tumor_analysis")
    summary.setdefault("warnings", [])
    summary.setdefault("evidence_sources", [])

    # Attach evidence bundle
    summary["evidence_bundle"] = _build_tumor_evidence_bundle(summary, step_results)
    summary["module_maturity"] = _build_tumor_module_maturity_assessment(summary)

    return summary


def get_tumor_module_contract() -> Dict[str, Any]:
    """Return the frozen tumor module maturity contract."""
    return {
        "schema_version": TUMOR_MODULE_MATURITY_SCHEMA_VERSION,
        "module": "tumor",
        "stable_entrypoints": list(TUMOR_STABLE_ENTRYPOINTS),
        "required_review_keys": sorted(TUMOR_REVIEW_SUMMARY_REQUIRED_KEYS),
        "expected_outputs": list(TUMOR_EXPECTED_OUTPUTS),
        "canonical_namespace": 'adata.uns["sclucid"]["tumor"]',
        "readiness_key": "readiness",
        "claim_boundary_key": "claim_boundary",
        "step_results_key": "step_results",
    }


def validate_tumor_module_completeness(
    adata: AnnData,
    *,
    require_ready: bool = False,
    raise_on_error: bool = False,
) -> Dict[str, Any]:
    """Validate that an AnnData object contains a benchmark-grade tumor result."""
    issues: List[str] = []
    warnings_list: List[str] = []
    tumor_ns = adata.uns.get("sclucid", {}).get("tumor", {})
    if not isinstance(tumor_ns, Mapping):
        issues.append('Missing or invalid adata.uns["sclucid"]["tumor"] namespace.')
        tumor_ns = {}

    review_summary = tumor_ns.get("review_summary")
    errors = validate_tumor_review_summary(review_summary or {})
    issues.extend(errors)

    if "step_results" not in tumor_ns:
        issues.append('Missing adata.uns["sclucid"]["tumor"]["step_results"].')

    readiness = (review_summary or {}).get("readiness", {})
    if require_ready and readiness.get("status") != "ready":
        issues.append(f"Tumor readiness is {readiness.get('status')!r}, expected 'ready'.")

    maturity = _build_tumor_module_maturity_assessment(review_summary or {})
    if maturity.get("status") == "incomplete":
        issues.extend(maturity.get("issues", []))
    elif maturity.get("status") == "review_required":
        warnings_list.extend(maturity.get("review_required", []))

    result = {
        "schema_version": TUMOR_MODULE_MATURITY_SCHEMA_VERSION,
        "module": "tumor",
        "valid": len(issues) == 0,
        "status": "valid" if not issues else "invalid",
        "issues": list(dict.fromkeys(str(i) for i in issues)),
        "warnings": list(dict.fromkeys(str(w) for w in warnings_list)),
        "maturity": maturity,
        "contract": get_tumor_module_contract(),
    }
    if result["issues"] and raise_on_error:
        raise ValueError("; ".join(result["issues"]))
    return result


def _build_tumor_evidence_bundle(
    summary: Mapping[str, Any],
    step_results: List[Any],
) -> Dict[str, Any]:
    readiness = summary.get("readiness", {}) if isinstance(summary, Mapping) else {}
    evidence_chain = [
        EvidenceItem(
            source="context",
            name="analysis_input_context",
            value={"step_results_present": len(step_results) > 0},
            rationale="Tumor stage depends on upstream analysis outputs.",
            related_keys=["step_results"],
        ),
        EvidenceItem(
            source="output_health",
            name="claim_boundary",
            value=summary.get("claim_boundary") if isinstance(summary, Mapping) else None,
            rationale="Claim boundary declares the evidentiary strength of tumor outputs.",
            related_keys=["claim_boundary"],
        ),
        EvidenceItem(
            source="output_health",
            name="tumor_readiness",
            value=readiness,
            rationale="Tumor readiness summarizes whether downstream interpretation is safe.",
            related_keys=["readiness"],
        ),
    ]
    raw_actions = summary.get("action_items") or []
    if isinstance(raw_actions, Mapping):
        raw_actions = raw_actions.values()
    actions = [ReviewAction(**item) for item in raw_actions if isinstance(item, Mapping)]
    bundle = EvidenceBundle(
        module="tumor",
        stage="run_tumor_analysis",
        status=str(readiness.get("status", "unknown")),
        confidence=readiness.get("score") if isinstance(readiness, Mapping) else None,
        context={"workflow_name": "tumor_analysis"},
        evidence_chain=evidence_chain,
        action_items=actions,
        related_review_keys=sorted(TUMOR_REVIEW_SUMMARY_REQUIRED_KEYS),
    )
    return model_to_dict(bundle)


def _build_tumor_module_maturity_assessment(
    summary: Mapping[str, Any],
) -> Dict[str, Any]:
    payload = _review_payload(summary)
    required = set(TUMOR_REVIEW_SUMMARY_REQUIRED_KEYS)
    required.discard("module_maturity")
    missing = sorted(required - set(payload.keys()))
    issues = [f"missing_required_key:{key}" for key in missing]
    review_required = []

    readiness = payload.get("readiness", {}) if isinstance(payload, Mapping) else {}
    if readiness.get("status") == "blocked":
        issues.extend(readiness.get("reasons", []))
    elif readiness.get("status") in ("degraded", "review_required"):
        review_required.extend(readiness.get("reasons", []))

    if issues:
        status = "incomplete"
    elif review_required:
        status = "review_required"
    else:
        status = "complete"

    return _json_safe(
        {
            "schema_version": TUMOR_MODULE_MATURITY_SCHEMA_VERSION,
            "module": "tumor",
            "status": status,
            "issues": issues,
            "review_required": review_required,
            "contract": get_tumor_module_contract(),
            "summary": (
                "Tumor review summary satisfies the benchmark module contract."
                if status == "complete"
                else "Tumor review summary is present but requires review."
                if status == "review_required"
                else "Tumor review summary does not satisfy the benchmark module contract."
            ),
        }
    )


def _json_safe(value: Any) -> Any:
    return sanitize_for_hdf5(value)


__all__ = [
    "TUMOR_MODULE_MATURITY_SCHEMA_VERSION",
    "TUMOR_REVIEW_SUMMARY_REQUIRED_KEYS",
    "TUMOR_TRACE_SCHEMA_VERSION",
    "enrich_tumor_review_summary",
    "get_tumor_module_contract",
    "validate_tumor_module_completeness",
    "validate_tumor_review_summary",
]
