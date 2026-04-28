"""Tests for shared scLucid evidence schema."""

from scLucid.utils.evidence import (
    EVIDENCE_SCHEMA_VERSION,
    DecisionRecord,
    EvidenceBundle,
    EvidenceItem,
    ReviewAction,
    model_to_dict,
)


def test_evidence_bundle_serializes_to_json_ready_dict():
    item = EvidenceItem(
        source="metric",
        name="retention_rate",
        value=0.91,
        confidence=0.8,
        rationale="High retention supports continuing downstream.",
        related_keys=["benchmark_summary.retention.retention_rate"],
    )
    decision = DecisionRecord(
        parameter="min_genes",
        recommended=200,
        applied=200,
        source="recommendation",
        confidence=0.8,
        evidence=[item],
    )
    action = ReviewAction(
        priority="optional",
        action="Archive QC summary.",
        rationale="No blocking issues detected.",
        evidence_keys=["qc_readiness"],
    )
    bundle = EvidenceBundle(
        module="qc",
        stage="run_standard_qc",
        status="ready",
        confidence=0.95,
        context={"sample_key": "sampleID"},
        decisions=[decision],
        evidence_chain=[item],
        action_items=[action],
    )

    payload = model_to_dict(bundle)

    assert payload["schema_version"] == EVIDENCE_SCHEMA_VERSION
    assert payload["module"] == "qc"
    assert payload["decisions"][0]["evidence"][0]["source"] == "metric"
    assert payload["action_items"][0]["priority"] == "optional"
