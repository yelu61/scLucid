"""Trust checks for built-in scLucid marker and geneset resources."""

from __future__ import annotations

import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .manager import Manager, get_marker_manager, load_gene_set_manager
from .resource_loader import load_json, load_toml, resource_exists

_GENE_SYMBOL_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
_ACTIVE_MARKER_RESOURCES = {
    "marker_registry_human": "registry_human",
    "marker_registry_mouse": "registry_mouse",
    "marker_tissue_human": "tissue_human",
    "marker_tumor_human": "tumor_human",
}
_ACTIVE_GENESET_RESOURCES = {
    "genesets_cancer_signatures": ("human", "cancer_signatures"),
    "genesets_cancer_hallmarks": ("human", "cancer_hallmarks"),
}
_MIN_MARKERS_BY_GRANULARITY = {
    "compartment": 1,
    "lineage": 3,
    "subtype": 3,
    "tissue_subtype": 3,
    "state": 3,
    "artifact": 3,
    "program": 5,
    "epithelial_support": 3,
    "tumor_type_hint": 3,
    "cancer_subtype": 3,
    "reference_anchor": 3,
}
_MAX_MARKERS_BY_GRANULARITY = {
    "compartment": 12,
    "lineage": 12,
    "subtype": 10,
    "tissue_subtype": 10,
    "state": 15,
    "artifact": 50,
    "program": 50,
    "epithelial_support": 15,
    "tumor_type_hint": 15,
    "cancer_subtype": 15,
    "reference_anchor": 12,
}


@dataclass(frozen=True)
class ResourceAuditIssue:
    """Structured issue emitted by the resource trust audit."""

    resource: str
    issue: str
    severity: str
    detail: Any = None

    def as_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable representation."""
        return {
            "resource": self.resource,
            "issue": self.issue,
            "severity": self.severity,
            "detail": self.detail,
        }


def load_resource_manifest() -> dict[str, Any]:
    """Load the marker resource manifest."""
    return load_toml("marker_resource_manifest.toml")


def load_reference_index() -> dict[str, dict[str, Any]]:
    """Load references keyed by ``source_id``."""
    data = load_toml("references.toml")
    return {
        str(item["source_id"]): dict(item)
        for item in data.get("references", [])
        if isinstance(item, dict) and item.get("source_id")
    }


def load_marker_curation_literature_index() -> list[dict[str, Any]]:
    """Load the local curation literature queue when present."""
    path = (
        Path(__file__).parents[3]
        / "docs"
        / "marker_resources"
        / "marker_curation_literature_index.jsonl"
    )
    if not path.is_file():
        raise FileNotFoundError(path)
    rows: list[dict[str, Any]] = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                rows.append(json.loads(stripped))
    return rows


def _resource_filename(resource_name: str) -> str:
    if resource_name.startswith("genesets_"):
        return f"{resource_name}.json"
    return f"{resource_name}.toml"


def _source_ids_from(value: Any) -> set[str]:
    source_ids: set[str] = set()
    if isinstance(value, dict):
        for key, nested in value.items():
            if key in {"source_ids", "framework_source_ids"}:
                if isinstance(nested, str):
                    source_ids.add(nested)
                elif isinstance(nested, list):
                    source_ids.update(str(item) for item in nested if item)
                continue
            source_ids.update(_source_ids_from(nested))
    elif isinstance(value, list):
        for item in value:
            source_ids.update(_source_ids_from(item))
    return source_ids


def _add_issue(
    issues: list[ResourceAuditIssue],
    resource: str,
    issue: str,
    severity: str,
    detail: Any = None,
) -> None:
    issues.append(ResourceAuditIssue(resource, issue, severity, detail))


def classify_literature_resource_utility(row: dict[str, Any]) -> dict[str, Any]:
    """
    Classify a curation-batch reference into resource-use categories.

    This is a conservative triage helper. It does not claim that markers have
    been extracted from the paper; it records where the paper should be reviewed
    and what kind of resource payload it can support.
    """
    title = str(row.get("title", "")).lower()
    source_type = str(row.get("source_type", "")).lower()
    recommendation = str(row.get("recommended_resource", "")).lower()
    batch_id = str(row.get("batch_id", ""))
    status = str(row.get("curation_status") or row.get("status") or "")

    target_resources: list[str] = []
    utilities: list[str] = []

    def add_target(name: str) -> None:
        if name not in target_resources:
            target_resources.append(name)

    def add_utility(name: str) -> None:
        if name not in utilities:
            utilities.append(name)

    if "marker_registry" in recommendation or "single_cell_atlas" in source_type:
        add_target("marker_registry_human")
        add_utility("marker_core")
    if "marker_tissue" in recommendation or batch_id == "01":
        add_target("marker_tissue_human")
        add_utility("tissue_context")
    if (
        "marker_tumor" in recommendation
        or "pan_cancer" in source_type
        or "disease_atlas" in source_type
        or "cancer" in title
        or "tumor" in title
        or batch_id in {"11", "12"}
    ):
        add_target("marker_tumor_human")
        add_utility("tumor_context")
    if "geneset" in recommendation or "signature" in title or "hallmark" in title:
        add_target("genesets_cancer_signatures")
        add_utility("geneset_scoring")
    if "hallmark" in title:
        add_target("genesets_cancer_hallmarks")
    if "naming" in recommendation or "nomenclature" in title or "guideline" in title:
        add_target("marker_aliases")
        add_utility("nomenclature_reference")
    if "validation" in recommendation or "framework" in recommendation:
        add_utility("validation_reference")
    if "computational" in source_type or "tool" in source_type or "database" in source_type:
        add_utility("benchmark_reference")
    if "reference_only" in recommendation:
        add_utility("reference_only")

    if not utilities:
        add_utility("context_reference")
    if not target_resources and "reference_only" not in utilities:
        add_target("references")

    extraction_status = (
        "queued_zotero_resolution"
        if status == "queued"
        else "triaged_not_extracted"
    )
    priority = "medium"
    if status == "queued":
        priority = "high"
    elif any(item in utilities for item in {"marker_core", "tumor_context", "geneset_scoring"}):
        priority = "high"
    elif utilities == ["reference_only"] or utilities == ["context_reference"]:
        priority = "low"

    return {
        "resource_utility": utilities,
        "target_resources": target_resources,
        "curation_priority": priority,
        "fulltext_review_required": "reference_only" not in utilities,
        "extraction_status": extraction_status,
    }


def audit_marker_entry_quality() -> dict[str, Any]:
    """
    Summarize marker-entry curation gaps without treating them as code errors.

    Contract failures belong in ``audit_marker_resources``. This function reports
    biological curation debt: thin marker lists, oversized lists, missing negative
    markers, and entries still waiting for review.
    """
    resources: dict[str, Any] = {}
    gaps: list[dict[str, Any]] = []

    def _effective_negative_markers(cell: Any) -> set[str]:
        markers = {str(gene) for gene in cell.negative_markers}
        parent = cell.parent
        while parent is not None:
            markers.update(str(gene) for gene in parent.negative_markers)
            parent = parent.parent
        return markers

    for resource_name, manager_name in _ACTIVE_MARKER_RESOURCES.items():
        mgr = Manager(manager_name, case_sensitive=True)
        by_granularity: dict[str, list[int]] = {}
        review_status = Counter(
            str(cell.metadata.get("review_status", "missing"))
            for cell in mgr.CELLS.values()
        )

        for name, cell in mgr.CELLS.items():
            granularity = str(cell.metadata.get("granularity", "missing"))
            kind = str(cell.metadata.get("kind", "missing"))
            marker_count = len(cell.markers)
            negative_count = len(_effective_negative_markers(cell))
            by_granularity.setdefault(granularity, []).append(marker_count)

            if cell.minor and marker_count == 0:
                continue
            if granularity.endswith("_collection") or (
                marker_count == 0 and kind.endswith("_collection")
            ):
                continue
            min_markers = _MIN_MARKERS_BY_GRANULARITY.get(granularity)
            max_markers = _MAX_MARKERS_BY_GRANULARITY.get(granularity)
            if min_markers is not None and marker_count < min_markers:
                gaps.append(
                    {
                        "resource": resource_name,
                        "entry": name,
                        "gap": "thin_marker_set",
                        "granularity": granularity,
                        "marker_count": marker_count,
                        "recommended_min": min_markers,
                    }
                )
            if max_markers is not None and marker_count > max_markers:
                gaps.append(
                    {
                        "resource": resource_name,
                        "entry": name,
                        "gap": "oversized_marker_set",
                        "granularity": granularity,
                        "marker_count": marker_count,
                        "recommended_max": max_markers,
                    }
                )
            if (
                kind in {"cell_type", "tumor_evidence", "cancer_context"}
                and granularity
                in {"lineage", "subtype", "tissue_subtype", "tumor_type_hint", "cancer_subtype"}
                and negative_count == 0
            ):
                gaps.append(
                    {
                        "resource": resource_name,
                        "entry": name,
                        "gap": "missing_effective_negative_markers",
                        "granularity": granularity,
                        "parent": cell.parent.name if cell.parent else None,
                    }
                )

        resources[resource_name] = {
            "review_status": dict(sorted(review_status.items())),
            "marker_counts_by_granularity": {
                granularity: {
                    "n": len(values),
                    "min": min(values),
                    "max": max(values),
                    "median": sorted(values)[len(values) // 2],
                }
                for granularity, values in sorted(by_granularity.items())
                if values
            },
        }

    gap_counts = Counter(str(gap["gap"]) for gap in gaps)
    return {
        "resources": resources,
        "gap_counts": dict(sorted(gap_counts.items())),
        "gaps": gaps,
        "priority_gaps": gaps[:50],
        "n_gaps": len(gaps),
    }


def audit_marker_resources(
    *,
    known_source_ids: set[str] | None = None,
) -> dict[str, Any]:
    """Audit built-in marker TOML resources and manager views."""
    known_source_ids = known_source_ids or set(load_reference_index())
    resources: dict[str, Any] = {}
    issues: list[ResourceAuditIssue] = []

    for resource_name, manager_name in _ACTIVE_MARKER_RESOURCES.items():
        mgr = Manager(manager_name, case_sensitive=True)
        contract_issues = mgr.validate_resource_contract(known_source_ids=known_source_ids)
        for item in contract_issues:
            _add_issue(
                issues,
                resource_name,
                str(item["issue"]),
                str(item["severity"]),
                {"entry": item["entry"], "detail": item.get("detail")},
            )
        resources[resource_name] = mgr.audit_summary(include_views=False)

    alias_data = load_toml("marker_aliases.toml")
    alias_source_ids = _source_ids_from(alias_data)
    unknown_alias_sources = sorted(alias_source_ids.difference(known_source_ids))
    if unknown_alias_sources:
        _add_issue(
            issues,
            "marker_aliases",
            "unknown_source_ids",
            "error",
            unknown_alias_sources,
        )
    for item in alias_data.get("gene_display_aliases", []):
        symbol = str(item.get("symbol", ""))
        if not _GENE_SYMBOL_RE.match(symbol):
            _add_issue(
                issues,
                "marker_aliases",
                "invalid_gene_symbol",
                "error",
                symbol,
            )

    view_summary = _audit_marker_views(issues)
    return {
        "resources": resources,
        "views": view_summary,
        "issues": [issue.as_dict() for issue in issues],
    }


def _audit_marker_views(issues: list[ResourceAuditIssue]) -> dict[str, int]:
    views = {
        "global_annotation": get_marker_manager(species="human", view="global_annotation"),
        "state_annotation": get_marker_manager(species="human", view="state_annotation"),
        "program_scoring": get_marker_manager(species="human", view="program_scoring"),
        "tumor_interpretation": get_marker_manager(
            species="human",
            cancer_type="Lung Cancer",
            view="tumor_interpretation",
        ),
        "doublet_detection": get_marker_manager(species="human", view="doublet_detection"),
        "plotting": get_marker_manager(species="human", view="plotting"),
    }
    view_summary = {name: len(mgr.CELLS) for name, mgr in views.items()}

    if any(
        cell.metadata.get("kind") in {"state", "functional_program", "tumor_evidence", "cancer_context"}
        for cell in views["global_annotation"].CELLS.values()
    ):
        _add_issue(
            issues,
            "marker_views",
            "global_annotation_contaminated",
            "error",
        )
    if "T cells" in views["program_scoring"].CELLS:
        _add_issue(
            issues,
            "marker_views",
            "program_scoring_contains_identity_marker",
            "error",
            "T cells",
        )
    if "LUSC" in views["global_annotation"].CELLS:
        _add_issue(
            issues,
            "marker_views",
            "tumor_entry_in_global_annotation",
            "error",
            "LUSC",
        )
    for view_name, count in view_summary.items():
        if count == 0:
            _add_issue(issues, "marker_views", "empty_view", "error", view_name)
    return view_summary


def audit_geneset_resources(
    *,
    known_source_ids: set[str] | None = None,
) -> dict[str, Any]:
    """Audit built-in geneset JSON resources used for scoring/enrichment."""
    known_source_ids = known_source_ids or set(load_reference_index())
    resources: dict[str, Any] = {}
    issues: list[ResourceAuditIssue] = []

    for resource_name, (species, geneset_name) in _ACTIVE_GENESET_RESOURCES.items():
        data = load_json(f"{resource_name}.json")
        payload = data.get(species, data)
        if not isinstance(payload, dict):
            _add_issue(issues, resource_name, "invalid_payload", "error")
            continue

        mgr = load_gene_set_manager(species=species, name=geneset_name, kind="geneset")
        categories = payload.get("_categories", {})
        missing_category_entries: dict[str, list[str]] = {}
        if isinstance(categories, dict):
            for category, names in categories.items():
                missing = [name for name in names if name not in mgr.CELLS]
                if missing:
                    missing_category_entries[str(category)] = missing
        if missing_category_entries:
            _add_issue(
                issues,
                resource_name,
                "category_references_missing_genesets",
                "error",
                missing_category_entries,
            )

        unknown_sources = sorted(_source_ids_from(data).difference(known_source_ids))
        if unknown_sources:
            _add_issue(
                issues,
                resource_name,
                "unknown_source_ids",
                "error",
                unknown_sources,
            )

        missing_usage = [
            name
            for name, cell in mgr.CELLS.items()
            if "use_for" not in cell.metadata or "not_for" not in cell.metadata
        ]
        if missing_usage:
            _add_issue(
                issues,
                resource_name,
                "missing_usage_metadata",
                "error",
                missing_usage[:20],
            )

        invalid_genes = {
            name: [gene for gene in cell.markers if not _GENE_SYMBOL_RE.match(gene)]
            for name, cell in mgr.CELLS.items()
        }
        invalid_genes = {name: genes for name, genes in invalid_genes.items() if genes}
        if invalid_genes:
            _add_issue(
                issues,
                resource_name,
                "invalid_gene_symbols",
                "error",
                invalid_genes,
            )

        resources[resource_name] = {
            "n_entries": len(mgr.CELLS),
            "categories": sorted(categories) if isinstance(categories, dict) else [],
        }

    return {
        "resources": resources,
        "issues": [issue.as_dict() for issue in issues],
    }


def audit_curation_index(
    *,
    known_source_ids: set[str] | None = None,
    expected_batches: int = 12,
) -> dict[str, Any]:
    """Audit the curated literature queue built from marker curation batches."""
    known_source_ids = known_source_ids or set(load_reference_index())
    issues: list[ResourceAuditIssue] = []
    try:
        rows = load_marker_curation_literature_index()
    except FileNotFoundError as exc:
        _add_issue(
            issues,
            "marker_curation_literature_index",
            "curation_index_missing",
            "warning",
            str(exc),
        )
        return {
            "n_rows": 0,
            "n_batches": 0,
            "status_counts": {},
            "issues": [issue.as_dict() for issue in issues],
        }

    batch_ids = {str(row.get("batch_id")) for row in rows if row.get("batch_id")}
    statuses = Counter(
        str(row.get("curation_status") or row.get("status") or "missing")
        for row in rows
    )
    utility_counts: Counter[str] = Counter()
    priority_counts: Counter[str] = Counter()
    missing_triage_fields: list[str] = []
    for row in rows:
        row_id = f"{row.get('batch_id', '?')}.{row.get('source_number', '?')}"
        if not row.get("resource_utility") or not row.get("target_resources"):
            missing_triage_fields.append(row_id)
            triage = classify_literature_resource_utility(row)
        else:
            triage = row
        for utility in triage.get("resource_utility", []):
            utility_counts[str(utility)] += 1
        priority_counts[str(triage.get("curation_priority", "missing"))] += 1

    if len(batch_ids) != expected_batches:
        _add_issue(
            issues,
            "marker_curation_literature_index",
            "missing_batch_coverage",
            "error",
            sorted(batch_ids),
        )

    unknown_sources = sorted(_source_ids_from(rows).difference(known_source_ids))
    if unknown_sources:
        _add_issue(
            issues,
            "marker_curation_literature_index",
            "unknown_source_ids",
            "error",
            unknown_sources,
        )

    queued = [
        row
        for row in rows
        if (row.get("curation_status") or row.get("status")) == "queued"
    ]
    if queued:
        _add_issue(
            issues,
            "marker_curation_literature_index",
            "queued_literature_requires_review",
            "warning",
            {
                "count": len(queued),
                "titles": [str(row.get("title", "")) for row in queued],
            },
        )
    if missing_triage_fields:
        _add_issue(
            issues,
            "marker_curation_literature_index",
            "missing_resource_utility_triage",
            "warning",
            {
                "count": len(missing_triage_fields),
                "examples": missing_triage_fields[:20],
            },
        )

    return {
        "n_rows": len(rows),
        "n_batches": len(batch_ids),
        "status_counts": dict(sorted(statuses.items())),
        "resource_utility_counts": dict(sorted(utility_counts.items())),
        "priority_counts": dict(sorted(priority_counts.items())),
        "issues": [issue.as_dict() for issue in issues],
    }


def audit_resource_manifest() -> dict[str, Any]:
    """Audit active manifest entries and package resource availability."""
    manifest = load_resource_manifest()
    issues: list[ResourceAuditIssue] = []
    active_resources = [
        resource
        for resource in manifest.get("resources", [])
        if isinstance(resource, dict) and resource.get("status") == "active"
    ]

    for resource in active_resources:
        name = str(resource.get("name", ""))
        filename = _resource_filename(name)
        if not resource_exists(filename):
            _add_issue(
                issues,
                "marker_resource_manifest",
                "active_resource_file_missing",
                "error",
                filename,
            )
        if not resource.get("workflow"):
            _add_issue(
                issues,
                "marker_resource_manifest",
                "missing_workflow_metadata",
                "error",
                name,
            )

    return {
        "version": manifest.get("metadata", {}).get("version"),
        "active_resources": [str(resource.get("name")) for resource in active_resources],
        "issues": [issue.as_dict() for issue in issues],
    }


def build_resource_trust_report(*, strict: bool = False) -> dict[str, Any]:
    """
    Build a structured trust report for resources used across scLucid modules.

    Non-strict mode permits review-queue warnings. Strict mode promotes warnings
    to a failing status so release workflows can demand full literature closure.
    """
    known_source_ids = set(load_reference_index())
    sections = {
        "manifest": audit_resource_manifest(),
        "markers": audit_marker_resources(known_source_ids=known_source_ids),
        "marker_quality": audit_marker_entry_quality(),
        "genesets": audit_geneset_resources(known_source_ids=known_source_ids),
        "curation_index": audit_curation_index(known_source_ids=known_source_ids),
    }
    issues = [
        issue
        for section in sections.values()
        for issue in section.get("issues", [])
    ]
    n_errors = sum(1 for issue in issues if issue.get("severity") == "error")
    n_warnings = sum(1 for issue in issues if issue.get("severity") == "warning")
    status = "pass"
    if n_errors or (strict and n_warnings):
        status = "fail"
    elif n_warnings:
        status = "warn"

    return {
        "status": status,
        "strict": strict,
        "n_errors": n_errors,
        "n_warnings": n_warnings,
        "known_source_ids": len(known_source_ids),
        "sections": sections,
        "issues": issues,
    }


def assert_trusted_resources(*, strict: bool = False) -> dict[str, Any]:
    """Raise ``AssertionError`` when the resource trust report is not acceptable."""
    report = build_resource_trust_report(strict=strict)
    if report["status"] == "fail":
        preview = report["issues"][:5]
        raise AssertionError(f"scLucid resource trust audit failed: {preview}")
    return report
