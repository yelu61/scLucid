# scLucid Documentation Map

This directory contains design, governance, and development documents for
scLucid. User-facing Sphinx documentation lives under `docs/source/`.

## Current Documents

| File | Purpose |
|------|---------|
| `SCLUCID_STRATEGIC_IMPLEMENTATION_PLAN.md` | Current strategic architecture, priorities, guardrails, support layers, and long-term roadmap. |
| `roadmap/README.md` | Phase-level execution playbook for building scLucid toward Nature Methods-level standards and Genome Biology / Nature Communications / Nature Computational Science submission readiness. |
| `BULK_SPATIAL_DESIGN.md` | Canonical design for `scLucid.tools.bulk` and `scLucid.tools.spatial`. |
| `OMICVERSE_FEATURE_AUDIT.md` | GPL-aware feature audit used to avoid direct OmicVerse code copying. |
| `DATA_USAGE_GUIDE.md` | Current local data fixture roles and handling rules. |
| `MARKER_RESOURCE_CURATION.md` | Marker resource organization and curation contract. |
| `MARKER_RESOURCE_QUALITY_SUMMARY.md` | Latest marker resource quality/trust summary snapshot. |
| `MARKER_NOMENCLATURE_CONTRACT.md` | Naming rules for marker-derived identities, states, programs, and tumor context. |
| `NAMING_CONVENTIONS.md` | Project code/API naming conventions. |
| `PLUGIN_DEVELOPMENT_GUIDE.md` | Extension guide for custom steps/plugins. |

## Archive

Older design drafts and marker curation batches are kept under `docs/archive/`
for provenance. They are not the current implementation plan.

| Archive Path | Contents |
|--------------|----------|
| `docs/archive/design_notes/` | Early strategy, QC, and analysis design notes superseded by current docs. |
| `docs/archive/marker_curation/` | Batch marker curation notes and historical coverage summaries. |

When adding new planning documentation, update either the strategic plan above
or the phase playbook below instead of creating another parallel roadmap.

## Phase Playbook

The `docs/roadmap/` directory expands the strategy into execution phases for
implementation and submission readiness:

- Phase 1: Core API and evidence contracts.
- Phase 2: QC evidence benchmark.
- Phase 3: Preprocess and analysis validation.
- Phase 4: Tumor interpretation case studies.
- Phase 5: Tumor ecosystem modeling.
- Phase 6: Knowledge and evidence infrastructure.
- Phase 7: Support evidence modules and R/Python parity.
- Phase 8: Release, manuscript, and submission package.
