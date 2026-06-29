> **⚠️ ARCHIVED / SUPERSEDED**
>
> This coverage summary is kept for provenance only. The live marker resource
> status is tracked in `docs/marker_curation_literature_index.jsonl`,
> `docs/marker_resource_quality_gaps.jsonl`, `docs/marker_curation_candidates.jsonl`,
> and `docs/MARKER_RESOURCE_CURATION.md`.

---

# Marker Curation Coverage Summary

Generated from `docs/marker_curation_batch_01-12.md`, Zotero local inventory, and
`src/scLucid/resources/references.toml`.

## Current Coverage

- Batch curation source papers indexed: 141
- Papers registered in `references.toml`: 137
- Papers matched to local Zotero items: 133
- Papers still queued for manual/Zotero follow-up: 4

## Remaining Queued Sources

| Batch | Title | DOI | Reason |
|---|---|---|---|
| 07 | Deciphering the performance of macrophages in tumour microenvironment | 10.1186/s13045-024-01559-0 | Not matched in Zotero inventory or DOI search |
| 10 | Conserved spatial subtypes and cellular neighborhoods of CAFs | 10.1016/j.ccell.2025.03.004 | DOI search pending due local API usage limit |
| 10 | Molecular features of CAF subtypes | 10.1158/1078-0432.CCR-20-4226 | DOI search pending due local API usage limit |
| 10 | Cross-tissue single-cell landscape of human monocytes and macrophages | 10.1016/j.immuni.2021.07.007 | Not matched in Zotero inventory or DOI search |

## Resource Use

- `marker_curation_literature_index.jsonl` is the machine-readable coverage queue.
- `references.toml` now contains stable `source_id` entries for all matched/registered sources.
- Marker and geneset entries should use these `source_id` values instead of free-text citations.
- Full-text review remains required before promoting `review_status` from `needs_review` to `reviewed`.

## High-Priority Follow-Up

1. Resolve the 4 remaining queued sources.
2. Use batch-specific source IDs to add `source_ids` to existing low-confidence marker entries.
3. Promote entries only after checking source tables/figures or full text in Zotero.
4. Continue expanding resources by lineage priority: fibroblast/CAF, myeloid/macrophage, neutrophil, DC, endothelial, tumor states.
