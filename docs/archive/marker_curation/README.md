# Historical Marker Curation Batches

This directory contains historical marker curation notes organized by batch.
These documents record the rationale, sources, and decisions made while curating
marker sets for scLucid. They are kept for provenance and are **not** the
current implementation plan.

## Current Status

The live marker resource status is tracked in the top-level docs:

- `docs/marker_curation_literature_index.jsonl` — indexed literature queue
- `docs/marker_resource_quality_gaps.jsonl` — detected quality gaps
- `docs/marker_curation_candidates.jsonl` — candidate entries awaiting review
- `docs/MARKER_RESOURCE_CURATION.md` — current curation contract
- `docs/MARKER_RESOURCE_QUALITY_SUMMARY.md` — latest quality snapshot

## Batches

| Batch | Topic |
|-------|-------|
| `marker_curation_batch_01_cross_tissue_atlas.md` | Cross-tissue fibroblast atlas |
| `marker_curation_batch_02_t_cell_papers.md` | T cell papers |
| `marker_curation_batch_03_b_cell_papers.md` | B cell papers |
| `marker_curation_batch_04_nk_cell_papers.md` | NK cell papers |
| `marker_curation_batch_05_myeloid_papers.md` | Myeloid papers |
| `marker_curation_batch_06_dc_cell_papers.md` | Dendritic cell papers |
| `marker_curation_batch_07_monocyte_macrophage_papers.md` | Monocyte / macrophage papers |
| `marker_curation_batch_08_neutrophil_papers.md` | Neutrophil papers |
| `marker_curation_batch_09_endothelial_papers.md` | Endothelial papers |
| `marker_curation_batch_10_fibroblast_papers.md` | Fibroblast / CAF papers |
| `marker_curation_batch_11_pan_cancer_tumor_cell_papers.md` | Pan-cancer tumor cell papers |
| `marker_curation_batch_12_cancer_type_atlas_papers.md` | Cancer-type atlas papers |

## Adding New Curation

New marker curation should follow the contract in
`docs/MARKER_RESOURCE_CURATION.md` and update the live JSONL queues rather than
adding new files to this archive.
