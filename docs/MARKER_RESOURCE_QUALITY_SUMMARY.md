# Marker Resource Quality Summary

Generated from `build_resource_trust_report()` and `audit_marker_entry_quality()`.

## Trust Status

- Trust status: `warn`
- Errors: 0
- Warnings: 1
- Known source IDs: 220

## Literature Utility Triage

- `benchmark_reference`: 5
- `geneset_scoring`: 86
- `marker_core`: 79
- `nomenclature_reference`: 10
- `reference_only`: 1
- `tissue_context`: 7
- `tumor_context`: 117
- `validation_reference`: 26

## Curation Priority

- `high`: 139
- `medium`: 2

## Marker Entry Quality Gaps

- Total quality gaps: 114
- `missing_effective_negative_markers`: 111
- `thin_marker_set`: 3

## Resource Review Status

- `marker_registry_human`: needs_review=64, reviewed=80, scaffold=26
- `marker_registry_mouse`: needs_review=58, reviewed=33, scaffold=9
- `marker_tissue_human`: needs_review=84
- `marker_tumor_human`: needs_review=42, scaffold=3

## First-Pass Priority Gaps

These are the first 50 gaps returned by the audit, not the full queue.

| Resource | Entry | Gap | Granularity | Detail |
|---|---|---|---|---|
| `marker_registry_human` | Stem Cell | `missing_effective_negative_markers` | `lineage` |  |
| `marker_registry_human` | Hematopoietic Stem Cell | `missing_effective_negative_markers` | `subtype` | parent=Stem Cell |
| `marker_registry_human` | Mesenchymal Stem Cell | `missing_effective_negative_markers` | `subtype` | parent=Stem Cell |
| `marker_registry_human` | RBC | `missing_effective_negative_markers` | `lineage` |  |
| `marker_registry_mouse` | Hepatocytes | `missing_effective_negative_markers` | `subtype` | parent=Specialized |
| `marker_registry_mouse` | Cardiomyocytes | `missing_effective_negative_markers` | `subtype` | parent=Specialized |
| `marker_registry_mouse` | Adipocytes | `missing_effective_negative_markers` | `subtype` | parent=Specialized |
| `marker_registry_mouse` | Erythrocytes | `missing_effective_negative_markers` | `subtype` | parent=Specialized |
| `marker_tissue_human` | Corneal EpC | `missing_effective_negative_markers` | `tissue_subtype` | parent=Cornea Tissue |
| `marker_tissue_human` | Basal EpC | `missing_effective_negative_markers` | `tissue_subtype` | parent=Cornea Tissue |
| `marker_tissue_human` | Superficial EpC | `missing_effective_negative_markers` | `tissue_subtype` | parent=Cornea Tissue |
| `marker_tissue_human` | Limbal Stem Cells | `missing_effective_negative_markers` | `tissue_subtype` | parent=Cornea Tissue |
| `marker_tissue_human` | Limbal EpC | `missing_effective_negative_markers` | `tissue_subtype` | parent=Cornea Tissue |
| `marker_tissue_human` | Progenitors | `thin_marker_set` | `tissue_subtype` | markers=2; min=3 |
| `marker_tissue_human` | Progenitors | `missing_effective_negative_markers` | `tissue_subtype` | parent=Cornea Tissue |
| `marker_tissue_human` | Corneal EndoC | `missing_effective_negative_markers` | `tissue_subtype` | parent=Cornea Tissue |
| `marker_tissue_human` | Corneal Stromal | `missing_effective_negative_markers` | `tissue_subtype` | parent=Cornea Tissue |
| `marker_tissue_human` | Keratocyte | `missing_effective_negative_markers` | `tissue_subtype` | parent=Cornea Tissue |
| `marker_tissue_human` | Limbal Stromal | `thin_marker_set` | `tissue_subtype` | markers=2; min=3 |
| `marker_tissue_human` | Limbal Stromal | `missing_effective_negative_markers` | `tissue_subtype` | parent=Cornea Tissue |
| `marker_tissue_human` | Scleral Cells | `missing_effective_negative_markers` | `tissue_subtype` | parent=Cornea Tissue |
| `marker_tissue_human` | Iridocytes | `missing_effective_negative_markers` | `tissue_subtype` | parent=Cornea Tissue |
| `marker_tissue_human` | Corneal Nerve | `missing_effective_negative_markers` | `tissue_subtype` | parent=Cornea Tissue |
| `marker_tissue_human` | Endocrine cells | `missing_effective_negative_markers` | `tissue_subtype` | parent=Pancreas Tissue |
| `marker_tissue_human` | Acinar cells | `missing_effective_negative_markers` | `tissue_subtype` | parent=Pancreas Tissue |
| `marker_tissue_human` | Ductal-like1 | `missing_effective_negative_markers` | `tissue_subtype` | parent=Pancreas Tissue |
| `marker_tissue_human` | Ductal-like2 | `missing_effective_negative_markers` | `tissue_subtype` | parent=Pancreas Tissue |
| `marker_tissue_human` | Alpha cell | `missing_effective_negative_markers` | `tissue_subtype` | parent=Pancreas Tissue |
| `marker_tissue_human` | Beta cell | `missing_effective_negative_markers` | `tissue_subtype` | parent=Pancreas Tissue |
| `marker_tissue_human` | Delta cell | `missing_effective_negative_markers` | `tissue_subtype` | parent=Pancreas Tissue |
| `marker_tissue_human` | PP cell | `missing_effective_negative_markers` | `tissue_subtype` | parent=Pancreas Tissue |
| `marker_tissue_human` | AT1 cell | `missing_effective_negative_markers` | `tissue_subtype` | parent=Lung Tissue |
| `marker_tissue_human` | AT2 cell | `missing_effective_negative_markers` | `tissue_subtype` | parent=Lung Tissue |
| `marker_tissue_human` | Ciliated cell | `missing_effective_negative_markers` | `tissue_subtype` | parent=Lung Tissue |
| `marker_tissue_human` | Club cell | `missing_effective_negative_markers` | `tissue_subtype` | parent=Lung Tissue |
| `marker_tissue_human` | Lung Goblet cell | `missing_effective_negative_markers` | `tissue_subtype` | parent=Lung Tissue |
| `marker_tissue_human` | Basal cell | `missing_effective_negative_markers` | `tissue_subtype` | parent=Lung Tissue |
| `marker_tissue_human` | Pulmonary EC | `missing_effective_negative_markers` | `tissue_subtype` | parent=Lung Tissue |
| `marker_tissue_human` | Pulmonary Fibroblast | `missing_effective_negative_markers` | `tissue_subtype` | parent=Lung Tissue |
| `marker_tissue_human` | Hepatocyte | `missing_effective_negative_markers` | `tissue_subtype` | parent=Liver Tissue |
| `marker_tissue_human` | Cholangiocyte | `missing_effective_negative_markers` | `tissue_subtype` | parent=Liver Tissue |
| `marker_tissue_human` | Kupffer cell | `missing_effective_negative_markers` | `tissue_subtype` | parent=Liver Tissue |
| `marker_tissue_human` | Hepatic stellate cell | `missing_effective_negative_markers` | `tissue_subtype` | parent=Liver Tissue |
| `marker_tissue_human` | Liver EC | `missing_effective_negative_markers` | `tissue_subtype` | parent=Liver Tissue |
| `marker_tissue_human` | Proximal tubule cell | `missing_effective_negative_markers` | `tissue_subtype` | parent=Kidney Tissue |
| `marker_tissue_human` | Distal tubule cell | `missing_effective_negative_markers` | `tissue_subtype` | parent=Kidney Tissue |
| `marker_tissue_human` | Podocyte | `missing_effective_negative_markers` | `tissue_subtype` | parent=Kidney Tissue |
| `marker_tissue_human` | Collecting duct cell | `missing_effective_negative_markers` | `tissue_subtype` | parent=Kidney Tissue |
| `marker_tissue_human` | Glomerular EC | `missing_effective_negative_markers` | `tissue_subtype` | parent=Kidney Tissue |
| `marker_tissue_human` | Mesangial cell | `missing_effective_negative_markers` | `tissue_subtype` | parent=Kidney Tissue |

## Practical Next Pass

1. Resolve the 4 queued literature rows in Zotero/local references.
2. Fix the remaining thin marker sets only after source-specific review.
3. Add tissue-level negative marker anchors to `marker_tissue_human.toml` or route tissue subtypes through lineage parents.
4. Add tumor-context negative markers for cancer subtype hints where normal epithelial or immune/stromal confusion is likely.
5. Promote `needs_review` entries only after checking the source table, figure, supplement, or indexed Zotero full text.
