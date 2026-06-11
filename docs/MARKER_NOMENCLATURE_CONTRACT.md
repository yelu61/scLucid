# scLucid Marker Nomenclature Contract

This contract generalizes the modular naming principles from `Guidelines for T cell nomenclature`
to all marker resources without copying T cell-specific subtype hierarchies into other lineages.

## Generalizable Principles

- Prefer stable lineage or subtype names for cell identity labels.
- Treat activation, exhaustion, cytotoxicity, hypoxia, cycling, stress, EMT, interferon response,
  metabolic rewiring, and secretion as states or programs, not primary identity labels.
- Use modular labels when a term combines identity and state:
  `Lineage | subtype/state/program modifier`.
- Keep paper-specific labels, abbreviations, and protein/display aliases in `marker_aliases.toml`,
  not in marker gene lists.
- Use official gene symbols in `markers` and `negative_markers`.
- Use `source_ids`, `source_collection`, `evidence_role`, `use_for`, and `not_for` to make
  provenance and routing explicit.

## Canonical Label Patterns

```text
Compartment:
  Immune
  Stromal
  Epithelial

Lineage:
  T cells
  B cells
  Macrophages
  Endothelial cells

Subtype:
  CD8+ T
  cDC1
  Plasma
  Tip-like endothelial

State:
  CD8+ T | exhausted-like
  Macrophages | inflammatory
  Fibroblasts | matrix-remodeling

Functional program:
  T cell activation
  Antigen presentation
  ITH hypoxia
  Stromal barrier

Tumor context:
  Malignant ITH stemness state
  LUAD
  Basal-like breast cancer
```

## Non-Generalizable T Cell Specifics

- Do not apply T cell memory/exhaustion substage terms to unrelated lineages.
- Do not promote a T cell state such as `Tex` into a stable cell subtype.
- Do not force every lineage to have naive/memory/effector-like levels.

## Routing Expectations

- `cell_type`: identity evidence; can enter `global_annotation` only when broad enough and reviewed.
- `state`: compact context-dependent state evidence; never enters `global_annotation`.
- `functional_program` and `geneset`: scoring/enrichment modules; never direct identity labels.
- `tumor_evidence`, `cancer_context`, `cancer_state`: tumor interpretation only; do not contaminate
  normal global annotation.

## Reference

- `SRC0132`: Guidelines for T cell nomenclature, Nature Reviews Immunology, 2026.
