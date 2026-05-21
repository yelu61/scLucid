# Marker Resource Curation

This document defines how marker resources are organized and extended. The goal is
to keep all marker-dependent workflows routed through `Manager` plus resource
files, while keeping biological curation auditable and easy to maintain.

## Resource Roles

- `marker_registry_human.toml` and `marker_registry_mouse.toml`: concise marker
  evidence for normal cell identity, cell state, artifact detection, and reusable
  functional programs.
- `marker_tissue_human.toml`: tissue-specific normal parenchymal and local
  subtype markers.
- `marker_tumor_human.toml`: tumor-context evidence, tumor type hints, malignant
  programs, cancer states, and diploid reference anchors.
- `genesets_*.json` and `*.gmt`: gene-set scoring and enrichment resources. These
  may overlap TOML files at the gene level, but they are not annotation marker
  registries.

## Registry Sections

Core registries use explicit top-level sections:

- `compartment`: broad compartments such as immune, stromal, epithelial, neural.
- `lineage`: reusable broad lineages that are not tied to one compartment tree.
- `subtype`: globally reusable subtype entries that should not live under a
  single lineage hierarchy.
- `state`: lineage-aware cell states; do not use these as primary global labels.
- `artifact`: QC, stress, ambient RNA, mitochondrial, ribosomal, or hemoglobin
  signatures that should warn annotation rather than define lineage.
- `functional_program`: module-scoring programs used for biological
  interpretation, not direct global cell type annotation.

Tumor resources use explicit top-level sections:

- `epithelial_support`: epithelial/tumor identity evidence and negative markers.
- `malignancy_program`: pan-cancer or tumor-context functional programs.
- `tumor_type_hint`: tissue or cancer-type marker hints.
- `cancer_state`: tumor cell states such as cycling, EMT/plasticity, hypoxia.
- `reference_anchor`: immune/stromal diploid reference anchors for CNV-assisted
  malignancy interpretation.

## Metadata Contract

Each entry should include or inherit the following fields when possible:

- `kind`: `cell_type`, `state`, `artifact`, `functional_program`,
  `tumor_evidence`, `cancer_context`, or `cancer_state`.
- `granularity`: `compartment`, `lineage`, `subtype`, `state`, `artifact`,
  `program`, `tumor_type_hint`, `epithelial_support`, or `reference_anchor`.
- `compartment`, `lineage`, `tissue`, `disease`, `cancer_type`: biological scope.
- `scope`: `all`, `lineage_restricted`, `tissue_specific`, `tumor_context`, or
  `cancer_type_specific`.
- `applies_to`: compatible lineages, tissues, or tumor types.
- `use_for_global_annotation`, `use_for_state_annotation`,
  `use_for_malignancy_interpretation`: routing flags used by `Manager` views.
- `evidence_tier`: `seed`, `curated_review`, `atlas_supported`, `consensus`, or
  `deprecated`.
- `source_type`: `expert_seed`, `review`, `single_cell_atlas`,
  `pan_cancer_atlas`, `disease_atlas`, or `literature`.
- `review_status`: `scaffold`, `needs_review`, `reviewed`, or `conflict`.

## Curation Rules

- Use official gene symbols only. Put protein aliases, antibodies, serum markers,
  and assay-only markers in metadata fields such as `alias_markers` or
  `excluded_non_gene_markers`.
- Separate cell identity from state and program evidence. For example, `T cells`
  is identity, `T cell exhaustion-like` is state, and `Cytotoxicity` is a
  functional program.
- Do not infer malignancy from epithelial identity alone. Malignancy requires
  tumor context, CNV evidence, malignant reference support, or cancer-specific
  programs.
- Treat ribosomal, mitochondrial, hemoglobin, dissociation-stress, and ambient RNA
  signatures as artifact/QC evidence unless a tissue context supports a biological
  interpretation.
- Prefer small, specific marker lists for TOML annotation entries. Put broad
  pathway modules in `genesets_*.json` or GMT resources.

## Curation Quality Standard

Marker entries should be curated as compact pieces of annotation evidence, not as
large pathway signatures. A useful TOML marker entry should satisfy most of the
following criteria:

- The marker set supports one biological claim: compartment, lineage, subtype,
  state, artifact, program, tumor context, or reference anchor.
- The positive markers are reasonably specific in single-cell RNA-seq data and
  are expected to be detectable at the RNA level.
- The entry includes negative markers when the label is often confused with
  another lineage, compartment, or artifact.
- The entry records source provenance, review status, scope, routing flags, and
  biological caveats.
- The entry can be used by `Manager.get_view(...)` without manual interpretation
  of the source paper.

Do not add an entry directly to the reviewed resource when:

- the source only reports protein-level, antibody, serum, pathology, or bulk
  tissue markers without RNA-level support;
- the gene list mostly contains housekeeping, ribosomal, mitochondrial,
  hemoglobin, interferon, cell-cycle, heat-shock, or dissociation-stress genes
  and the intended use is primary identity annotation;
- the source describes a broad pathway, hallmark, or score rather than a concise
  identity marker set;
- the marker set is tumor-type-specific but lacks enough context to separate
  malignant identity from normal tissue lineage;
- the LLM or curator cannot identify where the marker evidence appears in the
  source.

Use these practical list sizes unless there is a clear reason to deviate:

- `compartment` and broad `lineage`: 3-12 positive markers.
- `subtype` and tissue-specific subtype: 3-10 positive markers.
- `state`: 3-15 markers, with `use_for_global_annotation = false`.
- `artifact`: 3-20 markers, with annotation-warning semantics.
- `functional_program` and `malignancy_program`: 5-50 genes; if larger, move the
  entry to a gene-set resource.
- `reference_anchor`: compact immune/stromal markers with strong negative
  epithelial/tumor markers when possible.

## Evidence Tiers and Review Status

Use `evidence_tier` to describe the biological support and `review_status` to
describe whether scLucid maintainers have accepted the entry:

- `seed`: useful starter knowledge or expert seed; not yet literature-hardened.
- `curated_review`: extracted from a review and checked by a curator.
- `atlas_supported`: supported by one well-described single-cell atlas.
- `consensus`: repeatedly supported across multiple reviews, atlases, or
  datasets and stable across contexts.
- `deprecated`: kept only for compatibility or historical tracking.

Recommended `review_status` transitions:

- `scaffold`: automatically drafted, not yet inspected.
- `needs_review`: plausible but requires manual source and biology review.
- `reviewed`: curator-accepted for the stated scope and routing flags.
- `conflict`: biologically plausible but context-dependent or contradicted by
  another source; keep caveats in `notes`.

For pan-cancer atlas papers, default to `atlas_supported` and
`needs_review`. Promote to `consensus` only after the same marker relationship is
supported by multiple independent sources or is already canonical.

## Conflict and Specificity Policy

Many useful markers are not exclusive. Curate them as evidence with scope rather
than pretending they are universal identifiers.

- If a marker appears in multiple related lineages, keep it as a broad marker and
  add more specific subtype markers at the child level.
- If a marker is useful only after lineage restriction, set
  `scope = "lineage_restricted"` and fill `applies_to`.
- If a marker is useful only in one tissue or disease, route it to
  `marker_tissue_*` or `marker_tumor_*` rather than the global registry.
- If a marker frequently causes false-positive calls, add it to
  `negative_markers` for the confused label or document the caveat in `notes`.
- If a marker represents activation, exhaustion, cycling, hypoxia, EMT,
  interferon response, stress, or ribosomal/mitochondrial content, classify it as
  `state`, `artifact`, or `functional_program` rather than primary identity.

## Recommended Curation Workflow

1. Extract candidate marker evidence from a review, atlas table, figure legend,
   supplement, or author-provided annotation table.
2. Normalize names to official gene symbols for the target species.
3. Classify each candidate into one resource section and one biological role.
4. Separate identity markers from state/program/artifact markers.
5. Add scope, routing flags, source provenance, evidence tier, review status, and
   notes.
6. Run resource validation and `Manager` view tests before committing the marker
   change.
7. Promote entries from `needs_review` to `reviewed` only after manual inspection
   against the source and basic marker specificity checks.

## LLM Curation Prompt

Use this prompt when extracting markers from review articles or pan-cancer atlas
papers:

```text
You are a single-cell RNA-seq marker curation assistant for scLucid. Extract
marker evidence from the provided review, atlas, figure, table, or supplement and
format it for the scLucid marker resource schema. Your job is evidence
extraction and conservative classification, not biological guessing.

Inputs:
- Species:
- Tissue:
- Disease or cancer type:
- Source type: review | single_cell_atlas | pan_cancer_atlas | disease_atlas
- Source citation: title, year, DOI/PMID if available
- Text, table, or figure notes:
- Optional existing scLucid labels to merge with:

Rules:
1. Use only official gene symbols in `markers` and `negative_markers`.
2. Do not store protein aliases, antibody names, serum markers, or non-gene
   biomarkers as marker genes. Put them in metadata.
3. Classify each entry as exactly one of:
   compartment, lineage, subtype, state, artifact, functional_program,
   epithelial_support, malignancy_program, tumor_type_hint, cancer_state,
   reference_anchor.
4. Separate cell identity from cell state and functional programs.
5. If an entry is lineage-, tissue-, disease-, or cancer-type-specific, fill
   metadata fields `scope`, `applies_to`, `tissue`, `disease`, or `cancer_type`.
6. Do not interpret epithelial markers alone as malignant markers.
7. Flag broad, low-specificity, housekeeping, ribosomal, mitochondrial,
   hemoglobin, dissociation-stress, and ambient RNA signatures as artifact or
   program evidence, not primary cell identity.
8. Mark uncertain entries with `review_status = "needs_review"` and explain the
   caveat in `notes`.
9. Do not invent markers that are not present in the source text, table, figure,
   supplement, or provided notes.
10. Prefer concise marker sets. If a source gives a long pathway or signature,
    classify it as `functional_program` or `malignancy_program`, and set
    `recommended_resource = "geneset"` when it should not be a TOML marker entry.
11. For tumor epithelial markers, set `use_for_malignancy_interpretation = true`
    only when the entry is tumor-context evidence; do not set
    `use_for_global_annotation = true`.
12. For cell states, artifacts, and programs, set
    `use_for_global_annotation = false`.
13. Include negative markers when the source or domain knowledge indicates a
    common confusion, but mark them as `inferred_negative_markers = true` in
    metadata if they were not explicitly stated by the source.
14. Use `confidence = "high" | "medium" | "low"`:
    - high: canonical marker relationship or directly supported by a clear
      single-cell figure/table;
    - medium: plausible and source-supported but context-dependent;
    - low: weak, indirect, ambiguous, or needs human review.

Return JSON records with this structure:
[
  {
    "section": "",
    "name": "",
    "markers": [],
    "negative_markers": [],
    "metadata": {
      "kind": "",
      "granularity": "",
      "category": "",
      "species": "",
      "compartment": "",
      "lineage": "",
      "tissue": [],
      "disease": [],
      "cancer_type": [],
      "scope": "",
      "applies_to": [],
      "use_for_global_annotation": false,
      "use_for_state_annotation": false,
      "use_for_malignancy_interpretation": false,
      "evidence_tier": "",
      "source_type": "",
      "review_status": "",
      "confidence": "",
      "recommended_resource": "marker_toml",
      "inferred_negative_markers": false,
      "source": {
        "title": "",
        "year": "",
        "doi": "",
        "pmid": "",
        "figure_or_table": ""
      },
      "notes": ""
    }
  }
]
```

## LLM Review Prompt for Existing Entries

Use this second prompt after records have already been drafted. It is designed to
find routing mistakes, over-broad markers, and identity/state mixing.

```text
You are reviewing proposed scLucid marker resource entries. Evaluate whether each
entry is safe to use for the stated resource section and Manager routing flags.

For each entry, return:
- decision: accept | revise | reject
- main_issue: one sentence
- corrected_section, corrected_kind, corrected_granularity
- marker_edits: markers_to_keep, markers_to_remove, markers_to_move_to_geneset,
  suggested_negative_markers
- routing_edits: use_for_global_annotation, use_for_state_annotation,
  use_for_malignancy_interpretation
- evidence_tier and review_status recommendation
- notes for human curator

Review criteria:
1. Primary identity labels must not be dominated by state, cycling, stress,
   ribosomal, mitochondrial, hemoglobin, or broad pathway genes.
2. Tumor epithelial evidence must not be treated as malignancy by itself.
3. Tissue-specific markers should not be promoted to global annotation unless
   they are stable across tissues.
4. Pan-cancer programs belong in tumor interpretation or gene-set scoring, not
   broad cell identity annotation.
5. Negative markers should capture common confusions without excluding true
   biological subtypes.
6. If source provenance is missing or ambiguous, keep `review_status =
   "needs_review"` or reject the entry.
```
