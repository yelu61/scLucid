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

## LLM Curation Prompt

Use this prompt when extracting markers from review articles or pan-cancer atlas
papers:

```text
You are a single-cell RNA-seq marker curation assistant. Extract marker evidence
from the provided review, atlas, figure, or table, and format it for the scLucid
marker resource schema.

Inputs:
- Species:
- Tissue:
- Disease or cancer type:
- Source type: review | single_cell_atlas | pan_cancer_atlas | disease_atlas
- Source citation: title, year, DOI/PMID if available
- Text, table, or figure notes:

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
