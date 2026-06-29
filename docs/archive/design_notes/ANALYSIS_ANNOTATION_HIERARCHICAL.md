>
> **⚠️ ARCHIVED / SUPERSEDED**
>
> This document is kept for historical reference only. The current design,
> contracts, and implementation plans are maintained in the top-level `docs/`
> files and in `docs/source/`. Do not use this document as a source of truth.
> See `docs/README.md` for the documentation map.

---

# Hierarchical Annotation Workflow

This note defines the scLucid annotation boundary after adopting a
major-lineage-first workflow.

## Design

Global annotation should establish conservative major lineages and decide which
lineages deserve subset refinement. It should not force all subtype or state
labels from global cluster markers alone.

Recommended flow:

1. Run global clustering and major lineage annotation.
2. Build a subset refinement plan with `build_hierarchical_annotation_plan()`.
3. For selected lineages, run `run_subset_annotation_refinement()`.
4. Review subset clusters, markers, and optional marker-manager evidence.
5. Optionally write subset cluster or marker labels back to global `adata.obs`.

## Why Subset Reprocessing

Subset refinement should start from raw-count-like input, usually
`adata.layers["counts"]`, after global QC. It should rerun normalization, log
transformation, subset-specific HVG selection, PCA, neighbors, clustering, and
marker discovery inside the lineage.

Directly reclustering a slice of the global PCA/neighbors graph is only a quick
exploration mode. It is not the mature subset annotation path because global HVGs
and PCs are optimized for separating broad lineages, not subtype-level variation.

## Output Contract

`run_subset_annotation_refinement()` returns a dictionary:

```python
{
    "T cells": adata_t_subset,
    "Myeloid": adata_myeloid_subset,
}
```

The global object stores an audit summary in:

```text
adata.uns["sclucid"]["analysis"]["annotation"]["subset_annotation_refinement"]
```

By default, subset results are not mapped back to global `adata.obs`. Set
`write_back=True` only when the subset labels have been reviewed and should
become part of the global annotation table.

## LLM Annotation Bundle

`build_llm_annotation_bundle()` can now include `lineage_key`, `subtype_key`, and
`state_key`. This keeps data-driven LLM annotation constrained by the lineage
gate and reduces the risk of naming global clusters as over-specific subtypes.
