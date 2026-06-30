>
> **⚠️ ARCHIVED / SUPERSEDED**
>
> This document is kept for historical reference only. The current design,
> contracts, and implementation plans are maintained in the top-level `docs/`
> files and in `docs/user/`. Do not use this document as a source of truth.
> See `docs/README.md` for the documentation map.

---

# scLucid 核心价值定位与实施计划

## Historical Note

This document captures the earliest product positioning and implementation
sketches for scLucid. Many of the classes shown here (e.g.,
`IntelligentQCRecommender`, `BiologyAwareHVGSelector`,
`CancerPurityAwareAnnotation`) were exploratory design probes and do not
represent the current public API. The actual implementation is documented in
`docs/user/` and the top-level `README.md`.

Key ideas that survived into the current design:
- Data-driven, evidence-based QC recommendations.
- Conservative preprocessing that preserves tumor-relevant signals.
- Explicit inference semantics and auditable decisions under
  `adata.uns["sclucid"]`.

Ideas that were explicitly dropped or deferred:
- Project-level ambient RNA correction as part of the default path.
- Heavy R/Bioconductor execution branches in core QC/Preprocess.
- Fully automated therapy-response prediction without human review.
