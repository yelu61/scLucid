>
> **⚠️ ARCHIVED / SUPERSEDED**
>
> This document is kept for historical reference only. The current design,
> contracts, and implementation plans are maintained in the top-level `docs/`
> files and in `docs/user/`. Do not use this document as a source of truth.
> See `docs/README.md` for the documentation map.

---

# scLucid 功能完善与流程设计行动计划

## Historical Note

This document records an early sprint-style roadmap. The actual execution plan
and phase boundaries have been superseded by `docs/roadmap/` and
`docs/SCLUCID_STRATEGIC_IMPLEMENTATION_PLAN.md`.

Principles from this draft that remain in force:
- QC → Preprocess → Analysis is the stable workflow spine.
- Tumor interpretation consumes analysis outputs rather than pushing heavy
  corrections back into QC/Preprocess.
- Optional heavy dependencies (ScDblFinder, ambient RNA correction, custom
  R/Bioconductor paths) remain outside the default workflow.
