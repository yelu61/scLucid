>
> **⚠️ ARCHIVED / SUPERSEDED**
>
> This document is kept for historical reference only. The current design,
> contracts, and implementation plans are maintained in the top-level `docs/`
> files and in `docs/user/`. Do not use this document as a source of truth.
> See `docs/README.md` for the documentation map.
>
> For the current QC/Preprocess policy, see `docs/user/quickstart.md`,
> `docs/user/best_practices.md`, and `docs/user/qc_preprocess_maturity.md`.

---

# Intelligent QC Implementation Summary

## Historical Note

This document was written when `IntelligentQCRecommender` was first introduced.
The module still exists, but the current default QC path is lighter and does not
require intelligent recommendations. Heavy ambient RNA correction, ScDblFinder,
and custom R/Bioconductor execution branches are intentionally outside the
default workflow.

For the current QC/Preprocess policy, see:
- `docs/user/quickstart.md`
- `docs/user/best_practices.md`
- `docs/user/qc_preprocess_maturity.md`
