# Data Usage Guide

This repository may contain small local datasets for smoke tests, examples, and
real-data acceptance runs. Treat files under `data/` as development fixtures, not
as packaged public data.

## Current Local Data Layout

At the time of this guide, the local workspace contains:

| Path | Role | Notes |
|------|------|-------|
| `data/pbmc3k.h5ad` | Normal baseline | PBMC reference fixture for quick workflow checks. |
| `data/pbmc3k_raw.h5ad` | Raw PBMC baseline | Useful for QC/preprocess acceptance. May be untracked locally. |
| `data/processed/pbmc3k_prepared.h5ad` | Prepared PBMC fixture | Intermediate prepared object for examples or validation. |
| `data/lin2020.pdac.h5ad` | Tumor fixture | PDAC-oriented real-data acceptance candidate. |
| `data/schlesinger2020.pdac.h5ad` | Tumor fixture | PDAC-oriented real-data acceptance candidate. |

Older references to directories such as `data/pbmc3k/`,
`data/human_LUAD_GSE131907/`, or `data/mouse_melanoma_GSE119352/` are archived
design notes and are not the current local layout.

## Intended Uses

- **PBMC baseline**
  - quick smoke tests;
  - normal-tissue QC/preprocess validation;
  - Scanpy/Seurat parity comparisons.
- **PDAC/tumor fixtures**
  - tumor-aware QC and annotation checks;
  - malignancy/TME interpretation validation;
  - real-data acceptance gates.
- **Processed fixtures**
  - regression tests for downstream workflow contracts;
  - examples that should not repeat expensive preprocessing.

## Data Handling Rules

- Do not assume local data files are present in every checkout.
- Keep tests that require real data marked as integration or guarded by file
  existence checks.
- Do not commit large new datasets without an explicit project decision.
- Prefer scripts that write manifests describing input path, shape, species,
  tissue, cancer type, and preprocessing state.
- Keep generated outputs outside tracked source paths, for example under
  `results/`, `validation_outputs/`, or `/tmp`.

## Recommended Acceptance Data Roles

| Role | Dataset Type | Purpose |
|------|--------------|---------|
| Normal baseline | PBMC | Validate default QC/preprocess/analysis behavior on a familiar non-tumor dataset. |
| First tumor baseline | PDAC | Validate tumor-aware annotation, malignancy evidence, TME composition, and audit reports. |
| Second tumor type | TBD | Test generalization beyond PDAC before making broader claims. |
| Active project notebook | Project-specific | Validate user-facing ergonomics and failure modes in a real research workflow. |

## Example Pattern

```python
from pathlib import Path

import scanpy as sc

DATA_DIR = Path("data")
pbmc_path = DATA_DIR / "pbmc3k.h5ad"

if not pbmc_path.exists():
    raise FileNotFoundError(
        "PBMC fixture is not available in this checkout. "
        "Use an integration dataset path or skip this real-data check."
    )

adata = sc.read_h5ad(pbmc_path)
```

## Related Documents

- `docs/SCLUCID_STRATEGIC_IMPLEMENTATION_PLAN.md`
- `docs/source/workflow_hardening.rst`
- `docs/source/data_contracts.rst`
- `scripts/run_pbmc_golden_path.py`
- `scripts/run_analysis_acceptance.py`
