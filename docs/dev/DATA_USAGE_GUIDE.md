# Data Usage Guide

This repository may contain small local datasets for smoke tests, examples, and
real-data acceptance runs. Treat files under `data/` as development fixtures, not
as packaged public data.

## Current Local Data Layout

At the time of this guide, the local workspace contains the eight ready h5ad
datasets listed below. The authoritative inventory, benchmark roles, field
contract, and future dataset priorities are maintained in
`data/DATASETS.md`.

| Path | Role | Notes |
|------|------|-------|
| `data/pbmc3k.h5ad` | Normal baseline | PBMC reference fixture for quick workflow checks. |
| `data/lin2020.pdac.h5ad` | Tumor fixture | PDAC-oriented real-data acceptance candidate. |
| `data/schlesinger2020.pdac.h5ad` | Tumor fixture | PDAC-oriented real-data acceptance candidate. |
| `data/zilionis2019.nsclc.h5ad` | Tumor benchmark | NSCLC tumor/blood benchmark with author cell-type annotations. |
| `data/baron2016.pancreas.h5ad` | Normal donor benchmark | Multi-donor pancreas reference for preprocess and batch-diagnostic checks. |
| `data/lee2020.crc.h5ad` | Tumor benchmark | CRC tumor/normal benchmark with author cell-type and subtype annotations. |
| `data/kang2018.pbmc.h5ad` | Doublet/stimulation benchmark | PBMC demuxlet labels for doublet evidence plus ctrl/stim condition labels. |
| `data/cellbender_tiny.h5ad` | Ambient fixture | Tiny CellBender/heart10k-derived fixture for empty-droplet diagnostics. |

Optional local intermediates may also exist, for example
`data/pbmc3k_raw.h5ad` or `data/processed/pbmc3k_prepared.h5ad`. Treat those as
generated or checkout-specific files unless they are documented in
`data/DATASETS.md`.

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
| Second tumor type | NSCLC / CRC | Test generalization beyond PDAC before making broader claims. |
| Batch diagnostic | Normal multi-donor pancreas | Test conservative integration and batch-correction recommendations. |
| Doublet ground truth | Kang PBMC demuxlet | Validate doublet evidence and singlet/doublet/ambiguous handling. |
| Ambient diagnostic | CellBender tiny | Validate empty-droplet and ambient RNA diagnostic contracts. |
| Active project notebook | Project-specific | Validate user-facing ergonomics and failure modes in a real research workflow. |

## Validation Output Roles

Current Phase 2 QC evidence is consolidated by
`validation/qc/build_qc_evidence_package.py` into
`validation_outputs/qc_evidence_package/`.

| Output | Purpose |
|---|---|
| `qc_source_data.tsv` | Harmonized QC source data: workflow/ambient contract, threshold decision quality, tumor program fidelity, doublet evidence, and tumor reviewer narrative. |
| `qc_claim_scorecard.tsv` | Claim-level status for QC auditability, tumor-aware biological fidelity, doublet calibration, ambient diagnostic contract, and dataset coverage. |
| `qc_dataset_coverage.tsv` | Dataset availability, roles, and mapped Figure 2 panels. |
| `qc_evidence_report.md` | Compact human-readable evidence summary. |

The current package supports the claim that scLucid QC has a systematic
evidence framework and clear tumor/doublet strengths. It should not be used to
claim universal QC superiority: Kang demuxlet is an external genotype reference
for donor doublets, and CellBender tiny is an ambient contract fixture rather
than a performance benchmark.

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
- `docs/user/workflow_hardening.md`
- `docs/user/data_contracts.md`
- `docs/user/validation_scaffold.md`
- `docs/validation/qc_preprocess_evidence_pilot.md`
- `validation/README.md`
- `data/DATASETS.md`
- `scripts/run_pbmc_golden_path.py`
- `scripts/run_analysis_acceptance.py`
