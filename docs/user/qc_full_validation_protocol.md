# Full QC Validation Protocol

## Decision boundary

scLucid QC starts at the count matrix and its provenance. FASTQ alignment,
read mapping, and upstream quantification are recorded but are outside the
implemented performance claim. Missing raw or unfiltered droplets makes cell
calling and ambient heads `NOT_EVALUABLE` for that project; it never creates a
synthetic pass.

QC is released only when every required evidence head passes on its registered
dataset/endpoint combinations. There is no aggregate quality score and no
compensation between heads. A strong doublet result, for example, cannot offset
missing ambient or damaged-cell evidence.

The machine-readable contracts are:

- `validation/dataset_evidence_registry.json`: accession, availability,
  license boundary, metadata, truth type, endpoint, and limitation;
- `validation/qc_preprocess/acceptance_contract.json`: frozen QC head list,
  dataset/endpoint combinations, thresholds, and release semantics;
- `validation/evidence_run_index.json`: canonical executed artifacts only;
- `validation_outputs/current/qc_full_gate/`: generated head-level readiness.

## Complete QC evidence matrix

| Evidence head | Required truth and datasets | Primary validation target | Preregistered acceptance |
|---|---|---|---|
| Input contract | controlled truth suite plus all three real projects | counts/provenance integrity, sample/library resolution, read-only review | counts integrity 100%; missing multi-sample key fails closed 100%; review mutations 0 |
| Assay/profile selection | SampleQC, Lin PDAC, real projects | scRNA versus snRNA safety profile and absence of one universal mt cutoff | profile provenance 100%; locked high-quality false blocks 0; universal mt cutoff prohibited |
| Cell calling | microscopy labels E-MTAB-2600/PRJEB6455 and 10x HGMM | intact and low-RNA recovery versus empty droplets | true-cell recall >=95%; low-RNA recall >=90%; empty-droplet FDR <=1% |
| Ambient | 10x HGMM raw droplets | cross-species soup reduction without identity erasure | contamination reduction >=50%; native-marker and identity loss each <=2% |
| Catastrophic sample QC | Lin PDAC, SampleQC, real projects | failed-library detection before cell-level filtering | failure recall 100%; locked high-quality false-block rate 0 |
| Damaged-cell classification | E-MTAB-2600/PRJEB6455 microscopy, SampleQC, Lin, real projects | damaged-cell recall at fixed intact-cell protection | KEEP false removal <=2%; recall gain >=5 percentage points versus both registered baselines; uncertainty excluded from the primary binary endpoint |
| Doublets | Kang/demuxlet, Cell Hashing, hashed mouse kidney, 10x HGMM | cross-dataset method choice and calibration | AUPRC regret <=5%; ECE <=0.10 if probabilities exist; no automatic deletion without independent confirmation |
| Rare/tumor guardrail | SampleQC, Lin, real projects | rare lineage and tumor-program protection | rare false removal <=2%; maximum lineage retention gap <=5% |
| Selector generalization | SampleQC, Lin, real projects | data-specific choice versus expert global and per-sample MAD | KEEP false removal <=2%; absolute REMOVE recall gain >=5 points versus each baseline; grouped-bootstrap 95% CI lower bound >0 |
| Iterative review | SampleQC, Lin, real projects | whether a temporary map adds independently supported removals | round-2 KEEP false removal <=2%; cluster position alone cannot delete; no new evidence yields `REVIEW` or stop |
| Policy execution | controlled truth suite plus all real projects | immutable recommendation and exact/repeatable apply | review mutations 0; policy/apply agreement 100%; repeat agreement 100%; retained counts preservation 100% |
| DecisionCard UX | all three real projects | correct action with minimal editable configuration | config fields reduced >=70%; critical-action errors 0; manual workarounds 0; RunEvidence completion 100% |
| Scalability | controlled suite, 10x HGMM, real projects | bounded sparse execution and deterministic decisions | no unintended densification; repeat agreement 100%; failure rate 0; runtime and peak memory reported on reference hardware |

The microscopy-labelled data are registered under
[E-MTAB-2600](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2600)
with raw reads under [PRJEB6455](https://www.ebi.ac.uk/ena/browser/view/PRJEB6455). Their external
labels are useful for intact/damaged/empty discrimination, but their legacy
capture technology and imperfect microscopy phenotype are explicit
limitations. Cell Hashing provides orthogonal HTO negative/singlet/doublet
information ([GSE108313](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108313)).
SampleQC supplies multimodal and sample-shift stress tests
([PMC9912498](https://pmc.ncbi.nlm.nih.gov/articles/PMC9912498/)). The 10x HGMM
mixture supplies species-based cell-calling, ambient, and mixed-species
doublet evidence.

## Evidence hierarchy

The evidence statuses are intentionally not interchangeable:

- `PASS` and `PASS_BASELINE` satisfy an exact registered endpoint;
- `SIMULATION_PASS_NOT_EXTERNAL` verifies mechanism or engineering behavior but
  cannot satisfy external scientific validity;
- `CONTRACT_PASS_NOT_PERFORMANCE` verifies schema/execution only;
- `NOT_EVALUABLE` records missing required input without claiming failure or
  success;
- `REVIEW`, `BLOCKED`, `NOT_RUN`, and missing artifacts all block CORE release.

The controlled truth suite contains explicit intact, low-RNA, rare, damaged,
empty, doublet, and ambient components. Its exact identity is
`counts = native_counts + ambient_counts`. It is used to catch contract,
determinism, and directionality failures early. It does not replace microscopy,
HTO/demuxlet, species-mixture, or blinded real-project truth.

## Execution order

```bash
python validation/qc/generate_controlled_qc_truth_suite.py
python validation/qc/run_controlled_qc_contract_benchmark.py
python validation/build_dataset_evidence_registry_report.py
python validation/qc/build_full_qc_validation_report.py
python validation/qc_preprocess/run_locked_qc_acceptance.py \
  --pack-dir validation_outputs/current/qc_truth_pack \
  --output-dir validation_outputs/current/qc_locked \
  --sample-labels validation_outputs/current/qc_truth_labels/sample_labels.tsv \
  --cell-labels validation_outputs/current/qc_truth_labels/cell_labels.tsv
python validation/qc_preprocess/build_maturity_gate_report.py \
  --qc validation_outputs/current/qc_locked/locked_qc_acceptance.json \
  --preprocess validation_outputs/current/preprocess_mixology/mixology_preprocess_benchmark.json \
  --ux validation_outputs/current/real_project_ux/real_project_ux_acceptance.json \
  --portfolio validation_outputs/current/dataset_registry/dataset_evidence_readiness.json \
  --output-dir validation_outputs/current/maturity_gate
```

Generated data and reports live only under `validation_outputs/current/`.
Disposable exploration belongs under `validation_outputs/work/` and is never
read as release evidence.

## Current release interpretation

Engineering PASS on the controlled suite means the QC review/apply contract is
safe enough to continue validation. It does not make QC scientific status
`CORE`. Scientific release remains blocked until every exact public and
real-project endpoint above has a canonical passing artifact. Analysis and
Tumor remain frozen while the complete QC gate is blocked.
