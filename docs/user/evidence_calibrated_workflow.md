# Evidence-Calibrated QC And Preprocessing

## Product Boundary

scLucid is a decision and execution system for real multi-sample, especially
tumor, scRNA-seq projects. It is not a black-box biological analyst and it does
not claim that a complex method is better because it is newer or more
statistical.

QC and preprocessing are currently `REVIEW`, not `CORE`. Analysis and Tumor
remain callable for compatibility, but their feature development is frozen
until the locked scientific and real-project usability gates pass.

## The Four-Action Main Path

```python
import scLucid as scl

context = scl.ProjectContext(
    dataset_type="tumor_tissue",
    sample_key="sample",
    batch_key="technical_batch",  # only when this is truly a technical batch
    condition_key="condition",
    input_provenance="filtered_counts",
)

qc_review = scl.recommend_qc_policy(adata, context=context)
qc_result = scl.apply_qc_policy(adata, qc_review.policy)

pp_review = scl.recommend_preprocess_policy(
    qc_result.adata,
    context,
    consumer="exploration",
)
pp_result = scl.apply_preprocess_policy(qc_result.adata, pp_review.policy)
```

`recommend_*` is read-only. `apply_*` validates an input fingerprint before
execution and returns `RunEvidence`; the output AnnData is available as
`result.adata`.

## DecisionCard

The first screen is deliberately compact:

- state: `READY`, `REVIEW`, or `BLOCKED`
- recommendation and reason
- affected samples, cells, and lineages
- counterfactual difference from each candidate
- uncertainty and missing evidence
- one next action

There is no decorative total quality score. Cell quality, sample quality,
doublets, ambient RNA, stress/damage, and lineage sensitivity are separate
evidence heads.

## QC Contract

QC proceeds through four reviews:

1. input and experimental provenance
2. sample/library validity
3. cell-policy comparison
4. quick-map and lineage-sensitive review

A multi-sample dataset without a trustworthy sample key is `BLOCKED`. Missing
unfiltered droplets makes ambient/cell-calling evidence `NOT_EVALUABLE`, not
passed. A sample with concordant catastrophic evidence across gene complexity,
mitochondrial fraction, and top-gene dominance is blocked before cell-level
filtering.

The candidate families are:

- protocol-profiled expert global baseline
- per-sample robust MAD baseline
- a miQC-family joint genes-by-mitochondrial model
- a SampleQC-family robust multi-sample, multi-QC-population model

The current Python joint and multi-sample implementations are explicitly
labelled sensitivity proxies, not reimplementations of miQC or SampleQC. They
cannot cast automatic `REMOVE` votes. An unreliable mixture fit is excluded
from review decisions. This is intentional until the actual backends and
blinded benchmark have been validated. miQC motivates probabilistic joint QC,
while SampleQC motivates robust multi-sample and multi-cell-type structure:

- [miQC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8415599/)
- [SampleQC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9912498/)

Doublet calls remain review-only without suitable external labels because
benchmark performance varies by dataset and doublet type:
[doublet benchmark](https://pmc.ncbi.nlm.nih.gov/articles/PMC7897250/).
Ambient correction changes expression estimates and is not treated as cell
deletion: [DecontX](https://link.springer.com/article/10.1186/s13059-020-1950-6).

## Preprocessing Contract

The execution result preserves distinct downstream semantics:

| Space | Canonical location | Intended consumers |
|---|---|---|
| counts | `layers["counts"]` | pseudobulk and count models |
| normalized_full | `layers["normalized_full"]` and `raw` | markers, programs, interpretation |
| discovery_rep | `obsm["X_pca"]` | neighbors, clustering, exploration |
| integrated_rep | named `obsm` key, only when selected | integrated graph or mapping |

The standard baseline is library-size normalization plus log1p. Analytic
Pearson residuals are eligible only for UMI discovery space and cannot replace
count or full-gene interpretation space. Literature supports their count-model
basis but also shows that simple log-based transformations can remain
competitive, so they must win empirically:

- [analytic Pearson residuals](https://link.springer.com/article/10.1186/s13059-021-02451-7)
- [transformation comparison](https://www.nature.com/articles/s41592-023-01814-1)

Feature selection compares batch-aware HVG with multinomial deviance. The
default unsupervised space never forces protected markers into the HVG union;
hypothesis-driven feature sets are separate sensitivity analyses. See
[Townes et al.](https://link.springer.com/article/10.1186/s13059-019-1861-6).

Unintegrated PCA is always the baseline. If batch is confounded with condition
or protected biology, integration is `BLOCKED`. Otherwise integration still
requires a Pareto comparison of batch removal, biology conservation, rare
population retention, and stability. scIB demonstrates that performance is
task-dependent and evaluates these trade-offs with multiple metrics:
[scIB](https://www.nature.com/articles/s41592-021-01336-8).

## Locked Validation Gate

The machine-readable gate is
`validation/qc_preprocess/acceptance_contract.json`; calculations are in
`validation/qc_preprocess/locked_acceptance.py`.

The primary QC endpoint requires no more than 2% false removal among expert
`KEEP` cells, at least 5 percentage points more low-quality recall than both
global and per-sample MAD baselines, and a sample/library-grouped bootstrap 95%
CI lower bound above zero. `UNCERTAIN` labels are excluded from the primary
binary endpoint.

The preprocessing gate requires held-out regret no greater than 5% and biology
loss no greater than 2%. A complex method cannot be selected without Pareto
improvement over the simple unintegrated baseline.

Whole-pipeline validation is required because processing steps interact; an
isolated favorable metric is insufficient. See
[pipeComp](https://doi.org/10.1186/s13059-020-02136-7).

Until every locked gate and all three real-project usability runs pass:

- QC and Preprocess stay `REVIEW`
- Analysis and Tumor feature development stays frozen
- no universal superiority statement is permitted

## Compatibility

Direct module APIs and configuration objects remain available for one minor
release. Calling `recommend_qc_policy(adata)` without `context` returns the
legacy dictionary bundle. Supplying `context` activates the new DecisionCard
path. Existing scripts and documents are not deleted until their replacement,
migration mapping, and regression run are recorded.
