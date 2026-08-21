# Blinded QC validation protocol

scLucid cannot validate its own QC policy by treating its predictions as truth.
The locked validation therefore separates input-derived reviewer evidence from
all policy calls and original identifiers.

## Build the pack

The registered public mixology fixture is derived from the official
[`sc_mixology`](https://github.com/LuyiTian/sc_mixology) release. Its three
post-sample-QC objects are joined on their common genes while retaining
protocol and `cell_line_demuxlet` identity truth:

```bash
python validation/qc_preprocess/prepare_mixology_fixture.py \
  --source data/sincell_with_class.RData \
  --output data/public_mixology.h5ad
```

Because the source objects were already sample-QC filtered, mixology supports
controlled identity and preprocessing-fidelity endpoints; it is not treated as
low-quality-cell truth.

```bash
python validation/qc_preprocess/build_blinded_truth_pack.py \
  --output-dir validation_outputs/current/qc_truth_pack
```

The builder creates two directories:

- `reviewer/`: the only directory given to expert reviewers;
- `sealed/`: source paths, original sample/cell identifiers, sampling tiers,
  source fingerprints, and the unblinding map.

The primary cell endpoint uses a uniform random sample within every library.
A separate metric-tail challenge set improves coverage of unusual cells but is
not used to estimate the locked primary endpoint. Sampling never uses a
scLucid policy or candidate call.

## Reviewer instructions

Reviewers edit only `sample_labels.tsv` and `cell_labels.tsv`. Allowed labels:

- `KEEP`: sufficiently supported as usable for the intended QC stage;
- `REMOVE`: sufficiently supported as technically unusable or a catastrophic
  library failure;
- `UNCERTAIN`: evidence is inadequate or plausible biology cannot be separated
  from technical damage.

Every labeled row requires `reviewer_id`. A short rationale is strongly
recommended for `REMOVE` and disagreement cases. A suspected doublet alone
should normally remain `UNCERTAIN`, because doublet evidence is evaluated by a
separate endpoint. Dataset aliases and case order are deterministic, while
original identifiers and scLucid outputs remain sealed.

Labels should be frozen before any reviewer sees selector or baseline calls.
For two-reviewer annotation, reconcile disagreements without exposing method
predictions and preserve the original reviewer files as immutable sidecars.

The generated review workbook presents the same anonymous evidence with
filters, frozen headers, label dropdowns, and a live count of unfinished rows.
After review, convert it into label-only frozen tables:

```bash
python validation/qc_preprocess/import_review_workbook.py \
  --workbook scLucid_blinded_QC_review_v1.xlsx \
  --pack-dir validation_outputs/current/qc_truth_pack \
  --output-dir validation_outputs/current/qc_truth_labels
```

The import is rejected if case coverage, evidence hashes, allowed labels, or
reviewer identifiers do not match the frozen pack.

## Run the locked gate

```bash
python validation/qc_preprocess/run_locked_qc_acceptance.py \
  --pack-dir validation_outputs/current/qc_truth_pack \
  --output-dir validation_outputs/current/qc_locked \
  --sample-labels validation_outputs/current/qc_truth_labels/sample_labels.tsv \
  --cell-labels validation_outputs/current/qc_truth_labels/cell_labels.tsv \
  --doublet-summary validation_outputs/current/qc_doublet/doublet_benchmark_report_summary.json
```

The runner verifies complete labels, evidence hashes, exact case coverage, and
source-file SHA-256 fingerprints before unblinding. Missing labels produce
`BLOCKED`; completed labels that do not meet the preregistered thresholds
produce `FAIL`. Neither state supports a superiority claim.

The QC gate does not replace separate experimental doublet, ambient RNA, or
preprocessing validation. `PASS` is scoped only to the frozen datasets,
reviewed cases, tasks, and registered baselines.
