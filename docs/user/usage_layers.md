# Usage Layers

scLucid is designed as a three-layer product rather than a single API
style. The layers target different users and different moments in the
same research project.

## Layer 1: Workflow

The workflow layer is the first screen for most users. It runs the
supported QC -\> preprocessing -\> analysis path with conservative
defaults and records the decisions under `adata.uns["sclucid"]`.

Use this layer when:

- a user wants a complete first pass with minimal code
- the project needs a reproducible baseline before manual tuning
- the goal is to compare datasets with the same supported path

Primary entrypoints:

- `scLucid.plan_analysis`
- `scLucid.run_pipeline`
- `scLucid.review_run`
- `scLucid.qc.run_qc`
- `scLucid.preprocess.run_preprocessing`
- `scLucid.analysis.run_standard_analysis`

`scLucid.qc.run_standard_qc` remains available for compatibility, explicit
step/resume control, and legacy threshold-filtering behavior. It is not the
default reviewer-first QC entrypoint.

Canonical example:

- `examples/02_simple_api/qc_preprocess_review.py` for the recommended
  four-action QC/Preprocess review and apply path
- `examples/01_workflow/basic_pipeline.py` for compatibility pipeline use

Expected output:

- filtered and processed `AnnData`
- stage review summaries
- contract validation records
- optional sidecar reports and figures when `save_dir` is set
- a unified `adata.uns["sclucid"]["run_review"]` decision surface

## Layer 2: Simple API

The simple API layer exposes each workflow stage as composable steps. It
is for analysts who want to inspect a decision, replace one method, or
run a stage again with a different config while keeping the rest of the
package conventions.

Use this layer when:

- QC thresholds need explicit review before filtering
- preprocessing parameters need to be compared
- annotation evidence needs manual inspection
- a user wants more control without rewriting the full workflow

Primary entrypoints:

- `scLucid.recommend_qc_policy`
- `scLucid.apply_qc_policy`
- `scLucid.recommend_preprocess_policy`
- `scLucid.apply_preprocess_policy`
- `scLucid.qc.calculate_qc_metric`
- `scLucid.qc.run_qc_threshold_decision`
- `scLucid.qc.build_qc_decisions`
- `scLucid.qc.filter_cells`
- `scLucid.preprocess.normalize_data`
- `scLucid.preprocess.find_hvgs`
- `scLucid.preprocess.scale_data`
- `scLucid.preprocess.batch_correction`

`scLucid.qc.recommend_intelligent_qc` and
`scLucid.recommendation.run_intelligent_preprocessing` remain compatibility
or sensitivity APIs. They do not replace DecisionCard review or establish
scientific superiority.

Canonical examples:

- `examples/02_simple_api/qc_preprocess_review.py` for stage-level
  workflows with review-summary inspection
- `examples/02_simple_api/qc_step_by_step.py` and
  `examples/02_simple_api/preprocess_step_by_step.py` for low-level
  teaching examples that should be paired with manual review
  finalization before being promoted to project templates

Expected output:

- the same AnnData conventions as the workflow layer
- inspectable intermediate objects and tables
- reviewer-facing reports for the decisions the user chose manually
- a finalized review contract when the manual path becomes a handoff
  artifact

### Manual API Contract Rule

Manual API calls are first-class scLucid usage, not a fallback. The
requirement is that a manual path must not leave the package contract
behind. If a notebook or script performs QC, preprocessing, or analysis
manually, the final handoff must still write the same module-level keys
used by the workflow layer:

``` python
adata.uns["sclucid"][module]["workflow_config"]
adata.uns["sclucid"][module]["steps_executed"]
adata.uns["sclucid"][module]["review_summary"]
```

Use `scLucid.utils.finalize_manual_review_summary()` to normalize, validate,
store, and optionally export a manual stage into this contract. The advanced
notebooks remain narrative examples, but notebook-local finalizers should no
longer be copied into new projects.

## Layer 3: Advanced

The advanced layer is for real exploratory analysis where every decision
should be visible. It is usually a notebook or project script that uses
configs, review summaries, sidecar artifacts, and custom checkpoints
together.

Use this layer when:

- a manuscript workflow needs a complete audit trail
- tumor-specific assumptions need manual review
- multiple parameter choices must be compared before finalizing
- the user needs custom hooks, checkpoints, or project-specific metadata

Primary artifacts:

- `examples/03_advanced_notebooks/`
- project notebooks based on the same contracts
- golden-path scripts such as `scripts/run_pbmc_golden_path.py`

Recommended product-facing notebook sequence:

1.  `Step1A-QC_Audit.ipynb` -\> `Step1-sce_cleaned.h5ad`
2.  `Step1B-Preprocessing_Audit.ipynb` -\> `Step2-sce_preprocessed.h5ad`
3.  `Step2-Annotation_and_Malignancy.ipynb` -\>
    `Step3-sce_annotated.h5ad`
4.  `Step3-Standard_Downstream.ipynb`
5.  `Step4-Signature_and_Target_Analysis.ipynb`

Legacy unsplit notebooks may remain as references, but the split
Step1A/Step1B sequence is the canonical advanced handoff for QC and
preprocessing.

Expected output:

- final `.h5ad`
- review summaries for each stage
- figures for inspection
- machine-readable manifest
- explicit notes for user overrides and biological assumptions

## How The Layers Work Together

The three layers should not become three separate products. They should
share the same data contracts, config names, review summary envelope,
and output locations.

The frozen layer contract is available in code:

``` python
from scLucid.utils import get_api_layer_spec

print(get_api_layer_spec("workflow"))
```

A common project flow is:

1.  Run the workflow layer to get a baseline.
2.  Inspect `review_run()` and its prioritized next actions.
3.  Drop into the simple API layer for the stage that needs adjustment.
4.  Finalize the manual decisions into the same review contract.
5.  Promote the final decisions into an advanced notebook or golden-path
    script.

Examples should make this layering explicit:

- workflow examples demonstrate conservative defaults and review
  summaries
- simple API examples may expose low-level functions, but must state
  when a result is exploratory rather than publication-grade inference
- advanced downstream notebooks should use sample-level pseudobulk and
  compositional APIs for formal condition/proportion inference
- cell-level differential expression belongs to marker discovery and
  annotation evidence unless explicitly labelled as exploratory
- ambient RNA and empty-droplet examples are diagnostic-only unless an
  external correction result is registered in the QC namespace
- workflow examples should prefer one-call or stage-level workflow APIs
- simple API examples may call low-level functions, but should say
  whether they are teaching examples or contract-preserving handoff
  examples
- advanced notebooks may include rich manual review cells, but should
  keep package policy in `docs/` and use the same review-summary schema
  as the workflow layer

## Documentation Responsibilities

`docs/` should explain policy:

- which layer to start with
- which defaults are recommended
- how review summaries should be interpreted
- which features are stable versus experimental
- what contract a manual API path must write before handoff

`examples/` should show runnable usage:

- one short script per supported layer or scenario
- minimal assumptions near the top of the file
- no hidden package policy that is absent from docs
- explicit labels for canonical, teaching, legacy, or extension examples

`examples/03_advanced_notebooks/` should show full analysis narratives:

- richer intermediate plots
- parameter review sections
- real-data or project-style execution
- final outputs that can be inspected by a reviewer
- the same module maturity, compact summary, QC handoff, layer
  transition, and step-evidence contracts used by the workflow/simple
  API layers when a notebook implements QC or preprocessing manually

## Product Acceptance Bar

A layer is considered product-ready only when:

- the code path is covered by at least one runnable example
- the docs describe when to use it and when not to use it
- the output follows the scLucid AnnData contract
- failures produce actionable errors or review warnings
- a lightweight or golden-path test protects the expected behavior

An example is considered product-facing only when:

- it runs against the installed package without ad hoc path edits
- it does not imply unsupported plugin or registry behavior
- it writes or inspects the expected review contract
- it states optional heavy dependencies as opt-in enhancements, not
  defaults
