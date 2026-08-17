# Frequently Asked Questions

This page collects common questions about installing, using, and
contributing to scLucid. Sections will grow as questions accumulate; see
the `quickstart` and `best_practices` pages for orientation.

## Installation

**Q: Which install extra should I use?**

For a typical analysis: `pip install "sclucid[analysis,de]"`. The full
`[all]` extra pulls in heavy optional backends (scVI, scVelo,
infercnvpy, squidpy, etc.) and is only needed if you intend to use those
tools.

**Q: Why do I see "Could not import optional tools backend"?**

scLucid degrades gracefully when an optional backend's dependencies are
missing. The warning identifies which backend is unavailable — install
the corresponding extra (e.g. `pip install "sclucid[tools]"` for
squidpy, scVelo, infercnvpy) to enable it. As of v0.1, `squidpy` is
required for the `spatial` backend; missing optional dependencies are
now skipped silently rather than raising an `ImportWarning`.

## Workflow

**Q: When should I use `run_pipeline` versus the per-stage API?**

Start with `ProjectContext` and `plan_analysis()`, then use
`run_pipeline(plan=plan)` for a conservative first pass. It propagates
the project context through QC, preprocess, and analysis and records
contract validation results under `adata.uns["sclucid"]`. Call
`review_run()` to see prioritized next actions. Drop to the stage API
(`run_qc`, `run_preprocessing`, `run_standard_analysis`) only for the
stage that needs an override. `run_standard_qc` is the compatibility and
explicit step/resume-control path.

**Q: How do I tell scLucid that my data is tumor tissue?**

Create `ProjectContext(dataset_type="tumor_tissue", cancer_type=...)`
and pass it through `plan_analysis()` to `run_pipeline()`. The QC stage
then treats high-mitochondrial tumor states as biological-risk evidence
rather than a simple automatic filter, and later stages can adapt their
assumptions accordingly. See `workflow_hardening` for the tumor
acceptance path on PDAC.

**Q: How do I know what to do after the pipeline finishes?**

Use `review = review_run(adata)`. `review.to_frame()` shows the stage,
status, recommended/applied value, rationale, evidence, next action, and
rerun scope. Resolve `BLOCKED` rows first, then `REVIEW` rows. `READY`
means the declared handoff is valid; it does not make automated labels
or exploratory statistics final biological truth.

## Reproducibility

**Q: How do I save the full configuration of a run?**

Every workflow stage writes its effective configuration to
`adata.uns["sclucid"][stage]["workflow_config"]` and the inheritance
chain to `adata.uns["sclucid"][stage]["config_lineage"]`. Save the
`.h5ad` and the stage save directories together; both round-trip.

**Q: Are seeds set automatically?**

Yes. Each workflow accepts a `random_state` parameter (default 42) that
is threaded through HVG selection, PCA, neighbor graphs, UMAP, Leiden,
and recommendation samplers.
