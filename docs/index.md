# scLucid Documentation

scLucid is a diagnostic-first, audit-ready Python framework for tumor
single-cell interpretation. Its documentation is organized around stable user
workflows, reviewable data contracts, and implementation plans that are kept
separate from historical design notes.

## Start Here

- [Installation](user/installation.md)
- [Quickstart](user/quickstart.md)
- [Usage Layers](user/usage_layers.md)
- [Module Features And Stage Plan](user/module_features_and_plan.md)

## Current Documentation Contract

Use these layers in order when documents disagree:

1. Code, tests, and executed notebooks define what is implemented.
2. API pages and maintained workflow guides define the public user-facing
   contract.
3. [Module Features And Stage Plan](user/module_features_and_plan.md) and
   [Current Implementation Policy](CURRENT_IMPLEMENTATION_AND_DOCS_POLICY.md)
   summarize the current module design.
4. `docs/roadmap/` contains phase plans and submission-oriented execution
   goals.
5. `docs/dev/` and `docs/archive/` preserve audit history and provenance.

## Main Workflow Spine

```text
QC decision -> preprocess layer contract -> analysis inference policy -> tumor/annotation evidence
```

The core layer contract for preprocessing is:

```text
counts -> normalized -> raw -> HVG -> scaled -> PCA -> graph
```

Every mature workflow module should expose reviewer-facing summaries that make
recommendations, applied choices, confidence, review requirements, and
biological risks explicit.
