# scLucid Documentation

scLucid is a diagnostic-first, audit-ready Python framework for tumor
single-cell interpretation. Its documentation is organized around stable user
workflows, reviewable data contracts, and implementation plans that are kept
separate from historical design notes.

Its scientific bet is that single-cell software should not only execute accepted
best practices. It should also help analysts make context-aware, evidence-backed
decisions: what assumptions were made, what biological risk was protected, what
claim level is justified, and what evidence would change the interpretation.

## Start Here

- [Installation](user/installation.md)
- [Quickstart](user/quickstart.md)
- [Project Context](user/project_context.md)
- [Reviewing Results](user/reviewing_results.md)
- [Parameter Profiles](user/parameter_profiles.md)
- [Scientific Reasoning Contracts](user/scientific_reasoning_contracts.md)
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
ProjectContext -> analysis plan -> QC -> preprocess -> analysis -> tumor -> unified run review
```

The core layer contract for preprocessing is:

```text
counts -> normalized -> raw -> HVG -> scaled -> PCA -> graph
```

Every mature workflow module should expose reviewer-facing summaries that make
recommendations, applied choices, confidence, review requirements, and
biological risks explicit.

## Vision To Roadmap

The stable user contract remains the current workflow spine above. The roadmap
extends that spine toward three planned capabilities:

- context-aware parameter and interpretation guidance for tumor and other
  biomedical settings
- unified evidence, claim, and inference semantics across QC, preprocessing,
  analysis, annotation, tumor interpretation, and support evidence modules
- sensitivity and validation artifacts that make workflow claims reviewable and
  falsifiable instead of merely reproducible

The documentation-only
[Scientific Reasoning Contracts](user/scientific_reasoning_contracts.md)
define the planned `ReasoningBrief -> WorkflowReview -> RunEvidence` handoff
and keep optional expert-perspective skills upstream of the scLucid product.
They do not change the current Python API.
