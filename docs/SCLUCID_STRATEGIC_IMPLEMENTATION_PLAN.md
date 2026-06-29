# scLucid Strategic Implementation Plan

**Status**: Current strategic master plan  
**Updated**: 2026-06-12  
**Scope**: Product strategy, academic moat, long-term differentiation, and phased priorities for scLucid as a tumor single-cell interpretation system.

**Execution Playbook**: Phase-level execution details are maintained in `docs/roadmap/README.md`. This document defines the strategic architecture and priority system; roadmap documents define phase preparation, steps, acceptance standards, deliverables, and risk controls.

## 0. Strategic Positioning

scLucid should not become a Scanpy replacement, Seurat replacement, generic single-cell toolkit, general multi-omics platform, or OmicVerse-style large ecosystem. Its stronger and narrower role is:

> scLucid is a diagnostic-first, audit-ready, tumor-focused, evidence-driven interpretation system for tumor single-cell research.

## 1. Scientific Vision

The next decade of tumor single-cell biology will progressively move from
cataloging cell types toward interpreting cell states, cell communities,
ecosystems, and clinically meaningful phenotypes. scLucid is designed to evolve
along the same hierarchy:

> Gene -> Program -> Cell State -> Cell Community -> Ecotype -> Ecosystem -> Clinical Phenotype

| Interpretation Level | Current Coverage | Role In scLucid |
|----------------------|------------------|-----------------|
| Gene | Strong | Marker and gene-level evidence |
| Program | Active | Tumor and TME program scoring |
| Cell State | Active | Annotation, state interpretation, malignancy evidence |
| Cell Community | Emerging | Co-abundance and local ecosystem patterns |
| Ecotype | Planned / prototype | Sample-level archetypes and stability-tested patterns |
| Ecosystem | Planned | Integrated tumor microenvironment interpretation |
| Clinical Phenotype | Long-term | Response or outcome association with explicit limitations |

Short-term:

- Robust tumor cell-state interpretation.
- Audit-ready QC, preprocessing, analysis, annotation, and malignancy evidence.

Mid-term:

- Tumor ecosystem modeling.
- Stable sample-level ecotype and microenvironment archetype prototypes.

Long-term:

- Ecosystem-aware clinical interpretation.
- Patient/sample-level evidence summaries with explicit limitations.

This hierarchy is aspirational, not a claim that all levels are complete today.
The current core is `qc -> preprocess -> analysis`; the long-term academic moat
should come from connecting that core workflow to tumor state interpretation,
tumor ecosystem modeling, and a source-aware knowledge infrastructure.

## 2. Strategic Guardrails

### What scLucid Is

- A compact tumor single-cell workflow spine.
- A diagnostic and audit layer around routine analysis decisions.
- A tumor interpretation system that turns computational outputs into reviewable biological evidence.
- A framework for explicit inference semantics: exploratory cell-level, descriptive sample-level, sample-level inference, spatial evidence, and external evidence.
- A system that can consume validated bulk, spatial, clinical, R/Python, atlas, literature, and LLM evidence.

### What scLucid Is Not

- It is not a replacement for Scanpy or Seurat.
- It is not a general-purpose single-cell toolkit that must implement every method.
- It is not a general spatial transcriptomics platform.
- It is not a broad multi-omics platform.
- It is not a black-box automatic annotation engine.
- It is not an LLM-first biological claim generator.

### Strategic Rule

Every new feature must answer:

1. Does it improve tumor interpretation?
2. Does it strengthen evidence, auditability, or inference safety?
3. Can it be validated against data, literature, or a mature reference method?
4. Can it remain optional if it depends on heavy or specialized tooling?

If the answer is no, the feature should be delayed, kept as an external adapter, or rejected.

## 3. Target Architecture

```text
Core Workflow Layer
  QC -> Preprocess -> Analysis

Interpretation Layer
  Annotation -> Malignancy -> Tumor Programs -> TME States -> Therapy Signals

Tumor Ecosystem Layer
  Cell Communities -> Ecotypes -> Ecosystem Archetypes -> Patient Stratification

Knowledge Infrastructure
  Marker Manager, Gene Set Manager, Atlas Manager, Ontology Manager, Evidence Manager

Audit Layer
  Evidence Bundle, Review Summary, Inference Semantics, Reports, Limitations

Support Evidence Modules
  Bulk Evidence, Spatial Evidence, Clinical Evidence, Selected R/Python Parity

Long-Term Agent Interface Layer
  Dataset Review, Annotation Review, Tumor Interpretation, Tumor Board Reporting
```

The audit layer is not a separate optional feature. It is the shared contract that connects all other layers through `adata.uns["sclucid"]`, exported review summaries, evidence bundles, reports, and reproducibility metadata.

## 4. P0: Core Workflow, Audit Contract, And Validation

P0 is the current execution priority. It should receive most development attention until the core path is stable enough for real tumor projects and manuscript-grade benchmarks.

### Direction 1: Core Lucid Workflow

Goal:

Stabilize `qc -> preprocess -> analysis` as the compact single-cell workflow spine.

Scope:

- QC with adaptive thresholds, tumor-aware warnings, doublet evidence, and review summary.
- Preprocessing with explicit layer semantics, sparse-aware operations, HVG/PCA/neighbors/UMAP evidence, and batch-correction cautions.
- Analysis with clustering review, marker evidence, annotation evidence, consensus labels, post-hoc QC review, DE/proportion inference semantics, and optional malignancy bridge.

Acceptance standards:

- PBMC and PDAC golden paths can run reproducibly.
- Core outputs remain traceable in `adata.uns["sclucid"]`.
- Missing optional heavy dependencies do not break core import or core workflow.
- Major decisions produce reviewable evidence, not just output columns.

What not to do:

- Do not expand the core into every specialist method.
- Do not make bulk, spatial, clinical, or R parity part of mandatory execution.
- Do not over-optimize for UMAP appearance at the cost of evidence quality.

### Direction 2: Evidence And Audit Contract

Goal:

Make auditability the central product and academic differentiator.

Core objects:

- Review summary.
- Evidence bundle.
- Decision table.
- Inference semantics.
- Reproducibility manifest.
- Warning and limitation records.
- Method-specific evidence tables.

Acceptance standards:

- Every core module writes a structured review summary.
- Every major user-facing result has a source, confidence, limitation, and inference level.
- Reports can explain why a decision was made and what evidence supported it.
- Evidence structures are machine-readable enough to support future review agents.

Strategic value:

This is the layer that prevents scLucid from being perceived as another toolkit. It turns the package into a research review system.

### Direction 3: Layered Validation

Goal:

Prove scLucid is not merely runnable but reliable, interpretable, and safer in tumor single-cell workflows.

Validation levels:

1. **Execution parity**: scLucid produces reasonable Scanpy/Seurat-comparable baseline outputs.
2. **Decision quality**: recommendations and warnings are easier to inspect than ad hoc workflows.
3. **Biological concordance**: major lineages, tumor compartments, and markers match known biology.
4. **Inference safety**: outputs avoid common cell-level and sample-level overclaims.
5. **Usability and engineering**: fewer hidden parameters, clearer failure modes, reproducible artifacts.

Acceptance standards:

- README claims have validation material behind them.
- Benchmarks record where scLucid is better, equivalent, or weaker.
- At least PBMC, PDAC, and one second tumor type are used as acceptance gates.
- Validation produces figure-ready tables, not just notebooks.

## 5. P1: Tumor Interpretation, Ecosystem Modeling, And Scale

P1 is the main academic differentiation phase. It should start once the P0 core is sufficiently stable.

### Direction 4: Tumor State Interpretation

Goal:

Convert core analysis outputs into tumor-specific biological interpretation.

Scope:

- Malignant vs non-malignant evidence.
- CNV or CopyKAT/inferCNV-style evidence consumption.
- Tumor programs: proliferation, EMT, hypoxia, IFN response, antigen presentation, stress, apoptosis, stemness, metabolism.
- TME states: exhaustion, cytotoxicity, suppressive myeloid, inflammatory myeloid, CAF states, endothelial states, B/plasma states.
- Therapy signatures and response-associated programs where evidence is available.

Acceptance standards:

- Each interpretation output has input requirements, method notes, evidence sources, confidence, limitations, and inference semantics.
- Normal epithelial annotation is separated from malignant-cell interpretation.
- Tumor interpretation can consume stable analysis outputs without forcing all tumor algorithms into the core workflow.

### Direction 5: Tumor Ecosystem Modeling

Goal:

Move scLucid from cell-level annotation toward sample-level or region-level tumor ecosystem interpretation.

Definition:

Tumor ecosystem modeling aggregates malignant programs, TME states, cell type composition, CNV/malignancy confidence, optional spatial zones, optional bulk concordance, and optional clinical metadata into interpretable sample-level or region-level patterns.

Core concepts:

- Cell communities.
- Ecotypes.
- Microenvironment archetypes.
- Sample-level tumor states.
- Patient stratification.
- Clinical or therapy-associated ecosystem patterns.

Recommended architecture:

```text
Cell-level evidence
  annotations, programs, TME states, malignancy calls

Sample/region aggregation
  composition, program averages, state abundance, spatial zones

Stability and interpretation
  clustering, NMF, consensus, marker/program explanation

Ecosystem output
  ecotype assignment, confidence, signature, clinical association, review summary
```

Acceptance standards:

- Ecotype and ecosystem outputs are marked as exploratory unless validated at sample or cohort level.
- Results include stability analysis and explanation of dominant programs/states.
- At least one tumor case study demonstrates ecosystem summaries beyond standard UMAP plots.

Ecosystem acceptance datasets:

- **PDAC ecosystem**: CAF/myeloid/TME-rich case for malignant-stromal-immune ecosystem interpretation.
- **NSCLC ecosystem**: immune-infiltration and exhaustion-focused case for TME archetypes.
- **CRC ecosystem**: epithelial state, inflammation, and stromal/immune co-abundance case.

Ecosystem validation criteria:

- Ecotype stability under resampling and feature perturbation.
- Sample-level archetype reproducibility across cohorts or subsets.
- Program/state co-occurrence consistency.
- Patient-level separability when metadata support it.
- Concordance with literature or original-study interpretations.
- Optional support from spatial, bulk, or clinical evidence.

Strategic value:

This direction is the likely source of scLucid's long-term academic moat. Core workflow improvements are useful but easier to replicate; tumor ecosystem modeling plus audit evidence is much harder to copy.

### Direction 6: Engineering, Scale, And Productization

Goal:

Make scLucid reliable enough for real tumor projects, large datasets, and reproducible handoff.

Scope:

- Lightweight import and dependency hygiene.
- Optional extras for heavy methods.
- Sparse-aware computation.
- Memory and runtime benchmarks.
- Clear error messages.
- Manifest files and compact project outputs.
- HTML and markdown audit reports.
- Publication-oriented figure templates.

Acceptance standards:

- Large workflows avoid unnecessary `.toarray()`.
- Missing optional dependencies produce actionable messages and graceful fallback.
- Runtime and memory are benchmarked against sensible Scanpy baselines.
- Project outputs can be handed to another analyst and understood without hidden state.

## 6. P2: Knowledge Infrastructure, Selected Parity, Clinical Evidence, And Agent Interfaces

P2 contains long-term compounding assets. These should be designed early but implemented selectively, after the core and tumor interpretation layers are stable.

### Direction 7: Knowledge And Evidence Infrastructure

Goal:

Unify marker resources, gene sets, atlas references, literature knowledge, ontology, therapy knowledge, and LLM evidence into source-aware, versioned, auditable knowledge infrastructure.

Scope:

- Marker Manager: positive markers, negative markers, artifact markers, lineage/subtype/state/program/tumor views.
- Gene Set Manager: tumor programs, TME states, therapy signatures, pathway modules.
- Atlas Manager: cell type and tumor state references with version and context.
- Literature Evidence Manager: curated statements, citation links, scope notes.
- Ontology Manager: cell ontology, tumor ontology, tissue/cancer naming rules.
- Therapy Knowledge Manager: ICI response, resistance, targeted therapy signatures, if supported.
- LLM Evidence Manager: suggestions, summaries, and review aids, never ground truth.

Acceptance standards:

- Every knowledge item has source, version, species, tissue/cancer context, evidence level, limitations, and review date where possible.
- Knowledge resources are routable by task: annotation, program scoring, tumor interpretation, artifact review, therapy interpretation.
- LLM output is stored as evidence suggestion with provenance and limitations.

Strategic value:

This infrastructure makes scLucid interpretation durable. Without it, annotation, tumor programs, therapy signatures, ecotypes, and agent interfaces become ad hoc.

### Direction 8: Selected R/Python Parity

Goal:

Bridge mature R methods into Python only when they improve tumor interpretation and can be validated.

Priority:

- P0/P1 method evidence: scDblFinder parity, CNV/malignancy evidence, pseudobulk/sample-level DE, selected annotation/reference evidence, selected program scoring.
- Later optional evidence: communication methods, compositional analysis, selected deconvolution.
- Defer or avoid: broad trajectory parity, general spatial deep wrappers, and method ports that do not strengthen tumor interpretation.

Implementation modes:

1. Optional wrapper.
2. Clean-room reimplementation.
3. Parity adapter that consumes external results.

Acceptance standards:

- Every port or wrapper has a parity matrix.
- At least one synthetic and one real-data validation example are available.
- License boundaries are respected.
- Claims match the implementation level: wrapped, ported, approximated, or consumed.

### Direction 9: Clinical Evidence Integration

Goal:

Prepare scLucid to connect tumor ecosystem outputs with clinical context without overclaiming.

Scope:

- Patient metadata.
- Therapy response.
- Survival or outcome association.
- Clinical subtype.
- Pathology context.
- Sample-level tumor ecosystem features.

Acceptance standards:

- Clinical analysis remains optional and sample-level.
- Outputs distinguish association from prediction.
- Clinical evidence is never inferred from cell-level data alone.
- Missing clinical metadata does not weaken core workflow.

### Direction 10: Long-Term Agent Interface Layer

Goal:

Reserve a future interface layer for agents that summarize and critique evidence once the evidence contract is stable.

Potential agents:

- Dataset Review Agent: reviews QC/preprocess evidence and flags dataset risks.
- Annotation Review Agent: compares marker, reference, ontology, and artifact evidence.
- Tumor Interpretation Agent: summarizes malignancy, TME, programs, and ecosystem evidence.
- Tumor Board Report Agent: generates patient/sample-level interpretation summaries with explicit limitations.
- Hypothesis Generation Agent: proposes follow-up analyses, not final conclusions.

Guardrail:

Agents are interfaces over evidence, not primary analysis engines. They must consume evidence bundles and review summaries. They should summarize, critique, and route evidence; they must not generate unsupported biological claims.

Current priority:

Design only. Do not make agent features central until P0/P1 evidence and validation layers are mature.

## 7. Support Evidence Modules

Support layers are important, but they must not become the product center.

### Bulk Evidence

Role:

- Pseudobulk validation.
- Bulk deconvolution.
- Bulk-pseudobulk concordance.
- Sample-level DE support.

Boundary:

scLucid should not become a bulk RNA-seq platform. Bulk methods exist to strengthen tumor single-cell interpretation.

### Spatial Evidence

Role:

- Validate spatial organization of tumor programs and TME states.
- Support immune infiltration/exclusion interpretation.
- Support tumor-stroma boundary and niche evidence.
- Provide region-level evidence for ecosystem hypotheses.

Boundary:

scLucid should not become a general spatial transcriptomics platform. Spatial utilities should remain evidence modules for tumor interpretation.

Deferred:

- General spatial alignment platforms.
- Image segmentation platforms.
- Broad deep spatial method wrappers.
- Full spatial ecosystem replacement.

### Clinical Evidence

Role:

- Connect sample-level tumor ecosystem features with clinical metadata.
- Support response or outcome association when data are available.
- Enable translational interpretation reports.

Boundary:

Clinical outputs must remain association-level unless validated predictive modeling is explicitly performed with appropriate sample size and held-out testing.

## 8. Priority System

### P0: Current Priority

- Core workflow contract.
- Evidence and audit contract.
- Layered validation.
- QC/preprocess/analysis acceptance gates.
- Inference semantics.
- PBMC, PDAC, and second tumor golden paths.

### P1: Main Differentiation

- Tumor state interpretation.
- Tumor ecosystem modeling.
- Tumor-focused visualization templates.
- Engineering and scale.
- Knowledge infrastructure v1.
- Selected method parity only where it supports doublet, CNV, annotation, or program evidence.

### P2: Long-Term Assets

- Knowledge infrastructure expansion.
- Advanced R/Python parity.
- Clinical evidence integration.
- Spatial evidence expansion.
- Agent interface roadmap.
- More tumor types and larger validation cohorts.

### Deprioritize Or Avoid

- General spatial platform development.
- General multi-omics platform expansion.
- Broad trajectory/lineage parity unless required by a tumor case study.
- Cell-cell communication wrappers unless tied to tumor ecosystem evidence.
- Heavy deep spatial methods unless used as optional external evidence.
- LLM-first interpretation without evidence contract.

## 9. Milestones

### Milestone A: Core Lucid Workflow

Time estimate: 2-4 weeks.

Deliverables:

- PBMC golden path.
- PDAC golden path.
- QC/preprocess/analysis contracts.
- Audit report examples.
- README positioning and framework graphic.

Acceptance:

- Core workflow is reproducible.
- Audit summaries are complete.
- Optional heavy dependencies do not break import.

### Milestone B: Layered Validation

Time estimate: 4-8 weeks.

Deliverables:

- Scanpy baseline comparison.
- Seurat baseline comparison.
- Decision quality checklist.
- Inference safety suite.
- Biological concordance report.

Acceptance:

- At least 3 datasets complete validation.
- Results record better/equivalent/weaker relative performance.

### Milestone C: Tumor Interpretation System

Time estimate: 6-10 weeks.

Deliverables:

- Malignancy evidence table.
- Tumor program scoring.
- TME composition and state summaries.
- Therapy signature skeleton.
- Tumor-specific plots.

Acceptance:

- PDAC workflow can output a complete tumor interpretation report.
- A second tumor type demonstrates generalization.

### Milestone D: Tumor Ecosystem Modeling

Time estimate: 8-12 weeks.

Deliverables:

- Sample-level ecosystem feature matrix.
- Ecotype or microenvironment archetype prototype.
- Stability and confidence summaries.
- Patient-level visualization templates.

Acceptance:

- At least one case study produces ecosystem summaries with reviewable evidence and explicit exploratory status.

### Milestone E: Knowledge And Evidence Infrastructure

Time estimate: continuous, start after core contracts stabilize.

Deliverables:

- Source-aware marker and gene-set schema.
- Atlas/reference evidence adapters.
- Literature and ontology evidence schema.
- Therapy knowledge schema.
- LLM evidence guardrails.

Acceptance:

- Knowledge infrastructure can be routed by task and cited in review summaries.

### Milestone F: Release And Manuscript Package

Time estimate: 2-3 months after P0/P1 validation.

Deliverables:

- Release tag.
- Zenodo DOI.
- Documentation site.
- Reproducible benchmark scripts.
- Main figures and source data.
- Manuscript draft.

Acceptance:

- The work is ready for Genome Biology / Nature Communications / Nature Computational Science-level review, while being built to Nature Methods-style validation standards.

## 10. One-Sentence Strategy

scLucid should become the evidence-audited interpretation layer for tumor single-cell research: compact in workflow, conservative in inference, rich in tumor knowledge, and extensible through validated evidence modules.
