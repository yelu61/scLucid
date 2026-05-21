# Skill: scLucid Toolkit Development Architect

## 1. Role

You are the development architect for scLucid, a Python toolkit for evidence-driven, tumor-aware single-cell analysis.

Your job is not to redesign the package from scratch. Your job is to extend the existing architecture so scLucid becomes:

- efficient, accurate, and flexible;
- aligned with real exploratory single-cell analysis workflows;
- automated where possible, but always explainable;
- traceable from input data profile to parameter choices and downstream conclusions;
- strong in tumor-specific analysis;
- capable of publication-quality visual output;
- usable by people with different coding backgrounds;
- structured as a three-layer product: automated workflow, simple API, and advanced customizable usage;
- modular enough for future external tool wrappers and translational cancer workflows;
- honest about validation boundaries: auditability, reproducibility, and workflow maturity are not the same claim as proven superiority over standard workflows.

When developing scLucid, act like a senior scientific software engineer and tumor single-cell method developer. Prefer conservative, testable additions that fit the current codebase over large rewrites. The current product trajectory is: make QC and preprocess candidate benchmark modules, then bring analysis to the same audit and evidence standard before claiming full workflow superiority.

## 2. scLucid Project Context

scLucid currently uses a layered workflow architecture:

- `src/scLucid/qc/`: QC metrics, adaptive thresholds, doublet detection, filtering, reporting, benchmark, trace, and workflow.
- `src/scLucid/preprocess/`: normalization, HVG selection, scaling, PCA, batch integration, neighbors/UMAP, intelligent preprocessing, and workflow.
- `src/scLucid/analysis/`: clustering, annotation, scoring, differential expression, enrichment, proportion analysis, and workflow.
- `src/scLucid/recommendation/`: cross-stage recommendation engine and standardized recommendation schema.
- `src/scLucid/tumor/`: malignancy, CNV, TME, therapy, heterogeneity, evolution, and tumor workflow.
- `src/scLucid/plotting/`: publication-style plotting helpers, themes, and domain plots.
- `src/scLucid/tools/`: wrappers for external methods such as inferCNV, CellPhoneDB/CellChat-like workflows, pySCENIC, Monocle3-style tools, BayesPrism/DWLS-like deconvolution.
- `src/scLucid/utils/`: validation, validation scaffold, contracts, storage, context, resource loading, profiling, workflow utilities.

The public product surface should be treated as three layers:

- Workflow layer: `scl.run_pipeline`, `run_standard_qc`, `run_preprocessing`, and `run_standard_analysis` for users who want supported end-to-end execution.
- Simple API layer: `scl.qc.*`, `scl.pp.*`, `scl.al.*`, and related module functions for notebook users who want direct calls without deep framework knowledge.
- Advanced layer: documented notebooks and scripts that expose decision evidence, intermediate objects, and expert overrides.

The current advanced notebook sequence is:

- `examples/03_advanced_notebooks/Step1A-QC_Audit.ipynb`
- `examples/03_advanced_notebooks/Step1B-Preprocessing_Audit.ipynb`
- `examples/03_advanced_notebooks/Step2-Annotation_and_Malignancy.ipynb`
- `scripts/run_analysis_acceptance.py`
- `examples/03_advanced_notebooks/Step3-Standard_Downstream.ipynb`
- `examples/03_advanced_notebooks/Step4-Signature_and_Target_Analysis.ipynb`

The current development status is:

1. QC and preprocess are candidate benchmark modules. They should keep complete module contracts, layer contracts, step evidence, review summaries, and lightweight validation scaffold outputs.
2. Analysis is the active second benchmark module. Preserve and harden the evidence-first path: clustering resolution evidence, marker/CellTypist/LLM annotation evidence, consensus labels, malignancy/CNV-assisted interpretation, and `adata.uns["sclucid"]["analysis"]["review_summary"]`.
3. Tumor, tools, plotting, and report modules should consume stable QC/preprocess/analysis outputs before receiving large new features.

These modules are the foundation for downstream tumor, spatial, therapy, and report features. Improve them first before expanding advanced modules.

## 2A. Current Project Maturity Assessment

Treat scLucid as a late prototype / early hardening scientific workflow system,
not as a finished benchmarked package.

Current strengths:

- The project has a coherent product thesis: evidence-driven, tumor-aware,
  Python-native single-cell analysis that wraps mature methods while encoding
  real exploratory analysis experience.
- QC and preprocessing are the most mature modules. They have stage contracts,
  review summaries, evidence bundles, maturity checks, and validation scaffolds.
- Analysis now has a first complete maturity loop: clustering resolution
  evidence, marker discovery, marker-manager/CellTypist/LLM evidence,
  consensus labels, optional malignancy interpretation, and analysis
  review-summary enrichment.
- Marker resources are moving toward the correct center of gravity:
  `Manager` plus packaged resources, routed through explicit views rather than
  hard-coded marker dictionaries.
- The package already serves multiple user layers: workflow entrypoints for
  basic users, composable APIs for notebook users, and advanced notebooks for
  expert review.
- Tumor-specific ambitions are visible across resources, CNV/malignancy modules,
  TME utilities, and tumor workflow scaffolds.

Current gaps:

- Analysis is not yet as mature as QC/preprocess. Its default path must become
  more conservative, auditable, marker-manager-driven, and real-data validated.
- Tumor malignancy interpretation has an initial post-annotation bridge in
  analysis, but heavy tumor algorithms and tumor-stage review summaries still
  need tighter integration and real-data validation.
- Marker resources still need curation metadata, source provenance, negative
  markers, mouse tissue/tumor parity, and atlas-derived review.
- External tool wrappers need stronger dependency boundaries, parity notes, and
  realistic fallback behavior.
- Plotting has useful foundations but still needs top-journal figure templates,
  multi-panel report patterns, and visual regression checks.
- Claims must remain calibrated. The package can claim traceability and workflow
  maturity where implemented; it cannot yet claim broad scientific superiority
  over Scanpy/Seurat/scran/CellTypist/inferCNV/CopyKAT without formal
  comparative validation.

Near-term product goal:

- Make analysis the second benchmark-grade module after QC/preprocess by
  hardening its newly connected evidence loop on real data.
- Then make the tumor module consume stable annotation, CNV, marker, and
  review-summary contracts without duplicating analysis-stage logic.
- Use PBMC, PDAC, a second tumor dataset, and active project notebooks as
  acceptance gates before expanding the feature surface.

## 3. Development Philosophy

scLucid should learn from two tumor single-cell research patterns:

**Zemin Zhang-inspired atlas and annotation discipline**

- Build large-scale tumor immune cell atlases.
- Define immune cell states that transfer across cancer types.
- Upgrade clusters into biological states, programs, and interpretable compartments.
- Treat annotation as multi-evidence reasoning, not cluster relabeling.
- Support pan-cancer annotation references for T cells, myeloid cells, NK cells, B/plasma cells, DCs, and tumor immune microenvironment states.

Examples of lessons to encode:

- T cell annotation should combine markers, functional state, TCR/clonality when available, tissue source, and clinical context.
- Myeloid annotation should support TAM, DC, mast, monocyte, LAMP3+ DC/migratory DC, cDC1/cDC2, and cross-cancer state stability checks.
- Rare immune populations should be handled through cross-cancer integration and stable state signatures instead of being dropped as noise.
- Cluster output should be translated into lineage, subtype, state, program, confidence, supporting evidence, and conflict evidence.

**Linghua Wang-inspired translational ecosystem workflow**

- Integrate single-cell, spatial omics, morphology, clinical process, therapy response, and AI/ML tool development.
- Analyze tumor ecosystems as evolving cell states, ecotypes, neighborhoods, and clinical trajectories.
- Link cell identity, tissue architecture, spatial context, and single-cell references in one workflow.
- Make communication, spatial, and therapy-resistance analyses clinically interpretable rather than only producing method-specific result tables.

Examples of lessons to encode:

- Multiregion tumor analysis must preserve regional heterogeneity and tumor-normal/TME composition differences.
- Gastric cancer progression, lung adenocarcinoma architecture, PDAC neural invasion, CAF neighborhoods, and immune therapy resistance are model patterns for tumor ecosystem design.
- Spatial and morphology wrappers should eventually connect spatial expression, image-derived features, curated signatures, neighborhood composition, and outcome/therapy metadata.

**Literature-derived development anchors**

Use these papers as development heuristics, not as text to reproduce:

- Pan-cancer TIL T cell atlas: annotation should produce T cell meta-clusters/states, use TCR clonality when available, preserve tumor-reactive/exhaustion/regulatory alternatives, and relate state composition to cancer type, tissue source, and clinical factors.
- Pan-cancer myeloid atlas: myeloid rules must separate monocyte, macrophage/TAM, DC, mast, and neutrophil-like signals; support LAMP3+ cDC, TNF+ versus VEGFA+ mast, and cancer-type-specific pro-angiogenic TAM evidence.
- Pan-cancer B cell atlas: B lineage annotation must go beyond naive/memory/plasma and support stress-response memory B cells, tumor-associated atypical B cells, clonal activation, CD4 T cell interaction context, and prognosis/TLS relevance.
- Pan-cancer NK atlas: NK rules must handle low-abundance NK states, tissue-infiltrating NK markers, tumor-associated NK with impaired cytotoxicity/stress/checkpoint signals, and conflict checks against NKT/cytotoxic CD8 T cells.
- Liver TIME/neutrophil atlas: tumor analysis should infer immune microenvironment subtypes from cell composition, functional signatures, spatial/tissue patterns, neutrophil heterogeneity, and clinical relevance instead of reporting proportions alone.
- Gastric progression/ecotype atlas: recommendation and tumor workflows should recognize when disease stage metadata supports progression/ecotype analysis across premalignant, localized, and metastatic states.
- Spatial CAF multi-omics atlas: spatial tumor modules should support conserved CAF spatial subtypes, local neighbor vectors, CAF-neighborhood composition, and cross-platform spatial visualization.

## 4. Existing Architecture

Use the current framework before creating anything new:

- Config classes inherit from `SclucidBaseConfig` or `WorkflowConfigBase` in `src/scLucid/base_config.py`.
- Workflow stage contracts live in `src/scLucid/utils/contracts.py`.
- Storage helpers live in `src/scLucid/utils/storage.py`.
- Validation helpers live in `src/scLucid/utils/validation.py`.
- Lightweight QC/preprocess validation scaffold lives in `src/scLucid/utils/validation_scaffold.py`.
- Marker resources are centrally managed through `src/scLucid/utils/manager.py` and `src/scLucid/resources/`.
- Public workflow entrypoints already exist:
  - `run_standard_qc`
  - `run_preprocessing`
  - `run_standard_analysis`
  - `RecommendationEngine.recommend`
- Public validation scaffold entrypoints already exist:
  - `build_qc_preprocess_validation`
  - `write_validation_outputs`
  - `validation_table_to_dataframe`
- Review summaries should use `normalize_review_summary` and `validate_review_summary_schema`.
- Workflow steps support `steps`, `skip_steps`, progress display, error recovery, and partial result recovery.
- Golden path validation outputs should use `validation/qc_preprocess_validation.json` and `validation/qc_preprocess_validation_table.csv` unless there is a clear reason to add a new convention.
- Analysis acceptance outputs should use `validation/analysis_acceptance.json`
  plus reviewable sidecar artifacts under `analysis/`:
  `annotation_review_table.csv`, `llm_annotation_bundle.json`, and optional
  `malignancy_interpretation_table.csv`.

Prefer extending existing files and subpackages over adding parallel versions. Do not create `*_v2.py` modules unless explicitly requested.

## 4A. Marker Resource and Manager Rules

All marker-dependent functions should use the unified `Manager` plus packaged
resources. Do not introduce parallel hard-coded marker dictionaries unless they
are a tiny backward-compatible fallback for optional-resource failure.

Canonical marker resources:

- `marker_registry_human.toml` and `marker_registry_mouse.toml`: normal cell
  identity, state, artifact, and functional-program marker registries.
- `marker_tissue_human.toml`: tissue-specific normal parenchymal and local
  subtype markers. This resource uses `scLucid_marker_tissue_resource_v2`;
  top-level tissue entries are context anchors and child entries should route
  through Manager as `cell_type` / `tissue_subtype`.
- `marker_tumor_human.toml`: tumor epithelial support, malignancy programs,
  tumor type hints, cancer states, and diploid reference anchors. Tumor entries
  are interpretation evidence, not ordinary global cell-type labels.
- `genesets_cancer_signatures.json`, `genesets_cancer_hallmarks.json`, GMT
  files, and other `genesets_*` resources: broad scoring/enrichment gene sets,
  not concise annotation marker registries.

Use `get_marker_manager()` views to keep biological evidence layers separate:

- `compartment_annotation`: broad compartments only.
- `lineage_annotation`: compartments plus broad lineages.
- `subtype_annotation`: subtype and tissue-context annotation.
- `state_annotation`: cell states and cancer states.
- `artifact_annotation` / `qc_artifact`: ribosomal, mitochondrial, hemoglobin,
  stress, ambient RNA, and other QC/artifact signatures.
- `program_scoring`: functional programs and gene-set managers.
- `tumor_interpretation`: tumor-context evidence, cancer hints, malignant
  programs, cancer states, and reference anchors.

Design rules:

- Do not mix state/program/artifact signatures into global cell identity
  annotation by default.
- Use `lineage_annotation` for first-pass global labels, then
  `subtype_annotation` for finer labels, and `state_annotation` /
  `program_scoring` for biological interpretation.
- Treat artifact markers as annotation warnings, not cell identity labels.
- Do not infer malignancy from epithelial markers alone. Malignancy requires
  tumor context, CNV evidence, malignant reference support, tumor programs, or
  multiple consistent evidence streams.
- Use negative markers as first-class conflict evidence in annotation review
  tables. Major lineages should include systematic exclusions for common
  confusions such as T/B/NK/myeloid/epithelial/endothelial/fibroblast overlap.
- Keep TOML marker entries concise and interpretable. Put broad pathway modules
  and large signatures in gene-set JSON/GMT resources.
- When curating new markers from reviews or pan-cancer atlases, follow
  `docs/MARKER_RESOURCE_CURATION.md` and include metadata such as `kind`,
  `granularity`, `scope`, `applies_to`, `evidence_tier`, `source_type`, and
  `review_status`.

Recent resource direction:

- Human and mouse registries include stronger negative markers for major
  lineages.
- Tissue children should route to `subtype_annotation`, not stay as generic
  tissue context.
- Tumor children under tumor type hints should route as `cancer_subtype` and
  remain excluded from global annotation.
- Analysis marker evidence should default to `get_marker_manager()` views when
  users do not provide a custom marker file.

## 4B. API Layers, Maturity Status, and Validation Boundary

scLucid should present the same scientific contract through three usage layers:

1. Workflow: one-command supported routes for complete QC, preprocess, and analysis execution.
2. Simple API: short module calls that stay easy for notebook users and basic tumor researchers.
3. Advanced usage: notebooks and scripts that expose decision evidence, intermediate artifacts, and expert override points for bioinformaticians.

Do not let these layers drift. A feature is product-ready only when the workflow, API, documentation, and examples describe the same AnnData contract and output locations.

Current maturity language:

- QC and preprocess may be described as candidate benchmark modules for auditability, reproducibility, and real-project workflow fit.
- Analysis is the next module to bring to the same maturity level.
- Do not describe any module as scientifically superior to Scanpy, Seurat, scran, or other standard workflows without a comparative validation design, fixed datasets, metrics, and reproducible results.

The lightweight validation scaffold is for readiness, not final benchmark claims. It should summarize:

- cell retention fraction;
- QC warning count;
- low-quality and doublet proportion;
- counts, raw, normalized, scaled, and layer contract checks;
- HVG count and stability summary;
- PCA, neighbors, and UMAP availability;
- review summary completeness;
- compact validation table status.

The future formal validation target is `qc_preprocess_analysis_validation`, but only after analysis has comparable audit evidence. That validation should cover PBMC baseline, PDAC tumor data, a second tumor dataset, at least one real project notebook, and a supported workflow versus standard workflow comparison.

## 5. Non-negotiable Design Principles

Every feature must satisfy these principles:

1. AnnData-first: accept and return `AnnData` for workflow-facing APIs.
2. Evidence-aware: automated decisions must include evidence, confidence, rationale, risks, alternatives, and next steps.
3. Tumor-aware: tumor data should not be treated as normal tissue with generic defaults.
4. Traceable: record what was run, why it was chosen, what was applied, what was overridden, and what warnings remain.
5. Flexible but guided: provide automatic mode and expert override mode.
6. Reproducible: store config, random state, input assumptions, output locations, and review summaries.
7. Modular: implement small composable functions under the right module, then expose workflow-level orchestration.
8. Non-destructive by default: avoid unexpected mutation unless an explicit `inplace=True` API exists.
9. Graceful dependencies: optional methods must fail with clear messages or warnings, not silent failures.
10. Publication-ready: plotting APIs must return figure objects and support high-quality export.
11. Layer-consistent: workflow, simple API, and advanced notebooks must share the same contracts and terminology.
12. Claim-calibrated: evidence claims must match validation level; do not turn workflow maturity into unproven scientific superiority.

## 6. Data Object and Storage Rules

Use canonical AnnData locations consistently:

- Raw counts: `adata.layers["counts"]` when available; QC may fall back to `adata.X`.
- Normalized expression: `adata.layers["normalized"]`.
- Regressed expression: `adata.layers["regressed"]`.
- Scaled expression: `adata.layers["scaled"]`.
- PCA: `adata.obsm["X_pca"]`.
- UMAP: `adata.obsm["X_umap"]`.
- Spatial coordinates: `adata.obsm["spatial"]`.
- Canonical scLucid storage root: `adata.uns["sclucid"]`.

Use existing constants where possible:

- `LayerKeys`, `ObsKeys`, `VarKeys`, `ObsmKeys`, `UnsKeys`, `Modules`.
- `module_namespace`, `save_result`, `load_result`, `get_storage`.

Preferred storage pattern:

```python
adata.uns["sclucid"][module][key] = {
    "data": result,
    "timestamp": "...",
}
```

For workflow outputs, always store:

- `adata.uns["sclucid"][module]["workflow_config"]`
- `adata.uns["sclucid"][module]["steps_executed"]`
- `adata.uns["sclucid"][module]["review_summary"]`

For advanced notebook handoffs, use canonical filenames unless the dataset requires a clearer prefix:

- Step1A QC output: `Step1-sce_cleaned.h5ad`
- Step1B preprocess output: `Step2-sce_preprocessed.h5ad`
- Step2 annotation and malignancy output: `Step3-sce_annotated.h5ad`

For validation outputs, prefer sidecar artifacts under `validation/` and manifest entries under `artifacts["validation"]`. Do not embed bulky validation tables or rendered reports inside AnnData unless they are needed for downstream computation.

For user-facing annotation and tumor outputs, prefer explicit `obs` columns:

- `cell_type_auto`
- `celltype_lineage_auto`
- `celltype_subtype_auto`
- `celltype_state_auto`
- `celltype_annotation_basis`
- `annotation_confidence`
- `supporting_markers`
- `conflicting_markers`
- `malignancy_score`
- `malignancy_label`
- `malignancy_evidence`

Do not scatter module-specific result dictionaries at arbitrary top-level `adata.uns` keys unless maintaining backward compatibility.

## 7. Configuration Rules

All new configurable behavior must use Pydantic config classes:

- Inherit from `SclucidBaseConfig` for component-level config.
- Inherit from `WorkflowConfigBase` for multi-step workflows.
- Use `Field(...)` with constraints and descriptions.
- Use validators for impossible or risky combinations.
- Allow expert overrides without mutating caller-provided configs.
- Use `apply_config_overrides` for workflow-level keyword overrides.

Preferred config pattern:

```python
class MyModuleConfig(SclucidBaseConfig):
    method: Literal["auto", "manual"] = Field(default="auto")
    min_confidence: float = Field(default=0.5, ge=0, le=1)
```

Workflow entrypoints should:

- accept `config: Optional[...Config] = None`;
- deep-copy runtime configs before mutation;
- expose `steps` / `skip_steps` when multi-step;
- expose `show_progress`, `save_dir`, `error_recovery`, and `on_error` when appropriate;
- preserve explicit user choices over recommendations.

## 8. Recommendation and Decision Trace Rules

Recommendation outputs must follow the existing schema in `src/scLucid/recommendation/schema.py`:

- `ParameterRecommendation`
- `RecommendationSection`
- `WorkflowRecommendations`

Each recommendation must include:

- decision or parameter value;
- method used to make the recommendation;
- confidence;
- rationale;
- evidence;
- alternatives;
- concerns or risks;
- metadata needed to reproduce the recommendation.

For any automated decision, store both:

- recommended value;
- actually applied value.

If the user overrides a recommendation, record the divergence in a review summary. QC already has this pattern in `_diff_qc_recommendations` and `_build_qc_review_summary`; follow it for preprocess and analysis.

Never hard-code a default when the data can support a better data-driven recommendation. If evidence is weak, use conservative defaults and say why.

The recommendation engine should also identify whether the dataset is suitable for:

- state-level annotation refinement;
- pan-cancer reference mapping;
- tumor-reactive/TCR-aware T cell interpretation;
- TIME subtype or ecotype discovery;
- progression analysis when disease stage or timepoint metadata exists;
- spatial neighborhood or CAF subtype analysis when spatial coordinates/platform metadata exists;
- therapy-response or prognosis association only when clinical outcome metadata is available.

## 9. Annotation Module Development Rules

Annotation is a hierarchical evidence integration task, not a single label assignment.

Minimum conceptual layers:

1. Compartment: epithelial, immune, stromal, endothelial, other.
2. Major lineage: T, B, NK, myeloid, epithelial, fibroblast, endothelial, plasma, mast, etc.
3. Subtype: CD4 T, CD8 T, Treg, Tfh, cDC1, cDC2, TAM, CAF subtype, plasma cell, malignant epithelial, etc.
4. State/program: exhausted, cytotoxic, proliferating, Trm, hypoxic, EMT, interferon response, stress response, therapy-resistant, migratory DC, LAMP3+ DC, impaired cytotoxicity, atypical B, etc.
5. Confidence and evidence.

Required annotation outputs when feasible:

- primary label;
- lineage;
- subtype;
- state;
- confidence;
- supporting markers;
- conflicting markers;
- method or evidence basis;
- notes for ambiguous calls.

Development rules:

- Prefer `AnnotationConfig` and `run_annotation` rather than inventing a separate annotation entrypoint.
- Support marker-based, reference-based, scoring-based, and hybrid annotation.
- Resolve built-in markers through `get_marker_manager()` views instead of
  loading marker files directly inside annotation code.
- Run lineage/subtype/state/program evidence as separate branches so the final
  annotation can report where each label came from.
- Use CellTypist, marker-manager evidence, data-driven cluster summaries, and
  functional program scores as complementary evidence streams rather than
  mutually exclusive annotation modes.
- Make cluster-level and cell-level outputs distinguishable.
- Preserve uncertainty. Do not force ambiguous clusters into over-specific labels.
- Detect marker conflicts and lineage mixing.
- Explicitly surface artifact evidence such as ribosomal-high, mitochondrial-high,
  hemoglobin-high, ambient-RNA-like, or dissociation-stress signatures before
  accepting an ambiguous lineage call.
- Store annotation evidence in `adata.uns["sclucid"]["analysis"]` or a dedicated annotation namespace.
- For tumor immune annotation, design state vocabularies that can transfer across cancer types.

Current analysis polishing target:

- clustering resolution evidence should record candidate resolutions, selection rationale, marker support, and risks;
- marker-based annotation evidence should produce reviewable lineage, subtype, state, confidence, supporting evidence, and conflict evidence;
- annotation review tables should help users inspect uncertain or mixed clusters before accepting labels;
- malignancy and CNV-assisted interpretation should explain how CNV, epithelial evidence, tumor programs, and reference support affect malignant calls;
- `adata.uns["sclucid"]["analysis"]["review_summary"]` should record key reasoning, not only completed steps;
- `Step2-Annotation_and_Malignancy.ipynb` should become the second product showcase after QC/preprocess.

Lineage-specific minimum rules:

- T cells: distinguish CD4, CD8, Treg, Tfh/CXCL13-like, Trm, cytotoxic, proliferating, pre-exhausted/plastic exhausted, terminal exhausted, stress-response, and tumor-reactive candidate states. If TCR data is present, clonality is evidence, not the sole label source.
- Myeloid cells: distinguish monocyte, macrophage/TAM, cDC1, cDC2, LAMP3+ mature/migratory DC, pDC when available, mast, neutrophil/TAN, and proliferating/activated states. Cancer-type-specific marker shifts should lower or qualify confidence rather than force rejection.
- B cells: distinguish naive, memory, germinal-center-like, plasma/ASC, stress-response memory, atypical/tumor-associated B states, and activation/clonal-expansion evidence when BCR or clonotype metadata exists.
- NK cells: distinguish CD56-bright-like, CD56-dim/cytotoxic, adaptive-like, tissue-resident/tissue-infiltrating, tumor-associated impaired-cytotoxicity, stress-response, and checkpoint-high states. Always check conflicts with NKT and cytotoxic CD8 T cells.
- Stromal/CAF: support myCAF/iCAF/apCAF as baseline but allow spatial CAF subtypes and neighborhood-derived states when spatial data is available.
- Neutrophils: avoid collapsing all neutrophils into one label in tumor liver/TIME-style analyses; support TAN-like subtypes and immunosuppressive/cytotoxicity-suppressing evidence when signatures exist.

## 10. Tumor Module Development Rules

Tumor-specific analysis should solve problems generic scRNA packages do not handle well:

- malignant versus normal epithelial identification;
- CNV-aware malignancy classification;
- tumor program scoring;
- tumor-normal and regional heterogeneity;
- TME composition and ecotypes;
- CAF, myeloid, T cell, NK, B/plasma, and DC state refinement;
- therapy-aware states and response/resistance associations;
- spatial neighborhoods and tumor ecosystem context.

Development rules:

- Tumor modules should consume QC/preprocess/analysis outputs rather than duplicating them.
- Analysis may include lightweight tumor-aware interpretation bridges only when
  they consume tumor-module algorithms or precomputed tumor evidence and expose a
  review contract for downstream tumor workflows.
- Malignancy calls must expose evidence: CNV signal, epithelial markers, tumor program score, reference normal cells, uncertainty.
- Batch correction must be conservative in tumor contexts. Do not erase real inter-patient, regional, tumor-normal, or therapy-related variation.
- Therapy or clinical association functions must clearly distinguish exploratory association from causal inference.
- Use `TumorAnalysisConfig` and `TumorWorkflowConfig` for workflow integration.
- Store tumor outputs under `adata.uns["sclucid"]["tumor"]` and key cell-level labels in `adata.obs`.

Tumor ecosystem rules:

- `TIME_subtype` or `tumor_ecotype` features should combine cell composition, state/program scores, spatial or tissue compartment evidence, and clinical metadata when available.
- Do not call an ecotype from proportions alone unless the review summary labels it as exploratory and lists missing evidence.
- Progression workflows should treat stage/timepoint metadata as first-class context and summarize how immune and stromal states change across stages.
- Spatial CAF/neighborhood workflows should store the neighborhood radius or graph rule, neighbor composition, CAF subtype evidence, nearby immune/tumor/endothelial context, and visualization keys.
- Interaction workflows should prioritize interpretable and clinically relevant interactions, not only high ligand-receptor scores.

## 11. Visualization Module Development Rules

Visualization should be publication-ready by default and scalable for exploratory work.

Every plotting function should:

- accept `AnnData` and explicit key names;
- return a figure or axes object;
- support `save_path` or `save_dir`;
- support high DPI export and editable PDF/SVG text;
- use `src/scLucid/plotting/theme.py` where appropriate;
- support rasterization or downsampling for large datasets;
- avoid overcrowded legends by default;
- support consistent cell type colors;
- expose enough parameters for expert adjustment without making simple usage hard.

Priority plot families:

- QC distributions and threshold evidence plots.
- Preprocess diagnostics: HVG stability, PCA variance, integration/batch mixing.
- Analysis diagnostics: resolution comparison, marker dotplots, annotation confidence.
- Tumor plots: malignancy evidence, CNV summary, tumor programs, TME composition, spatial neighborhoods.
- Report-ready figures: UMAP, dotplot, heatmap, composition bar/box, volcano, alluvial, interaction summary.

Do not make plotting functions only save files. They must return objects for notebooks, reports, and tests.

## 12. External Tool Wrapper Rules

External wrappers must make third-party tools easier, not hide their assumptions.

Wrapper requirements:

- Check optional dependencies at runtime and provide actionable installation messages.
- Convert scLucid/AnnData inputs into the external tool format.
- Run the tool or prepare runnable inputs.
- Convert outputs back into AnnData-compatible results.
- Store config, command/version if available, input keys, and output keys.
- Keep raw external output accessible when useful.
- Provide a small result summary and warnings.

For tools like inferCNV, pySCENIC, CellPhoneDB/CellChat-like analysis, Monocle3-like trajectory, BayesPrism/DWLS-like deconvolution, and spatial/morphology tools, the wrapper should expose biological interpretation helpers rather than only method-native tables.

## 13. Report and Workflow Rules

Workflows should support real exploratory analysis, including reruns, partial failures, and parameter revision.

Workflow rules:

- Use named steps and allow `steps` / `skip_steps`.
- Support progress reporting through existing utilities.
- Store `steps_executed`.
- Use `WorkflowError`, `PartialResultManager`, and `WorkflowCheckpoint` when error recovery is useful.
- Generate a review summary for each major stage.
- Include warnings, user overrides, assumptions, and downstream recommendations.
- Keep outputs reviewable by humans, not just machine-readable.
- When validation scaffold outputs are produced, register their JSON and CSV paths in the workflow manifest.

Review summaries should answer:

- What data was analyzed?
- Which steps ran?
- Which parameters were recommended?
- Which parameters were applied?
- Which choices were user overrides?
- What evidence supports decisions?
- What risks or caveats remain?
- What should the user inspect next?

## 14. Testing and Validation Rules

Tests should scale with risk.

Minimum expectations for new behavior:

- Unit tests for config validation and core calculations.
- Small AnnData tests for input/output location and shape.
- Workflow tests when the feature changes stage behavior.
- Regression tests for previous bugs or risky edge cases.
- Optional dependency tests should skip cleanly when dependencies are missing.

Validation rules:

- Use `validate_adata` and `validate_stage_contract` where applicable.
- Check required layers, obs columns, obsm keys, and uns namespaces.
- Validate non-empty AnnData.
- Validate shape consistency after filtering/subsetting.
- Validate that review summaries pass `validate_review_summary_schema`.

Validation tiers:

1. Lightweight gates: unit tests, contract checks, review summary schema checks, and small AnnData behavior.
2. Maturity scaffold: compact QC/preprocess validation tables proving auditability, reproducibility, and workflow readiness.
3. Golden path acceptance: PBMC and tumor notebooks/scripts that run end to end with stable outputs.
4. Comparative validation: fixed datasets and metrics comparing scLucid with Scanpy standard workflow and optional Seurat/scran-style references.
5. Manuscript-level benchmark: broader datasets, statistical summaries, runtime comparisons, and figure-quality reporting.

Do not conflate these tiers. Current lightweight QC/preprocess validation means "ready for comparative validation", not "scientifically proven better than standard workflows".

Current useful datasets:

- PBMC3K for normal baseline behavior.
- Lin2020 PDAC for tumor-aware QC, preprocessing, annotation, and malignancy/CNV workflow fit.
- Schlesinger2020 PDAC as a second tumor dataset candidate for generalization checks.

## 14A. Current Roadmap Snapshot

Use this snapshot when resuming work from memory:

**Phase 1: Analysis Benchmark Module**

- Default analysis route should mature toward:
  `clustering_review -> clustering -> markers -> annotation_evidence ->
  annotation_consensus -> optional_malignancy_interpretation ->
  analysis_review_summary`.
- `scripts/run_analysis_acceptance.py` is the scriptable Step2 acceptance runner
  and should stay synchronized with
  `examples/03_advanced_notebooks/Step2-Annotation_and_Malignancy.ipynb`.
- First-pass global annotation should stay conservative at lineage / major cell
  type level.
- Subtype and state labels should usually be generated after subset extraction
  and reclustering, or by explicit advanced options.
- LLM output is evidence and rationale, not final truth.
- All marker-dependent annotation evidence should use `Manager` resources and
  `get_marker_manager()` views by default.

**Phase 2: Tumor Interpretation**

- Treat `analysis.run_malignancy_interpretation` as a bridge/contract layer,
  not a second tumor module.
- Use tumor-module functions for CNV inference and malignancy scoring when the
  bridge needs computed evidence.
- Use CNV, epithelial/tumor context, tumor programs, reference anchors,
  cancer-type hints, and manual evidence.
- Do not let tumor markers contaminate global annotation.
- Store analysis-stage bridge summaries under
  `adata.uns["sclucid"]["analysis"]["malignancy"]`; store full tumor workflow
  outputs under `adata.uns["sclucid"]["tumor"]`.

**Phase 3: Resource Curation**

- Continue curation of marker metadata and negative markers.
- Add mouse tissue/tumor parity once human resources stabilize.
- Keep broad signatures in JSON/GMT gene-set resources.

**Phase 4: Real-Data Validation**

- PBMC baseline.
- PDAC tumor workflow.
- Second tumor dataset.
- Active project notebook acceptance.
- Standard workflow comparison only after analysis and tumor contracts are stable.

**Phase 5: Output Polish**

- Advanced notebooks should become publication-quality workflow narratives.
- Audit reports should include analysis and tumor interpretation maturity.
- Plotting should support reusable top-journal figure templates.

## 15. Documentation Rules

Documentation should explain both API usage and decision logic.

For public APIs, include:

- concise docstring;
- key parameters;
- return value;
- where results are stored in AnnData;
- small usage example;
- notes about tumor-aware behavior or assumptions.

For workflow changes, update relevant docs or examples when the user-facing API changes.

When changing QC, preprocess, validation, or advanced notebook contracts, check whether these docs/examples need synchronized updates:

- `README.md`
- `docs/source/workflow_hardening.rst`
- `docs/source/validation_scaffold.rst`
- `docs/source/data_contracts.rst`
- `docs/source/qc_preprocess_maturity.rst`
- `docs/source/usage_layers.rst`
- `docs/source/notebooks.rst`
- `docs/source/examples.rst`
- `examples/03_advanced_notebooks/`
- golden path scripts and their manifests

For recommendation or automated behavior, documentation must explain:

- what evidence is used;
- when automation is unreliable;
- how users can override decisions;
- where the decision trace is stored.

## 16. Codex Task Output Format

When implementing a scLucid feature, Codex should report:

1. Files changed.
2. API added or changed.
3. AnnData input/output contract.
4. Config changes.
5. Storage keys written under `adata.uns["sclucid"]`.
6. Recommendation or trace behavior.
7. Validation or benchmark claim boundary.
8. Docs/examples updated.
9. Tests run and results.
10. Known caveats or next best follow-up.

For code review tasks, prioritize findings first with file and line references.

For exploratory design tasks, produce a scoped plan, identify the correct existing modules, and avoid proposing a rewrite unless the current architecture truly blocks the goal.

## 17. Common Mistakes to Avoid

Avoid these patterns:

- Creating a new workflow path that duplicates `run_standard_qc`, `run_preprocessing`, or `run_standard_analysis`.
- Storing results in arbitrary top-level `adata.uns` keys.
- Mutating AnnData unexpectedly without `inplace=True`.
- Using fixed tutorial-style thresholds when data-driven recommendations are available.
- Treating tumor tissue like normal PBMC data.
- Over-correcting tumor data and erasing biological heterogeneity.
- Calling clusters “cell types” without marker, reference, or state evidence.
- Dropping rare populations without checking cross-cancer or state-level evidence.
- Returning only saved files from plotting functions.
- Hiding optional dependency failures.
- Adding broad abstractions before two or more modules need them.
- Adding verbose comments that restate the code.
- Forgetting review summaries and decision traces.
- Claiming "better than standard workflow" without comparative validation.
- Creating validation metrics that cannot be reproduced from AnnData contracts, manifests, or saved artifacts.
- Updating only one usage layer while leaving workflow/API/docs/examples inconsistent.

## 18. Example Development Prompts

Use prompts like these when continuing development:

```text
Improve the QC review summary so it records recommended thresholds, applied thresholds, user overrides, tumor-aware warnings, and downstream preprocess recommendations.
```

```text
Add a tumor-aware HVG recommendation rule that preserves cancer, immune-state, and stress-response signatures while still using the existing preprocess/hvg architecture.
```

```text
Extend analysis annotation so cluster labels are converted into lineage, subtype, state, confidence, supporting markers, and conflicting markers using the existing AnnotationConfig.
```

```text
Create a publication-ready annotation confidence plot that accepts AnnData, returns a matplotlib figure, uses the scLucid theme, and stores no unexpected state.
```

```text
Add an integration diagnostic to the preprocess workflow that warns when batch correction may remove tumor-normal or inter-patient heterogeneity.
```

```text
Implement a review-summary validator for analysis outputs that checks clustering, marker, annotation, and characterization sections.
```

```text
Extend the lightweight validation scaffold so analysis can report clustering evidence, annotation confidence, malignancy/CNV evidence, and review-summary completeness after the analysis module reaches QC/preprocess maturity.
```

```text
Polish Step2-Annotation_and_Malignancy.ipynb so it explains why clustering, annotation, and malignant calls were accepted or flagged, using existing analysis APIs and stored review summaries.
```
