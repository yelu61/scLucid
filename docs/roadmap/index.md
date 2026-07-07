# Roadmap

The roadmap is the phase-level execution playbook for building scLucid toward
submission-grade evidence and reproducible release assets. It records where the
project is going; the current user-facing module contract lives in
[Module Features And Stage Plan](../user/module_features_and_plan.md).

## Phases

- [Phase 1: Core API And Contracts](PHASE_1_CORE_API_AND_CONTRACTS.md)
- [Phase 2: QC Evidence Benchmark](PHASE_2_QC_EVIDENCE_BENCHMARK.md)
- [Phase 3: Preprocess Analysis Validation](PHASE_3_PREPROCESS_ANALYSIS_VALIDATION.md)
- [Phase 4: Tumor Interpretation Case Studies](PHASE_4_TUMOR_INTERPRETATION_CASE_STUDIES.md)
- [Phase 5: Tumor Ecosystem Modeling](PHASE_5_TUMOR_ECOSYSTEM_MODELING.md)
- [Phase 6: Knowledge Evidence Infrastructure](PHASE_6_KNOWLEDGE_EVIDENCE_INFRASTRUCTURE.md)
- [Phase 7: Support Evidence And R/Python Parity](PHASE_7_SUPPORT_EVIDENCE_AND_R_PYTHON_PARITY.md)
- [Phase 8: Release Manuscript And Submission](PHASE_8_RELEASE_MANUSCRIPT_AND_SUBMISSION.md)

## Vision Tasks Across Phases

These tasks translate scLucid's differentiation vision into the existing phase
plan. They should be implemented inside the phase files above rather than as a
separate roadmap.

| Vision task | Primary phase | Supporting phases | Expected evidence |
|----|----|----|----|
| Context-aware decisions | Phase 1 | Phase 2, Phase 3, Phase 4 | Review summaries record dataset context, analysis goal, recommended choice, applied choice, rationale, risk, and limitation. |
| QC as biological-risk attribution | Phase 2 | Phase 3, Phase 4 | Benchmarks compare QC policies by retention, marker/program fidelity, tumor purity, immune/TME composition, and review-required cells. |
| Unified claim and evidence semantics | Phase 1 | Phase 3, Phase 4, Phase 6, Phase 7 | A consistent vocabulary for `claim_level`, `inference_level`, `evidence_level`, confidence, source, and limitation appears in module outputs and docs. |
| Tumor interpretation case studies | Phase 4 | Phase 2, Phase 3, Phase 5 | PDAC and at least one second tumor route show malignancy, CNV, TME, program, and therapy evidence with explicit limits. |
| Ecosystem-level interpretation | Phase 5 | Phase 4, Phase 6, Phase 7 | Sample-level ecosystem feature matrices, archetype/ecotype prototypes, stability reports, and support evidence summaries. |
| Knowledge and evidence routing | Phase 6 | Phase 4, Phase 5, Phase 7 | Marker, gene set, atlas, literature, ontology, therapy, and LLM-suggestion evidence use a shared schema and provenance record. |
| Productized review surfaces | Phase 8 | All phases | Audit reports, Methods-ready summaries, figure source data, and release assets expose the evidence contract without overclaiming implementation maturity. |
