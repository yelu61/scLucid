# Full Preprocess Validation Protocol

## Release question

The Preprocess release question is not whether one mixology run can be
completed. It is whether scLucid can preserve count semantics, select a
consumer-appropriate transformation and feature space, build a stable
neighborhood representation, decline unsafe integration, preserve rare and
tumor biology, and execute the reviewed policy reproducibly in real projects.

The post-QC authoritative count matrix is the scope boundary. FASTQ processing,
cell calling, ambient correction, and cell removal belong to the QC validation
contract. Formal downstream inference belongs to Analysis. A Preprocess `PASS`
therefore supports only the construction and selection of declared expression
or neighborhood representations; it does not validate a biological mechanism
or a condition effect.

The machine-readable contract is split between:

- `validation/dataset_evidence_registry.json`, which owns endpoint estimands,
  metrics, thresholds, dataset accessions, truth types, and metadata;
- `validation/qc_preprocess/acceptance_contract.json`, which owns the exact
  endpoint-by-dataset release portfolio and cross-head rules;
- `validation/evidence_run_index.json`, which binds executed artifacts to one
  dataset and one endpoint.

No aggregate Preprocess score is permitted. Every required head must pass its
own exact dataset portfolio.

## Twelve independent evidence heads

| Head | Primary validation target | Locked acceptance boundary | Exact required datasets |
|---|---|---|---|
| `pp_input_representation_contract` | Preserve authoritative counts and keep `counts`, `normalized_full`, `discovery_rep`, and `integrated_rep` semantically distinct | exact count preservation and representation provenance 100%; zero consumer-contract violations; review mutation count 0 | controlled Preprocess contract fixture; PBMC3k engineering fixture; scMixology; three real projects |
| `pp_normalization_selection` | Select shifted-log, Pearson residual, or another eligible normalization for an unseen protocol or dataset | held-out utility regret <=5%; biology loss <=2%; shifted-log baseline always evaluated; zero silent fallback; every compared candidate assumption-eligible | scMixology; Mereu; Ding; three real projects |
| `pp_feature_selection` | Select unsupervised features without circularly forcing protected markers into the discovery graph | utility regret <=5%; rare-class recall and program-preservation loss each <=2%; all three preregistered feature-count variants pass; zero protected-marker injection | scMixology; Mereu; Kinker; three real projects |
| `pp_selector_regret` | Generalize the full candidate selector beyond one protocol | held-out regret <=5%; biology loss <=2%; selected candidate consistent across all three preregistered parameter variants; all eligible candidates evaluated | scMixology; Mereu; Ding; three real projects |
| `pp_graph_stability` | Preserve neighborhood conclusions across seeds, feature counts, PCs, and neighborhood sizes | graph-seed and neighbor-identity loss each <=2% relative to the simple baseline; partition-stability regret <=5%; conclusion consistency 100% | scMixology; scIB pancreas; three real projects |
| `pp_integration_need_confounding` | Distinguish no-integration, review, and fail-closed cases | unintegrated baseline coverage 100%; false `READY` under registered confounding 0%; unnecessary-integration recommendation <=5%; condition-structure loss <=2%; Cramer's V >=0.7 blocks automatic integration | scDesign3 controlled scenarios; scMixology; scIB pancreas; three real projects |
| `pp_integration_pareto` | Improve batch behavior without trading away biology, rare populations, programs, or graph stability | unintegrated baseline mandatory; complex method must Pareto-improve; biology, rare-population, program, and graph-stability loss each <=2% | scMixology; Mereu; Ding; scIB pancreas; three real projects |
| `pp_identity_preservation` | Preserve known identity and low-frequency populations across protocol or batch shifts | utility regret <=5%; rare-class recall loss <=2%; absolute rare-population abundance bias <=2% | scMixology; Mereu; scIB pancreas; 10x HGMM; three real projects |
| `pp_tumor_structure_preservation` | Preserve patient, lineage, and tumor-program structure relative to the simple unintegrated baseline | lineage-retention, tumor-program-correlation, and patient-structure loss each <=2% | Lin PDAC; Moncada PDAC; Zilionis NSCLC; Lee CRC; Kinker pan-cancer lines; three real projects |
| `pp_policy_execution` | Apply only the explicitly reviewed policy and reproduce its outputs | review mutation count 0; policy/apply agreement, repeat-run agreement, count preservation, and representation provenance each 100% | controlled Preprocess contract fixture; PBMC3k engineering fixture; scMixology; three real projects |
| `pp_decisioncard_ux` | Make the correct consumer-specific action usable without notebook workarounds | user-edited fields reduced >=70%; zero critical-action errors, manual workarounds, or project-specific patches; RunEvidence completion 100% | three real projects |
| `pp_scalability` | Retain sparse, bounded, reproducible execution on declared hardware | zero unintended dense expansion and failures; repeated policy agreement 100%; runtime and peak memory reported; no hardware-independent wall-time claim | controlled Preprocess contract fixture; PBMC3k engineering fixture; 10x HGMM; Kinker; three real projects |

The thresholds are guardrails for the registered tasks, not universal constants
for all tissues or consumers. Any post-hoc threshold change requires a contract
version bump, a written reason, and rerunning every affected binding.

## Representation and consumer contract

The four spaces are not interchangeable:

| Representation | Required consumer boundary |
|---|---|
| `counts` | authoritative untransformed counts for pseudobulk and count models |
| `normalized_full` | full-gene normalized expression for markers, programs, and interpretable plots |
| `discovery_rep` | unsupervised feature space for PCA, neighbors, and clustering |
| `integrated_rep` | optional graph or latent representation for an explicitly reviewed integration consumer |

An integrated graph or latent space is not corrected expression. It cannot be
silently used for marker testing or pseudobulk differential expression.
Pearson residuals require UMI-count assumptions. scran is eligible only when
its declared dependency is available. No method may silently fall back to a
different implementation.

Protected markers are reviewed in `normalized_full`. They are not forced into
the default unsupervised feature union. A hypothesis-driven feature set is a
labelled sensitivity analysis and cannot replace the unsupervised result.

## Evaluation order and leakage control

Validation follows the dependency chain:

1. lock source hashes, metadata, counts, software versions, seeds, and the
   intended downstream consumer;
2. verify the representation contract before comparing scientific utility;
3. fit normalization, feature selection, PCA, and selector decisions only on
   training protocols, batches, or datasets;
4. evaluate the held-out protocol or dataset once using identity, rare-class,
   program, depth-dependence, and stability endpoints;
5. evaluate the unintegrated baseline before any integration method;
6. block automatic integration when batch and condition are not identifiable;
7. test reasonable seeds, feature counts, PC counts, and neighborhood sizes;
8. execute the selected policy through the public apply path and record
   `RunEvidence`;
9. repeat the maintained four-action workflow in 202604JJH, 202507LPJ, and
   202603AK112 without project-specific code patches.

Cells are not independent validation replicates. The experimental unit is the
held-out protocol, dataset, library, patient, or project declared by the
endpoint. Confidence intervals and sensitivity summaries must respect that
unit.

## `NOT_EVALUABLE` is not a pass

`NOT_EVALUABLE` is allowed only when a head or candidate lacks a required input
or declared dependency. The evidence record must name:

- the affected head;
- the missing requirement;
- the affected candidates;
- the resulting claim limitation;
- the single next action that would make the question evaluable.

Examples include a count-based normalization candidate without authoritative
counts, or scran without its declared dependency. A confounded batch/condition
design is different: the biological correction may be non-identifiable, but
the validation scenario is evaluable and should pass only if scLucid returns
the required fail-closed decision. `NOT_EVALUABLE`, simulation-only, and
contract-only statuses never count as release `PASS`.

## Dataset roles and truth limits

- scMixology, Mereu, Ding, and 10x HGMM provide controlled identities and
  protocol or species perturbations. They are strong for technical
  generalization but weak models of a primary tumor ecosystem.
- scIB pancreas provides a multi-study integration stress test with
  silver-standard annotations; it is not experimental truth.
- scDesign3 provides registered positive and negative confounding controls. A
  simulation can validate mechanism sanity and fail-closed behavior, not
  external superiority.
- Kinker provides controlled cell-line identity and recurrent program
  structure, but cell lines do not reproduce fragile TME populations.
- Lin, Moncada, Zilionis, Lee, and the three real projects test patient,
  lineage, and tumor-program preservation. Their expert or author labels are
  review or silver standards, not complete cell-level ground truth.
- The deterministic controlled Preprocess fixture is the provenance source for
  the maintained 120-cell by 160-gene public review/apply contract benchmark.
  Its `CONTRACT_PASS_NOT_PERFORMANCE` result is engineering evidence, not
  scientific or production-scale evidence.
- PBMC3k is a separate engineering fixture. It may pass representation,
  determinism, and memory contracts but cannot support scientific superiority.

## Current conclusion

The deterministic controlled contract benchmark passes all 15 registered
public review/apply checks and is recorded as
`CONTRACT_PASS_NOT_PERFORMANCE`. This establishes the local API and
representation-contract mechanism on its 120-cell by 160-gene synthetic
fixture. It does not establish PBMC3k execution, production scale, scientific
utility, or real-project usability, and therefore cannot pass those endpoint
portfolios by itself.

The controlled scMixology run supports two local bindings. The selected simple
unintegrated policy has regret about 0.49% on the registered held-out protocol
task and is stable across the three registered HVG/PC variants, so
`scmixology_gse118767 × pp_selector_regret` is `PASS`. The same run retains the
unintegrated candidate because Harmony did not satisfy every registered loss
guardrail, so `scmixology_gse118767 × pp_integration_pareto` is
`PASS_BASELINE`. Neither result passes the corresponding external dataset
portfolio.

It does **not** pass normalization selection on external technologies, feature
selection for tumor programs, graph stability on primary tissue, integration
need/confounding, the full multi-dataset Pareto integration head, tumor
preservation, policy execution in real projects, UX, or scalability. Until all
12 head portfolios pass,
Preprocess remains `REVIEW`; the allowed claim is a controlled mixology result,
not general readiness or superiority over traditional preprocessing.
