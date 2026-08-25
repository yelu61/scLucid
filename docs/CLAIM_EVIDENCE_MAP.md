# scLucid Claim–Evidence Map

**Status**: working draft, 2026-08-24. Owner: maintainer.
**Purpose**: define what scLucid claims, at what strength, and what evidence
must back each claim — so that "done" is finite and checkable, the paper
skeleton is derivable, and development effort is not spread over surfaces that
carry no claim.

This document is the bridge between three things that were previously loose:

- the **positioning** ("tumor data has systematic blind spots in standard
  pipelines; analyze like an expert"),
- the **validation machinery** (`validation/qc_preprocess/acceptance_contract.json`,
  `validation/evidence_run_index.json`, `validation/dataset_evidence_registry.json`),
- the **paper** (each locked endpoint maps to a figure panel).

It does not replace the strategic plan or the full validation protocols; it
assigns every module a claim tier and every tier a definition of done.

---

## Part 0 — The three claim tiers

Every module, feature, and catalog scenario belongs to exactly one tier.
The tier defines what "完善" means for it. Nothing else is required.

| Tier | Name | Claim made | Definition of done | Evidence required |
|------|------|-----------|--------------------|-------------------|
| **T1** | Core decision points | "Naive defaults systematically fail here on tumor data; scLucid's risk-attributed decision is measurably better or safer." | Locked endpoint passes on the full preregistered dataset portfolio, incl. grouped-bootstrap CI | Truth-bearing or reference datasets, pre-registered numeric thresholds, ≥1 real-project instance |
| **T2** | Supporting analyses | "Correct, safe, auditable." No superiority claim. | Inference contracts enforced (replicate-aware, fail-closed), engineering tests green | Contract tests, simulation-level sanity; explicitly NOT external superiority benchmarks |
| **T3** | Convenience / compatibility surface | Nothing. It works and is labeled. | Smoke tests, API stability | None beyond tests |

**Rules that keep the scope finite:**

1. **No T1 claim without a locked endpoint.** If a module has no endpoint in
   `acceptance_contract.json`, it is not T1, and no benchmark is owed for it.
2. **T2 is explicitly finished without benchmarks.** DE, proportion, clustering,
   annotation are done when inference safety is enforced. Resist adding
   accuracy competitions here.
3. **T3 must never grow a claim.** Plotting, `tools/` R-ports, IO helpers are
   compatibility surfaces. They are documented as support evidence only.
4. **`NOT_EVALUABLE` / `CONTRACT_PASS_NOT_PERFORMANCE` /
   `SIMULATION_PASS_NOT_EXTERNAL` never count as T1 pass** (per the full
   validation protocols).

### Module × tier assignment (current)

| Module | Tier | Claim (if T1) | Locked endpoints | Current status |
|--------|------|---------------|------------------|----------------|
| QC decision layer (`qc/policy/reviewer.py`) | **T1** | Risk-attributed filtering beats fixed-threshold defaults on tumor data | `qc_selector_superiority`, `qc_damage_classification`, `qc_catastrophic_sample_detection`, `qc_rare_population_preservation`, `qc_doublet_calibration`, `qc_ambient_correction`, `qc_cell_calling`, `qc_profile_selection` | 2 FAIL, 3 BLOCKED, 1 REVIEW, engineering suite PASS only |
| Preprocess decision layer (`preprocess/policy.py`) | **T1** | Consumer-specific representation choice with bounded regret; integration only when needed and Pareto-improving | `pp_selector_regret`, `pp_integration_need_confounding`, `pp_integration_pareto`, `pp_tumor_structure_preservation`, `pp_graph_stability`, `pp_identity_preservation`, `pp_normalization_selection`, `pp_feature_selection` | Mixology selector PASS, integration PASS_BASELINE; 10 of 12 heads not passed |
| Analysis inference contracts | **T2** | — | (contract only) `an_pseudobulk_de`, `an_composition`, `an_annotation_accuracy` | Fail-closed slice validated 2026-08-17; external benchmarks not run |
| Annotation (multi-evidence consensus) | **T2** | — | `an_annotation_accuracy` (regret ≤5%, optional) | Contract path exists; accuracy benchmark not run |
| Tumor module | **T2→frozen** | — (application claims deferred to the application chapter, Part C) | none locked | Feature development frozen until P0 passes |
| `tools/` R-ports (pyMonocle3, pyCellChat, pyBayesPrism, pyDWLS, …) | **T3** | — | none | Compatibility/support evidence only; candidates for de-emphasis |
| Plotting, IO, utils | **T3** | — | none | Tests only |

> TODO(owner): confirm the tumor module's tier. Recommendation: keep tumor
> interpretation T2 (correctness) in the tool paper, and let biological
> novelty live exclusively in the application chapter (Part C, Fig 5).

---

## Part 1 — Tumor decision-scenario catalog

One row = one **decision scenario where a naive default can go wrong on tumor
data**. This catalog is the paper's biological motivation, the differentiated
feature list, and the concrete meaning of "analyze like an expert".

A scenario is **established (立得住)** when all three legs hold:

- **Lit** — at least one literature anchor that the failure mode is real
- **Pub** — demonstrated on at least one public dataset with truth or reference
- **Own** — observed in at least one own real project (202604JJH / 202507LPJ /
  202603AK112 or another clinical/mouse-model dataset)

Two legs = candidate, stated as hypothesis. One leg = do not claim.

**Ecosystem lens (motivation, not a new claim).** Every scenario below is a
cell-level decision, but its harm is ultimately measured at the ecosystem
level: TME composition estimates, malignant program scores, and
patient/ecotype structure are the downstream victims of each naive default.
Standard-pipeline blind spots rarely fail at cell counts — they fail at
ecosystem conclusions. The "Ecosystem stakes" column records this for each
row. Scope guard: ecosystem readouts are **evaluation rulers** for the
existing T1 claims (candidate secondary endpoints), not new claims; locking
any of them as an endpoint still follows the Part 3 rules.

| # | Scenario | Ecosystem stakes | Naive default failure | Expert judgment encoded | scLucid evidence head / API | Locked endpoint & threshold | Lit | Pub | Own | Status |
|---|----------|------------------|----------------------|-------------------------|-----------------------------|------------------------------|-----|-----|-----|--------|
| S1 | **Stress/high-MT injury**: dissociation stress and tumor epithelial physiology raise MT%; viable stressed cells are not dead | Removing stressed/fragile cells systematically biases TME composition estimates — stressed epithelial/malignant fractions are undercounted | Fixed 5–10% MT filter removes viable stressed tumor/epithelial cells and erases stress programs | Per-sample MT distribution + joint RNA–MT modeling; MT high ≠ dead | `qc_profile_selection`, `qc_damage_classification`; miQC-family candidate (currently sensitivity proxy) | `qc_selector_superiority`: KEEP false removal ≤2%, recall gain ≥5 pt vs both MAD baselines, 95% CI lower bound >0 | van den Brink 2017 (dissociation stress); Hippen 2021 (miQC) | E-MTAB-2600 (damage labels) — currently **FAIL** | TODO(owner): which project showed this | **FAIL — top repair target** |
| S2 | **Catastrophic sample failure**: an entire library is degraded and quietly poisons integration | One failed library can turn ecotype clustering into batch/quality clustering | Sample enters pooling/integration; batch structure is artifact | Sample-level gate before cell-level filtering | `qc_catastrophic_sample_detection` | detection rate 1.0; locked high-quality false-block rate 0.0 | SampleQC 2022 (multivariate sample QC) | lin2020 PDAC (GSM4679533 must be blocked) — currently **BLOCKED** | TODO(owner) | **BLOCKED** |
| S3 | **Doublets vs. real transition/rare states**: EMT-like, cycling, or fragile lineages look doublet-like | Killing transition/rare states erases EMT-like or cycling compartments from the ecosystem composition vector | Fixed-rate doublet removal kills rare/transitional populations | Calibrate per sample; check marker coherence before removal; rare-population guard | `qc_doublet_calibration`, `qc_rare_population_preservation`; scDblFinder parity | Kang demuxlet + cell-hashing (GSE108313) + HGMM; AUPRC/calibration thresholds per full protocol | Kang 2018 (demuxlet); Germain 2022 (scDblFinder) | Kang **REVIEW** (AUPRC 0.605); HGMM **PASS** | TODO(owner) | **REVIEW — calibration gap** |
| S4 | **Ambient RNA contamination**: necrotic tumor tissue has high ambient; marker reads are contaminated (e.g. HBB/Ig signal in epithelium) | Contaminated counts blur tumor/stroma/immune boundaries; marker-based annotation and program scores are polluted | Marker interpretation and DE on contaminated counts | Estimate contamination fraction; correct conservatively or flag; never silent | `qc_ambient_correction` | HGMM mixture truth; contamination reduction ≥50% (full protocol) | Young & Behjati 2020 (SoupX); Fleming 2023 (CellBender) | HGMM **FAIL**; CellBender tiny fixture = plumbing only | TODO(owner) | **FAIL — top repair target** |
| S5 | **Low-RNA biology vs. damaged cells**: neutrophils/platelets have genuinely low RNA | Dropping low-RNA biology (neutrophils/platelets) leaves the TME immune compartment systematically incomplete | `nFeature` floor removes real low-RNA cell types | Distinguish empty / damaged / low-RNA-biology before filtering | `qc_cell_calling`, `qc_damage_classification` | E-MTAB-2600 microscopy labels; true-cell recall ≥95%, empty FDR ≤1% | Ilicic 2016 (microscopy QC labels) | **BLOCKED** (no scLucid call output; no low-RNA truth subset) | TODO(owner) | **BLOCKED** |
| S6 | **Integration necessity**: not every dataset needs integration; overcorrection erases tumor-vs-normal / patient biology | Overcorrection erases real patient-to-patient ecosystem differences; ecotypes reflect batch, not biology | Always integrate (Harmony/scVI by default) | Test batch–biology confounding first; Cramér's V ≥0.7 blocks auto-integration; complex method must Pareto-dominate unintegrated baseline | `pp_integration_need_confounding`, `pp_integration_pareto` | Mixology: PASS_BASELINE (Harmony did not dominate); multi-dataset Pareto portfolio per contract | Luecken 2022 (scIB benchmark) | Mixology **PASS_BASELINE**; multi-dataset portfolio not run | TODO(owner) | **Partial — extend beyond Mixology** |
| S7 | **Preprocessing choice regret**: tutorial defaults copied across datasets | A wrong representation warps the neighborhood graph that all program and composition readouts are built on | Same normalization/HVG/PC settings everywhere | Held-out protocol evaluation; regret-bounded selection | `pp_selector_regret` | Mixology leave-one-protocol-out: regret ≤5%, biology loss ≤2% — **PASS** (regret ≈0.49%) | Tian 2019 (scMixology); Mereu 2020; Ding 2020 | **PASS** (Mixology); Mereu/Ding/scIB not acquired | — | **PASS on one dataset — portfolio incomplete** |
| S8 | **Tumor structure preservation**: preprocessing/integration must not distort malignant programs or TME structure | The ecosystem feature matrix is computed directly on the structure preprocessing may have destroyed | Method ranked on batch-mixing metrics only | Identity/program/structure retention audited alongside mixing | `pp_tumor_structure_preservation`, `pp_identity_preservation` | lin2020/moncada2020 PDAC, zilionis2019 NSCLC, lee2020 CRC; biology/program/graph losses ≤2% | TODO(owner): cite tumor-atlas integration caution literature | Not run | TODO(owner) | **Not run** |
| S9 | **Graph/clustering instability**: neighborhoods and clusters unstable under resampling or parameter jitter | Unstable graphs make ecotype assignments irreproducible | One clustering presented as the structure | Stability quantified; unstable regions flagged for review | `pp_graph_stability` | Per full protocol thresholds; primary tissue datasets | TODO(owner) | Not run | — | **Not run** |
| S10 | **Cell-level exploration ≠ sample-level conclusion** (pseudo-replication) | Cell-level conclusions make sample-level ecosystem–clinical associations statistically invalid | Cell-level DE p-values read as biological evidence | Replicate-aware pseudobulk; fail closed without replicates; inference level labeled | `an_pseudobulk_de`, analysis inference contracts | Squair compendium empirical FDR ≤0.06, sign concordance ≥95%; muscat simulations | Squair 2021 (pseudo-replication); Crowell 2020 (muscat) | **CONTRACT_PASS_NOT_PERFORMANCE** (contract only); external benchmarks not acquired | ✓ (design intent from own projects) | **Contract pass — performance unverified** |

### Gaps this catalog exposes (decision-relevant)

1. **Own-project leg is unwritten everywhere.** The three real projects exist
   in the contract but no per-scenario incident log exists. This is the
   cheapest high-value gap to close: for each S1–S10, record one concrete
   incident (or "not observed") from JJH/LPJ/AK112 or mouse-model data.
2. **FAIL endpoints concentrate in QC (S1, S4, S5).** These are exactly the
   scenarios with the strongest literature anchors — fixing them is both the
   gate blocker and the paper's core evidence.
3. **S8/S9 have no literature anchor or dataset run yet** — decide whether
   they stay in the catalog as hypotheses or are scoped out of v1 of the map.

---

## Part 2 — Claim → endpoint → figure map (paper skeleton)

| Paper element | Claim | Backing endpoints | Dataset portfolio | Status |
|---------------|-------|-------------------|-------------------|--------|
| **Fig 1** Framework | Tumor scRNA-seq analysis is a chain of auditable biological-risk decisions | (design figure; no endpoint) | — | Draftable now |
| **Fig 2** QC blind spots | S1–S5 failure modes are real and scLucid's decisions are safer | `qc_selector_superiority`, `qc_damage_classification`, `qc_ambient_correction`, `qc_cell_calling`, `qc_catastrophic_sample_detection`, `qc_doublet_calibration` | E-MTAB-2600, HGMM, Kang, cell-hashing, SampleQC metrics, lin2020, real_project_panel | 2 FAIL / 2 BLOCKED / 1 REVIEW — **the gate** |
| **Fig 3** Preprocess decisions | Bounded-regret representation choice; integration only when justified | `pp_selector_regret`, `pp_integration_need_confounding`, `pp_integration_pareto`, `pp_tumor_structure_preservation`, `pp_graph_stability` | Mixology, Mereu, Ding, scIB, PDAC×2, NSCLC, CRC | 1 PASS / 1 PASS_BASELINE / rest not run |
| **Fig 4** Real-project interception | REVIEW/BLOCKED decisions caught concrete errors in real tumor projects | UX acceptance: config-field reduction ≥70%, zero project-specific patches, run evidence complete | 202604JJH, 202507LPJ, 202603AK112 | **BLOCKED — records incomplete** |
| **Fig 5** (optional, tier-up) Application | Novel biological finding in one cancer type or pan-cancer public cohort | separate workstream; reuses A outputs | own clinical cohort or public pan-cancer | Not started — must not block Fig 2–4 |

**Publication decoupling rule**: Figures 2–4 constitute the tool/method paper.
Figure 5 is the tier-raiser and a separate workstream. Do not let Fig 5
thinking delay the locked gate.

**Ecosystem readouts as secondary endpoints (optional, not locked).** For
Fig 2–3, ecosystem-distortion metrics may strengthen existing claims as
downstream harm measures — e.g. TME composition bias before/after S1-type
removal, or patient ecosystem-difference retention before/after S6-type
integration. These sharpen the ruler; they do not add claims, and none become
locked endpoints without following Part 3.

---

## Part 3 — What this map forbids (scope guards)

1. No new T1 claims without a locked endpoint and portfolio entry.
2. No accuracy/superiority benchmarks for T2 modules (DE/proportion/annotation
   beyond the already-locked contract endpoints).
3. No claims at all for T3 (`tools/`, plotting, IO).
4. Simulation results establish engineering only; they are never external
   scientific evidence (already encoded in the contract's status semantics).
5. The catalog grows only by the three-leg rule; a scenario without a path to
   all three legs is a hypothesis, not a feature.

## Part 4 — Tumor research-question map (application layer)

This part maps **current tumor single-cell research trends** to the biological
questions scLucid is *designed to help answer* — mainly around ecosystems and
drug sensitivity/mechanism. Two disciplines apply:

1. **Voice discipline.** Nothing here is claimable today. In the tool paper
   these rows appear only as motivation ("designed to answer"); the verb
   "answers" is reserved for the application paper (Fig 5) after the locked
   gate passes and per-row validation exists.
2. **Status honesty.** Rows marked *skeleton* have real but frozen,
   unvalidated code (T2 at best). Rows marked *gap* have no implementation.
   The strongest trend rows are currently gaps — that asymmetry is deliberate
   information, not an oversight.

| # | Research trend | Biological question | scLucid asset / gap | Status today | Path to claim | Paper |
|---|---------------|---------------------|---------------------|--------------|---------------|-------|
| R1 | TME ecotypes / multicellular communities | How many ecosystem types exist in this cohort, and how do they associate with outcome or treatment response? | `microenvironment/ecosystem.py`, composition analysis; sample-level ecosystem feature matrix and ecotype prototype (Phase 5) **not built** | Skeleton; core piece missing | Build Phase 5 matrix → cross-cohort ecotype stability → outcome association in own/public cohort | Fig 5 (anchor lit: Bagaev 2021) |
| R2 | Malignant cell states / meta-programs | Which malignant programs recur, and how plastic are they? | `find_tumor` / malignancy scoring, `heterogeneity`, program scoring | Skeleton | Program calling validated against CNV/references; cross-cohort recurrence audit | Fig 5 (anchor lit: Gavish 2023) |
| R3 | Clonal evolution under therapy | Which clones/states are selected by treatment, and where does resistance originate? | `evolution`, CNV (inferCNV-style) | Skeleton | Clone calling validated (e.g. Numbat/CopyKAT cross-check); pre/post cohort | Fig 5 |
| R4 | Drug sensitivity / resistance mechanism | Which cell states predict response or resistance (incl. ICB-resistance programs)? | `therapy/prediction.py`, `therapy/resistance.py`; bulk bridge pyBayesPrism/pyDWLS | Skeleton | State→response association replicated across cohorts; bulk↔sc bridge benchmarked | Fig 5 |
| R5 | Cell–cell interaction / druggable axes | Which ligand–receptor axes drive malignant phenotypes and are targetable? | `microenvironment/interaction.py`, pyCellChat | Skeleton | Interaction calls vs. known axes; spatial or perturbation corroboration | Fig 5 |
| R6 | **Neoadjuvant / longitudinal treatment designs** (pre / post / on-treatment biopsies) | Which clones and states are selected under therapy? What distinguishes responder vs non-responder baseline ecosystems? | **Gap** — no longitudinal paired-design contract; `ProjectContext` would need treatment/timepoint fields | **Not built; highest-value gap** | Longitudinal contract → own AK112/LPJ cohorts (ivonescimab PD-1/VEGF context is squarely this trend) | **Fig 5 main line candidate** (anchor lit: Yost 2019) |
| R7 | Immune repertoire integration (TCR/BCR + scRNA) | Does response come from pre-existing clone expansion or new recruitment? What is the fate of exhausted clones? | **Gap** — no repertoire support | Not built | Mount as external evidence (e.g. scRepertoire-style inputs), not in-house build | Fig 5 support |
| R8 | Spatial integration (scRNA as reference) | How are ecotypes organized in space? Where is the tumor–immune interface / invasive margin? | `tools/` spatial support layer — reference-provider role only (README boundary) | Skeleton, bounded | Deconvolution reference quality audit; collaborator spatial data | Fig 5 support |
| R9 | Drug-tolerant persisters / adaptive non-genetic resistance | Which states survive therapy and seed relapse? Reversible or fixed? | `therapy/resistance.py` + `evolution` skeleton; needs state-persistence tracking | Skeleton | State-persistence metric; longitudinal (R6) dependence | Fig 5 (drug-mechanism core) |
| R10 | Premalignancy / field cancerization | When does the microenvironment permit tumor emergence? What is already altered in adjacent-normal tissue? | No dedicated module; QC/analysis layers apply directly | Gap (no build needed for a case study) | Cohort with precursor/adjacent-normal sampling | Fig 5 alternative scenario |
| R11 | Therapeutic target prioritization | Which malignant program or TME interaction axis is druggable? | `therapy/target.py` | Skeleton | Prioritization recall vs. known targets | Fig 5 translational close |

**Explicitly out of scope (drift guards):**

- **Foundation models / virtual-cell predictors** (scGPT/scFoundation/GEARS
  class): conflicts with the lightweight-core, ecosystem-aware-not-replacing
  positioning; may be mounted as external annotation evidence at most.
- **Building a spatial platform**: scLucid provides scRNA references to
  spatial tools; it is not a spatial analysis platform (README boundary).
- **Wet-lab / assay method development**: outside the tool boundary.

**Literature status**: anchor citations above are established references;
representative 2024–2025 exemplars per row still need a literature
verification pass before any of these rows is cited in a manuscript.

## Open owner decisions (blocking finalization of this map)

- [ ] Tumor module tier confirmation (Part 0 table)
- [ ] Fill the **Own** leg for S1–S8 from real projects (or mark not-observed)
- [ ] Literature anchors for S8, S9 (or scope them out)
- [ ] Decide whether S10 is promoted to T1 (it has the strongest prior
      literature of all rows) or remains a T2 contract claim
- [ ] Pick the Fig 5 application direction (own clinical cohort vs public
      pan-cancer) — non-blocking for the gate
- [ ] Ratify Part 4: 11 research-question rows + 3 out-of-scope lines; confirm
      R6 (longitudinal treatment response, AK112/LPJ) as the Fig 5 main line
- [ ] Run the 2024–2025 literature verification pass for Part 4 rows
