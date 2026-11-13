# scLucid

### A question-driven framework for tumor single-cell research

[![PyPI version](https://badge.fury.io/py/sclucid.svg)](https://badge.fury.io/py/sclucid)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/yelu61/scLucid/actions/workflows/build.yml/badge.svg)](https://github.com/yelu61/scLucid/actions)

**Start with a research question. Connect study design, analysis and evidence.**

scLucid is a Python research framework that starts from tumor single-cell data
and connects biological and clinical questions to study design, reusable analysis
tasks and reviewable evidence. Its intended use is to help researchers decide
what their data can answer, select compatible analyses, inspect sample-level
results and identify the evidence needed for stronger conclusions. Spatial,
bulk, clinical and perturbation evidence can contribute to the same question
when their measurement and study designs are compatible.

**中文定位：** scLucid 是以单细胞数据为起点、面向肿瘤研究问题的分析与证据组织框架。它旨在帮助研究者将生物学或临床现象转化为可检验问题，判断数据与设计能回答什么，组合分析任务，并形成可追溯、可复核且解释边界明确的患者／动物级证据。

This is the product direction. Research-entry cards and cross-source evidence
rules are currently documentation-guided; they are not an automated question
translator or a completed multimodal workflow.

> **Development status:** this page presents the current framework direction
> and local development progress. The published GitHub code and detailed docs
> are an earlier snapshot. The Q1 example, Q2/Q4 additions and revised task
> cards described below have not yet been synchronized to GitHub; cloning this
> repository does not currently provide those additions. This README update
> is not a code release or a scientific-validation milestone.

[![scLucid v7: six research entries, conditional task selection and evidence returned to the original question](docs/assets/sclucid_question_driven_framework.png)](docs/assets/sclucid_question_driven_framework.png)

*Framework v7: R01–R06 define overlapping research contexts; Q1–Q9 supply shared,
bounded analytical evidence. The R02 example illustrates conditional selection,
not a mandatory pipeline. Optional Q4 supplies frozen program definitions for
Q1 scores; freezing is not independent validation. Q5 additionally requires
appropriate genetic/lineage evidence and a separately supported implementation.
Task roles and development maturity are separate. All plots are schematic.
Click the image to enlarge it, or open the
[vector source](docs/assets/sclucid_question_driven_framework.svg).*

**[Get started](#get-started)** · **[Research questions](#research-questions)** ·
**[Worked Q1 example](#run-a-research-task)** · **[Documentation](docs/README.md)** ·
**[Roadmap](docs/roadmap/index.md)**

## Research questions

Researchers can start from patient heterogeneity; treatment response, residual
disease and recurrence; ecosystem/spatial organization; progression and
metastasis; intervention priorities; or clinical measurement and stratification.
The local research-entry cards translate these overlapping contexts into design
requirements and task choices. R01/R02 have detailed local documentation;
R03–R06 are overview entries. Those detailed cards remain pending publication.
The definitions and bounded example below summarize their role.

**Research question → testable subquestions and design eligibility → selected Q
tasks → evidence synthesis answering the original question.** Q1–Q9 retain their
IDs as reusable analytical tasks. Their outputs return to the research entry
for interpretation, including unanswered parts and next validation.
The overview image shows six research contexts and one bounded R02 mapping
example. Entry-documentation status does not change task implementation maturity.

### Research entries and Q tasks have different responsibilities

| Layer | What it specifies | What counts as a useful result |
| --- | --- | --- |
| **Research entry (R01–R06)** | A biological/clinical concern, population, intended use, competing explanations and the evidence needed to answer it | An evidence-backed answer to the instantiated research question, with unresolved alternatives and a next step |
| **Analysis task (Q1–Q9)** | A bounded analytical question, eligible design/input, reusable methods, output meaning and interpretation limits | Correct, reviewable estimates or summaries for the specified task/version/design branch; unsupported or non-estimable results remain visible |

**中文边界：研究入口定义“在什么研究情境下，要回答什么、需要哪些证据”；Q 任务定义“在明确的数据与设计条件下，能够产出哪一种分析证据”。** Q 是分析任务的组织与验证单元；研究场景是选择、组合和评价这些任务的情境。场景中的研究问题仍须具体化，不能用一个场景名称替代研究设计。

The relationship is **conditional and many-to-many**. One entry may use one or
several Q tasks; the same Q can support several entries. A task is not one
algorithm, and an entry is not a second implementation of the tasks. Passing
each selected Q does not establish that the original research question is
answered. “Clinical/scenario layer” refers to R entries here; C1–C3 are reserved
for candidate paper contributions, not entry IDs.

For example, **“Which states are associated with residual disease or differences
in treatment response?”** must first be decomposed:

| Testable subquestion | Task choice and condition | Boundary |
| --- | --- | --- |
| Do captured-cell composition or fixed within-population scores differ between eligible conditions? | Q1; the local example supports independent or genuine paired two-condition designs | A pre/post difference alone does not identify a treatment effect or resistance mechanism |
| Are new program definitions needed, and do candidate programs recur across patients? | Q4 prototype only when discovery is needed; frozen program definitions supply Q1 scores | Freezing is not independent validation; discovery patients are not independent validation patients |
| Do the measured state and microenvironment features co-vary across samples? | Optional Q2 prototype with compatible unit/feature definitions | Sample association does not establish spatial contact or signaling |
| Are state differences related to genetic clones or lineage labels? | Conditional Q5 with appropriate evidence; not a complete executable task | Expression snapshots alone do not establish state conversion |

Without response labels, the response-association part remains unanswered;
without valid pairing, a within-patient change claim is unavailable. The final
review states which parts were supported, contradicted or not assessed, rather
than reporting only a list of completed tasks. Full design branches and maturity are maintained in local entry/task cards
pending publication; the image and this overview do not add executable tasks.

### Reusable Q task catalog

| Task | Question it is designed to address |
| --- | --- |
| **Q1 · Composition and state change** | Do conditions or real follow-up visits differ in cellular composition, within-population programs, or both? |
| **Q2 · Ecosystem relationships** | Which tumor and microenvironment features co-vary across samples, and which relationships replicate? |
| **Q3 · Tissue references** | How do tumor, adjacent and healthy tissue differ under a compatible reference design? |
| **Q4 · Recurrent programs** | Which expression programs recur across patients and independent studies? |
| **Q5 · Clone–state relationships** | How do genetic clones, receptor clonotypes or lineage labels relate to expression states? |
| **Q6 · Model fidelity and transport** | Which human disease features are represented in a model, and which model findings transfer to humans? |
| **Q7 · Target prioritization** | Which targets merit experiments given patient coverage and normal-tissue expression risk? |
| **Q8 · Perturbation response** | How might a specified intervention affect a defined population, and which predictions merit experimental testing? |
| **Q9 · Clinical feature transport** | Do frozen single-cell features transfer to independent clinical cohorts and retain their intended outcome relevance? |

**Current scope at a glance — this catalog does not imply nine completed workflows:**

- **Q1 · Worked example in local development:** independent or genuinely paired two-condition
  comparisons, with scoped CRC application and reuse evidence.
- **Q2 / Q4 · Bounded local prototypes:** ecosystem relationships and recurrent
  programs, respectively; independent scientific validation remains unfinished.
- **Q3 / Q5 / Q6 · Later or external-evidence tasks:** reference comparisons,
  clone–state relationships and model fidelity require further task development.
- **Q7–Q9 · Planned:** no complete executable tasks declared. Virtual-cell
  models are candidate methods within Q8, subject to context-specific evaluation.

Detailed maturity, supported branches and stop rules are maintained in the
local task cards pending publication. Q1's scoped execution and reuse also
remain distinct from independent scientific validation.

## How the framework works

1. **Define the research question and its subquestions.** State the intended
   use, competing explanations and required evidence. Declare biological units,
   specimens, conditions, real pairing/time and the identifiable comparison.
2. **Prepare eligible inputs.** Review counts, QC, expression representations
   and cell identities. Preserve provenance when using externally prepared data.
3. **Select and combine tasks.** Reuse methods appropriate to the design. For
   example, Q4 can define program features for Q1, while Q2 can examine their
   sample-level relationships with the microenvironment.
4. **Return the evidence to the research question.** Read effects, uncertainty, coverage, sensitivity
   and non-estimable results alongside linked source tables and figures.
   State the supported answer, unresolved alternatives and missing evidence.
   Add independent functional, spatial or clinical evidence where required.

Tasks share QC, preprocessing, annotation, statistical analysis, tumor-program
and visualization capabilities. Mature ecosystem tools can supply optional
methods and evidence. Task composition depends on compatible inputs and
measurement definitions; it does not require every task or modality.

## What scLucid brings to a project

| Principle | Practical consequence |
| --- | --- |
| **Question-led analysis** | Study design and intended outputs guide which methods to use. |
| **Composition and state stay distinct** | Population abundance and within-population expression answer separate questions and retain their own denominators. |
| **Biological units stay explicit** | Patients or animals, specimens and cells have different roles; more cells do not create more independent replicates. |
| **Reviewable decisions** | QC and preprocessing expose a recommendation before an explicit policy is applied, with parameters, warnings and provenance retained. |
| **Evidence-linked interpretation** | Task outputs connect estimates and figures to source rows, sensitivities and unresolved evidence gaps. |
| **Reusable foundations** | Stepwise examples and existing APIs support project-specific analyses; optional dependencies are added when needed. |

These are design commitments and scoped implementation features. Comparative
benefits require their own evaluation; no universal performance advantage is
claimed.

## Get started

### Install

For an editable installation of the published GitHub snapshot:

```bash
git clone https://github.com/yelu61/scLucid.git
cd scLucid
pip install -e ".[analysis,de]"
```

For a published package installation:

```bash
pip install sclucid
```

The Python import is **`import scLucid as scl`** (capital `L`), matching the
source package directory `src/scLucid/`.

Optional extras cover annotation and clustering (`analysis`), differential
expression (`de`), experimental cNMF programs (`programs`, local development), integration,
CNV, velocity, bulk and spatial support. Install only those required by the
chosen workflow; see the [installation guide](docs/user/installation.md).
Uncommitted local work is available only in its local editable checkout, and
published packages may differ from the current development documentation.

The `velocity` extra supports optional RNA-velocity methods. Their possible
role in candidate state-transition analyses combining Q1/Q4/Q5
is a design direction, not a validated transition workflow. Velocity or
pseudotime alone does not establish lineage, reversibility or resistance.

### Prepare a dataset

Start from an AnnData object with verified count provenance and real sample
metadata. In this example, `sample`, `patient` and `condition` are existing
columns in `adata.obs`; replace them with the corresponding project columns.
The input must retain raw counts in `layers["counts"]`.

```python
import scanpy as sc
import scLucid as scl

adata = sc.read_h5ad("data/cohort_counts.h5ad")
context = scl.ProjectContext(
    dataset_type="tumor_tissue",
    species="human",
    sample_key="sample",
    experimental_unit_key="patient",
    condition_key="condition",
    input_provenance="filtered_counts",
)

# Review the recommendation, reason and required next action.
qc_review = scl.recommend_qc_policy(adata, context=context)
print(qc_review.status, qc_review.reason, qc_review.next_action)
```

After reviewing the recommendation and resolving any blocking issue, explicitly
apply the policy and review the next stage:

```python
qc_result = scl.apply_qc_policy(adata, qc_review.policy)

pp_review = scl.recommend_preprocess_policy(
    qc_result.adata, context, consumer="exploration"
)
print(pp_review.status, pp_review.reason, pp_review.next_action)

# Apply after reviewing this stage's policy and unresolved issues.
pp_result = scl.apply_preprocess_policy(qc_result.adata, pp_review.policy)
```

Declare `paired_key` when genuine pairing exists and `batch_key` only for a
known technical batch. `recommend_*` leaves its input unchanged; `apply_*`
checks the reviewed input fingerprint and returns the output in `result.adata`.
The exploration representation above is preparation for inspection; each
downstream statistical method still needs its appropriate expression input.
See the [review/apply workflow](docs/user/evidence_calibrated_workflow.md) and
[data contracts](docs/user/data_contracts.md).

### Run a research task

**Local development preview:** this example and its task card are not yet in
the published GitHub snapshot. The command and output layout below describe
the development checkout containing `examples/02_key_api/q1_composition_state.py`.

The Q1 worked example separates
captured-cell composition from fixed within-population score changes. It is a
stepwise Python script driven by a JSON study specification; the command below
runs those same steps from the terminal.

To continue from `pp_result.adata`, first review the cell-population annotations
and compute the prespecified program scores using an appropriate expression
representation. Store the population labels and per-cell scores in the prepared object’s
`.obs`, then save that annotated object as H5AD. In the JSON specification, point
`input_path` to this file, map the metadata and score columns, and record the
score definitions and provenance. Preprocessing alone does not supply these
annotations or scores, and the Q1 script does not infer them.

In that development checkout, read the Q1 task card in
`docs/CLAIM_EVIDENCE_MAP.md` before creating the project specification.
Set `SCLUCID_Q1_OUT` to a new or empty
directory; existing results are protected from overwrite.

```bash
SCLUCID_Q1_SPEC=/path/to/q1_spec.json \
SCLUCID_Q1_OUT=/path/to/new_results \
python examples/02_key_api/q1_composition_state.py
```

A completed run writes the following linked outputs directly to
`SCLUCID_Q1_OUT`:

```text
new_results/
├── sample_design.csv               # Specimens, biological units and conditions
├── population_features.csv         # Sample/population abundance, scores and eligibility
├── abundance_state_effects.csv     # Separate abundance/state rows, intervals and test status
├── sensitivity_full.csv            # Full-data comparison for the sensitivity analysis
├── sensitivity_omissions.csv       # Every leave-one-unit-out result, including failures
├── sensitivity_summary.csv         # Effect shifts and influential omissions
├── specification.json              # Study specification used for this run
├── figure_legends.txt              # Design, denominators and interpretation limits
├── sample_observations_*.{pdf,svg,png}
├── effects_*.{pdf,svg,png}
└── run.json                        # Execution record, provenance and output hashes
```

Here `*` denotes a contrast/population/score page identifier; each figure page
is saved in all three formats. Abundance and state estimates are separate rows
in one effects table. Coverage and eligibility are retained in the source
tables and figure-page records. `run.json` records execution completion
separately from scientific performance (`NOT_EVALUATED`).

The current entrypoint supports
independent or genuinely paired two-condition comparisons. Adjusted,
multi-visit and treatment-benefit models are outside that entrypoint's scope.
For composition alone, use the proportion/comparison APIs without inventing a
state score.

For the currently published interactive examples, see the
[notebook guide](docs/user/notebooks.md). Existing plan/pipeline
and stage APIs remain available with their own execution contracts; see
[usage layers](docs/user/usage_layers.md).

## Evidence and development boundaries

scLucid is under active development. A runnable example, an engineering check
and independent scientific validation establish different things. QC and
preprocessing remain **`REVIEW`** under their existing scientific and
real-project usability gates; scoped task execution does not promote the
whole framework to validated status.

- **Associations retain their meaning.** Co-variation does not establish
  functional communication, spatial proximity or causality.
- **Expression retains its limits.** Cell identity alone does not establish
  malignancy; an RNA program or target score does not establish drug response,
  a therapeutic window or clinical utility.
- **Additional evidence is conditional.** Spatial, lineage, protein,
  perturbation and clinical-outcome data extend a task only when available,
  compatible and relevant to the claim.
- **Generalization is evaluated separately.** Human and mouse studies,
  discovery and independent evaluation, and measurement transfer and outcome
  prediction require their own evidence.

Current implementation and evaluation records are maintained in the
[claim–evidence map](docs/CLAIM_EVIDENCE_MAP.md),
[maturity guide](docs/user/qc_preprocess_maturity.md) and
[roadmap](docs/roadmap/index.md).

## Documentation

These links lead to the currently published documentation snapshot. Its
positioning, roadmap and task coverage predate the local development described
above and will be synchronized with the corresponding code.

| Start here | Reference |
| --- | --- |
| Learn the framework | [Documentation home](docs/README.md) · [Positioning and design](docs/SCLUCID_STRATEGIC_IMPLEMENTATION_PLAN.md) |
| Choose a question | [Q1–Q9 overview and scope](#research-questions) |
| Prepare and inspect data | [Quickstart](docs/user/quickstart.md) · [Project context](docs/user/project_context.md) · [Reviewing results](docs/user/reviewing_results.md) |
| Work through an analysis | [Examples](docs/user/examples.md) · [Notebooks](docs/user/notebooks.md) · [Best practices](docs/user/best_practices.md) |
| Use specific functions | [API reference](docs/api/) · [Data contracts](docs/user/data_contracts.md) |
| Check scope and progress | [Roadmap](docs/roadmap/index.md) · [Documentation policy](docs/CURRENT_IMPLEMENTATION_AND_DOCS_POLICY.md) |

## Contributing

Contributions that improve task usability, evidence quality, reproducibility or
documentation are welcome. See the [contributing guide](docs/CONTRIBUTING.md)
and [issue tracker](https://github.com/yelu61/scLucid/issues).

## Citation and license

Until a formal methods paper is available, cite the repository and report the
package version or commit used for your analysis. A manuscript citation will be
added when available.

scLucid is released under the [MIT License](LICENSE).
