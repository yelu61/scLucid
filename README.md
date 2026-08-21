# scLucid: Lucid Tumor Single-Cell Interpretation System

[![PyPI version](https://badge.fury.io/py/sclucid.svg)](https://badge.fury.io/py/sclucid)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/yelu61/scLucid/actions/workflows/build.yml/badge.svg)](https://github.com/yelu61/scLucid/actions)

**scLucid** is an evidence-calibrated review, strategy-comparison, and
reproducible execution system for real multi-sample, especially tumor,
single-cell RNA-seq projects. It is not a large automatic analysis package and
does not replace biological judgment. Its job is to expose where common
strategies can fail, compare their consequences, and return an explicit
`READY / REVIEW / BLOCKED` decision with one next action.

The current product focus is QC and preprocessing. Analysis and tumor modules
remain available for compatibility, but new feature development is frozen until
the locked QC/preprocessing scientific and real-project usability gates pass.
Until then, scLucid does **not** claim universal superiority over traditional
QC or preprocessing workflows.

## 60-Second Quickstart

The standard path has four explicit actions. Review does not modify the input;
execution requires the reviewed policy:

```python
import scLucid as scl

adata = scl.read_10x("path/to/filtered_feature_bc_matrix/", species="human")
context = scl.ProjectContext(
    dataset_type="tumor_tissue",
    sample_key="sample",
    input_provenance="filtered_counts",
)

qc_review = scl.recommend_qc_policy(adata, context=context)
qc_result = scl.apply_qc_policy(adata, qc_review.policy)

pp_review = scl.recommend_preprocess_policy(
    qc_result.adata,
    context,
    consumer="exploration",
)
pp_result = scl.apply_preprocess_policy(qc_result.adata, pp_review.policy)
```

`scl.read_10x` handles both Cell Ranger directories and `.h5` files, copies the
counts to `layers["counts"]` automatically, and attaches your dataset context
(species / tissue / cancer type) so downstream stages pick it up without extra
arguments.

The compatibility workflows (`run_qc`, `run_preprocessing`, and `run_pipeline`)
remain available for one minor release. See the
[evidence-calibrated workflow](docs/user/evidence_calibrated_workflow.md) for
the new contract and migration boundary.

## Why scLucid

| Principle | What It Means |
|-----------|---------------|
| **Decision-first** | Review is read-only; a fingerprinted policy is applied only after an explicit decision. |
| **Audit-ready by default** | Decisions, parameters, warnings, contracts, reviewer tables, and review summaries are stored under `adata.uns["sclucid"]` and can be exported to an HTML audit report. |
| **Explicit inference semantics** | Results distinguish exploratory, descriptive, and sample-level inferences so exploratory signals are not overstated. |
| **Tumor project orientation** | Sample failure, fragile lineages, stress, ambient RNA, doublets, and integration confounding are separate evidence heads. |
| **Ecosystem-aware, not ecosystem-replacing** | Mature Python/R tools can be wrapped or validated when useful, but scLucid keeps a lightweight core and records method-specific evidence. |
| **Context-aware decisions** | Dataset context, tumor biology, and analysis intent are treated as part of the decision record rather than hidden analyst assumptions. |

## Installation

The toolkit is modular. Install the lightweight core and add extras as needed.

```bash
# Standard installation (core QC, preprocessing, analysis, and plotting)
pip install sclucid

# Optional analysis extras (CellTypist, cosg, etc.)
pip install "sclucid[analysis]"

# Optional advanced tools (scVelo, infercnvpy, scVI-related workflows)
pip install "sclucid[tools]"

# Optional bulk RNA-seq evidence modules
pip install "sclucid[bulk]"

# Optional spatial transcriptomics evidence modules
pip install "sclucid[spatial]"

# Install everything
pip install "sclucid[all]"
```

Install the latest development version from GitHub:

```bash
pip install "git+https://github.com/yelu61/scLucid.git"
```

For editable development:

```bash
git clone https://github.com/yelu61/scLucid.git
cd scLucid
pip install -e ".[all]"
```

## A 5-Minute Example

```python
import scLucid as scl

# Load data (Cell Ranger dir, .h5, or .h5ad)
adata = scl.read_10x(
    "data/pbmc3k/filtered_feature_bc_matrix/",
    species="human",
    tissue="PBMC",
)

# Run the supported core workflow
context = scl.ProjectContext(
    dataset_type="pbmc_or_blood",
    sample_key="sample",
    study_objective="broad cell atlas",
)
plan = scl.plan_analysis(adata, context=context)
adata_final = scl.run_pipeline(
    adata,
    plan=plan,
    show_progress=True,
)

# Read the prioritized decisions before interpreting results
review = scl.review_run(adata_final)
print(review.to_frame())
print(review.show_next_actions())

# Export an auditable HTML report
scl.export_audit_report(adata_final, "results/audit_report.html")

# Visualize results with a publication-ready style
from scLucid import FONT_NATURE
scl.set_figure_params(dpi=300, font_style=FONT_NATURE)
scl.pl.plot_embedding(adata_final, color_by="cell_type_auto", show=False)

import matplotlib.pyplot as plt
plt.savefig("results.pdf", dpi=600, bbox_inches="tight")
```

## Development Boundary

scLucid is built as an evidence-driven tumor single-cell workflow system, not as
a claim that every automated result is superior to Scanpy, Seurat, scran,
inferCNV, CopyKAT, CellTypist, or other mature tools.

The stable promise is narrower: make assumptions visible, preserve biological
caution, label inference boundaries, and keep enough review evidence for another
analyst to audit the workflow.

Portable scientific-reasoning and researcher-perspective skills remain
upstream of the package and evolve independently. They may provide a neutral
`ReasoningBrief`, but they are not scLucid dependencies and their identities or
opinions are never treated as execution evidence. See the
[scientific reasoning contracts](docs/user/scientific_reasoning_contracts.md)
for the draft boundary.

Current module maturity, implementation plans, validation status, and roadmap
details live in `docs/`, especially `docs/README.md`, `docs/user/`,
`docs/api/`, and `docs/roadmap/`.

## Core Differentiators

scLucid's long-term value is not a larger list of wrappers. Its core
competition is interpretability and scientific restraint around tumor
single-cell workflows:

- **Reviewable decisions**: QC, preprocessing, and analysis steps leave behind
  structured review records that explain what was recommended, what was applied,
  and what still needs human review.
- **Tumor-aware caution**: high mitochondrial content, stress programs,
  doublets, batch correction, malignancy calls, and TME signals are treated as
  biological-risk decisions rather than simple automatic filters.
- **Layer and inference contracts**: expression layers, embeddings,
  annotations, and DE/proportion results are labeled with their intended use so
  exploratory cell-level findings are not mistaken for sample-level claims.
- **Evidence bridges**: single-cell results are designed to connect with bulk,
  spatial, marker, CNV, and external reference evidence without hiding the
  confidence boundary.
- **Falsifiable workflow records**: outputs should make it possible to ask what
  would change under a different QC, preprocessing, annotation, or inference
  choice, instead of presenting one automated path as final.

Implementation details, field names, design plans, and roadmaps live in
`docs/`, especially `docs/README.md`, `docs/api/`, and
`docs/user/best_practices.md`.

## Choose Your Analysis Mode

| Goal | Layer | Entry Point | Best For |
|------|-------|-------------|----------|
| Guided first pass | **Workflow** | `scl.plan_analysis()`, `scl.run_pipeline()`, `scl.review_run()` | Beginners, standard projects, reproducible baselines |
| One stage | **Stage workflow** | `scl.run_qc()`, `scl.pp.run_preprocessing()`, `scl.analysis.run_standard_analysis()` | Analysts rerunning one decision boundary |
| Composable steps | **Simple API** | `scl.qc.calculate_qc_metric()`, `scl.pp.normalize_data()`, etc. | Analysts who need parameter control |
| Full transparency | **Advanced** | `examples/03_advanced_notebooks/Step1A-QC_Audit.ipynb` | Review-grade audits |

## Documentation

* **Quick Start**: [docs/user/quickstart.md](docs/user/quickstart.md)
* **Installation Guide**: [docs/user/installation.md](docs/user/installation.md)
* **Best Practices**: [docs/user/best_practices.md](docs/user/best_practices.md)
* **Project Context**: [docs/user/project_context.md](docs/user/project_context.md)
* **Reviewing Results**: [docs/user/reviewing_results.md](docs/user/reviewing_results.md)
* **Parameter Profiles**: [docs/user/parameter_profiles.md](docs/user/parameter_profiles.md)
* **Core Data Contracts**: [docs/user/data_contracts.md](docs/user/data_contracts.md)
* **API Reference**: [docs/api/](docs/api/)
* **Strategic Plan**: [docs/SCLUCID_STRATEGIC_IMPLEMENTATION_PLAN.md](docs/SCLUCID_STRATEGIC_IMPLEMENTATION_PLAN.md)
* **Roadmap**: [docs/roadmap/](docs/roadmap/)
* **Plugin Development**: [docs/dev/PLUGIN_DEVELOPMENT_GUIDE.md](docs/dev/PLUGIN_DEVELOPMENT_GUIDE.md)
* **R Parity Matrix**: [docs/user/r_parity.md](docs/user/r_parity.md)
* **PBMC Golden Path**: [scripts/run_pbmc_golden_path.py](scripts/run_pbmc_golden_path.py)
* **Analysis Acceptance Runner**: [scripts/run_analysis_acceptance.py](scripts/run_analysis_acceptance.py)

## What scLucid Is Not

- It is not a claim that every scLucid result is automatically better than
  Scanpy, Seurat, or specialist tools.
- It is not a broad multi-omics platform that tries to cover every published
  method.
- It is not a general spatial transcriptomics platform.
- It is not a black-box automated annotation engine.

scLucid's bet is narrower and more useful: make tumor single-cell workflows
clearer, safer to interpret, easier to audit, and easier to connect with bulk,
spatial, clinical, and mature ecosystem evidence.

## Contributing

We welcome contributions from the community! Please check out our
[Contributing Guidelines](CONTRIBUTING.md) and the issue tracker.

## License

`scLucid` is licensed under the MIT License.

## How to Cite

If you use `scLucid` in your research before a formal methods paper is available,
please cite the GitHub repository and include the package version used in your
analysis. A manuscript citation will be added once available.
