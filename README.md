# scLucid: Lucid Tumor Single-Cell Interpretation System

[![PyPI version](https://badge.fury.io/py/sclucid.svg)](https://badge.fury.io/py/sclucid)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/yelu61/scLucid/actions/workflows/build.yml/badge.svg)](https://github.com/yelu61/scLucid/actions)

**scLucid** is a diagnostic-first, audit-ready Python framework for tumor
single-cell interpretation. The goal is to make single-cell analysis
**lucid**: clear in its assumptions, explicit in its inference boundaries, and
reviewable from raw AnnData to biological interpretation.

scLucid is not trying to replace Scanpy, Seurat, or the broader single-cell
ecosystem. It builds a tumor-focused research system around them: adaptive QC,
conservative preprocessing, evidence-based annotation, malignancy and TME
interpretation, explicit inference semantics, and optional bulk/spatial evidence
modules.

## 60-Second Quickstart

From a Cell Ranger output to a clustered, annotated AnnData with a shareable
HTML audit trail in four lines:

```python
import scLucid as scl

adata = scl.read_10x("path/to/filtered_feature_bc_matrix/", species="human")
adata = scl.run_pipeline(adata, dataset_type="pbmc_or_blood")
scl.export_audit_report(adata, "report.html")
```

`scl.read_10x` handles both Cell Ranger directories and `.h5` files, copies the
counts to `layers["counts"]` automatically, and attaches your dataset context
(species / tissue / cancer type) so downstream stages pick it up without extra
arguments.

## Why scLucid

| Principle | What It Means |
|-----------|---------------|
| **Diagnostic-first** | QC, preprocessing, DE, proportion, bulk, and spatial utilities are paired with checks and warnings before results are trusted. |
| **Audit-ready by default** | Decisions, parameters, warnings, contracts, reviewer tables, and review summaries are stored under `adata.uns["sclucid"]` and can be exported to an HTML audit report. |
| **Explicit inference semantics** | Results distinguish exploratory, descriptive, and sample-level inferences so exploratory signals are not overstated. |
| **Tumor ecosystem orientation** | Annotation, CNV/malignancy evidence, TME composition, therapy signatures, cell communities, and ecotype-style concepts are first-class design targets. |
| **Ecosystem-aware, not ecosystem-replacing** | Mature Python/R tools can be wrapped or validated when useful, but scLucid keeps a lightweight core and records method-specific evidence. |

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
adata_final = scl.run_pipeline(
    adata,
    stages=["qc", "preprocess", "analysis"],
    dataset_type="pbmc_or_blood",
    show_progress=True,
)

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

Implementation details, field names, design plans, and roadmaps live in
`docs/`, especially `docs/README.md`, `docs/api/`, and
`docs/user/best_practices.md`.

## Choose Your Analysis Mode

| Goal | Layer | Entry Point | Best For |
|------|-------|-------------|----------|
| One-line analysis | **Workflow** | `scl.run_pipeline()` | Beginners, standard projects, reproducible pipelines |
| Composable steps | **Simple API** | `scl.qc.calculate_qc_metric()`, `scl.pp.normalize_data()`, etc. | Analysts who need parameter control |
| Full transparency | **Advanced** | `examples/03_advanced_notebooks/Step1A-QC_Audit.ipynb` | Review-grade audits |

## Documentation

* **Quick Start**: [docs/user/quickstart.md](docs/user/quickstart.md)
* **Installation Guide**: [docs/user/installation.md](docs/user/installation.md)
* **Best Practices**: [docs/user/best_practices.md](docs/user/best_practices.md)
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
