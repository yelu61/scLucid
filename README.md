# scLucid: A Comprehensive System for Single-Cell Analysis

[![PyPI version](https://badge.fury.io/py/sclucid.svg)](https://badge.fury.io/py/sclucid)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/yelu61/scLucid/actions/workflows/build.yml/badge.svg)](https://github.com/yelu61/scLucid/actions)

**scLucid** is a powerful and flexible Python toolkit for the analysis of single-cell RNA-sequencing data. It is designed to be more than just a wrapper; it's a complete **analysis system** that guides researchers from raw data to deep biological insights.

The toolkit's philosophy is to balance ease-of-use for standard workflows with deep customizability for advanced, exploratory research. It achieves this through a modular architecture, high-level workflow functions, and a unique, biology-aware marker management system.

### 60-Second Quickstart

From a Cell Ranger output to a clustered, annotated AnnData with a shareable HTML audit trail — four lines:

```python
import scLucid as scl

adata = scl.read_10x("path/to/filtered_feature_bc_matrix/", species="human")
adata = scl.run_pipeline(adata, dataset_type="pbmc_or_blood")
scl.export_audit_report(adata, "report.html")
```

`scl.read_10x` handles both Cell Ranger directories and `.h5` files, copies the counts to `layers["counts"]` automatically, and attaches your dataset context (species / tissue / cancer type) so downstream stages pick it up without extra arguments. For an existing `.h5ad` file, use `scl.read_h5ad` instead.

### Project Status

scLucid is in active development and is best described as an **evidence-driven
single-cell workflow system in late prototype / early hardening stage**.

The package has moved beyond a collection of wrappers. It already contains
stable workflow entrypoints, AnnData contracts, review summaries, marker-resource
routing, recommendation scaffolds, real-data golden-path scripts, and targeted
tests for QC, preprocessing, analysis, tumor utilities, resources, and reporting.
The strongest modules today are QC and preprocessing: both are close to
benchmark-module maturity for auditability, reproducibility, and workflow fit.
Analysis is the active module being raised to the same standard. It now has an
evidence-first closed loop for clustering-resolution review, marker discovery,
marker-manager/CellTypist/LLM annotation evidence, consensus labels, optional
malignancy interpretation, and analysis review-summary maturity contracts.

The current development boundary is important: scLucid can already provide a
traceable, biologically informed workflow, but it should not yet claim broad
scientific superiority over Scanpy, Seurat, scran, inferCNV, CopyKAT, CellTypist,
or other mature tools. The next stage is real-data workflow hardening: run the
same auditable path on PBMC, PDAC, and active tumor projects; compare outputs to
standard workflows; and turn those results into documented acceptance criteria.

### Current Maturity Assessment

| Area | Current Level | What Works Now | Main Gaps |
|------|---------------|----------------|-----------|
| QC | Candidate benchmark module | Adaptive thresholds, tumor-aware warnings, doublet heuristics, review summaries, benchmark scaffolds | Broader real-data benchmarks and clearer user-facing threshold narratives |
| Preprocessing | Candidate benchmark module | Layer contracts, normalization/HVG/PCA/neighbors/UMAP evidence, batch-correction cautions, maturity contract | Larger multi-sample validation, stronger batch-correction recommendation evidence |
| Analysis | Second benchmark module in active hardening | `clustering_review -> markers -> annotation_evidence -> annotation_consensus -> malignancy_interpretation`, manager-routed marker resources, review-summary contract | Real-data acceptance runs, richer CellTypist/reference evidence, better human-facing review tables |
| Marker Resources | Strong architectural direction | Unified `Manager`, human/mouse registry resources, tissue/tumor marker views, artifact/program/tumor routing, curation SOP | Source provenance at scale, mouse tissue/tumor parity, atlas-derived marker review |
| Tumor Module | Feature-rich but needs integration hardening | CNV, malignancy scoring/classification, TME, therapy, heterogeneity, workflow scaffolds | Consume stable analysis outputs more tightly, store tumor-stage review summaries, validate on tumor datasets |
| Plotting | Useful foundation | Publication-style themes and domain plots | Top-journal figure templates, richer multi-panel reports, visual regression checks |
| Tools / R Parity | Broad wrapper coverage | Python-facing wrappers for mature ecosystem methods | Dependency isolation, parity matrices, realistic fallbacks, method-specific validation |
| Documentation / Examples | Good skeleton | Three usage layers, advanced notebooks, golden-path scripts | Keep docs synchronized with maturity contracts and real-data acceptance results |

### Development Roadmap

The roadmap is intentionally staged so the package matures from traceable
execution to evidence-backed biological usefulness.

**Phase 1 — Finish Analysis As The Second Benchmark Module**

- Harden the evidence-first `run_standard_analysis` path:
  clustering-resolution evidence, marker discovery, marker-manager annotation
  evidence, optional reference/CellTypist evidence, optional data-driven LLM
  suggestion bundles, consensus labels, optional malignancy interpretation, and
  review summary.
- Keep first-pass annotation conservative: lineage / major cell type first;
  subtype and state annotation should be driven by subset reclustering or explicit
  user request.
- Route all marker-dependent analysis through `get_marker_manager()` views:
  `lineage_annotation`, `subtype_annotation`, `state_annotation`,
  `artifact_annotation`, `program_scoring`, and `tumor_interpretation`.
- Treat LLM output as annotation evidence, not ground truth.

**Phase 2 — Tumor-Aware Interpretation Contract**

- Keep `analysis.run_malignancy_interpretation` as a lightweight bridge that
  consumes final annotation, tumor marker evidence, optional CNV scores, optional
  malignancy signatures, and user-provided cancer context.
- Keep heavy tumor-specific algorithms in `scLucid.tumor`: CNV inference,
  malignancy scoring/classification, TME, therapy, heterogeneity, and ecosystem
  workflows.
- Separate normal epithelial annotation from malignant-cell interpretation.
- Support multiple evidence backends: lightweight CNV score, inferCNV-style
  output, CopyKAT-like calls, malignancy signatures, and manual evidence.
- Store malignant/non-malignant/suspect/unresolved calls with confidence,
  reasons, and review requirements.

**Phase 3 — Resource Curation And Validation**

- Continue upgrading marker resources from “readable” to “routable,
  reviewable, source-aware”.
- Add mouse tissue/tumor marker parity after the human route stabilizes.
- Add resource validation tests for required metadata, marker symbol hygiene,
  view routing, negative markers, artifact exclusion, and tumor evidence
  isolation.
- Curate immune and tumor-state markers from pan-cancer atlases while keeping
  broad pathway signatures in gene-set JSON/GMT resources rather than concise
  annotation TOML files.

**Phase 4 — Real-Data Acceptance Gates**

- Maintain PBMC as the normal baseline.
- Maintain PDAC as the first tumor acceptance workflow.
- Add at least one second tumor type and one active research project notebook.
- Record acceptance criteria for cell retention, preprocessing readiness,
  cluster interpretability, annotation confidence, marker consistency,
  malignancy evidence, and report completeness.

**Phase 5 — Publication Output And User Experience**

- Convert advanced notebooks into polished, reproducible workflow narratives.
- Expand audit reports to include analysis and tumor interpretation maturity.
- Add top-journal figure templates and visual regression checks for important
  plotting functions.
- Keep beginner workflow, simple API, and advanced expert routes synchronized.

### Key Features

* **🧪 End-to-End Workflows**: High-level functions like `run_standard_qc` and `run_preprocessing` to go from raw data to a clustered UMAP with just a few lines of code.
* **🧠 Intelligent QC**: Advanced doublet detection using a combination of `scrublet` and a novel, marker co-expression heuristic.
* **🧬 Biology-Aware Analysis**: A unified `Manager`/resource system routes curated markers and gene sets into compartment, lineage, subtype, state, artifact, program-scoring, and tumor-interpretation views.
* **🔬 Multi-Evidence Annotation**: A complete suite of tools to annotate cell types using automated methods (`CellTypist`), gene scoring, and evidence-gathering functions (`characterize_clusters`).
* **🔧 Advanced Tools Module**: Seamlessly integrated wrappers for specialized analyses, including:
    * RNA Velocity (`scVelo`)
    * CNV Inference (`infercnvpy`, `CopyKAT`)
    * Trajectory Inference (`PAGA`, `Monocle3`)
    * Cell-Cell Communication (`CellChat`, `CellPhoneDB`)
    * Bulk Deconvolution (`BayesPrism`, `DWLS`)
* **📊 Publication-Quality Visualizations**: A rich plotting library to generate stunning and informative figures for every step of the analysis.
* **🎨 Academic Journal Font Styles**: Pre-configured font styles for top journals - Nature (Arial), Cell (Helvetica), and Traditional (Times New Roman).
* **🔄 Reproducible Science**: A configuration-driven approach using **Pydantic** ensures automatic validation, type safety, and reproducibility with JSON serialization.
* **📝 Auditable Reports**: `scl.export_audit_report(adata, "report.html")` renders every recommendation rationale, applied threshold, configuration lineage, and contract validation result into one self-contained HTML page — review-ready out of the box.
* **🔌 Extensible Plugin Architecture**: Abstract base classes and factory pattern allow you to create custom analysis plugins without modifying core code. See [Plugin Development Guide](docs/PLUGIN_DEVELOPMENT_GUIDE.md) for details.

### Choose Your Analysis Mode

scLucid offers **three user-facing layers** designed for different levels of control and expertise:

| Your Goal | Recommended Layer | Entry Point | Best For |
|-----------|-------------------|-------------|----------|
| **One-line analysis** — load data and run the full pipeline | **Workflow** | `scl.run_pipeline()` | Beginners, standard projects, reproducible pipelines |
| **Composable steps** — inspect or replace individual stages | **Simple API** | `scl.qc.calculate_qc_metric()`, `scl.pp.normalize_data()`, etc. | Analysts who need parameter control |
| **Full transparency** — every threshold, diagnostic, and override visible | **Advanced** | `examples/03_advanced_notebooks/Step1A-QC_Audit.ipynb` | Real exploratory projects, review-grade audits |

> **💡 How to choose**: If you just want results, use **Workflow**. If you need to tweak parameters, use **Simple API**. If you are doing research where every decision must be auditable, use **Advanced**.

**Examples for each layer:**
- **Workflow**: `examples/01_workflow/basic_pipeline.py`
- **Simple API**: `examples/02_simple_api/qc_step_by_step.py`
- **Advanced**: `examples/03_advanced_notebooks/Step1A-QC_Audit.ipynb` -> `Step1B-Preprocessing_Audit.ipynb` -> `Step2-Annotation_and_Malignancy.ipynb`

### Installation

The toolkit is modular. You can install the lightweight core and add extras as needed.

```bash
# Standard Installation (Core QC, Preprocessing, Analysis, and Plotting)
pip install sclucid

# To include additional analysis packages (CellTypist, cosg, etc.)
pip install "sclucid[analysis]"

# To include all advanced tools (scVelo, rpy2, infercnvpy, etc.)
pip install "sclucid[tools]"

# To install everything
pip install "sclucid[all]"
```

You can also install the latest development version directly from GitHub:
```bash
pip install "git+https://github.com/yelu61/scLucid.git"
```

For developers, clone the repository and install in editable mode:
```bash
git clone https://github.com/yelu61/scLucid.git
cd scLucid
pip install -e ".[all]"
```

On the local development machine, the maintained single-cell environment can run
the lightweight gates directly:

```bash
MAMBA_EXE=/opt/homebrew/bin/mamba \
SCLUCID_TEST_ENV_PATH=/Users/luye/micromamba/envs/scrna-env \
scripts/run_test_gates.sh
```

The first real-data workflow gate is the PBMC golden path:

```bash
/Users/luye/micromamba/envs/scrna-env/bin/python \
  scripts/run_pbmc_golden_path.py \
  --n-cells 300 \
  --output-dir results/golden/pbmc3k_subset \
  --overwrite
```

### Quick Start: A 5-Minute Analysis

Here is a minimal example of a complete workflow.

```python
import scLucid as scl

# --- 1. Load data (Cell Ranger dir, .h5, or .h5ad) ---
adata = scl.read_10x(
    "data/pbmc3k/filtered_feature_bc_matrix/",
    species="human",
    tissue="PBMC",
)

# --- 2. Run the supported core workflow ---
adata_final = scl.run_pipeline(
    adata,
    stages=["qc", "preprocess", "analysis"],
    dataset_type="pbmc_or_blood",
    show_progress=True,
)

# --- 3. Export an auditable HTML report ---
# Every threshold, parameter source, and warning is rendered into one file.
scl.export_audit_report(adata_final, "results/audit_report.html")

# --- 4. Visualize Final Results ---
# Set publication-ready font style for your target journal
from scLucid import FONT_NATURE, FONT_CELL, FONT_TRADITIONAL
scl.set_figure_params(dpi=300, font_style=FONT_NATURE)  # For Nature/Science

scl.pl.plot_embedding(adata_final, color_by="cell_type_auto", show=False)

# Save with embedded fonts for publication
import matplotlib.pyplot as plt
plt.savefig("results.pdf", dpi=600, bbox_inches="tight")
```

### Documentation

For detailed tutorials, how-to guides, and the full API reference:

* **Plugin Development**: [Plugin Development Guide](docs/PLUGIN_DEVELOPMENT_GUIDE.md) - Create custom analysis plugins
* **R Parity Matrix**: [docs/source/r_parity.rst](docs/source/r_parity.rst) - What each R-package port (BayesPrism, Monocle3, CellChat, DWLS) covers vs the R original
* **Naming Conventions**: [Naming Conventions](docs/NAMING_CONVENTIONS.md) - Code style guidelines
* **Local Documentation Source**: [docs/source/](docs/source/) - Sphinx documentation sources for installation, quickstart, API references, and best practices
* **Core Data Contracts**: [docs/source/data_contracts.rst](docs/source/data_contracts.rst) - Stable AnnData and review-summary conventions shared across workflow stages
* **Workflow Hardening Plan**: [docs/source/workflow_hardening.rst](docs/source/workflow_hardening.rst) - Real-data vertical-slice plan for PBMC, PDAC, and active project validation
* **PBMC Golden Path**: [scripts/run_pbmc_golden_path.py](scripts/run_pbmc_golden_path.py) - Runnable real-data baseline that emits a manifest, final `.h5ad`, and inspection figures
* **Analysis Acceptance Runner**: [scripts/run_analysis_acceptance.py](scripts/run_analysis_acceptance.py) - Runnable Step2 analysis hardening path for clustering review, annotation evidence, consensus labels, and optional malignancy interpretation

For quick examples, see the `examples/` directory.

### Contributing

We welcome contributions from the community! If you'd like to contribute, please check out our [Contributing Guidelines](CONTRIBUTING.md) and the issue tracker.

### License

`scLucid` is licensed under the MIT License.

### How to Cite

If you use `scLucid` in your research before a formal methods paper is available,
please cite the GitHub repository and include the package version used in your
analysis. A manuscript citation will be added once available.
