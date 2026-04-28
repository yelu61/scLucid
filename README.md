# scLucid: A Comprehensive System for Single-Cell Analysis

[![PyPI version](https://badge.fury.io/py/sclucid.svg)](https://badge.fury.io/py/sclucid)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/yelu61/scLucid/actions/workflows/build.yml/badge.svg)](https://github.com/yelu61/scLucid/actions)

**scLucid** is a powerful and flexible Python toolkit for the analysis of single-cell RNA-sequencing data. It is designed to be more than just a wrapper; it's a complete **analysis system** that guides researchers from raw data to deep biological insights.

The toolkit's philosophy is to balance ease-of-use for standard workflows with deep customizability for advanced, exploratory research. It achieves this through a modular architecture, high-level workflow functions, and a unique, biology-aware marker management system.

### Project Status

scLucid is in active development. The core package already has stable workflow
entrypoints, AnnData contracts, review summaries, and lightweight CI gates. The
next development stage is real-data workflow hardening: running the same
traceable path on PBMC and PDAC datasets, then using active research projects as
acceptance tests for biological plausibility and usability.

### Key Features

* **🧪 End-to-End Workflows**: High-level functions like `run_standard_qc` and `run_preprocessing` to go from raw data to a clustered UMAP with just a few lines of code.
* **🧠 Intelligent QC**: Advanced doublet detection using a combination of `scrublet` and a novel, marker co-expression heuristic.
* **🧬 Biology-Aware Analysis**: A powerful `MarkerManager` system that integrates expert-curated gene sets for annotation, QC, and clustering evaluation.
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
* **🔌 Extensible Plugin Architecture**: Abstract base classes and factory pattern allow you to create custom analysis plugins without modifying core code. See [Plugin Development Guide](docs/PLUGIN_DEVELOPMENT_GUIDE.md) for details.

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

### Quick Start: A 5-Minute Analysis

Here is a minimal example of a complete workflow.

```python
import scanpy as sc
import scLucid as scl

# --- 1. Load Data ---
adata = sc.read_h5ad("data/pbmc3k.h5ad")
adata.layers["counts"] = adata.X.copy()

# --- 2. Run the supported core workflow ---
adata_final = scl.run_pipeline(
    adata,
    stages=["qc", "preprocess", "analysis"],
    dataset_type="pbmc_or_blood",
    species="human",
    show_progress=True,
)

# --- 3. Visualize Final Results ---
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
* **Naming Conventions**: [Naming Conventions](docs/NAMING_CONVENTIONS.md) - Code style guidelines
* **OpenSpec Specifications**: [openspec/specs/](openspec/specs/) - Technical specifications
* **Local Documentation Source**: [docs/source/](docs/source/) - Sphinx documentation sources for installation, quickstart, API references, and best practices
* **Core Data Contracts**: [docs/source/data_contracts.rst](docs/source/data_contracts.rst) - Stable AnnData and review-summary conventions shared across workflow stages
* **Workflow Hardening Plan**: [docs/source/workflow_hardening.rst](docs/source/workflow_hardening.rst) - Real-data vertical-slice plan for PBMC, PDAC, and active project validation

For quick examples, see the `examples/` directory.

### Contributing

We welcome contributions from the community! If you'd like to contribute, please check out our [Contributing Guidelines](CONTRIBUTING.md) and the issue tracker.

### License

`scLucid` is licensed under the MIT License.

### How to Cite

If you use `scLucid` in your research before a formal methods paper is available,
please cite the GitHub repository and include the package version used in your
analysis. A manuscript citation will be added once available.
