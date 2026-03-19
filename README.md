# scLucid: A Comprehensive System for Single-Cell Analysis

[![PyPI version](https://badge.fury.io/py/sclucid.svg)](https://badge.fury.io/py/sclucid)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/yelu61/scLucid/actions/workflows/build.yml/badge.svg)](https://github.com/yelu61/scLucid/actions)

**scLucid** is a powerful and flexible Python toolkit for the analysis of single-cell RNA-sequencing data. It is designed to be more than just a wrapper; it's a complete **analysis system** that guides researchers from raw data to deep biological insights.

The toolkit's philosophy is to balance ease-of-use for standard workflows with deep customizability for advanced, exploratory research. It achieves this through a modular architecture, high-level workflow functions, and a unique, biology-aware marker management system.

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
pip install scLucid

# To include additional analysis packages (CellTypist, cosg, etc.)
pip install scLucid[analysis]

# To include all advanced tools (scVelo, rpy2, infercnvpy, etc.)
pip install scLucid[tools]

# To install everything
pip install scLucid[all]
```

You can also install the latest development version directly from GitHub:
```bash
pip install git+[https://github.com/yelu61/scLucid.git](https://github.com/yelu61/scLucid.git)
```

For developers, clone the repository and install in editable mode:
```bash
git clone [https://github.com/yelu61/scLucid.git](https://github.com/yelu61/scLucid.git)
cd scLucid
pip install -e .[all]
```

### Quick Start: A 5-Minute Analysis

Here is a minimal example of a complete workflow.

```python
import scLucid as scl # It is common to create a short alias
from scLucid.preprocess import WorkflowConfig
from scLucid.analysis import ClusteringConfig, AnnotationConfig

# --- 1. Load Data (assuming adata is already loaded) ---
# adata = scl.utils.load_10x_data(...)

# --- 2. Run Quality Control ---
adata_qc = scl.qc.run_standard_qc(adata, species="human")

# --- 3. Run Preprocessing ---
# This single function handles normalization, HVG selection, scaling, PCA, and Harmony integration
prep_config = WorkflowConfig(
    integration={"method": "harmony", "batch_key": "sampleID"}
)
adata_prep = scl.preprocess.run_preprocessing(adata_qc, config=prep_config)

# --- 4. Run Clustering & Annotation ---
# Define a clustering strategy
cluster_config = ClusteringConfig(resolution=0.8, use_rep="X_harmony")
adata_clustered = scl.analysis.cluster_cells(adata_prep, config=cluster_config)

# Define an annotation strategy
anno_config = AnnotationConfig(
    cluster_key="leiden_res0.8",
    marker_species="human",
    final_method="combined" # Use both scoring and enrichment
)
adata_final = scl.analysis.run_annotation(adata_clustered, config=anno_config)

# --- 5. Visualize Final Results ---
# Set publication-ready font style for your target journal
from scLucid import FONT_NATURE, FONT_CELL, FONT_TRADITIONAL
scl.set_figure_params(dpi=300, font_style=FONT_NATURE)  # For Nature/Science

scl.utils.plot_embedding(adata_final, color_by="cell_type")

# Save with embedded fonts for publication
import matplotlib.pyplot as plt
plt.savefig("results.pdf", dpi=600, bbox_inches="tight")
```

### Documentation

For detailed tutorials, how-to guides, and the full API reference:

* **Plugin Development**: [Plugin Development Guide](docs/PLUGIN_DEVELOPMENT_GUIDE.md) - Create custom analysis plugins
* **Naming Conventions**: [Naming Conventions](docs/NAMING_CONVENTIONS.md) - Code style guidelines
* **OpenSpec Specifications**: [openspec/specs/](openspec/specs/) - Technical specifications
* **Full Documentation**: [**Your Documentation URL Here**](docs/) (Under construction)

For quick examples, see the `examples/` directory.

### Contributing

We welcome contributions from the community! If you'd like to contribute, please check out our [Contributing Guidelines](CONTRIBUTING.md) and the issue tracker.

### License

`scLucid` is licensed under the MIT License.

### How to Cite

If you use `scLucid` in your research, please cite our paper:

> *Your Name, et al. (2025). scLucid: A comprehensive and flexible system for single-cell analysis. Journal Name.*