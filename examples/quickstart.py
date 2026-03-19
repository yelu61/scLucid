"""
Quick Start Example - scLucid in 5 Minutes

This script demonstrates the complete scLucid pipeline:
QC → Preprocessing → Clustering → Annotation
"""

import scanpy as sc
import scLucid
from scLucid.qc import run_standard_qc, QCWorkflowConfig
from scLucid.preprocess import run_preprocessing, WorkflowConfig
from scLucid.analysis import cluster_cells, annotate_clusters, ClusteringConfig
from scLucid.utils.manager import get_marker_manager

# 1. Load data
print("Loading data...")
adata = sc.read_h5ad("data/pbmc_raw.h5ad")
adata.layers["counts"] = adata.X.copy()  # Backup raw counts

# 2. Quality Control
print("Running QC...")

# Option A: Simple dictionary-based configuration (recommended for quick start)
qc_config = QCWorkflowConfig.from_simple_dict({
    "thresholds_min_genes": 200,
    "thresholds_pc_mt": 20.0,
    "doublet_method": "scrublet",
    "species": "human",
    "save_dir": "results/qc"
})

# Option B: Ultra-quick configuration for standard analyses
# qc_config = QCWorkflowConfig.quick(min_genes=200, pc_mt=20.0, species="human")

# Option C: Full nested configuration (for advanced users)
# from scLucid.qc import QCThresholds, DoubletConfig, MarkingConfig
# qc_config = QCWorkflowConfig(
#     marking_config=MarkingConfig(thresholds=QCThresholds(min_genes=200, pc_mt=20.0)),
#     doublet_config=DoubletConfig(method="scrublet"),
#     save_dir="results/qc"
# )

adata = run_standard_qc(adata, config=qc_config)
print(f"Cells after QC: {adata.n_obs}")

# 3. Preprocessing
print("Running preprocessing...")

# Option A: Simple dictionary-based configuration (recommended for quick start)
pp_config = WorkflowConfig.from_simple_dict({
    "normalization_target_sum": 1e4,
    "hvg_n_top_genes": 2000,
    "graph_n_pcs": 50,
    "save_dir": "results/preprocess"
})

# Option B: Ultra-quick configuration for standard analyses
# pp_config = WorkflowConfig.quick(n_top_genes=2000, run_regression=False)

# Option C: Full nested configuration (for advanced users)
# from scLucid.preprocess import NormalizationConfig, HVGConfig
# pp_config = WorkflowConfig(
#     normalization=NormalizationConfig(target_sum=1e4),
#     hvg=HVGConfig(n_top_genes=2000),
#     save_dir="results/preprocess"
# )

# Example: Skip regression and use all genes
# adata = run_preprocessing(
#     adata,
#     config=pp_config,
#     skip_steps=["regression", "subset_hvg"]
# )

adata = run_preprocessing(
    adata,
    config=pp_config,
    show_progress=True,  # Show progress bar for each step
    # error_recovery=True,  # Enable error recovery (optional)
    # recovery_save_dir="./recovery",  # Save partial results on error
)
print("Preprocessing complete!")

# 4. Clustering
print("Clustering cells...")
cluster_config = ClusteringConfig(method="leiden", resolution=1.0)
adata = cluster_cells(adata, config=cluster_config)

# 5. Annotation
print("Annotating cell types...")
mgr = get_marker_manager(species="human", tissue="Blood")
adata = annotate_clusters(adata, marker_config=mgr)

# 6. Save and visualize
print("Saving results...")
adata.write("results/pbmc_analyzed.h5ad")
sc.pl.umap(adata, color=["leiden", "cell_type"], save="_final.pdf")

print("✅ Analysis complete!")
print(f"Results saved to: results/")
