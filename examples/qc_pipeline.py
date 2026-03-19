"""
Quality Control Pipeline Example

Demonstrates comprehensive QC workflow with doublet detection and filtering.
"""

import scanpy as sc
import scLucid
from scLucid.qc import (
    calculate_qc_metric,
    predict_doublets,
    filter_cells,
    QCThresholds,
    DoubletConfig
)
from scLucid.config import set_config

# Configure
set_config(verbosity=1, n_jobs=4)

# Load data
adata = sc.read_h5ad("data/raw_counts.h5ad")
adata.layers["counts"] = adata.X.copy()

# Step 1: Calculate QC metrics
print("Calculating QC metrics...")
adata = calculate_qc_metric(
    adata,
    species="human",
    config=None,  # Use defaults
    save_dir="results/qc/metrics"
)

# Step 2: Detect doublets
print("Detecting doublets...")
doublet_config = DoubletConfig(
    method="scrublet",
    expected_doublet_rate=0.06
)
adata = predict_doublets(adata, config=doublet_config)

# Step 3: Define thresholds
thresholds = QCThresholds(
    min_genes=200,
    max_genes=5000,
    pc_mt=20.0,
    nmads=5.0
)

# Step 4: Filter cells
print("Filtering cells...")
adata_filtered = filter_cells(adata, thresholds)

# Step 5: Visualize
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Before vs After
sc.pl.violin(adata, ['n_genes_by_counts', 'pct_counts_mt'], ax=axes[0], show=False)
axes[0].set_title('Before QC')

sc.pl.violin(adata_filtered, ['n_genes_by_counts', 'pct_counts_mt'], ax=axes[1], show=False)
axes[1].set_title('After QC')

sc.pl.violin(adata, ['doublet_score'], ax=axes[2], show=False)
axes[2].set_title('Doublet Scores')

plt.tight_layout()
plt.savefig('results/qc/qc_comparison.pdf', dpi=300)

# Save
adata_filtered.write("results/qc_filtered.h5ad")
print(f"✅ QC complete!")
print(f"Original cells: {adata.n_obs}")
print(f"Filtered cells: {adata_filtered.n_obs}")
print(f"Removed: {(1 - adata_filtered.n_obs/adata.n_obs)*100:.1f}%")
