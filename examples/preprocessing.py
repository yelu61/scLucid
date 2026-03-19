"""
Preprocessing Pipeline Example

Demonstrates normalization, HVG selection, scaling, PCA, and batch correction.
"""

import scanpy as sc
import scLucid
from scLucid.preprocess import (
    normalize_data,
    find_hvgs,
    select_hvg_sets,
    scale_data,
    batch_correction,
    NormalizationConfig,
    HVGConfig,
    ScalingConfig,
    IntegrationConfig
)
from scLucid.config import set_config

# Configure
set_config(verbosity=1, n_jobs=4)

# Load QC-filtered data
adata = sc.read_h5ad("results/qc_filtered.h5ad")

# Step 1: Normalization
print("Normalizing data...")
norm_config = NormalizationConfig(
    target_sum=1e4,
    output_layer="normalized"
)
adata = normalize_data(adata, config=norm_config, save_dir="results/preprocess")

# Step 2: Find HVGs
print("Finding highly variable genes...")
hvg_config = HVGConfig(
    method="seurat",
    n_top_genes=2000,
    min_mean=0.0125,
    max_mean=3,
    span=0.3
)
adata = find_hvgs(
    adata,
    config=hvg_config,
    input_layer="normalized",
    force=True,
    save_dir="results/preprocess/hvg"
)

# Step 3: Subset to HVGs
print("Subsetting to HVGs...")
adata = select_hvg_sets(
    adata,
    hvg_keys=["highly_variable"],
    mode="direct",
    subset=True,
    keep_raw=False
)

# Step 4: Scale data
print("Scaling data...")
scale_config = ScalingConfig(
    scale_zero_center=True,
    max_value=10
)
adata = scale_data(adata, config=scale_config)

# Step 5: PCA
print("Running PCA...")
sc.tl.pca(adata, n_comps=50)
sc.pl.pca_variance_ratio(adata, log=True, save="_preprocess.pdf")

# Step 6: Batch correction (if needed)
if "batch" in adata.obs.columns and adata.obs["batch"].nunique() > 1:
    print("Running batch correction...")
    integration_config = IntegrationConfig(
        method="harmony",
        batch_key="batch"
    )
    adata = batch_correction(
        adata,
        config=integration_config,
        save_dir="results/preprocess/integration"
    )
else:
    print("No batch correction needed (single sample)")

# Step 7: Neighbors and UMAP
print("Computing neighbors and UMAP...")
sc.pp.neighbors(adata, n_pcs=50, n_neighbors=15)
sc.tl.umap(adata)

# Visualize
sc.pl.umap(adata, color=["batch"] if "batch" in adata.obs else ["phase"],
            save="_preprocess_umap.pdf", dpi=300)

# Save
adata.write("results/preprocessed.h5ad")
print("✅ Preprocessing complete!")
print(f"HVGs: {adata.var['highly_variable'].sum()}")
print(f"PCs: {adata.obsm['X_pca'].shape[1]}")
print(f"UMAP: {adata.obsm['X_umap'].shape}")
