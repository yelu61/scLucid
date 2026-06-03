"""
Preprocessing Pipeline Example

Demonstrates the light-dependency default path: normalization, HVG selection,
scaling, PCA, neighbors, and UMAP. Batch correction is shown as an explicit
opt-in enhancement.

Before running this script, QC should already have removed clear low-quality
cells. Keep raw counts in ``layers["counts"]``; optional methods such as scran,
Harmony, scVI/scANVI, and BBKNN are enhancements, not defaults.
"""

import scanpy as sc
import scLucid
from scLucid.preprocess import (
    normalize_data,
    find_hvgs,
    select_hvg_sets,
    scale_data,
    batch_correction,
    run_preprocessing,
    NormalizationConfig,
    HVGConfig,
    ScalingConfig,
    IntegrationConfig,
    PreprocessingWorkflowConfig,
)
from scLucid.preprocess.intelligent import run_intelligent_preprocessing
from scLucid.config import set_config

# Configure
set_config(verbosity=1, n_jobs=4)

# Load QC-filtered data
adata = sc.read_h5ad("results/qc_filtered.h5ad")

# ---------------------------------------------------------------------------
# Option A: Standard default path (one-liner)
# ---------------------------------------------------------------------------
# config = PreprocessingWorkflowConfig.default()
# adata = run_preprocessing(adata, config=config, save_dir="results/preprocess")

# ---------------------------------------------------------------------------
# Option B: Intelligent preprocessing with data-driven parameter selection
# ---------------------------------------------------------------------------
# adata, strategy = run_intelligent_preprocessing(
#     adata,
#     batch_key="batch",
#     save_dir="results/preprocess",
# )
# # Review summary is stored in adata.uns["sclucid"]["preprocess"]["intelligent_review_summary"]
# # and written to results/preprocess/preprocess_review_summary.json / .md

# ---------------------------------------------------------------------------
# Option C: Manual step-by-step (shown below)
# ---------------------------------------------------------------------------

# Step 0: Gene filtering
print("Filtering low-detection genes...")
sc.pp.filter_genes(adata, min_cells=3)

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
    method="scanpy",
    n_top_genes=2000,
    flavor="auto",
    span=0.3,
    exclude_gene_types=["mitochondrial", "ribosomal"],
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
    max_value=10
)
adata = scale_data(adata, config=scale_config)

# Step 5: PCA
print("Running PCA...")
sc.tl.pca(adata, n_comps=50)
sc.pl.pca_variance_ratio(adata, log=True, save="_preprocess.pdf")

# Step 6: Optional batch correction (enable only after diagnostics)
RUN_BATCH_CORRECTION = False
if RUN_BATCH_CORRECTION and "batch" in adata.obs.columns and adata.obs["batch"].nunique() > 1:
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
    print("Skipping batch correction in the light default path")
    if "batch" in adata.obs.columns and adata.obs["batch"].nunique() > 1:
        medians = adata.obs.groupby("batch", observed=False)["total_counts"].median()
        if len(medians) > 1 and medians.max() / max(medians.min(), 1) > 2:
            print(
                "Review note: sample-level sequencing depth differs by >2x; "
                "inspect PCA/UMAP before enabling integration."
            )

# Step 7: Neighbors and UMAP
print("Computing neighbors and UMAP...")
sc.pp.neighbors(adata, n_pcs=50, n_neighbors=15)
sc.tl.umap(adata)

# Visualize
sc.pl.umap(adata, color=["batch"] if "batch" in adata.obs else ["phase"],
            save="_preprocess_umap.pdf", dpi=300)

# Save
adata.write("results/preprocessed.h5ad")
print("Preprocessing complete!")
print(f"HVGs: {adata.var['highly_variable'].sum()}")
print(f"PCs: {adata.obsm['X_pca'].shape[1]}")
print(f"UMAP: {adata.obsm['X_umap'].shape}")
