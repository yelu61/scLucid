import scanpy as sc
import scRNA # Your amazing toolkit!

# 1. Load Data
adata = sc.read_10x_mtx(...)
adata.layers['counts'] = adata.X.copy() # IMPORTANT: Backup raw counts

# --- QC Module ---
# 2. Calculate QC metrics
adata = scRNA.qc.calculate_qc_metric(adata, sample_key='batch')

# 3. Identify problematic cells
adata = scRNA.qc.is_low_quality_cell(adata, sample_key='batch')
adata = scRNA.qc.is_doublet(adata, sample_key='batch')

# 4. Filter cells based on the annotations
adata = scRNA.qc.filter_cells(
    adata,
    filter_by_outliers=True,
    filter_by_mt=True,
    filter_by_doublets=True
)

# --- Preprocess Module ---
# 5. Normalize the filtered counts
adata = scRNA.preprocess.normalize_data(
    adata,
    layer='counts', # Use the raw counts as input
    output_layer='log1p_norm'
)

# 6. Annotate Highly Variable Genes
adata = scRNA.preprocess.annotate_hvg(
    adata,
    method='scanpy',
    layer='log1p_norm',
    n_top_genes=2000
)

# 7. Select HVGs (subset the object) and Scale the data
# Note: select_hvg now returns a new, subsetted AnnData object
adata_hvg = scRNA.preprocess.select_hvg(
    adata,
    hvg_key='highly_variable_scanpy',
    subset=True
)

# Now, scale the log-normalized data within the HVG-subsetted object
adata_hvg = scRNA.preprocess.scale_data(
    adata_hvg,
    layer='log1p_norm', # This layer was carried over during subsetting
    output_layer='scaled'
)

# 8. Run PCA
# PCA uses the scaled data in .X by default
sc.pp.pca(adata_hvg, n_comps=50)

# 9. Perform Batch Correction
# This adds the final 'X_harmony' embedding
adata_hvg = scRNA.preprocess.batch_correction(
    adata_hvg,
    batch_key='batch',
    method='harmony',
    use_rep='X_pca' # Use the PCA as input
)

# Your data is now ready for downstream analysis!
print("Preprocessing complete!")
print("Corrected embedding is in adata_hvg.obsm['X_harmony']")