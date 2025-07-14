import scanpy as sc

import scRNA  # Your amazing toolkit!

# 1. Load Data
adata = sc.read_10x_mtx(...)
adata.layers["counts"] = adata.X.copy()  # IMPORTANT: Backup raw counts

# --- QC Module ---
# 2. Calculate QC metrics
adata = scRNA.qc.calculate_qc_metric(adata, sample_key="batch")

# 3. Identify problematic cells
adata = scRNA.qc.is_low_quality_cell(adata, sample_key="batch")
adata = scRNA.qc.is_doublet(adata, sample_key="batch")

# 4. Filter cells based on the annotations
adata = scRNA.qc.filter_cells(
    adata, filter_by_outliers=True, filter_by_mt=True, filter_by_doublets=True
)

# --- Preprocess Module ---
# 5. Normalize the filtered counts
adata = scRNA.preprocess.normalize_data(
    adata,
    layer="counts",  # Use the raw counts as input
    output_layer="log1p_norm",
)

# 6. Annotate Highly Variable Genes
adata = scRNA.preprocess.annotate_hvg(
    adata, method="scanpy", layer="log1p_norm", n_top_genes=2000
)

# 7. Select HVGs (subset the object) and Scale the data
# Note: select_hvg now returns a new, subsetted AnnData object
adata_hvg = scRNA.preprocess.select_hvg(
    adata, hvg_key="highly_variable_scanpy", subset=True
)

# Now, scale the log-normalized data within the HVG-subsetted object
adata_hvg = scRNA.preprocess.scale_data(
    adata_hvg,
    layer="log1p_norm",  # This layer was carried over during subsetting
    output_layer="scaled",
)

# 8. Run PCA
# PCA uses the scaled data in .X by default
sc.pp.pca(adata_hvg, n_comps=50)

# 9. Perform Batch Correction
# This adds the final 'X_harmony' embedding
adata_hvg = scRNA.preprocess.batch_correction(
    adata_hvg,
    batch_key="batch",
    method="harmony",
    use_rep="X_pca",  # Use the PCA as input
)

# Your data is now ready for downstream analysis!
print("Preprocessing complete!")
print("Corrected embedding is in adata_hvg.obsm['X_harmony']")

# --- Analysis Module ---
# Assume 'adata' is the preprocessed object from the previous guide.
# It has PCA and a corrected embedding (e.g., 'X_harmony').

# 1. Find Optimal Clustering Resolution
# This is a powerful first step that uses your marker genes to find biologically meaningful clusters.
adata = scRNA.analysis.find_resolution(
    adata,
    metric="marker_separation",
    marker_config="path/to/your/manager_human.toml",
    use_rep="X_harmony",  # Use the batch-corrected embedding
)
# The result is in adata.obs['leiden_optimal'] or adata.obs['louvain_optimal']

# 2. Find Marker Genes for the new clusters
markers_df = scRNA.analysis.find_markers(adata, groupby="leiden_optimal")

# 3. Filter for high-quality markers
filtered_markers = scRNA.analysis.filter_markers(adata, min_log2fc=1, max_padj=0.05)
# You can now inspect the 'filtered_markers' DataFrame.

# 4. Annotate Clusters
# First, score all cells for cell type identity
adata = scRNA.analysis.score_cell_types(
    adata, marker_config="path/to/your/manager_human.toml"
)

# Now, use those scores to annotate the clusters
adata = scRNA.analysis.annotate_clusters(
    adata,
    cluster_key="leiden_optimal",
    marker_config="path/to/your/manager_human.toml",
    method="max_score",  # Use the max_score method
    key_added="cell_type_annotation",
)

# --- Visualization ---
# 5. Visualize the final annotated UMAP
scRNA.analysis.plot_embedding(
    adata, color_by="cell_type_annotation", title="Final Annotated Cell Types"
)

# 6. Visualize the composition of your clusters
scRNA.analysis.plot_composition(
    adata,
    group_key="batch",  # See how batches are distributed across cell types
    stack_key="cell_type_annotation",
)

# 7. Visualize a heatmap of the top marker genes
scRNA.analysis.plot_marker_heatmap(
    adata, markers_df=filtered_markers, groupby="leiden_optimal", n_genes=5
)
