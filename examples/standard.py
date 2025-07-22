import scanpy as sc

import scRNA  # Your amazing toolkit!

# 1. Load Data
adata = sc.read_10x_mtx(...)
adata.layers["counts"] = adata.X.copy()  # IMPORTANT: Backup raw counts

# --- QC Module ---
# 2. Calculate QC metrics
adata = scRNA.qc.calculate_qc_metric(adata, sample_key="batch")

# 3. Identify problematic cells
adata = scRNA.qc.is_doublet(adata, sample_key="batch")
adata = scRNA.qc.is_low_quality_cell(adata, sample_key="batch")

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

# 7. Regress out cell cycle (optional, good to do before scaling)
adata = scRNA.preprocess.score_cell_cycle(adata)
adata = scRNA.preprocess.regress_out(
    adata,
    keys=["S_score", "G2M_score"],
    layer="log1p_norm",
    output_layer="log1p_norm_regressed" # Store in a new layer
)

# 8. Select HVGs (subset the object) and Scale the data
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

# 9. Run PCA
# PCA uses the scaled data in .X by default
sc.pp.pca(adata_hvg, n_comps=50)

# 10. Perform Batch Correction
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

# 假设 adata_hvg 是预处理好的 AnnData 对象
# ...

# --- 阶段一: 聚类 ---
# 1. 优化参数 (假如之前已完成，这里直接使用)
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30, use_rep="X_harmony")
sc.tl.umap(adata_hvg)

# 2. 寻找最优分辨率并聚类
# 使用marker分离度作为指标，需要提供marker配置文件
adata_hvg = scRNA.analysis.find_resolution(
    adata_hvg,
    metric="marker_separation",
    marker_config="./manager_human.toml"
)

# 3. 可视化聚类结果
scRNA.analysis.plot_embedding(adata_hvg, color_by="leiden_optimal", title="Optimal Leiden Clustering")

# --- 阶段二: 细胞注释 ---
# 4. 基于打分进行初步注释
adata_hvg = scRNA.analysis.score_cell_types(adata_hvg, marker_config="./manager_human.toml")
adata_hvg = scRNA.analysis.annotate_clusters(
    adata_hvg,
    cluster_key="leiden_optimal",
    marker_config="./manager_human.toml",
    method="max_score"
)

# 5. 通过DE分析寻找数据驱动的marker
markers_df = scRNA.analysis.find_markers(adata_hvg, groupby="leiden_optimal")
filtered_markers = scRNA.analysis.filter_markers(markers_df, min_log2fc=1.0)

# 6. 人工校验
# 查看自动注释结果
scRNA.analysis.plot_embedding(adata_hvg, color_by="leiden_optimal_annotated", title="Automated Annotation (Max Score)")

# 验证B细胞(假设是cluster '3')的marker
scRNA.analysis.plot_embedding(adata_hvg, color=["MS4A1", "CD79A"], title="B Cell Markers")

# 查看各cluster的top marker热图
scRNA.analysis.plot_marker_heatmap(adata_hvg, markers_df=filtered_markers, groupby="leiden_optimal", n_genes=5)

# 假设经过校验，您创建了一个最终的手动注释列 adata_hvg.obs['cell_type_manual']

# --- 阶段三: 下游分析 ---
# 7. 查看最终注释的细胞类型组成
scRNA.analysis.plot_composition(adata_hvg, group_key="sample", stack_key="cell_type_manual")