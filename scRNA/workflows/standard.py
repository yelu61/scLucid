"""
Workflow module for single-cell RNA-seq data analysis.
"""

import scanpy as sc
import anndata as ad
from typing import Optional, List, Union, Dict
from . import qc, norm, hvg, integrate
from .config import load_config
# workflow.py 修改
def standard_analysis_pipeline(
    adata,
    sample_key="sampleID",
    n_top_genes=2000,
    normalization_method="standard",
    batch_correction_method="harmony",
    clustering_resolution=0.8,
    output_dir=None,
    **kwargs
):
    """
    执行标准单细胞RNA-seq分析流程
    
    完整流程: QC -> 标准化 -> 批次校正 -> HVG选择 -> PCA -> 聚类 -> 可视化
    """
    import scanpy as sc
    import os
    from . import qc, preprocess
    
    # 设置输出目录
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    # 步骤1: 质量控制
    print("步骤1: 执行质量控制")
    adata = qc.calculate_qc_metric(adata, sample_key=sample_key)
    adata = qc.is_low_quality_cell(adata, sample_key=sample_key)
    adata = qc.is_doublet(adata, sample_key=sample_key)
    adata_filtered = qc.filter_low_quality_cells(
        adata, filter_outliers=True, filter_mt=True, filter_doublets=True
    )
    
    if output_dir:
        adata_filtered.write(os.path.join(output_dir, "01_filtered_data.h5ad"))
    
    # 步骤2: 标准化和对数变换
    print("步骤2: 数据标准化")
    adata_filtered = preprocess.normalize_data(
        adata_filtered, method=normalization_method, plot=True,
        save_dir=output_dir
    )
    
    # 步骤3: 细胞周期评分
    print("步骤3: 细胞周期评分")
    adata_filtered = preprocess.score_cell_cycle(
        adata_filtered, species="auto", save_dir=output_dir
    )
    
    # 步骤4: 高变基因选择
    print("步骤4: 识别高变基因")
    adata_filtered = preprocess.annotate_hvg(
        adata_filtered, method="scanpy", n_top_genes_scanpy=n_top_genes
    )
    adata_filtered = preprocess.select_hvg(
        adata_filtered, subset=True
    )
    
    if output_dir:
        adata_filtered.write(os.path.join(output_dir, "02_normalized_hvg.h5ad"))
    
    # 步骤5: 缩放数据
    print("步骤5: 数据缩放")
    adata_filtered = preprocess.scale_data(
        adata_filtered, layer="log1p_norm", vars_to_regress=["S_score", "G2M_score"]
    )
    
    # 步骤6: 降维
    print("步骤6: 主成分分析")
    sc.pp.pca(adata_filtered, n_comps=50)
    
    # 步骤7: 批次校正(如果需要)
    if len(adata_filtered.obs[sample_key].unique()) > 1:
        print("步骤7: 批次校正")
        adata_filtered = preprocess.batch_correction(
            adata_filtered, batch_key=sample_key, method=batch_correction_method,
            plot=True, save_dir=output_dir
        )
    
    # 步骤8: 构建KNN图和UMAP嵌入
    print("步骤8: 构建KNN图和UMAP嵌入")
    sc.pp.neighbors(adata_filtered)
    sc.tl.umap(adata_filtered)
    
    # 步骤9: 聚类
    print("步骤9: 执行聚类")
    sc.tl.leiden(adata_filtered, resolution=clustering_resolution)
    
    # 步骤10: 生成结果可视化
    print("步骤10: 生成结果可视化")
    if output_dir:
        sc.pl.umap(
            adata_filtered, 
            color=["leiden", "phase", sample_key], 
            save=os.path.join(output_dir, "umap_clusters.pdf")
        )
        
        # 保存最终结果
        adata_filtered.write(os.path.join(output_dir, "final_analysis.h5ad"))
    
    # 返回分析结果
    return adata_filtered

# 其他自定义工作流...