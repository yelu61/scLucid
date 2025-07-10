# 添加新模块 annotation.py
"""
细胞类型注释模块，用于自动化细胞类型识别
"""

import scanpy as sc
import pandas as pd
import numpy as np
from typing import Dict, List, Optional

def score_gene_groups(
    adata: sc.AnnData,
    markers: Dict[str, List[str]],
    ctrl_size: int = 50,
    score_name: str = "cell_type_scores"
) -> sc.AnnData:
    """
    根据已知标记基因对细胞类型进行评分
    
    Args:
        adata: AnnData对象
        markers: 细胞类型标记基因字典，{细胞类型: [标记基因列表]}
        ctrl_size: 每个标记基因集的控制基因数量
        score_name: 结果评分存储前缀名称
        
    Returns:
        添加了细胞类型评分的AnnData对象
    """
    # 检查标记基因是否存在于数据中
    for cell_type, genes in markers.items():
        genes_found = [gene for gene in genes if gene in adata.var_names]
        if len(genes_found) == 0:
            print(f"警告: 未找到细胞类型 '{cell_type}' 的任何标记基因")
        elif len(genes_found) < len(genes):
            print(f"警告: 细胞类型 '{cell_type}' 的部分标记基因未找到, 仅使用 {len(genes_found)}/{len(genes)} 个基因")
            markers[cell_type] = genes_found
    
    # 为每种细胞类型计算评分
    for cell_type, genes in markers.items():
        if len(genes) > 0:
            score_key = f"{score_name}_{cell_type}"
            sc.tl.score_genes(adata, genes, ctrl_size=ctrl_size, score_name=score_key)
    
    return adata

def classify_cells(
    adata: sc.AnnData,
    score_name: str = "cell_type_scores",
    min_score: float = 0.1,
    category_name: str = "predicted_cell_type"
) -> sc.AnnData:
    """
    基于得分将细胞分类为最可能的细胞类型
    
    Args:
        adata: 包含细胞类型评分的AnnData对象
        score_name: 评分的前缀名称
        min_score: 分配细胞类型的最小评分阈值
        category_name: 分类结果的列名
        
    Returns:
        添加了预测细胞类型的AnnData对象
    """
    # 获取所有评分列
    score_cols = [col for col in adata.obs.columns if col.startswith(score_name)]
    
    if len(score_cols) == 0:
        raise ValueError(f"未找到以'{score_name}'开头的评分列")
    
    # 提取实际细胞类型名称
    cell_types = [col.replace(f"{score_name}_", "") for col in score_cols]
    
    # 找出每个细胞的最高评分及其对应的细胞类型
    cell_type_scores = adata.obs[score_cols]
    max_score_idx = np.argmax(cell_type_scores.values, axis=1)
    max_scores = np.max(cell_type_scores.values, axis=1)
    
    # 分配细胞类型，如果最高分低于阈值则标记为"未知"
    predicted_types = np.array(["未知"] * adata.n_obs, dtype=object)
    for i, (idx, score) in enumerate(zip(max_score_idx, max_scores)):
        if score >= min_score:
            predicted_types[i] = cell_types[idx]
    
    # 添加结果到adata.obs
    adata.obs[category_name] = pd.Categorical(predicted_types)
    
    # 打印分类统计信息
    type_counts = adata.obs[category_name].value_counts()
    print(f"细胞类型分类结果:")
    for cell_type, count in type_counts.items():
        print(f"  {cell_type}: {count} 细胞 ({count/adata.n_obs:.2%})")
    
    return adata

def annotate_clusters(
    adata: sc.AnnData,
    cluster_key: str = "leiden",
    cell_type_key: str = "predicted_cell_type",
    output_key: str = "cluster_annotation",
    min_fraction: float = 0.5
) -> sc.AnnData:
    """
    基于预测的细胞类型注释聚类
    
    Args:
        adata: AnnData对象
        cluster_key: 聚类标识的obs列名
        cell_type_key: 细胞类型预测结果的obs列名
        output_key: 输出聚类注释的obs列名
        min_fraction: 将聚类标记为特定类型所需的最小细胞比例
        
    Returns:
        添加了聚类注释的AnnData对象
    """
    if cluster_key not in adata.obs:
        raise ValueError(f"聚类键'{cluster_key}'不在adata.obs中")
    
    if cell_type_key not in adata.obs:
        raise ValueError(f"细胞类型键'{cell_type_key}'不在adata.obs中")
    
    # 获取每个聚类中的细胞类型分布
    cluster_annotations = {}
    
    for cluster in adata.obs[cluster_key].unique():
        cluster_mask = adata.obs[cluster_key] == cluster
        cell_types = adata.obs.loc[cluster_mask, cell_type_key]
        
        # 计算细胞类型计数和比例
        type_counts = cell_types.value_counts()
        type_fractions = type_counts / type_counts.sum()
        
        # 找出最多的细胞类型
        most_common_type = type_counts.idxmax()
        fraction = type_fractions[most_common_type]
        
        # 如果比例足够高，则注释为该类型
        if fraction >= min_fraction:
            annotation = f"{most_common_type} ({fraction:.1%})"
        else:
            annotation = f"Mixed (top: {most_common_type}, {fraction:.1%})"
        
        cluster_annotations[cluster] = annotation
    
    # 将注释添加到adata.obs
    adata.obs[output_key] = adata.obs[cluster_key].map(cluster_annotations)
    
    return adata