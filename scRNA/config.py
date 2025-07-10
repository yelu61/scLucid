# preprocess/config.py
"""单细胞RNA-seq数据预处理配置模块"""

import os
from typing import Dict, Any, Optional

# 默认参数
DEFAULT_PARAMS = {
    # 通用参数
    "batch_key": "sampleID",     # 批次/样本标识符
    "random_seed": 42,           # 随机种子
    
    # 层名称
    "layer_raw": "counts",       # 原始计数层
    "layer_norm": "log1p_norm",  # 标准化层
    "layer_scale": "scaled",     # 缩放层
    
    # 标准化参数
    "target_sum": 1e4,           # 标准化目标总和
    
    # 高变基因参数
    "n_top_genes": 2000,         # HVG数量
    "min_mean": 0.0125,          # 最小平均表达
    "max_mean": 3,               # 最大平均表达
    "min_disp": 0.5,             # 最小离散度
    
    # 降维参数
    "n_pcs": 50,                 # PCA维度
    "n_neighbors": 15,           # KNN邻居数
    
    # 聚类参数
    "resolution": 0.8,           # 聚类分辨率
    
    # 质量控制参数
    "min_genes": 200,            # 最小基因数
    "max_genes": 5000,           # 最大基因数
    "min_cells": 3,              # 最小细胞数
    "max_mt_percent": 20,        # 最大线粒体百分比
}

def get_param(name: str, user_params: Optional[Dict[str, Any]] = None) -> Any:
    """获取参数值，优先使用用户指定的值"""
    if user_params and name in user_params:
        return user_params[name]
    if name in DEFAULT_PARAMS:
        return DEFAULT_PARAMS[name]
    raise ValueError(f"未知参数: {name}")

def load_config(config_file: str) -> Dict[str, Any]:
    """从配置文件加载参数"""
    import json
    
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"配置文件不存在: {config_file}")
    
    with open(config_file, 'r') as f:
        user_params = json.load(f)
    
    return user_params