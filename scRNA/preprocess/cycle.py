"""
Utility functions for preprocessing single-cell RNA-seq data.
"""

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import sparse
from typing import Optional, Tuple, List, Dict, Union

def check_adata(
    adata: sc.AnnData,
    required_obs: Optional[List[str]] = None,
    required_var: Optional[List[str]] = None,
    required_layers: Optional[List[str]] = None,
    required_obsm: Optional[List[str]] = None,
) -> bool:
    """
    Check if AnnData object contains required attributes.
    
    Args:
        adata: AnnData object
        required_obs: Required columns in adata.obs
        required_var: Required columns in adata.var
        required_layers: Required layers in adata.layers
        required_obsm: Required matrices in adata.obsm
        
    Returns:
        True if all requirements are met, otherwise raises ValueError
    """
    if required_obs:
        missing_obs = [key for key in required_obs if key not in adata.obs]
        if missing_obs:
            raise ValueError(f"Missing required columns in adata.obs: {', '.join(missing_obs)}")
    
    if required_var:
        missing_var = [key for key in required_var if key not in adata.var]
        if missing_var:
            raise ValueError(f"Missing required columns in adata.var: {', '.join(missing_var)}")
    
    if required_layers:
        missing_layers = [key for key in required_layers if key not in adata.layers]
        if missing_layers:
            raise ValueError(f"Missing required layers in adata.layers: {', '.join(missing_layers)}")
    
    if required_obsm:
        missing_obsm = [key for key in required_obsm if key not in adata.obsm]
        if missing_obsm:
            raise ValueError(f"Missing required matrices in adata.obsm: {', '.join(missing_obsm)}")
    
    return True

def validate_adata(adata):
    """验证AnnData对象的有效性和完整性"""
    import anndata as ad
    
    if not isinstance(adata, ad.AnnData):
        raise TypeError("输入必须是AnnData对象")
    
    if adata.n_obs == 0:
        raise ValueError("AnnData对象不包含任何细胞")
    
    if adata.n_vars == 0:
        raise ValueError("AnnData对象不包含任何基因")
    
    # 检查计数数据的基本特性
    if not sparse.issparse(adata.X) and not isinstance(adata.X, np.ndarray):
        raise TypeError("adata.X必须是稀疏矩阵或numpy数组")
    
    # 检查表达值是否为负
    if isinstance(adata.X, np.ndarray):
        if np.any(adata.X < 0):
            print("警告: 表达矩阵包含负值，可能不是原始计数数据")
    else:
        if sparse.issparse(adata.X) and np.any(adata.X.data < 0):
            print("警告: 表达矩阵包含负值，可能不是原始计数数据")
    
    return True

def plot_preprocessing_summary(
    adata: sc.AnnData,
    save_dir: Optional[str] = None,
    file_name: str = "preprocessing_summary.png",
):
    """
    Generate a summary of preprocessing steps.
    
    Args:
        adata: AnnData object
        save_dir: Directory to save plot
        file_name: File name for saved plot
    """
    # Collect information about the data
    n_cells = adata.n_obs
    n_genes = adata.n_vars
    
    # Check which preprocessing steps have been done
    has_normalized = "log1p_norm" in adata.layers
    has_hvg = "highly_variable" in adata.var or "highly_variable_scanpy" in adata.var
    has_scaled = "scaled" in adata.layers
    has_integrated = any(key.startswith("X_") and key != "X_pca" for key in adata.obsm)
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("Preprocessing Summary", fontsize=16)
    
    # Plot 1: Cell count stats
    if "sampleID" in adata.obs:
        sample_counts = adata.obs["sampleID"].value_counts()
        sample_counts.plot(kind="bar", ax=axes[0, 0])
        axes[0, 0].set_title("Cells per Sample")
        axes[0, 0].set_ylabel("Number of Cells")
        axes[0, 0].tick_params(axis='x', rotation=45)
    else:
        axes[0, 0].text(0.5, 0.5, f"Total Cells: {n_cells}\nTotal Genes: {n_genes}", 
                     ha='center', va='center', fontsize=12)
        axes[0, 0].set_title("Dataset Size")
        axes[0, 0].axis('off')
    
    # Plot 2: Gene expression distribution
    if has_normalized:
        if "log1p_norm" in adata.layers:
            gene_means = np.mean(adata.layers["log1p_norm"], axis=0)
            axes[0, 1].hist(gene_means, bins=50)
            axes[0, 1].set_title("Gene Expression Distribution (Normalized)")
            axes[0, 1].set_xlabel("Mean Expression")
            axes[0, 1].set_ylabel("Frequency")
    else:
        axes[0, 1].text(0.5, 0.5, "Normalization not performed", 
                     ha='center', va='center', fontsize=12)
        axes[0, 1].set_title("Gene Expression")
        axes[0, 1].axis('off')
    
    # Plot 3: HVG stats
    if has_hvg:
        hvg_key = "highly_variable" if "highly_variable" in adata.var else "highly_variable_scanpy"
        n_hvgs = sum(adata.var[hvg_key])
        
        if "dispersions" in adata.var:
            # Create scatter plot of dispersion vs mean
            axes[1, 0].scatter(
                adata.var["means"],
                adata.var["dispersions"],
                c=adata.var[hvg_key].map({True: "red", False: "gray"}),
                alpha=0.5
            )
            axes[1, 0].set_title(f"Highly Variable Genes ({n_hvgs})")
            axes[1, 0].set_xlabel("Mean")
            axes[1, 0].set_ylabel("Dispersion")
            axes[1, 0].set_xscale("log")
            axes[1, 0].set_yscale("log")
        else:
            # Just show number of HVGs
            axes[1, 0].text(0.5, 0.5, f"Number of HVGs: {n_hvgs}", 
                         ha='center', va='center', fontsize=12)
            axes[1, 0].set_title("Highly Variable Genes")
            axes[1, 0].axis('off')
    else:
        axes[1, 0].text(0.5, 0.5, "HVG selection not performed", 
                     ha='center', va='center', fontsize=12)
        axes[1, 0].set_title("Highly Variable Genes")
        axes[1, 0].axis('off')
    
    # Plot 4: Preprocessing steps
    steps = ["Normalized", "HVG Selection", "Scaled", "Integrated"]
    completed = [has_normalized, has_hvg, has_scaled, has_integrated]
    
    axes[1, 1].bar(steps, [int(c) for c in completed], color=["green" if c else "red" for c in completed])
    axes[1, 1].set_title("Preprocessing Steps Completed")
    axes[1, 1].set_ylim(0, 1.5)
    axes[1, 1].set_yticks([0, 1])
    axes[1, 1].set_yticklabels(["No", "Yes"])
    axes[1, 1].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    
    if save_dir:
        import os
        os.makedirs(save_dir, exist_ok=True)
        plt.savefig(os.path.join(save_dir, file_name), dpi=300)
    
    plt.show()

import logging

def setup_logger(name, level=logging.INFO, log_file=None):
    """设置日志记录器"""
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # 控制台处理程序
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    logger.addHandler(console)
    
    # 文件处理程序(如果指定)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger
def progress_monitor(iterable, desc="Processing", total=None):
    """使用tqdm创建进度条"""
    try:
        from tqdm.notebook import tqdm
        return tqdm(iterable, desc=desc, total=total)
    except ImportError:
        try:
            from tqdm import tqdm
            return tqdm(iterable, desc=desc, total=total)
        except ImportError:
            print("建议安装tqdm包以显示进度条")
            return iterable

def memory_efficient_operation(func):
    """装饰器用于内存高效操作"""
    import gc
    import functools
    
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # 执行操作前收集垃圾
        gc.collect()
        
        # 调用原始函数
        result = func(*args, **kwargs)
        
        # 操作后再次收集垃圾
        gc.collect()
        
        return result
    
    return wrapper