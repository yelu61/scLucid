"""
数据加载工具函数

提供便捷的数据加载接口，用于测试、示例和文档。
支持物种差异处理（人源 vs 鼠源）。
"""

from pathlib import Path
from typing import Dict

import scanpy as sc
from anndata import AnnData

# 数据目录
DATA_DIR = Path(__file__).parent.parent / "data"


def load_pbmc3k() -> AnnData:
    """
    加载PBMC3K数据集（正常组织，人源）

    Returns:
    -------
    adata : AnnData
        PBMC数据，已计算QC指标
    """
    from ..qc import calculate_qc_metrics

    adata = sc.read_h5ad(DATA_DIR / "pbmc3k" / "pbmc3k_raw.h5ad")
    adata = calculate_qc_metrics(adata)

    # 添加元数据
    adata.obs["sample_id"] = "pbmc"
    adata.obs["tissue_type"] = "normal"
    adata.obs["species"] = "human"
    adata.obs["batch"] = "pbmc_batch"  # 单批次

    return adata


def load_luad() -> AnnData:
    """
    加载LUAD数据集（肺腺癌，人源，肿瘤组织）

    Returns:
    -------
    adata : AnnData
        LUAD数据，已计算QC指标
    """
    from ..qc import calculate_qc_metrics

    # 尝试不同的文件名
    luad_dir = DATA_DIR / "human_LUAD_GSE131907"
    luad_files = list(luad_dir.glob("*.h5ad"))

    if not luad_files:
        raise FileNotFoundError(f"LUAD数据集未找到: {luad_dir}")

    adata = sc.read_h5ad(luad_files[0])
    adata = calculate_qc_metrics(adata)

    # 添加元数据
    adata.obs["sample_id"] = "luad"
    adata.obs["tissue_type"] = "lung_tumor"
    adata.obs["species"] = "human"

    # 检查批次信息
    if "batch" not in adata.obs.columns:
        adata.obs["batch"] = "luad_batch"  # 如果没有批次信息

    return adata


def load_melanoma() -> AnnData:
    """
    加载黑色素瘤数据集（鼠源，多批次，肿瘤组织）

    Returns:
    -------
    adata : AnnData
        黑色素瘤数据，已计算QC指标
    """
    from ..qc import calculate_qc_metrics

    # 尝试不同的文件名
    melanoma_dir = DATA_DIR / "mouse_melanoma_GSE119352"
    melanoma_files = list(melanoma_dir.glob("*.h5ad"))

    if not melanoma_files:
        raise FileNotFoundError(f"黑色素瘤数据集未找到: {melanoma_dir}")

    adata = sc.read_h5ad(melanoma_files[0])
    adata = calculate_qc_metrics(adata)

    # 添加元数据
    adata.obs["sample_id"] = "melanoma"
    adata.obs["tissue_type"] = "melanoma"
    adata.obs["species"] = "mouse"  # 鼠源

    # 多批次（如果没有批次信息，创建模拟批次）
    if "batch" not in adata.obs.columns:
        n_batches = 3
        adata.obs["batch"] = [f"melanoma_batch_{i % n_batches}" for i in range(adata.n_obs)]

    return adata


def load_all_datasets() -> Dict[str, AnnData]:
    """
    加载所有数据集

    Returns:
    -------
    datasets : dict
        {dataset_name: AnnData}

    Note:
    ----
    数据集特征：
    - PBMC: 人源，正常组织，单批次
    - LUAD: 人源，肺腺癌，肿瘤组织
    - 黑色素瘤: 鼠源，多批次，肿瘤组织

    Example:
    -------
    >>> from scLucid.utils.data_loader import load_all_datasets
    >>>
    >>> # 加载所有数据集
    >>> datasets = load_all_datasets()
    >>>
    >>> # 查看每个数据集的特征
    >>> for name, adata in datasets.items():
    ...     print(f"{name}:")
    ...     print(f"  物种: {adata.obs['species'].iloc[0]}")
    ...     print(f"  组织类型: {adata.obs['tissue_type'].iloc[0]}")
    ...     print(f"  细胞数: {adata.n_obs}")
    ...     print(f"  批次数: {adata.obs['batch'].nunique()}")
    """
    datasets = {}

    try:
        datasets["PBMC"] = load_pbmc3k()
    except Exception as e:
        print(f"⚠ PBMC加载失败: {e}")

    try:
        datasets["LUAD"] = load_luad()
    except Exception as e:
        print(f"⚠ LUAD加载失败: {e}")

    try:
        datasets["Melanoma"] = load_melanoma()
    except Exception as e:
        print(f"⚠ 黑色素瘤加载失败: {e}")

    if len(datasets) == 0:
        raise ValueError("没有成功加载任何数据集")

    print(f"✓ 成功加载 {len(datasets)} 个数据集")

    return datasets


def get_dataset_info(adata: AnnData) -> Dict[str, any]:
    """
    获取数据集的基本信息

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix

    Returns:
    -------
    info : dict
        数据集信息字典

    Example:
    -------
    >>> from scLucid.utils.data_loader import load_pbmc3k, get_dataset_info
    >>>
    >>> adata = load_pbmc3k()
    >>> info = get_dataset_info(adata)
    >>>
    >>> print(f"数据集: {info['sample_id']}")
    >>> print(f"物种: {info['species']}")
    >>> print(f"组织类型: {info['tissue_type']}")
    >>> print(f"细胞数: {info['n_cells']}")
    >>> print(f"基因数: {info['n_genes']}")
    >>> print(f"批次数: {info['n_batches']}")
    """
    info = {
        "sample_id": adata.obs["sample_id"].iloc[0] if "sample_id" in adata.obs else "unknown",
        "species": adata.obs["species"].iloc[0] if "species" in adata.obs else "unknown",
        "tissue_type": (
            adata.obs["tissue_type"].iloc[0] if "tissue_type" in adata.obs else "unknown"
        ),
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "n_batches": adata.obs["batch"].nunique() if "batch" in adata.obs else 1,
    }

    # QC指标
    if "n_genes" in adata.obs.columns:
        info["median_n_genes"] = adata.obs["n_genes"].median()
        info["mean_n_genes"] = adata.obs["n_genes"].mean()

    if "pct_counts_mt" in adata.obs.columns:
        info["median_mt"] = adata.obs["pct_counts_mt"].median()
        info["mean_mt"] = adata.obs["pct_counts_mt"].mean()

    return info


def print_dataset_summary(datasets: Dict[str, AnnData]):
    """
    打印所有数据集的摘要信息

    Parameters
    ----------
    datasets : dict
        {dataset_name: AnnData}

    Example:
    -------
    >>> from scLucid.utils.data_loader import load_all_datasets, print_dataset_summary
    >>>
    >>> datasets = load_all_datasets()
    >>> print_dataset_summary(datasets)
    """
    print("\n" + "=" * 70)
    print("数据集摘要")
    print("=" * 70)

    for name, adata in datasets.items():
        info = get_dataset_info(adata)

        print(f"\n{name}:")
        print(f"  物种: {info['species']}")
        print(f"  组织类型: {info['tissue_type']}")
        print(f"  细胞数: {info['n_cells']:,}")
        print(f"  基因数: {info['n_genes']:,}")
        print(f"  批次数: {info['n_batches']}")

        if "median_n_genes" in info:
            print(f"  中位基因数: {info['median_n_genes']:.0f}")
            print(f"  平均基因数: {info['mean_n_genes']:.0f}")

        if "median_mt" in info:
            print(f"  中位线粒体%: {info['median_mt']:.1f}%")
            print(f"  平均线粒体%: {info['mean_mt']:.1f}%")

    print("\n" + "=" * 70)


def filter_by_species(datasets: Dict[str, AnnData], species: str) -> Dict[str, AnnData]:
    """
    按物种过滤数据集

    Parameters
    ----------
    datasets : dict
        {dataset_name: AnnData}
    species : str
        物种 ('human' 或 'mouse')

    Returns:
    -------
    filtered : dict
        只包含指定物种的数据集

    Example:
    -------
    >>> from scLucid.utils.data_loader import load_all_datasets, filter_by_species
    >>>
    >>> datasets = load_all_datasets()
    >>> human_datasets = filter_by_species(datasets, 'human')
    >>>
    >>> print(f"人源数据集: {list(human_datasets.keys())}")
    >>> # ['PBMC', 'LUAD']
    """
    filtered = {
        name: adata for name, adata in datasets.items() if adata.obs["species"].iloc[0] == species
    }

    return filtered


def filter_by_tissue_type(datasets: Dict[str, AnnData], tissue_type: str) -> Dict[str, AnnData]:
    """
    按组织类型过滤数据集

    Parameters
    ----------
    datasets : dict
        {dataset_name: AnnData}
    tissue_type : str
        组织类型 ('normal' 或 'tumor')

    Returns:
    -------
    filtered : dict
        只包含指定组织类型的数据集

    Example:
    -------
    >>> from scLucid.utils.data_loader import load_all_datasets, filter_by_tissue_type
    >>>
    >>> datasets = load_all_datasets()
    >>> tumor_datasets = filter_by_tissue_type(datasets, 'tumor')
    >>>
    >>> print(f"肿瘤数据集: {list(tumor_datasets.keys())}")
    >>> # ['LUAD', 'Melanoma']
    """
    filtered = {}

    for name, adata in datasets.items():
        ttype = adata.obs["tissue_type"].iloc[0]
        # 检查是否包含关键词
        if tissue_type.lower() in ttype.lower():
            filtered[name] = adata

    return filtered


__all__ = [
    "load_pbmc3k",
    "load_luad",
    "load_melanoma",
    "load_all_datasets",
    "get_dataset_info",
    "print_dataset_summary",
    "filter_by_species",
    "filter_by_tissue_type",
]
