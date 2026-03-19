"""
细胞比例分析方法选择和推荐模块。

提供三种细胞比例差异分析方法：
1. Pseudo-bulk: 聚合到样本级别 + 传统统计检验
2. scCODA: 贝叶斯组成数据分析（处理批次效应）
3. Milo: 基于邻域的细胞水平分析（保留空间异质性）

使用方法推荐函数根据数据特征自动选择最合适的方法。
"""

from __future__ import annotations

import logging
from enum import Enum
from typing import Literal, Optional, Union

import pandas as pd
from anndata import AnnData
from pydantic import Field

from .config import MethodSelectionConfig

log = logging.getLogger(__name__)


class ProportionMethod(str, Enum):
    """
    细胞比例分析方法枚举。

    Attributes
    ----------
    PSEUDOBULK : str
        Pseudo-bulk 方法：聚合到样本级别，使用传统统计检验（DESeq2, t-test, Wilcoxon等）
        适用场景：样本量充足（每组 N≥5），无明显批次效应
        优势：成熟稳定，统计功效高，易于解释
    SCCODA : str
        scCODA 方法：贝叶斯组成数据分析
        适用场景：样本量少（每组 N<5），存在批次效应
        优势：处理批次效应，提供可信区间，多条件友好
    MILO : str
        Milo 方法：基于邻域的细胞水平差异分析
        适用场景：需要检测亚群水平变化，关注空间分布
        优势：保留单细胞分辨率，检测局部变化
    """
    PSEUDOBULK = "pseudobulk"
    SCCODA = "sccoda"
    MILO = "milo"

    @classmethod
    def get_description(cls) -> dict:
        """获取各方法的描述信息."""
        return {
            cls.PSEUDOBULK: {
                "name": "Pseudo-bulk Analysis",
                "description": "聚合到样本级别，使用传统统计检验",
                "best_for": "样本充足（N≥5/组），无批次效应",
                "output": "(prop_df, stat_df) 元组",
                "ref": "Love et al., 2014 (DESeq2)"
            },
            cls.SCCODA: {
                "name": "scCODA",
                "description": "贝叶斯组成数据分析",
                "best_for": "样本少（N<5/组），有批次效应",
                "output": "AnnData (结果在 .uns)",
                "ref": "Büttner et al., 2021"
            },
            cls.MILO: {
                "name": "Milo",
                "description": "基于邻域的细胞水平分析",
                "best_for": "需要检测亚群变化，关注空间分布",
                "output": "AnnData (结果在 .uns)",
                "ref": "Dan et al., 2022"
            }
        }


def recommend_method(
    adata: AnnData,
    sample_col: str = "sample_id",
    condition_col: str = "condition",
    config: Optional[MethodSelectionConfig] = None
) -> ProportionMethod:
    """
    根据数据特征推荐最合适的细胞比例分析方法。

    该函数分析数据特征（样本量、细胞类型数量、批次效应等），
    自动推荐最合适的分析方法。

    Parameters
    ----------
    adata : AnnData
        单细胞数据对象
    sample_col : str
        样本 ID 列名
    condition_col : str
        条件列名
    config : MethodSelectionConfig, optional
        方法选择配置，如果为 None 则自动推断

    Returns
    -------
    ProportionMethod
        推荐的分析方法

    Examples
    --------
    >>> from scLucid.analysis.proportion_methods import recommend_method
    >>>
    >>> # 示例 1: 自动推荐
    >>> method = recommend_method(adata, sample_col="sample", condition_col="condition")
    >>> print(f"推荐方法: {method.value}")
    推荐方法: sccoda

    >>> # 示例 2: 提供配置
    >>> from scLucid.analysis.proportion_methods import MethodSelectionConfig
    >>> config = MethodSelectionConfig(
    ...     n_samples_per_group=3,
    ...     has_batch_effect=True
    ... )
    >>> method = recommend_method(adata, config=config)

    See Also
    --------
    analyze_celltype_proportion : 使用推荐方法进行分析

    Notes
    -----
    **推荐逻辑**:

    1. **优先级 1: 空间分辨率需求**
       - 如果 spatial_resolution=True → 推荐 Milo
       - 理由：只有 Milo 能保留空间异质性信息

    2. **优先级 2: 批次效应 + 小样本**
       - 如果 has_batch_effect=True 或 n_samples_per_group < 5 → 推荐 scCODA
       - 理由：scCODA 专为处理批次效应和小样本设计

    3. **默认: Pseudo-bulk**
       - 大样本（N≥5），无批次效应 → 推荐 Pseudo-bulk
       - 理由：成熟稳定，统计功效高，易于解释

    **每种方法的优势**:

    - **Pseudo-bulk**:
      * 文献接受度高，结果可重复
      * DESeq2 对低丰度细胞类型稳健
      * 统计功效高于单细胞水平方法
      * 适合验证性研究

    - **scCODA**:
      * 贝叶斯框架处理不确定性
      * 自动处理批次效应
      * 提供可信区间
      * 适合探索性研究

    - **Milo**:
      * 检测亚群水平变化
      * 保留空间分布模式
      * 无需预先定义细胞类型
      * 适合发现新亚群
    """
    if config is None:
        # 自动推断数据特征
        log.info("自动分析数据特征...")

        # 1. 计算样本量
        samples_per_group = adata.obs.groupby(condition_col)[sample_col].nunique()
        n_samples = samples_per_group.min()

        # 2. 计算细胞类型数量
        n_celltypes = adata.obs[condition_col].nunique() if condition_col in adata.obs else 10

        # 3. 检测批次效应（简化判断）
        has_batch = 'batch' in adata.obs.columns and adata.obs['batch'].nunique() > 1

        log.info(f"数据特征: {n_samples} 样本/组, {n_celltypes} 细胞类型, 批次效应={has_batch}")

        # 自动配置
        config = MethodSelectionConfig(
            n_samples_per_group=n_samples,
            n_celltypes=n_celltypes,
            has_batch_effect=has_batch
        )

    # 推荐逻辑
    if config.spatial_resolution:
        log.info("推荐 Milo 方法（需要空间分辨率）")
        return ProportionMethod.MILO

    if config.has_batch_effect or config.n_samples_per_group < 5:
        log.info(f"推荐 scCODA 方法（批次效应={config.has_batch_effect}, "
                f"样本量={config.n_samples_per_group} < 5）")
        return ProportionMethod.SCCODA

    log.info(f"推荐 Pseudo-bulk 方法（大样本={config.n_samples_per_group} ≥ 5, 无批次效应）")
    return ProportionMethod.PSEUDOBULK


def compare_methods(
    adata: AnnData,
    methods: list[ProportionMethod] = None,
    sample_col: str = "sample_id",
    condition_col: str = "condition",
    celltype_col: str = "cell_type"
) -> pd.DataFrame:
    """
    比较不同方法的适用性评分。

    为每种方法计算适用性评分，帮助用户选择最合适的方法。

    Parameters
    ----------
    adata : AnnData
        单细胞数据对象
    methods : list[ProportionMethod], optional
        要比较的方法列表，默认比较所有方法
    sample_col : str
        样本 ID 列名
    condition_col : str
        条件列名
    celltype_col : str
        细胞类型列名

    Returns
    -------
    pd.DataFrame
        方法比较表，包含各方法的适用性评分

    Examples
    --------
    >>> from scLucid.analysis.proportion_methods import compare_methods
    >>>
    >>> comparison = compare_methods(adata)
    >>> print(comparison[['method', 'overall_score', 'recommendation']])
    """
    if methods is None:
        methods = list(ProportionMethod)

    # 分析数据特征
    samples_per_group = adata.obs.groupby(condition_col)[sample_col].nunique()
    n_samples = samples_per_group.min()
    n_celltypes = adata.obs[celltype_col].nunique()
    has_batch = 'batch' in adata.obs.columns and adata.obs['batch'].nunique() > 1

    results = []

    for method in methods:
        scores = {}

        if method == ProportionMethod.PSEUDOBULK:
            # Pseudo-bulk 评分
            scores['样本充足度'] = min(n_samples / 10, 1.0)  # 10+ 样本满分
            scores['批次稳健性'] = 0.0 if has_batch else 1.0  # 不适合批次效应
            scores['细胞类型友好'] = min(n_celltypes / 20, 1.0)
            scores['计算效率'] = 1.0  # 最快
            scores['文献接受度'] = 1.0  # 最高
            scores['亚群分辨率'] = 0.0  # 无

        elif method == ProportionMethod.SCCODA:
            # scCODA 评分
            scores['样本充足度'] = 0.8 if n_samples < 5 else 0.6
            scores['批次稳健性'] = 1.0  # 最适合批次效应
            scores['细胞类型友好'] = 0.9
            scores['计算效率'] = 0.5  # MCMC 较慢
            scores['文献接受度'] = 0.7
            scores['亚群分辨率'] = 0.3

        elif method == ProportionMethod.MILO:
            # Milo 评分
            scores['样本充足度'] = min(n_samples / 5, 1.0)
            scores['批次稳健性'] = 0.5
            scores['细胞类型友好'] = 0.7
            scores['计算效率'] = 0.3  # 最慢
            scores['文献接受度'] = 0.6
            scores['亚群分辨率'] = 1.0  # 最高

        # 计算总分
        scores['overall_score'] = sum(scores.values()) / len(scores)

        # 推荐建议
        if scores['overall_score'] >= 0.7:
            scores['recommendation'] = '✅ 强烈推荐'
        elif scores['overall_score'] >= 0.5:
            scores['recommendation'] = '⚠️  可用'
        else:
            scores['recommendation'] = '❌ 不推荐'

        results.append({
            'method': method.value,
            **scores
        })

    df = pd.DataFrame(results)
    df = df.sort_values('overall_score', ascending=False)

    return df


# 便捷别名
MethodRecommender = recommend_method
MethodComparator = compare_methods
