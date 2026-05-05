"""
统一的细胞比例分析接口。

提供统一的 `analyze_celltype_proportion()` 函数，自动分发到：
- Pseudo-bulk 方法 (analysis/proportion.py)
- scCODA 方法 (tools/sccoda.py)

使用 `recommend_method()` 自动选择最合适的方法。
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Literal, Optional, Union

import pandas as pd
from anndata import AnnData

from .config import ProportionConfig
from .methods import ProportionMethod, recommend_method
from .pseudobulk import celltype_proportion_analysis as pb_analysis

log = logging.getLogger(__name__)


def analyze_celltype_proportion(
    adata: AnnData,
    method: Union[ProportionMethod, str, None] = None,
    config: Optional[ProportionConfig] = None,
    sample_col: str = "sample_id",
    condition_col: str = "condition",
    celltype_col: str = "cell_type",
    out_dir: Optional[Union[str, Path]] = None,
    return_type: Literal["auto", "tuple", "anndata"] = "auto",
    **kwargs,
) -> Union[tuple[pd.DataFrame, pd.DataFrame], AnnData]:
    """
    统一的细胞比例分析接口。

    这是分析细胞类型比例变化的统一入口函数。根据指定的方法（或自动推荐），
    调用相应的分析流程：Pseudo-bulk 或 scCODA。

    Parameters
    ----------
    adata : AnnData
        单细胞数据对象，必须包含样本、条件和细胞类型列
    method : {None, 'pseudobulk', 'sccoda'} or ProportionMethod, optional
        分析方法选择：
        - None (默认): 自动推荐最合适的方法
        - 'pseudobulk': Pseudo-bulk + 传统统计检验
        - 'sccoda': 贝叶斯组成数据分析
    config : ProportionConfig, optional
        分析配置对象。如果为 None，使用默认配置
    sample_col : str
        样本 ID 列名 (default: 'sample_id')
    condition_col : str
        条件列名 (default: 'condition')
    celltype_col : str
        细胞类型列名 (default: 'cell_type')
    out_dir : str or Path, optional
        输出目录，用于保存结果和图表
    return_type : {'auto', 'tuple', 'anndata'}
        返回类型：
        - 'auto': 根据方法自动选择（默认）
        - 'tuple': 总是返回 (prop_df, stat_df) 元组
        - 'anndata': 总是返回 AnnData 对象
    **kwargs
        传递给具体分析方法的额外参数

    Returns:
    -------
    根据方法和 return_type 参数：
        - **Pseudo-bulk**: (prop_df, stat_df) 元组
        - **scCODA**: AnnData（结果存储在 adata.uns['sclucid']['sccoda']）

    Raises:
    ------
    ValueError
        如果指定的方法无效
    ValueError
        如果指定的方法无效

    Examples:
    --------
    **示例 1: 自动推荐方法**
    >>> from scLucid.analysis import analyze_celltype_proportion
    >>>
    >>> # 自动选择最合适的方法
    >>> result = analyze_celltype_proportion(
    ...     adata,
    ...     sample_col="sample",
    ...     condition_col="condition"
    ... )

    **示例 2: 指定 Pseudo-bulk 方法**
    >>> from scLucid.analysis import ProportionConfig, analyze_celltype_proportion
    >>>
    >>> config = ProportionConfig(
    ...     test_method='wilcoxon',
    ...     plot_types=['bar', 'box', 'volcano']
    ... )
    >>>
    >>> prop_df, stat_df = analyze_celltype_proportion(
    ...     adata,
    ...     method='pseudobulk',
    ...     config=config,
    ...     out_dir='./results'
    ... )

    **示例 3: 使用 scCODA 方法**
    >>> result = analyze_celltype_proportion(
    ...     adata,
    ...     method='sccoda',
    ...     reference_cell_type='T_cells',
    ...     reference_level='control'
    ... )

    **示例 4: 先推荐再分析**
    >>> from scLucid.analysis.proportion_methods import recommend_method
    >>>
    >>> # 步骤 1: 推荐方法
    >>> method = recommend_method(
    ...     adata,
    ...     sample_col='sample',
    ...     condition_col='condition'
    ... )
    >>> print(f"推荐方法: {method}")
    >>>
    >>> # 步骤 2: 使用推荐方法
    >>> result = analyze_celltype_proportion(
    ...     adata,
    ...     method=method,
    ...     config=config
    ... )

    See Also:
    --------
    recommend_method : 根据数据特征推荐分析方法
    ProportionConfig : 分析配置类
    scLucid.analysis.proportion.run_sccoda : scCODA 详细接口

    Notes:
    -----
    **方法选择指南**:

    ┌─────────────┬────────────────┬─────────────────┬──────────────┐
    │ 方法        │ 样本量/组      │ 批次效应        │ 适用场景     │
    ├─────────────┼────────────────┼─────────────────┼──────────────┤
    │ Pseudo-bulk │ N ≥ 5         │ 无             │ 标准分析     │
    │ scCODA      │ N < 5 或任意  │ 有             │ 批次校正     │
    └─────────────┴────────────────┴─────────────────┴──────────────┘

    **输出格式**:

    - **Pseudo-bulk** 返回元组:
        * prop_df: 样本 × 细胞类型的比例矩阵
        * stat_df: 统计检验结果（p-value, 效应量等）

    - **scCODA** 返回 AnnData:
        * 结果存储在 adata.uns['sclucid']['sccoda']
        * 可使用对应的 summarise 函数提取结果

    **性能考虑**:

    - Pseudo-bulk: 最快（秒级）
    - scCODA: 中等（分钟级，MCMC 采样）

    References:
    ----------
    .. [1] Love, M.I. et al. (2014) "Moderated estimation of fold change and
       dispersion for RNA-seq data with DESeq2." Genome Biol.
    .. [2] Büttner, M. et al. (2021) "scCODA: Bayesian composition analysis
       of single-cell data." Nat. Methods.
    """
    # 1. 方法选择
    if method is None:
        log.info("未指定方法，自动推荐...")
        method = recommend_method(adata, sample_col=sample_col, condition_col=condition_col)
    elif isinstance(method, str):
        method = ProportionMethod(method.lower())

    log.info(f"使用方法: {method.value}")

    # 2. 设置输出目录
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

    # 3. 根据方法分发
    if method == ProportionMethod.PSEUDOBULK:
        # 使用 Pseudo-bulk 方法
        log.info("运行 Pseudo-bulk 分析...")

        # 确保配置存在
        if config is None:
            from .config import ProportionConfig as DefaultConfig

            config = DefaultConfig(
                celltype_col=celltype_col,
                sample_col=sample_col,
                condition_col=condition_col,
                out_dir=str(out_dir) if out_dir else None,
            )
        else:
            # 更新配置中的列名
            config.celltype_col = celltype_col
            config.sample_col = sample_col
            config.condition_col = condition_col
            if out_dir:
                config.out_dir = str(out_dir)

        # 运行 Pseudo-bulk 分析
        prop_df, stat_df = pb_analysis(adata, config)

        # 根据返回类型要求返回
        if return_type == "anndata":
            # 将结果存储到 AnnData 中
            adata.uns.setdefault("sclucid", {})["proportion"] = {
                "method": "pseudobulk",
                "prop_df": prop_df,
                "stat_df": stat_df,
                "config": config.model_dump(),
            }
            return adata
        else:
            return prop_df, stat_df

    elif method == ProportionMethod.SCCODA:
        # 使用 scCODA 方法
        log.info("运行 scCODA 分析...")

        try:
            from .sccoda import run_sccoda
        except ImportError:
            raise ImportError("scCODA 未安装。请安装: pip install scCODA")

        # 提取 scCODA 特定参数
        sccoda_kwargs = {
            "cell_type_col": celltype_col,
            "sample_col": sample_col,
            "condition_col": condition_col,
            "out_dir": out_dir,
        }
        sccoda_kwargs.update(kwargs)

        # 运行 scCODA
        result = run_sccoda(adata, **sccoda_kwargs)

        return result

    else:
        raise ValueError(
            f"未知的方法: {method}\n" f"可用方法: {[m.value for m in ProportionMethod]}"
        )


# 便捷别名
analyze_proportion = analyze_celltype_proportion
celltype_proportion_analysis = analyze_celltype_proportion


def analyze_all_methods(
    adata: AnnData,
    methods: list[ProportionMethod] = None,
    sample_col: str = "sample_id",
    condition_col: str = "condition",
    celltype_col: str = "cell_type",
    out_dir: Optional[Union[str, Path]] = None,
    compare: bool = True,
) -> dict[str, Union[tuple, AnnData]]:
    """
    使用多种方法分析细胞比例，便于结果比较。

    Parameters
    ----------
    adata : AnnData
        单细胞数据对象
    methods : list[ProportionMethod], optional
        要运行的方法列表，默认运行所有已实现的方法
    sample_col : str
        样本 ID 列名
    condition_col : str
        条件列名
    celltype_col : str
        细胞类型列名
    out_dir : str or Path, optional
        输出目录，各方法结果保存在子目录中
    compare : bool
        是否生成方法比较报告

    Returns:
    -------
    dict
        字典，键为方法名，值为对应方法的结果

    Examples:
    --------
    >>> from scLucid.analysis import analyze_all_methods
    >>>
    >>> results = analyze_all_methods(
    ...     adata,
    ...     methods=['pseudobulk', 'sccoda'],
    ...     out_dir='./comparison'
    ... )
    >>>
    >>> # 访问各方法结果
    >>> prop_df, stat_df = results['pseudobulk']
    >>> adata_sccoda = results['sccoda']

    Notes:
    -----
    该函数便于：
    - 方法验证：比较不同方法的结果一致性
    - 方法选择：查看哪个方法在你的数据上表现最好
    - 综合分析：结合多种方法的优势

    警告：运行多种方法会显著增加计算时间，尤其是 scCODA（MCMC 采样）。
    """
    if methods is None:
        # 运行所有已实现的方法
        methods = [ProportionMethod.PSEUDOBULK, ProportionMethod.SCCODA]

    log.info(f"使用 {len(methods)} 种方法进行分析: {[m.value for m in methods]}")

    results = {}
    out_dir = Path(out_dir) if out_dir else None

    for method in methods:
        log.info(f"\n{'='*60}")
        log.info(f"运行 {method.value} 方法")
        log.info(f"{'='*60}")

        # 为每个方法创建子目录
        method_out_dir = out_dir / method.value if out_dir else None

        try:
            result = analyze_celltype_proportion(
                adata,
                method=method,
                sample_col=sample_col,
                condition_col=condition_col,
                celltype_col=celltype_col,
                out_dir=method_out_dir,
                return_type="auto",
            )
            results[method.value] = result

        except Exception as e:
            log.error(f"{method.value} 分析失败: {e}")
            results[method.value] = None

    # 生成比较报告
    if compare:
        log.info("\n生成方法比较报告...")
        from .proportion_methods import compare_methods

        comparison = compare_methods(
            adata,
            methods=methods,
            sample_col=sample_col,
            condition_col=condition_col,
            celltype_col=celltype_col,
        )

        if out_dir:
            comparison.to_csv(out_dir / "method_comparison.csv", index=False)
            log.info(f"比较报告保存至: {out_dir / 'method_comparison.csv'}")

        # 打印比较表
        print("\n" + "=" * 80)
        print("方法适用性比较")
        print("=" * 80)
        print(comparison[["method", "overall_score", "recommendation"]].to_string(index=False))

    return results
