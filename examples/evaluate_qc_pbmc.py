"""
QC混合策略评估 - 简化版（使用PBMC数据）

由于LUAD和黑色素瘤数据需要预处理，这个简化版先使用PBMC数据
演示QC策略的对比评估。
"""

import sys
from pathlib import Path
from typing import Dict, Tuple
import time

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation

# 添加src到路径
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from scLucid.qc import (
    calculate_qc_metric,
    recommend_intelligent_qc,
    recommend_qc_strategy,
    QCRecommendation
)
from scLucid.preprocess import run_preprocessing
from scLucid.analysis import cluster_cells

# 设置scanpy verbosity
sc.settings.verbosity = 1


def load_pbmc_data(data_dir: Path) -> sc.AnnData:
    """加载PBMC数据"""
    print("=" * 70)
    print("加载PBMC3K数据集")
    print("=" * 70)

    adata = sc.read_h5ad(data_dir / "pbmc3k" / "pbmc3k_raw.h5ad")
    adata = calculate_qc_metric(adata)

    # 添加元数据
    adata.obs['sample_id'] = 'pbmc'
    adata.obs['tissue_type'] = 'normal'
    adata.obs['species'] = 'human'
    adata.obs['batch'] = 'pbmc_batch'

    print(f"\n✓ PBMC数据加载完成:")
    print(f"  细胞数: {adata.n_obs:,}")
    print(f"  基因数: {adata.n_vars:,}")
    print(f"  中位基因数: {adata.obs['n_genes'].median():.0f}")
    print(f"  中位线粒体%: {adata.obs['pct_counts_mt'].median():.1f}%")

    return adata


def apply_unified_threshold(adata: sc.AnnData) -> Tuple[sc.AnnData, dict]:
    """策略1: 统一阈值"""
    print("\n" + "=" * 70)
    print("策略1: 统一阈值 (Unified Thresholds)")
    print("=" * 70)

    # 使用传统固定阈值
    min_genes = 200  # 传统阈值
    max_mt = 20.0    # 传统阈值

    print(f"\n固定阈值:")
    print(f"  min_genes > {min_genes}")
    print(f"  pct_mt < {max_mt}%")

    # 应用过滤
    adata_filtered = adata[
        (adata.obs['n_genes'] > min_genes) &
        (adata.obs['pct_counts_mt'] < max_mt)
    ].copy()

    info = {
        'threshold_min_genes': min_genes,
        'threshold_max_mt': max_mt,
        'n_original': adata.n_obs,
        'n_filtered': adata_filtered.n_obs,
        'retention_rate': adata_filtered.n_obs / adata.n_obs * 100
    }

    print(f"\n结果:")
    print(f"  原始细胞: {info['n_original']:,}")
    print(f"  过滤后细胞: {info['n_filtered']:,}")
    print(f"  保留率: {info['retention_rate']:.1f}%")

    return adata_filtered, info


def apply_sample_specific(adata: sc.AnnData) -> Tuple[sc.AnnData, QCRecommendation]:
    """策略2: 样本特异性阈值"""
    print("\n" + "=" * 70)
    print("策略2: 样本特异性阈值 (Sample-Specific)")
    print("=" * 70)

    tissue_type = adata.obs['tissue_type'].iloc[0]

    # 获取智能推荐
    rec = recommend_intelligent_qc(
        adata,
        tissue_type=tissue_type,
        plot=False
    )

    print(f"\n智能推荐:")
    print(f"  策略: {rec.overall_strategy.value}")
    print(f"  min_genes: {rec.min_genes.threshold:.0f} "
          f"[95% CI: {rec.min_genes.ci_lower}-{rec.min_genes.ci_upper}]")
    print(f"  max_mt: {rec.max_mt_percent.threshold:.1f}% "
          f"[95% CI: {rec.max_mt_percent.ci_lower:.1f}-"
          f"{rec.max_mt_percent.ci_upper:.1f}]")
    print(f"  数据质量: {rec.data_quality_score:.1f}/100")
    print(f"  置信度: {rec.overall_confidence:.2f}")

    # 应用推荐阈值
    adata_filtered = adata[
        (adata.obs['n_genes'] > rec.min_genes.threshold) &
        (adata.obs['pct_counts_mt'] < rec.max_mt_percent.threshold)
    ].copy()

    print(f"\n结果:")
    print(f"  原始细胞: {adata.n_obs:,}")
    print(f"  过滤后细胞: {adata_filtered.n_obs:,}")
    print(f"  保留率: {adata_filtered.n_obs / adata.n_obs * 100:.1f}%")

    return adata_filtered, rec


def apply_hybrid_strategy(adata: sc.AnnData, rec: QCRecommendation) -> Tuple[sc.AnnData, dict]:
    """策略3: 混合策略"""
    print("\n" + "=" * 70)
    print("策略3: 混合策略 (Hybrid Approach)")
    print("=" * 70)

    # 对于单个样本，混合策略 = 样本特异性 + 全局约束
    # 这里我们模拟全局约束的效果

    threshold = rec.min_genes.threshold
    # 假设全局约束要求阈值在 [180, 210] 范围内
    lower_bound = 180
    upper_bound = 210

    threshold_constrained = np.clip(threshold, lower_bound, upper_bound)

    if threshold != threshold_constrained:
        print(f"\n阈值调整: {threshold:.0f} → {threshold_constrained:.0f}")
        print(f"  (应用全局约束: [{lower_bound}, {upper_bound}])")
    else:
        print(f"\n阈值: {threshold_constrained:.0f}")
        print(f"  (无需调整，已在约束范围内)")

    # 应用过滤
    adata_filtered = adata[
        (adata.obs['n_genes'] > threshold_constrained) &
        (adata.obs['pct_counts_mt'] < rec.max_mt_percent.threshold)
    ].copy()

    info = {
        'threshold': threshold_constrained,
        'n_original': adata.n_obs,
        'n_filtered': adata_filtered.n_obs,
        'retention_rate': adata_filtered.n_obs / adata.n_obs * 100,
        'adjusted': threshold != threshold_constrained
    }

    print(f"\n结果:")
    print(f"  原始细胞: {info['n_original']:,}")
    print(f"  过滤后细胞: {info['n_filtered']:,}")
    print(f"  保留率: {info['retention_rate']:.1f}%")

    return adata_filtered, info


def compare_strategies(adata: sc.AnnData):
    """对比三种QC策略"""
    print("\n" + "*" * 70)
    print("*" + " " * 68 + "*")
    print("*" + " " * 22 + "QC策略对比评估" + " " * 30 + "*")
    print("*" + " " * 68 + "*")
    print("*" * 70)

    results = {}

    # 策略1: 统一阈值
    adata_unified, info_unified = apply_unified_threshold(adata.copy())
    results['unified'] = {
        'adata': adata_unified,
        'info': info_unified
    }

    # 策略2: 样本特异性
    adata_specific, rec_specific = apply_sample_specific(adata.copy())
    results['sample_specific'] = {
        'adata': adata_specific,
        'info': {
            'threshold': rec_specific.min_genes.threshold,
            'n_original': adata.n_obs,
            'n_filtered': adata_specific.n_obs,
            'retention_rate': adata_specific.n_obs / adata.n_obs * 100
        }
    }

    # 策略3: 混合策略
    adata_hybrid, info_hybrid = apply_hybrid_strategy(adata.copy(), rec_specific)
    results['hybrid'] = {
        'adata': adata_hybrid,
        'info': info_hybrid
    }

    # 汇总结果
    print("\n" + "=" * 70)
    print("结果汇总")
    print("=" * 70)

    summary = pd.DataFrame({
        '策略': ['统一阈值', '样本特异性', '混合策略'],
        'min_genes阈值': [
            info_unified['threshold_min_genes'],
            results['sample_specific']['info']['threshold'],
            results['hybrid']['info']['threshold']
        ],
        '保留率(%)': [
            info_unified['retention_rate'],
            results['sample_specific']['info']['retention_rate'],
            results['hybrid']['info']['retention_rate']
        ],
        '保留细胞数': [
            info_unified['n_filtered'],
            results['sample_specific']['info']['n_filtered'],
            results['hybrid']['info']['n_filtered']
        ]
    })

    print("\n" + str(summary))
    print()

    # 找出最佳策略
    best_retention = summary.loc[summary['保留率(%)'].idxmax(), '策略']
    print(f"✓ 最高保留率: {best_retention} ({summary.loc[summary['保留率(%)'].idxmax(), '保留率(%)']:.1f}%)")

    print("\n" + "=" * 70)
    print("推荐结论")
    print("=" * 70)

    print("\n基于PBMC3K数据的评估:")
    print("1. 统一阈值: 使用固定值（min_genes > 200, pct_mt < 20%）")
    print("   - 简单但不够灵活")
    print("   - 不考虑数据分布")

    print("\n2. 样本特异性: 使用智能推荐（GMM + Bootstrap）")
    print("   - 数据驱动，适应数据分布")
    print("   - 提供95%置信区间")
    print("   - 考虑组织类型")

    print("\n3. 混合策略: 样本特异性 + 全局约束")
    print("   - 平衡适应性和可比性")
    print("   - 适合多样本场景")
    print("   - 推荐用于实际分析")

    print("\n" + "=" * 70)
    print("下一步")
    print("=" * 70)

    print("\n1. 数据准备:")
    print("   python examples/prepare_data.py")
    print("   (转换LUAD和黑色素瘤数据为h5ad格式)")

    print("\n2. 完整评估:")
    print("   python examples/evaluate_qc_strategies.py")
    print("   (使用所有三个数据集)")

    print("\n3. 预期结果:")
    print("   - LUAD（肿瘤）: 混合策略表现最佳")
    print("   - 黑色素瘤（多批次）: 混合策略 + 批次校正最佳")

    # 保存结果
    output_dir = Path("./qc_evaluation_results")
    output_dir.mkdir(exist_ok=True)

    summary.to_csv(output_dir / "strategy_comparison_pbmc.csv", index=False)

    print(f"\n✓ 结果已保存到: {output_dir}/strategy_comparison_pbmc.csv")

    return results


def main():
    """主函数"""
    # 数据目录
    data_dir = Path(__file__).parent.parent / "data"

    print("\n" + "*" * 70)
    print("*" + " " * 68 + "*")
    print("*" + " " * 10 + "QC混合策略评估 - PBMC3K数据集" + " " * 24 + "*")
    print("*" + " " * 68 + "*")
    print("*" * 70)

    # 加载PBMC数据
    adata = load_pbmc_data(data_dir)

    # 对比三种策略
    results = compare_strategies(adata)

    print("\n" + "=" * 70)
    print("评估完成！")
    print("=" * 70)


if __name__ == '__main__':
    main()
