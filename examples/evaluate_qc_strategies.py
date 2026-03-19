"""
QC混合策略全面评估脚本

使用data/文件夹中的三个数据集评估三种QC策略：
1. PBMC3K (人源，正常组织，单批次)
2. LUAD (人源，肺腺癌，肿瘤组织)
3. 黑色素瘤 (鼠源，多批次，肿瘤组织)

评估维度：
- 细胞保留率
- 聚类质量
- 批次效应残留
- 数据质量改善
- 计算效率

作者：scLucid Development Team
日期：2025-02-08
"""

import sys
from pathlib import Path
from typing import Dict, List, Tuple, Any
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


# =============================================================================
# 数据加载
# =============================================================================

def load_datasets(data_dir: Path) -> Dict[str, sc.AnnData]:
    """
    加载所有数据集

    Returns
    -------
    datasets : dict
        {dataset_name: AnnData}
    """
    print("=" * 70)
    print("加载数据集")
    print("=" * 70)

    datasets = {}

    # 1. PBMC3K (人源，正常组织，单批次)
    print("\n1. 加载PBMC3K数据集...")
    pbmc = sc.read_h5ad(data_dir / "pbmc3k" / "pbmc3k_raw.h5ad")
    pbmc = calculate_qc_metric(pbmc)
    pbmc.obs['sample_id'] = 'pbmc'
    pbmc.obs['tissue_type'] = 'normal'
    pbmc.obs['species'] = 'human'
    pbmc.obs['batch'] = 'pbmc_batch'  # 单批次
    datasets['PBMC'] = pbmc
    print(f"  ✓ PBMC: {pbmc.n_obs} 细胞, {pbmc.n_vars} 基因")

    # 2. LUAD (人源，肺腺癌，肿瘤组织)
    print("\n2. 加载LUAD数据集...")
    try:
        # 尝试不同的文件名
        luad_files = list((data_dir / "human_LUAD_GSE131907").glob("*.h5ad"))
        if luad_files:
            luad = sc.read_h5ad(luad_files[0])
            luad = calculate_qc_metric(luad)
            luad.obs['sample_id'] = 'luad'
            luad.obs['tissue_type'] = 'lung_tumor'
            luad.obs['species'] = 'human'
            # 假设有多批次
            if 'batch' not in luad.obs.columns:
                luad.obs['batch'] = 'luad_batch'  # 如果没有批次信息
            datasets['LUAD'] = luad
            print(f"  ✓ LUAD: {luad.n_obs} 细胞, {luad.n_vars} 基因")
        else:
            print(f"  ⚠ LUAD数据集未找到，跳过")
    except Exception as e:
        print(f"  ⚠ LUAD加载失败: {e}")

    # 3. 黑色素瘤 (鼠源，多批次，肿瘤组织)
    print("\n3. 加载黑色素瘤数据集...")
    try:
        melanoma_files = list((data_dir / "mouse_melanoma_GSE119352").glob("*.h5ad"))
        if melanoma_files:
            melanoma = sc.read_h5ad(melanoma_files[0])
            melanoma = calculate_qc_metric(melanoma)
            melanoma.obs['sample_id'] = 'melanoma'
            melanoma.obs['tissue_type'] = 'melanoma'
            melanoma.obs['species'] = 'mouse'
            # 多批次
            if 'batch' not in melanoma.obs.columns:
                # 创建模拟批次（用于演示）
                n_batches = 3
                melanoma.obs['batch'] = [
                    f'melanoma_batch_{i % n_batches}'
                    for i in range(melanoma.n_obs)
                ]
            datasets['Melanoma'] = melanoma
            print(f"  ✓ 黑色素瘤: {melanoma.n_obs} 细胞, {melanoma.n_vars} 基因")
        else:
            print(f"  ⚠ 黑色素瘤数据集未找到，跳过")
    except Exception as e:
        print(f"  ⚠ 黑色素瘤加载失败: {e}")

    print(f"\n总计: 加载了 {len(datasets)} 个数据集")

    return datasets


# =============================================================================
# QC策略实现
# =============================================================================

def apply_unified_thresholds(
    datasets: Dict[str, sc.AnnData]
) -> Dict[str, sc.AnnData]:
    """
    策略1: 统一阈值（传统方法）

    计算所有样本的全局QC指标分布，选择统一阈值
    """
    print("\n" + "=" * 70)
    print("策略1: 统一阈值 (Unified Thresholds)")
    print("=" * 70)

    # 计算全局阈值
    all_genes = []
    all_mt = []

    for name, adata in datasets.items():
        all_genes.extend(adata.obs['n_genes'].values)
        all_mt.extend(adata.obs['pct_counts_mt'].values)

    # 使用固定百分位数
    global_min_genes = np.percentile(all_genes, 10)  # 第10百分位数
    global_max_mt = np.percentile(all_mt, 90)  # 第90百分位数

    print(f"\n全局阈值:")
    print(f"  min_genes > {global_min_genes:.0f}")
    print(f"  pct_mt < {global_max_mt:.1f}%")

    # 应用统一阈值
    filtered = {}
    for name, adata in datasets.items():
        adata_filtered = adata[
            (adata.obs['n_genes'] > global_min_genes) &
            (adata.obs['pct_counts_mt'] < global_max_mt)
        ].copy()

        filtered[name] = adata_filtered

        retention_rate = len(adata_filtered) / len(adata) * 100
        print(f"\n{name}:")
        print(f"  原始: {len(adata)} 细胞")
        print(f"  过滤后: {len(adata_filtered)} 细胞")
        print(f"  保留率: {retention_rate:.1f}%")

    return filtered


def apply_sample_specific_thresholds(
    datasets: Dict[str, sc.AnnData]
) -> Tuple[Dict[str, sc.AnnData], Dict[str, QCRecommendation]]:
    """
    策略2: 样本特异性阈值

    每个样本独立计算QC指标，选择自己的最优阈值
    """
    print("\n" + "=" * 70)
    print("策略2: 样本特异性阈值 (Sample-Specific)")
    print("=" * 70)

    filtered = {}
    recommendations = {}

    for name, adata in datasets.items():
        tissue_type = adata.obs['tissue_type'].iloc[0]
        species = adata.obs['species'].iloc[0]

        print(f"\n{name}:")
        print(f"  组织类型: {tissue_type}")
        print(f"  物种: {species}")

        # 获取智能推荐
        rec = recommend_intelligent_qc(
            adata,
            tissue_type=tissue_type,
            plot=False
        )

        recommendations[name] = rec

        # 显示推荐
        print(f"  推荐策略: {rec.overall_strategy.value}")
        print(f"  min_genes: {rec.min_genes.threshold:.0f} "
              f"[95% CI: {rec.min_genes.ci_lower}-{rec.min_genes.ci_upper}]")
        print(f"  max_mt: {rec.max_mt_percent.threshold:.1f}% "
              f"[95% CI: {rec.max_mt_percent.ci_lower:.1f}-"
              f"{rec.max_mt_percent.ci_upper:.1f}]")
        print(f"  数据质量: {rec.data_quality_score:.1f}/100")

        # 应用阈值
        adata_filtered = adata[
            (adata.obs['n_genes'] > rec.min_genes.threshold) &
            (adata.obs['pct_counts_mt'] < rec.max_mt_percent.threshold)
        ].copy()

        filtered[name] = adata_filtered

        retention_rate = len(adata_filtered) / len(adata) * 100
        print(f"  保留率: {retention_rate:.1f}%")

    return filtered, recommendations


def apply_hybrid_strategy(
    datasets: Dict[str, sc.AnnData],
    recommendations: Dict[str, QCRecommendation]
) -> Dict[str, sc.AnnData]:
    """
    策略3: 混合策略（推荐）

    使用intelligent_qc获取样本特异性推荐，
    然后应用全局约束（median ± 3×MAD）
    """
    print("\n" + "=" * 70)
    print("策略3: 混合策略 (Hybrid Approach)")
    print("=" * 70)

    # 计算全局约束
    thresholds = [rec.min_genes.threshold for rec in recommendations.values()]
    global_median = np.median(thresholds)
    global_mad = median_abs_deviation(thresholds)

    lower_bound = global_median - 3 * global_mad
    upper_bound = global_median + 3 * global_mad

    print(f"\n全局约束:")
    print(f"  中位数: {global_median:.0f}")
    print(f"  MAD: {global_mad:.0f}")
    print(f"  约束范围: [{lower_bound:.0f}, {upper_bound:.0f}]")

    # 应用混合策略
    filtered = {}
    for name, adata in datasets.items():
        rec = recommendations[name]
        original_threshold = rec.min_genes.threshold

        # 应用约束
        threshold = np.clip(original_threshold, lower_bound, upper_bound)

        # 检查是否调整
        adjusted = ""
        if original_threshold < lower_bound:
            adjusted = f" (向上调整: {original_threshold:.0f} → {threshold:.0f})"
        elif original_threshold > upper_bound:
            adjusted = f" (向下调整: {original_threshold:.0f} → {threshold:.0f})"

        print(f"\n{name}:")
        print(f"  阈值: {threshold:.0f}{adjusted}")

        # 过滤
        adata_filtered = adata[
            (adata.obs['n_genes'] > threshold) &
            (adata.obs['pct_counts_mt'] < rec.max_mt_percent.threshold)
        ].copy()

        filtered[name] = adata_filtered

        retention_rate = len(adata_filtered) / len(adata) * 100
        print(f"  保留率: {retention_rate:.1f}%")

    return filtered


# =============================================================================
# 评估指标
# =============================================================================

def evaluate_cell_retention(
    original: Dict[str, sc.AnnData],
    filtered: Dict[str, sc.AnnData],
    strategy_name: str
) -> pd.DataFrame:
    """
    评估细胞保留率
    """
    print(f"\n{'='*70}")
    print(f"细胞保留率评估 - {strategy_name}")
    print(f"{'='*70}")

    results = []

    for name in original.keys():
        n_original = original[name].n_obs
        n_filtered = filtered[name].n_obs
        retention_rate = n_filtered / n_original * 100

        mean_genes = filtered[name].obs['n_genes'].mean()
        median_mt = filtered[name].obs['pct_counts_mt'].median()

        results.append({
            'dataset': name,
            'n_original': n_original,
            'n_filtered': n_filtered,
            'retention_rate': retention_rate,
            'mean_genes': mean_genes,
            'median_mt': median_mt
        })

        print(f"\n{name}:")
        print(f"  保留率: {retention_rate:.1f}%")
        print(f"  平均基因数: {mean_genes:.0f}")
        print(f"  中位线粒体%: {median_mt:.1f}%")

    df = pd.DataFrame(results)
    df.set_index('dataset', inplace=True)

    return df


def evaluate_clustering_quality(
    filtered: Dict[str, sc.AnnData],
    strategy_name: str
) -> pd.DataFrame:
    """
    评估聚类质量

    注意：这会运行完整的preprocessing和clustering，可能需要一些时间
    """
    print(f"\n{'='*70}")
    print(f"聚类质量评估 - {strategy_name}")
    print(f"{'='*70}")

    results = []

    for name, adata in filtered.items():
        print(f"\n{name}:")

        try:
            # Preprocessing
            print("  运行preprocessing...")
            adata_pp = run_preprocessing(adata.copy())

            # Clustering
            print("  运行clustering...")
            adata_pp = cluster_cells(adata_pp, resolution=0.8)

            # 计算聚类质量指标
            from sklearn.metrics import silhouette_score, davies_bouldin_score

            silhouette = silhouette_score(
                adata_pp.obsm['X_pca'],
                adata_pp.obs['leiden']
            )

            davies_bouldin = davies_bouldin_score(
                adata_pp.obsm['X_pca'],
                adata_pp.obs['leiden']
            )

            n_clusters = adata_pp.obs['leiden'].nunique()

            results.append({
                'dataset': name,
                'silhouette_score': silhouette,
                'davies_bouldin_score': davies_bouldin,
                'n_clusters': n_clusters
            })

            print(f"  Silhouette score: {silhouette:.3f} (越高越好)")
            print(f"  Davies-Bouldin score: {davies_bouldin:.3f} (越低越好)")
            print(f"  聚类数: {n_clusters}")

        except Exception as e:
            print(f"  ⚠ 评估失败: {e}")
            results.append({
                'dataset': name,
                'silhouette_score': np.nan,
                'davies_bouldin_score': np.nan,
                'n_clusters': np.nan
            })

    df = pd.DataFrame(results)
    df.set_index('dataset', inplace=True)

    return df


def evaluate_batch_effect(
    filtered: Dict[str, sc.AnnData],
    strategy_name: str
) -> pd.DataFrame:
    """
    评估批次效应残留

    注意：需要多个批次的数据
    """
    print(f"\n{'='*70}")
    print(f"批次效应评估 - {strategy_name}")
    print(f"{'='*70}")

    results = []

    for name, adata in filtered.items():
        print(f"\n{name}:")

        # 检查是否有批次信息
        if 'batch' not in adata.obs.columns:
            print(f"  ⚠ 无批次信息，跳过")
            results.append({
                'dataset': name,
                'n_batches': 1,
                'batch_mixing': 'N/A'
            })
            continue

        n_batches = adata.obs['batch'].nunique()

        # 简单的批次混合评估
        # 如果批次混合得好，每个batch在UMAP上的分布应该是均匀的
        # 这里我们使用一个简化的指标

        try:
            # Preprocessing for visualization
            adata_pp = run_preprocessing(adata.copy())

            # 计算每个batch的细胞在近邻中的比例
            # 这是一个简化的批次混合评估
            from sklearn.neighbors import NearestNeighbors

            nbrs = NearestNeighbors(n_neighbors=30)
            nbrs.fit(adata_pp.obsm['X_pca'])
            distances, indices = nbrs.kneighbors(adata_pp.obsm['X_pca'])

            # 计算batch混合度
            batch_mixing_scores = []
            for i in range(len(adata_pp)):
                neighbor_batches = adata_pp.obs['batch'].iloc[indices[i]]
                own_batch = adata_pp.obs['batch'].iloc[i]
                mixing_score = (neighbor_batches != own_batch).sum() / 30
                batch_mixing_scores.append(mixing_score)

            avg_mixing = np.mean(batch_mixing_scores)

            results.append({
                'dataset': name,
                'n_batches': n_batches,
                'batch_mixing': avg_mixing
            })

            print(f"  批次数: {n_batches}")
            print(f"  批次混合度: {avg_mixing:.3f} (越高越好，1=完全混合)")

        except Exception as e:
            print(f"  ⚠ 评估失败: {e}")
            results.append({
                'dataset': name,
                'n_batches': n_batches,
                'batch_mixing': np.nan
            })

    df = pd.DataFrame(results)
    df.set_index('dataset', inplace=True)

    return df


# =============================================================================
# 主函数
# =============================================================================

def main():
    """主函数：运行完整的QC策略对比"""

    print("\n" + "*" * 70)
    print("*" + " " * 68 + "*")
    print("*" + " " * 15 + "QC混合策略全面评估" + " " * 31 + "*")
    print("*" + " " * 68 + "*")
    print("*" * 70)

    # 数据目录
    data_dir = Path(__file__).parent.parent / "data"

    if not data_dir.exists():
        print(f"\n❌ 数据目录不存在: {data_dir}")
        return

    # ========================================
    # 1. 加载数据
    # ========================================
    datasets = load_datasets(data_dir)

    if len(datasets) == 0:
        print("\n❌ 没有可用的数据集")
        return

    # ========================================
    # 2. 应用三种QC策略
    # ========================================

    all_results = {}

    # 策略1: 统一阈值
    filtered_unified = apply_unified_thresholds(datasets)
    all_results['unified'] = filtered_unified

    # 策略2: 样本特异性
    filtered_specific, recommendations = apply_sample_specific_thresholds(datasets)
    all_results['sample_specific'] = filtered_specific

    # 策略3: 混合策略
    filtered_hybrid = apply_hybrid_strategy(datasets, recommendations)
    all_results['hybrid'] = filtered_hybrid

    # ========================================
    # 3. 评估所有策略
    # ========================================

    print("\n" + "=" * 70)
    print("开始评估所有策略")
    print("=" * 70)

    evaluation_results = {}

    for strategy_name, filtered_datasets in all_results.items():
        # 细胞保留率
        retention_df = evaluate_cell_retention(
            datasets,
            filtered_datasets,
            strategy_name
        )

        # 聚类质量（可选，耗时较长）
        print(f"\n是否运行聚类质量评估？(y/n): ", end='')
        # response = input()  # 取消注释以启用交互
        response = 'n'  # 默认跳过，节省时间

        if response.lower() == 'y':
            clustering_df = evaluate_clustering_quality(
                filtered_datasets,
                strategy_name
            )
        else:
            print("  跳过聚类质量评估（节省时间）")
            clustering_df = None

        # 批次效应（可选）
        batch_df = evaluate_batch_effect(
            filtered_datasets,
            strategy_name
        )

        evaluation_results[strategy_name] = {
            'retention': retention_df,
            'clustering': clustering_df,
            'batch_effect': batch_df
        }

    # ========================================
    # 4. 汇总结果
    # ========================================

    print("\n" + "=" * 70)
    print("结果汇总")
    print("=" * 70)

    # 细胞保留率对比
    print("\n" + "-" * 70)
    print("细胞保留率对比")
    print("-" * 70)

    retention_comparison = pd.DataFrame({
        strategy: eval_res['retention']['retention_rate']
        for strategy, eval_res in evaluation_results.items()
    })
    print(retention_comparison)

    # 计算平均保留率
    print("\n平均保留率:")
    for strategy in retention_comparison.columns:
        avg_rate = retention_comparison[strategy].mean()
        print(f"  {strategy}: {avg_rate:.1f}%")

    # ========================================
    # 5. 保存结果
    # ========================================

    output_dir = Path("./qc_evaluation_results")
    output_dir.mkdir(exist_ok=True)

    # 保存所有评估结果
    with open(output_dir / "evaluation_summary.txt", 'w') as f:
        f.write("QC混合策略评估结果\n")
        f.write("=" * 70 + "\n\n")

        f.write("细胞保留率对比:\n")
        retention_comparison.to_csv(f)

        if any(eval_res['clustering'] is not None
               for eval_res in evaluation_results.values()):
            f.write("\n\n聚类质量对比:\n")
            for strategy, eval_res in evaluation_results.items():
                if eval_res['clustering'] is not None:
                    f.write(f"\n{strategy}:\n")
                    eval_res['clustering'].to_csv(f)

        f.write("\n\n批次效应对比:\n")
        for strategy, eval_res in evaluation_results.items():
            f.write(f"\n{strategy}:\n")
            eval_res['batch_effect'].to_csv(f)

    print(f"\n✓ 结果已保存到: {output_dir}")

    # ========================================
    # 6. 最终建议
    # ========================================

    print("\n" + "=" * 70)
    print("最终建议")
    print("=" * 70)

    # 使用决策树获取推荐
    strategy_rec, rationale = recommend_qc_strategy(
        datasets,
        batch_key='batch',
        tissue_key='tissue_type'
    )

    print("\n基于数据特征的自动推荐:")
    for line in rationale:
        print(line)

    print("\n" + "=" * 70)
    print("评估完成！")
    print("=" * 70)


if __name__ == '__main__':
    main()
