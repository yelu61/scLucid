"""
QC策略决策树 - 帮助用户选择合适的QC策略

这个模块提供了一个决策树框架，根据数据特征推荐最合适的QC策略。
"""

from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from anndata import AnnData


class QCStrategyDecisionTree:
    """
    QC策略决策树

    根据数据特征自动推荐最合适的QC策略：
    - 统一阈值 (Unified)
    - 样本特异性 (Sample-Specific)
    - 混合策略 (Hybrid)
    """

    def __init__(self):
        """初始化决策树"""
        self.decision_factors = {}
        self.recommendation = None
        self.rationale = []

    def assess_data_characteristics(
        self,
        samples_dict: Dict[str, AnnData],
        batch_key: str = 'batch',
        tissue_key: str = 'tissue_type'
    ) -> Dict:
        """
        评估数据特征

        Parameters
        ----------
        samples_dict : dict
            {sample_id: AnnData}
        batch_key : str
            批次信息的列名
        tissue_key : str
            组织类型的列名

        Returns
        -------
        characteristics : dict
            包含所有评估特征的字典
        """
        characteristics = {}

        # 1. 批次数目
        if batch_key and batch_key in list(samples_dict.values())[0].obs.columns:
            all_batches = []
            for adata in samples_dict.values():
                all_batches.extend(adata.obs[batch_key].unique())
            n_batches = len(set(all_batches))
        else:
            n_batches = len(samples_dict)

        characteristics['n_batches'] = n_batches

        # 2. 组织类型数目
        if tissue_key and tissue_key in list(samples_dict.values())[0].obs.columns:
            all_tissues = []
            for adata in samples_dict.values():
                all_tissues.extend(adata.obs[tissue_key].unique())
            n_tissues = len(set(all_tissues))
        else:
            n_tissues = 1

        characteristics['n_tissues'] = n_tissues

        # 3. 计算样本间QC指标差异
        all_medians = {'n_genes': [], 'n_counts': [], 'pct_counts_mt': []}

        for sample_id, adata in samples_dict.items():
            for metric in all_medians.keys():
                if metric in adata.obs.columns:
                    all_medians[metric].append(adata.obs[metric].median())

        # 计算变异系数 (CV)
        cvs = {}
        for metric, values in all_medians.items():
            if len(values) > 0 and np.mean(values) > 0:
                cvs[metric] = np.std(values) / np.mean(values)

        characteristics['cvs'] = cvs

        # 4. 平均CV (跨所有指标)
        if cvs:
            characteristics['mean_cv'] = np.mean(list(cvs.values()))
        else:
            characteristics['mean_cv'] = 0.0

        # 5. 样本量
        characteristics['n_samples'] = len(samples_dict)

        # 6. 总细胞数
        total_cells = sum(adata.n_obs for adata in samples_dict.values())
        characteristics['total_cells'] = total_cells

        # 7. 每个样本平均细胞数
        characteristics['mean_cells_per_sample'] = total_cells / len(samples_dict)

        self.decision_factors = characteristics
        return characteristics

    def recommend_strategy(
        self,
        samples_dict: Dict[str, AnnData],
        batch_key: str = 'batch',
        tissue_key: str = 'tissue_type'
    ) -> Tuple[str, List[str]]:
        """
        推荐QC策略

        Parameters
        ----------
        samples_dict : dict
            {sample_id: AnnData}
        batch_key : str
            批次信息的列名
        tissue_key : str
            组织类型的列名

        Returns
        -------
        strategy : str
            推荐的策略: 'unified', 'sample_specific', 'hybrid'
        rationale : list
            推荐理由
        """
        # 评估数据特征
        characteristics = self.assess_data_characteristics(
            samples_dict, batch_key, tissue_key
        )

        rationale = []

        # 决策树逻辑
        # ========================================

        # 决策1: 批次数目
        if characteristics['n_batches'] == 1:
            rationale.append(
                f"✓ 单批次数据 (n_batches={characteristics['n_batches']})"
            )
            single_batch = True
        else:
            rationale.append(
                f"⚠ 多批次数据 (n_batches={characteristics['n_batches']})"
            )
            single_batch = False

        # 决策2: 组织类型
        if characteristics['n_tissues'] == 1:
            rationale.append(
                f"✓ 单一组织类型 (n_tissues={characteristics['n_tissues']})"
            )
            single_tissue = True
        else:
            rationale.append(
                f"⚠ 多组织类型 (n_tissues={characteristics['n_tissues']})"
            )
            single_tissue = False

        # 决策3: 样本间变异
        mean_cv = characteristics['mean_cv']
        if mean_cv < 0.1:
            rationale.append(
                f"✓ 样本间变异小 (CV={mean_cv:.2f} < 0.1)"
            )
            low_variation = True
        elif mean_cv < 0.3:
            rationale.append(
                f"⚠ 样本间变异中等 (CV={mean_cv:.2f})"
            )
            low_variation = False
        else:
            rationale.append(
                f"❌ 样本间变异大 (CV={mean_cv:.2f} > 0.3)"
            )
            low_variation = False

        # 最终决策
        # ========================================
        if single_batch and single_tissue and low_variation:
            # 场景1: 简单场景
            strategy = 'unified'
            rationale.append("\n🎯 推荐策略: 统一阈值 (Unified Thresholds)")
            rationale.append("   原因: 单批次、单组织、低变异")

        elif not single_batch and not single_tissue and not low_variation:
            # 场景3: 极端异质性
            strategy = 'sample_specific'
            rationale.append("\n🎯 推荐策略: 样本特异性 + 批次校正")
            rationale.append("   原因: 多批次、多组织、高变异")
            rationale.append("   提示: 需要配合Harmony/BBKNN/scVI")

        else:
            # 场景2: 中等复杂度（最常见）
            strategy = 'hybrid'
            rationale.append("\n🎯 推荐策略: 混合策略 (Hybrid Approach)")
            rationale.append("   原因: 平衡全局一致性和局部适应性")
            rationale.append("   方法: intelligent_qc + 全局约束 (median ± 3×MAD)")

        self.recommendation = strategy
        self.rationale = rationale

        return strategy, rationale

    def print_recommendation(self):
        """打印推荐结果"""
        if not self.rationale:
            print("请先调用 recommend_strategy() 方法")
            return

        print("\n" + "=" * 70)
        print("QC策略推荐结果")
        print("=" * 70)

        for line in self.rationale:
            print(line)

        # 打印详细建议
        print("\n" + "=" * 70)
        print("详细建议")
        print("=" * 70)

        if self.recommendation == 'unified':
            print("""
1. 计算所有样本的QC指标分布
2. 选择统一阈值:
   - min_genes: 第5-10百分位数
   - max_mt_percent: 第85-90百分位数
3. 应用到所有样本

示例代码:
   from scLucid.qc import run_standard_qc
   adata_filtered = run_standard_qc(adata, min_genes=200, max_mt_percent=20)
""")

        elif self.recommendation == 'sample_specific':
            print("""
1. 每个样本独立QC
2. 使用intelligent_qc获取样本特异性推荐
3. 合并样本后进行批次校正

示例代码:
   from scLucid.qc import recommend_intelligent_qc
   import harmonypy as hm

   # 样本特异性QC
   for sample_id, adata in samples.items():
       rec = recommend_intelligent_qc(adata)
       adata = adata[adata.obs['n_genes'] > rec.min_genes.threshold]

   # 合并和批次校正
   adata_combined = ad.concat(samples, join='outer')
   adata_combined.obsm['X_pca_harmony'] = hm.run_harmony(...)
""")

        elif self.recommendation == 'hybrid':
            print("""
1. 使用intelligent_qc获取每个样本的推荐
2. 计算全局约束 (median ± 3×MAD)
3. 应用约束和过滤

示例代码:
   from scLucid.qc import recommend_intelligent_qc
   import numpy as np
   from scipy.stats import median_abs_deviation

   # 1. 获取推荐
   recommendations = {}
   thresholds = []
   for sample_id, adata in samples.items():
       rec = recommend_intelligent_qc(adata, tissue_type='lung_tumor')
       recommendations[sample_id] = rec
       thresholds.append(rec.min_genes.threshold)

   # 2. 全局约束
   global_median = np.median(thresholds)
   global_mad = median_abs_deviation(thresholds)
   lower_bound = global_median - 3 * global_mad
   upper_bound = global_median + 3 * global_mad

   # 3. 应用约束和过滤
   for sample_id, adata in samples.items():
       rec = recommendations[sample_id]
       threshold = np.clip(rec.min_genes.threshold, lower_bound, upper_bound)
       adata_filtered = adata[adata.obs['n_genes'] > threshold]
""")

        print("=" * 70)


# 便捷函数
def recommend_qc_strategy(
    samples_dict: Dict[str, AnnData],
    batch_key: str = 'batch',
    tissue_key: str = 'tissue_type'
) -> Tuple[str, List[str]]:
    """
    便捷函数：推荐QC策略

    Parameters
    ----------
    samples_dict : dict
        {sample_id: AnnData}
    batch_key : str
        批次信息的列名
    tissue_key : str
        组织类型的列名

    Returns
    -------
    strategy : str
        推荐的策略: 'unified', 'sample_specific', 'hybrid'
    rationale : list
        推荐理由

    Examples
    --------
    >>> from scLucid.qc import recommend_qc_strategy
    >>>
    >>> # 假设有多个样本
    >>> samples = {'sample1': adata1, 'sample2': adata2}
    >>>
    >>> # 获取推荐
    >>> strategy, rationale = recommend_qc_strategy(
    ...     samples,
    ...     batch_key='batch',
    ...     tissue_key='tissue_type'
    ... )
    >>>
    >>> print(f"推荐策略: {strategy}")
    >>> for line in rationale:
    ...     print(line)
    """
    decision_tree = QCStrategyDecisionTree()
    strategy, rationale = decision_tree.recommend_strategy(
        samples_dict, batch_key, tissue_key
    )
    decision_tree.print_recommendation()

    return strategy, rationale


# 使用示例
if __name__ == '__main__':
    print("""
QC策略决策树 - 使用示例
====================================

# 示例1: 简单场景（单批次、单组织）
samples = {'rep1': adata1, 'rep2': adata2, 'rep3': adata3}
strategy, rationale = recommend_qc_strategy(samples)
# 预期: 'unified'

# 示例2: 复杂场景（多批次、肿瘤-正常）
samples = {
    'tumor_01': adata_tumor1,
    'tumor_02': adata_tumor2,
    'normal_01': adata_normal1,
    'normal_02': adata_normal2,
}
strategy, rationale = recommend_qc_strategy(
    samples,
    batch_key='batch',
    tissue_key='tissue_type'
)
# 预期: 'hybrid' 或 'sample_specific'

# 示例3: 极端异质性
samples = {
    'batch1_tumor': adata1,
    'batch2_normal': adata2,
    'batch3_tumor': adata3,
}
strategy, rationale = recommend_qc_strategy(
    samples,
    batch_key='batch',
    tissue_key='tissue_type'
)
# 预期: 'sample_specific'
""")
