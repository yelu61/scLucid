# QC策略对比分析：统一阈值 vs 样本特异性阈值

## 📋 问题背景

**当前策略**：所有样本统一阈值
- 计算所有样本的QC指标分布
- 选取一个统一的阈值标准
- 过滤所有细胞

**替代策略**：每个样本使用自己的阈值
- 每个样本独立计算QC指标
- 每个样本选择自己的最优阈值
- 分别过滤

---

## 🔬 方法学对比分析

### 策略A：统一阈值（Unified Thresholds）

#### 实现方式
```python
# 方法1：全局分布
import numpy as np

# 合并所有样本
all_genes = np.concatenate([
    sample1.obs['n_genes'],
    sample2.obs['n_genes'],
    sample3.obs['n_genes']
])

# 基于全局分布选择阈值
global_threshold = np.percentile(all_genes, 5)  # 第5百分位数

# 应用到所有样本
for sample in [sample1, sample2, sample3]:
    sample = sample[sample.obs['n_genes'] > global_threshold]
```

#### 优点 ✅
1. **跨样本可比性**
   - 所有样本使用相同标准
   - 避免样本间的系统性偏差
   - 下游分析更可靠

2. **保守策略**
   - 不会过度过滤任何一个样本
   - 保留所有样本的"正常"细胞

3. **简单可重复**
   - 阈值选择透明
   - 易于文档化
   - 便于重现

4. **适合批次效应较小的场景**
   - 如果所有样本用相同协议制备
   - 技术变异小于生物变异

#### 缺点 ❌
1. **忽略样本特异性**
   - 不同组织类型有不同特性（肿瘤 vs 正常）
   - 不同制备方案（新鲜 vs 冷冻）
   - 不同测序深度

2. **可能欠过滤或过过滤**
   - 对低质量样本：阈值太宽松，保留低质量细胞
   - 对高质量样本：阈值太严格，丢失有价值的细胞

3. **不适合异质性数据**
   - 肿瘤-正常混合
   - 多时间点
   - 多批次

---

### 策略B：样本特异性阈值（Sample-Specific Thresholds）

#### 实现方式
```python
# 方法2：样本特定阈值
from scLucid.qc import recommend_intelligent_qc

samples = {
    'sample1': sample1,
    'sample2': sample2,
    'sample3': sample3
}

for sample_id, adata in samples.items():
    # 获取样本特异性推荐
    recommendation = recommend_intelligent_qc(
        adata,
        tissue_type=adata.obs['tissue_type'].iloc[0]
    )

    # 使用样本特异性阈值
    threshold = recommendation.min_genes.threshold

    print(f"{sample_id}: min_genes = {threshold} "
          f"[95% CI: {recommendation.min_genes.ci_lower}-"
          f"{recommendation.min_genes.ci_upper}]")

    adata = adata[adata.obs['n_genes'] > threshold]
```

#### 优点 ✅
1. **适应样本特性**
   - 考虑组织类型（肿瘤 vs 正常）
   - 考虑制备方案
   - 考虑测序深度

2. **数据驱动**
   - 基于每个样本的数据分布
   - 更客观
   - 有置信区间

3. **适合异质性数据**
   - 处理批次效应
   - 处理不同组织类型
   - 处理不同时间点

4. **提高敏感性和特异性**
   - 低质量样本：更严格的过滤
   - 高质量样本：保留更多细胞

#### 缺点 ❌
1. **引入样本间偏差**
   - 不同标准可能引入系统性差异
   - 需要后续批次校正

2. **复杂度增加**
   - 需要为每个样本单独评估
   - 需要验证每个样本的阈值

3. **可能过拟合**
   - 某些样本可能选择极端阈值
   - 需要正则化或约束

---

## 🎯 推荐策略：**混合方法**（Hybrid Approach）

### 核心思想
**全局指导 + 局部适应**

```python
from scLucid.qc import recommend_intelligent_qc, IntelligentQCRecommender

def hybrid_qc_strategy(
    samples_dict,
    global_constraint=True,
    tissue_aware=True
):
    """
    混合QC策略：全局约束 + 样本自适应

    Parameters
    ----------
    samples_dict : dict
        {sample_id: AnnData}
    global_constraint : bool
        是否应用全局约束
    tissue_aware : bool
        是否考虑组织类型
    """

    recommendations = {}
    all_thresholds = []

    # 第1步：获取每个样本的推荐阈值
    for sample_id, adata in samples_dict.items():
        tissue_type = adata.obs.get('tissue_type', ['unknown'])[0]

        rec = recommend_intelligent_qc(
            adata,
            tissue_type=tissue_type,
            plot=False
        )

        recommendations[sample_id] = rec
        all_thresholds.append(rec.min_genes.threshold)

    # 第2步：应用全局约束（如果启用）
    if global_constraint:
        # 计算全局阈值范围
        global_median = np.median(all_thresholds)
        global_mad = median_abs_deviation(all_thresholds)

        # 约束：样本阈值不能偏离全局中位数太远
        lower_bound = global_median - 3 * global_mad
        upper_bound = global_median + 3 * global_mad

        # 应用约束
        for sample_id, rec in recommendations.items():
            original_threshold = rec.min_genes.threshold

            # 如果样本阈值超出范围，调整到边界
            if rec.min_genes.threshold < lower_bound:
                rec.min_genes.threshold = lower_bound
                print(f"{sample_id}: threshold adjusted upward "
                      f"({original_threshold:.0f} -> {lower_bound:.0f})")

            elif rec.min_genes.threshold > upper_bound:
                rec.min_genes.threshold = upper_bound
                print(f"{sample_id}: threshold adjusted downward "
                      f"({original_threshold:.0f} -> {upper_bound:.0f})")

    # 第3步：应用过滤
    filtered_samples = {}
    for sample_id, adata in samples_dict.items():
        rec = recommendations[sample_id]

        # 使用调整后的阈值
        adata_filtered = adata[
            adata.obs['n_genes'] > rec.min_genes.threshold
        ].copy()

        filtered_samples[sample_id] = adata_filtered

        # 报告
        print(f"\n{sample_id}:")
        print(f"  Original: {len(adata)} cells")
        print(f"  Filtered: {len(adata_filtered)} cells")
        print(f"  Retained: {len(adata_filtered)/len(adata):.1%}")
        print(f"  Threshold: {rec.min_genes.threshold:.0f} "
              f"[95% CI: {rec.min_genes.ci_lower}-{rec.min_genes.ci_upper}]")

    return filtered_samples, recommendations
```

---

## 📊 策略对比表

| 维度 | 统一阈值 | 样本特异性 | 混合策略（推荐） |
|------|---------|-----------|----------------|
| **跨样本可比性** | ✅ 高 | ❌ 低 | ✅ 中高 |
| **适应性** | ❌ 低 | ✅ 高 | ✅ 高 |
| **鲁棒性** | ✅ 高 | ⚠️ 中 | ✅ 高 |
| **复杂度** | ✅ 低 | ❌ 高 | ⚠️ 中 |
| **适合批次效应** | ❌ 否 | ✅ 是 | ✅ 是 |
| **适合异质性** | ❌ 否 | ✅ 是 | ✅ 是 |
| **可重复性** | ✅ 高 | ⚠️ 中 | ✅ 中高 |
| **推荐使用场景** | 简单项目 | 复杂多批次 | **通用** |

---

## 🎯 场景推荐

### 场景1：同批次、同组织类型（简单）
**推荐：统一阈值**

```python
# 所有样本来自同一批次、同一组织类型
# 例如：PBMC的3个重复样本

# 使用全局阈值
threshold = np.percentile(all_genes, 5)
for sample in samples:
    sample = sample[sample.obs['n_genes'] > threshold]
```

**原因**：
- 样本间技术变异小
- 统一阈值保证可比性
- 简单可靠

---

### 场景2：多批次、不同组织类型（复杂）
**推荐：混合策略**

```python
# 样本来自不同批次、不同组织类型
# 例如：5个肿瘤样本 + 5个正常样本

filtered, recs = hybrid_qc_strategy(
    samples_dict,
    global_constraint=True,  # 全局约束
    tissue_aware=True  # 组织感知
)
```

**原因**：
- 批次效应明显
- 肿瘤vs正常组织特性不同
- 需要平衡适应性和可比性

---

### 场景3：极端异质性（如肿瘤-正常混合）
**推荐：样本特异性 + 批次校正**

```python
# 每个样本独立QC
# 然后使用Harmony/BBKNN进行批次校正

for sample_id, adata in samples.items():
    # 样本特异性QC
    rec = recommend_intelligent_qc(
        adata,
        tissue_type=adata.obs['tissue_type'][0]
    )
    adata = adata[adata.obs['n_genes'] > rec.min_genes.threshold]

# 合并样本
adata_combined = ad.concat(samples, join='outer')

# 批次校正（关键！）
import scvi
import harmonypy as hm

# 使用Harmony
adata_combined.obsm['X_pca_harmony'] = hm.run_harmony(
    adata_combined.obsm['X_pca'],
    adata_combined.obs,
    'batch'
)
```

**原因**：
- 样本间差异太大
- 需要样本特异性处理
- 批次校正可以去除技术偏差

---

## 🔬 scLucid的推荐实现

### 基于intelligent_qc的混合策略

```python
from scLucid.qc import IntelligentQCRecommender

class UnifiedQCRecommender:
    """
    统一QC推荐器：平衡全局一致性和样本特异性
    """

    def __init__(self, constraint_type='median_mad'):
        """
        Parameters
        ----------
        constraint_type : str
            'none': 无约束（纯样本特异性）
            'median_mad': 基于中位数和MAD的约束（推荐）
            'quantile': 基于分位数的约束
        """
        self.constraint_type = constraint_type

    def recommend_unified(
        self,
        samples_dict,
        tissue_key=None,
        plot=True
    ):
        """
        为多个样本生成统一的QC推荐

        Returns
        -------
        recommendations : dict
            {sample_id: QCRecommendation}
        summary : pd.DataFrame
            所有样本的推荐摘要
        """
        recommendations = {}
        sample_thresholds = []

        # 第1步：获取每个样本的推荐
        for sample_id, adata in samples_dict.items():
            tissue_type = 'unknown'
            if tissue_key and tissue_key in adata.obs:
                tissue_type = adata.obs[tissue_type].iloc[0]

            recommender = IntelligentQCRecommender(strategy='auto')
            rec = recommender.recommend(
                adata,
                tissue_type=tissue_type,
                plot=False
            )

            recommendations[sample_id] = rec
            sample_thresholds.append(rec.min_genes.threshold)

        # 第2步：计算全局约束
        if self.constraint_type == 'median_mad':
            global_median = np.median(sample_thresholds)
            global_mad = median_abs_deviation(sample_thresholds)

            lower_bound = global_median - 3 * global_mad
            upper_bound = global_median + 3 * global_mad

            print(f"\n全局约束:")
            print(f"  中位数: {global_median:.0f}")
            print(f"  MAD: {global_mad:.0f}")
            print(f"  范围: [{lower_bound:.0f}, {upper_bound:.0f}]")

            # 第3步：应用约束
            for sample_id, rec in recommendations.items():
                original = rec.min_genes.threshold

                if rec.min_genes.threshold < lower_bound:
                    rec.min_genes.threshold = lower_bound
                    rec.min_genes.ci_lower = lower_bound
                elif rec.min_genes.threshold > upper_bound:
                    rec.min_genes.threshold = upper_bound
                    rec.min_genes.ci_upper = upper_bound

                if original != rec.min_genes.threshold:
                    print(f"  {sample_id}: {original:.0f} -> "
                          f"{rec.min_genes.threshold:.0f}")

        # 第4步：生成摘要
        summary_data = []
        for sample_id, rec in recommendations.items():
            summary_data.append({
                'sample_id': sample_id,
                'min_genes': rec.min_genes.threshold,
                'ci_lower': rec.min_genes.ci_lower,
                'ci_upper': rec.min_genes.ci_upper,
                'strategy': rec.overall_strategy.value,
                'confidence': rec.overall_confidence
            })

        summary = pd.DataFrame(summary_data)

        if plot:
            self._plot_thresholds(summary, sample_thresholds)

        return recommendations, summary

    def _plot_thresholds(self, summary, original_thresholds):
        """可视化阈值分布和约束"""
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # 左图：阈值分布
        axes[0].hist(original_thresholds, bins=10, alpha=0.7, edgecolor='black')
        axes[0].axvline(np.median(original_thresholds), color='red',
                       linestyle='--', linewidth=2, label='全局中位数')
        axes[0].set_xlabel('min_genes 阈值')
        axes[0].set_ylabel('样本数')
        axes[0].set_title('各样本推荐阈值分布')
        axes[0].legend()

        # 右图：样本对比
        x_pos = np.arange(len(summary))
        axes[1].bar(x_pos, summary['min_genes'],
                   yerr=[summary['min_genes'] - summary['ci_lower'],
                         summary['ci_upper'] - summary['min_genes']],
                   alpha=0.7, capsize=5, edgecolor='black')
        axes[1].axhline(np.median(original_thresholds), color='red',
                       linestyle='--', linewidth=2, label='全局中位数')
        axes[1].set_xticks(x_pos)
        axes[1].set_xticklabels(summary['sample_id'], rotation=45)
        axes[1].set_ylabel('min_genes 阈值')
        axes[1].set_title('各样本阈值对比（含95% CI）')
        axes[1].legend()

        plt.tight_layout()
        plt.savefig('qc_thresholds_comparison.pdf', dpi=300, bbox_inches='tight')
        plt.close()

        print("\n图表已保存: qc_thresholds_comparison.pdf")
```

---

## 📝 使用示例

```python
# 假设有多个样本
samples = {
    'tumor_01': adata_tumor1,
    'tumor_02': adata_tumor2,
    'normal_01': adata_normal1,
    'normal_02': adata_normal2,
}

# 使用混合策略
unified_rec = UnifiedQCRecommender(constraint_type='median_mad')
recommendations, summary = unified_rec.recommend_unified(
    samples,
    tissue_key='tissue_type',
    plot=True
)

# 应用过滤
filtered_samples = {}
for sample_id, adata in samples.items():
    rec = recommendations[sample_id]

    adata_filtered = adata[
        (adata.obs['n_genes'] > rec.min_genes.threshold) &
        (adata.obs['pct_counts_mt'] < rec.max_mt_percent.threshold)
    ].copy()

    filtered_samples[sample_id] = adata_filtered

# 合并样本
adata_combined = ad.concat(filtered_samples, join='outer')

# 后续批次校正
# ...
```

---

## 🎯 总结和建议

### 1. **统一阈值**（简单场景）
- ✅ 适合：同批次、同组织、少量样本
- ❌ 不适合：多批次、异质性大
- 📊 风险：可能欠过滤或过过滤

### 2. **样本特异性阈值**（极端异质）
- ✅ 适合：极端异质性、不同组织类型
- ❌ 不适合：需要跨样本可比性的分析
- 📊 风险：引入系统性偏差

### 3. **混合策略**（推荐！⭐⭐⭐⭐⭐）
- ✅ 适合：大多数实际场景
- ✅ 平衡全局一致性和局部适应性
- ✅ 提供可视化和证据
- 📊 **scLucid的推荐策略**

### 实施建议

1. **默认使用混合策略**
   - 使用`intelligent_qc`获取样本特异性推荐
   - 应用全局约束（median ± 3×MAD）
   - 可视化和记录所有阈值

2. **极端情况使用样本特异性 + 批次校正**
   - 样本间差异太大
   - 必须配合批次校正（Harmony, BBKNN, scVI）

3. **简单情况使用统一阈值**
   - 同批次、同组织
   - 技术变异小

### 关键原则

**"在保证可比性的前提下，最大化适应性"**

- 可比性：通过全局约束保证
- 适应性：通过样本特异性推荐实现
- 平衡点：混合策略

---

## 🔬 未来改进方向

1. **自适应约束强度**
   - 根据样本间差异自动调整约束强度
   - 批次效应大 → 弱约束
   - 批次效应小 → 强约束

2. **基于组织的分组约束**
   - 肿瘤样本组内约束
   - 正常样本组内约束
   - 而非全局统一约束

3. **动态阈值优化**
   - 使用下游分析质量优化阈值
   - 例如：聚类质量、注释质量

4. **Bayesian分层模型**
   - 更优雅地处理全局和局部信息
   - 自动borrowing strength

---

**这是scLucid可以发表的另一个方法学创新！** 🚀
