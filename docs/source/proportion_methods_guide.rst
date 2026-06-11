# 细胞比例分析方法选择指南

本文档说明如何在三种细胞比例分析方法之间选择和使用：
- **Pseudo-bulk / CLR**: 聚合到 biological sample 级别，优先使用 CLR 转换后的组成数据检验
- **Pseudo-bulk / covariate-aware logCPM**: 当存在 batch、patient、paired design 时，使用样本级线性模型
- **scCODA**: 可选贝叶斯组成数据分析 backend
- **Milo**: 基于邻域的细胞水平分析，当前仍为计划接口

## 快速开始

.. warning::

   细胞比例是 compositional data。原始比例上的 ``t-test`` / ``wilcoxon``
   在 scLucid 中仅保留为 legacy exploratory 路径。正式解释应优先使用
   sample-level CLR、DESeq2-style 计数模型或带协变量的样本级模型。

### 自动推荐方法

最简单的方式是让系统自动推荐最合适的方法：

```python
from scLucid.analysis import analyze_celltype_proportion

# 自动推荐并分析
result = analyze_celltype_proportion(
    adata,
    sample_col="sample_id",
    condition_col="condition"
)

# 系统会打印推荐的原因，例如：
# INFO: 推荐方法: sccoda
# INFO: 原因: 样本量=3 < 5, 批次效应=True
```

### 查看方法推荐但不运行

```python
from scLucid.analysis import recommend_method

# 获取推荐方法
method = recommend_method(
    adata,
    sample_col="sample",
    condition_col="condition"
)

print(f"推荐方法: {method.value}")
# 输出: 'pseudobulk' 或 'sccoda'（Milo 仍为计划接口）
```

### 比较所有方法的适用性

```python
from scLucid.analysis import compare_methods

# 生成方法比较表
comparison = compare_methods(adata)

print(comparison[['method', 'overall_score', 'recommendation']])
# 输出:
#          method  overall_score  recommendation
#      pseudobulk           0.80  ✅ 强烈推荐
#          sccoda           0.65  ⚠️  可用
#            milo           0.45  ❌ 不推荐
```

---

## 方法详解

### 1. Pseudo-bulk 方法

**原理**：聚合到 biological sample 级别，默认使用 compositional-aware
CLR 检验；在有 batch/patient 等设计因素时，可使用
``linear_model_logcpm`` 进行样本级协变量建模。

``composition_transform(method="clr")`` 会先判断输入是否已经是闭合比例
（行和约等于 1）、百分比（行和约等于 100）、原始 counts，或 0-1 范围但
未闭合的子组成。除已闭合比例外，其余非负输入都会先按样本行闭合后再做
CLR；负值或非有限值会报错，避免静默产生错误的组成数据结果。

**优势**：
- ✅ 成熟稳定，文献广泛接受
- ✅ 统计功效高（样本级聚合）
- ✅ 避免把 cell 当作独立 biological replicate
- ✅ 支持 CLR、paired/batch-aware 路径和 FDR
- ✅ 易于解释

**劣势**：
- ❌ 忽略细胞间异质性
- ❌ 丢失单细胞分辨率信息
- ❌ 原始比例检验不适合正式 compositional inference

**适用场景**：
- 每组至少有 biological replicates
- 有明确 sample、condition、可选 batch/patient 元数据
- 细胞类型注释完整
- 关注细胞类型水平变化

**使用示例**：

```python
from scLucid.analysis import analyze_celltype_proportion, ProportionConfig

# 配置分析
config = ProportionConfig(
    test_method='clr-t-test',  # 推荐：sample-level CLR 检验
    plot_types=['bar', 'box', 'volcano'],  # 可视化
    out_dir='./results'
)

# 运行分析
prop_df, stat_df = analyze_celltype_proportion(
    adata,
    method='pseudobulk',
    config=config
)

# 查看结果
print(stat_df[stat_df['padj'] < 0.05])  # 显著的细胞类型

# 结果表会包含 inference_level / compositional_data_warning 等审计字段。
```

### Pseudo-bulk 条件 DE：带 batch/patient 的样本级模型

当目标是基因层面的 condition DE，而不是细胞比例变化，请使用
``run_pseudobulk_de``。如果存在 batch 或 patient blocking，优先使用
纯 Python 的 ``linear_model_logcpm`` 路径：

```python
from scLucid.analysis import PseudobulkDEConfig, run_pseudobulk_de

config = PseudobulkDEConfig(
    sample_col="sampleID",
    condition_key="condition",
    groupby="cell_type_final",
    group_names=["T cell", "B cell"],
    contrasts=[("control", "treated")],
    method="linear_model_logcpm",
    design_covariates=["batch"],
    block_col="patient_id",   # 可选；会作为 categorical covariate 进入模型
    min_cells_per_sample=10,
)

de_df = run_pseudobulk_de(adata, config)

# 默认 robust_cov_type="HC3"，使用异方差稳健标准误。
# 如需复现普通 OLS 标准误，可显式设置 robust_cov_type="nonrobust"。

# 正式结果应检查：
# - inference_level == "sample_level"
# - valid_for_publication_inference == True
# - design_formula / design_covariates / block_col / covariance_type
```

每组只有一个 biological sample 时，scLucid 默认返回
``descriptive_single_sample`` effect-size-only 结果，不给正式 p 值。
强制 cell-level fallback 会被标记为 ``exploratory_cell_level``。

---

### 2. scCODA 方法

**原理**：贝叶斯组成数据模型，专门用于 compositional abundance
变化。scLucid 中该路径是可选 backend；使用前需确认依赖、reference
cell type 和采样参数。

**优势**：
- ✅ 处理批次效应
- ✅ 适合小样本（N<5）
- ✅ 提供可信区间
- ✅ 多条件比较友好

**劣势**：
- ❌ MCMC 采样较慢
- ❌ 贝叶斯模型调参复杂
- ❌ 文献接受度较新

**适用场景**：
- 每组 N < 5 样本
- 存在批次效应
- 需要贝叶斯可信区间
- 多条件比较

**使用示例**：

```python
from scLucid.analysis import analyze_celltype_proportion

# 运行 scCODA 分析
adata_result = analyze_celltype_proportion(
    adata,
    method='sccoda',
    reference_cell_type='T_cells',  # 参考细胞类型
    reference_level='control',       # 参考条件
    n_samples=25000,                  # MCMC 采样数
    out_dir='./results'
)

# 结果存储在 adata.uns
sccoda_results = adata_result.uns['sclucid']['sccoda']

# 查看显著变化
print(sccoda_results['final_results'])
```

**scCODA 专属功能**：

```python
from scLucid.tools import (
    run_sccoda,
    summarize_sccoda,
    plot_sccoda_proportion_with_significance
)

# 运行分析
adata = run_sccoda(
    adata,
    cell_type_col='cell_type',
    sample_col='sample_id',
    condition_col='condition'
)

# 汇总结果
summary = summarize_sccoda(adata)

# 绘图
plot_sccoda_proportion_with_significance(
    adata,
    condition='condition',
    save_path='./sccoda_plot.pdf'
)
```

---

### 3. Milo 方法

**状态**：⚠️ **尚未实现**（计划中）。不要把该接口作为当前
production 分析路径。

**原理**：在 UMAP/PCA 空间定义邻域，检验邻域细胞组成变化

**优势**：
- ✅ 保留单细胞分辨率
- ✅ 检测亚群水平变化
- ✅ 无需预先定义细胞类型
- ✅ 可视化空间模式

**劣势**：
- ❌ 计算复杂度高
- ❌ 需要调参（邻域大小）
- ❌ 结果解释较复杂

**适用场景**：
- 需要检测亚群变化
- 细胞类型注释可能不完整
- 关注空间分布模式
- 发现新的细胞亚群

**未来使用示例（计划）**：

```python
# 待实现
from scLucid.analysis import analyze_celltype_proportion

adata_result = analyze_celltype_proportion(
    adata,
    method='milo',
    n_neighbors=30,     # 邻域大小
    n_pcs=30,           # PCs 数量
    alpha=0.1,          # 显著性阈值
    out_dir='./results'
)
```

---

## 方法对比表

| 特性 | Pseudo-bulk | scCODA | Milo |
|------|-------------|--------|------|
| **样本要求** | biological replicates | 小样本可选 | N ≥ 3/组 |
| **批次效应** | ✅ ``linear_model_logcpm`` 可显式建模 | ✅ 可建模 | ⚠️ 部分处理 |
| **空间分辨率** | ❌ 无 | ❌ 无 | ✅ 保留 |
| **计算速度** | ⚡⚡⚡ 快 | ⚡⚡ 中等 | ⚡ 慢 |
| **统计功效** | ✅ 高 | ⚠️ 中等 | ⚠️ 中等 |
| **结果解释** | ✅ 简单 | ⚠️ 复杂 | ⚠️ 复杂 |
| **文献接受度** | ✅ 高 | ⚠️ 中等 | ⚠️ 中等 |
| **成熟度** | ✅ 成熟 | ⚠️ 较新 | ⚠️ 较新 |

---

## 工作流程示例

### 工作流 1: 自动化分析（推荐）

适合：不确定使用哪种方法，希望系统自动选择

```python
from scLucid.analysis import analyze_celltype_proportion, recommend_method

# 步骤 1: 查看推荐
method = recommend_method(adata)
print(f"推荐方法: {method.value}")

# 步骤 2: 运行分析（使用推荐方法）
result = analyze_celltype_proportion(adata)

# 步骤 3: 提取结果
if isinstance(result, tuple):
    prop_df, stat_df = result
    sig_celltypes = stat_df[stat_df['padj'] < 0.05]
else:  # scCODA 返回 AnnData
    sccoda_results = result.uns['sclucid']['sccoda']
```

### 工作流 2: 方法比较验证

适合：需要验证结果一致性，选择最合适的方法

```python
from scLucid.analysis import analyze_all_methods

# 运行所有方法并比较
results = analyze_all_methods(
    adata,
    methods=['pseudobulk', 'sccoda'],
    out_dir='./comparison',
    compare=True
)

# 查看比较报告
# 保存于: ./comparison/method_comparison.csv

# 比较结果
pb_prop, pb_stat = results['pseudobulk']
adata_sccoda = results['sccoda']

# 比较 p-values
import matplotlib.pyplot as plt
plt.scatter(pb_stat['pval'], sccoda_stat['pval'])
plt.xlabel('Pseudo-bulk p-value')
plt.ylabel('scCODA p-value')
plt.savefig('./comparison/pval_correlation.pdf')
```

### 工作流 3: 大样本标准分析

适合：大样本量（N≥10/组），批次/配对信息清楚

```python
from scLucid.analysis import analyze_celltype_proportion, ProportionConfig

config = ProportionConfig(
    test_method='clr-t-test',
    plot_types=['bar', 'box', 'heatmap', 'volcano']
)

prop_df, stat_df = analyze_celltype_proportion(
    adata,
    method='pseudobulk',
    config=config
)
```

### 工作流 4: 小样本批次校正

适合：小样本量（N<5/组），存在批次效应

```python
from scLucid.analysis import analyze_celltype_proportion

adata_result = analyze_celltype_proportion(
    adata,
    method='sccoda',
    reference_cell_type='T_cells',
    reference_level='control',
    # scCODA 特定参数
    n_samples=25000,
    n_burnin=5000
)
```

---

## 常见问题

### Q1: 我应该选择哪种方法？

**A**: 使用自动推荐功能：

```python
from scLucid.analysis import recommend_method, compare_methods

# 快速推荐
method = recommend_method(adata)

# 详细比较
comparison = compare_methods(adata)
print(comparison[['method', 'overall_score', 'recommendation']])
```

### Q2: 可以同时使用多种方法吗？

**A**: 可以！使用 `analyze_all_methods()`：

```python
from scLucid.analysis import analyze_all_methods

results = analyze_all_methods(
    adata,
    methods=['pseudobulk', 'sccoda'],
    out_dir='./comparison',
    compare=True
)
```

### Q3: Pseudo-bulk 和 scCODA 结果不一致怎么办？

**A**: 这是正常的，它们检测不同类型的变化：

- **Pseudo-bulk/CLR**: 检测 sample-level 组成变化
- **scCODA**: 检测组成变化（考虑组成性质）

建议：
1. 查看两种方法都显著的变化（高置信度）
2. 使用生物学知识判断哪种结果更合理
3. 考虑样本量和批次效应的影响

### Q4: Milo 什么时候会实现？

**A**: Milo 在开发路线图上。当前建议：

- **短期**：使用 Pseudo-bulk + 手动亚群分析
- **中期**：使用聚类 + Milo（需单独实现）
- **长期**：集成到 scLucid 统一接口

你可以通过 GitHub Issues 请求 Milo 功能优先级提升。

---

## 参考资料

1. **DESeq2**: Love et al., *Genome Biology* 2014
2. **scCODA**: Büttner et al., *Nature Methods* 2021
3. **Milo**: Dan et al., *Nature Methods* 2022
4. **scCODA 教程**: https://github.com/theislab/scCODA
5. **Milo 教程**: https://github.com/MarioniLab/milo

---

## 联系与反馈

如有问题或建议，请：
- 提交 GitHub Issue
- 查看 scLucid 文档
- 联系开发团队
