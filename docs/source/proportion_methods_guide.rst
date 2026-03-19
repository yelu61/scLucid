# 细胞比例分析方法选择指南

本文档说明如何在三种细胞比例分析方法之间选择和使用：
- **Pseudo-bulk**: 聚合到样本级别 + 传统统计
- **scCODA**: 贝叶斯组成数据分析
- **Milo**: 基于邻域的细胞水平分析

## 快速开始

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
# 输出: 'sccoda', 'pseudobulk', 或 'milo'
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

**原理**：聚合到样本级别，使用 bulk RNA-seq 统计方法（DESeq2, t-test, Wilcoxon等）

**优势**：
- ✅ 成熟稳定，文献广泛接受
- ✅ 统计功效高（样本级聚合）
- ✅ DESeq2 对低丰度细胞类型稳健
- ✅ 易于解释

**劣势**：
- ❌ 忽略细胞间异质性
- ❌ 丢失单细胞分辨率信息

**适用场景**：
- 每组 N ≥ 5 样本
- 无明显批次效应
- 细胞类型注释完整
- 关注细胞类型水平变化

**使用示例**：

```python
from scLucid.analysis import analyze_celltype_proportion, ProportionConfig

# 配置分析
config = ProportionConfig(
    test_method='wilcoxon',  # 统计方法
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
```

---

### 2. scCODA 方法

**原理**：贝叶斯层次模型，专门为单细胞数据设计

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

**状态**：⚠️ **尚未实现**（计划中）

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
| **样本要求** | N ≥ 5/组 | N < 5/组 | N ≥ 3/组 |
| **批次效应** | ❌ 不适合 | ✅ 专门处理 | ⚠️ 部分处理 |
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

适合：大样本量（N≥10/组），无批次效应

```python
from scLucid.analysis import analyze_celltype_proportion, ProportionConfig

config = ProportionConfig(
    test_method='deseq2',  # DESeq2 对大样本效果好
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

- **Pseudo-bulk**: 检测平均比例变化
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
