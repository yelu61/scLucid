# scLucid 功能完善与流程设计行动计划

## 📋 当前后续工作优先级

基于当前架构状态（已完成模块化、抽象基类、命名规范），建议按以下优先级推进：

---

## 🎯 第一阶段：核心流程验证（1-2周）

### 目标
确保 QC → Preprocess → Analysis 主流程可以无缝衔接

### 任务清单

#### 1.1 数据流验证
- [ ] **输入输出接口标准化**
  - 检查每个模块的 AnnData 输入要求
  - 统一 `.obs` 和 `.uns` 的键名规范
  - 确保数据格式一致性

- [ ] **配置传递链路**
  ```
  QCWorkflowConfig → PreprocessingWorkflowConfig → AnalysisWorkflowConfig
  ```
  - 验证配置类之间的参数传递
  - 添加配置继承/复用机制

#### 1.2 端到端流程测试
创建测试脚本 `tests/integration/test_full_pipeline.py`：

```python
def test_complete_pipeline():
    """测试完整分析流程"""
    # 1. 创建测试数据
    adata = create_test_data()

    # 2. QC
    from scLucid.qc import run_standard_qc, QCWorkflowConfig
    qc_config = QCWorkflowConfig(...)
    adata = run_standard_qc(adata, qc_config)
    assert 'leiden' not in adata.obs  # QC后没有聚类

    # 3. Preprocessing
    from scLucid.preprocess import run_preprocessing
    adata = run_preprocessing(adata)
    assert 'X_pca' in adata.obsm  # 检查PCA
    assert 'X_umap' in adata.obsm  # 检查UMAP

    # 4. Analysis
    from scLucid.analysis import run_standard_analysis
    adata = run_standard_analysis(adata)
    assert 'cell_type' in adata.obs  # 检查注释

    print("✓ 完整流程测试通过")
```

#### 1.3 错误处理完善
- [ ] 添加输入验证函数
- [ ] 统一错误消息格式
- [ ] 提供清晰的错误提示

---

## 🔧 第二阶段：功能完善（2-3周）

### 2.1 QC 模块完善

#### 优先级 P0（必须）
- [ ] **基础指标完整性**
  ```python
  # 检查必需指标
  required_metrics = [
      'n_genes', 'n_counts', 'pct_counts_mt',
      'doublet_score', 'cell_cycle_phase'
  ]
  ```

- [ ] **过滤功能健壮性**
  - 边界值处理（空数据、极端值）
  - 过滤前后统计对比
  - 过滤报告生成

#### 优先级 P1（重要）
- [ ] **双细胞检测整合**
  - 统一 scrublet + heuristic 双方法结果
  - 提供置信度评分
  - 可视化双细胞分布

### 2.2 Preprocess 模块完善

#### 优先级 P0
- [ ] **HVG 选择优化**
  - 多种 HVG 选择方法的统一接口
  - HVG 稳定性评估
  - 自动推荐 HVG 数量

- [ ] **批次校正方法验证**
  - 测试所有批次校正方法（harmony, scanorama, combat, bbknn）
  - 提供方法选择建议
  - 批次效应评估

#### 优先级 P1
- [ ] **数据标准化**
  - 统一 normalization 参数
  - log1p 标准化验证
  - scaling 方法选择

### 2.3 Analysis 模块完善

#### 优先级 P0
- [ ] **聚类方法验证**
  - 测试所有聚类方法（leiden, louvain, kmeans, hdbscan）
  - 参数优化建议
  - 聚类质量评估

- [ ] **注释方法整合**
  - 统一 marker-based 和 CellTypist 结果
  - 置信度整合
  - 注释质量评分

#### 优先级 P1
- [ ] **差异表达完整性**
  - 确保 compare_groups, compare_conditions 可用
  - 结果格式统一
  - 可视化完整性

### 2.4 Plotting 模块完善

#### 优先级 P0
- [ ] **基础图表完整性**
  ```python
  # 必需图表类型
  required_plots = [
      'embedding',      # UMAP/t-SNE
      'violin',         # 表达分布
      'heatmap',        # 表达热图
      'volcano',        # 差异表达
      'dotplot'         # 点图
  ]
  ```

- [ ] **参数命名统一**
  ```python
  # 统一使用
  plot_xxx(adata, groupby='leiden', color='cell_type')
  ```

#### 优先级 P1
- [ ] **图表保存功能**
  - 支持多种格式（PNG, PDF, SVG）
  - DPI 和尺寸配置
  - 批量保存

---

## 🎨 第三阶段：高级功能与优化（3-4周）

### 3.1 工作流优化

#### 高级工作流函数
- [ ] **半自动分析模式**
  ```python
  # 交互式选择参数
  run_interactive_analysis(adata)
  ```

- [ ] **流程检查点**
  - 保存中间结果
  - 从断点恢复
  - 缓存机制

#### 配置模板
- [ ] **预设配置模板**
  ```python
  templates = {
      'pbmc_3k': {...},
      'large_dataset': {...},
      'spatial_data': {...}
  }
  ```

### 3.2 性能优化

- [ ] **并行处理**
  - 多样本并行分析
  - 多线程计算

- [ ] **内存优化**
  - 大数据集处理
  - 数据类型优化

### 3.3 用户体验

- [ ] **进度显示**
  - rich 进度条
  - 详细日志

- [ ] **报告生成**
  - HTML 分析报告
  - PDF 报告

---

## 📊 具体实施建议

### 推荐方案：敏捷迭代

**Sprint 1（本周）**：核心流程验证
- 创建端到端测试
- 修复发现的bug
- 文档化数据流

**Sprint 2（下周）**：QC + Preprocess 完善
- 补充缺失功能
- 优化参数传递
- 添加单元测试

**Sprint 3（第三周）**：Analysis 完善
- 注释流程优化
- 聚类方法验证
- 差异表达完整性

**Sprint 4（第四周）**：Plotting + 优化
- 图表完整性检查
- 性能优化
- 用户体验改进

---

## 🔍 每个功能完善的标准

### 完善度检查清单

对于每个功能，应该满足：

#### 代码质量
- [ ] 有完整的 docstring
- [ ] 类型注解完整
- [ ] 参数验证
- [ ] 错误处理
- [ ] 日志记录

#### 功能完整性
- [ ] 输入输出明确
- [ ] 边界情况处理
- [ ] 默认值合理
- [ ] 可配置参数充分

#### 测试覆盖
- [ ] 单元测试
- [ ] 集成测试
- [ ] 边界测试
- [ ] 性能测试（如需要）

#### 文档
- [ ] API 文档
- [ ] 使用示例
- [ ] 参数说明
- [ ] 返回值说明

---

## 🎯 建议的起点

### 推荐从 **QC → Preprocess → Analysis** 顺序开始

**原因**：
1. QC 是基础，必须稳健
2. Preprocess 承上启下
3. Analysis 在稳固基础上更容易完善

### 第一个任务建议

**创建 `tests/integration/test_data_flow.py`**

```python
"""测试数据在模块间的流转"""

def test_qc_to_preprocess():
    """测试 QC 到 Preprocess 的数据流"""
    # 模拟数据
    adata = create_test_adata()

    # QC
    from scLucid.qc import run_standard_qc
    adata_qc = run_standard_qc(adata.copy())

    # 验证 QC 输出符合 Preprocess 期望
    assert 'X_pca' not in adata_qc.obsm  # QC不做PCA
    assert adata_qc.obs['pct_counts_mt'] is not None

    # Preprocessing
    from scLucid.preprocess import run_preprocessing
    adata_pp = run_preprocessing(adata_qc)

    # 验证输出
    assert 'X_pca' in adata_pp.obsm
    assert 'X_umap' in adata_pp.obsm
    assert 'leiden' in adata_pp.obs
```

---

## 💡 需要您的决策

请告诉我：

1. **您更倾向于哪种方式？**
   - A) 系统化地完善（按模块顺序）
   - B) 从核心流程开始（端到端）
   - C) 针对特定需求（您有具体功能要优先完善吗？）

2. **当前最关心的问题是什么？**
   - 功能完整性（某些功能缺失或不够完善）
   - 流程设计（各模块如何协同工作）
   - 性能优化（处理大数据集）
   - 用户体验（易用性、文档）
   - 其他？

3. **您希望我先做什么？**
   - 创建测试框架
   - 完善某个特定模块
   - 设计标准流程
   - 其他？

告诉我您的想法，我会制定详细的实施计划并开始执行！🚀
