# scLucid 核心价值定位与实施计划

## 🎯 工具定位

### 不同于其他工具
```
传统工具（Seurat, Scanpy）:
  └─ 静态阈值：n_genes > 200（硬编码）
  └─ 固定流程：按步骤调用函数
  └─ 通用工具：适用于所有数据

scLucid（您的目标）:
  └─ 智能阈值：基于数据分布自动推荐
  └─ 自适应流程：根据数据特征调整策略
  └─ 证据驱动：每个决策都有统计/生物学依据
  └─ 肿瘤优化：专门针对肿瘤分析场景
```

### 核心创新点（论文卖点）

1. **自适应QC (Adaptive QC)**
   - 数据驱动的阈值推荐
   - 考虑批次、组织类型、数据质量
   - 提供置信度和解释

2. **智能注释策略**
   - 基于肿瘤纯度调整注释
   - Marker互斥分析
   - 多证据整合

3. **探索性分析友好**
   - 保留多种可能性的中间结果
   - 可追溯的决策链
   - 交互式参数调整

---

## 📋 系统化实施计划

### 🏗️ Phase 1: QC模块智能化（Week 1-2）

#### 目标
从"硬阈值"转向"数据驱动的智能QC"

#### 1.1 自适应阈值系统

**已有基础**：
- `qc/adaptive_threshold.py` - GMM, KDE, DBSCAN阈值学习
- `qc/doublet.py` - heuristic marker co-expression

**需要完善**：

```python
class IntelligentQCRecommender:
    """
    智能QC推荐系统

    功能：
    1. 分析数据分布（而非固定阈值）
    2. 考虑批次效应
    3. 提供置信区间
    4. 可视化建议依据
    """

    def recommend_thresholds(
        self,
        adata: AnnData,
        tissue_type: str = "unknown",
        strategy: str = "auto"
    ) -> QCThresholdRecommendation:
        """
        推荐QC阈值

        Returns:
        --------
        QCThresholdRecommendation with:
        - min_genes: 推荐值 + 95% CI
        - max_mt_percent: 推荐值 + 依据
        - doublet_threshold: 自适应推荐
        - confidence: 整体置信度评分
        - evidence: 支持证据（统计检验、图表）
        """
```

**实施步骤**：
1. 创建 `qc/intelligent_qc.py`
2. 整合现有的 `adaptive_threshold.py`
3. 添加肿瘤特异性规则：
   - 肿瘤样本的细胞周期异常
   - 肿瘤纯度低的预期
   - 双细胞（肿瘤-正常混合）的特殊处理

#### 1.2 质量评估与置信度

```python
class DataQualityAssessment:
    """
    数据质量评估

    评估维度：
    - 技术重复性（相关性）
    - 线粒体污染程度
    - 细胞周期分布
    - 批次效应强度
    - 双细胞污染估计
    """

    def assess(self, adata: AnnData) -> QualityReport:
        """
        返回质量报告：
        - overall_score: 0-100
        - quality_flags: 质量问题列表
        - recommendations: 改进建议
        - confidence: 质量评估的置信度
        """
```

---

### 🔬 Phase 2: Preprocess 模块优化（Week 3-4）

#### 目标
确保预处理步骤不丢失重要的生物学信号，特别是肿瘤相关信号

#### 2.1 HVG 选择优化

**创新点**：基于生物学知识的 HVG 选择

```python
class BiologyAwareHVGSelector:
    """
    生物学感知的HVG选择

    特点：
    - 不仅考虑方差，还考虑肿瘤相关基因
    - 保留低表达但关键的肿瘤标志物
    - 细胞类型特异性 vs 肿瘤通用性权衡
    """

    def select_hvgs(
        self,
        adata: AnnData,
        mode: str = "auto",  # tumor_aware, standard, conserved
        cancer_type: str = None
    ) -> List[str]:
        """
        选择HVGs

        tumor_aware模式：
        1. 标准方差筛选
        2. 添加已知癌症基因（如 KRAS, TP53, EGFR）
        3. 添加低表达但关键的标志物
        4. 评估HVG在不同细胞类型中的表达模式
        """
```

#### 2.2 批次效应智能处理

```python
class AdaptiveBatchCorrection:
    """
    自适应批次校正

    特点：
    - 自动评估是否需要批次校正
    - 选择合适的方法（harmony vs BBKNN vs scanorama）
    - 保留肿瘤-正常差异（不过度校正）
    """

    def should_correct(
        self,
        adata: AnnData,
        cancer_purity_threshold: float = 0.3
    ) -> BatchCorrectionRecommendation:
        """
        决策树：
        1. 评估批次混合度
        2. 评估肿瘤纯度分布
        3. 如果肿瘤纯度变化大，谨慎校正
        4. 推荐校正方法（或不校正）
        """
```

---

### 🧬 Phase 3: Analysis 模块 - 肿瘤特异性（Week 5-7）

#### 这是核心创新部分！

#### 3.1 肿瘤纯度感知的注释

```python
class CancerPurityAwareAnnotation:
    """
    肿瘤纯度感知的注释系统

    核心思想：
    - 高纯度肿瘤样本 → 严格注释
    - 低纯度/混合样本 → 混合细胞类型标签
    - 根据纯度自动调整注释策略
    """

    def annotate_with_purity(
        self,
        adata: AnnData,
        reference: AnnData = None,
        method: str = "adaptive"
    ) -> AnnData:
        """
        自适应注释流程

        步骤：
        1. 估计每个cluster的肿瘤纯度
        2. 高纯度cluster (>90%):
           - 使用严格的marker匹配
           - 细分亚型
        3. 低纯度cluster (<70%):
           - 允许混合标签
           - 标注为混合细胞类型
        4. 中等纯度:
           - 同时给出主要和次要细胞类型
        """

    def estimate_cancer_purity(
        self,
        adata: AnnData,
        cluster: str,
        cancer_markers: List[str]
    ) -> float:
        """
        估计肿瘤纯度

        方法：
        - 基于已知的癌症标志物
        - 考虑 Marker co-expression patterns
        - 使用infercnv结果（如果有）
        - 估计置信区间
        """
```

#### 3.2 Marker互斥分析系统

```python
class MarkerMutualExclusivityAnalyzer:
    """
    Marker互斥分析

    目的：
    - 检测生物学上不应该共表达的marker
    - 识别双细胞或异常状态
    - 指导注释策略调整

    应用：
    - 上皮 + 免疫marker共存 → 可能是双细胞或CMT（癌肉瘤样）
    - 多个lineage marker共存 → 未分化/干性
    - 肿瘤marker + 正常marker → 肿瘤微环境
    """

    def analyze_mutual_exclusivity(
        self,
        adata: AnnData,
        marker_sets: Dict[str, List[str]]
    ) -> MutualExclusivityReport:
        """
        分析marker共表达模式

        返回：
        - 互斥性评分（0-1，0=完全互斥，1=高共表达）
        - 异常cluster列表
        - 生物学解释
        - 建议的处理策略
        """
```

#### 3.3 自适应注释策略

```python
class AdaptiveAnnotationStrategy:
    """
    自适应注释策略

    根据数据特征自动选择注释方法：

    数据特征 → 策略
    ─────────────────────
    高纯度，清晰marker  → marker_based（严格）
    低纯度，混合      → semi_supervised（聚类）
    批次效应强         → batch_correction + annotation
    样本量小（<5）    → transfer_learning
    样本量大（>10）   → supervised（CellTypist）
    有参考 atlas     → label_transfer
    """

    def recommend_strategy(
        self,
        adata: AnnData,
        context: AnalysisContext
    ) -> AnnotationStrategy:
        """
        推荐最佳注释策略

        返回：
        - recommended_method: CellTypist vs marker-based vs ensemble
        - parameters: 推荐的参数
        - rationale: 推荐理由
        - expected_performance: 预期性能
        """
```

---

### 📊 Phase 4: 证据驱动与可追溯性（Week 8-9）

#### 4.1 决策记录系统

```python
from datetime import datetime
from typing import Any, Dict, List

class AnalysisDecisionLogger:
    """
    分析决策记录系统

    记录每个重要决策的：
    - 输入数据特征
    - 决策依据（统计检验、可视化）
    - 选择的参数
    - 决策结果
    - 置信度
    """

    def log_decision(
        self,
        decision_type: str,  # "threshold_selection", "annotation_method", etc.
        inputs: Dict[str, Any],
        evidence: List[Evidence],
        decision: Any,
        confidence: float
    ):
        """
        记录决策

        Evidence类型：
        - StatisticalTest: p-value, effect size
        - Visualization: plot path
        - PriorKnowledge: 文献引用
        - DataPattern: 数据分布描述
        """

    def generate_report(self) -> DecisionReport:
        """
        生成决策报告

        包含：
        - 决策树（为什么选择这个参数）
        - 证据链
        - 可视化支持
        - 可重现性信息
        """
```

#### 4.2 置信度评分系统

```python
class ConfidenceScorer:
    """
    置信度评分系统

    为每个分析结果提供多维度的置信度：

    1. Statistical Confidence
       - p-value, effect size
       - 样本量是否足够

    2. Technical Confidence
       - 数据质量
       - 批次效应控制
       - 方法的适用性

    3. Biological Confidence
       - 已知marker支持
       - 文献一致性
       - 生物学合理性

    4. Reproducibility Confidence
       - 随机种子影响
       - 参数敏感性
       - 方法稳定性
    """

    def score_confidence(
        self,
        result: AnalysisResult
    ) -> ConfidenceScore:
        """
        综合置信度评分

        返回：
        - overall_score: 0-100
        - dimension_scores: 各维度得分
        - concerns: 潜在问题列表
        - recommendations: 如何提高置信度
        """
```

---

### 🎨 Phase 5: 肿瘤特异性功能（Week 10-12）

#### 5.1 肿瘤微环境分析

```python
class TumorMicroenvironmentAnalyzer:
    """
    肿瘤微环境分析

    分析维度：
    1. 免疫浸润评估
    - T细胞亚群比例
    - 耗竭vs效应比例
    - 免疫检查点表达

    2. 基质细胞
    - 成纤维细胞激活状态
    - 血管生成相关
    - 细胞外基质重塑

    3. 细胞-细胞通讯
    - 肿瘤-免疫相互作用
    - 配体-受体对分析
    """

    def analyze_tme(
        self,
        adata: AnnData,
        cancer_type: str = None
    ) -> TMEReport:
        """
        分析肿瘤微环境

        返回：
        - immune_infiltration_score
        - stromal_activation_score
        - cell_cell_communication
        - therapeutic_targets
        """
```

#### 5.2 CNV与克隆结构

```python
class ClonalArchitectureAnalyzer:
    """
    克隆结构分析（基于infercnv）

    分析：
    1. 染色体异常模式
    - 大规模gain/loss
    - 局部扩增（如HER2, EGFR）
    - 染色体不稳定性

    2. 推断克隆结构
    - 主克隆 vs 亚克隆
    - 进化关系

    3. 治疗相关标记
    - 靶向标志
    - 耐药标志
    """

    def infer_clonal_architecture(
        self,
        adata: AnnData
    ) -> ClonalArchitecture:
        """
        推断克隆架构

        返回：
        - clonal_expansions: 扩增事件
        - subclonal_mutations: 亚克隆
        - therapeutic_targets
        """
```

#### 5.3 治疗响应预测

```python
class TherapyResponsePredictor:
    """
    治疗响应预测

    预测：
    - 免疫检查点阻断响应
    - 靶向治疗响应
    - 化疗敏感性
    - 耐药可能标志物

    基于：
    - 基因表达特征
    - 免疫浸润模式
    - TMB（Tumor Mutational Burden）
    - 既往文献
    """
```

---

## 📝 Phase 6: 文档与重现性（Week 13-14）

### 6.1 智能分析报告生成

```python
class IntelligentAnalysisReportGenerator:
    """
    智能分析报告生成器

    生成的报告包含：

    1. 数据质量总结
    - 质量评分
    - 技术指标
    - 可视化

    2. 分析流程追溯
    - 每个步骤的决策依据
    - 参数选择理由
    - 可视化支持

    3. 结果总结
    - 主要发现
    - 置信度评估
    - 局限性讨论

    4. 肿瘤特异性发现
    - 肿瘤纯度估计
    - 微环境特征
    - 治疗相关标志物

    5. 方法学创新
    - 使用的智能算法
    - 与现有工具的对比
    - 可重现性说明
    """

    def generate_report(
        self,
        adata: AnnData,
        analysis_history: AnalysisHistory
    ) -> Report:
        """
        生成HTML/PDF报告
        """
```

---

## 🎯 论文核心价值主张

### 与 Seurat/Scanpy 的对比

| 维度 | Seurat/Scanpy | scLucid（您的工具） |
|------|---------------|-------------------|
| **阈值选择** | 固定阈值（用户设定） | 数据驱动推荐 + 置信区间 |
| **流程适应性** | 固定流程 | 自适应流程选择 |
| **注释策略** | 单一方法 | 纯度感知 + 多策略整合 |
| **肿瘤支持** | 通用工具 | 肿瘤特异性功能 |
| **决策依据** | 隐式 | 显式证据链 |
| **可重现性** | 参数依赖 | 证据驱动，可追溯 |

### 论文可以发表的要点

1. **方法学创新**
   - 自适应阈值学习系统
   - 肿瘤纯度感知注释
   - Marker互斥分析算法
   - 证据驱动决策框架

2. **肿瘤特异性应用**
   - 肿瘤微环境分析
   - 克隆结构推断
   - 治疗响应预测
   - 混合样本处理

3. **探索性分析支持**
   - 保留多种可能性的中间结果
   - 交互式决策
   - 可追溯的决策链

---

## 📅 实施时间表

### Month 1-2: 核心智能算法
- Week 1-2: QC 智能化
- Week 3-4: Preprocess 优化
- Week 5-7: Analysis 肿瘤特异性 ⭐
- Week 8-9: 证据驱动系统

### Month 3: 集成与测试
- Week 10-11: 肿瘤特异性功能
- Week 12: 文档和示例
- Week 13-14: 完整流程测试

### Month 4-6: 论文撰写与发表
- 数据收集和案例研究
- 方法论论文撰写
- 配套分析工具
- 投稿到 Nature Methods / Genome Biology / Cell Systems

---

## 🚀 第一步：从哪里开始？

**本周开始任务**：创建 `qc/intelligent_qc.py`

这个文件将整合现有的 `adaptive_threshold.py`，创建一个智能QC推荐系统。

需要开始吗？我可以立即开始实现这个核心创新功能！🎯
