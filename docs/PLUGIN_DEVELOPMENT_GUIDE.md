# 抽象基类插件开发 - 可视化指南

## 🤔 什么是插件开发？

### 传统方式 vs 插件方式

**传统方式（❌ 问题）**：
```python
# 想要添加自定义QC方法？必须修改核心代码
def my_custom_qc(adata, threshold):
    # ... 你的代码 ...
    pass

# 需要修改 scLucid/qc/filtering.py
# 如果出错了，可能破坏整个包
# 代码审查困难，维护成本高
```

**插件方式（✅ 解决）**：
```python
# 创建独立的插件文件，不修改核心代码
from scLucid import AnalysisStep

class MyCustomQC(AnalysisStep):
    def validate_input(self, adata):
        return True

    def run(self, adata, **kwargs):
        # ... 你的代码 ...
        return adata

# 注册到工厂
from scLucid import AnalysisStepFactory
AnalysisStepFactory.register('my_qc', MyCustomQC)

# 使用
qc = AnalysisStepFactory.create('my_qc')
qc.run(adata)
```

---

## 🎯 核心概念：抽象基类（ABC）

### 什么是抽象基类？

**抽象基类**就像一个"模板"或"合同"，定义了：
- 必须实现哪些方法
- 方法的签名是什么
- 如何与系统其他部分交互

### 可视化类比

```
📋 抽象基类 = 建筑蓝图
                         │
         ┌───────────────┼───────────────┐
         │               │               │
    具体实现1      具体实现2      具体实现3
    (MyQC)       (YourQC)      (CustomQC)
         │               │               │
         └───────────────┴───────────────┘
                         │
                    所有实现都遵循
                    相同的"蓝图"
```

---

## 🔧 scLucid 的抽象基类体系

### 6大抽象基类

```
scLucid 基础架构
│
├─ AnalysisStep          → 所有分析步骤的基础
│   ├─ QC Filter        → QC 过滤器
│   ├─ Cell Annotator   → 细胞注释器
│   ├─ Scoring Method   → 评分方法
│   ├─ Plotting Backend → 绘图后端
│   └─ Proportion Method → 比例分析方法
│
└─ AnalysisStepFactory  → 工厂模式
    ├─ register()        → 注册插件
    ├─ create()          → 创建实例
    └─ list_steps()      → 列出所有插件
```

---

## 💡 实际使用场景

### 场景1：创建自定义QC过滤器

```python
from scLucid import AnalysisStep, SclucidBaseConfig
from anndata import AnnData
from pydantic import Field

# 步骤1: 定义配置
class StrictQCConfig(SclucidBaseConfig):
    min_genes: int = Field(default=500, ge=0)
    max_mt: float = Field(default=10.0, ge=0, le=100)

# 步骤2: 创建插件类（继承 AnalysisStep）
class StrictQCFilter(AnalysisStep):
    def __init__(self, config: StrictQCConfig):
        self.config = config

    def validate_input(self, adata: AnnData) -> bool:
        """检查输入是否有效"""
        if adata.n_obs == 0:
            raise ValueError("没有细胞")
        return True

    def run(self, adata: AnnData, **kwargs) -> AnnData:
        """执行QC过滤"""
        # 你的自定义逻辑
        import scanpy as sc
        sc.pp.filter_cells(adata, min_genes=self.config.min_genes)
        return adata

    def get_summary(self) -> dict:
        """返回摘要"""
        return {"min_genes": self.config.min_genes}

# 步骤3: 注册插件
from scLucid import AnalysisStepFactory
AnalysisStepFactory.register('strict_qc', StrictQCFilter)

# 步骤4: 使用插件
qc = AnalysisStepFactory.create('strict_qc')
qc.run(adata)
```

### 场景2：组合多个插件

```python
# 定义工作流
steps = [
    'strict_qc',           # 你的自定义QC
    'clustering',          # 标准聚类
    'my_annotator'         # 你的自定义注释器
]

# 配置每一步
configs = {
    'strict_qc': {'config': StrictQCConfig(min_genes=500)},
    'clustering': {'resolution': 0.8},
    'my_annotator': {...}
}

# 运行自定义工作流
from scLucid.analysis import run_custom_analysis
adata = run_custom_analysis(adata, steps=steps, step_configs=configs)
```

---

## 🏗️ 工厂模式（Factory Pattern）

### 为什么需要工厂？

```
直接创建：
    qc1 = StrictQCFilter(config1)
    qc2 = AnotherQC(config2)
    ❌ 需要知道所有具体的类名

工厂创建：
    qc1 = AnalysisStepFactory.create('strict_qc', config=config1)
    qc2 = AnalysisStepFactory.create('another_qc', config=config2)
    ✅ 只需要知道插件名称，更灵活
```

### 工厂的优势

1. **解耦**: 使用者不需要知道具体类名
2. **动态**: 运行时决定创建哪个插件
3. **可扩展**: 添加新插件不需要修改使用代码
4. **可测试**: 容易mock和测试

---

## 📦 插件开发的好处

### 对开发者

| 优势 | 说明 |
|------|------|
| **独立开发** | 在自己的文件中开发，不影响核心代码 |
| **快速迭代** | 不需要等待代码审查和合并 |
| **版本控制** | 可以有自己的版本号和发布周期 |
| **选择性使用** | 其他用户可以选择是否使用你的插件 |

### 对用户

| 优势 | 说明 |
|------|------|
| **更多选择** | 可以使用社区开发的多种算法 |
| **定制化** | 可以为自己实验室定制特定分析流程 |
| **实验性** | 可以尝试实验性功能而不影响稳定性 |

### 对维护者

| 优势 | 说明 |
|------|------|
| **核心稳定** | 核心代码修改减少，更稳定 |
| **降低负担** | 不需要维护所有可能的算法变体 |
| **社区贡献** | 社区可以贡献和共享插件 |

---

## 🎓 实际例子

### 例子1：数据库注释器

```python
class DatabaseAnnotator(AnalysisStep):
    """从外部数据库注释细胞类型"""

    def validate_input(self, adata):
        return 'X_pca' in adata.obsm

    def run(self, adata, database_path, **kwargs):
        # 连接数据库
        # 查询最相似的细胞类型
        # 返回注释结果
        return adata

# 注册
AnalysisStepFactory.register('db_annotator', DatabaseAnnotator)
```

### 例子2：机器学习评分器

```python
class MLScorer(AnalysisStep):
    """使用机器学习模型评分"""

    def validate_input(self, adata):
        return adata.shape[1] > 0

    def run(self, adata, model_path, **kwargs):
        # 加载训练好的模型
        # 对每个细胞评分
        # 存储在 adata.obs
        return adata

# 注册
AnalysisStepFactory.register('ml_scorer', MLScorer)
```

### 例子3：可视化后端

```python
class InteractivePlotter(AnalysisStep):
    """创建交互式可视化"""

    def validate_input(self, adata):
        return adata.n_obs > 0

    def run(self, adata, **kwargs):
        # 创建 plotly 图表
        # 添加交互控件
        # 返回 HTML
        return adata

    def get_summary(self):
        return {"type": "interactive"}

# 注册
AnalysisStepFactory.register('interactive_plot', InteractivePlotter)
```

---

## 📚 学习资源

### 在项目内
- `src/scLucid/base_interfaces.py` - 抽象基类定义
- `examples/plugin_development_example.py` - 完整示例代码
- `docs/NAMING_CONVENTIONS.md` - 命名规范

### 外部资源
- Python ABC (Abstract Base Classes): https://docs.python.org/3/library/abc.html
- Factory Pattern: https://refactoring.guru/design-patterns/factory-method
- Plugin Architecture: https://en.wikipedia.org/wiki/Plugin_(computing)

---

## ❓ 常见问题

### Q1: 我必须使用插件吗？
**A**: 不必须！插件是可选的。你可以继续使用现有的函数和类。插件只是提供了一个扩展的选项。

### Q2: 插件会影响性能吗？
**A**: 插件是普通的Python类，性能与核心函数相当。工厂创建的开销极小（<1ms）。

### Q3: 如何分享我的插件？
**A**: 你可以：
- 发布到 PyPI 作为独立包
- 提供代码片段给用户
- 提交 PR 到 scLucid 仓库（如果通用性强）

### Q4: 插件会破坏核心代码吗？
**A**: 不会。插件是独立运行的，不会修改核心代码。

### Q5: 我可以同时使用核心功能和插件吗？
**A**: 可以！你可以在同一个分析流程中混合使用核心功能和自定义插件。

---

## 🚀 快速开始

1. **学习抽象基类**: 查看 `src/scLucid/base_interfaces.py`
2. **运行示例**: `python examples/plugin_development_example.py`
3. **创建第一个插件**: 继承 `AnalysisStep`，实现必需方法
4. **注册和使用**: 使用 `AnalysisStepFactory`
5. **分享**: 如果有用，考虑分享给社区

祝你开发出强大的插件！🎉
