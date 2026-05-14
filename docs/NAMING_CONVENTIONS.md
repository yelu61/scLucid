# scLucid 命名规范

本文档定义了 scLucid 项目的统一命名规范，确保代码库的一致性和可读性。

## 函数命名规范

### 公共 API 函数

所有公共函数应使用描述性名称，遵循以下模式：

| 前缀 | 用途 | 示例 |
|------|------|------|
| `run_*` | 运行完整的工作流或分析流程 | `run_preprocessing()`, `run_annotation()`, `analyze_celltype_proportion()` |
| `calculate_*` | 计算指标或分数 | `calculate_qc_metrics()`, `calculate_signature_matrix()` |
| `find_*` | 查找或识别特征 | `find_markers()`, `find_hvgs()` |
| `get_*` | 获取数据或配置 | `get_marker_manager()`, `get_summary()` |
| `plot_*` | 绘制图表 | `plot_embedding()`, `plot_volcano()` |
| `score_*` | 评分计算 | `score_cell_types()`, `score_by_gene_sets()` |
| `compare_*` | 比较分析 | `compare_groups()`, `compare_conditions()` |
| `filter_*` | 过滤操作 | `filter_cells()`, `filter_markers()` |
| `predict_*` | 预测操作 | `predict_doublets()` |
| `suggest_*` | 建议参数 | `suggest_qc_thresholds()`, `suggest_hvg_choice()` |

### 私有函数

私有函数（内部使用）应以单下划线 `_` 开头：

```python
def _validate_input(adata):
    """内部验证函数"""
    pass

def _calculate_metric(values):
    """内部计算函数"""
    pass
```

## 类命名规范

### 配置类

所有配置类应以 `Config` 结尾：

```python
class QCWorkflowConfig(SclucidBaseConfig):
    """QC 工作流配置"""
    pass

class ClusteringConfig(SclucidBaseConfig):
    """聚类配置"""
    pass
```

### 管理器类

管理器类应以 `Manager` 结尾：

```python
class MarkerManager:
    """Marker 基因管理器"""
    pass

class CacheManager:
    """缓存管理器"""
    pass
```

### 分析器类

分析器、预测器等应以描述性名称命名：

```python
class CellAnnotator:
    """细胞类型注释器"""
    pass

class DoubletPredictor:
    """双细胞预测器"""
    pass
```

### 抽象基类

抽象基类应以描述性名称命名，不使用特殊前缀：

```python
class AnalysisStep(ABC):
    """分析步骤抽象基类"""
    pass

class QCFilter(ABC):
    """QC 过滤器抽象基类"""
    pass
```

## 变量命名规范

### 常量

常量应使用全大写字母和下划线：

```python
DEFAULT_RESOLUTION = 0.8
MAX_MARKERS = 100
```

### 普通变量

普通变量应使用小写字母和下划线：

```python
n_cells = 1000
cluster_labels = adata.obs['leiden']
```

### 布尔变量

布尔变量应以 `is_`, `has_`, `can_` 等前缀开头：

```python
is_normalized = True
has_batch_effect = False
can_parallelize = True
```

## 模块和包命名规范

### 模块文件

模块文件应使用小写字母和下划线：

```
qc/filtering.py
preprocess/normalize.py
analysis/clustering.py
```

### 子包

子包应使用小写字母，避免使用 `_v2` 等版本标识符：

```
analysis/proportion/        # ✓ 良好
analysis/differential_expression/  # ✓ 良好
qc/workflow_v2.py           # ✗ 避免
```

## 最佳实践

### 1. 使用动词开头

函数名应以动词开头，清晰表达操作：

```python
# ✓ 良好
find_markers(adata)
filter_cells(adata)
plot_volcano(results)

# ✗ 避免
markers(adata)
cells_filter(adata)
volcano_plot(results)
```

### 2. 保持简洁但描述性

函数名应简洁但足够描述：

```python
# ✓ 良好
run_standard_qc(adata)

# ✗ 太长
run_standard_quality_control_workflow(adata)

# ✗ 太短
qc(adata)
```

### 3. 避免缩写

避免使用不明确的缩写：

```python
# ✓ 良好
calculate_n_genes(adata)
plot_umap(adata)

# ✗ 避免缩写
calc_n_genes(adata)
plt_umap(adata)
```

### 4. 统一术语

在整个代码库中使用一致的术语：

| 概念 | 统一术语 | 避免 |
|------|----------|------|
| 细胞类型 | `cell_type` | `celltype`, `CellType`, `celltype` |
| 样本ID | `sample_id` | `sample`, `SampleID`, `sampleID` |
| 标记基因 | `markers` | `marker_genes`, `marker` |
| 双细胞 | `doublet` | `doublets`（复数） |

## 示例

### 符合规范的代码示例

```python
from scLucid.base_interfaces import AnalysisStep
from scLucid.base_config import SclucidBaseConfig

class MyAnalysisConfig(SclucidBaseConfig):
    """自定义分析配置"""
    threshold: float = 0.5
    max_iter: int = 100

class MyAnalyzer(AnalysisStep):
    """自定义分析器"""

    def __init__(self, config: MyAnalysisConfig):
        self.config = config

    def validate_input(self, adata):
        """验证输入"""
        if adata.n_obs == 0:
            raise ValueError("No cells in data")
        return True

    def run(self, adata, **kwargs):
        """运行分析"""
        # 实现分析逻辑
        return adata

    def get_summary(self):
        """获取结果摘要"""
        return {"config": self.config.model_dump()}

# 公共API函数
def run_my_analysis(adata, config=None):
    """运行自定义分析"""
    if config is None:
        config = MyAnalysisConfig()

    analyzer = MyAnalyzer(config)
    analyzer.validate_input(adata)
    return analyzer.run(adata)
```

## 检查清单

在提交代码前，请确认：

- [ ] 公共函数使用规范的前缀（run_, calculate_, find_, get_, plot_ 等）
- [ ] 私有函数以 `_` 开头
- [ ] 类名使用大写字母和下划线分隔
- [ ] 配置类以 `Config` 结尾
- [ ] 管理器类以 `Manager` 结尾
- [ ] 常量使用全大写字母
- [ ] 布尔变量以 `is_`, `has_`, `can_` 开头
- [ ] 避免使用不明确的缩写
- [ ] 使用一致的术语（cell_type, sample_id, markers）
- [ ] 模块文件使用小写字母和下划线
