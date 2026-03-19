# Data文件夹使用指南

## 📁 数据集概览

`data/` 文件夹包含3个精心挑选的数据集，用于：

✅ **测试** (tests/)
✅ **示例** (examples/)
✅ **文档** (docs/notebooks/)

---

## 📊 数据集详情

### 1. PBMC3K (正常组织)
- **路径**: `data/pbmc3k/`
- **文件**: `pbmc3k_raw.h5ad`
- **来源**: 10x Genomics PBMC datasets
- **类型**: 外周血单核细胞（正常组织）
- **细胞数**: ~2,700
- **基因数**: ~2,000
- **特点**:
  - ✅ 单批次、同组织
  - ✅ 高质量数据
  - ✅ 适合作为基线

**适用场景**：
- 测试统一阈值策略
- 基础功能测试
- 快速示例（小数据集）

---

### 2. LUAD (肺腺癌 - 肿瘤组织)
- **路径**: `data/human_LUAD_GSE131907/`
- **来源**: GEO GSE131907
- **类型**: 肺腺癌（肿瘤组织）
- **细胞数**: 待确认
- **特点**:
  - ✅ 肿瘤组织
  - ✅ 肿瘤微环境
  - ✅ 更高的线粒体含量
  - ✅ 肿瘤-正常混合

**适用场景**：
- 测试肿瘤感知策略
- 验证intelligent_qc的肿瘤特异性
- 肿瘤研究示例

---

### 3. 黑色素瘤 (多批次、肿瘤)
- **路径**: `data/mouse_melanoma_GSE119352/`
- **来源**: GEO GSE119352
- **类型**: 小鼠黑色素瘤（肿瘤组织）
- **细胞数**: 待确认
- **特点**:
  - ✅ 多批次
  - ✅ 肿瘤组织
  - ✅ 极端异质性
  - ✅ 批次效应明显

**适用场景**：
- 测试混合策略
- 测试样本特异性策略
- 批次校正示例
- 极端场景测试

---

## 🎯 使用场景矩阵

| 数据集 | 正常/肿瘤 | 批次数 | 主要用途 |
|--------|----------|--------|---------|
| **PBMC3K** | 正常 | 1 | 基线测试、快速示例 |
| **LUAD** | 肿瘤 | 待确认 | 肿瘤感知策略测试 |
| **黑色素瘤** | 肿瘤 | 多 | 混合策略、批次校正 |

---

## 📝 使用示例

### 1. 在测试中使用

#### tests/qc/test_intelligent_qc.py
```python
import pytest
from pathlib import Path
from scLucid.qc import recommend_intelligent_qc

# 数据路径
DATA_DIR = Path(__file__).parent.parent.parent / "data"

@pytest.fixture
def pbmc_data():
    """PBMC数据（正常组织）"""
    import scanpy as sc
    adata = sc.read_h5ad(DATA_DIR / "pbmc3k" / "pbmc3k_raw.h5ad")
    return adata

@pytest.fixture
def luad_data():
    """LUAD数据（肿瘤组织）"""
    import scanpy as sc
    adata = sc.read_h5ad(DATA_DIR / "human_LUAD_GSE131907" / "luad.h5ad")
    return adata

def test_intelligent_qc_normal_tissue(pbmc_data):
    """测试正常组织的智能QC"""
    rec = recommend_intelligent_qc(pbmc_data, tissue_type="normal")

    # 正常组织应该使用更严格的MT阈值
    assert rec.max_mt_percent.threshold < 20.0
    assert rec.overall_strategy in ['standard', 'conservative']

def test_intelligent_qc_tumor_tissue(luad_data):
    """测试肿瘤组织的智能QC"""
    rec = recommend_intelligent_qc(luad_data, tissue_type="lung_tumor")

    # 肿瘤组织应该有更高的MT阈值
    assert rec.max_mt_percent.threshold > 15.0  # 允许更高的MT
    assert rec.overall_strategy in ['tumor_aware', 'auto']
```

#### tests/integration/test_strategies_comparison.py
```python
"""对比三种QC策略"""
from scLucid.qc import recommend_qc_strategy, recommend_intelligent_qc

def test_strategies_comparison():
    """使用三个数据集对比策略"""
    import scanpy as sc
    from pathlib import Path

    DATA_DIR = Path(__file__).parent.parent.parent / "data"

    # 加载三个数据集
    pbmc = sc.read_h5ad(DATA_DIR / "pbmc3k" / "pbmc3k_raw.h5ad")
    luad = sc.read_h5ad(DATA_DIR / "human_LUAD_GSE131907" / "luad.h5ad")
    melanoma = sc.read_h5ad(DATA_DIR / "mouse_melanoma_GSE119352" / "melanoma.h5ad")

    # 添加组织类型标签
    pbmc.obs['tissue_type'] = 'normal'
    luad.obs['tissue_type'] = 'lung_tumor'
    melanoma.obs['tissue_type'] = 'melanoma'

    samples = {
        'pbmc': pbmc,
        'luad': luad,
        'melanoma': melanoma
    }

    # 获取策略推荐
    strategy, rationale = recommend_qc_strategy(
        samples,
        tissue_key='tissue_type'
    )

    # 验证推荐
    assert strategy in ['unified', 'sample_specific', 'hybrid']
```

---

### 2. 在示例中使用

#### examples/compare_qc_strategies.py
```python
"""
对比三种QC策略的完整示例

使用 data/ 文件夹中的三个数据集演示：
1. 统一阈值策略
2. 样本特异性策略
3. 混合策略（推荐）
"""

import scanpy as sc
from pathlib import Path
from scLucid.qc import (
    recommend_intelligent_qc,
    recommend_qc_strategy,
    calculate_qc_metrics
)

# 数据路径
DATA_DIR = Path(__file__).parent.parent / "data"

# 加载数据
print("加载数据集...")
pbmc = sc.read_h5ad(DATA_DIR / "pbmc3k" / "pbmc3k_raw.h5ad")
luad = sc.read_h5ad(DATA_DIR / "human_LUAD_GSE131907" / "luad.h5ad")
melanoma = sc.read_h5ad(DATA_DIR / "mouse_melanoma_GSE119352" / "melanoma.h5ad")

# 计算QC指标
pbmc = calculate_qc_metrics(pbmc)
luad = calculate_qc_metrics(luad)
melanoma = calculate_qc_metrics(melanoma)

# 添加组织类型
pbmc.obs['tissue_type'] = 'normal'
luad.obs['tissue_type'] = 'lung_tumor'
melanoma.obs['tissue_type'] = 'melanoma'

samples = {
    'PBMC (正常)': pbmc,
    'LUAD (肺腺癌)': luad,
    '黑色素瘤 (多批次)': melanoma
}

# ============================================
# 1. 自动策略推荐
# ============================================
print("\n" + "=" * 70)
print("步骤1: 自动策略推荐")
print("=" * 70)

strategy, rationale = recommend_qc_strategy(
    samples,
    tissue_key='tissue_type'
)

for line in rationale:
    print(line)

# ============================================
# 2. 智能QC推荐
# ============================================
print("\n" + "=" * 70)
print("步骤2: 获取智能QC推荐")
print("=" * 70)

recommendations = {}
for sample_name, adata in samples.items():
    tissue_type = adata.obs['tissue_type'].iloc[0]

    rec = recommend_intelligent_qc(
        adata,
        tissue_type=tissue_type,
        plot=False
    )

    recommendations[sample_name] = rec

    print(f"\n{sample_name}:")
    print(f"  策略: {rec.overall_strategy.value}")
    print(f"  min_genes: {rec.min_genes.threshold} "
          f"[95% CI: {rec.min_genes.ci_lower}-{rec.min_genes.ci_upper}]")
    print(f"  max_mt_percent: {rec.max_mt_percent.threshold:.1f}% "
          f"[95% CI: {rec.max_mt_percent.ci_lower:.1f}-"
          f"{rec.max_mt_percent.ci_upper:.1f}]")
    print(f"  数据质量: {rec.data_quality_score:.1f}/100")

# ============================================
# 3. 应用混合策略
# ============================================
print("\n" + "=" * 70)
print("步骤3: 应用混合策略")
print("=" * 70)

import numpy as np
from scipy.stats import median_abs_deviation

# 计算全局约束
thresholds = [rec.min_genes.threshold for rec in recommendations.values()]
global_median = np.median(thresholds)
global_mad = median_abs_deviation(thresholds)
lower_bound = global_median - 3 * global_mad
upper_bound = global_median + 3 * global_mad

print(f"\n全局约束:")
print(f"  中位数: {global_median:.0f}")
print(f"  MAD: {global_mad:.0f}")
print(f"  范围: [{lower_bound:.0f}, {upper_bound:.0f}]")

# 应用混合策略
filtered_samples = {}
for sample_name, adata in samples.items():
    rec = recommendations[sample_name]
    threshold = np.clip(rec.min_genes.threshold, lower_bound, upper_bound)

    adata_filtered = adata[
        (adata.obs['n_genes'] > threshold) &
        (adata.obs['pct_counts_mt'] < rec.max_mt_percent.threshold)
    ].copy()

    filtered_samples[sample_name] = adata_filtered

    adjusted = " ✓" if (rec.min_genes.threshold != threshold) else ""
    print(f"\n{sample_name}:")
    print(f"  原始: {len(adata)} 细胞")
    print(f"  过滤后: {len(adata_filtered)} 细胞")
    print(f"  保留率: {len(adata_filtered)/len(adata):.1%}")
    print(f"  阈值调整: {rec.min_genes.threshold:.0f} → {threshold:.0f}{adjusted}")

print("\n" + "=" * 70)
print("✓ 混合策略应用完成")
print("=" * 70)
```

---

### 3. 在文档中使用

#### docs/notebooks/05_intelligent_qc.ipynb
```python
# notebook cell 1
{
    "cell_type": "markdown",
    "metadata": {},
    "source": [
        "# 智能QC示例\n",
        "\n",
        "本notebook演示如何使用scLucid的智能QC功能。\n",
        "\n",
        "## 数据集\n",
        "\n",
        "我们使用三个真实数据集：\n",
        "1. **PBMC3K**: 正常组织，单批次\n",
        "2. **LUAD**: 肺腺癌，肿瘤组织\n",
        "3. **黑色素瘤**: 多批次，肿瘤组织\n",
        "\n",
        "这些数据集位于 `data/` 文件夹中。"
    ]
}

# notebook cell 2
{
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": [
        "import scanpy as sc\n",
        "from pathlib import Path\n",
        "from scLucid.qc import (\n",
        "    recommend_intelligent_qc,\n",
        "    recommend_qc_strategy,\n",
        "    calculate_qc_metrics\n",
        ")\n",
        "\n",
        "# 数据路径\n",
        "DATA_DIR = Path('../../data')\n",
        "\n",
        "# 加载PBMC数据\n",
        "pbmc = sc.read_h5ad(DATA_DIR / 'pbmc3k' / 'pbmc3k_raw.h5ad')\n",
        "pbmc = calculate_qc_metrics(pbmc)\n",
        "\n",
        "# 获取智能推荐\n",
        "rec = recommend_intelligent_qc(pbmc, tissue_type='normal')\n",
        "\n",
        "print(f\"推荐阈值: min_genes > {rec.min_genes.threshold}\")\n",
        "print(f\"置信区间: [{rec.min_genes.ci_lower}, {rec.min_genes.ci_upper}]\")\n",
        "print(f\"数据质量: {rec.data_quality_score:.1f}/100\")"
    ]
}
```

---

## 🔄 数据加载工具函数

创建 `src/scLucid/utils/data_loader.py`:

```python
"""
数据加载工具函数

提供便捷的数据加载接口，用于测试、示例和文档。
"""

from pathlib import Path
from typing import Dict, Optional
import scanpy as sc
from anndata import AnnData

# 数据目录
DATA_DIR = Path(__file__).parent.parent.parent / "data"


def load_pbmc3k() -> AnnData:
    """
    加载PBMC3K数据集（正常组织）

    Returns
    -------
    adata : AnnData
        PBMC数据，已计算QC指标
    """
    from ..qc import calculate_qc_metrics

    adata = sc.read_h5ad(DATA_DIR / "pbmc3k" / "pbmc3k_raw.h5ad")
    adata = calculate_qc_metrics(adata)
    adata.obs['tissue_type'] = 'normal'

    return adata


def load_luad() -> AnnData:
    """
    加载LUAD数据集（肺腺癌，肿瘤组织）

    Returns
    -------
    adata : AnnData
        LUAD数据，已计算QC指标
    """
    from ..qc import calculate_qc_metrics

    adata = sc.read_h5ad(DATA_DIR / "human_LUAD_GSE131907" / "luad.h5ad")
    adata = calculate_qc_metrics(adata)
    adata.obs['tissue_type'] = 'lung_tumor'

    return adata


def load_melanoma() -> AnnData:
    """
    加载黑色素瘤数据集（多批次，肿瘤组织）

    Returns
    -------
    adata : AnnData
        黑色素瘤数据，已计算QC指标
    """
    from ..qc import calculate_qc_metrics

    adata = sc.read_h5ad(DATA_DIR / "mouse_melanoma_GSE119352" / "melanoma.h5ad")
    adata = calculate_qc_metrics(adata)
    adata.obs['tissue_type'] = 'melanoma'

    return adata


def load_all_datasets() -> Dict[str, AnnData]:
    """
    加载所有数据集

    Returns
    -------
    datasets : dict
        {dataset_name: AnnData}
    """
    return {
        'pbmc': load_pbmc3k(),
        'luad': load_luad(),
        'melanoma': load_melanoma()
    }


# 在tests/fixtures/data_loader.py中也可以使用
```

---

## 📊 数据集统计表

创建 `data/DATA_SUMMARY.md`:

```markdown
# Data文件夹数据集汇总

## 数据集统计

| 数据集 | 细胞数 | 基因数 | 组织类型 | 批次数 | 来源 | 文件大小 |
|--------|-------|-------|---------|-------|------|---------|
| PBMC3K | ~2,700 | ~2,000 | 正常（PBMC） | 1 | 10x Genomics | ~20MB |
| LUAD | TBD | TBD | 肺腺癌 | TBD | GEO GSE131907 | TBD |
| 黑色素瘤 | TBD | TBD | 黑色素瘤 | 多 | GEO GSE119352 | TBD |

## 使用场景

| 数据集 | 主要用途 | 适用测试 |
|--------|---------|---------|
| PBMC3K | 基线测试、快速示例 | • 统一阈值策略<br>• 基础功能测试<br>• 单样本分析 |
| LUAD | 肿瘤感知策略 | • 肿瘤vs正常对比<br>• MT阈值适应性<br>• 肿瘤微环境 |
| 黑色素瘤 | 混合策略、批次校正 | • 多批次处理<br>• 极端异质性<br>• 批次效应评估 |

## 引用

如果使用这些数据集，请引用：

- **PBMC3K**: 10x Genomics (https://support.10xgenomics.com/single-cell-gene-expression/datasets)
- **LUAD**: GSE131907 (原论文引用)
- **黑色素瘤**: GSE119352 (原论文引用)
```

---

## 🧪 测试覆盖率目标

### 数据集使用目标

| 数据集 | 单元测试 | 集成测试 | 示例 | 文档 |
|--------|---------|---------|------|------|
| PBMC3K | ✅ 高 | ✅ 高 | ✅ 是 | ✅ 是 |
| LUAD | ✅ 中 | ✅ 高 | ✅ 是 | ✅ 是 |
| 黑色素瘤 | ✅ 中 | ✅ 高 | ✅ 是 | ✅ 是 |

### 测试场景覆盖

| 场景 | PBMC3K | LUAD | 黑色素瘤 |
|------|--------|------|----------|
| 基础QC | ✅ | ✅ | ✅ |
| 智能QC | ✅ | ✅ | ✅ |
| 肿瘤感知 | ❌ | ✅ | ✅ |
| 多批次 | ❌ | ❌ | ✅ |
| 策略对比 | ❌ | ✅ | ✅ |

---

## 📝 最佳实践

### 1. 数据加载

```python
# ✅ 推荐：使用工具函数
from scLucid.utils.data_loader import load_pbmc3k, load_luad, load_melanoma
pbmc = load_pbmc3k()

# ❌ 不推荐：硬编码路径
adata = sc.read_h5ad('../../data/pbmc3k/pbmc3k_raw.h5ad')
```

### 2. 测试fixture

```python
# ✅ 推荐：使用fixture
@pytest.fixture
def pbmc_data():
    from scLucid.utils.data_loader import load_pbmc3k
    return load_pbmc3k()

def test_qc(pbmc_data):
    # 使用 pbmc_data
    pass

# ❌ 不推荐：在每个测试中加载数据
def test_qc():
    pbmc = load_pbmc3k()  # 重复加载
    pass
```

### 3. 示例代码

```python
# ✅ 推荐：提供完整context
# Load data
adata = load_pbmc3k()

# Get recommendation
rec = recommend_intelligent_qc(adata, tissue_type='normal')

# Print results
print(f"Threshold: {rec.min_genes.threshold}")

# ❌ 不推荐：跳过加载步骤
rec = recommend_intelligent_qc(adata, tissue_type='normal')  # adata从哪来？
```

---

**总结：data文件夹的三个数据集完美覆盖了QC策略验证的所有场景！** 🎯
