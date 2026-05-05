"""
抽象基类插件开发示例

这个示例展示如何使用抽象基类创建自定义分析插件。
"""

from scLucid.base_interfaces import AnalysisStep, AnalysisStepFactory
from scLucid.base_config import SclucidBaseConfig
from anndata import AnnData
from typing import Optional
from pydantic import Field

# ========================================
# 示例1: 创建自定义QC步骤
# ========================================

class CustomQCConfig(SclucidBaseConfig):
    """自定义QC配置"""
    min_genes: int = Field(default=200, ge=0, description="最小基因数")
    max_mt_percent: float = Field(default=20.0, ge=0, le=100, description="最大线粒体百分比")

class HighStringencyQC(AnalysisStep):
    """
    高严紧度QC过滤器

    这是对标准QC的自定义实现，提供了更严格的过滤标准。
    """

    def __init__(self, config: Optional[CustomQCConfig] = None):
        self.config = config or CustomQCConfig()
        self._results = None

    def validate_input(self, adata: AnnData) -> bool:
        """验证输入数据"""
        if not isinstance(adata, AnnData):
            raise ValueError("输入必须是 AnnData 对象")

        if adata.n_obs == 0:
            raise ValueError("数据中没有细胞")

        return True

    def run(self, adata: AnnData, **kwargs) -> AnnData:
        """执行高严紧度QC过滤"""
        import scanpy as sc

        # 应用过滤
        sc.pp.filter_cells(adata, min_genes=self.config.min_genes)
        adata.obs['pct_counts_mt'] = adata.obs['pct_counts_mt'].astype(float)

        # 过滤高线粒体含量的细胞
        adata = adata[adata.obs['pct_counts_mt'] < self.config.max_mt_percent, :].copy()

        # 存储结果
        self._results = {
            'n_cells_before': adata.n_obs,
            'min_genes': self.config.min_genes,
            'max_mt_percent': self.config.max_mt_percent
        }

        return adata

    def get_summary(self) -> dict:
        """获取过滤摘要"""
        return self._results or {"status": "not_run"}

# ========================================
# 示例2: 创建自定义注释方法
# ========================================

class CustomAnnotatorConfig(SclucidBaseConfig):
    """自定义注释器配置"""
    reference_path: str = Field(description="参考数据路径")
    similarity_threshold: float = Field(default=0.8, ge=0, le=1)

class MyCustomAnnotator(AnalysisStep):
    """
    自定义细胞类型注释器

    使用自己的算法或模型进行细胞类型注释。
    """

    def __init__(self, config: Optional[CustomAnnotatorConfig] = None):
        self.config = config or CustomAnnotatorConfig()
        self._results = None

    def validate_input(self, adata: AnnData) -> bool:
        """验证输入"""
        if 'X_pca' not in adata.obsm:
            raise ValueError("需要先进行PCA降维")
        return True

    def run(self, adata: AnnData, **kwargs) -> AnnData:
        """运行自定义注释算法"""
        # 这里可以实现任何自定义的注释逻辑
        # 例如：使用自己的机器学习模型、数据库查找等

        # 示例：随机注释（实际使用时替换为真实算法）
        import numpy as np
        cell_types = ['T cells', 'B cells', 'NK cells', 'Monocytes']
        predicted = np.random.choice(cell_types, size=adata.n_obs)

        adata.obs['custom_cell_type'] = predicted
        adata.obs['custom_cell_type_confidence'] = np.random.random(adata.n_obs)

        self._results = {
            'method': 'custom_annotation',
            'n_cell_types': len(cell_types)
        }

        return adata

    def get_summary(self) -> dict:
        return self._results

# ========================================
# 示例3: 注册和使用插件
# ========================================

def register_custom_plugins():
    """注册自定义插件到工厂"""

    # 注册自定义QC步骤
    AnalysisStepFactory.register(
        name='high_stringency_qc',
        step_class=HighStringencyQC
    )

    # 注册自定义注释器
    AnalysisStepFactory.register(
        name='my_annotator',
        step_class=MyCustomAnnotator
    )

    print("✓ 自定义插件已注册")

# ========================================
# 示例4: 使用插件
# ========================================

def example_usage():
    """展示如何使用插件"""

    # 1. 直接使用插件类
    from scLucid import AnnData
    import scanpy as sc

    # 创建测试数据
    adata = sc.datasets.pbmc3k()

    # 使用自定义QC
    qc_filter = HighStringencyQC()
    qc_filter.validate_input(adata)
    adata_filtered = qc_filter.run(adata)
    print(f"过滤后细胞数: {adata_filtered.n_obs}")

    # 2. 通过工厂使用插件
    register_custom_plugins()

    qc_step = AnalysisStepFactory.create('high_stringency_qc')
    annotator = AnalysisStepFactory.create('my_annotator')

    # 3. 组合多个插件
    adata = qc_step.run(adata)
    adata = annotator.run(adata)

    return adata

# ========================================
# 示例5: 插件与工作流集成
# ========================================

def create_custom_workflow():
    """创建包含自定义步骤的工作流"""

    from scLucid.analysis import run_custom_analysis

    # 定义自定义工作流步骤
    steps = ['high_stringency_qc', 'clustering', 'my_annotator']

    # 配置每一步
    step_configs = {
        'high_stringency_qc': {
            'config': CustomQCConfig(min_genes=300, max_mt_percent=15.0)
        },
        'clustering': {'resolution': 0.8},
        'my_annotator': {'reference_path': '/path/to/reference.h5ad'}
    }

    # 运行自定义工作流
    # adata = run_custom_analysis(adata, steps=steps, step_configs=step_configs)

if __name__ == '__main__':
    print("插件开发示例")
    print("=" * 70)
    print("\n这个文件展示了如何使用抽象基类创建自定义分析插件。")
    print("\n关键概念：")
    print("1. 继承抽象基类（如 AnalysisStep）")
    print("2. 实现必需方法：validate_input(), run()")
    print("3. 使用工厂注册：AnalysisStepFactory.register()")
    print("4. 通过工厂创建实例：AnalysisStepFactory.create()")
    print("\n这样做的优势：")
    print("- 不需要修改核心代码")
    print("- 插件可独立开发和测试")
    print("- 插件可动态加载和卸载")
    print("- 保持核心代码的稳定性")
