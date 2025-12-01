"""
BayesPrism Python Implementation
贝叶斯细胞比例重建推断统计边缘化方法的Python实现
"""

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import nbinom, dirichlet
from scipy.special import digamma, gammaln
from sklearn.decomposition import NMF
import warnings
from typing import Optional, Union, Tuple, List, Dict
from dataclasses import dataclass
import matplotlib.pyplot as plt
import seaborn as sns
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

# ==================== 核心数据类 ====================

@dataclass
class PrismConfig:
    """BayesPrism配置参数"""
    n_iter: int = 100  # Gibbs采样迭代次数
    n_chains: int = 4  # 马尔科夫链数量
    burnin: int = 50  # 预热迭代数
    update_bulk: bool = True  # 是否更新bulk表达
    pseudo_min: float = 1e-8  # 最小伪计数
    key: Optional[str] = None  # 肿瘤样本关键字
    outlier_cutoff: float = 0.01  # 异常值截断
    max_outlier: int = 10  # 最大异常基因数
    gibbs_control: Dict = None  # Gibbs采样控制参数
    
    def __post_init__(self):
        if self.gibbs_control is None:
            self.gibbs_control = {
                'chain_length': self.n_iter,
                'burn_in': self.burnin,
                'thinning': 1
            }


class BayesPrismReference:
    """
    BayesPrism参考数据类
    处理scRNA-seq参考数据
    """
    
    def __init__(
        self,
        reference: Union[pd.DataFrame, sparse.spmatrix],
        cell_type_labels: pd.Series,
        cell_state_labels: Optional[pd.Series] = None,
        input_type: str = "count.matrix",
        pseudo_min: float = 1e-8
    ):
        """
        初始化参考数据
        
        Parameters:
        -----------
        reference : 参考表达矩阵 (genes x cells)
        cell_type_labels : 细胞类型标签
        cell_state_labels : 细胞状态标签（可选）
        input_type : "count.matrix" 或 "GEP"（基因表达谱）
        pseudo_min : 伪计数最小值
        """
        self.reference = reference
        self.cell_type_labels = cell_type_labels
        self.cell_state_labels = cell_state_labels if cell_state_labels is not None else cell_type_labels
        self.input_type = input_type
        self.pseudo_min = pseudo_min
        
        # 转换为稀疏矩阵（如果需要）
        if isinstance(reference, pd.DataFrame):
            self.reference_matrix = sparse.csr_matrix(reference.values)
            self.gene_names = reference.index.tolist()
            self.cell_names = reference.columns.tolist()
        else:
            self.reference_matrix = reference
            self.gene_names = None
            self.cell_names = None
        
        # 生成参考谱
        self.cell_types = cell_type_labels.unique().tolist()
        self.cell_states = cell_state_labels.unique().tolist()
        self._generate_reference_profile()
    
    def _generate_reference_profile(self):
        """生成细胞类型特异性表达谱"""
        n_genes = self.reference_matrix.shape[0]
        n_cell_types = len(self.cell_types)
        
        # 初始化参考谱矩阵
        self.phi = np.zeros((n_genes, n_cell_types))
        
        for i, cell_type in enumerate(self.cell_types):
            # 获取该细胞类型的细胞索引
            cell_idx = self.cell_type_labels == cell_type
            
            if sparse.issparse(self.reference_matrix):
                type_expr = self.reference_matrix[:, cell_idx].toarray()
            else:
                type_expr = self.reference_matrix[:, cell_idx]
            
            # 计算平均表达（添加伪计数）
            total_counts = type_expr.sum(axis=0, keepdims=True)
            total_counts[total_counts == 0] = 1  # 避免除零
            
            # 归一化
            normalized = type_expr / total_counts
            mean_expr = normalized.mean(axis=1)
            
            # 添加伪计数避免零值
            mean_expr[mean_expr == 0] = self.pseudo_min
            
            self.phi[:, i] = mean_expr
        
        # 归一化参考谱（每列和为1）
        self.phi = self.phi / self.phi.sum(axis=0, keepdims=True)
        
    def get_cell_state_profile(self) -> Dict:
        """获取细胞状态特异性表达谱"""
        n_genes = self.reference_matrix.shape[0]
        state_profiles = {}
        
        for state in self.cell_states:
            cell_idx = self.cell_state_labels == state
            
            if sparse.issparse(self.reference_matrix):
                state_expr = self.reference_matrix[:, cell_idx].toarray()
            else:
                state_expr = self.reference_matrix[:, cell_idx]
            
            total_counts = state_expr.sum(axis=0, keepdims=True)
            total_counts[total_counts == 0] = 1
            
            normalized = state_expr / total_counts
            mean_expr = normalized.mean(axis=1)
            mean_expr[mean_expr == 0] = self.pseudo_min
            
            state_profiles[state] = mean_expr
        
        return state_profiles


class BayesPrism:
    """
    BayesPrism主类
    实现完整的贝叶斯反卷积算法
    """
    
    def __init__(
        self,
        reference: BayesPrismReference,
        mixture: pd.DataFrame,
        config: Optional[PrismConfig] = None
    ):
        """
        初始化BayesPrism对象
        
        Parameters:
        -----------
        reference : BayesPrismReference对象
        mixture : bulk RNA-seq混合表达矩阵 (genes x samples)
        config : 配置参数
        """
        self.reference = reference
        self.mixture = mixture
        self.config = config if config is not None else PrismConfig()
        
        # 对齐基因
        self._align_genes()
        
        # 初始化结果容器
        self.theta_initial = None  # θ0: 初始细胞类型比例
        self.theta_updated = None  # θf: 更新后的细胞类型比例
        self.Z = None  # Z: 细胞类型特异性表达
        self.posterior_samples = None  # 后验样本
        
    def _align_genes(self):
        """对齐参考和混合数据的基因"""
        ref_genes = set(self.reference.gene_names) if self.reference.gene_names else set()
        mix_genes = set(self.mixture.index)
        
        common_genes = sorted(list(ref_genes.intersection(mix_genes)))
        
        if len(common_genes) == 0:
            raise ValueError("参考和混合数据没有共同基因！")
        
        print(f"对齐后共有 {len(common_genes)} 个基因")
        
        # 重新排序
        if self.reference.gene_names:
            gene_idx = [self.reference.gene_names.index(g) for g in common_genes]
            self.reference.reference_matrix = self.reference.reference_matrix[gene_idx, :]
            self.reference.phi = self.reference.phi[gene_idx, :]
            self.reference.gene_names = common_genes
        
        self.mixture = self.mixture.loc[common_genes, :]
        self.aligned_genes = common_genes
        
    def cleanup_genes(
        self,
        remove_ribo: bool = True,
        remove_mito: bool = True,
        remove_sex: bool = True,
        min_expression: int = 0,
        min_cells: int = 0
    ):
        """
        清理基因
        
        Parameters:
        -----------
        remove_ribo : 是否移除核糖体基因
        remove_mito : 是否移除线粒体基因
        remove_sex : 是否移除性染色体基因
        min_expression : 最小表达量
        min_cells : 最小表达细胞数
        """
        genes_to_keep = []
        
        for gene in self.aligned_genes:
            gene_upper = gene.upper()
            
            # 检查核糖体基因
            if remove_ribo and (gene_upper.startswith('RPS') or gene_upper.startswith('RPL')):
                continue
            
            # 检查线粒体基因
            if remove_mito and gene_upper.startswith('MT-'):
                continue
            
            # 检查性染色体基因
            if remove_sex and (gene_upper.startswith('XIST') or gene_upper.startswith('Y')):
                continue
            
            genes_to_keep.append(gene)
        
        # 更新基因列表
        self._update_gene_subset(genes_to_keep)
        print(f"清理后剩余 {len(genes_to_keep)} 个基因")
    
    def _update_gene_subset(self, genes: List[str]):
        """更新基因子集"""
        if self.reference.gene_names:
            gene_idx = [self.reference.gene_names.index(g) for g in genes]
            self.reference.reference_matrix = self.reference.reference_matrix[gene_idx, :]
            self.reference.phi = self.reference.phi[gene_idx, :]
            self.reference.gene_names = genes
        
        self.mixture = self.mixture.loc[genes, :]
        self.aligned_genes = genes
    
    def select_markers(
        self,
        n_markers: int = 500,
        method: str = "t-test"
    ) -> List[str]:
        """
        选择标记基因
        
        Parameters:
        -----------
        n_markers : 每个细胞类型的标记基因数
        method : 差异表达方法
        
        Returns:
        --------
        marker_genes : 标记基因列表
        """
        marker_genes = set()
        n_cell_types = len(self.reference.cell_types)
        
        for i, cell_type in enumerate(self.reference.cell_types):
            # 获取该细胞类型的表达
            type_expr = self.reference.phi[:, i]
            
            # 获取其他细胞类型的平均表达
            other_expr = np.delete(self.reference.phi, i, axis=1).mean(axis=1)
            
            # 计算fold change
            fold_change = np.log2(type_expr + 1e-10) - np.log2(other_expr + 1e-10)
            
            # 选择top基因
            top_idx = np.argsort(fold_change)[-n_markers:]
            top_genes = [self.aligned_genes[idx] for idx in top_idx]
            
            marker_genes.update(top_genes)
        
        marker_genes = sorted(list(marker_genes))
        print(f"选择了 {len(marker_genes)} 个标记基因")
        
        return marker_genes
    
    def run_deconvolution(
        self,
        n_cores: int = 1,
        verbose: bool = True
    ):
        """
        运行反卷积
        
        Parameters:
        -----------
        n_cores : 并行核心数
        verbose : 是否显示进度
        """
        n_samples = self.mixture.shape[1]
        n_cell_types = len(self.reference.cell_types)
        n_genes = len(self.aligned_genes)
        
        # 初始化结果
        self.theta_initial = np.zeros((n_cell_types, n_samples))
        self.theta_updated = np.zeros((n_cell_types, n_samples))
        self.Z = np.zeros((n_genes, n_cell_types, n_samples))
        
        # 准备混合数据
        mixture_array = self.mixture.values
        
        if verbose:
            print("开始反卷积...")
            iterator = tqdm(range(n_samples), desc="处理样本")
        else:
            iterator = range(n_samples)
        
        # 对每个样本进行反卷积
        if n_cores > 1:
            with ProcessPoolExecutor(max_workers=n_cores) as executor:
                results = list(executor.map(
                    self._deconvolve_sample,
                    [(mixture_array[:, i], i) for i in range(n_samples)]
                ))
        else:
            results = [self._deconvolve_sample((mixture_array[:, i], i)) for i in iterator]
        
        # 整理结果
        for i, (theta0, theta_f, Z_sample) in enumerate(results):
            self.theta_initial[:, i] = theta0
            self.theta_updated[:, i] = theta_f
            self.Z[:, :, i] = Z_sample
        
        if verbose:
            print("反卷积完成！")
    
    def _deconvolve_sample(self, args: Tuple) -> Tuple:
        """
        对单个样本进行反卷积
        
        Parameters:
        -----------
        args : (mixture_vector, sample_index)
        
        Returns:
        --------
        theta0 : 初始细胞类型比例
        theta_f : 更新后的细胞类型比例
        Z : 细胞类型特异性表达
        """
        mixture_vector, sample_idx = args
        
        # 初始估计 (使用NNLS或其他方法)
        theta0 = self._initial_estimate(mixture_vector)
        
        # Gibbs采样
        theta_samples, Z_samples = self._gibbs_sampling(mixture_vector, theta0)
        
        # 计算后验均值
        theta_f = theta_samples.mean(axis=0)
        Z_mean = Z_samples.mean(axis=0)
        
        return theta0, theta_f, Z_mean
    
    def _initial_estimate(self, mixture: np.ndarray) -> np.ndarray:
        """
        初始估计细胞类型比例
        使用非负最小二乘法(NNLS)
        """
        from scipy.optimize import nnls
        
        # 归一化
        mixture_norm = mixture / (mixture.sum() + 1e-10)
        phi_norm = self.reference.phi
        
        # NNLS求解
        theta, _ = nnls(phi_norm, mixture_norm)
        
        # 归一化为概率
        theta = theta / (theta.sum() + 1e-10)
        
        return theta
    
    def _gibbs_sampling(
        self,
        mixture: np.ndarray,
        theta_init: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Gibbs采样估计后验分布
        
        Parameters:
        -----------
        mixture : 混合表达向量
        theta_init : 初始细胞类型比例
        
        Returns:
        --------
        theta_samples : θ的后验样本
        Z_samples : Z的后验样本
        """
        n_genes = len(mixture)
        n_cell_types = len(self.reference.cell_types)
        n_iter = self.config.n_iter
        burnin = self.config.burnin
        
        # 初始化
        theta = theta_init.copy()
        Z = np.zeros((n_genes, n_cell_types))
        
        # 存储样本
        theta_samples = []
        Z_samples = []
        
        # 总reads数
        N = mixture.sum()
        
        for iteration in range(n_iter + burnin):
            # Step 1: 采样 Z | θ, X
            Z = self._sample_Z(mixture, theta, self.reference.phi)
            
            # Step 2: 采样 θ | Z
            theta = self._sample_theta(Z)
            
            # 保存样本（跳过burnin）
            if iteration >= burnin:
                theta_samples.append(theta.copy())
                Z_samples.append(Z.copy())
        
        return np.array(theta_samples), np.array(Z_samples)
    
    def _sample_Z(
        self,
        mixture: np.ndarray,
        theta: np.ndarray,
        phi: np.ndarray
    ) -> np.ndarray:
        """
        采样细胞类型特异性表达 Z
        使用多项分布
        """
        n_genes, n_cell_types = phi.shape
        Z = np.zeros((n_genes, n_cell_types))
        
        for g in range(n_genes):
            if mixture[g] == 0:
                continue
            
            # 计算每个细胞类型的概率
            probs = phi[g, :] * theta
            probs = probs / (probs.sum() + 1e-10)
            
            # 多项采样
            Z[g, :] = np.random.multinomial(int(mixture[g]), probs)
        
        return Z
    
    def _sample_theta(self, Z: np.ndarray) -> np.ndarray:
        """
        采样细胞类型比例 θ
        使用Dirichlet分布
        """
        # 计算每个细胞类型的总reads
        cell_type_counts = Z.sum(axis=0)
        
        # Dirichlet先验参数（均匀先验）
        alpha = np.ones(len(cell_type_counts))
        
        # 后验参数
        alpha_posterior = alpha + cell_type_counts
        
        # 从Dirichlet分布采样
        theta = np.random.dirichlet(alpha_posterior)
        
        return theta
    
    def get_fraction(self, updated: bool = True) -> pd.DataFrame:
        """
        获取细胞类型比例
        
        Parameters:
        -----------
        updated : 是否使用更新后的θ
        
        Returns:
        --------
        fraction_df : 细胞类型比例DataFrame
        """
        theta = self.theta_updated if updated else self.theta_initial
        
        fraction_df = pd.DataFrame(
            theta.T,
            index=self.mixture.columns,
            columns=self.reference.cell_types
        )
        
        return fraction_df
    
    def get_expression(self, cell_type: str = None) -> pd.DataFrame:
        """
        获取细胞类型特异性表达
        
        Parameters:
        -----------
        cell_type : 细胞类型名称（None表示所有）
        
        Returns:
        --------
        expression_df : 表达DataFrame
        """
        if cell_type is None:
            # 返回所有细胞类型
            result = {}
            for i, ct in enumerate(self.reference.cell_types):
                result[ct] = pd.DataFrame(
                    self.Z[:, i, :],
                    index=self.aligned_genes,
                    columns=self.mixture.columns
                )
            return result
        else:
            ct_idx = self.reference.cell_types.index(cell_type)
            return pd.DataFrame(
                self.Z[:, ct_idx, :],
                index=self.aligned_genes,
                columns=self.mixture.columns
            )
    
    def compute_cv(self) -> pd.DataFrame:
        """
        计算细胞类型比例的变异系数(CV)
        量化后验分布的不确定性
        """
        # 这里需要后验样本，简化实现
        # 实际应该从Gibbs采样保存的样本计算
        cv = np.std(self.theta_updated, axis=1) / (np.mean(self.theta_updated, axis=1) + 1e-10)
        
        cv_df = pd.DataFrame({
            'cell_type': self.reference.cell_types,
            'CV': cv
        })
        
        return cv_df
    
    def plot_fraction(self, figsize: Tuple = (12, 6)):
        """绘制细胞类型比例热图"""
        fraction_df = self.get_fraction(updated=True)
        
        plt.figure(figsize=figsize)
        sns.heatmap(
            fraction_df.T,
            cmap='YlOrRd',
            cbar_kws={'label': 'Fraction'},
            xticklabels=True,
            yticklabels=True
        )
        plt.title('Cell Type Fractions')
        plt.xlabel('Samples')
        plt.ylabel('Cell Types')
        plt.tight_layout()
        plt.show()
    
    def plot_correlation(self, figsize: Tuple = (10, 8)):
        """绘制细胞类型之间的相关性"""
        fraction_df = self.get_fraction(updated=True)
        corr = fraction_df.corr()
        
        plt.figure(figsize=figsize)
        sns.heatmap(
            corr,
            annot=True,
            fmt='.2f',
            cmap='coolwarm',
            center=0,
            square=True
        )
        plt.title('Cell Type Correlation')
        plt.tight_layout()
        plt.show()


# ==================== 嵌入学习模块 ====================

class BayesPrismEmbedding:
    """
    BayesPrism嵌入学习模块
    使用EM算法近似肿瘤表达
    """
    
    def __init__(
        self,
        prism: BayesPrism,
        tumor_key: str,
        n_programs: int = 5
    ):
        """
        初始化嵌入学习
        
        Parameters:
        -----------
        prism : BayesPrism对象（已完成反卷积）
        tumor_key : 肿瘤细胞类型标识
        n_programs : 恶性基因程序数量
        """
        self.prism = prism
        self.tumor_key = tumor_key
        self.n_programs = n_programs
        
        # 获取肿瘤表达
        self.tumor_expr = self._extract_tumor_expression()
        
        # 初始化结果
        self.W = None  # 基因程序矩阵
        self.H = None  # 程序系数矩阵
    
    def _extract_tumor_expression(self) -> pd.DataFrame:
        """提取肿瘤特异性表达"""
        tumor_Z = self.prism.get_expression(self.tumor_key)
        return tumor_Z
    
    def run_nmf(
        self,
        max_iter: int = 200,
        tol: float = 1e-4,
        verbose: bool = True
    ):
        """
        运行非负矩阵分解(NMF)学习基因程序
        
        Parameters:
        -----------
        max_iter : 最大迭代次数
        tol : 收敛容差
        verbose : 是否显示进度
        """
        if verbose:
            print(f"学习 {self.n_programs} 个恶性基因程序...")
        
        # 使用NMF分解
        model = NMF(
            n_components=self.n_programs,
            init='random',
            max_iter=max_iter,
            tol=tol,
            random_state=42
        )
        
        # 拟合
        self.W = model.fit_transform(self.tumor_expr.values)
        self.H = model.components_
        
        # 转换为DataFrame
        self.W_df = pd.DataFrame(
            self.W,
            index=self.tumor_expr.index,
            columns=[f'Program_{i+1}' for i in range(self.n_programs)]
        )
        
        self.H_df = pd.DataFrame(
            self.H,
            index=[f'Program_{i+1}' for i in range(self.n_programs)],
            columns=self.tumor_expr.columns
        )
        
        if verbose:
            print("基因程序学习完成！")
    
    def get_gene_programs(self) -> pd.DataFrame:
        """获取基因程序矩阵"""
        return self.W_df
    
    def get_program_usage(self) -> pd.DataFrame:
        """获取程序使用系数"""
        return self.H_df
    
    def plot_gene_programs(self, top_n: int = 20, figsize: Tuple = (15, 10)):
        """
        可视化基因程序
        
        Parameters:
        -----------
        top_n : 每个程序显示的top基因数
        figsize : 图形大小
        """
        fig, axes = plt.subplots(
            (self.n_programs + 1) // 2, 2,
            figsize=figsize
        )
        axes = axes.flatten()
        
        for i in range(self.n_programs):
            program = self.W_df[f'Program_{i+1}']
            top_genes = program.nlargest(top_n)
            
            axes[i].barh(range(top_n), top_genes.values[::-1])
            axes[i].set_yticks(range(top_n))
            axes[i].set_yticklabels(top_genes.index[::-1])
            axes[i].set_xlabel('Weight')
            axes[i].set_title(f'Program {i+1}')
        
        plt.tight_layout()
        plt.show()
    
    def plot_program_usage(self, figsize: Tuple = (12, 6)):
        """可视化程序使用情况"""
        plt.figure(figsize=figsize)
        sns.heatmap(
            self.H_df,
            cmap='viridis',
            cbar_kws={'label': 'Usage'},
            xticklabels=True,
            yticklabels=True
        )
        plt.title('Malignant Gene Program Usage')
        plt.xlabel('Samples')
        plt.ylabel('Programs')
        plt.tight_layout()
        plt.show()


# ==================== 工具函数 ====================

def cleanup_genes(
    gene_names: List[str],
    remove_ribo: bool = True,
    remove_mito: bool = True,
    remove_sex: bool = True
) -> List[str]:
    """
    清理基因列表
    
    Parameters:
    -----------
    gene_names : 基因名称列表
    remove_ribo : 是否移除核糖体基因
    remove_mito : 是否移除线粒体基因
    remove_sex : 是否移除性染色体基因
    
    Returns:
    --------
    cleaned_genes : 清理后的基因列表
    """
    cleaned_genes = []
    
    for gene in gene_names:
        gene_upper = gene.upper()
        
        if remove_ribo and (gene_upper.startswith('RPS') or gene_upper.startswith('RPL')):
            continue
        
        if remove_mito and gene_upper.startswith('MT-'):
            continue
        
        if remove_sex and (gene_upper.startswith('XIST') or 'CHROMOSOME_Y' in gene_upper):
            continue
        
        cleaned_genes.append(gene)
    
    return cleaned_genes


def find_outlier_genes(
    mixture: pd.DataFrame,
    reference_genes: List[str],
    cutoff: float = 0.01
) -> List[str]:
    """
    查找异常基因
    
    Parameters:
    -----------
    mixture : 混合表达矩阵
    reference_genes : 参考基因列表
    cutoff : 异常值截断阈值
    
    Returns:
    --------
    outlier_genes : 异常基因列表
    """
    outliers = []
    
    for gene in mixture.index:
        if gene not in reference_genes:
            # 检查是否在混合中高表达
            expr_ratio = (mixture.loc[gene] > mixture.max() * cutoff).sum() / len(mixture.columns)
            if expr_ratio > 0.1:
                outliers.append(gene)
    
    return outliers


def compute_correlation(
    deconv_result: pd.DataFrame,
    true_fraction: pd.DataFrame
) -> pd.DataFrame:
    """
    计算反卷积结果与真实比例的相关性
    
    Parameters:
    -----------
    deconv_result : 反卷积结果
    true_fraction : 真实细胞类型比例
    
    Returns:
    --------
    correlation_df : 相关性DataFrame
    """
    from scipy.stats import pearsonr, spearmanr
    
    results = []
    
    for cell_type in deconv_result.columns:
        if cell_type in true_fraction.columns:
            pred = deconv_result[cell_type].values
            true = true_fraction[cell_type].values
            
            pearson_r, pearson_p = pearsonr(pred, true)
            spearman_r, spearman_p = spearmanr(pred, true)
            
            results.append({
                'cell_type': cell_type,
                'pearson_r': pearson_r,
                'pearson_p': pearson_p,
                'spearman_r': spearman_r,
                'spearman_p': spearman_p
            })
    
    return pd.DataFrame(results)


# ==================== 示例使用代码 ====================

def example_usage():
    """
    示例：如何使用BayesPrism Python实现
    """
    
    # 1. 准备数据
    # 加载scRNA-seq参考数据（示例）
    # reference_counts: genes x cells的表达矩阵
    # cell_type_labels: 细胞类型标签
    
    print("=" * 60)
    print("BayesPrism Python实现示例")
    print("=" * 60)
    
    # 模拟数据（实际使用时替换为真实数据）
    n_genes = 2000
    n_cells = 500
    n_samples = 50
    n_cell_types = 5
    
    # 模拟参考数据
    np.random.seed(42)
    reference_counts = sparse.random(n_genes, n_cells, density=0.1, format='csr')
    reference_counts.data = np.random.poisson(5, size=reference_counts.data.shape)
    
    cell_type_labels = pd.Series(
        np.random.choice(['T_cell', 'B_cell', 'Macrophage', 'Fibroblast', 'Tumor'],
                        size=n_cells)
    )
    
    gene_names = [f'Gene_{i}' for i in range(n_genes)]
    
    # 转换为DataFrame（可选）
    reference_df = pd.DataFrame(
        reference_counts.toarray(),
        index=gene_names,
        columns=[f'Cell_{i}' for i in range(n_cells)]
    )
    
    # 2. 创建参考对象
    print("\n步骤1: 创建参考数据对象")
    ref = BayesPrismReference(
        reference=reference_df,
        cell_type_labels=cell_type_labels,
        input_type="count.matrix"
    )
    print(f"参考包含 {len(ref.cell_types)} 个细胞类型")
    
    # 3. 准备混合数据
    print("\n步骤2: 准备bulk混合数据")
    # 模拟混合数据
    mixture_data = np.random.poisson(
        10,
        size=(n_genes, n_samples)
    )
    mixture_df = pd.DataFrame(
        mixture_data,
        index=gene_names,
        columns=[f'Sample_{i}' for i in range(n_samples)]
    )
    
    # 4. 创建BayesPrism对象
    print("\n步骤3: 创建BayesPrism对象")
    config = PrismConfig(
        n_iter=50,  # 示例用较少迭代
        burnin=25,
        n_chains=2
    )
    
    bp = BayesPrism(
        reference=ref,
        mixture=mixture_df,
        config=config
    )
    
    # 5. 清理基因
    print("\n步骤4: 清理基因")
    bp.cleanup_genes(
        remove_ribo=True,
        remove_mito=True,
        remove_sex=True
    )
    
    # 6. 运行反卷积
    print("\n步骤5: 运行反卷积")
    bp.run_deconvolution(n_cores=1, verbose=True)
    
    # 7. 获取结果
    print("\n步骤6: 获取结果")
    fraction_df = bp.get_fraction(updated=True)
    print("\n细胞类型比例（前5个样本）:")
    print(fraction_df.head())
    
    # 8. 计算CV
    cv_df = bp.compute_cv()
    print("\n细胞类型比例的变异系数:")
    print(cv_df)
    
    # 9. 可视化
    print("\n步骤7: 可视化结果")
    # bp.plot_fraction()  # 取消注释以显示
    # bp.plot_correlation()
    
    # 10. 嵌入学习（如果有肿瘤细胞）
    print("\n步骤8: 嵌入学习")
    if 'Tumor' in ref.cell_types:
        embedding = BayesPrismEmbedding(
            prism=bp,
            tumor_key='Tumor',
            n_programs=3
        )
        embedding.run_nmf(verbose=True)
        
        programs = embedding.get_gene_programs()
        print("\n基因程序（前5个基因）:")
        print(programs.head())
        
        # embedding.plot_gene_programs()  # 取消注释以显示
        # embedding.plot_program_usage()
    
    print("\n" + "=" * 60)
    print("示例完成!")
    print("=" * 60)


if __name__ == "__main__":
    # 运行示例
    example_usage()