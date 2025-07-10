import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
from pathlib import Path

# 导入自定义包
import scRNA
from scRNA.analysis import Manager, annotate_clusters, score_cell_types, marker_guided_clustering
from scRNA.analysis import find_markers, marker_enrichment_analysis
from scRNA.analysis import plot_marker_expression, plot_cell_type_composition, umap_with_annotated_clusters
from scRNA.analysis.cluster import merge_clusters

def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='单细胞RNA-seq分析流程')
    parser.add_argument('--input', type=str, help='输入数据路径 (.h5ad 格式)')
    parser.add_argument('--output', type=str, default='results', help='输出目录')
    parser.add_argument('--sample_key', type=str, default='sampleID', help='样本信息所在的列名')
    parser.add_argument('--species', type=str, default='human', choices=['human', 'mouse'], help='物种类型')
    parser.add_argument('--marker_config', type=str, help='Marker配置文件路径')
    parser.add_argument('--resolution', type=float, default=0.8, help='Leiden聚类分辨率')
    parser.add_argument('--skip_qc', action='store_true', help='跳过质量控制步骤')
    parser.add_argument('--skip_pp', action='store_true', help='跳过预处理步骤')
    parser.add_argument('--skip_batch', action='store_true', help='跳过批次校正步骤')
    parser.add_argument('--skip_analysis', action='store_true', help='跳过高级分析步骤')
    
    return parser.parse_args()

def main():
    """主函数：执行完整分析流程"""
    args = parse_args()
    
    # 创建输出目录
    results_dir = args.output
    os.makedirs(results_dir, exist_ok=True)
    
    # 设置随机种子确保结果可重现
    np.random.seed(42)
    
    # 加载数据
    print("=== 加载数据 ===")
    if args.input:
        adata = sc.read_h5ad(args.input)
        print(f"从{args.input}加载数据: {adata.shape[0]}个细胞, {adata.shape[1]}个基因")
    else:
        # 使用示例数据集
        adata = sc.datasets.pbmc3k()
        print(f"加载示例数据集: {adata.shape[0]}个细胞, {adata.shape[1]}个基因")
        
        # 为示例数据添加批次信息
        if args.sample_key not in adata.obs:
            adata.obs[args.sample_key] = np.random.choice(['sample1', 'sample2', 'sample3'], size=adata.n_obs)
            adata.obs[args.sample_key] = adata.obs[args.sample_key].astype('category')
            print(f"已添加模拟样本信息 ({args.sample_key})")

    # 选择适当的marker配置文件
    if args.marker_config:
        marker_config = args.marker_config
    else:
        # 默认配置文件
        if args.species == 'human':
            marker_config = os.path.join('marker_configs', 'manager_human.toml')
        else:
            marker_config = os.path.join('marker_configs', 'manager_mouse.toml')
        
        # 确保配置文件存在
        if not os.path.exists(marker_config):
            print(f"警告: 默认marker配置文件'{marker_config}'不存在，将创建临时配置文件")
            # 创建临时配置文件目录
            os.makedirs(os.path.dirname(marker_config), exist_ok=True)
            # 将样例配置写入文件
            with open(marker_config, 'w') as f:
                if args.species == 'human':
                    f.write(human_marker_sample_config)
                else:
                    f.write(mouse_marker_sample_config)
    
    print(f"使用marker配置: {marker_config}")

    # 执行质量控制
    if not args.skip_qc:
        print("\n=== 质量控制 ===")
        # 计算QC指标
        adata = scRNA.qc.calculate_qc_metric(
            adata, 
            sample_key=args.sample_key,
            plot_violin=True,
            plot_scatter=True,
            save_dir=os.path.join(results_dir, "qc_metrics")
        )

        # 标记低质量细胞
        adata = scRNA.qc.is_low_quality_cell(
            adata,
            sample_key=args.sample_key,
            min_genes=200,
            pc_mt=10,
            plot_outliers=True,
            save_dir=os.path.join(results_dir, "qc_metrics")
        )

        # 检测双细胞
        adata = scRNA.qc.is_doublet(
            adata,
            sample_key=args.sample_key,
            rate=0.05,
            plot_umap=True,
            save_dir=os.path.join(results_dir, "doublets")
        )

        # 过滤低质量细胞和双细胞
        adata_filtered = scRNA.qc.filter_low_quality_cells(
            adata,
            filter_low_genes=True,
            filter_outliers=True,
            filter_mt=True,
            filter_doublets=True
        )

        print(f"质量控制后的细胞数: {adata_filtered.n_obs}")
    else:
        print("\n跳过质量控制")
        adata_filtered = adata.copy()

    # 执行预处理
    if not args.skip_pp:
        print("\n=== 预处理 ===")
        # 数据标准化
        adata_filtered = scRNA.pp.normalize_data(
            adata_filtered,
            method="standard",
            target_sum=1e4,
            log_transform=True,
            plot=True,
            save_dir=os.path.join(results_dir, "normalization")
        )

        # 细胞周期评分
        adata_filtered = scRNA.pp.score_cell_cycle(
            adata_filtered,
            species=args.species,
            plot=True,
            save_dir=os.path.join(results_dir, "cell_cycle")
        )

        # 回归掉细胞周期效应
        adata_filtered = scRNA.pp.regress_out(
            adata_filtered,
            keys=["S_score", "G2M_score"],
            layer="log1p_norm",
            output_layer="cell_cycle_regressed"
        )

        # 标记高变异基因
        adata_filtered = scRNA.pp.annotate_hvg(
            adata_filtered,
            method="scanpy",
            n_top_genes_scanpy=2000,
            layer="log1p_norm",
            batch_key=args.sample_key
        )

        # 选择高变异基因
        adata_filtered = scRNA.pp.select_hvg(
            adata_filtered,
            method="scanpy",
            subset=True,
            n_top_genes=2000
        )

        print(f"选择的高变异基因数: {adata_filtered.n_vars}")

        # 数据缩放
        adata_filtered = scRNA.pp.scale_data(
            adata_filtered,
            max_value=10.0,
            layer="log1p_norm",
            output_layer="scaled"
        )

        # 执行PCA
        print("计算PCA...")
        sc.pp.pca(adata_filtered, n_comps=50)
    else:
        print("\n跳过预处理")
        # 确保有基本的预处理结果
        if "log1p_norm" not in adata_filtered.layers:
            print("执行基本的标准化...")
            sc.pp.normalize_total(adata_filtered, target_sum=1e4)
            sc.pp.log1p(adata_filtered)
            adata_filtered.layers["log1p_norm"] = adata_filtered.X.copy()
        
        if "X_pca" not in adata_filtered.obsm:
            print("执行基本的PCA...")
            sc.pp.pca(adata_filtered, n_comps=50)

    # 执行批次校正
    if not args.skip_batch and args.sample_key in adata_filtered.obs:
        print("\n=== 批次校正 ===")
        # 检查样本数量
        n_samples = len(adata_filtered.obs[args.sample_key].unique())
        if n_samples > 1:
            adata_filtered = scRNA.pp.batch_correction(
                adata_filtered,
                batch_key=args.sample_key,
                method="harmony",
                n_pcs=30,
                plot=True,
                save_dir=os.path.join(results_dir, "batch_correction")
            )
            # 使用批次校正后的表示进行邻居计算
            sc.pp.neighbors(adata_filtered, use_rep="X_harmony", n_neighbors=15)
        else:
            print(f"只有一个样本，跳过批次校正")
            sc.pp.neighbors(adata_filtered, n_neighbors=15)
    else:
        print("\n跳过批次校正")
        # 使用PCA表示计算邻居
        if "neighbors" not in adata_filtered.uns:
            sc.pp.neighbors(adata_filtered, n_neighbors=15)

    # 计算UMAP降维
    print("\n=== 降维与基础聚类 ===")
    sc.tl.umap(adata_filtered)
    
    # 使用marker引导的聚类找到最佳分辨率
    if not args.skip_analysis:
        # 使用预先选定的分辨率范围尝试marker-guided聚类
        print("\n执行marker引导的聚类...")
        try:
            adata_filtered = marker_guided_clustering(
                adata_filtered,
                marker_config=marker_config,
                resolution_range=(0.1, 1.5, 10),
                metric="marker_separation",
                clustering_method="leiden",
                key_added="leiden_optimized",
                plot=True
            )
            # 使用优化的聚类结果
            cluster_key = "leiden_optimized"
        except Exception as e:
            print(f"Marker引导聚类失败: {str(e)}")
            print("使用标准Leiden聚类...")
            # 退回到标准聚类
            sc.tl.leiden(adata_filtered, resolution=args.resolution)
            cluster_key = "leiden"
    else:
        # 使用标准Leiden聚类
        print("\n执行标准Leiden聚类...")
        sc.tl.leiden(adata_filtered, resolution=args.resolution)
        cluster_key = "leiden"

    # 可视化基础结果
    print("\n绘制基础UMAP结果...")
    sc.settings.set_figure_params(figsize=(12, 10))
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    sc.pl.umap(adata_filtered, color=cluster_key, ax=axes[0, 0], show=False, title="Clusters")
    if args.sample_key in adata_filtered.obs:
        sc.pl.umap(adata_filtered, color=args.sample_key, ax=axes[0, 1], show=False, title="Samples")
    else:
        axes[0, 1].text(0.5, 0.5, "No sample information", ha='center', va='center')
        axes[0, 1].set_title("Samples")
        axes[0, 1].axis('off')
    
    if "phase" in adata_filtered.obs:
        sc.pl.umap(adata_filtered, color="phase", ax=axes[1, 0], show=False, title="Cell Cycle Phase")
    else:
        axes[1, 0].text(0.5, 0.5, "No cell cycle information", ha='center', va='center')
        axes[1, 0].set_title("Cell Cycle Phase")
        axes[1, 0].axis('off')
    
    sc.pl.umap(adata_filtered, color="n_genes_by_counts", ax=axes[1, 1], show=False, title="Gene Count")

    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, "basic_umap.png"), dpi=300)
    plt.close()

    # 执行高级分析
    if not args.skip_analysis:
        print("\n=== 高级分析 ===")
        
        # 1. 基于已知的marker基因为聚类打分
        print("\n计算细胞类型得分...")
        adata_filtered = score_cell_types(
            adata_filtered,
            marker_config=marker_config,
            key_added="cell_type_scores",
            min_genes=2,
            method="scanpy"
        )

        # 查看前几个得分列
        score_cols = [col for col in adata_filtered.obs.columns if col.endswith("_score")]
        print(f"生成了 {len(score_cols)} 个细胞类型得分")

        # 2. 使用多种方法进行自动细胞类型注释
        print("\n执行聚类注释...")
        print("方法1: 基于相关性的注释")
        try:
            adata_filtered = annotate_clusters(
                adata_filtered,
                cluster_key=cluster_key,
                marker_config=marker_config,
                method="correlation",
                key_added="annotation_correlation"
            )
            
            print("方法2: 基于最大得分的注释")
            adata_filtered = annotate_clusters(
                adata_filtered,
                cluster_key=cluster_key,
                marker_config=marker_config,
                method="max_score",
                key_added="annotation_max_score"
            )
            
            print("方法3: 基于标记富集的注释")
            adata_filtered = annotate_clusters(
                adata_filtered,
                cluster_key=cluster_key,
                marker_config=marker_config,
                method="marker_enrichment",
                key_added="annotation_enrichment"
            )
            
            # 可视化不同注释方法的结果
            fig, axes = plt.subplots(1, 3, figsize=(18, 6))
            
            sc.pl.umap(adata_filtered, color="annotation_correlation", ax=axes[0], 
                      show=False, title="Correlation-based")
            sc.pl.umap(adata_filtered, color="annotation_max_score", ax=axes[1], 
                      show=False, title="Max-score-based")
            sc.pl.umap(adata_filtered, color="annotation_enrichment", ax=axes[2], 
                      show=False, title="Enrichment-based")
            
            plt.tight_layout()
            plt.savefig(os.path.join(results_dir, "annotation_methods_comparison.png"), dpi=300)
            plt.close()
            
            # 选择默认注释方法
            annotation_key = "annotation_correlation"
        except Exception as e:
            print(f"自动注释失败: {str(e)}")
            annotation_key = None

        # 3. 对每个聚类找到marker基因
        print("\n寻找聚类marker基因...")
        adata_filtered = find_markers(
            adata_filtered,
            groupby=cluster_key,
            method="wilcoxon",
            pts=True,
            min_fold_change=1.5,
            max_pval=0.05,
            filter_genes=True,
            plot=True,
            key_added=f"{cluster_key}_markers"
        )

        # 4. 使用Manager展示细胞类型信息
        print("\n显示细胞类型marker信息...")
        mgr = Manager(marker_config)
        mgr.intersect_with(adata_filtered)
        mgr.show_clusters()
        
        # 5. 可视化marker基因表达
        print("\n可视化marker基因表达...")
        try:
            # 选择前几个主要细胞类型
            major_types = [ct for ct, cell in mgr.CELLS.items() if cell.level == "major"][:3]
            if major_types:
                plot_marker_expression(
                    adata_filtered,
                    marker_config=marker_config,
                    cell_types=major_types,
                    basis="umap",
                    n_markers=2,
                    figsize=(12, 8),
                    ncols=3,
                    save=os.path.join(results_dir, "marker_expression.png")
                )
        except Exception as e:
            print(f"绘制marker表达失败: {str(e)}")

        # 6. 可视化聚类中的细胞类型组成
        if annotation_key:
            print("\n可视化细胞类型组成...")
            try:
                plot_cell_type_composition(
                    adata_filtered,
                    cluster_key=cluster_key,
                    annotation_key=annotation_key,
                    marker_config=marker_config,
                    normalize=True,
                    plot_type="bar",
                    save=os.path.join(results_dir, "cell_type_composition_bar.png")
                )
                
                plot_cell_type_composition(
                    adata_filtered,
                    cluster_key=cluster_key,
                    annotation_key=annotation_key,
                    marker_config=marker_config,
                    normalize=True,
                    plot_type="heatmap",
                    save=os.path.join(results_dir, "cell_type_composition_heatmap.png")
                )
            except Exception as e:
                print(f"绘制细胞类型组成失败: {str(e)}")
        
        # 7. 合并相似聚类
        print("\n合并相似聚类...")
        try:
            adata_filtered = merge_clusters(
                adata_filtered,
                cluster_key=cluster_key,
                marker_config=marker_config,
                similarity_threshold=0.7,
                method="marker_overlap",
                key_added=f"{cluster_key}_merged"
            )
            
            # 可视化合并结果
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
            sc.pl.umap(adata_filtered, color=cluster_key, ax=axes[0], 
                      show=False, title="Original clusters")
            sc.pl.umap(adata_filtered, color=f"{cluster_key}_merged", ax=axes[1], 
                      show=False, title="Merged clusters")
            plt.tight_layout()
            plt.savefig(os.path.join(results_dir, "cluster_merging.png"), dpi=300)
            plt.close()
        except Exception as e:
            print(f"合并聚类失败: {str(e)}")
        
        # 8. 富集分析
        if f"{cluster_key}_markers_dataframe" in adata_filtered.uns:
            print("\n执行marker富集分析...")
            try:
                # 获取前100个差异表达基因
                de_genes = adata_filtered.uns[f"{cluster_key}_markers_dataframe"]["names"].head(100).tolist()
                
                enrichment_results = marker_enrichment_analysis(
                    adata_filtered,
                    de_genes=de_genes,
                    marker_config=marker_config,
                    method="hypergeometric",
                    min_genes=2,
                    plot=True,
                    n_top=10
                )
                
                # 保存富集结果
                enrichment_results.to_csv(os.path.join(results_dir, "enrichment_results.csv"))
            except Exception as e:
                print(f"富集分析失败: {str(e)}")
        
        # 9. 在UMAP上标注聚类
        if annotation_key:
            print("\n创建带注释的UMAP可视化...")
            try:
                umap_with_annotated_clusters(
                    adata_filtered,
                    cluster_key=cluster_key,
                    annotation_key=annotation_key,
                    color_by="annotation",
                    text_size=10,
                    save=os.path.join(results_dir, "annotated_umap.png")
                )
            except Exception as e:
                print(f"创建带注释的UMAP失败: {str(e)}")

    # 保存处理后的数据
    final_h5ad = os.path.join(results_dir, "analyzed_data.h5ad")
    print(f"\n保存分析后的数据到: {final_h5ad}")
    adata_filtered.write(final_h5ad)
    
    print("\n分析完成! 结果保存到:", results_dir)


# 样例marker配置
human_marker_sample_config = """
[["T cells"]]
name = "T"
color = "#1ba169"
markers = ["CD3D", "CD3E", "CD3G", "CD8A", "TRBC2"]

[["B cells"]]
name = "B"
color = "#7ccbc5"
markers = ["CD79A", "CD79B", "MS4A1", "IGHM", "VPREB3"]

[["NK cells"]]
name = "NK"
color = "#014431"
markers = ["GNLY", "NKG7", "FGFBP2", "FCGR3A", "CX3CR1", "KLRB1", "NCR1"]

[["Myeloid cells"]]
name = "Myeloid"
color = "#1151b4"
markers = ["CD68", "CD14", "LYZ", "CD1E", "IL3RA"]
"""

mouse_marker_sample_config = """
[["T cells"]]
name = "T"
color = "#1ba169"
markers = ["Cd3d", "Cd3e", "Cd3g", "Cd8a", "Trbc2"]

[["B cells"]]
name = "B"
color = "#7ccbc5"
markers = ["Cd79a", "Cd79b", "Ms4a1", "Ighm", "Vpreb3"]

[["NK cells"]]
name = "NK"
color = "#014431"
markers = ["Klrb1", "Ncr1", "Nkg7"]

[["Myeloid cells"]]
name = "Myeloid"
color = "#1151b4"
markers = ["Cd68", "Cd14", "Lyz2", "Il3ra"]
"""

if __name__ == "__main__":
    main()