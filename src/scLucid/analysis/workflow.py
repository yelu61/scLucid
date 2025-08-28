"""
workflow.py
------------

Recommended main workflow for single-cell RNA-seq analysis.
Covers clustering, marker analysis, enrichment, annotation (auto/manual/AI), scoring, visualization, and evaluation.

Assumes all functions/modules are imported from analysis package:
- de_enrichment.py
- annotation.py
- scoring.py
"""

import scanpy as sc
from analysis.de_enrichment import (
    find_markers, filter_markers, run_enrichment, summarize_markers_and_enrichment,
    visualize_markers
)
from analysis.annotation import (
    score_cell_types, annotate_clusters, run_celltypist, transfer_labels,
    summarize_annotation_evidence, apply_annotation_mapping, evaluate_annotation
)
from analysis.scoring import (
    score_by_gene_sets, compare_scores, batch_compare_scores, plot_score_comparison
)

# 1. 预处理、降维、聚类（可用scanpy主流程）
def preprocess_and_cluster(adata, n_pcs=30, n_neighbors=15, resolution=1.0, key_added="leiden"):
    sc.pp.normalize_total(adata); sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata); adata = adata[:, adata.var.highly_variable]
    sc.pp.pca(adata, n_comps=n_pcs); sc.pp.neighbors(adata, n_neighbors=n_neighbors)
    sc.tl.umap(adata); sc.tl.leiden(adata, resolution=resolution, key_added=key_added)
    return adata

# 2. 差异分析与marker筛选
def marker_analysis(adata, cluster_key="leiden", method="wilcoxon"):
    markers_df = find_markers(adata, groupby=cluster_key, method=method)
    filtered_df = filter_markers(adata, key="rank_genes_groups")
    return filtered_df

# 3. 富集分析
def enrichment_analysis(adata, cluster_key="leiden", markers_df=None):
    enrichment_dict = run_enrichment(adata, groupby=cluster_key)
    return enrichment_dict

# 4. 导出摘要供AI/人工命名
def export_annotation_summary(adata, markers_df, enrichment_dict, out_file="annotation_summary.md"):
    return summarize_markers_and_enrichment(adata, markers_df, enrichment_dict, summary_file=out_file)

# 5. 应用人工/AI命名mapping
def apply_manual_annotation(adata, cluster_key, mapping_file, key_added="cell_type_ai"):
    apply_annotation_mapping(adata, cluster_key, mapping_file, key_added=key_added)

# 6. 自动注释（如需）
def auto_annotation(adata, marker_config, cluster_key="leiden", method="combined"):
    score_cell_types(adata, marker_config)
    annotate_clusters(adata, cluster_key, marker_config, method=method, key_added="cell_type_auto")

# 7. 注释评估
def annotation_evaluation(adata, cluster_key, annotation_key, marker_config, plot=True):
    return evaluate_annotation(adata, cluster_key, annotation_key, marker_config, plot=plot)

# 8. 分数/通路打分及对比
def score_and_compare(adata, gene_sets, groupby="leiden"):
    score_by_gene_sets(adata, gene_sets)
    for score in gene_sets:
        plot_score_comparison(adata, f"{score}_score", groupby=groupby)

# 9. 可视化markers
def bulk_marker_visualization(adata, markers_df, cluster_key="leiden"):
    visualize_markers(adata, markers_df, groupby=cluster_key, plot_type="dotplot")

# 10. 全流程执行示例（入口）
def main_workflow(adata, marker_config, gene_sets, mapping_file=None):
    adata = preprocess_and_cluster(adata)
    filtered_markers = marker_analysis(adata)
    enrichment_dict = enrichment_analysis(adata, markers_df=filtered_markers)
    export_annotation_summary(adata, filtered_markers, enrichment_dict)
    if mapping_file:
        apply_manual_annotation(adata, "leiden", mapping_file)
        annotation_evaluation(adata, "leiden", "cell_type_ai", marker_config)
    else:
        auto_annotation(adata, marker_config, "leiden")
        annotation_evaluation(adata, "leiden", "cell_type_auto", marker_config)
    score_and_compare(adata, gene_sets)
    bulk_marker_visualization(adata, filtered_markers)

# 11. 可选：团队调用
# main_workflow(adata, marker_config="markers.toml", gene_sets={"Tcell": ["CD3D", "CD3E"]}, mapping_file="mapping.csv")