import os
import pandas as pd
import anndata
import subprocess
import numpy as np

try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    R_AVAILABLE = True
except ImportError:
    R_AVAILABLE = False

def run_monocle3(
    adata: anndata.AnnData, 
    out_dir: str = "monocle3_tmp",
    r_path: str = "Rscript",
    root_cells: list = None
) -> anndata.AnnData:
    """
    在R中自动调用monocle3分析，并将pseudotime结果写回AnnData.obs
    """
    if not R_AVAILABLE:
        raise ImportError("请先安装rpy2，并确保R和monocle3已安装！")

    os.makedirs(out_dir, exist_ok=True)
    counts_csv = os.path.join(out_dir, "counts.csv")
    meta_csv = os.path.join(out_dir, "meta.csv")
    pseudo_csv = os.path.join(out_dir, "pseudotime.csv")

    # 1. 导出表达矩阵和meta
    pd.DataFrame(
        adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
        index=adata.obs_names,
        columns=adata.var_names
    ).T.to_csv(counts_csv)  # monocle3要求: gene x cell
    adata.obs.to_csv(meta_csv)

    # 2. 写R脚本
    r_script = os.path.join(out_dir, "run_monocle3.R")
    with open(r_script, "w") as f:
        f.write(f"""
library(monocle3)
library(Matrix)
library(data.table)
counts <- as.matrix(fread("{counts_csv}", data.table=FALSE), rownames=1)
meta <- read.csv("{meta_csv}", row.names=1)
cds <- new_cell_data_set(counts, cell_metadata=meta)
cds <- preprocess_cds(cds)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
""")
        if root_cells is not None:
            f.write(f'cds <- order_cells(cds, root_cells=c("{",".join(root_cells)}"))\n')
        else:
            f.write('cds <- order_cells(cds)\n')
        f.write(f'write.csv(pseudotime(cds), file="{pseudo_csv}")\n')

    # 3. 调用R
    cmd = f'{r_path} {r_script}'
    print(f"Running monocle3 R script: {cmd}")
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError("R monocle3分析失败，请检查R环境和依赖包。")

    # 4. 读取结果回AnnData
    pseudo = pd.read_csv(pseudo_csv, index_col=0)
    adata.obs["pseudotime"] = pseudo.loc[adata.obs_names, "x"].values
    return adata

def run_cellchat(
    adata: anndata.AnnData,
    groupby: str = "celltype",
    species: str = "Human",  # or "Mouse"
    out_dir: str = "cellchat_tmp",
    r_path: str = "Rscript",
    project_name: str = "CellChatProject"
) -> anndata.AnnData:
    """
    调用R的CellChat包进行细胞间通讯分析，并将通信强度矩阵写入AnnData.uns

    参数:
        adata: AnnData对象
        groupby: 细胞类型分组字段
        species: "Human"或"Mouse"
        out_dir: 临时目录
        r_path: Rscript路径（如环境变量已配置可用默认）
        project_name: CellChat工程名称
    返回值:
        结果写入adata.uns["cellchat"]
    """
    os.makedirs(out_dir, exist_ok=True)
    expr_csv = os.path.join(out_dir, "expr.csv")
    meta_csv = os.path.join(out_dir, "meta.csv")
    comm_csv = os.path.join(out_dir, "comm.csv")
    group_csv = os.path.join(out_dir, "group.csv")

    # 1. 导出表达矩阵（gene x cell）和meta信息
    expr_df = pd.DataFrame(
        adata.X.T.toarray() if hasattr(adata.X, "toarray") else adata.X.T,
        index=adata.var_names,
        columns=adata.obs_names
    )
    expr_df.to_csv(expr_csv)
    meta_df = adata.obs[[groupby]].copy()
    meta_df[groupby] = meta_df[groupby].astype(str)
    meta_df.to_csv(meta_csv)

    # 2. 生成R脚本
    r_script = os.path.join(out_dir, "run_cellchat.R")
    with open(r_script, "w") as f:
        f.write(f"""
library(CellChat)
library(patchwork)
library(data.table)
library(ggplot2)

data.input <- as.matrix(fread("{expr_csv}", data.table=FALSE), rownames=1)
meta <- read.csv("{meta_csv}", row.names=1)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "{groupby}")
cellchat@DB <- CellChatDB.{species}
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells=10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# 导出聚合后的通信强度（细胞类型-by-细胞类型矩阵）
write.csv(cellchat@net$weight, file="{comm_csv}")

# 导出分组信息
write.csv(meta, file="{group_csv}")

# 可选: 保存CellChat对象
saveRDS(cellchat, file="{out_dir}/cellchat_object.rds")
""")

    # 3. 调用R
    cmd = f'{r_path} {r_script}'
    print(f"Running CellChat R script: {cmd}")
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError("R CellChat分析失败，请检查R环境和依赖包。")

    # 4. 读取通信矩阵结果写入 AnnData.uns
    comm_df = pd.read_csv(comm_csv, index_col=0)
    adata.uns["cellchat"] = {
        "communication_matrix": comm_df,
        "groupby": groupby,
        "species": species
    }
    return adata


def run_scenic(
    adata: anndata.AnnData,
    species: str = "hs",  # "hs" for human, "mm" for mouse
    out_dir: str = "scenic_tmp",
    r_path: str = "Rscript",
    scenic_db_dir: str = "/path/to/cistarget_databases",
    project_name: str = "SCENIC"
) -> anndata.AnnData:
    """
    调用R的SCENIC包进行转录因子活性推断，并将AUC矩阵写入AnnData.obsm

    参数:
        adata: AnnData对象
        species: "hs" (human) 或 "mm" (mouse)
        out_dir: 临时目录
        r_path: Rscript路径
        scenic_db_dir: SCENIC数据库目录
        project_name: SCENIC工程名称
    返回值:
        结果写入adata.obsm["AUC"]
    """
    os.makedirs(out_dir, exist_ok=True)
    expr_csv = os.path.join(out_dir, "expr.csv")
    auc_csv = os.path.join(out_dir, "auc.csv")

    # 1. 导出表达矩阵（gene x cell）
    expr_df = pd.DataFrame(
        adata.X.T.toarray() if hasattr(adata.X, "toarray") else adata.X.T,
        index=adata.var_names,
        columns=adata.obs_names
    )
    expr_df.to_csv(expr_csv)

    # 2. 生成 R 脚本
    r_script = os.path.join(out_dir, "run_scenic.R")
    with open(r_script, "w") as f:
        f.write(f"""
library(SCENIC)
library(data.table)
library(SingleCellExperiment)

exprMat <- as.matrix(fread("{expr_csv}", data.table=FALSE), rownames=1)
cellInfo <- data.frame(CellType=rep("cell", ncol(exprMat)))
rownames(cellInfo) <- colnames(exprMat)
org <- "{species}"
dbDir <- "{scenic_db_dir}"
dbs <- list.files(dbDir, full.names=TRUE, pattern="{species}")
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle="{project_name}", nCores=4)

# SCENIC主流程
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
runGenie3(exprMat_filtered, scenicOptions)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

# AUC矩阵导出
aucell <- loadInt(scenicOptions, "aucell")
auc <- getAUC(aucell)
write.csv(auc, file="{auc_csv}")
""")

    # 3. 调用R
    cmd = f'{r_path} {r_script}'
    print(f"Running SCENIC R script: {cmd}")
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError("R SCENIC分析失败，请检查R环境和依赖包。")

    # 4. 读取AUC结果写入 AnnData.obsm
    auc_df = pd.read_csv(auc_csv, index_col=0)
    auc_df = auc_df.T  # 转为cell x regulon
    adata.obsm["AUC"] = auc_df.loc[adata.obs_names].values  # 顺序对齐
    adata.uns["AUC_regulon_names"] = list(auc_df.columns)
    return adata


def run_cellphonedb(
    adata: anndata.AnnData,
    groupby: str = "celltype",
    out_dir: str = "cellphonedb_tmp",
    cellphonedb_cmd: str = "cellphonedb",
    project_name: str = "cpdb_project"
) -> anndata.AnnData:
    """
    Run CellPhoneDB analysis and import key results back into AnnData.uns.

    Args:
        adata: AnnData object.
        groupby: Column in .obs for cluster/celltype.
        out_dir: Directory for temporary files and results.
        cellphonedb_cmd: Command for CellPhoneDB (default: 'cellphonedb').
        project_name: Used for output naming.

    Returns:
        adata: AnnData object with CellPhoneDB results in .uns['cellphonedb'].
    """
    os.makedirs(out_dir, exist_ok=True)
    meta_file = os.path.join(out_dir, "meta.txt")
    counts_file = os.path.join(out_dir, "counts.txt")

    # 1. Export meta.txt (cell & cluster/celltype)
    meta = adata.obs[[groupby]].copy()
    meta.reset_index(inplace=True)
    meta.columns = ["Cell", "cell_type"]
    meta.to_csv(meta_file, sep="\t", index=False)

    # 2. Export counts.txt (gene x cell, raw counts, as required by CPDB)
    # CPDB expects genes as rows, cells as columns, integers
    if hasattr(adata.raw, "X"):
        X = adata.raw.X
        var_names = adata.raw.var_names
    else:
        X = adata.X
        var_names = adata.var_names

    counts = pd.DataFrame(
        X.T.toarray() if hasattr(X, "toarray") else X.T,
        index=var_names,
        columns=adata.obs_names
    )
    counts = counts.round().astype(int)
    counts.to_csv(counts_file, sep="\t")

    # 3. Run CellPhoneDB (subprocess call)
    cwd = os.getcwd()
    os.chdir(out_dir)
    try:
        # (1) Initialize (optional)
        # subprocess.run([cellphonedb_cmd, "database", "generate"], check=True)
        # (2) Analysis
        result = subprocess.run([
            cellphonedb_cmd, "method", "statistical_analysis",
            "meta.txt", "counts.txt",
            "--output-path", ".", "--project-name", project_name,
            "--threads", "4"
        ], check=True, capture_output=True, text=True)
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(e.stdout)
        print(e.stderr)
        raise RuntimeError("CellPhoneDB analysis failed. See output above for details.")
    finally:
        os.chdir(cwd)

    # 4. Load results into AnnData.uns
    result_files = {
        "means": os.path.join(out_dir, "means.txt"),
        "pvalues": os.path.join(out_dir, "pvalues.txt"),
        "deconvoluted": os.path.join(out_dir, "deconvoluted.txt"),
        "significant_means": os.path.join(out_dir, "significant_means.txt"),
    }
    adata.uns["cellphonedb"] = {}
    for key, path in result_files.items():
        if os.path.exists(path):
            adata.uns["cellphonedb"][key] = pd.read_csv(path, sep="\t", index_col=0)
    return adata

import os
import pandas as pd
import anndata
import subprocess

def run_slingshot(
    adata: anndata.AnnData,
    embedding_key: str = "X_umap",  # 可以是PCA、UMAP、TSNE等
    cluster_key: str = "leiden",    # 用于分群，决定拟时序的起点/分支
    start_cluster: str = None,      # 可选，指定起始cluster
    out_dir: str = "slingshot_tmp",
    r_path: str = "Rscript"
) -> anndata.AnnData:
    """
    调用R的slingshot包进行拟时序分析，将pseudotime结果写入AnnData.obs

    参数:
        adata: AnnData对象
        embedding_key: obsm中的降维坐标key
        cluster_key: obs中的分群信息key
        start_cluster: 拟时序起点（可选）
        out_dir: 临时目录
        r_path: Rscript路径
    返回:
        adata.obs['slingshot_pseudotime'] 增加pseudotime
        adata.obs['slingshot_lineage']    增加主分支/轨迹编号
    """
    os.makedirs(out_dir, exist_ok=True)
    embed_csv = os.path.join(out_dir, "embedding.csv")
    cluster_csv = os.path.join(out_dir, "cluster.csv")
    pseudo_csv = os.path.join(out_dir, "pseudotime.csv")
    lineage_csv = os.path.join(out_dir, "lineage.csv")

    # 1. 导出降维坐标和分群信息
    pd.DataFrame(
        adata.obsm[embedding_key],
        index=adata.obs_names,
        columns=[f"Dim_{i+1}" for i in range(adata.obsm[embedding_key].shape[1])]
    ).to_csv(embed_csv)
    pd.DataFrame({
        "cluster": adata.obs[cluster_key].astype(str)
    }, index=adata.obs_names).to_csv(cluster_csv)

    # 2. 生成R脚本
    r_script = os.path.join(out_dir, "run_slingshot.R")
    with open(r_script, "w") as f:
        f.write(f"""
library(slingshot)
library(data.table)
embed <- as.matrix(fread("{embed_csv}", data.table=FALSE, row.names=1))
cluster <- read.csv("{cluster_csv}", row.names=1)$cluster
""")
        if start_cluster is not None:
            f.write(f'ss <- slingshot(embed, cluster, start.clus="{start_cluster}")\n')
        else:
            f.write('ss <- slingshot(embed, cluster)\n')
        f.write(f"""
# pseudotime 和 lineage（主分支）导出
pt <- slingPseudotime(ss)
write.csv(pt, file="{pseudo_csv}")
lg <- slingClusterLabels(ss)
write.csv(lg, file="{lineage_csv}")
""")

    # 3. 调用R
    cmd = f'{r_path} {r_script}'
    print(f"Running Slingshot R script: {cmd}")
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError("R slingshot分析失败，请检查R和依赖包。")

    # 4. 读取结果写回 AnnData
    pt = pd.read_csv(pseudo_csv, index_col=0)
    lg = pd.read_csv(lineage_csv, index_col=0)
    # 默认取第一条主分支
    adata.obs["slingshot_pseudotime"] = pt.iloc[:, 0].reindex(adata.obs_names)
    adata.obs["slingshot_lineage"] = lg.iloc[:, 0].reindex(adata.obs_names)
    return adata