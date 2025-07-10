import numpy as np
import pandas as pd
import scanpy as sc
import random
#import torch
from scipy.stats import median_abs_deviation
import os
from typing import List, Dict

def identify_outliers(adata: sc.AnnData, metric: str, nmads: int) -> pd.Series:
    """
    Identify outliers based on the given metric and number of median absolute deviations.

    Args:
        adata (AnnData): AnnData object to check for outliers.
        metric (str): The metric to use for outlier detection. Must be a valid column in adata.obs.
        nmads (int): Number of median absolute deviations for outlier detection.

    Returns:
        outliers (pandas.Series): Boolean mask indicating if a cell is an outlier or not.
    """
    if metric not in adata.obs.columns:
        raise ValueError(f"Invalid metric '{metric}'. Must be a column in adata.obs.")

    values = adata.obs[metric]
    median = np.median(values)
    mad = median_abs_deviation(values)
    outliers = [abs(value - median) > nmads * mad for value in values]
    return pd.Series(outliers, index=adata.obs_names)


def merge_data(all_samples: List[str], raw_dir: str, group_dict: Dict[str, str]) -> sc.AnnData:
    """
    Merge single-cell data from multiple samples into one AnnData object.

    Args:
        all_samples (List[str]): List of sample names.
        raw_dir (str): Directory containing raw data files.
        group_dict (Dict[str, str]): Dictionary mapping sample names to group labels.

    Returns:
        adata (AnnData): Merged AnnData object.
    """
    if not all_samples:
        raise ValueError("No samples provided.")

    if not os.path.isdir(raw_dir):
        raise ValueError(f"Invalid raw data directory: {raw_dir}")

    adata_list = []
    for sample in all_samples:
        sample_path = os.path.join(raw_dir, sample)
        if os.path.isdir(sample_path):
            adata = sc.read_10x_mtx(
                sample_path,
                var_names="gene_symbols",
                make_unique=True,
                cache=False,
                gex_only=True,
                prefix=None,
            )
            adata.obs["sampleID"] = sample
            adata.obs["group"] = group_dict.get(sample, "Unknown")
            adata_list.append(adata)

    if not adata_list:
        raise ValueError("No valid sample data found in the provided directory.")

    adata = sc.concat(adata_list, join="outer", index_unique="-")
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    return adata

def seed_everything(seed):
    np.random.seed(seed)
    random.seed(seed)
    #torch.manual_seed(seed)
    #torch.cuda.manual_seed_all(seed)
    #torch.backends.cudnn.deterministic = True
    #torch.backends.cudnn.benchmark = False
    

def save_cell_counts_to_csv(adata, obs_column, column_name, output_file):
    # 根据指定的 obs 列统计细胞数量
    cell_counts = adata.obs[obs_column].value_counts()
    # 将 Series 转换为 DataFrame
    cell_counts_df = cell_counts.to_frame(name=column_name)
    # 重置索引并添加列名
    cell_counts_df = cell_counts_df.reset_index()
    cell_counts_df.columns = [obs_column, column_name]
    # 保存为 CSV 文件
    cell_counts_df.to_csv(output_file, index=False)
    print(f"细胞数量统计已保存至 {output_file}")
