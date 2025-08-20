
import logging
import pandas as pd
# TCGAbiolinks is a powerful package for this
try:
    import TCGAbiolinks
except ImportError:
    TCGAbiolinks = None

log = logging.getLogger(__name__)
__all__ = ["load_tcga_data"]

def load_tcga_data(
    project: str, # e.g., "TCGA-LUAD"
    data_category: str = "Transcriptome Profiling",
    data_type: str = "Gene Expression Quantification",
    workflow_type: str = "STAR - Counts",
    download: bool = True,
    cache_dir: str = "./tcga_data"
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Downloads and prepares TCGA RNA-seq count data and clinical metadata.

    Args:
        project: The TCGA project ID (e.g., 'TCGA-BRCA').
        ... (other parameters for GDC query) ...
        download: Whether to download the data if not cached.
        cache_dir: Directory to store downloaded data.

    Returns:
        A tuple of (counts_df, clinical_df).
    """
    if TCGAbiolinks is None:
        raise ImportError("Please install TCGAbiolinks: `pip install TCGAbiolinks`")
    
    log.info(f"Loading data for TCGA project: {project}")
    
    # 1. Query and download data from GDC
    query = TCGAbiolinks.GDCquery(
        project=project,
        data.category=data_category,
        data.type=data_type,
        workflow.type=workflow_type
    )
    
    if download:
        TCGAbiolinks.GDCdownload(query, directory=cache_dir)
    
    # 2. Prepare expression data
    data = TCGAbiolinks.GDCprepare(query, directory=cache_dir)
    counts_df = TCGAbiolinks.assay(data, "unstranded") # Get the raw counts
    
    # 3. Download clinical data
    clinical_df = TCGAbiolinks.GDCquery_clinic(project=project, type="clinical")
    
    # Clean up clinical data index to match counts data
    clinical_df.set_index("submitter_id", inplace=True)
    counts_df.columns = [ "-".join(x.split("-")[0:3]) for x in counts_df.columns ]
    
    # Align data
    common_samples = clinical_df.index.intersection(counts_df.columns)
    counts_df = counts_df[common_samples]
    clinical_df = clinical_df.loc[common_samples]

    log.info(f"Loaded {len(common_samples)} samples with both expression and clinical data.")
    return counts_df, clinical_df