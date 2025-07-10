"""
Workflow module for single-cell RNA-seq data analysis.
"""

from typing import Optional, List, Dict, Union
import scanpy as sc
import anndata as ad
from . import qc, norm, hvg, integrate, config
from .config import load_config

class scRNAWorkflow:
    
    def __init__(self, config_path: Optional[str] = None):
        """Initialize Workflow Manager"""
        self.config = config.load_config(config_path)
        self.results = {}  # Storing analysis results
        
    def run_workflow(
        self, 
        adata: ad.AnnData,
        steps: List[str] = ["qc", "norm", "hvg", "dim_reduction", "clustering"],
        sample_key: str = "sampleID",
    ) -> ad.AnnData:
        """Run a standard workflow"""
        self.adata = adata.copy()
        
        for step in steps:
            if step == "qc":
                self._run_qc(sample_key)
            elif step == "norm":
                self._run_normalization()
            elif step == "hvg":
                self._run_hvg(sample_key)
            elif step == "dim_reduction":
                self._run_dim_reduction()
            elif step == "clustering":
                self._run_clustering()
            # Other steps...
        
        return self.adata
    
    def _run_qc(self, sample_key: str):
        """Run quality control steps"""
        qc_config = self.config["qc"]
        self.adata = qc.calculate_qc_metric(
            self.adata, 
            sample_key=sample_key,
            plot_violin=qc_config.get("plot_violin", True)
        )
        self.adata = qc.is_low_quality_cell(
            self.adata,
            sample_key=sample_key,
            min_genes=qc_config.get("min_genes", 200),
            pc_mt=qc_config.get("pc_mt", 20)
        )
        # Filter low quality cells
        self.adata = self.adata[~self.adata.obs.filter(regex="outlier").any(axis=1)].copy()
    
    # Implementing other steps...