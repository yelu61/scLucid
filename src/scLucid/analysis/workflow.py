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

from analysis.config import AnalysisWorkflowConfig, ClusteringConfig, DifferentialConfig, AnnotationConfig
from analysis.pipeline import run_full_analysis # Assuming you create a new pipeline module

def setup_workflow_config() -> AnalysisWorkflowConfig:
    """Creates a complete configuration for the analysis."""
    
    # Define each step's configuration
    clustering_cfg = ClusteringConfig(resolution=1.2, key_added="leiden_1.2")
    
    de_cfg = DifferentialConfig(
        groupby="leiden_1.2", 
        method="wilcoxon"
    )

    annotation_cfg = AnnotationConfig(
        cluster_key="leiden_1.2",
        marker_species="human",
        final_method="combined",
        key_added="cell_type_auto"
    )

    # Combine into a master workflow config
    workflow_config = AnalysisWorkflowConfig(
        clustering=clustering_cfg,
        de=de_cfg,
        annotation=annotation_cfg
    )
    return workflow_config

def main():
    # 1. Load data
    adata = sc.read_h5ad("path/to/your/data.h5ad")
    
    # 2. Get the configuration
    config = setup_workflow_config()
    
    # 3. Run the entire analysis with one function call
    adata = run_full_analysis(adata, config)

    # 4. Save results
    adata.write_h5ad("path/to/results.h5ad")

if __name__ == "__main__":
    main()