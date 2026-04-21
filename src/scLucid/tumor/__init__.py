"""
scLucid Tumor Module - Comprehensive toolkit for cancer single-cell analysis

This module provides specialized tools for tumor research, including:
- Copy number variation (CNV) analysis
- Tumor microenvironment characterization
- Malignant cell identification
- Tumor evolution tracking
- Therapy response prediction
- Tumor heterogeneity quantification

Example:
    >>> from scLucid.tumor import cnv, microenvironment, malignancy
    >>>
    >>> # CNV analysis
    >>> adata = cnv.infer_cnv(adata, reference_cells='immune')
    >>>
    >>> # Identify malignant cells
    >>> adata = malignancy.score_malignancy(adata)
    >>>
    >>> # Deconvolve TME
    >>> proportions = microenvironment.deconvolve_tme(adata)
"""

__version__ = "0.1.0"
__author__ = "scLucid"

# CNV Analysis
from .cnv.infercnv import (
    infer_cnv,
    find_tumor_cells,
    identify_clones,
    calculate_cnv_score,
)
from .cnv.clone_analysis import (
    infer_clonal_phylogeny,
    calculate_clonal_diversity,
    identify_clones as identify_clones_from_cnv,
    CloneAnalyzer,
)
from .cnv.cnv_signature import (
    extract_cnv_signatures,
    assign_cnv_signature,
    CNVSigExtractor,
)

# Tumor Microenvironment
from .microenvironment.deconvolution import (
    deconvolve_tme,
    estimate_stromal_content,
    analyze_immune_infiltration,
)
from .microenvironment.interaction import (
    analyze_cell_interactions,
    find_dominant_interactions,
    score_immune_interactions,
)
from .microenvironment.ecosystem import (
    analyze_ecosystem_composition,
    calculate_tumor_microenvironment_score,
    compare_ecosystems,
)

# Malignancy Analysis
from .malignancy.scoring import (
    score_malignancy,
    calculate_proliferation_index,
    estimate_metastatic_potential,
)
from .malignancy.classification import (
    classify_malignant_cells,
    score_malignancy_potential,
)
from .malignancy.stemness import (
    calculate_stemness_score,
    identify_cancer_stem_cells,
    compare_stemness_between_groups,
)

# Tumor Evolution
from .evolution.phylogeny import (
    build_phylogenetic_tree,
    root_tree,
)
from .evolution.trajectory import (
    analyze_tumor_progression,
    identify_transition_states,
    align_progression_trajectories,
)
from .evolution.metastasis import (
    predict_metastasis_risk,
    analyze_dissemination,
    compare_primary_vs_metastasis,
)

# Therapy Response
from .therapy.resistance import (
    identify_resistance_mechanisms,
    score_drug_resistance,
    compare_resistance_between_groups,
)
from .therapy.prediction import (
    predict_therapy_response,
    stratify_patients,
    evaluate_biomarker,
)
from .therapy.target import (
    discover_therapeutic_targets,
    prioritize_druggable_genes,
    suggest_targeted_therapies,
)

# Tumor Heterogeneity
from .heterogeneity.diversity import (
    calculate_diversity_indices,
    estimate_intratumoral_heterogeneity,
    calculate_transcriptional_diversity,
)
from .heterogeneity.regional import (
    analyze_regional_heterogeneity,
    identify_spatial_patterns,
    calculate_regional_expression_differences,
)
from .heterogeneity.temporal import (
    track_temporal_dynamics,
    analyze_treatment_response_trajectory,
    detect_clonal_sweep,
)

# Utilities
from .utils.markers import (
    get_tumor_markers,
    get_immune_markers,
    get_stromal_markers,
)
from .utils.signatures import (
    load_hallmark_signatures,
    calculate_signature_scores,
)
from .utils.databases import (
    query_cancer_gene_census,
    get_drug_targets,
)

__all__ = [
    # CNV
    "infer_cnv",
    "find_tumor_cells",
    "identify_clones",
    "calculate_cnv_score",
    "infer_clonal_phylogeny",
    "calculate_clonal_diversity",
    "identify_clones_from_cnv",
    "CloneAnalyzer",
    "extract_cnv_signatures",
    "assign_cnv_signature",
    "CNVSigExtractor",

    # Microenvironment
    "deconvolve_tme",
    "estimate_stromal_content",
    "analyze_immune_infiltration",
    "analyze_cell_interactions",
    "find_dominant_interactions",
    "score_immune_interactions",
    "analyze_ecosystem_composition",
    "calculate_tumor_microenvironment_score",
    "compare_ecosystems",

    # Malignancy
    "score_malignancy",
    "calculate_proliferation_index",
    "estimate_metastatic_potential",
    "classify_malignant_cells",
    "score_malignancy_potential",
    "calculate_stemness_score",
    "identify_cancer_stem_cells",
    "compare_stemness_between_groups",

    # Evolution
    "build_phylogenetic_tree",
    "root_tree",
    "analyze_tumor_progression",
    "identify_transition_states",
    "align_progression_trajectories",
    "predict_metastasis_risk",
    "analyze_dissemination",
    "compare_primary_vs_metastasis",

    # Therapy
    "identify_resistance_mechanisms",
    "score_drug_resistance",
    "compare_resistance_between_groups",
    "predict_therapy_response",
    "stratify_patients",
    "evaluate_biomarker",
    "discover_therapeutic_targets",
    "prioritize_druggable_genes",
    "suggest_targeted_therapies",

    # Heterogeneity
    "calculate_diversity_indices",
    "estimate_intratumoral_heterogeneity",
    "calculate_transcriptional_diversity",
    "analyze_regional_heterogeneity",
    "identify_spatial_patterns",
    "calculate_regional_expression_differences",
    "track_temporal_dynamics",
    "analyze_treatment_response_trajectory",
    "detect_clonal_sweep",

    # Utils
    "get_tumor_markers",
    "get_immune_markers",
    "get_stromal_markers",
    "load_hallmark_signatures",
    "calculate_signature_scores",
    "query_cancer_gene_census",
    "get_drug_targets",
]
