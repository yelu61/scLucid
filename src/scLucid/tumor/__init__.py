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
    construct_phylogeny,
    calculate_clonal_diversity,
    track_clonal_evolution,
)
from .cnv.cnv_signature import (
    extract_cnv_signatures,
    classify_cnv_pattern,
)

# Tumor Microenvironment
from .microenvironment.deconvolution import (
    deconvolve_tme,
    estimate_stromal_content,
    analyze_immune_infiltration,
)
from .microenvironment.interaction import (
    analyze_cell_interactions,
    calculate_communication_strength,
)
from .microenvironment.ecosystem import (
    score_tme_ecosystem,
    identify_ecosystem_types,
)

# Malignancy Analysis
from .malignancy.scoring import (
    score_malignancy,
    calculate_proliferation_index,
    estimate_metastatic_potential,
)
from .malignancy.classification import (
    classify_malignant_status,
    distinguish_tumor_normal,
)
from .malignancy.stemness import (
    calculate_stemness_index,
    identify_cancer_stem_cells,
)

# Tumor Evolution
from .evolution.phylogeny import (
    build_phylogenetic_tree,
    root_tree,
)
from .evolution.trajectory import (
    analyze_tumor_progression,
    identify_transition_states,
)
from .evolution.metastasis import (
    predict_metastasis_risk,
    analyze_dissemination,
)

# Therapy Response
from .therapy.resistance import (
    identify_resistance_mechanisms,
    score_drug_resistance,
)
from .therapy.prediction import (
    predict_therapy_response,
    stratify_patients,
)
from .therapy.target import (
    discover_therapeutic_targets,
    prioritize_druggable_genes,
)

# Tumor Heterogeneity
from .heterogeneity.diversity import (
    calculate_diversity_indices,
    estimate_intratumoral_heterogeneity,
)
from .heterogeneity.regional import (
    analyze_regional_heterogeneity,
    identify_spatial_patterns,
)
from .heterogeneity.temporal import (
    track_temporal_dynamics,
    analyze_treatment_response_trajectory,
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
    "construct_phylogeny",
    "calculate_clonal_diversity",
    "track_clonal_evolution",
    "extract_cnv_signatures",
    "classify_cnv_pattern",

    # Microenvironment
    "deconvolve_tme",
    "estimate_stromal_content",
    "analyze_immune_infiltration",
    "analyze_cell_interactions",
    "calculate_communication_strength",
    "score_tme_ecosystem",
    "identify_ecosystem_types",

    # Malignancy
    "score_malignancy",
    "calculate_proliferation_index",
    "estimate_metastatic_potential",
    "classify_malignant_status",
    "distinguish_tumor_normal",
    "calculate_stemness_index",
    "identify_cancer_stem_cells",

    # Evolution
    "build_phylogenetic_tree",
    "root_tree",
    "analyze_tumor_progression",
    "identify_transition_states",
    "predict_metastasis_risk",
    "analyze_dissemination",

    # Therapy
    "identify_resistance_mechanisms",
    "score_drug_resistance",
    "predict_therapy_response",
    "stratify_patients",
    "discover_therapeutic_targets",
    "prioritize_druggable_genes",

    # Heterogeneity
    "calculate_diversity_indices",
    "estimate_intratumoral_heterogeneity",
    "analyze_regional_heterogeneity",
    "identify_spatial_patterns",
    "track_temporal_dynamics",
    "analyze_treatment_response_trajectory",

    # Utils
    "get_tumor_markers",
    "get_immune_markers",
    "get_stromal_markers",
    "load_hallmark_signatures",
    "calculate_signature_scores",
    "query_cancer_gene_census",
    "get_drug_targets",
]
