"""Pydantic configurations for bulk RNA-seq analysis."""

from __future__ import annotations

from typing import List, Literal, Optional

from pydantic import Field

from ...base_config import SclucidBaseConfig


class BulkDiagnosticsConfig(SclucidBaseConfig):
    """Configuration for bulk data quality diagnostics."""

    model_config = SclucidBaseConfig.model_config

    min_samples_total: int = Field(default=2, ge=2)
    min_samples_per_condition: int = Field(default=2, ge=1)
    require_replicates: bool = Field(default=True)
    max_library_size_cv: Optional[float] = Field(default=None, ge=0)
    max_zero_gene_fraction: Optional[float] = Field(default=None, ge=0, le=1)


class BulkNormalizationConfig(SclucidBaseConfig):
    """Configuration for bulk RNA-seq normalization."""

    model_config = SclucidBaseConfig.model_config

    method: Literal["CPM", "TPM", "FPKM", "RPKM", "deseq2_size_factors"] = Field(
        default="CPM"
    )
    gene_length_col: Optional[str] = Field(
        default=None,
        description="Column in adata.var with gene lengths (bp) for TPM/FPKM/RPKM.",
    )
    target_sum: float = Field(default=1e6, gt=0)
    pseudocount: float = Field(default=0.0, ge=0)


class BulkDEConfig(SclucidBaseConfig):
    """Configuration for bulk differential expression."""

    model_config = SclucidBaseConfig.model_config

    condition_col: str = Field(description="Column in adata.obs defining conditions.")
    condition1: str = Field(description="Reference condition.")
    condition2: str = Field(description="Test condition.")
    method: Literal["ttest", "welch", "pydeseq2", "limma"] = Field(default="welch")
    covariates: List[str] = Field(default_factory=list)
    min_counts_per_gene: int = Field(default=10, ge=0)
    min_samples_expressing: int = Field(default=2, ge=1)
    p_adjust_method: Literal["fdr_bh", "bonferroni"] = Field(default="fdr_bh")
    alpha: float = Field(default=0.05, gt=0, lt=1)
    fallback_to_descriptive: bool = Field(
        default=False,
        description="If True and replicates are insufficient, return effect sizes without p-values.",
    )


class BulkTraitAssociationConfig(SclucidBaseConfig):
    """Configuration for continuous trait association on bulk samples."""

    model_config = SclucidBaseConfig.model_config

    trait_col: str = Field(description="Column in adata.obs with continuous trait.")
    method: Literal["pearson", "spearman", "ols"] = Field(default="spearman")
    covariates: List[str] = Field(default_factory=list)
    min_samples: int = Field(default=10, ge=5)


class BulkDeconvolutionConfig(SclucidBaseConfig):
    """Configuration for bulk deconvolution."""

    model_config = SclucidBaseConfig.model_config

    method: Literal["BayesPrism", "DWLS"] = Field(default="BayesPrism")
    cell_type_key: str = Field(description="Column in reference obs with cell type labels.")
    sample_key: str = Field(default="sampleID")
    min_common_genes: int = Field(default=100, ge=10)


class BulkAbundanceConfig(SclucidBaseConfig):
    """Configuration for differential abundance of deconvolved proportions."""

    model_config = SclucidBaseConfig.model_config

    group_col: str = Field(description="Column in metadata defining groups.")
    group1: str = Field(description="Reference group.")
    group2: str = Field(description="Test group.")
    method: Literal["ttest", "wilcoxon"] = Field(default="wilcoxon")
    p_adjust_method: Literal["fdr_bh", "bonferroni"] = Field(default="fdr_bh")


class BulkClinicalAssociationConfig(SclucidBaseConfig):
    """Configuration for correlating abundance with a clinical variable."""

    model_config = SclucidBaseConfig.model_config

    clinical_variable: str = Field(description="Column in metadata with continuous variable.")
    method: Literal["pearson", "spearman"] = Field(default="spearman")
    min_samples: int = Field(default=10, ge=5)
