"""Dataset context helpers shared by scLucid workflows."""

from __future__ import annotations

from typing import Any, Dict, Literal, Optional, Union

from anndata import AnnData
from pydantic import ConfigDict, Field

from ..base_config import SclucidBaseConfig

DatasetType = Literal[
    "unknown",
    "pbmc_or_blood",
    "normal_tissue",
    "tumor_tissue",
    "cell_line",
    "organoid",
    "spatial",
]

_DATASET_TYPE_ALIASES = {
    "unknown": "unknown",
    "auto": "unknown",
    "pbmc": "pbmc_or_blood",
    "blood": "pbmc_or_blood",
    "immune": "pbmc_or_blood",
    "pbmc_or_blood": "pbmc_or_blood",
    "normal": "normal_tissue",
    "normal_tissue": "normal_tissue",
    "tissue": "normal_tissue",
    "tumor": "tumor_tissue",
    "tumour": "tumor_tissue",
    "cancer": "tumor_tissue",
    "tumor_tissue": "tumor_tissue",
    "tumour_tissue": "tumor_tissue",
    "cellline": "cell_line",
    "cell_line": "cell_line",
    "cell line": "cell_line",
    "organoid": "organoid",
    "spatial": "spatial",
    "visium": "spatial",
}

_MULTI_SAMPLE_HINTS = {"multi_sample", "multisample", "multi sample", "multiple_samples"}


def _normalize_token(value: Optional[str]) -> str:
    if value is None:
        return ""
    token = str(value).strip().lower().replace("-", "_")
    return " ".join(token.split())


def is_multi_sample_hint(value: Optional[str]) -> bool:
    """Return True when a string describes sample structure rather than biology."""
    token = _normalize_token(value)
    return token in _MULTI_SAMPLE_HINTS or token.replace(" ", "_") in _MULTI_SAMPLE_HINTS


def normalize_dataset_type(value: Optional[str]) -> DatasetType:
    """Normalize user-provided biological dataset type strings into canonical values."""
    token = _normalize_token(value)
    if not token:
        return "unknown"
    compact = token.replace(" ", "_")

    if compact in _DATASET_TYPE_ALIASES:
        return _DATASET_TYPE_ALIASES[compact]  # type: ignore[return-value]
    if token in _DATASET_TYPE_ALIASES:
        return _DATASET_TYPE_ALIASES[token]  # type: ignore[return-value]
    if is_multi_sample_hint(token):
        return "unknown"
    if "tumor" in token or "tumour" in token or "cancer" in token:
        return "tumor_tissue"
    if "pbmc" in token or "blood" in token:
        return "pbmc_or_blood"
    if "cell" in token and "line" in token:
        return "cell_line"
    if "organoid" in token:
        return "organoid"
    if "spatial" in token or "visium" in token:
        return "spatial"
    if "normal" in token:
        return "normal_tissue"
    return "unknown"


class AnalysisContext(SclucidBaseConfig):
    """Shared dataset context used to tune defaults without splitting workflows."""

    model_config = ConfigDict(extra="ignore")

    dataset_type: DatasetType = Field(default="unknown")
    species: str = Field(default="human")
    tissue: Optional[str] = Field(default=None)
    tissue_type: str = Field(default="unknown")
    cancer_type: Optional[str] = Field(default=None)
    sample_key: Optional[str] = Field(default=None)
    batch_key: Optional[str] = Field(default=None)
    condition_key: Optional[str] = Field(default=None)
    cell_type_key: Optional[str] = Field(default=None)
    is_spatial: bool = Field(default=False)
    is_multi_sample: bool = Field(default=False)
    n_cells: Optional[int] = Field(default=None)
    n_genes: Optional[int] = Field(default=None)
    notes: list[str] = Field(default_factory=list)

    @property
    def qc_tissue_type(self) -> str:
        """Return the legacy tissue_type value expected by QC recommenders."""
        if self.dataset_type == "tumor_tissue":
            return self.tissue_type if self.tissue_type != "unknown" else "tumor"
        if self.dataset_type == "spatial" and self._looks_tumor_context():
            return self.tissue_type if self.tissue_type != "unknown" else "tumor"
        if self.dataset_type in {"pbmc_or_blood", "normal_tissue", "cell_line", "organoid"}:
            return self.tissue_type if self.tissue_type != "unknown" else self.dataset_type
        return self.tissue_type or "unknown"

    @property
    def enables_tumor_module(self) -> bool:
        """Whether tumor-specific analysis is appropriate by default."""
        return self.dataset_type == "tumor_tissue" or (
            self.dataset_type in {"unknown", "spatial"} and self._looks_tumor_context()
        )

    def _looks_tumor_context(self) -> bool:
        if self.cancer_type:
            return True
        values = [self.tissue_type, self.tissue]
        return any(
            value is not None
            and any(token in str(value).lower() for token in ("tumor", "tumour", "cancer"))
            for value in values
        )

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serializable representation."""
        return self.model_dump()


DatasetProfile = AnalysisContext


def _first_existing_obs_key(adata: AnnData, candidates: list[str]) -> Optional[str]:
    for key in candidates:
        if key in adata.obs.columns:
            return key
    return None


def _single_obs_value(adata: AnnData, key: Optional[str]) -> Optional[str]:
    if key is None or key not in adata.obs.columns or adata.n_obs == 0:
        return None
    values = adata.obs[key].dropna().astype(str).unique().tolist()
    return values[0] if len(values) == 1 else None


def infer_analysis_context(
    adata: AnnData,
    *,
    context: Optional[Union[AnalysisContext, Dict[str, Any]]] = None,
    dataset_type: Optional[str] = None,
    species: str = "human",
    tissue: Optional[str] = None,
    tissue_type: str = "unknown",
    cancer_type: Optional[str] = None,
    sample_key: Optional[str] = None,
    batch_key: Optional[str] = None,
    condition_key: Optional[str] = None,
    cell_type_key: Optional[str] = None,
) -> AnalysisContext:
    """Infer a conservative analysis context from explicit hints and AnnData metadata."""
    if isinstance(context, AnalysisContext):
        base = context.model_copy(deep=True)
    elif isinstance(context, dict):
        base = AnalysisContext.model_validate(context)
    else:
        base = AnalysisContext()

    explicit_dataset_type = dataset_type or base.dataset_type
    explicit_multi_sample_hint = is_multi_sample_hint(dataset_type)
    obs_tissue_type = _single_obs_value(
        adata, _first_existing_obs_key(adata, ["dataset_type", "tissue_type", "sample_type"])
    )
    resolved_tissue_type = (
        tissue_type
        if tissue_type and tissue_type != "unknown"
        else base.tissue_type
        if base.tissue_type != "unknown"
        else obs_tissue_type
        or "unknown"
    )
    resolved_dataset_type = normalize_dataset_type(
        explicit_dataset_type if explicit_dataset_type != "unknown" else resolved_tissue_type
    )

    resolved_sample_key = sample_key or base.sample_key or _first_existing_obs_key(
        adata, ["sampleID", "sample", "Sample", "orig.ident", "orig_ident", "donor", "patient"]
    )
    resolved_batch_key = batch_key or base.batch_key or _first_existing_obs_key(
        adata, ["batch", "Batch", "sampleID", "sample", "orig.ident"]
    )
    resolved_condition_key = condition_key or base.condition_key or _first_existing_obs_key(
        adata, ["condition", "group", "treatment", "response", "disease", "phenotype"]
    )
    resolved_cell_type_key = cell_type_key or base.cell_type_key or _first_existing_obs_key(
        adata, ["cell_type_auto", "cell_type", "celltype", "annotation", "CellType"]
    )
    is_spatial = base.is_spatial or resolved_dataset_type == "spatial" or "spatial" in adata.obsm

    sample_n = (
        int(adata.obs[resolved_sample_key].nunique())
        if resolved_sample_key and resolved_sample_key in adata.obs.columns
        else 1
    )
    batch_n = (
        int(adata.obs[resolved_batch_key].nunique())
        if resolved_batch_key and resolved_batch_key in adata.obs.columns
        else 1
    )
    is_multi_sample = base.is_multi_sample or explicit_multi_sample_hint or max(sample_n, batch_n) > 1

    notes = list(base.notes)
    if explicit_multi_sample_hint:
        notes.append("Multi-sample was treated as sample structure, not dataset_type.")
    if resolved_dataset_type == "unknown":
        notes.append("Dataset type was not explicit; using conservative defaults.")
    if is_spatial and resolved_dataset_type != "spatial":
        notes.append("Spatial coordinates detected in adata.obsm['spatial'].")
    if is_multi_sample:
        notes.append("Multiple samples or batches detected.")

    return AnalysisContext(
        dataset_type=resolved_dataset_type,
        species=species or base.species,
        tissue=tissue if tissue is not None else base.tissue,
        tissue_type=resolved_tissue_type,
        cancer_type=cancer_type if cancer_type is not None else base.cancer_type,
        sample_key=resolved_sample_key,
        batch_key=resolved_batch_key,
        condition_key=resolved_condition_key,
        cell_type_key=resolved_cell_type_key,
        is_spatial=is_spatial,
        is_multi_sample=is_multi_sample,
        n_cells=int(adata.n_obs),
        n_genes=int(adata.n_vars),
        notes=list(dict.fromkeys(notes)),
    )


infer_dataset_profile = infer_analysis_context

__all__ = [
    "AnalysisContext",
    "DatasetProfile",
    "DatasetType",
    "infer_analysis_context",
    "infer_dataset_profile",
    "is_multi_sample_hint",
    "normalize_dataset_type",
]
