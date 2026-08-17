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

CELL_TYPE_KEY_CANDIDATES = (
    "cell_type",
    "cell_type_auto",
    "cell_type_final",
    "celltype",
    "cell_type_major",
    "annotation",
    "cell_annotation",
)

CELL_LINEAGE_KEY_CANDIDATES = (
    "celltype_lineage",
    "celltype_lineage_auto",
    "cell_lineage",
    "lineage",
)


def _normalize_token(value: Optional[str]) -> str:
    if value is None:
        return ""
    token = str(value).strip().lower().replace("-", "_")
    return " ".join(token.split())


def is_multi_sample_hint(value: Optional[str]) -> bool:
    """Return True when a string describes sample structure rather than biology."""
    token = _normalize_token(value)
    return token in _MULTI_SAMPLE_HINTS or token.replace(" ", "_") in _MULTI_SAMPLE_HINTS


def is_tumor_context(*values: Optional[str]) -> bool:
    """Return True when any provided context string indicates tumor biology."""
    tumor_tokens = ("tumor", "tumour", "cancer", "malignan")
    for value in values:
        token = _normalize_token(value)
        if token and any(marker in token for marker in tumor_tokens):
            return True
    return False


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
    if is_tumor_context(token):
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
    """Shared project context used to tune defaults and expose design assumptions.

    The original dataset-oriented fields remain stable.  The additional study
    design fields make the biological replicate and paired structure explicit
    before scLucid recommends integration or sample-level inference.
    """

    model_config = ConfigDict(extra="ignore")

    dataset_type: DatasetType = Field(default="unknown")
    species: str = Field(default="human")
    assay: str = Field(default="scrna")
    tissue: Optional[str] = Field(default=None)
    tissue_type: str = Field(default="unknown")
    cancer_type: Optional[str] = Field(default=None)
    study_objective: Optional[str] = Field(default=None)
    sample_key: Optional[str] = Field(default=None)
    batch_key: Optional[str] = Field(default=None)
    condition_key: Optional[str] = Field(default=None)
    experimental_unit_key: Optional[str] = Field(default=None)
    paired_key: Optional[str] = Field(default=None)
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
        return is_tumor_context(self.tissue_type, self.tissue)

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serializable representation."""
        return self.model_dump()


DatasetProfile = AnalysisContext
ProjectContext = AnalysisContext


def _first_existing_obs_key(adata: AnnData, candidates: list[str]) -> Optional[str]:
    for key in candidates:
        if key in adata.obs.columns:
            return key
    return None


def resolve_obs_key(
    adata: AnnData,
    candidates: tuple[str, ...] | list[str],
    *,
    preferred: Optional[str] = None,
) -> Optional[str]:
    """Return the first available obs key using a canonical candidate order."""
    ordered: list[str] = []
    if preferred:
        ordered.append(preferred)
    ordered.extend(str(key) for key in candidates if key and key not in ordered)
    return _first_existing_obs_key(adata, ordered)


def resolve_cell_type_key(adata: AnnData, preferred: Optional[str] = None) -> Optional[str]:
    """Resolve the canonical downstream cell-type annotation key."""
    return resolve_obs_key(adata, CELL_TYPE_KEY_CANDIDATES, preferred=preferred)


def resolve_cell_lineage_key(adata: AnnData, preferred: Optional[str] = None) -> Optional[str]:
    """Resolve the canonical downstream cell-lineage annotation key."""
    return resolve_obs_key(adata, CELL_LINEAGE_KEY_CANDIDATES, preferred=preferred)


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
    species: Optional[str] = None,
    assay: Optional[str] = None,
    tissue: Optional[str] = None,
    tissue_type: str = "unknown",
    cancer_type: Optional[str] = None,
    study_objective: Optional[str] = None,
    sample_key: Optional[str] = None,
    batch_key: Optional[str] = None,
    condition_key: Optional[str] = None,
    experimental_unit_key: Optional[str] = None,
    paired_key: Optional[str] = None,
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
        adata,
        [
            "batch",
            "Batch",
            "library_batch",
            "sequencing_batch",
            "seq_batch",
            "run_id",
        ],
    )
    resolved_condition_key = condition_key or base.condition_key or _first_existing_obs_key(
        adata, ["condition", "group", "treatment", "response", "disease", "phenotype"]
    )
    resolved_experimental_unit_key = (
        experimental_unit_key
        or base.experimental_unit_key
        or _first_existing_obs_key(
            adata,
            ["patient_id", "patient", "donor_id", "donor", "subject_id", "subject"],
        )
        or resolved_sample_key
    )
    resolved_paired_key = paired_key or base.paired_key
    if (
        resolved_paired_key is None
        and resolved_condition_key
        and resolved_experimental_unit_key
        and resolved_condition_key in adata.obs.columns
        and resolved_experimental_unit_key in adata.obs.columns
    ):
        condition_counts = adata.obs.groupby(
            resolved_experimental_unit_key, observed=True
        )[resolved_condition_key].nunique()
        if bool((condition_counts > 1).any()):
            resolved_paired_key = resolved_experimental_unit_key
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
    if resolved_experimental_unit_key:
        notes.append(
            f"Experimental unit defaults to obs[{resolved_experimental_unit_key!r}]; "
            "verify this before formal inference."
        )
    if resolved_paired_key:
        notes.append(
            f"Repeated conditions were detected within obs[{resolved_paired_key!r}]."
        )

    return AnalysisContext(
        dataset_type=resolved_dataset_type,
        species=species or base.species,
        assay=assay or base.assay,
        tissue=tissue if tissue is not None else base.tissue,
        tissue_type=resolved_tissue_type,
        cancer_type=cancer_type if cancer_type is not None else base.cancer_type,
        study_objective=(
            study_objective if study_objective is not None else base.study_objective
        ),
        sample_key=resolved_sample_key,
        batch_key=resolved_batch_key,
        condition_key=resolved_condition_key,
        experimental_unit_key=resolved_experimental_unit_key,
        paired_key=resolved_paired_key,
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
    "ProjectContext",
    "infer_analysis_context",
    "infer_dataset_profile",
    "is_multi_sample_hint",
    "is_tumor_context",
    "normalize_dataset_type",
]
