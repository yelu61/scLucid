"""Dataset roles for scLucid validation benchmarks.

This registry is intentionally small and explicit. It maps local h5ad files to
the evidence claims they can support, so benchmark scripts can stay aligned with
``data/DATASETS.md`` and the roadmap phase documents.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DatasetSpec:
    key: str
    path: Path
    tissue: str
    modality_role: str
    qc_roles: tuple[str, ...]
    preprocess_roles: tuple[str, ...]
    figure2_panels: tuple[str, ...]
    figure3_panels: tuple[str, ...]
    required_obs: tuple[str, ...] = ("sample", "condition", "dataset")
    annotation_obs: tuple[str, ...] = ("cell_type", "cell_subtype")
    benchmark_notes: str = ""


DATASETS: tuple[DatasetSpec, ...] = (
    DatasetSpec(
        key="pbmc3k",
        path=Path("data/pbmc3k.h5ad"),
        tissue="PBMC",
        modality_role="normal_baseline",
        qc_roles=("fixed_threshold_baseline", "lightweight_smoke"),
        preprocess_roles=("layer_contract", "scanpy_parity", "hvg_baseline"),
        figure2_panels=("2A", "2B"),
        figure3_panels=("3A", "3B"),
        required_obs=("dataset",),
        annotation_obs=(),
        benchmark_notes="Fast baseline; legacy fixture lacks reliable sample/condition/cell_type annotations. Cell-type retention benchmarks degrade gracefully.",
    ),
    DatasetSpec(
        key="lin2020.pdac",
        path=Path("data/lin2020.pdac.h5ad"),
        tissue="PDAC",
        modality_role="tumor_baseline",
        qc_roles=("tumor_aware_qc", "sample_retention_bias"),
        preprocess_roles=("tumor_marker_preservation", "golden_path"),
        figure2_panels=("2B", "2C"),
        figure3_panels=("3A", "3B", "3C"),
        required_obs=("sample", "dataset"),
        benchmark_notes="Legacy PDAC fixture. cell_type/cell_subtype annotations are missing/unreliable; benchmarks that require them degrade to sample-level retention only.",
    ),
    DatasetSpec(
        key="schlesinger2020.pdac",
        path=Path("data/schlesinger2020.pdac.h5ad"),
        tissue="PDAC",
        modality_role="tumor_generalization",
        qc_roles=("tumor_aware_qc", "single_sample_behavior"),
        preprocess_roles=("tumor_marker_preservation",),
        figure2_panels=("2B", "2C"),
        figure3_panels=("3B",),
        required_obs=("dataset",),
        benchmark_notes="Legacy PDAC fixture. cell_type/cell_subtype annotations are missing/unreliable; benchmarks that require them degrade to sample-level retention only.",
    ),
    DatasetSpec(
        key="zilionis2019.nsclc",
        path=Path("data/zilionis2019.nsclc.h5ad"),
        tissue="NSCLC",
        modality_role="second_tumor_type",
        qc_roles=("tumor_aware_qc", "paired_tumor_blood", "cell_type_retention"),
        preprocess_roles=("marker_preservation", "tme_state_preservation"),
        figure2_panels=("2B", "2C"),
        figure3_panels=("3B", "3C"),
        benchmark_notes="Strong tumor-aware QC case with paired blood controls.",
    ),
    DatasetSpec(
        key="lee2020.crc",
        path=Path("data/lee2020.crc.h5ad"),
        tissue="CRC",
        modality_role="second_tumor_type",
        qc_roles=("tumor_aware_qc", "sample_retention_bias", "cell_type_retention"),
        preprocess_roles=("marker_preservation", "patient_integration_diagnostic"),
        figure2_panels=("2B", "2C"),
        figure3_panels=("3B", "3C"),
        benchmark_notes="Large tumor/normal cohort for retention bias and marker fidelity.",
    ),
    DatasetSpec(
        key="baron2016.pancreas",
        path=Path("data/baron2016.pancreas.h5ad"),
        tissue="pancreas",
        modality_role="normal_batch_diagnostic",
        qc_roles=("normal_reference",),
        preprocess_roles=("batch_diagnostic", "donor_structure", "hvg_baseline"),
        figure2_panels=("2B",),
        figure3_panels=("3A", "3B", "3C"),
        required_obs=("sample", "donor", "condition", "dataset"),
        benchmark_notes="Normal multi-donor dataset for opt-in integration diagnostics.",
    ),
    DatasetSpec(
        key="kang2018.pbmc",
        path=Path("data/kang2018.pbmc.h5ad"),
        tissue="PBMC",
        modality_role="doublet_ground_truth",
        qc_roles=("doublet_ground_truth", "perturbation_qc", "donor_multiplexing"),
        preprocess_roles=("stimulation_batch_diagnostic", "donor_structure"),
        figure2_panels=("2D",),
        figure3_panels=("3C",),
        required_obs=(
            "sample",
            "condition",
            "donor",
            "demuxlet_multiplets",
            "doublet_ground_truth",
            "dataset",
        ),
        benchmark_notes="Use singlet/doublet for metrics; report ambs separately.",
    ),
    DatasetSpec(
        key="cellbender_tiny",
        path=Path("data/cellbender_tiny.h5ad"),
        tissue="heart",
        modality_role="ambient_tiny_fixture",
        qc_roles=("ambient_rna", "empty_droplet"),
        preprocess_roles=(),
        figure2_panels=("2A", "2D"),
        figure3_panels=(),
        required_obs=(
            "sample",
            "condition",
            "likely_cell",
            "likely_empty_droplet",
            "dataset",
        ),
        annotation_obs=(),
        benchmark_notes="Tiny fixture for diagnostic contracts, not biological claims.",
    ),
)
