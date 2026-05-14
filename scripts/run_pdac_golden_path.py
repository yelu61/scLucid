#!/usr/bin/env python
"""Run the PDAC tumor golden-path workflow for scLucid.

This script runs a traceable tumor baseline workflow on the Lin 2020 PDAC dataset,
writes reviewable artifacts, and emits a JSON manifest. It exercises:
  QC -> Preprocessing -> Analysis -> Tumor (malignancy + CNV)
"""

from __future__ import annotations

import argparse
import json
import shutil
import time
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

import scLucid as scl
from scLucid.analysis import (
    AnalysisWorkflowConfig,
    AnnotationConfig,
    ClusteringConfig,
    DifferentialConfig,
)
from scLucid.preprocess import WorkflowConfig
from scLucid.qc import QCWorkflowConfig
from scLucid.tumor.config import TumorAnalysisConfig
from scLucid.utils import (
    build_qc_preprocess_validation,
    validate_all_stage_contracts,
    write_validation_outputs,
)


DEFAULT_DATA_PATH = Path("data/lin2020.pdac.h5ad")
DEFAULT_OUTPUT_DIR = Path("results/golden/pdac")


def _package_version(package: str) -> str:
    try:
        return version(package)
    except PackageNotFoundError:
        return "unknown"


def load_pdac_data(
    data_path: Path,
    *,
    n_cells: int | None = None,
    random_state: int = 42,
):
    """Load PDAC data, optionally subsample, and ensure raw counts are in place."""
    adata = sc.read_h5ad(str(data_path))
    original_shape = {"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars)}

    if n_cells is not None and adata.n_obs > n_cells:
        rng = np.random.default_rng(random_state)
        idx = rng.choice(adata.n_obs, size=n_cells, replace=False)
        adata = adata[idx].copy()

    # Ensure raw counts are available
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    else:
        adata.layers["counts"] = adata.X.copy()

    # Add required metadata if missing
    if "species" not in adata.obs.columns:
        adata.obs["species"] = "human"
    if "sampleID" not in adata.obs.columns:
        # Use orig.ident or study as fallback
        if "orig.ident" in adata.obs.columns:
            adata.obs["sampleID"] = adata.obs["orig.ident"]
        else:
            adata.obs["sampleID"] = "pdac"

    return adata, original_shape


def make_hdf5_safe(value: Any) -> Any:
    """Convert nested workflow metadata to forms AnnData can write to HDF5."""
    if isinstance(value, pd.DataFrame):
        return {
            "format": "dataframe_records",
            "index": {str(index): str(item) for index, item in enumerate(value.index)},
            "columns": {str(index): str(item) for index, item in enumerate(value.columns)},
            "records": {
                str(index): make_hdf5_safe(record)
                for index, record in enumerate(value.to_dict(orient="records"))
            },
        }
    if isinstance(value, pd.Series):
        return {
            "format": "series",
            "index": {str(index): str(item) for index, item in enumerate(value.index)},
            "values": {
                str(index): make_hdf5_safe(item)
                for index, item in enumerate(value.tolist())
            },
        }
    if isinstance(value, np.ndarray):
        if value.dtype == object:
            return make_hdf5_safe(value.tolist())
        return value
    if isinstance(value, tuple):
        return {str(index): make_hdf5_safe(item) for index, item in enumerate(value)}
    if isinstance(value, list):
        return {str(index): make_hdf5_safe(item) for index, item in enumerate(value)}
    if isinstance(value, dict):
        return {str(key): make_hdf5_safe(item) for key, item in value.items()}
    if isinstance(value, (str, int, float, bool, np.number, np.bool_)):
        return value
    if value is None:
        return "null"
    return str(value)


def compact_sclucid_metadata(value: Any) -> dict[str, Any]:
    """Create compact provenance instead of embedding heavy workflow evidence."""
    if not isinstance(value, dict):
        return {"type": type(value).__name__}

    compact: dict[str, Any] = {
        "metadata_policy": (
            "Full workflow evidence is kept in review artifacts and the run "
            "manifest; the final h5ad stores review summaries and stage keys."
        )
    }
    for stage, payload in value.items():
        if isinstance(payload, dict):
            stage_payload: dict[str, Any] = {
                "_stage_keys": sorted(map(str, payload.keys()))
            }
            if "review_summary" in payload:
                stage_payload["review_summary"] = payload["review_summary"]
            compact[str(stage)] = stage_payload
        else:
            compact[str(stage)] = {"type": type(payload).__name__}

    return compact


def prepare_adata_for_write(adata) -> None:
    """Clean known workflow metadata that can otherwise block ``.h5ad`` writes."""
    if "sclucid" in adata.uns:
        adata.uns["sclucid"] = make_hdf5_safe(
            compact_sclucid_metadata(adata.uns["sclucid"])
        )


def build_configs(
    output_dir: Path,
    *,
    n_top_genes: int = 2000,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    include_annotation: bool = True,
) -> tuple[QCWorkflowConfig, WorkflowConfig, AnalysisWorkflowConfig, TumorAnalysisConfig]:
    """Build configs for the PDAC tumor golden path."""
    # Tumor-aware QC: more permissive mitochondrial threshold
    qc_config = QCWorkflowConfig(
        save_dir=str(output_dir / "qc"),
        species="human",
        tissue_type="tumor_tissue",
        use_parallel=False,
        n_jobs=1,
    )
    qc_config.metrics_reporting_config.plot_violin = False
    qc_config.metrics_reporting_config.plot_scatter = False
    qc_config.metrics_reporting_config.plot_top_genes = False
    qc_config.metrics_reporting_config.show_plots = False
    qc_config.metrics_reporting_config.export_stats = True
    qc_config.metrics_reporting_config.export_xlsx = False
    qc_config.marking_config.plot_outliers = False
    qc_config.marking_config.show_plots = False
    qc_config.doublet_config.plot_summary = False
    qc_config.doublet_config.plot_bar = False
    qc_config.doublet_config.plot_scatter = False
    qc_config.doublet_config.plot_upset = False
    qc_config.doublet_config.show_plots = False
    qc_config.doublet_config.export_stats = True

    preprocess_config = WorkflowConfig.quick(
        n_top_genes=n_top_genes,
        run_regression=False,
        run_integration=False,
        save_dir=str(output_dir / "preprocess"),
        n_jobs=1,
    )
    preprocess_config.normalization.plot = False
    preprocess_config.normalization.report = False
    preprocess_config.plot = False
    preprocess_config.graph.n_pcs = n_pcs
    preprocess_config.graph.n_neighbors = n_neighbors

    analysis_config = AnalysisWorkflowConfig(
        save_dir=str(output_dir / "analysis"), n_jobs=1
    )
    analysis_config.clustering = ClusteringConfig(
        method="leiden",
        resolution=0.8,  # Slightly higher for tumor heterogeneity
        use_rep="X_pca",
        key_added="leiden_clusters",
        plot=False,
    )
    analysis_config.de = DifferentialConfig(
        groupby="leiden_clusters",
        method="wilcoxon",
        use_raw=True,
        key_added="rank_genes_groups",
    )
    analysis_config.annotation = (
        AnnotationConfig(
            cluster_key="leiden_clusters",
            marker_species="human",
            run_celltypist=False,
            run_scoring=True,
            final_method="combined",
            key_added="cell_type_auto",
        )
        if include_annotation
        else None
    )

    tumor_config = TumorAnalysisConfig(
        run_malignancy=True,
        malignancy_method="cnv",
        run_tme=False,  # Skip TME for now — needs cell_type_auto
        run_cnv=True,
        run_therapy=False,
    )

    return qc_config, preprocess_config, analysis_config, tumor_config


def save_embedding_figures(adata, figures_dir: Path) -> list[str]:
    """Save a small inspection figure set and return artifact paths."""
    figures_dir.mkdir(parents=True, exist_ok=True)
    artifacts: list[str] = []

    figure_specs = [
        ("leiden_clusters", "umap_leiden_clusters.png"),
        ("cell_type_auto", "umap_cell_type_auto.png"),
        ("sampleID", "umap_sampleID.png"),
    ]
    for color_by, filename in figure_specs:
        if color_by not in adata.obs.columns or "X_umap" not in adata.obsm:
            continue
        try:
            fig = scl.pl.plot_embedding(
                adata,
                color_by=color_by,
                basis="umap",
                show=False,
                legend_style="legend",
                show_labels=color_by != "sampleID",
            )
            path = figures_dir / filename
            fig.savefig(path, dpi=200, bbox_inches="tight")
            plt.close(fig)
            artifacts.append(str(path))
        except Exception as exc:
            print(f"Warning: could not save {filename}: {exc}")

    return artifacts


def run_tumor_stage(adata, tumor_config: TumorAnalysisConfig):
    """Run tumor-specific analysis on pre-processed data."""
    import logging

    log = logging.getLogger(__name__)
    executed_steps = []
    stage_warnings = []

    # Malignancy scoring
    if tumor_config.run_malignancy:
        try:
            from scLucid.tumor.malignancy.scoring import score_malignancy

            log.info("Tumor stage: scoring malignancy")
            adata = score_malignancy(adata, key_added="malignancy")
            executed_steps.append("malignancy_scoring")
        except Exception as exc:
            msg = f"Malignancy scoring failed: {exc}"
            log.warning(msg)
            stage_warnings.append(msg)

    # Malignancy classification
    if tumor_config.run_malignancy:
        try:
            from scLucid.tumor.malignancy.classification import classify_malignant_cells

            log.info("Tumor stage: classifying malignant cells")
            classify_malignant_cells(adata, method="cnv", key_added="is_malignant")
            executed_steps.append("malignancy_classification")
        except Exception as exc:
            msg = f"Malignancy classification failed: {exc}"
            log.warning(msg)
            stage_warnings.append(msg)

    # CNV inference
    if tumor_config.run_cnv:
        try:
            from scLucid.tumor.cnv.infercnv import infer_cnv

            log.info("Tumor stage: inferring CNV")
            # Try to use cell_type_auto to identify reference normal cells
            ref_key = None
            if "cell_type_auto" in adata.obs.columns:
                ref_key = "cell_type_auto"
                # Find cells annotated as normal/immune/stromal
                normal_types = ["B cells", "T cells", "NK cells", "Monocytes",
                                "Macrophages", "Dendritic cells", "Fibroblasts",
                                "Endothelial cells"]
                obs_lower = adata.obs[ref_key].astype(str).str.lower()
                ref_mask = obs_lower.str.contains(
                    "b cell|t cell|nk |monocyte|macrophage|dendritic|fibroblast|endothelial|immune|stromal"
                )
                if ref_mask.any():
                    ref_values = adata.obs.loc[ref_mask, ref_key].unique().tolist()
                    infer_cnv(adata, reference_cells=ref_values, reference_key=ref_key, key_added="cnv")
                else:
                    infer_cnv(adata, key_added="cnv")
            else:
                infer_cnv(adata, key_added="cnv")
            executed_steps.append("cnv_inference")
        except Exception as exc:
            msg = f"CNV inference failed: {exc}"
            log.warning(msg)
            stage_warnings.append(msg)

    return adata, executed_steps, stage_warnings


def write_manifest(
    *,
    manifest_path: Path,
    data_path: Path,
    output_dir: Path,
    original_shape: dict[str, int],
    input_shape: dict[str, int],
    adata,
    elapsed_seconds: float,
    figure_artifacts: list[str],
    tumor_steps: list[str],
    tumor_warnings: list[str],
    subset_n_cells: int | None,
    random_state: int,
) -> dict[str, Any]:
    """Write a machine-readable run manifest."""
    contract_results = validate_all_stage_contracts(adata, when="output")
    obs_summary = {
        "n_clusters": (
            int(adata.obs["leiden_clusters"].nunique())
            if "leiden_clusters" in adata.obs.columns
            else None
        ),
        "n_cell_types": (
            int(adata.obs["cell_type_auto"].nunique())
            if "cell_type_auto" in adata.obs.columns
            else None
        ),
        "samples": (
            sorted(map(str, adata.obs["sampleID"].unique().tolist()))
            if "sampleID" in adata.obs.columns
            else []
        ),
        "n_malignant": (
            int(adata.obs["is_malignant"].sum())
            if "is_malignant" in adata.obs.columns
            else None
        ),
    }
    manifest: dict[str, Any] = {
        "workflow": "pdac_golden_path",
        "schema_version": "1.0",
        "data_path": str(data_path),
        "output_dir": str(output_dir),
        "subset_n_cells": subset_n_cells,
        "random_state": random_state,
        "original_shape": original_shape,
        "input_shape": input_shape,
        "final_shape": {"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars)},
        "retention_fraction": float(adata.n_obs / input_shape["n_cells"]),
        "obs_summary": obs_summary,
        "versions": {
            "sclucid": getattr(scl, "__version__", "unknown"),
            "scanpy": _package_version("scanpy"),
            "anndata": _package_version("anndata"),
            "numpy": _package_version("numpy"),
        },
        "contracts": {
            stage: result.to_dict() for stage, result in contract_results.items()
        },
        "tumor": {
            "steps_executed": tumor_steps,
            "warnings": tumor_warnings,
        },
        "artifacts": {
            "final_h5ad": str(output_dir / "pdac_golden_final.h5ad"),
            "manifest": str(manifest_path),
            "figures": figure_artifacts,
            "qc_dir": str(output_dir / "qc"),
            "preprocess_dir": str(output_dir / "preprocess"),
            "analysis_dir": str(output_dir / "analysis"),
        },
        "elapsed_seconds": elapsed_seconds,
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, default=str), encoding="utf-8"
    )
    return manifest


def run_pdac_golden_path(
    *,
    data_path: Path = DEFAULT_DATA_PATH,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    n_cells: int | None = None,
    random_state: int = 42,
    n_top_genes: int = 2000,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    include_annotation: bool = True,
    overwrite: bool = False,
    show_progress: bool = False,
) -> dict[str, Any]:
    """Run the PDAC golden path and return the manifest dictionary."""
    data_path = data_path.resolve()
    output_dir = output_dir.resolve()

    if output_dir.exists() and overwrite:
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    started = time.perf_counter()
    adata, original_shape = load_pdac_data(
        data_path, n_cells=n_cells, random_state=random_state
    )
    input_shape = {"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars)}

    qc_config, preprocess_config, analysis_config, tumor_config = build_configs(
        output_dir,
        n_top_genes=n_top_genes,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        include_annotation=include_annotation,
    )

    analysis_skip_steps = [] if include_annotation else ["annotation"]

    # Run core pipeline
    adata = scl.run_pipeline(
        adata,
        stages=["qc", "preprocess", "analysis"],
        dataset_type="tumor_tissue",
        species="human",
        qc_config=qc_config,
        preprocess_config=preprocess_config,
        analysis_config=analysis_config,
        show_progress=show_progress,
        analysis_skip_steps=analysis_skip_steps,
    )

    # Run tumor stage
    adata, tumor_steps, tumor_warnings = run_tumor_stage(adata, tumor_config)

    figure_artifacts = save_embedding_figures(adata, output_dir / "figures")
    elapsed = time.perf_counter() - started
    validation = build_qc_preprocess_validation(
        adata,
        run_manifest={
            "workflow": "pdac_golden_path",
            "input_shape": input_shape,
            "retention_fraction": float(adata.n_obs / input_shape["n_cells"]),
        },
        dataset_role="pdac_tumor",
        workflow_name="pdac_golden_path",
    )
    validation_artifacts = write_validation_outputs(validation, output_dir / "validation")
    prepare_adata_for_write(adata)
    final_h5ad = output_dir / "pdac_golden_final.h5ad"
    adata.write(final_h5ad)

    manifest = write_manifest(
        manifest_path=output_dir / "manifest.json",
        data_path=data_path,
        output_dir=output_dir,
        original_shape=original_shape,
        input_shape=input_shape,
        adata=adata,
        elapsed_seconds=elapsed,
        figure_artifacts=figure_artifacts,
        tumor_steps=tumor_steps,
        tumor_warnings=tumor_warnings,
        subset_n_cells=n_cells,
        random_state=random_state,
    )
    manifest["validation"] = validation
    manifest["artifacts"]["validation"] = validation_artifacts
    Path(manifest["artifacts"]["manifest"]).write_text(
        json.dumps(manifest, indent=2, default=str), encoding="utf-8"
    )

    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-path", type=Path, default=DEFAULT_DATA_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--n-cells", type=int, default=None)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--skip-annotation", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--show-progress", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = run_pdac_golden_path(
        data_path=args.data_path,
        output_dir=args.output_dir,
        n_cells=args.n_cells,
        random_state=args.random_state,
        n_top_genes=args.n_top_genes,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        include_annotation=not args.skip_annotation,
        overwrite=args.overwrite,
        show_progress=args.show_progress,
    )
    print(json.dumps(manifest, indent=2, default=str))


if __name__ == "__main__":
    main()
