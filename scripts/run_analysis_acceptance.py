#!/usr/bin/env python
"""Run the evidence-first analysis acceptance workflow.

This runner is the scriptable counterpart of
``examples/03_advanced_notebooks/Step2-Annotation_and_Malignancy.ipynb``. It
expects a preprocessed AnnData object and writes a compact set of artifacts that
make clustering, annotation, and optional malignancy interpretation reviewable.
"""

from __future__ import annotations

import argparse
import json
import shutil
import time
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

import scLucid as scl
from scLucid.analysis import (
    AnalysisWorkflowConfig,
    AnnotationConfig,
    ClusteringConfig,
    DifferentialConfig,
    run_standard_analysis,
    summarize_analysis_review_summary,
    validate_analysis_module_completeness,
)
from scLucid.utils import sanitize_for_hdf5

DEFAULT_INPUT_PATH = Path("data/processed/Step2-sce_preprocessed.h5ad")
DEFAULT_OUTPUT_DIR = Path("results/analysis_acceptance")
DEFAULT_RESOLUTIONS = (0.5, 0.6, 0.8)


def _package_version(package: str) -> str:
    try:
        return version(package)
    except PackageNotFoundError:
        return "unknown"


def make_json_safe(value: Any) -> Any:
    """Convert common scientific Python objects to JSON-safe containers."""
    if isinstance(value, pd.DataFrame):
        return value.to_dict(orient="records")
    if isinstance(value, pd.Series):
        return value.to_dict()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, dict):
        return {str(k): make_json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [make_json_safe(v) for v in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def prepare_adata_for_write(adata: AnnData) -> None:
    """Make known scLucid metadata writable without losing review summaries."""
    if "sclucid" in adata.uns:
        adata.uns["sclucid"] = sanitize_for_hdf5(adata.uns["sclucid"])


def load_preprocessed_adata(
    input_path: Path,
    *,
    n_cells: int | None = None,
    random_state: int = 42,
) -> tuple[AnnData, dict[str, int]]:
    """Load a preprocessed object and optionally subsample cells."""
    adata = sc.read_h5ad(str(input_path))
    original_shape = {"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars)}
    if n_cells is not None and adata.n_obs > n_cells:
        rng = np.random.default_rng(random_state)
        idx = rng.choice(adata.n_obs, size=int(n_cells), replace=False)
        adata = adata[idx].copy()
    if adata.raw is None:
        adata.raw = adata
    return adata, original_shape


def build_analysis_acceptance_config(
    output_dir: Path,
    *,
    resolutions: tuple[float, ...] = DEFAULT_RESOLUTIONS,
    marker_config: str | None = None,
    species: str = "human",
    tissue: str | None = None,
    run_celltypist: bool = False,
    celltypist_model: str = "Immune_All_Low.pkl",
    run_malignancy: bool = False,
    run_cnv: bool = False,
    cancer_type: str | None = None,
) -> AnalysisWorkflowConfig:
    """Build the standard Step2 evidence-first analysis config."""
    annotation_methods = ["marker_manager", "data_driven"]
    if run_celltypist:
        annotation_methods.insert(0, "celltypist")

    return AnalysisWorkflowConfig(
        save_dir=str(output_dir / "analysis"),
        n_jobs=1,
        clustering=ClusteringConfig(
            method="leiden",
            resolution=resolutions[0],
            use_rep="X_pca",
            key_added="leiden_clusters",
            plot=False,
        ),
        de=DifferentialConfig(
            groupby="leiden_clusters",
            method="wilcoxon",
            use_raw=True,
            key_added="rank_genes_groups",
        ),
        annotation=AnnotationConfig(
            cluster_key="leiden_clusters",
            marker_species=species,
            marker_tissue=tissue,
            lineage_marker_config=marker_config,
            run_celltypist=run_celltypist,
            celltypist_model=celltypist_model,
            run_scoring=False,
            final_method="celltypist" if run_celltypist else "combined",
            key_added="cell_type_auto",
            lineage_key="celltype_lineage_auto",
        ),
        run_clustering_review=True,
        candidate_resolutions=list(resolutions),
        use_recommended_resolution=True,
        run_annotation_evidence=True,
        annotation_methods=tuple(annotation_methods),
        final_annotation_strategy="consensus",
        annotation_level="lineage",
        run_malignancy_interpretation=run_malignancy,
        run_cnv_for_malignancy=run_cnv,
        run_malignancy_score=run_malignancy,
        malignancy_cancer_type=cancer_type,
        characterize=False,
    )


def export_analysis_artifacts(adata: AnnData, output_dir: Path) -> dict[str, str]:
    """Write reviewable analysis artifacts and return their paths."""
    output_dir.mkdir(parents=True, exist_ok=True)
    analysis_ns = adata.uns.get("sclucid", {}).get("analysis", {})
    annotation_ns = analysis_ns.get("annotation", {}) if isinstance(analysis_ns, dict) else {}
    malignancy_ns = analysis_ns.get("malignancy", {}) if isinstance(analysis_ns, dict) else {}

    artifacts: dict[str, str] = {}

    review_table = annotation_ns.get("annotation_review_table")
    if isinstance(review_table, pd.DataFrame):
        path = output_dir / "annotation_review_table.csv"
        review_table.to_csv(path, index=False)
        artifacts["annotation_review_table"] = str(path)

    marker_table = annotation_ns.get("marker_annotation_evidence")
    if isinstance(marker_table, pd.DataFrame):
        path = output_dir / "marker_annotation_evidence.csv"
        marker_table.to_csv(path, index=False)
        artifacts["marker_annotation_evidence"] = str(path)

    llm_bundle = annotation_ns.get("llm_annotation_bundle")
    if isinstance(llm_bundle, dict):
        path = output_dir / "llm_annotation_bundle.json"
        path.write_text(json.dumps(make_json_safe(llm_bundle), indent=2), encoding="utf-8")
        artifacts["llm_annotation_bundle"] = str(path)

    malignancy_table = malignancy_ns.get("malignancy_interpretation_table")
    if isinstance(malignancy_table, pd.DataFrame):
        path = output_dir / "malignancy_interpretation_table.csv"
        malignancy_table.to_csv(path, index=False)
        artifacts["malignancy_interpretation_table"] = str(path)

    review_summary = analysis_ns.get("review_summary")
    if isinstance(review_summary, dict):
        path = output_dir / "analysis_review_summary_compact.json"
        compact = summarize_analysis_review_summary(review_summary)
        path.write_text(json.dumps(make_json_safe(compact), indent=2), encoding="utf-8")
        artifacts["analysis_review_summary_compact"] = str(path)

    return artifacts


def build_analysis_acceptance(
    adata: AnnData,
    *,
    workflow_name: str,
    input_shape: dict[str, int],
    original_shape: dict[str, int],
    elapsed_seconds: float,
) -> dict[str, Any]:
    """Build a compact analysis acceptance report."""
    analysis_ns = adata.uns.get("sclucid", {}).get("analysis", {})
    review_summary = analysis_ns.get("review_summary", {}) if isinstance(analysis_ns, dict) else {}
    validation = validate_analysis_module_completeness(adata)
    compact = summarize_analysis_review_summary(review_summary) if review_summary else {}
    malignancy = (
        review_summary.get("malignancy_interpretation_summary", {})
        if isinstance(review_summary, dict)
        else {}
    )
    metrics = {
        "n_clusters": compact.get("n_clusters"),
        "review_table_rows": compact.get("review_table_rows"),
        "needs_review_clusters": compact.get("needs_review_clusters"),
        "n_final_labels": compact.get("n_final_labels"),
        "mean_annotation_confidence": compact.get("mean_confidence"),
        "malignancy_enabled": malignancy.get("enabled"),
        "n_malignant": malignancy.get("n_malignant"),
        "n_suspect_malignant": malignancy.get("n_suspect_malignant"),
        "analysis_valid": validation["valid"],
        "readiness_status": compact.get("readiness_status"),
        "maturity_status": compact.get("maturity_status"),
    }
    blocking_failures = []
    if not validation["valid"]:
        blocking_failures.extend(validation.get("issues", []))
    if not metrics["n_clusters"] or int(metrics["n_clusters"]) < 2:
        blocking_failures.append("analysis_acceptance.n_clusters<2")
    if not metrics["review_table_rows"]:
        blocking_failures.append("analysis_acceptance.annotation_review_table_missing")
    if not metrics["n_final_labels"]:
        blocking_failures.append("analysis_acceptance.final_labels_missing")

    return {
        "schema_version": "0.1",
        "scope": "analysis_acceptance",
        "workflow_name": workflow_name,
        "ready_for_real_data_review": not blocking_failures,
        "claim_boundary": (
            "This lightweight acceptance report checks traceability and workflow "
            "readiness. It does not claim scientific superiority over standard "
            "single-cell workflows."
        ),
        "original_shape": original_shape,
        "input_shape": input_shape,
        "final_shape": {"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars)},
        "metrics": metrics,
        "blocking_failures": list(dict.fromkeys(map(str, blocking_failures))),
        "warnings": validation.get("warnings", []),
        "elapsed_seconds": float(elapsed_seconds),
        "versions": {
            "sclucid": getattr(scl, "__version__", "unknown"),
            "scanpy": _package_version("scanpy"),
            "anndata": _package_version("anndata"),
            "numpy": _package_version("numpy"),
        },
    }


def run_analysis_acceptance(
    *,
    input_path: Path = DEFAULT_INPUT_PATH,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    n_cells: int | None = None,
    resolutions: tuple[float, ...] = DEFAULT_RESOLUTIONS,
    marker_config: str | None = None,
    species: str = "human",
    tissue: str | None = None,
    run_celltypist: bool = False,
    celltypist_model: str = "Immune_All_Low.pkl",
    run_malignancy: bool = False,
    run_cnv: bool = False,
    cancer_type: str | None = None,
    random_state: int = 42,
    overwrite: bool = False,
    show_progress: bool = False,
    write_h5ad: bool = True,
) -> dict[str, Any]:
    """Run analysis acceptance and return the machine-readable manifest."""
    input_path = input_path.resolve()
    output_dir = output_dir.resolve()
    if output_dir.exists() and overwrite:
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    started = time.perf_counter()
    adata, original_shape = load_preprocessed_adata(
        input_path,
        n_cells=n_cells,
        random_state=random_state,
    )
    input_shape = {"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars)}
    config = build_analysis_acceptance_config(
        output_dir,
        resolutions=resolutions,
        marker_config=marker_config,
        species=species,
        tissue=tissue,
        run_celltypist=run_celltypist,
        celltypist_model=celltypist_model,
        run_malignancy=run_malignancy,
        run_cnv=run_cnv,
        cancer_type=cancer_type,
    )
    steps = [
        "clustering_review",
        "clustering",
        "markers",
        "annotation_evidence",
        "annotation_consensus",
    ]
    if run_celltypist:
        steps.insert(3, "annotation")
    if run_malignancy:
        steps.append("malignancy_interpretation")

    adata = run_standard_analysis(
        adata,
        config=config,
        steps=steps,
        show_progress=show_progress,
    )
    elapsed = time.perf_counter() - started

    analysis_artifacts = export_analysis_artifacts(adata, output_dir / "analysis")
    acceptance = build_analysis_acceptance(
        adata,
        workflow_name="analysis_acceptance",
        input_shape=input_shape,
        original_shape=original_shape,
        elapsed_seconds=elapsed,
    )
    validation_dir = output_dir / "validation"
    validation_dir.mkdir(parents=True, exist_ok=True)
    acceptance_path = validation_dir / "analysis_acceptance.json"
    acceptance_path.write_text(
        json.dumps(make_json_safe(acceptance), indent=2),
        encoding="utf-8",
    )

    final_h5ad = output_dir / "Step3-sce_annotated.h5ad"
    if write_h5ad:
        prepare_adata_for_write(adata)
        adata.write(final_h5ad, compression="gzip")

    manifest = {
        "workflow": "analysis_acceptance",
        "schema_version": "1.0",
        "input_path": str(input_path),
        "output_dir": str(output_dir),
        "subset_n_cells": n_cells,
        "random_state": random_state,
        "resolutions": list(resolutions),
        "steps": steps,
        "acceptance": acceptance,
        "artifacts": {
            "acceptance_json": str(acceptance_path),
            "analysis_dir": str(output_dir / "analysis"),
            "analysis": analysis_artifacts,
            "final_h5ad": str(final_h5ad) if write_h5ad else None,
        },
    }
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(make_json_safe(manifest), indent=2), encoding="utf-8")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-path", type=Path, default=DEFAULT_INPUT_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--n-cells", type=int, default=None)
    parser.add_argument(
        "--resolutions",
        type=float,
        nargs="+",
        default=list(DEFAULT_RESOLUTIONS),
        help="Candidate Leiden resolutions for clustering review.",
    )
    parser.add_argument("--marker-config", type=str, default=None)
    parser.add_argument("--species", type=str, default="human")
    parser.add_argument("--tissue", type=str, default=None)
    parser.add_argument("--run-celltypist", action="store_true")
    parser.add_argument("--celltypist-model", type=str, default="Immune_All_Low.pkl")
    parser.add_argument("--run-malignancy", action="store_true")
    parser.add_argument("--run-cnv", action="store_true")
    parser.add_argument("--cancer-type", type=str, default=None)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--show-progress", action="store_true")
    parser.add_argument("--skip-h5ad", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = run_analysis_acceptance(
        input_path=args.input_path,
        output_dir=args.output_dir,
        n_cells=args.n_cells,
        resolutions=tuple(args.resolutions),
        marker_config=args.marker_config,
        species=args.species,
        tissue=args.tissue,
        run_celltypist=args.run_celltypist,
        celltypist_model=args.celltypist_model,
        run_malignancy=args.run_malignancy,
        run_cnv=args.run_cnv,
        cancer_type=args.cancer_type,
        random_state=args.random_state,
        overwrite=args.overwrite,
        show_progress=args.show_progress,
        write_h5ad=not args.skip_h5ad,
    )
    print(json.dumps(make_json_safe(manifest), indent=2))


if __name__ == "__main__":
    main()
