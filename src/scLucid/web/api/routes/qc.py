"""
QC analysis API routes for Dash backend.

Provides endpoints for:
- Running QC workflows
- Calculating QC metrics
- Generating QC reports
- Data callbacks for Dash interactive components
"""

from typing import Optional

import numpy as np
import pandas as pd
from fastapi import APIRouter, HTTPException

from ....base_config import SclucidBaseConfig

router = APIRouter()


class QCConfig(SclucidBaseConfig):
    """QC configuration model for web API."""

    sample_key: str = "sampleID"
    species: str = "human"
    min_genes: Optional[int] = None
    max_mt_percent: Optional[float] = None
    min_counts: Optional[int] = None
    doublet_detection: bool = True
    cell_cycle_scoring: bool = False


@router.post("/metrics")
async def calculate_qc_metrics(project_id: str):
    """
    Calculate QC metrics for a project.

    Args:
        project_id: Project identifier

    Returns:
        Calculated metrics summary for Plotly.js visualization
    """
    import scLucid as scl
    from scLucid.web.services.data_manager import get_project_data

    adata = get_project_data(project_id)
    if adata is None:
        raise HTTPException(status_code=404, detail="Project not found")

    # Calculate metrics
    scl.qc.calculate_qc_metric(adata, sample_key="sampleID")

    # Return metrics for Plotly.js
    return {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "metrics": {
            "n_genes_by_counts": adata.obs["n_genes_by_counts"].tolist(),
            "total_counts": adata.obs["total_counts"].tolist(),
            "pct_counts_mt": adata.obs.get("pct_counts_mt", []).tolist(),
            "log1p_total_counts": adata.obs.get("log1p_total_counts", []).tolist(),
        },
        "statistics": {
            "n_genes_by_counts": {
                "mean": float(adata.obs["n_genes_by_counts"].mean()),
                "median": float(adata.obs["n_genes_by_counts"].median()),
            },
            "total_counts": {
                "mean": float(adata.obs["total_counts"].mean()),
                "median": float(adata.obs["total_counts"].median()),
            },
        },
    }


@router.post("/filter-preview")
async def preview_qc_filtering(
    project_id: str,
    min_genes: int = 200,
    max_mt_percent: float = 20.0,
    min_counts: int = 1000,
):
    """
    Preview QC filtering results for interactive Plotly.js charts.

    Args:
        project_id: Project identifier
        min_genes: Minimum genes threshold
        max_mt_percent: Maximum mitochondrial percentage
        min_counts: Minimum counts threshold

    Returns:
        Filtered data and statistics for Plotly.js visualization
    """
    import scLucid as scl
    from scLucid.web.services.data_manager import get_project_data

    adata = get_project_data(project_id)
    if adata is None:
        raise HTTPException(status_code=404, detail="Project not found")

    # Calculate metrics if not already done
    if "n_genes_by_counts" not in adata.obs:
        scl.qc.calculate_qc_metric(adata, sample_key="sampleID")

    # Apply filters
    mask = pd.Series(True, index=adata.obs_names)
    mask &= adata.obs["n_genes_by_counts"] >= min_genes
    mask &= adata.obs["pct_counts_mt"] <= max_mt_percent
    mask &= adata.obs["total_counts"] >= min_counts

    n_total = len(adata)
    n_retained = mask.sum()
    n_removed = n_total - n_retained

    # Prepare data for Plotly.js scatter plot
    scatter_data = {
        "x": adata.obs["n_genes_by_counts"].tolist(),
        "y": adata.obs["pct_counts_mt"].tolist(),
        "retained": mask.tolist(),
        "hovertext": [
            f"Cell: {i}<br>Genes: {g}<br>MT%: {mt:.1f}"
            for i, g, mt in zip(
                adata.obs_names,
                adata.obs["n_genes_by_counts"],
                adata.obs["pct_counts_mt"],
            )
        ],
    }

    return {
        "scatter_data": scatter_data,
        "summary": {
            "total": n_total,
            "retained": int(n_retained),
            "removed": int(n_removed),
            "retained_pct": float(n_retained / n_total * 100),
        },
        "thresholds": {
            "min_genes": min_genes,
            "max_mt_percent": max_mt_percent,
            "min_counts": min_counts,
        },
    }


@router.get("/violin-data/{project_id}")
async def get_violin_plot_data(project_id: str, metric: str = "n_genes_by_counts"):
    """
    Get data for Plotly.js violin plot.

    Args:
        project_id: Project identifier
        metric: Metric name to plot

    Returns:
        Violin plot data formatted for Plotly.js
    """
    from scLucid.web.services.data_manager import get_project_data

    adata = get_project_data(project_id)
    if adata is None:
        raise HTTPException(status_code=404, detail="Project not found")

    if metric not in adata.obs:
        raise HTTPException(status_code=400, detail=f"Metric {metric} not found")

    # Prepare data for Plotly.js violin plot
    metric_data = adata.obs[metric].dropna().tolist()

    return {
        "metric": metric,
        "data": metric_data,
        "stats": {
            "min": float(np.min(metric_data)),
            "q1": float(np.percentile(metric_data, 25)),
            "median": float(np.median(metric_data)),
            "q3": float(np.percentile(metric_data, 75)),
            "max": float(np.max(metric_data)),
        },
    }


@router.get("/histogram-data/{project_id}")
async def get_histogram_data(project_id: str, metric: str = "n_genes_by_counts", nbins: int = 50):
    """
    Get histogram data for Plotly.js.

    Args:
        project_id: Project identifier
        metric: Metric name to plot
        nbins: Number of bins

    Returns:
        Histogram data formatted for Plotly.js
    """
    from scLucid.web.services.data_manager import get_project_data

    adata = get_project_data(project_id)
    if adata is None:
        raise HTTPException(status_code=404, detail="Project not found")

    if metric not in adata.obs:
        raise HTTPException(status_code=400, detail=f"Metric {metric} not found")

    metric_data = adata.obs[metric].dropna().values

    # Calculate histogram
    hist, bin_edges = np.histogram(metric_data, bins=nbins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    return {
        "metric": metric,
        "x": bin_centers.tolist(),
        "y": hist.tolist(),
        "nbins": nbins,
    }


@router.post("/apply-filter")
async def apply_qc_filter(project_id: str, config: QCConfig):
    """
    Apply QC filtering and return filtered data.

    Args:
        project_id: Project identifier
        config: QC configuration

    Returns:
        Filtered AnnData data summary
    """
    import scLucid as scl
    from scLucid.web.services.data_manager import get_project_data, update_project_data

    adata = get_project_data(project_id)
    if adata is None:
        raise HTTPException(status_code=404, detail="Project not found")

    # Run QC
    adata_qc = scl.qc.run_standard_qc(adata, config=config)

    # Update project data
    update_project_data(project_id, adata_qc)

    return {
        "status": "success",
        "original_n_obs": adata.n_obs,
        "filtered_n_obs": adata_qc.n_obs,
        "cells_removed": adata.n_obs - adata_qc.n_obs,
        "retention_rate": float(adata_qc.n_obs / adata.n_obs * 100),
    }


@router.get("/available-metrics/{project_id}")
async def get_available_metrics(project_id: str):
    """
    Get list of available QC metrics for a project.

    Args:
        project_id: Project identifier

    Returns:
        List of available metric names
    """
    from scLucid.web.services.data_manager import get_project_data

    adata = get_project_data(project_id)
    if adata is None:
        raise HTTPException(status_code=404, detail="Project not found")

    qc_metrics = []
    for col in adata.obs.columns:
        if col.startswith(("log1p_", "pct_", "n_", "total_")) or col in [
            "phase",
            "S_score",
            "G2M_score",
        ]:
            qc_metrics.append(col)

    return {"metrics": qc_metrics}
