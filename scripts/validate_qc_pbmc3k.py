"""Validate QC on pbmc3k.h5ad (non-tumor baseline, 4 samples)."""

import json
from pathlib import Path

import scanpy as sc

from scLucid.qc.workflow import run_standard_qc
from scLucid.qc.config import QCWorkflowConfig, MetricsReportingConfig, MarkingConfig, DoubletConfig

DATA_PATH = Path("data/pbmc3k.h5ad")
OUT_DIR = Path("validation_outputs/qc_pbmc3k")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main():
    print(f"Loading {DATA_PATH} ...")
    adata = sc.read_h5ad(DATA_PATH)
    print(f"  shape: {adata.shape}")
    print(f"  samples: {adata.obs['sampleID'].nunique()}")

    config = QCWorkflowConfig(
        save_dir=str(OUT_DIR / "results"),
        use_recommendations=True,
        threshold_mode="hierarchical",
        metrics_reporting_config=MetricsReportingConfig(show_plots=False),
        marking_config=MarkingConfig(show_plots=False, plot_outliers=False),
        doublet_config=DoubletConfig(run_algorithm=False, use_heuristics=False, show_plots=False),
        filter_config={"criteria_to_filter": ["outlier_min_genes", "outlier_mt"]},
    )

    adata_filtered = run_standard_qc(adata, config=config)

    qc_uns = adata_filtered.uns.get("sclucid", {}).get("qc", {})
    context = qc_uns.get("context", {}).get("data", {})
    recommendation = qc_uns.get("recommendation", {}).get("data")
    sample_thresholds = qc_uns.get("sample_thresholds", {}).get("data", {})
    filtering_summary = qc_uns.get("filtering_summary", {}).get("data", {})
    warnings = qc_uns.get("warnings", {}).get("data", [])

    summary = {
        "dataset": "pbmc3k.h5ad",
        "shape_before": list(adata.shape),
        "shape_after": list(adata_filtered.shape),
        "context": context,
        "recommendation_strategy": recommendation.get("overall_strategy") if recommendation else None,
        "sample_thresholds": sample_thresholds,
        "filtering_summary": filtering_summary,
        "warnings": warnings,
    }

    json_path = OUT_DIR / "execution_summary.json"
    json_path.write_text(json.dumps(summary, indent=2, default=str))
    print(f"Wrote {json_path}")

    md_lines = [
        "# QC Validation: pbmc3k.h5ad",
        "",
        f"- **Dataset**: pbmc3k.h5ad",
        f"- **Cells before/after**: {adata.n_obs} / {adata_filtered.n_obs}",
        f"- **Threshold mode**: {context.get('threshold_mode')}",
        f"- **Strategy**: {summary['recommendation_strategy']}",
        "",
        "## Sample thresholds",
        "",
        "```json",
        json.dumps(sample_thresholds, indent=2, default=str),
        "```",
        "",
        "## Filtering summary",
        "",
        "```json",
        json.dumps(filtering_summary, indent=2, default=str),
        "```",
        "",
        "## Warnings",
        "",
    ]
    for w in warnings:
        md_lines.append(f"- {w}")
    md_lines.append("")

    md_path = OUT_DIR / "validation_notes.md"
    md_path.write_text("\n".join(md_lines))
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
