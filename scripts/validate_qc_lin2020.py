"""Validate QC on lin2020.pdac.h5ad (multi-sample tumor) comparing pooled vs hierarchical."""

import json
from pathlib import Path

import scanpy as sc

from scLucid.qc.workflow import run_standard_qc
from scLucid.qc.config import QCWorkflowConfig, MetricsReportingConfig, MarkingConfig, DoubletConfig

DATA_PATH = Path("data/lin2020.pdac.h5ad")
OUT_DIR = Path("validation_outputs/qc_lin2020")
OUT_DIR.mkdir(parents=True, exist_ok=True)

SUSPECT_SAMPLES = ["GSM4679533", "GSM4679535"]


def run_mode(adata, mode: str):
    print(f"\n=== Running {mode} mode ===")
    config = QCWorkflowConfig(
        save_dir=str(OUT_DIR / f"results_{mode}"),
        tissue_type="tumor",
        use_recommendations=True,
        threshold_mode=mode,
        metrics_reporting_config=MetricsReportingConfig(show_plots=False),
        marking_config=MarkingConfig(show_plots=False, plot_outliers=False),
        doublet_config=DoubletConfig(run_algorithm=False, use_heuristics=False, show_plots=False),
        filter_config={"criteria_to_filter": ["outlier_min_genes", "outlier_mt"]},
    )
    adata_f = run_standard_qc(adata, config=config)
    qc_uns = adata_f.uns.get("sclucid", {}).get("qc", {})
    return {
        "shape_after": list(adata_f.shape),
        "context": qc_uns.get("context", {}).get("data", {}),
        "recommendation": qc_uns.get("recommendation", {}).get("data"),
        "sample_thresholds": qc_uns.get("sample_thresholds", {}).get("data", {}),
        "filtering_summary": qc_uns.get("filtering_summary", {}).get("data", {}),
        "warnings": qc_uns.get("warnings", {}).get("data", []),
    }


def main():
    print(f"Loading {DATA_PATH} ...")
    adata = sc.read_h5ad(DATA_PATH)
    print(f"  shape: {adata.shape}")
    print(f"  samples: {adata.obs['sampleID'].nunique()}")
    for s in adata.obs["sampleID"].unique():
        print(f"    - {s}: {(adata.obs['sampleID'] == s).sum()} cells")

    results_pooled = run_mode(adata, "pooled")
    results_hierarchical = run_mode(adata, "hierarchical")

    # Compare suspect samples
    comparison = {"suspect_samples": {}}
    pooled_st = results_pooled["sample_thresholds"]
    hier_st = results_hierarchical["sample_thresholds"]
    for sample in SUSPECT_SAMPLES:
        comparison["suspect_samples"][sample] = {
            "pooled": pooled_st.get(sample, {}),
            "hierarchical": hier_st.get(sample, {}),
        }

    summary = {
        "dataset": "lin2020.pdac.h5ad",
        "shape_before": list(adata.shape),
        "pooled": results_pooled,
        "hierarchical": results_hierarchical,
        "comparison": comparison,
    }

    json_path = OUT_DIR / "execution_summary.json"
    json_path.write_text(json.dumps(summary, indent=2, default=str))
    print(f"\nWrote {json_path}")

    md_lines = [
        "# QC Validation: lin2020.pdac.h5ad",
        "",
        f"- **Dataset**: lin2020.pdac.h5ad",
        f"- **Cells before**: {adata.n_obs}",
        f"- **Cells after (pooled)**: {results_pooled['shape_after'][0]}",
        f"- **Cells after (hierarchical)**: {results_hierarchical['shape_after'][0]}",
        "",
        "## Suspect sample comparison",
        "",
    ]
    for sample in SUSPECT_SAMPLES:
        md_lines.append(f"### {sample}")
        md_lines.append("```json")
        md_lines.append(json.dumps(comparison["suspect_samples"][sample], indent=2, default=str))
        md_lines.append("```")
        md_lines.append("")

    md_lines.extend([
        "## Pooled warnings",
        "",
    ])
    for w in results_pooled["warnings"]:
        md_lines.append(f"- {w}")
    md_lines.append("")

    md_lines.extend([
        "## Hierarchical warnings",
        "",
    ])
    for w in results_hierarchical["warnings"]:
        md_lines.append(f"- {w}")
    md_lines.append("")

    md_path = OUT_DIR / "validation_notes.md"
    md_path.write_text("\n".join(md_lines))
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
