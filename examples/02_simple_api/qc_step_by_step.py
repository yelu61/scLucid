"""
QC Step-by-Step — Simple API Layer Example

Demonstrates composable QC functions for analysts who want to inspect
or replace individual steps. Each step can be run independently, inspected,
and reconfigured before proceeding to the next.

Use this when:
- You need to review QC thresholds before filtering
- You want to try multiple doublet detection methods
- You need sample-aware adaptive thresholds for multi-sample data
- You want to export QC reports for reviewers
"""

from pathlib import Path

import scanpy as sc

import scLucid as scl


DATA_PATH = Path("data/pbmc3k.h5ad")

# Load raw data
adata = sc.read_h5ad(DATA_PATH)
if "counts" in adata.layers:
    adata.X = adata.layers["counts"].copy()
else:
    adata.layers["counts"] = adata.X.copy()

if "sampleID" not in adata.obs.columns:
    adata.obs["sampleID"] = "pbmc3k"

print(f"Input: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# ---------------------------------------------------------------------------
# Step 1: Calculate QC metrics
# ---------------------------------------------------------------------------
print("\n--- Step 1: QC Metrics ---")
adata = scl.qc.calculate_qc_metric(
    adata,
    sample_key="sampleID",
    calculate_cell_cycle=True,
    cell_cycle_species="human",
)
print(f"Metrics added: n_genes_by_counts, total_counts, pct_counts_mt, phase, ...")

# ---------------------------------------------------------------------------
# Step 2: Intelligent threshold recommendations
# ---------------------------------------------------------------------------
print("\n--- Step 2: Intelligent Recommendations ---")
rec = scl.qc.recommend_intelligent_qc(adata, tissue_type="pbmc_or_blood")
print(f"Recommended min_genes: {rec.min_genes.threshold}")
print(f"Recommended max_mt:    {rec.max_mt_percent.threshold}")
print(f"Overall confidence:    {rec.overall_confidence:.2f}")

# ---------------------------------------------------------------------------
# Step 3: Doublet detection (ensemble)
# ---------------------------------------------------------------------------
print("\n--- Step 3: Doublet Detection ---")
doublet_rates = scl.qc.generate_doublet_rates(adata, sample_key="sampleID")
print(f"Expected doublet rates: {doublet_rates}")

adata = scl.qc.predict_doublets(
    adata,
    config=scl.qc.DoubletConfig(
        method="scrublet",
        expected_doublet_rate=doublet_rates,
        use_heuristics=True,
        merge_strategy="weighted_average",
    ),
    sample_key="sampleID",
)
n_doublets = int(adata.obs["predicted_doublet"].sum())
print(f"Predicted doublets: {n_doublets} ({n_doublets/adata.n_obs*100:.1f}%)")

# ---------------------------------------------------------------------------
# Step 4: Suggest global thresholds (MAD-based)
# ---------------------------------------------------------------------------
print("\n--- Step 4: Global Thresholds ---")
threshold_table, suggested = scl.qc.suggest_qc_thresholds(adata, method="mad")
print("Suggested thresholds:")
print(threshold_table.to_string())

# ---------------------------------------------------------------------------
# Step 5: Adaptive per-sample marking
# ---------------------------------------------------------------------------
print("\n--- Step 5: Adaptive Marking ---")
adata = scl.qc.mark_low_quality_cells_adaptive(
    adata,
    batch_key="sampleID",
    metrics=["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    method="hierarchical",
)
n_adaptive = int(adata.obs["outlier_n_genes_by_counts_adaptive"].sum())
print(f"Adaptive outliers: {n_adaptive} cells")

# ---------------------------------------------------------------------------
# Step 6: Unified marking (combine all flags)
# ---------------------------------------------------------------------------
print("\n--- Step 6: Unified Marking ---")
adata = scl.qc.mark_low_quality_cell(
    adata,
    config=scl.qc.MarkingConfig(
        thresholds=scl.qc.QCThresholds(
            min_genes=suggested.min_genes,
            pc_mt=suggested.pc_mt,
        ),
    ),
)

# ---------------------------------------------------------------------------
# Step 7: Filter cells
# ---------------------------------------------------------------------------
print("\n--- Step 7: Filtering ---")
adata_filtered = scl.qc.filter_cells(
    adata,
    config=scl.qc.FilterConfig(
        criteria_to_filter=[
            "outlier_min_genes",
            "outlier_mt",
            "predicted_doublet",
        ],
        combination_logic="any",
    ),
    copy=True,
)
print(f"Before: {adata.n_obs:,} cells → After: {adata_filtered.n_obs:,} cells")
print(f"Retention: {adata_filtered.n_obs/adata.n_obs*100:.1f}%")

# ---------------------------------------------------------------------------
# Step 8: Export QC report
# ---------------------------------------------------------------------------
print("\n--- Step 8: Report ---")
scl.qc.generate_qc_report(
    adata_filtered,
    save_dir="results/qc_report",
    sample_key="sampleID",
    adata_before=adata,
)
print("QC report saved to results/qc_report/")

print("\n✅ QC complete!")
print(f"Final: {adata_filtered.n_obs:,} cells × {adata_filtered.n_vars:,} genes")
