"""
QC Step-by-Step -- Simple API Layer Example

Demonstrates composable QC functions for analysts who want to inspect
or replace individual steps. Each step can be run independently, inspected,
and reconfigured before proceeding to the next.

Use this when:
- You need to review QC thresholds before filtering
- You want to try multiple doublet detection methods
- You need sample-aware adaptive thresholds for multi-sample data
- You want to export QC reports for reviewers

The default stance is conservative: ambient RNA correction, CellBender/SoupX,
ScDblFinder, and R bridges are not run here. Treat them as project-specific
upstream or optional analyses when the dataset and environment justify them.
"""

from __future__ import annotations

from pathlib import Path

import scanpy as sc

import scLucid as scl

DATA_PATH = Path("data/pbmc3k.h5ad")


def main() -> None:
    # Load raw data
    adata = sc.read_h5ad(DATA_PATH)
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    else:
        adata.layers["counts"] = adata.X.copy()

    if "sampleID" not in adata.obs.columns:
        adata.obs["sampleID"] = "pbmc3k"

    print(f"Input: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

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
    print("Metrics added: n_genes_by_counts, total_counts, pct_counts_mt, phase, ...")

    # ---------------------------------------------------------------------------
    # Step 2: Read-only policy review
    # ---------------------------------------------------------------------------
    print("\n--- Step 2: Evidence-Calibrated Policy Review ---")
    qc_review = scl.recommend_qc_policy(
        adata,
        scl.ProjectContext(
            dataset_type="pbmc_or_blood",
            species="human",
            sample_key="sampleID",
            input_provenance="filtered_counts",
        ),
    )
    print(f"Status: {qc_review.status}")
    print(f"Reason: {qc_review.reason}")
    print(f"Candidate impacts: {qc_review.comparison}")
    print(f"Next action: {qc_review.next_action}")

    # ---------------------------------------------------------------------------
    # Step 3: Ambient RNA / empty-droplet diagnostics
    # ---------------------------------------------------------------------------
    print("\n--- Step 3: Ambient RNA Diagnostics ---")
    ambient = scl.qc.diagnose_ambient_rna(adata, layer="counts", top_n_genes=10)
    print(
        "Ambient RNA:",
        f"available={ambient.get('available')}",
        f"diagnostic_only={ambient.get('diagnostic_only')}",
        f"risk={ambient.get('risk_level')}",
    )
    empty = scl.qc.ambient.diagnose_empty_droplets(
        adata,
        layer="counts",
        min_barcodes=min(100, max(20, adata.n_obs // 2)),
        top_n_genes=10,
    )
    print(
        "Empty droplets:",
        f"available={empty.get('available')}",
        f"diagnostic_only={empty.get('diagnostic_only')}",
        f"risk={empty.get('risk_level')}",
    )
    print("These diagnostics do not apply CellBender/SoupX/scAR correction.")

    # ---------------------------------------------------------------------------
    # Step 4: Doublet detection (ensemble)
    # ---------------------------------------------------------------------------
    print("\n--- Step 4: Doublet Detection ---")
    doublet_rates = scl.qc.generate_doublet_rates(adata, sample_key="sampleID")
    print(f"Expected doublet rates: {doublet_rates}")

    adata = scl.qc.predict_doublets(
        adata,
        config=scl.qc.DoubletConfig(
            method="scrublet",
            expected_doublet_rate=doublet_rates,
            use_heuristics=True,
            merge_strategy="weighted_average",
            # Note: ignore_coexpression_pairs defaults to empty.
            # For tumor/EMT studies, explicitly set:
            #   ignore_coexpression_pairs=(("Epithelial", "Mesenchymal"),)
        ),
        sample_key="sampleID",
    )
    n_doublets = int(adata.obs["predicted_doublet"].sum())
    print(f"Predicted doublets: {n_doublets} ({n_doublets/adata.n_obs*100:.1f}%)")
    for col in ["algorithm_doublet_score", "heterotypic_doublet_risk", "homotypic_doublet_risk"]:
        if col in adata.obs:
            print(f"{col}: mean={adata.obs[col].mean():.3f}")

    # ---------------------------------------------------------------------------
    # Step 5: Resolve threshold decisions without filtering
    # ---------------------------------------------------------------------------
    print("\n--- Step 5: Threshold Decisions ---")
    threshold_decision = scl.qc.decide_qc_thresholds(
        adata,
        threshold_method="mad",
        threshold_policy="mad_then_intelligent",
    )
    print("Resolved thresholds:")
    print(threshold_decision["resolved_thresholds"].to_dict())

    # ---------------------------------------------------------------------------
    # Step 6: Apply threshold decisions as evidence labels
    # ---------------------------------------------------------------------------
    print("\n--- Step 6: Threshold Evidence ---")
    threshold_application = scl.qc.apply_qc_threshold_decision(
        adata,
        resolved_thresholds=threshold_decision["resolved_thresholds"],
        sample_key="sampleID",
        filter_cells_result=False,
    )
    adata = threshold_application["adata"]
    print(f"Applied threshold evidence: {threshold_application['decision_record']}")

    # ---------------------------------------------------------------------------
    # Step 7: Unified QC decisions
    # ---------------------------------------------------------------------------
    print("\n--- Step 7: Unified Decisions ---")
    decision_summary = scl.qc.build_qc_decisions(
        adata,
        tissue_type="pbmc_or_blood",
        policy="conservative",
        sample_key="sampleID",
    )
    adata.obs["qc_remove"] = adata.obs["qc_decision"].astype(str) == "remove"
    print(f"QC decisions: {decision_summary['decision_counts']}")

    # ---------------------------------------------------------------------------
    # Step 8: Filter cells
    # ---------------------------------------------------------------------------
    print("\n--- Step 8: Filtering ---")
    adata_filtered = scl.qc.filter_cells(
        adata,
        config=scl.qc.FilterConfig(
            criteria_to_filter=[
                "qc_remove",
            ],
            combination_logic="any",
        ),
        copy=True,
    )
    print(f"Before: {adata.n_obs:,} cells -> After: {adata_filtered.n_obs:,} cells")
    print(f"Retention: {adata_filtered.n_obs/adata.n_obs*100:.1f}%")
    print("Review doublet-heavy clusters before treating all doublet calls as hard failures.")

    # ---------------------------------------------------------------------------
    # Step 9: Export QC report
    # ---------------------------------------------------------------------------
    print("\n--- Step 9: Report ---")
    scl.qc.generate_qc_report(
        adata_filtered,
        save_dir="results/qc_report",
        sample_key="sampleID",
        adata_before=adata,
    )
    print("QC report saved to results/qc_report/")

    # Save for downstream preprocessing example
    output_path = Path("results/qc_filtered.h5ad")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata_filtered.write_h5ad(output_path)
    print("\nQC complete!")
    print(f"Final: {adata_filtered.n_obs:,} cells x {adata_filtered.n_vars:,} genes")
    print(f"Saved to: {output_path}")


if __name__ == "__main__":
    main()
