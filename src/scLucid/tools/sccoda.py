"""
Compositional data analysis using scCODA.

This module provides a wrapper for scCODA to identify significant changes in
cell type proportions between different experimental conditions.
"""

import logging
import os
from typing import List, Optional

import anndata
import matplotlib.pyplot as plt
import pandas as pd
import sccoda.util.comp_ana as sccoda_plot
from sccoda.model import compositional_data as com_data
from sccoda.model import scCODAModel

log = logging.getLogger(__name__)


def run_sccoda(
    adata: anndata.AnnData,
    cell_type_col: str,
    sample_col: str,
    condition_col: str,
    reference_level: Optional[str] = None,
    reference_cell_type: Optional[str] = None,
    min_cells: int = 10,
    pseudo_count: float = 0.5,
    n_samples: int = 25000,
    n_burnin: int = 5000,
    credible_interval: float = 0.95,
    cell_types_to_exclude: Optional[List[str]] = None,
    out_dir: Optional[str] = None,
    key_added: str = "sccoda",
    plot: bool = True,
    copy: bool = False,
) -> anndata.AnnData:
    """
    Run scCODA analysis to identify compositional changes in cell types.
    ... (docstring remains the same) ...
    """
    if copy:
        adata = adata.copy()

    # --- 1. Parameter and Data Validation ---
    log.info(f"Starting scCODA analysis on {adata.n_obs} cells")
    for col in [cell_type_col, sample_col, condition_col]:
        if col not in adata.obs:
            raise ValueError(f"Column '{col}' not found in adata.obs")

    # --- 2. Data Preparation (Efficiently) ---
    log.info("Preparing data for scCODA analysis")
    obs_df = adata.obs[[sample_col, condition_col, cell_type_col]].copy()
    if cell_types_to_exclude:
        obs_df = obs_df[~obs_df[cell_type_col].isin(cell_types_to_exclude)]
        log.info(f"Excluded {len(cell_types_to_exclude)} cell types")

    count_df = pd.crosstab(obs_df[sample_col], obs_df[cell_type_col])
    if pseudo_count > 0:
        count_df[count_df < min_cells] += pseudo_count

    # Create covariate DataFrame
    covariate_df = (
        obs_df[[sample_col, condition_col]].drop_duplicates().set_index(sample_col)
    )

    # Align data
    count_df, covariate_df = count_df.align(covariate_df, axis=0, join="inner")

    # --- 3. Set up scCODA Model ---
    if reference_cell_type is None:
        reference_cell_type = count_df.columns[-1]
        log.info(f"Using last cell type as reference: {reference_cell_type}")

    if reference_level is None:
        reference_level = covariate_df[condition_col].unique()[0]
        log.info(f"Using first condition as reference level: {reference_level}")

    data_coda = com_data.from_pandas(
        count_df, covariate_df, covariate_columns=[condition_col]
    )

    # --- 4. Run scCODA Model ---
    log.info(
        f"Running scCODA model with {n_samples} samples and {n_burnin} burn-in iterations"
    )
    model = scCODAModel.from_formula(
        f"~{condition_col}", data_coda, reference_cell_type=reference_cell_type
    )
    result = model.sample_hmc(n_draws=n_samples, num_burnin=n_burnin)

    # --- 5. Process and Store Results ---
    log.info("Processing scCODA results")
    final_effects = result.credible_effects()

    adata.uns[key_added] = {
        "summary": result.summary_df,
        "final_effects": final_effects,
        "credible_interval": credible_interval,
        "reference_level": reference_level,
        "reference_cell_type": reference_cell_type,
    }

    significant_effects = final_effects[final_effects["Final"]].shape[0]
    log.info(f"Found {significant_effects} significant cell type changes")

    if significant_effects > 0:
        sig_results = final_effects[final_effects["Final"]].copy()
        sig_results["effect"] = sig_results.apply(
            lambda row: f"{row['Covariate']}[{row['Contrast']}] "
            f"(logFC: {row['log-fold change']:.2f})",
            axis=1,
        )
        effect_map = sig_results.set_index("Cell Type")["effect"].to_dict()

        adata.obs[f"{key_added}_effect"] = (
            adata.obs[cell_type_col].map(effect_map).fillna("Not significant")
        )
        log.info(f"Added compositional effects to adata.obs['{key_added}_effect']")

    # --- 6. Plotting ---
    if plot:
        log.info("Generating visualizations")
        if out_dir is None:
            out_dir = os.path.join(os.getcwd(), "sccoda_results")
        os.makedirs(out_dir, exist_ok=True)

        try:
            sccoda_plot.plot_credible_intervals(result).savefig(
                os.path.join(out_dir, "credible_intervals.png")
            )
            plt.close()
            sccoda_plot.plot_effects(result, ["log-fold change"]).savefig(
                os.path.join(out_dir, "effects_plot.png")
            )
            plt.close()
            log.info(f"Saved visualizations to {out_dir}")
        except Exception as e:
            log.warning(f"Error during plotting: {e}")

    return adata
