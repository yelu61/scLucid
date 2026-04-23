"""
Compositional data analysis using scCODA.

This module provides a wrapper for scCODA to identify significant changes in
cell type proportions between different experimental conditions. It also provides
integrated visualization, including overlaying significance from scCODA on cell type proportion plots.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata
import matplotlib.pyplot as plt
import pandas as pd
import sccoda.util.comp_ana as sccoda_plot
import seaborn as sns
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
    out_dir: Optional[Path] = None,
    key_added: str = "sccoda",
    plot: bool = True,
    copy: bool = False,
) -> anndata.AnnData:
    """
    Run scCODA analysis to identify compositional changes in cell types.
    Results are stored in adata.uns['sclucid']['sccoda'][key_added].

    Args:
        adata: AnnData object.
        cell_type_col: Column in .obs for cell type annotation.
        sample_col: Column in .obs for sample ID.
        condition_col: Column in .obs for experimental condition.
        reference_level: Reference group for the condition variable.
        reference_cell_type: Reference cell type to anchor the model.
        min_cells: Minimum number of cells to avoid low-count errors.
        pseudo_count: Value to add to low-counts for stability.
        n_samples: Number of HMC samples for the model.
        n_burnin: Number of burn-in steps for HMC.
        credible_interval: Credible interval for credible effects.
        cell_types_to_exclude: List of cell types to exclude from the analysis.
        out_dir: Output directory for results (figures, tables).
        key_added: Key under adata.uns['sclucid']['sccoda'] for storing results.
        plot: Whether to generate and save plots.
        copy: Whether to operate on a copy of adata.

    Returns:
        AnnData object with results stored in .uns.
    """
    if copy:
        adata = adata.copy()

    adata.uns.setdefault("sclucid", {}).setdefault("sccoda", {})

    # 1. Data validation
    log.info(f"Starting scCODA analysis on {adata.n_obs} cells")
    for col in [cell_type_col, sample_col, condition_col]:
        if col not in adata.obs:
            raise ValueError(f"Column '{col}' not found in adata.obs")

    # 2. Data preparation
    log.info("Preparing data for scCODA analysis")
    obs_df = adata.obs[[sample_col, condition_col, cell_type_col]].copy()
    if cell_types_to_exclude:
        obs_df = obs_df[~obs_df[cell_type_col].isin(cell_types_to_exclude)]
        log.info(f"Excluded {len(cell_types_to_exclude)} cell types: {cell_types_to_exclude}")

    count_df = pd.crosstab(obs_df[sample_col], obs_df[cell_type_col])
    if pseudo_count > 0:
        count_df[count_df < min_cells] += pseudo_count

    # Create covariate DataFrame
    covariate_df = obs_df[[sample_col, condition_col]].drop_duplicates().set_index(sample_col)

    # Align data
    count_df, covariate_df = count_df.align(covariate_df, axis=0, join="inner")

    # 3. scCODA model setup
    if reference_cell_type is None:
        reference_cell_type = count_df.columns[-1]
        log.info(f"Using last cell type as reference: {reference_cell_type}")

    if reference_level is None:
        reference_level = covariate_df[condition_col].unique()[0]
        log.info(f"Using first condition as reference level: {reference_level}")

    data_coda = com_data.from_pandas(count_df, covariate_df, covariate_columns=[condition_col])

    # 4. Run the scCODA model
    log.info(f"Running scCODA model with {n_samples} samples and {n_burnin} burn-in iterations")
    model = scCODAModel.from_formula(
        f"~{condition_col}", data_coda, reference_cell_type=reference_cell_type
    )
    result = model.sample_hmc(n_draws=n_samples, num_burnin=n_burnin)

    # 5. Process and store results
    log.info("Processing scCODA results")
    final_effects = result.credible_effects()

    # Structured results
    adata.uns["sclucid"]["sccoda"][key_added] = {
        "summary": result.summary_df,
        "final_effects": final_effects,
        "credible_interval": credible_interval,
        "reference_level": reference_level,
        "reference_cell_type": reference_cell_type,
        "params": {
            "cell_type_col": cell_type_col,
            "sample_col": sample_col,
            "condition_col": condition_col,
            "min_cells": min_cells,
            "pseudo_count": pseudo_count,
            "n_samples": n_samples,
            "n_burnin": n_burnin,
            "cell_types_to_exclude": cell_types_to_exclude,
        },
    }

    significant_effects = final_effects[final_effects["Final"]].shape[0]
    log.info(f"Found {significant_effects} significant cell type changes")

    if significant_effects > 0:
        sig_results = final_effects[final_effects["Final"]].copy()
        sig_results["effect"] = sig_results.apply(
            lambda row: f"{row['Covariate']}[{row['Contrast']}] (logFC: {row['log-fold change']:.2f})",
            axis=1,
        )
        effect_map = sig_results.set_index("Cell Type")["effect"].to_dict()

        adata.obs[f"{key_added}_effect"] = (
            adata.obs[cell_type_col].map(effect_map).fillna("Not significant")
        )
        log.info(f"Added compositional effects to adata.obs['{key_added}_effect']")

    # 6. Plotting and export
    if plot:
        log.info("Generating visualizations")
        if out_dir is None:
            out_dir = Path.cwd() / "sccoda_results"
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        try:
            sccoda_plot.plot_credible_intervals(result).savefig(
                out_dir / f"{key_added}_credible_intervals.png"
            )
            plt.close()
            sccoda_plot.plot_effects(result, ["log-fold change"]).savefig(
                out_dir / f"{key_added}_effects_plot.png"
            )
            plt.close()
            # Export main tables
            result.summary_df.to_csv(out_dir / f"{key_added}_summary.csv", index=False)
            final_effects.to_csv(out_dir / f"{key_added}_final_effects.csv", index=False)
            log.info(f"Saved visualizations and tables to {out_dir}")
        except Exception as e:
            log.warning(f"Error during plotting or export: {e}")

    return adata


def run_sccoda_batch(
    adatas: List[anndata.AnnData],
    cell_type_col: str,
    sample_col: str,
    condition_col: str,
    out_dir: Path,
    sample_ids: Optional[List[str]] = None,
    **kwargs,
) -> Dict[str, Optional[anndata.AnnData]]:
    """
    Batch run scCODA for a list of AnnData objects.
    Returns dict[sample_id] = AnnData with scCODA results or None if failed.
    """
    results = {}
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if sample_ids is None:
        sample_ids = [f"sample{i + 1}" for i in range(len(adatas))]
    for adata, sid in zip(adatas, sample_ids):
        sample_dir = out_dir / sid
        try:
            results[sid] = run_sccoda(
                adata=adata,
                cell_type_col=cell_type_col,
                sample_col=sample_col,
                condition_col=condition_col,
                out_dir=sample_dir,
                **kwargs,
            )
            log.info(f"scCODA completed for {sid}")
        except Exception as e:
            log.error(f"scCODA failed for {sid}: {e}")
            results[sid] = None
    return results


def summarize_sccoda(
    adata: anndata.AnnData,
    key_added: str = "sccoda",
    save_dir: Optional[Path] = None,
    top_n: int = 10,
):
    """
    Summarize and export main scCODA result tables and top effects.
    """
    sccoda = adata.uns.get("sclucid", {}).get("sccoda", {}).get(key_added, {})
    summary = sccoda.get("summary")
    effects = sccoda.get("final_effects")
    params = sccoda.get("params", {})

    if save_dir:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
    # Export tables
    if summary is not None and save_dir:
        summary.to_csv(save_dir / f"{key_added}_summary.csv", index=False)
    if effects is not None and save_dir:
        effects.to_csv(save_dir / f"{key_added}_final_effects.csv", index=False)
        sig = effects[effects["Final"]].copy()
        sig.head(top_n).to_csv(save_dir / f"{key_added}_top_effects.csv", index=False)
    # Export params
    if params and save_dir:
        with open(save_dir / f"{key_added}_params.txt", "w") as f:
            for k, v in params.items():
                f.write(f"{k}: {v}\n")
    log.info(f"scCODA summary exported to: {save_dir or Path.cwd()}")


def plot_sccoda_proportion_with_significance(
    adata: anndata.AnnData,
    celltype_col: str,
    sample_col: str,
    condition_col: str,
    sccoda_key: str = "sccoda",
    out_path: Optional[Path] = None,
    figsize: Tuple[int, int] = (10, 6),
    star_threshold: float = 0.05,
):
    """
    Visualize cell type proportions per sample/condition, overlaying scCODA statistical results.
    Significant cell types (by scCODA) are marked with a star.

    Args:
        adata: AnnData object with .obs and .uns['sclucid']['sccoda'][sccoda_key] results.
        celltype_col: Column in .obs for cell type.
        sample_col: Column in .obs for sample.
        condition_col: Column in .obs for condition/group.
        sccoda_key: Key in .uns['sclucid']['sccoda'] for scCODA run.
        out_path: If provided, saves the plot to this path.
        figsize: Figure size.
        star_threshold: p-value threshold for significance.
    """
    # 1. Compute cell type proportions
    count_df = adata.obs.groupby([sample_col, celltype_col]).size().unstack(fill_value=0)
    prop_df = count_df.div(count_df.sum(axis=1), axis=0)
    cond_map = adata.obs.drop_duplicates(sample_col)[[sample_col, condition_col]].set_index(
        sample_col
    )[condition_col]
    df_long = prop_df.reset_index().melt(
        id_vars=sample_col, var_name="celltype", value_name="proportion"
    )
    df_long["condition"] = df_long[sample_col].map(cond_map)

    # 2. Identify significant cell types from scCODA
    sccoda = adata.uns.get("sclucid", {}).get("sccoda", {}).get(sccoda_key, {})
    effects = sccoda.get("final_effects", pd.DataFrame())
    if "p-value" in effects.columns:
        significant_cts = effects.loc[
            (effects["Final"]) & (effects["p-value"] < star_threshold), "Cell Type"
        ].tolist()
    else:
        significant_cts = (
            effects.loc[effects["Final"], "Cell Type"].tolist() if "Final" in effects else []
        )

    # 3. Plot
    plt.figure(figsize=figsize)
    ax = sns.boxplot(data=df_long, x="celltype", y="proportion", hue="condition")
    plt.xticks(rotation=45)
    ymax = df_long["proportion"].max()
    for i, ct in enumerate(ax.get_xticklabels()):
        label = ct.get_text()
        if label in significant_cts:
            ax.annotate(
                "*",
                xy=(i, ymax * 1.05),
                ha="center",
                va="bottom",
                color="red",
                fontsize=18,
            )
    plt.tight_layout()
    if out_path:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_path)
        plt.close()
    else:
        plt.show()
        plt.close()
