"""Tumor-focused bulk RNA-seq utilities.

These functions interpret bulk RNA-seq through the lens of tumor scRNA-seq
analysis: TME composition, immune infiltration, therapy response, and tumor
purity.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData

from .abundance import run_bulk_abundance_test
from .clinical import correlate_abundance_with_clinical
from .config import BulkAbundanceConfig, BulkClinicalAssociationConfig
from .deconvolution import deconvolve_bulk


_COMPARTMENTS = ["malignant", "immune", "stromal", "other"]


def _get_compartment_map(compartment_map: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    """Return a cell-type-label -> compartment map.

    If ``compartment_map`` is None, import the default map from the tumor
    microenvironment profiler so labels stay consistent across scLucid.
    """
    if compartment_map is not None:
        return {str(k).lower().replace("-", " ").replace("_", " "): v for k, v in compartment_map.items()}

    from ...tumor.microenvironment.deconvolution import TMEProfiler

    return dict(TMEProfiler._DEFAULT_COMPARTMENT_MAP)


def _normalize_label(label: str) -> str:
    return str(label).lower().replace("-", " ").replace("_", " ").strip()


def deconvolve_tumor_tme(
    adata_ref: AnnData,
    bulk_data: pd.DataFrame,
    cell_type_key: str,
    compartment_map: Optional[Dict[str, str]] = None,
    *,
    method: str = "DWLS",
    key_added: str = "tumor_tme",
    sample_key: str = "sampleID",
    **method_kwargs,
) -> AnnData:
    """Deconvolve bulk RNA-seq into tumor-TME compartments.

    This wraps ``deconvolve_bulk`` and aggregates inferred cell-type
    proportions into canonical compartments: ``malignant``, ``immune``,
    ``stromal``, and ``other``. It is designed for tumor scRNA-seq references
    where cell types have been annotated.

    Parameters
    ----------
    adata_ref
        Single-cell reference with cell type labels.
    bulk_data
        Bulk RNA-seq matrix (genes x samples).
    cell_type_key
        Column in ``adata_ref.obs`` with cell type labels.
    compartment_map
        Mapping from cell-type labels to compartments. If None, uses the
        ``TMEProfiler`` default compartment map.
    method
        Deconvolution backend ("DWLS" or "BayesPrism").
    key_added
        Key under ``adata_ref.uns["sclucid"]["tools"][key_added]``.
    sample_key
        Sample ID column for BayesPrism.
    **method_kwargs
        Passed to ``deconvolve_bulk``.

    Returns
    -------
    AnnData
        ``adata_ref`` with compartment proportions and inference tags stored.
    """
    adata_ref = deconvolve_bulk(
        adata_ref,
        bulk_data,
        cell_type_key=cell_type_key,
        method=method,
        key_added=key_added,
        sample_key=sample_key,
        **method_kwargs,
    )

    proportions = adata_ref.uns["sclucid"]["tools"][key_added]["proportions"]
    cmap = _get_compartment_map(compartment_map)

    compartment_rows = []
    unmapped: List[str] = []
    for sample_id, row in proportions.iterrows():
        bucket = {c: 0.0 for c in _COMPARTMENTS}
        for cell_type, value in row.items():
            norm = _normalize_label(cell_type)
            compartment = cmap.get(norm)
            if compartment is None:
                unmapped.append(str(cell_type))
                compartment = "other"
            bucket[compartment] += float(value)
        compartment_rows.append(pd.Series(bucket, name=sample_id))

    compartment_df = pd.DataFrame(compartment_rows)
    purity = compartment_df.get("malignant", 0.0)
    immune_score = compartment_df.get("immune", 0.0)
    stromal_score = compartment_df.get("stromal", 0.0)

    n_samples = bulk_data.shape[1]
    inference_level = "sample_level" if n_samples >= 3 else "exploratory_sample_level"
    valid_for_publication = n_samples >= 3

    adata_ref.uns["sclucid"]["tools"][key_added].update(
        {
            "compartment_proportions": compartment_df,
            "compartment_map": cmap,
            "unmapped_cell_types": sorted(set(unmapped)),
            "tumor_purity": purity,
            "immune_score": immune_score,
            "stromal_score": stromal_score,
            "inference_level": inference_level,
            "valid_for_publication_inference": valid_for_publication,
            "result_warning": (
                None
                if valid_for_publication
                else "Few bulk samples; TME proportions are exploratory."
            ),
        }
    )

    return adata_ref


def estimate_tumor_purity_from_bulk(
    bulk_data: pd.DataFrame,
    *,
    method: str = "expression",
    adata_ref: Optional[AnnData] = None,
    cell_type_key: Optional[str] = None,
    compartment_map: Optional[Dict[str, str]] = None,
    key_added: str = "tumor_purity",
) -> pd.DataFrame:
    """Estimate tumor purity per bulk sample.

    Parameters
    ----------
    bulk_data
        Bulk RNA-seq matrix (genes x samples).
    method
        - ``"expression"``: clean-room purity score based on a small set of
          epithelial/tumor marker genes and immune/stromal reference genes.
        - ``"tme"``: requires ``adata_ref`` and ``cell_type_key``; uses
          ``deconvolve_tumor_tme`` and returns the malignant compartment.
    adata_ref
        Single-cell reference for method="tme".
    cell_type_key
        Cell type column for method="tme".
    compartment_map
        Optional compartment map for method="tme".
    key_added
        Column name prefix for returned DataFrame.

    Returns
    -------
    pd.DataFrame
        One row per sample with purity estimate and inference tags.
    """
    if method == "tme":
        if adata_ref is None or cell_type_key is None:
            raise ValueError("method='tme' requires adata_ref and cell_type_key")
        adata_ref = deconvolve_tumor_tme(
            adata_ref,
            bulk_data,
            cell_type_key=cell_type_key,
            compartment_map=compartment_map,
            key_added=key_added,
        )
        purity = adata_ref.uns["sclucid"]["tools"][key_added]["tumor_purity"]
        result = pd.DataFrame({f"{key_added}": purity})
    elif method == "expression":
        # Clean-room expression-based estimate: ratio of tumor-like signal to
        # total tumor+immune+stromal signal. Uses a tiny built-in marker set
        # as a fallback when no sc reference is available.
        tumor_markers = ["EPCAM", "KRT5", "KRT8", "KRT18", "KRT19"]
        immune_markers = ["PTPRC", "CD3E", "CD68", "CD79A"]
        stromal_markers = ["COL1A1", "COL1A2", "ACTA2", "VIM"]

        def _mean_marker_score(df: pd.DataFrame, markers: List[str]) -> pd.Series:
            present = [g for g in markers if g in df.index]
            if not present:
                return pd.Series(0.0, index=df.columns)
            return df.loc[present].mean(axis=0)

        tumor_score = _mean_marker_score(bulk_data, tumor_markers)
        immune_score = _mean_marker_score(bulk_data, immune_markers)
        stromal_score = _mean_marker_score(bulk_data, stromal_markers)
        denom = tumor_score + immune_score + stromal_score + 1e-10
        result = pd.DataFrame({f"{key_added}": tumor_score / denom})
    else:
        raise ValueError(f"Unknown purity method: {method}")

    n_samples = bulk_data.shape[1]
    result["inference_level"] = "exploratory_sample_level"
    result["valid_for_publication_inference"] = False
    result["method"] = method
    result["result_warning"] = (
        "Expression-based purity estimates are approximate and should be "
        "validated with orthogonal methods (e.g., pathology, CNV)."
    )
    result["n_samples"] = n_samples
    return result


# ---------------------------------------------------------------------------
# Built-in immune signatures for bulk_immune_landscape
# ---------------------------------------------------------------------------

_DEFAULT_IMMUNE_SIGNATURES: Dict[str, List[str]] = {
    "T_cells": ["CD3D", "CD3E", "CD3G", "TRAC"],
    "Cytotoxic": ["GZMA", "GZMB", "PRF1", "GNLY"],
    "Exhausted_T": ["PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT"],
    "Macrophage": ["CD68", "CD14", "CD163", "CSF1R"],
    "Dendritic": ["CD1C", "CLEC9A", "XCR1", "BATF3"],
    "B_plasma": ["CD79A", "CD79B", "MZB1", "JCHAIN"],
    "TLS_like": ["CD19", "CD20", "CXCL13", "CCL19", "CCL21"],
    "Interferon_gamma": ["IFNG", "IRF1", "STAT1", "CXCL9", "CXCL10"],
}


def bulk_immune_landscape(
    bulk_data: pd.DataFrame,
    signature_dict: Optional[Dict[str, List[str]]] = None,
    *,
    method: str = "mean_expression",
    key_added: str = "bulk_immune_landscape",
) -> pd.DataFrame:
    """Score immune signatures in bulk RNA-seq samples.

    Parameters
    ----------
    bulk_data
        Bulk RNA-seq matrix (genes x samples).
    signature_dict
        Custom signature dictionary. If None, uses built-in immune signatures.
    method
        ``"mean_expression"`` (clean-room) or ``"ssgsea_like"`` (rank-based
        enrichment approximation, optional quality).
    key_added
        Prefix for returned columns.

    Returns
    -------
    pd.DataFrame
        Samples x signatures with inference tags.
    """
    signatures = signature_dict or _DEFAULT_IMMUNE_SIGNATURES

    if method == "mean_expression":
        scores = {}
        for name, markers in signatures.items():
            present = [g for g in markers if g in bulk_data.index]
            if not present:
                scores[f"{key_added}_{name}"] = np.nan
                continue
            scores[f"{key_added}_{name}"] = bulk_data.loc[present].mean(axis=0).values
        score_df = pd.DataFrame(scores, index=bulk_data.columns)
    elif method == "ssgsea_like":
        # Clean-room rank-based enrichment score approximation.
        from scipy.stats import rankdata

        scores = {}
        for sample in bulk_data.columns:
            vals = bulk_data[sample].values
            ranks = rankdata(vals, method="average")
            gene_ranks = pd.Series(ranks, index=bulk_data.index)
            sample_scores = {}
            for name, markers in signatures.items():
                present = [g for g in markers if g in gene_ranks.index]
                if not present:
                    sample_scores[f"{key_added}_{name}"] = np.nan
                    continue
                # Enrichment score: mean rank of signature genes normalized to [0, 1]
                sample_scores[f"{key_added}_{name}"] = gene_ranks[present].mean() / len(gene_ranks)
            scores[sample] = sample_scores
        score_df = pd.DataFrame(scores).T
    else:
        raise ValueError(f"Unknown immune landscape method: {method}")

    score_df["inference_level"] = "exploratory_sample_level"
    score_df["valid_for_publication_inference"] = False
    score_df["method"] = method
    score_df["result_warning"] = (
        "Bulk immune signatures are enrichment approximations; interpret in "
        "conjunction with single-cell or IHC evidence."
    )
    return score_df


def associate_tme_with_response(
    proportions_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    response_col: str,
    *,
    method: str = "mannwhitney",
    config: Optional[BulkAbundanceConfig] = None,
    group_col: Optional[str] = None,
    key_added: str = "tme_response_association",
) -> pd.DataFrame:
    """Associate TME compartment proportions with a response variable.

    Parameters
    ----------
    proportions_df
        Samples x compartments or cell types DataFrame.
    metadata_df
        Sample metadata indexed by sample ID.
    response_col
        Column in metadata. If numeric, treated as continuous. If categorical
        binary, used for group comparison.
    method
        ``"mannwhitney"`` (binary), ``"logistic"`` (binary), or
        ``"correlation"`` (continuous).
    config
        Optional abundance-test config for binary group comparisons.
    group_col
        Optional binary group column. If provided with ``method='mannwhitney'``,
        compare the two groups directly.
    key_added
        Audit key under ``adata.uns["sclucid"]["tools"][key_added]`` if an
        AnnData is supplied in future; currently unused.

    Returns
    -------
    pd.DataFrame
        Per-compartment association results with inference tags.
    """
    data = proportions_df.join(metadata_df, how="inner")
    response = data[response_col]
    is_binary = response.dtype == object or response.nunique() <= 2

    results = []
    for compartment in proportions_df.columns:
        vals = data[compartment].dropna()
        if len(vals) < 5:
            continue

        if method in ("mannwhitney", "logistic") and is_binary:
            if group_col is not None:
                groups = data[group_col].dropna().unique()
                if len(groups) != 2:
                    continue
                g1 = vals[data[group_col] == groups[0]]
                g2 = vals[data[group_col] == groups[1]]
            else:
                # Use response_col itself as binary group
                groups = response.dropna().unique()
                if len(groups) != 2:
                    continue
                g1 = vals[response == groups[0]]
                g2 = vals[response == groups[1]]

            if len(g1) < 2 or len(g2) < 2:
                continue

            if method == "mannwhitney":
                from scipy.stats import mannwhitneyu

                stat, pval = mannwhitneyu(g1, g2, alternative="two-sided")
                results.append(
                    {
                        "feature": compartment,
                        "method": "mannwhitney",
                        "statistic": stat,
                        "pval": pval,
                        "effect_size": float(g2.mean() - g1.mean()),
                        "n_group1": len(g1),
                        "n_group2": len(g2),
                    }
                )
            else:  # logistic
                try:
                    from sklearn.linear_model import LogisticRegression
                except ImportError as exc:
                    raise ImportError(
                        "method='logistic' requires scikit-learn. Install with: pip install scikit-learn"
                    ) from exc
                y = (response == groups[1]).astype(int).loc[vals.index].values
                x = vals.values.reshape(-1, 1)
                model = LogisticRegression(solver="lbfgs")
                model.fit(x, y)
                coef = float(model.coef_[0][0])
                # Approximate p-value via Wald test
                pred = model.predict_proba(x)[:, 1]
                se = np.sqrt(np.sum((x.ravel() ** 2) * pred * (1 - pred)))
                z = coef / max(se, 1e-10)
                from scipy.stats import norm

                pval = 2 * (1 - norm.cdf(abs(z)))
                results.append(
                    {
                        "feature": compartment,
                        "method": "logistic",
                        "statistic": z,
                        "pval": pval,
                        "effect_size": coef,
                        "n_group1": int((y == 0).sum()),
                        "n_group2": int((y == 1).sum()),
                    }
                )
        elif method == "correlation" or not is_binary:
            from scipy.stats import pearsonr, spearmanr

            mask = response.notna() & vals.notna()
            if mask.sum() < 5:
                continue
            corr, pval = spearmanr(vals[mask], response[mask])
            results.append(
                {
                    "feature": compartment,
                    "method": "spearman",
                    "statistic": corr,
                    "pval": pval,
                    "effect_size": corr,
                    "n_samples": int(mask.sum()),
                }
            )

    result = pd.DataFrame(results)
    if result.empty:
        return result

    from statsmodels.stats.multitest import multipletests

    _, qvals, _, _ = multipletests(result["pval"].values, method="fdr_bh")
    result["pvals_adj"] = qvals

    replicate_met = result.get("n_group1", pd.Series([0])).ge(3) & result.get("n_group2", pd.Series([0])).ge(3)
    if "n_samples" in result.columns:
        replicate_met = result["n_samples"].ge(10)

    result["inference_level"] = replicate_met.map(lambda x: "sample_level" if x else "exploratory_trait_association")
    result["valid_for_publication_inference"] = replicate_met
    result["replicate_requirement_met"] = replicate_met
    result["diagnostic_status"] = "passed"
    result["result_warning"] = replicate_met.map(
        lambda x: (
            None
            if x
            else "Insufficient samples for formal response association; treat as exploratory."
        )
    )
    return result.sort_values("pval")
