"""Infer CNV from single-cell RNA-seq data.

This module implements a CopyKAT-like CNV inference pipeline that:
1. Uses raw counts as input (with automatic detection and layer fallback)
2. Integrates genomic coordinates for chromosome-aware analysis
3. Filters artifact genes (mitochondrial, ribosomal, cell-cycle)
4. Applies chromosome-level smoothing with uniform_filter1d
5. Robustly scales against reference cells using MAD
6. Computes per-chromosome CNV scores and global scores
7. Predicts aneuploid / diploid cells with conservative thresholds
8. Provides comprehensive visualization functions

The pipeline is designed to be self-contained and does not require
external CNV packages. It is suitable for quick tumor-cell identification
inside analysis workflows.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.ndimage import uniform_filter1d
from scipy.sparse import issparse

log = logging.getLogger(__name__)

# -----------------------------------------------------------------------------
# Artifact gene patterns
# -----------------------------------------------------------------------------
_MT_PATTERNS = (r"^(mt-|MT-|Mt-)",)
_RIBO_PATTERNS = (r"^(Rpl|Rps|RPL|RPS)",)

# Default cell-cycle genes (can be overridden)
_DEFAULT_CELL_CYCLE_GENES: set = set()

# Default standard chromosomes (human)
_DEFAULT_STANDARD_CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]

# Default standard chromosomes (mouse)
_MOUSE_STANDARD_CHROMS = [str(i) for i in range(1, 20)] + ["X", "Y"]


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
def _resolve_expression_matrix(
    adata: AnnData,
    layer: Optional[str] = None,
    require_counts: bool = True,
) -> Tuple[np.ndarray, AnnData]:
    """Return expression matrix and a copy of adata with the chosen matrix set to X.

    If ``layer`` is provided, use that layer. Otherwise use ``adata.X``.
    When ``require_counts=True``, validates that the data looks like raw counts
    (integer values, max > 20) and attempts to fall back to known count layers.
    """
    if layer is not None:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in adata.layers")
        X = adata.layers[layer]
        work = adata.copy()
        work.X = X
        return _to_dense_if_sparse(work.X), work

    X = adata.X
    work = adata.copy()

    if require_counts:
        sample = _to_dense_if_sparse(X[: min(100, X.shape[0]), : min(100, X.shape[1])])
        max_val = float(np.max(sample))
        has_decimals = np.any(sample % 1 != 0)

        if max_val < 20 or has_decimals:
            log.warning("Data appears normalized/log-transformed. Searching for count layers...")
            for candidate in ("counts", "raw_counts", "raw"):
                if candidate in adata.layers:
                    log.info(f"Using '{candidate}' layer for CNV inference")
                    work.X = adata.layers[candidate].copy()
                    return _to_dense_if_sparse(work.X), work
            raise ValueError(
                "Raw count data required for CNV inference. "
                "Provide a count layer via the 'layer' argument."
            )

    return _to_dense_if_sparse(work.X), work


def _to_dense_if_sparse(X: Any) -> np.ndarray:
    if issparse(X):
        return X.toarray()
    return np.asarray(X)


def _mad(x: np.ndarray, axis: int = 0, eps: float = 1e-6) -> np.ndarray:
    """Median absolute deviation along an axis."""
    med = np.median(x, axis=axis, keepdims=True)
    mad_val = np.median(np.abs(x - med), axis=axis, keepdims=True)
    return np.squeeze(mad_val, axis=axis) + eps


def _match_gene_symbols(var: pd.DataFrame, var_names: pd.Index) -> pd.Series:
    """Extract gene symbols from var, falling back to var_names."""
    symbol_cols = [c for c in ("gene_name", "gene_symbol", "symbol", "external_gene_name") if c in var.columns]
    if symbol_cols:
        return var[symbol_cols[0]].astype(str)
    return pd.Series(var_names, index=var_names).astype(str)


def _attach_genomic_coordinates(
    adata: AnnData,
    gene_info: Optional[Union[pd.DataFrame, str, Path]] = None,
) -> AnnData:
    """Attach chromosome / start / end to adata.var from a gene info table.

    Parameters
    ----------
    adata : AnnData
        Must have gene names aligned with ``gene_info``.
    gene_info : DataFrame or path, optional
        Expected columns after optional renaming:
        ``chromosome`` (or ``chromosome_name``),
        ``start`` (or ``start_position``),
        ``end`` (or ``end_position``),
        index or ``gene_name`` column matching adata.var_names.
    """
    if gene_info is None:
        if not all(c in adata.var.columns for c in ("chromosome", "start", "end")):
            log.warning(
                "Genomic coordinates not available. CNV inference will run without "
                "chromosome-aware smoothing and ordering."
            )
        return adata

    if isinstance(gene_info, (str, Path)):
        gene_info = pd.read_csv(gene_info, sep="\t")

    gene_info = gene_info.copy()

    # Normalize column names
    rename_map = {
        "chromosome_name": "chromosome",
        "start_position": "start",
        "end_position": "end",
        "external_gene_name": "gene_name",
        "gene_symbol": "gene_name",
        "symbol": "gene_name",
    }
    gene_info = gene_info.rename(columns={k: v for k, v in rename_map.items() if k in gene_info.columns})

    # Filter to protein coding when biotype available
    if "gene_biotype" in gene_info.columns:
        gene_info = gene_info[gene_info["gene_biotype"] == "protein_coding"].copy()

    # Deduplicate
    if "gene_name" in gene_info.columns:
        gene_info = gene_info.drop_duplicates(subset=["gene_name"], keep="first")
        gene_info = gene_info.set_index("gene_name")
    else:
        # Assume index is gene names
        gene_info = gene_info.loc[~gene_info.index.duplicated(keep="first")]

    for col in ("chromosome", "start", "end"):
        if col in gene_info.columns:
            adata.var[col] = gene_info.reindex(adata.var_names)[col]

    return adata


# -----------------------------------------------------------------------------
# CNVAnalyzer
# -----------------------------------------------------------------------------
class CNVAnalyzer:
    """Analyze copy number variations from scRNA-seq data (CopyKAT-like).

    Parameters
    ----------
    gene_order : pd.DataFrame, optional
        DataFrame with columns ``chromosome``, ``start``, ``end`` indexed by gene name.
        If None, coordinates must already be present in ``adata.var``.
    window_size : int
        Smoothing window size for chromosome-level uniform filtering.
    min_cells_per_gene : int, optional
        Minimum number of cells expressing a gene to retain it.
        Defaults to 5% of total cells.
    min_counts_per_gene : int, optional
        Minimum total counts per gene. Defaults to ``min_cells_per_gene``.
    filter_mt : bool
        Remove mitochondrial genes before CNV inference.
    filter_ribo : bool
        Remove ribosomal genes before CNV inference.
    filter_cc : bool
        Remove cell-cycle genes before CNV inference.
    cc_genes : set, optional
        Custom cell-cycle gene set. If None and ``filter_cc=True``,
        an empty set is used (no genes removed unless provided).
    layer : str, optional
        Layer containing raw counts. If None, ``adata.X`` is used
        (with automatic count validation).
    target_sum : float
        Target sum for ``sc.pp.normalize_total``.
    standard_chroms : list, optional
        List of valid chromosome names. Defaults to human 1-22, X, Y.
    """

    def __init__(
        self,
        gene_order: Optional[pd.DataFrame] = None,
        window_size: int = 101,
        min_cells_per_gene: Optional[int] = None,
        min_counts_per_gene: Optional[int] = None,
        filter_mt: bool = True,
        filter_ribo: bool = True,
        filter_cc: bool = False,
        cc_genes: Optional[set] = None,
        layer: Optional[str] = None,
        target_sum: float = 1e4,
        standard_chroms: Optional[List[str]] = None,
    ):
        self.gene_order = gene_order
        self.window_size = window_size
        self.min_cells_per_gene = min_cells_per_gene
        self.min_counts_per_gene = min_counts_per_gene
        self.filter_mt = filter_mt
        self.filter_ribo = filter_ribo
        self.filter_cc = filter_cc
        self.cc_genes = cc_genes or _DEFAULT_CELL_CYCLE_GENES.copy()
        self.layer = layer
        self.target_sum = target_sum
        self.standard_chroms = standard_chroms or _DEFAULT_STANDARD_CHROMS

        # Fitted attributes
        self.cnv_matrix_: Optional[np.ndarray] = None
        self.tumor_scores_: Optional[pd.Series] = None
        self.chromosome_scores_: Optional[pd.DataFrame] = None
        self.adata_work_: Optional[AnnData] = None
        self.reference_mask_: Optional[np.ndarray] = None
        self.threshold_: Optional[float] = None
        self.Z_: Optional[np.ndarray] = None

    # -------------------------------------------------------------------------
    # Internal pipeline steps
    # -------------------------------------------------------------------------
    def _prepare_data(self, adata: AnnData) -> AnnData:
        """Resolve expression matrix and attach genomic coordinates."""
        _, work = _resolve_expression_matrix(adata, layer=self.layer, require_counts=True)
        work = _attach_genomic_coordinates(work, gene_info=self.gene_order)
        return work

    def _filter_genes(self, adata: AnnData) -> AnnData:
        """Filter genes by expression and genomic validity."""
        X = adata.X
        if issparse(X):
            gene_counts = np.array(X.sum(axis=0)).flatten()
            cells_per_gene = np.array((X > 0).sum(axis=0)).flatten()
        else:
            gene_counts = np.asarray(X.sum(axis=0)).flatten()
            cells_per_gene = np.asarray((X > 0).sum(axis=0)).flatten()

        min_cells = self.min_cells_per_gene or max(1, int(0.05 * adata.n_obs))
        min_counts = self.min_counts_per_gene or min_cells

        keep = (cells_per_gene >= min_cells) & (gene_counts >= min_counts)

        has_coords = "chromosome" in adata.var.columns
        if has_coords:
            coord_ok = adata.var[["chromosome", "start", "end"]].notna().all(axis=1).values
            valid_chrom_mask = adata.var["chromosome"].astype(str).isin(self.standard_chroms).values
            keep = keep & coord_ok & valid_chrom_mask

        adata = adata[:, keep].copy()

        # Drop chromosomes with too few genes
        if has_coords and "chromosome" in adata.var.columns:
            chr_counts = adata.var["chromosome"].value_counts()
            valid_chroms = chr_counts[chr_counts >= 30].index.tolist()
            adata = adata[:, adata.var["chromosome"].isin(valid_chroms)].copy()

        return adata

    def _filter_artifact_genes(self, adata: AnnData) -> AnnData:
        """Remove mitochondrial, ribosomal, and optionally cell-cycle genes."""
        gene_symbols = _match_gene_symbols(adata.var, adata.var_names)

        masks: List[np.ndarray] = []
        if self.filter_mt:
            mt_mask = gene_symbols.str.match(_MT_PATTERNS[0], na=False)
            masks.append(np.asarray(mt_mask))
            log.info(f"Filtering {int(mt_mask.sum())} mitochondrial genes")
        if self.filter_ribo:
            ribo_mask = gene_symbols.str.match(_RIBO_PATTERNS[0], na=False)
            masks.append(np.asarray(ribo_mask))
            log.info(f"Filtering {int(ribo_mask.sum())} ribosomal genes")
        if self.filter_cc and self.cc_genes:
            cc_mask = gene_symbols.isin(self.cc_genes).values
            masks.append(np.asarray(cc_mask))
            log.info(f"Filtering {int(cc_mask.sum())} cell-cycle genes")

        if not masks:
            return adata

        drop_mask = masks[0]
        for m in masks[1:]:
            drop_mask = drop_mask | m

        adata = adata[:, ~drop_mask].copy()
        log.info(f"Retained {adata.n_vars} genes after artifact filtering")
        return adata

    def _order_genes_by_chromosome(self, adata: AnnData) -> AnnData:
        """Sort genes by chromosome and genomic position."""
        if "chromosome" not in adata.var.columns or "start" not in adata.var.columns:
            log.warning("Genomic coordinates unavailable; skipping chromosome ordering")
            return adata
        adata.var["chromosome"] = pd.Categorical(
            adata.var["chromosome"].astype(str),
            categories=self.standard_chroms,
            ordered=True,
        )
        adata.var["start"] = pd.to_numeric(adata.var["start"], errors="coerce")
        order = adata.var.sort_values(["chromosome", "start"]).index
        return adata[:, order].copy()

    def _normalize_expression(self, adata: AnnData) -> np.ndarray:
        """Log1p CPM normalization. Returns dense matrix."""
        sc = _optional_scanpy()
        if sc is None:
            # Fallback manual normalization if scanpy unavailable
            X = _to_dense_if_sparse(adata.X)
            lib_sizes = X.sum(axis=1, keepdims=True)
            lib_sizes = np.where(lib_sizes == 0, 1, lib_sizes)
            X_norm = X / lib_sizes * self.target_sum
            return np.log1p(X_norm)

        # Use scanpy for normalization
        work = adata.copy()
        sc.pp.normalize_total(work, target_sum=self.target_sum)
        sc.pp.log1p(work)
        return _to_dense_if_sparse(work.X)

    def _smooth_by_chromosome(self, X: np.ndarray, chromosomes: np.ndarray) -> np.ndarray:
        """Apply uniform_filter1d per chromosome, or globally if no coordinates."""
        if chromosomes is None or len(chromosomes) == 0:
            log.warning("No chromosome information; applying global smoothing")
            return uniform_filter1d(
                X.astype(np.float32), size=self.window_size, axis=1, mode="nearest"
            )

        unique_chroms = pd.unique(chromosomes)
        X_smooth = np.zeros_like(X, dtype=np.float32)

        for chrom in unique_chroms:
            idx = np.where(chromosomes == chrom)[0]
            if len(idx) < 30:
                X_smooth[:, idx] = X[:, idx].astype(np.float32)
                continue
            X_chr = X[:, idx].astype(np.float32)
            X_smooth[:, idx] = uniform_filter1d(
                X_chr, size=self.window_size, axis=1, mode="nearest"
            )

        return X_smooth

    def _robust_scale_to_reference(
        self, X_smooth: np.ndarray, ref_mask: np.ndarray
    ) -> np.ndarray:
        """Center and scale using reference-cell median and MAD."""
        if ref_mask.sum() < 50:
            raise ValueError(
                f"Too few reference cells ({int(ref_mask.sum())}). "
                "Provide at least 50 reference cells."
            )

        X_ref = X_smooth[ref_mask, :]
        ref_mean = X_ref.mean(axis=0)
        ref_mad = _mad(X_ref, axis=0)

        Z = (X_smooth - ref_mean) / ref_mad
        Z = np.clip(Z, -10, 10)
        return Z

    def _compute_chromosome_scores(
        self, Z: np.ndarray, chromosomes: Optional[np.ndarray]
    ) -> pd.DataFrame:
        """Compute per-chromosome median absolute Z scores."""
        index = self.adata_work_.obs_names if self.adata_work_ is not None else None

        if chromosomes is None or len(chromosomes) == 0:
            # Global score only when coordinates unavailable
            return pd.DataFrame(
                {"global_score": np.median(np.abs(Z), axis=1)},
                index=index,
            )

        unique_chroms = pd.unique(chromosomes)
        scores: Dict[str, np.ndarray] = {}

        for chrom in unique_chroms:
            idx = np.where(chromosomes == chrom)[0]
            if len(idx) < 30:
                continue
            scores[f"chr{chrom}_score"] = np.median(np.abs(Z[:, idx]), axis=1)

        return pd.DataFrame(scores, index=index)

    # -------------------------------------------------------------------------
    # Public API
    # -------------------------------------------------------------------------
    def fit(
        self,
        adata: AnnData,
        reference_cells: Optional[Union[str, List[str]]] = None,
        reference_key: str = "cell_type",
    ) -> "CNVAnalyzer":
        """Infer CNV from expression data using the CopyKAT-like pipeline.

        Parameters
        ----------
        adata : AnnData
            Single-cell expression data with raw counts.
        reference_cells : str or list, optional
            Cell type labels to use as diploid reference. If None,
            all cells are used as reference (not recommended).
        reference_key : str
            Column in ``adata.obs`` containing cell type labels.

        Returns
        -------
        CNVAnalyzer
            Fitted analyzer with results stored in attributes.
        """
        log.info("Starting CNV inference (CopyKAT-like pipeline)")

        # Step 1: prepare data
        work = self._prepare_data(adata)
        log.info(f"Initial: {work.n_obs} cells x {work.n_vars} genes")

        # Step 2: filter genes
        work = self._filter_genes(work)
        log.info(f"After genomic filtering: {work.n_vars} genes")

        # Step 3: filter artifacts
        work = self._filter_artifact_genes(work)

        # Step 4: order by chromosome
        work = self._order_genes_by_chromosome(work)

        # Step 5: normalize
        X_norm = self._normalize_expression(work)
        log.info("Normalized expression (log1p CPM)")

        # Step 6: chromosome-level smoothing
        chromosomes = (
            work.var["chromosome"].astype(str).values
            if "chromosome" in work.var.columns
            else None
        )
        X_smooth = self._smooth_by_chromosome(X_norm, chromosomes)
        log.info(f"Smoothed with window_size={self.window_size}")

        # Step 7: reference mask
        if reference_cells is not None:
            if isinstance(reference_cells, str):
                ref_mask = (work.obs[reference_key] == reference_cells).values
            else:
                ref_mask = work.obs[reference_key].isin(reference_cells).values
        else:
            log.warning("No reference_cells provided; using all cells as reference")
            ref_mask = np.ones(work.n_obs, dtype=bool)

        self.reference_mask_ = ref_mask
        n_ref = int(ref_mask.sum())
        log.info(f"Reference cells: {n_ref} / {work.n_obs} ({n_ref / work.n_obs * 100:.1f}%)")

        # Step 8: robust scaling
        Z = self._robust_scale_to_reference(X_smooth, ref_mask)
        self.Z_ = Z
        log.info("Robust scaling complete (centered to reference, MAD-scaled)")

        # Step 9: compute scores
        self.chromosome_scores_ = self._compute_chromosome_scores(Z, chromosomes)
        self.tumor_scores_ = self.chromosome_scores_.median(axis=1)

        # Step 10: store CNV matrix
        self.cnv_matrix_ = Z
        self.adata_work_ = work

        # Step 11: compute extreme fraction
        extreme_frac = (np.abs(Z) > 2.0).mean(axis=1)
        work.obs["cnv_score"] = self.tumor_scores_.values
        work.obs["cnv_extreme_frac"] = extreme_frac

        log.info("CNV inference complete")
        return self

    def predict_aneuploid(
        self,
        threshold_mad: float = 3.0,
        use_gmm: bool = True,
        key_added: str = "predicted_class",
    ) -> pd.Series:
        """Predict aneuploid cells based on CNV scores.

        Uses a conservative threshold based on reference-cell CNV score
        distribution (median + ``threshold_mad`` * MAD). Optionally refines
        with a 2-component Gaussian Mixture Model on
        (cnv_score, cnv_extreme_frac).

        Parameters
        ----------
        threshold_mad : float
            Number of MADs above reference median for threshold-based call.
        use_gmm : bool
            Whether to run GMM refinement when sklearn is available.
        key_added : str
            Column name for predictions in ``adata_work_.obs``.

        Returns
        -------
        pd.Series
            Predicted class (``'aneuploid'`` or ``'diploid'``).
        """
        if self.tumor_scores_ is None:
            raise ValueError("Must call fit() first")

        ref_scores = self.tumor_scores_.values[self.reference_mask_]
        ref_med = float(np.median(ref_scores))
        ref_mad_score = float(np.median(np.abs(ref_scores - ref_med))) + 1e-6

        thr = ref_med + threshold_mad * ref_mad_score
        self.threshold_ = thr

        pred_thr = np.where(self.tumor_scores_.values > thr, "aneuploid", "diploid")

        pred_gmm = None
        if use_gmm:
            try:
                from sklearn.mixture import GaussianMixture

                feats = np.column_stack([
                    self.tumor_scores_.values,
                    self.adata_work_.obs["cnv_extreme_frac"].values,
                ])
                gmm = GaussianMixture(n_components=2, random_state=42)
                gmm_labels = gmm.fit_predict(feats)
                means = [
                    self.tumor_scores_.values[gmm_labels == i].mean()
                    for i in (0, 1)
                ]
                aneu_comp = int(np.argmax(means))
                pred_gmm = np.where(gmm_labels == aneu_comp, "aneuploid", "diploid")
            except Exception as exc:
                log.warning(f"GMM refinement skipped: {exc}")

        # Prefer GMM if available; otherwise threshold
        final = pd.Series(
            pred_gmm if pred_gmm is not None else pred_thr,
            index=self.adata_work_.obs_names,
            name=key_added,
        )

        if self.adata_work_ is not None:
            self.adata_work_.obs[f"{key_added}_thr"] = pd.Categorical(
                pred_thr, categories=["diploid", "aneuploid"]
            )
            if pred_gmm is not None:
                self.adata_work_.obs[f"{key_added}_gmm"] = pd.Categorical(
                    pred_gmm, categories=["diploid", "aneuploid"]
                )
            self.adata_work_.obs[key_added] = final

        n_aneu = int((final == "aneuploid").sum())
        log.info(
            f"Aneuploid prediction: threshold={thr:.4f}, "
            f"aneuploid={n_aneu} / {len(final)} ({n_aneu / len(final) * 100:.1f}%)"
        )
        return final

    def predict_tumor_cells(
        self,
        threshold: Optional[float] = None,
    ) -> pd.Series:
        """Predict tumor cells based on CNV scores.

        This is a convenience wrapper that calls ``predict_aneuploid``
        when the analyzer has been fitted. If ``threshold`` is None,
        the MAD-based threshold from ``predict_aneuploid`` is used.

        Parameters
        ----------
        threshold : float, optional
            Direct threshold on CNV score. If None, uses the auto-derived
            threshold from ``predict_aneuploid``.

        Returns
        -------
        pd.Series
            Boolean series indicating tumor (aneuploid) cells.
        """
        if self.tumor_scores_ is None:
            raise ValueError("Must call fit() first")

        if threshold is not None:
            return pd.Series(
                self.tumor_scores_.values > threshold,
                index=self.adata_work_.obs_names if self.adata_work_ is not None else None,
            )

        # Use aneuploid prediction if available
        if self.adata_work_ is not None and "predicted_class" in self.adata_work_.obs.columns:
            return pd.Series(
                self.adata_work_.obs["predicted_class"].astype(str) == "aneuploid",
                index=self.adata_work_.obs_names,
            )

        # Fallback: run prediction with defaults
        pred = self.predict_aneuploid()
        return pd.Series(pred == "aneuploid", index=pred.index)

    def get_results(self) -> Dict[str, Any]:
        """Return a dictionary of all computed results."""
        if self.adata_work_ is None:
            raise ValueError("Must call fit() first")

        results = {
            "cnv_matrix": self.cnv_matrix_,
            "cnv_score": self.tumor_scores_,
            "chromosome_scores": self.chromosome_scores_,
            "reference_mask": self.reference_mask_,
            "threshold": self.threshold_,
            "Z": self.Z_,
        }
        return results


# -----------------------------------------------------------------------------
# Convenience functions (backward-compatible)
# -----------------------------------------------------------------------------
def infer_cnv(
    adata: AnnData,
    reference_cells: Optional[Union[str, List[str]]] = None,
    reference_key: str = "cell_type",
    gene_order: Optional[Union[pd.DataFrame, str, Path]] = None,
    layer: Optional[str] = None,
    key_added: str = "cnv",
    window_size: int = 101,
    filter_mt: bool = True,
    filter_ribo: bool = True,
    filter_cc: bool = False,
    cc_genes: Optional[set] = None,
    standard_chroms: Optional[List[str]] = None,
    predict_aneuploid: bool = True,
    threshold_mad: float = 3.0,
    use_gmm: bool = True,
    copy: bool = False,
) -> AnnData:
    """Infer CNV from single-cell expression data (CopyKAT-like pipeline).

    Parameters
    ----------
    adata : AnnData
        Single-cell expression data. Raw counts are required.
    reference_cells : str or list, optional
        Cell types to use as normal (diploid) reference.
    reference_key : str
        Column in ``adata.obs`` with cell type labels.
    gene_order : DataFrame or path, optional
        Gene position information. If None, must be in ``adata.var``.
    layer : str, optional
        Layer containing raw counts. If None, uses ``adata.X``
        with automatic count validation.
    key_added : str
        Prefix for storing results.
    window_size : int
        Smoothing window size.
    filter_mt, filter_ribo, filter_cc : bool
        Whether to filter artifact genes.
    cc_genes : set, optional
        Custom cell-cycle genes.
    standard_chroms : list, optional
        Valid chromosome names (default human 1-22, X, Y).
    predict_aneuploid : bool
        Whether to run aneuploid prediction after CNV inference.
    threshold_mad : float
        MAD multiplier for threshold-based aneuploid call.
    use_gmm : bool
        Whether to refine with GMM.
    copy : bool
        Return a copy of adata.

    Returns
    -------
    AnnData
        Annotated data with CNV information stored in ``obsm``, ``obs``,
        and ``uns``.
    """
    if copy:
        adata = adata.copy()

    analyzer = CNVAnalyzer(
        gene_order=gene_order,
        window_size=window_size,
        filter_mt=filter_mt,
        filter_ribo=filter_ribo,
        filter_cc=filter_cc,
        cc_genes=cc_genes,
        layer=layer,
        standard_chroms=standard_chroms,
    )
    analyzer.fit(adata, reference_cells=reference_cells, reference_key=reference_key)

    work = analyzer.adata_work_

    # Store results back to original adata (aligning by obs_names)
    # CNV matrix: full Z score matrix (cells x genes)
    adata.obsm[f"X_{key_added}"] = analyzer.cnv_matrix_
    adata.obs[f"{key_added}_score"] = analyzer.tumor_scores_.reindex(adata.obs_names)
    adata.obs[f"{key_added}_extreme_frac"] = work.obs["cnv_extreme_frac"].reindex(adata.obs_names)

    # Chromosome-level scores
    for col in analyzer.chromosome_scores_.columns:
        adata.obs[col] = analyzer.chromosome_scores_[col].reindex(adata.obs_names)

    # Aneuploid predictions
    if predict_aneuploid:
        analyzer.predict_aneuploid(threshold_mad=threshold_mad, use_gmm=use_gmm)
        for suffix in ("predicted_class", "predicted_class_thr", "predicted_class_gmm"):
            col = suffix
            if col in work.obs.columns:
                adata.obs[f"{key_added}_{suffix}"] = work.obs[col].reindex(adata.obs_names)

    # Store parameters
    adata.uns[f"{key_added}_params"] = {
        "window_size": window_size,
        "reference_key": reference_key,
        "reference_cells": list(reference_cells) if isinstance(reference_cells, list) else reference_cells,
        "filter_mt": filter_mt,
        "filter_ribo": filter_ribo,
        "filter_cc": filter_cc,
        "threshold_mad": threshold_mad,
        "use_gmm": use_gmm,
    }

    # Summary statistics
    if predict_aneuploid and f"{key_added}_predicted_class" in adata.obs.columns:
        calls = adata.obs[f"{key_added}_predicted_class"].astype(str)
        adata.uns[f"{key_added}_summary"] = {
            "n_aneuploid": int((calls == "aneuploid").sum()),
            "n_diploid": int((calls == "diploid").sum()),
            "aneuploid_fraction": float((calls == "aneuploid").mean()),
            "mean_cnv_score": float(adata.obs[f"{key_added}_score"].mean()),
            "threshold": analyzer.threshold_,
        }

    log.info(f"CNV inference complete. Results stored with prefix '{key_added}'")
    return adata


def find_tumor_cells(
    adata: AnnData,
    method: str = "cnv_score",
    threshold: float = 0.5,
    key: str = "cnv",
) -> pd.Series:
    """Identify tumor cells based on CNV analysis.

    Parameters
    ----------
    adata : AnnData
        Data with CNV information.
    method : str
        Method for tumor identification:
        - ``"cnv_score"``: direct threshold on CNV score
        - ``"predicted_class"``: use aneuploid prediction
        - ``"clustering"``: K-means on CNV matrix
    threshold : float
        Threshold for ``"cnv_score"`` method.
    key : str
        Key prefix for CNV data in adata.

    Returns
    -------
    pd.Series
        Boolean series indicating tumor cells.
    """
    if method == "cnv_score":
        score_col = f"{key}_score"
        if score_col not in adata.obs.columns:
            raise KeyError(f"'{score_col}' not found in adata.obs")
        return adata.obs[score_col] > threshold

    elif method == "predicted_class":
        class_col = f"{key}_predicted_class"
        if class_col not in adata.obs.columns:
            raise KeyError(f"'{class_col}' not found in adata.obs. Run infer_cnv with predict_aneuploid=True.")
        return adata.obs[class_col].astype(str) == "aneuploid"

    elif method == "clustering":
        from sklearn.cluster import KMeans

        cnv_data = adata.obsm[f"X_{key}"]
        kmeans = KMeans(n_clusters=2, random_state=42, n_init=10).fit(cnv_data)
        cluster_means = [np.abs(cnv_data[kmeans.labels_ == i]).mean() for i in range(2)]
        tumor_cluster = np.argmax(cluster_means)
        return pd.Series(kmeans.labels_ == tumor_cluster, index=adata.obs_names)

    else:
        raise ValueError(f"Unknown method: {method}")


def identify_clones(
    adata: AnnData,
    cnv_key: str = "cnv",
    n_clusters: int = 5,
    method: str = "hierarchical",
) -> pd.Series:
    """Identify tumor clones based on CNV patterns.

    Parameters
    ----------
    adata : AnnData
        Data with CNV information.
    cnv_key : str
        Key for CNV data.
    n_clusters : int
        Number of expected clones.
    method : str
        Clustering method ("hierarchical", "kmeans", "leiden").

    Returns
    -------
    pd.Series
        Clone assignments for each cell.
    """
    cnv_data = adata.obsm[f"X_{cnv_key}"]

    if method == "hierarchical":
        from sklearn.cluster import AgglomerativeClustering

        clusterer = AgglomerativeClustering(n_clusters=n_clusters)
    elif method == "kmeans":
        from sklearn.cluster import KMeans

        clusterer = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    elif method == "leiden":
        import scanpy as sc

        adata_cnv = AnnData(X=cnv_data)
        sc.pp.neighbors(adata_cnv, n_neighbors=15)
        sc.tl.leiden(adata_cnv, resolution=1.0)
        return adata_cnv.obs["leiden"]
    else:
        raise ValueError(f"Unknown method: {method}")

    labels = clusterer.fit_predict(cnv_data)
    return pd.Series(labels, index=adata.obs_names, name="clone")


def calculate_cnv_score(
    adata: AnnData,
    cnv_key: str = "cnv",
    method: str = "mean_absolute",
) -> pd.Series:
    """Calculate overall CNV burden score.

    Parameters
    ----------
    adata : AnnData
        Data with CNV information.
    cnv_key : str
        Key for CNV data.
    method : str
        Scoring method:
        - ``"mean_absolute"``: mean absolute Z score
        - ``"median_absolute"``: median absolute Z score (default from pipeline)
        - ``"variance"``: variance of Z scores
        - ``"extreme_frac"``: fraction of |Z| > 2
        - ``"gini"``: Gini coefficient of absolute Z scores

    Returns
    -------
    pd.Series
        CNV scores per cell.
    """
    cnv_data = adata.obsm[f"X_{cnv_key}"]

    if method == "mean_absolute":
        scores = np.abs(cnv_data).mean(axis=1)
    elif method == "median_absolute":
        scores = np.median(np.abs(cnv_data), axis=1)
    elif method == "variance":
        scores = np.var(cnv_data, axis=1)
    elif method == "extreme_frac":
        scores = (np.abs(cnv_data) > 2.0).mean(axis=1)
    elif method == "gini":
        scores = np.array([_gini_coefficient(np.abs(x)) for x in cnv_data])
    else:
        raise ValueError(f"Unknown method: {method}")

    return pd.Series(scores, index=adata.obs_names, name="cnv_burden")


def _gini_coefficient(x: np.ndarray) -> float:
    """Calculate Gini coefficient for inequality."""
    sorted_x = np.sort(x)
    n = len(x)
    cumsum = np.cumsum(sorted_x)
    if cumsum[-1] == 0:
        return 0.0
    return float((n + 1 - 2 * np.sum(cumsum) / cumsum[-1]) / n)


def _optional_scanpy() -> Optional[Any]:
    """Return scanpy module if available, else None."""
    try:
        import scanpy

        return scanpy
    except ImportError:
        return None


# -----------------------------------------------------------------------------
# Visualization functions
# -----------------------------------------------------------------------------
def plot_cnv_distribution(
    adata: AnnData,
    cnv_key: str = "cnv",
    groupby: Optional[str] = None,
    reference_label: Optional[str] = None,
    figsize: Tuple[float, float] = (14, 5),
    save: Optional[str] = None,
    show: bool = True,
) -> Any:
    """Plot CNV score distribution by reference vs non-reference and by predicted class.

    Parameters
    ----------
    adata : AnnData
        Data with CNV inference results.
    cnv_key : str
        Key prefix for CNV columns.
    groupby : str, optional
        Column for stratifying the reference group. If provided,
        ``reference_label`` should be a value in this column.
    reference_label : str, optional
        Value indicating reference cells in ``groupby`` or
        ``adata.obs[f"{cnv_key}_is_reference"]``.
    figsize : tuple
        Figure size.
    save : str, optional
        Path to save figure.
    show : bool
        Whether to call ``plt.show()``.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    score_col = f"{cnv_key}_score"
    class_col = f"{cnv_key}_predicted_class"

    if score_col not in adata.obs.columns:
        raise KeyError(f"'{score_col}' not found in adata.obs")

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    scores = adata.obs[score_col].values

    # Panel 1: reference vs non-reference
    if groupby is not None and reference_label is not None:
        ref_mask = adata.obs[groupby] == reference_label
    elif f"{cnv_key}_is_reference" in adata.obs.columns:
        ref_mask = adata.obs[f"{cnv_key}_is_reference"].astype(bool)
    else:
        ref_mask = None

    if ref_mask is not None and ref_mask.sum() > 0:
        axes[0].hist(
            scores[ref_mask], bins=50, alpha=0.6, label="Reference",
            color="#3498db", density=True,
        )
        axes[0].hist(
            scores[~ref_mask], bins=50, alpha=0.6, label="Non-reference",
            color="#e74c3c", density=True,
        )
    else:
        axes[0].hist(scores, bins=50, alpha=0.6, color="#3498db", density=True)

    # Threshold line
    summary = adata.uns.get(f"{cnv_key}_summary", {})
    thr = summary.get("threshold")
    if thr is not None:
        axes[0].axvline(thr, color="black", linestyle="--", linewidth=2, label=f"Thr={thr:.3f}")
        axes[1].axvline(thr, color="black", linestyle="--", linewidth=2, label=f"Thr={thr:.3f}")

    axes[0].set_xlabel("CNV Score")
    axes[0].set_ylabel("Density")
    axes[0].set_title("CNV Score: Reference vs Others")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    # Panel 2: by predicted class
    if class_col in adata.obs.columns:
        for pred_class, color in [("diploid", "#2ecc71"), ("aneuploid", "#e74c3c")]:
            mask = adata.obs[class_col].astype(str) == pred_class
            if mask.sum() == 0:
                continue
            axes[1].hist(
                scores[mask], bins=50, alpha=0.6,
                label=f"{pred_class} (n={int(mask.sum())})",
                density=True, color=color,
            )
    else:
        axes[1].hist(scores, bins=50, alpha=0.6, color="#3498db", density=True)

    axes[1].set_xlabel("CNV Score")
    axes[1].set_ylabel("Density")
    axes[1].set_title("CNV Score by Predicted Class")
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    plt.tight_layout()

    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")

    if show:
        plt.show()

    return fig


def plot_cnv_heatmap(
    adata: AnnData,
    cnv_key: str = "cnv",
    n_show: int = 800,
    random_state: int = 42,
    sort_by: str = "cnv_score",
    vmin: float = -3,
    vmax: float = 3,
    figsize: Tuple[float, float] = (18, 7),
    cmap: str = "bwr",
    save: Optional[str] = None,
    show: bool = True,
) -> Any:
    """Plot CNV Z-score heatmap along the genome.

    Parameters
    ----------
    adata : AnnData
        Data with CNV matrix in ``obsm[f"X_{cnv_key}"]``.
    cnv_key : str
        Key prefix.
    n_show : int
        Number of cells to sample for visualization.
    random_state : int
        Random seed for sampling.
    sort_by : str
        How to sort cells ("cnv_score", "reference_first", "predicted_class").
    vmin, vmax : float
        Color scale limits.
    figsize : tuple
        Figure size.
    cmap : str
        Colormap.
    save : str, optional
        Path to save figure.
    show : bool
        Whether to call ``plt.show()``.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    cnv_mat = adata.obsm.get(f"X_{cnv_key}")
    if cnv_mat is None:
        raise KeyError(f"'X_{cnv_key}' not found in adata.obsm")

    cnv_mat = np.asarray(cnv_mat)
    n_show = min(n_show, adata.n_obs)

    rng = np.random.default_rng(random_state)
    show_cells = rng.choice(adata.obs_names, size=n_show, replace=False)
    cell_idx = adata.obs_names.get_indexer(show_cells)

    # Sort cells
    show_obs = adata.obs.loc[show_cells].copy()
    if sort_by == "reference_first" and f"{cnv_key}_is_reference" in show_obs.columns:
        show_obs["_is_ref"] = show_obs[f"{cnv_key}_is_reference"].astype(bool)
        show_obs = show_obs.sort_values(["_is_ref", f"{cnv_key}_score"], ascending=[False, True])
    elif sort_by == "predicted_class" and f"{cnv_key}_predicted_class" in show_obs.columns:
        show_obs["_is_aneu"] = show_obs[f"{cnv_key}_predicted_class"].astype(str) == "aneuploid"
        show_obs = show_obs.sort_values(["_is_aneu", f"{cnv_key}_score"], ascending=[False, True])
    else:
        show_obs = show_obs.sort_values(f"{cnv_key}_score", ascending=True)

    show_cells = show_obs.index.values
    cell_idx = adata.obs_names.get_indexer(show_cells)
    Z_show = cnv_mat[cell_idx, :]

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        Z_show, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest"
    )
    ax.set_xlabel("Genes ordered by chromosome & position")
    ax.set_ylabel("Cells (sorted)")
    ax.set_title("CNV Z-score Heatmap")
    cbar = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    cbar.set_label("Z")

    # Chromosome boundaries
    if "chromosome" in adata.var.columns:
        chrom = np.array(adata.var["chromosome"].astype(str))
        boundaries = np.where(chrom[:-1] != chrom[1:])[0] + 0.5
        for b in boundaries:
            ax.axvline(b, color="k", linewidth=0.3, alpha=0.4)

    plt.tight_layout()

    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")

    if show:
        plt.show()

    return fig


def plot_per_chromosome_scores(
    adata: AnnData,
    cnv_key: str = "cnv",
    class_col: Optional[str] = None,
    figsize: Tuple[float, float] = (16, 5),
    save: Optional[str] = None,
    show: bool = True,
) -> Any:
    """Plot per-chromosome CNV scores as split violin plots.

    Parameters
    ----------
    adata : AnnData
        Data with per-chromosome score columns (``chr*_score``).
    cnv_key : str
        Key prefix.
    class_col : str, optional
        Column for grouping (defaults to ``f"{cnv_key}_predicted_class"``).
    figsize : tuple
        Figure size.
    save : str, optional
        Path to save figure.
    show : bool
        Whether to call ``plt.show()``.

    Returns
    -------
    matplotlib.figure.Figure
    """
    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("seaborn is required for plot_per_chromosome_scores")

    import matplotlib.pyplot as plt

    if class_col is None:
        class_col = f"{cnv_key}_predicted_class"

    chr_score_cols = [c for c in adata.obs.columns if c.startswith("chr") and c.endswith("_score")]
    if not chr_score_cols:
        raise ValueError("No per-chromosome score columns found (expected 'chr*_score')")

    df = adata.obs[chr_score_cols + ([class_col] if class_col in adata.obs.columns else [])].copy()

    id_vars = [class_col] if class_col in df.columns else []
    long = df.melt(id_vars=id_vars, var_name="chrom", value_name="cnv_chr_score")
    long["chrom"] = (
        long["chrom"]
        .str.replace("chr", "", regex=False)
        .str.replace("_score", "", regex=False)
    )

    # Standard chromosome ordering
    all_chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
    present = [c for c in all_chroms if c in long["chrom"].unique()]
    long["chrom"] = pd.Categorical(long["chrom"], categories=present, ordered=True)
    long = long.sort_values("chrom")

    fig, ax = plt.subplots(figsize=figsize)

    if class_col in long.columns:
        present_classes = [c for c in ("diploid", "aneuploid") if c in long[class_col].unique()]
        palette = {c: col for c, col in zip(
            ("diploid", "aneuploid"), ("#2ecc71", "#e74c3c")
        ) if c in present_classes}
        sns.violinplot(
            data=long, x="chrom", y="cnv_chr_score", hue=class_col,
            split=True, inner="quartile", palette=palette, ax=ax, cut=0,
            order=present,
        )
    else:
        sns.violinplot(
            data=long, x="chrom", y="cnv_chr_score", ax=ax, cut=0,
            order=present,
        )

    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Median(|Z|) per chromosome")
    ax.set_title("Per-chromosome CNV score by predicted class")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")

    if show:
        plt.show()

    return fig


def plot_aneuploid_proportion(
    adata: AnnData,
    cnv_key: str = "cnv",
    groupby: Optional[Union[str, List[str]]] = None,
    figsize: Tuple[float, float] = (16, 6),
    save: Optional[str] = None,
    show: bool = True,
) -> Any:
    """Plot aneuploid cell proportion by group (sample, cluster, etc.).

    Parameters
    ----------
    adata : AnnData
        Data with predicted class.
    cnv_key : str
        Key prefix.
    groupby : str or list, optional
        Column(s) to group by. Defaults to ``["sampleID", "leiden_clusters"]``
        if available.
    figsize : tuple
        Figure size.
    save : str, optional
        Path to save figure.
    show : bool
        Whether to call ``plt.show()``.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    class_col = f"{cnv_key}_predicted_class"
    if class_col not in adata.obs.columns:
        raise KeyError(f"'{class_col}' not found. Run infer_cnv with predict_aneuploid=True.")

    if groupby is None:
        candidates = [c for c in ("sampleID", "sample", "leiden_clusters", "cell_type") if c in adata.obs.columns]
        groupby = candidates[:2] if len(candidates) >= 2 else candidates[:1]
    elif isinstance(groupby, str):
        groupby = [groupby]

    valid_groups = [g for g in groupby if g in adata.obs.columns]
    if not valid_groups:
        raise ValueError("No valid groupby columns found in adata.obs")

    n_panels = len(valid_groups)
    fig, axes = plt.subplots(1, n_panels, figsize=figsize)
    if n_panels == 1:
        axes = [axes]

    for ax, group_col in zip(axes, valid_groups):
        tab = (
            adata.obs.groupby(group_col)[class_col]
            .value_counts(normalize=False)
            .unstack(fill_value=0)
        )
        tab["total"] = tab.sum(axis=1)
        if "aneuploid" not in tab.columns:
            tab["aneuploid"] = 0
        tab["aneuploid_pct"] = tab["aneuploid"] / tab["total"] * 100
        tab = tab.sort_values("aneuploid_pct", ascending=True)

        y_pos = range(len(tab))
        ax.barh(y_pos, tab["aneuploid_pct"], color="#e74c3c", height=0.6)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(tab.index.astype(str))
        ax.set_xlabel("Aneuploid %")
        ax.set_title(f"Aneuploid proportion by {group_col}")
        max_val = tab["aneuploid_pct"].max()
        ax.set_xlim(0, max_val * 1.15 if max_val > 0 else 100)
        ax.grid(axis="x", alpha=0.3)

        for i, v in enumerate(tab["aneuploid_pct"]):
            ax.text(
                v + (max_val * 0.01), i, f"{v:.1f}%",
                va="center", color="black", fontsize=9,
            )

    plt.tight_layout()

    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")

    if show:
        plt.show()

    return fig
