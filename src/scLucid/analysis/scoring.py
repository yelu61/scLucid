"""
Enhanced scoring module with comprehensive gene set management and visualization.
"""

import logging
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData
from scipy.stats import mannwhitneyu, ttest_ind, zscore

from ..markers.manager import _get_marker_path, _load_marker_file
from ..utils import sanitize_for_hdf5

log = logging.getLogger(__name__)

__all__ = [
    "FunctionalSignatureManager",
    "score_by_gene_sets",
    "calculate_signature_matrix",
    "plot_signature_heatmap",
    "plot_delta_heatmap",
    "batch_plot_delta_heatmap",
    "plot_score_violin_with_stats",
    "batch_compare_scores",
]


# ===================== Gene Set Manager =====================


class FunctionalSignatureManager:
    def __init__(self, species: str = "human", custom_signatures: Optional[str] = None):
        self.species = species.lower()
        self.signatures: Dict[str, List[str]] = {}
        self.categories: Dict[str, List[str]] = {}  # ✅ 新增
        self._load_builtin_signatures()
        if custom_signatures:
            self._load_custom_signatures(custom_signatures)

    def _load_builtin_signatures(self):
        """Load built-in functional signatures from resources."""
        try:
            marker_path = _get_marker_path("functional_signatures")
            data = _load_marker_file(marker_path)

            if self.species in data:
                species_data = data[self.species]

                # ✅ Load categories if available
                if "_categories" in species_data:
                    self.categories = species_data["_categories"]
                    # Remove from signatures dict
                    species_data = {
                        k: v for k, v in species_data.items() if k != "_categories"
                    }

                self.signatures = species_data
                log.info(
                    f"Loaded {len(self.signatures)} built-in signatures for {self.species}"
                )
                if self.categories:
                    log.info(f"Loaded {len(self.categories)} signature categories")
            else:
                raise ValueError(f"Species '{self.species}' not found")
        except FileNotFoundError:
            log.warning("Built-in functional_signatures.json not found.")

    def get_category(self, category_name: str) -> Dict[str, List[str]]:
        """
        Get all signatures in a category.

        Parameters
        ----------
        category_name : str
            Category name (e.g., 'Immune_Function', 'Metabolism')

        Returns
        -------
        dict
            Dictionary of {signature_name: genes} for that category

        Examples
        --------
        >>> manager = FunctionalSignatureManager(species='human')
        >>> immune_sigs = manager.get_category('Immune_Function')
        >>> print(list(immune_sigs.keys()))
        ['Cytotoxicity', 'Exhausted', 'Pro_inflammatory', ...]
        """
        if category_name not in self.categories:
            raise KeyError(
                f"Category '{category_name}' not found. "
                f"Available: {list(self.categories.keys())}"
            )

        sig_names = self.categories[category_name]
        return {
            name: self.signatures[name] for name in sig_names if name in self.signatures
        }

    def list_categories(self) -> List[str]:
        """List all available categories."""
        return list(self.categories.keys())

    def _load_custom_signatures(self, custom_path: str):
        """Load custom signatures from user-provided file."""
        try:
            custom_data = _load_marker_file(Path(custom_path))

            # If custom file has species structure
            if self.species in custom_data:
                custom_data = custom_data[self.species]

            # Merge with existing signatures
            self.signatures.update(custom_data)
            log.info(f"Loaded custom signatures from {custom_path}")
        except Exception as e:
            log.warning(f"Failed to load custom signatures: {e}")

    def get_signature(self, name: str) -> List[str]:
        """Get a specific signature by name."""
        if name not in self.signatures:
            raise KeyError(
                f"Signature '{name}' not found. Available: {list(self.signatures.keys())}"
            )
        return self.signatures[name]

    def get_all_signatures(self) -> Dict[str, List[str]]:
        """Get all available signatures."""
        return self.signatures.copy()

    def add_signature(self, name: str, genes: List[str]):
        """Add a custom signature."""
        self.signatures[name] = genes
        log.info(f"Added custom signature '{name}' with {len(genes)} genes")

    def list_signatures(self) -> List[str]:
        """List all available signature names."""
        return list(self.signatures.keys())

    def convert_to_species(self, target_species: str) -> "FunctionalSignatureManager":
        """
        Convert gene names to another species using gProfiler.
        """
        try:
            from gprofiler import GProfiler
        except ImportError:
            raise ImportError(
                "gprofiler-official required for species conversion. "
                "Install with: pip install gprofiler-official"
            )
        
        org_map = {"human": "hsapiens", "mouse": "mmusculus"}
        
        gp = GProfiler(return_dataframe=True)
        
        # Collect all genes
        all_genes = list({gene for genes in self.signatures.values() for gene in genes})
        
        # Convert
        log.info(f"Converting {len(all_genes)} genes from {self.species} to {target_species}...")
        ortho = gp.orth(
            query=all_genes,
            organism=org_map.get(self.species, self.species),
            target=org_map.get(target_species, target_species),
        )
        
        mapping = (
            ortho[ortho["ortholog.name"] != "n.s."][["incoming.name", "ortholog.name"]]
            .set_index("incoming.name")["ortholog.name"]
            .to_dict()
        )
        
        # Create new manager with converted genes
        new_manager = FunctionalSignatureManager.__new__(FunctionalSignatureManager)
        new_manager.species = target_species
        new_manager.signatures = {}
        new_manager.categories = {}  # ✅ 初始化 categories
        
        # Convert signatures
        for name, genes in self.signatures.items():
            converted = list({mapping[g] for g in genes if g in mapping})
            if converted:
                new_manager.add_signature(name, converted)
            else:
                log.warning(f"No genes converted for signature '{name}'")
        
        # ✅ 尝试保留分类信息（如果签名名称未改变）
        for category, sig_names in self.categories.items():
            converted_sigs = [s for s in sig_names if s in new_manager.signatures]
            if converted_sigs:
                new_manager.categories[category] = converted_sigs
        
        log.info(f"Conversion complete. {len(new_manager.signatures)} signatures available.")
        if new_manager.categories:
            log.info(f"Preserved {len(new_manager.categories)} categories.")
        
        return new_manager


# ===================== Helper Functions =====================


def _ensure_scoring_namespace(adata: AnnData) -> dict:
    """Ensure the scoring namespace exists in adata.uns and return it."""
    return (
        adata.uns.setdefault("sclucid", {})
        .setdefault("analysis", {})
        .setdefault("scoring", {})
    )


def _cohens_d(x: np.ndarray, y: np.ndarray) -> Optional[float]:
    """Compute Cohen's d for two independent samples (unequal n)."""
    try:
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        nx, ny = len(x), len(y)
        if nx < 2 or ny < 2:
            return None
        sx, sy = x.std(ddof=1), y.std(ddof=1)
        if sx == 0 and sy == 0:
            return 0.0
        sp = np.sqrt(((nx - 1) * sx**2 + (ny - 1) * sy**2) / (nx + ny - 2))
        if sp == 0:
            return None
        return (x.mean() - y.mean()) / sp
    except Exception:
        return None


def _validate_score_column(
    adata: AnnData,
    score_key: str,
    expected_source: str = "score_genes"
) -> bool:
    """
    Validate if a score column exists and was generated correctly.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    score_key : str
        Score column name
    expected_source : str
        Expected generation method
    
    Returns
    -------
    bool
        True if valid, False otherwise
    """
    if score_key not in adata.obs.columns:
        return False
    
    # Check if it's numeric
    if not pd.api.types.is_numeric_dtype(adata.obs[score_key]):
        log.warning(f"Score '{score_key}' is not numeric")
        return False
    
    # Check for suspicious values
    score_data = adata.obs[score_key].dropna()
    if len(score_data) == 0:
        log.warning(f"Score '{score_key}' has no valid values")
        return False
    
    # Check if values are reasonable (scores from score_genes are typically -5 to 5)
    if expected_source == "score_genes":
        if score_data.min() < -10 or score_data.max() > 10:
            log.warning(
                f"Score '{score_key}' has unusual range [{score_data.min():.2f}, {score_data.max():.2f}]. "
                f"Expected range for score_genes is typically [-5, 5]."
            )
    
    return True


# ===================== Core Scoring Functions =====================


def score_by_gene_sets(
    adata: AnnData,
    gene_sets: Union[Dict[str, List[str]], FunctionalSignatureManager],
    layer: Optional[str] = "log1p_norm",
    use_raw: bool = False,
    ctrl_size: int = 50,
    score_name_suffix: str = "_score",
    preserve_missing: bool = True,
    min_genes_required: int = 1,
    **kwargs,
) -> AnnData:
    """
    Score cells for multiple gene sets and add columns to adata.obs.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    gene_sets : dict or FunctionalSignatureManager
        Dictionary of {set_name: [genes]} or a signature manager
    layer : str, optional
        Expression layer when use_raw=False
    use_raw : bool
        If True, score on adata.raw
    ctrl_size : int
        Control gene set size for scanpy.tl.score_genes
    score_name_suffix : str
        Suffix for new obs columns
    preserve_missing : bool
        If True, allow partial matches; if False, skip sets without all genes
    min_genes_required : int
        Minimum number of found genes to attempt scoring
    **kwargs
        Additional arguments passed to scanpy.tl.score_genes

    Returns
    -------
    AnnData
        Modified AnnData with scores added to .obs
    """
    if ctrl_size < 1:
        raise ValueError(f"ctrl_size must be >= 1, got {ctrl_size}")
    if min_genes_required < 1:
        raise ValueError(f"min_genes_required must be >= 1, got {min_genes_required}")
    if not score_name_suffix:
        log.warning("score_name_suffix is empty, scores will overwrite each other!")

    ns = _ensure_scoring_namespace(adata)

    # Handle FunctionalSignatureManager input
    if isinstance(gene_sets, FunctionalSignatureManager):
        gene_sets = gene_sets.get_all_signatures()

    if use_raw:
        if adata.raw is None:
            raise ValueError("adata.raw is not set, but use_raw=True.")
        source_adata = adata.raw
        target_layer = None
    else:
        source_adata = adata
        target_layer = layer
        if target_layer is None:
            log.warning("layer=None and use_raw=False; using adata.X for scoring.")
        elif target_layer not in adata.layers:
            raise ValueError(f"Layer '{target_layer}' not found in adata.layers.")

    total_sets = len(gene_sets)
    scored_count = 0
    skipped_sets: List[str] = []
    per_set_stats: Dict[str, Dict[str, int]] = {}

    log.info(
        f"Scoring cells for {total_sets} gene sets (use_raw={use_raw}, layer={target_layer})..."
    )

    for set_name, genes in gene_sets.items():
        genes = [g for g in genes if isinstance(g, str) and len(g) > 0]
        if len(genes) == 0:
            log.warning(f"Gene set '{set_name}' is empty. Skipping.")
            skipped_sets.append(set_name)
            per_set_stats[set_name] = {"n_input": 0, "n_found": 0, "scored": 0}
            continue

        found = [g for g in genes if g in source_adata.var_names]
        n_input, n_found = len(genes), len(found)
        per_set_stats[set_name] = {"n_input": n_input, "n_found": n_found, "scored": 0}

        if n_found < min_genes_required:
            log.warning(
                f"Gene set '{set_name}': found {n_found}/{n_input} (<{min_genes_required}). Skipping."
            )
            skipped_sets.append(set_name)
            continue

        if not preserve_missing and n_found < n_input:
            log.warning(
                f"Gene set '{set_name}': partial match {n_found}/{n_input} and preserve_missing=False. Skipping."
            )
            skipped_sets.append(set_name)
            continue

        score_name = f"{set_name}{score_name_suffix}"
        try:
            sc.tl.score_genes(
                adata,
                found,
                score_name=score_name,
                use_raw=use_raw,
                layer=target_layer,
                ctrl_size=ctrl_size,
                **kwargs,
            )
            per_set_stats[set_name]["scored"] = 1
            scored_count += 1
        except Exception as e:
            log.warning(f"Failed to score set '{set_name}': {e}")
            skipped_sets.append(set_name)

    ns["gene_set_scoring"] = sanitize_for_hdf5(
        {
            "n_sets_input": total_sets,
            "n_sets_scored": scored_count,
            "n_sets_skipped": len(skipped_sets),
            "skipped_sets": skipped_sets,
            "per_set_stats": per_set_stats,
            "params": {
                "use_raw": use_raw,
                "layer": layer,
                "ctrl_size": ctrl_size,
                "score_name_suffix": score_name_suffix,
                "preserve_missing": preserve_missing,
                "min_genes_required": min_genes_required,
                "scanpy_version": getattr(sc, "__version__", "unknown"),
            },
        }
    )
    log.info(
        f"Completed scoring: {scored_count}/{total_sets} sets scored, {len(skipped_sets)} skipped."
    )
    return adata


# ===================== Matrix-based Analysis =====================


def calculate_signature_matrix(
    adata: AnnData,
    gene_sets: Union[Dict[str, List[str]], FunctionalSignatureManager],
    groupby: str,
    subset_cells: Optional[List[str]] = None,
    use_raw: bool = True,
    ctrl_size: int = 50,
    z_score: bool = True,
) -> pd.DataFrame:
    """
    Calculate mean signature scores per group with optional z-score normalization.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    gene_sets : dict or FunctionalSignatureManager
        Gene sets to score
    groupby : str
        Column in adata.obs to group by
    subset_cells : list, optional
        Subset of cell types/groups to include
    use_raw : bool
        Use adata.raw for scoring
    ctrl_size : int
        Control size for scoring
    z_score : bool
        Apply z-score normalization row-wise

    Returns
    -------
    pd.DataFrame
        Matrix of signature scores (signatures × groups)
    """
    # Handle FunctionalSignatureManager input
    if isinstance(gene_sets, FunctionalSignatureManager):
        gene_sets = gene_sets.get_all_signatures()
    source_var_names = (
        adata.raw.var_names if (use_raw and adata.raw is not None) else adata.var_names
    )

    # Filter valid gene sets
    valid_sets = {
        k: [g for g in v if g in source_var_names] for k, v in gene_sets.items()
    }
    valid_sets = {k: v for k, v in valid_sets.items() if v}

    if not valid_sets:
        raise ValueError("None of the signature genes found in adata.var_names")

    # Subset cells if requested
    if subset_cells is not None:
        adata_ = adata[adata.obs[groupby].isin(subset_cells)].copy()
    else:
        adata_ = adata.copy()

    # Calculate scores
    for sig, genes in valid_sets.items():
        sc.tl.score_genes(
            adata_,
            gene_list=genes,
            score_name=sig,
            ctrl_size=min(len(genes), ctrl_size),
            use_raw=use_raw,
        )

    # Aggregate by group
    df = adata_.obs.groupby(groupby)[list(valid_sets.keys())].mean().T

    # Apply z-score normalization if requested
    if z_score:
        df = df.apply(zscore, axis=1)

    return df


# ===================== Visualization Functions =====================


def plot_signature_heatmap(
    matrix: pd.DataFrame,
    cmap: str = "RdBu_r",
    center: float = 0,
    title: str = "Functional Signature Heatmap",
    figsize: Optional[Tuple[float, float]] = None,
    annot: bool = True,
    fmt: str = ".2f",
    save: Optional[str] = None,
    **kwargs,
) -> plt.Axes:
    """
    Plot a heatmap of signature scores.

    Parameters
    ----------
    matrix : pd.DataFrame
        Signature matrix (signatures × groups)
    cmap : str
        Colormap name
    center : float
        Value to center the colormap at
    title : str
        Plot title
    figsize : tuple, optional
        Figure size (auto-calculated if None)
    annot : bool
        Show values in cells
    fmt : str
        Format string for annotations
    save : str, optional
        Path to save figure
    **kwargs
        Additional arguments for sns.heatmap

    Returns
    -------
    plt.Axes
        Matplotlib axes object
    """
    if figsize is None:
        figsize = (max(8, matrix.shape[1] * 0.9), max(4, matrix.shape[0] * 0.6))

    plt.figure(figsize=figsize)
    ax = sns.heatmap(
        matrix,
        cmap=cmap,
        center=center,
        linewidths=1.5,
        linecolor="white",
        annot=annot,
        fmt=fmt,
        cbar_kws={"label": "Z-score"},
        **kwargs,
    )
    ax.set_title(title, fontsize=16, pad=20)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    plt.tight_layout()

    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")
        log.info(f"Saved heatmap to {save}")

    plt.show()
    return ax


def plot_delta_heatmap(
    adata: AnnData,
    gene_sets: Union[Dict[str, List[str]], FunctionalSignatureManager],
    groupby: str,
    compare_group: str,
    ref_group: str,
    target_group: str,
    cmap: str = "RdBu_r",
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (8, 6),
    use_raw: bool = True,
    ctrl_size: int = 50,
    save: Optional[str] = None,
    **kwargs,
) -> Tuple[plt.Axes, pd.DataFrame]:
    """
    Plot a delta heatmap showing score differences between two groups.

    This is particularly useful for comparing functional remodeling across
    different conditions (e.g., treatment vs control).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    gene_sets : dict or FunctionalSignatureManager
        Gene sets to score
    groupby : str
        Column for cell subsets (e.g., 'celltype')
    compare_group : str
        Column for conditions to compare (e.g., 'treatment')
    ref_group : str
        Reference group name (e.g., 'Control')
    target_group : str
        Target group name (e.g., 'Treated')
    cmap : str
        Colormap (diverging recommended)
    title : str, optional
        Plot title (auto-generated if None)
    figsize : tuple
        Figure size
    use_raw : bool
        Use adata.raw for scoring
    trl_size : int
        Control gene set size for scanpy.tl.score_genes (default: 50)
    save : str, optional
        Path to save figure
    **kwargs
        Additional arguments for sns.heatmap

    Returns
    -------
    tuple
        (matplotlib axes, delta DataFrame)
    """
    # Handle FunctionalSignatureManager input
    if isinstance(gene_sets, FunctionalSignatureManager):
        gene_sets = gene_sets.get_all_signatures()

    source_var_names = (
        adata.raw.var_names if (use_raw and adata.raw is not None) else adata.var_names
    )

    # Calculate scores if not already present
    valid_sigs = []
    for sig, genes in gene_sets.items():
        genes_in = [g for g in genes if g in source_var_names]
        if len(genes_in) > 0:
            if sig not in adata.obs.columns:
                sc.tl.score_genes(
                    adata,
                    gene_list=genes_in,
                    score_name=sig,
                    use_raw=use_raw,
                    ctrl_size=min(len(genes_in), ctrl_size),
                )
            else:
                if not _validate_score_column(adata, sig, "score_genes"):
                    log.warning(
                        f"Existing column '{sig}' may not be a valid score. "
                        f"Recalculating..."
                    )
                    sc.tl.score_genes(
                    adata,
                    gene_list=genes_in,
                    score_name=sig,
                    use_raw=use_raw,
                    ctrl_size=min(len(genes_in), ctrl_size),
                )
            valid_sigs.append(sig)

    if not valid_sigs:
        raise ValueError("No valid signatures found.")

    # Extract data
    obs_data = adata.obs[[groupby, compare_group] + valid_sigs].copy()

    # Calculate means for each group
    mean_ref = (
        obs_data[obs_data[compare_group] == ref_group]
        .groupby(groupby)[valid_sigs]
        .mean()
    )
    mean_target = (
        obs_data[obs_data[compare_group] == target_group]
        .groupby(groupby)[valid_sigs]
        .mean()
    )

    # Align subclusters
    common_subclusters = mean_ref.index.intersection(mean_target.index)
    mean_ref = mean_ref.loc[common_subclusters]
    mean_target = mean_target.loc[common_subclusters]

    # Calculate delta
    delta_df = mean_target - mean_ref

    # Plot
    plt.figure(figsize=figsize)
    max_val = np.max(np.abs(delta_df.values))

    ax = sns.heatmap(
        delta_df,
        cmap=cmap,
        center=0,
        vmin=-max_val,
        vmax=max_val,
        annot=True,
        fmt=".2f",
        linewidths=1,
        linecolor="white",
        cbar_kws={"label": f"Δ Score ({target_group} - {ref_group})"},
        square=True,
        **kwargs,
    )

    if title is None:
        title = f"Functional Remodeling: {target_group} vs {ref_group}"

    ax.set_title(title, fontsize=14, pad=20)
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")
        log.info(f"Saved delta heatmap to {save}")

    plt.show()
    return ax, delta_df


def batch_plot_delta_heatmap(
    adata: AnnData,
    gene_sets: Union[Dict[str, List[str]], FunctionalSignatureManager],
    groupby: str,
    compare_group: str,
    condition_pairs: List[Tuple[str, str]],
    ncols: int = 2,
    figsize_per_plot: Tuple[float, float] = (6, 5),
    use_raw: bool = True,
    ctrl_size: int = 50,
    save_prefix: Optional[str] = None,
    **kwargs,
) -> List[Tuple[pd.DataFrame, Dict]]:
    """
    Plot multiple delta heatmaps for different condition pairs.

    Parameters
    ----------
    condition_pairs : list of tuples
        List of (ref_group, target_group) pairs to compare
    ncols : int
        Number of columns in the subplot grid
    figsize_per_plot : tuple
        Size of each individual subplot
    save_prefix : str, optional
        Prefix for saved figures

    Returns
    -------
    list of tuples
        List of (delta_df, stats_dict) for each comparison

    Examples
    --------
    >>> results = batch_plot_delta_heatmap(
    ...     adata,
    ...     sig_manager,
    ...     groupby='celltype',
    ...     compare_group='treatment',
    ...     condition_pairs=[('Control', 'Drug_A'), ('Control', 'Drug_B')],
    ...     save_prefix='delta_analysis'
    ... )
    """
    # Handle FunctionalSignatureManager input
    if isinstance(gene_sets, FunctionalSignatureManager):
        gene_sets = gene_sets.get_all_signatures()
    
    n_pairs = len(condition_pairs)
    nrows = (n_pairs + ncols - 1) // ncols
    
    # Calculate figure size
    total_width = figsize_per_plot[0] * ncols
    total_height = figsize_per_plot[1] * nrows
    
    fig, axes = plt.subplots(
        nrows, ncols, 
        figsize=(total_width, total_height),
        squeeze=False  # ✅ 确保返回 2D 数组
    )
    
    results = []
    
    # Pre-calculate all delta matrices (避免重复计算)
    source_var_names = (
        adata.raw.var_names if (use_raw and adata.raw is not None) else adata.var_names
    )
    
    # Calculate scores once for all comparisons
    valid_sigs = []
    for sig, genes in gene_sets.items():
        genes_in = [g for g in genes if g in source_var_names]
        if len(genes_in) > 0:
            if sig not in adata.obs.columns:
                sc.tl.score_genes(
                    adata,
                    gene_list=genes_in,
                    score_name=sig,
                    use_raw=use_raw,
                    ctrl_size=min(len(genes_in), ctrl_size),
                )
            valid_sigs.append(sig)
    
    if not valid_sigs:
        raise ValueError("No valid signatures found.")
    
    obs_data = adata.obs[[groupby, compare_group] + valid_sigs].copy()
    
    # Plot each comparison
    for idx, (ref, target) in enumerate(condition_pairs):
        row = idx // ncols
        col = idx % ncols
        ax = axes[row, col]
        
        # Calculate delta for this pair
        mean_ref = (
            obs_data[obs_data[compare_group] == ref]
            .groupby(groupby)[valid_sigs]
            .mean()
        )
        mean_target = (
            obs_data[obs_data[compare_group] == target]
            .groupby(groupby)[valid_sigs]
            .mean()
        )
        
        common_groups = mean_ref.index.intersection(mean_target.index)
        if len(common_groups) == 0:
            log.warning(f"No common groups for {ref} vs {target}")
            ax.text(0.5, 0.5, f"No data\n{ref} vs {target}", 
                   ha='center', va='center', transform=ax.transAxes)
            ax.axis('off')
            continue
        
        mean_ref = mean_ref.loc[common_groups]
        mean_target = mean_target.loc[common_groups]
        delta_df = mean_target - mean_ref
        
        # Plot on subplot
        max_val = np.max(np.abs(delta_df.values))
        
        sns.heatmap(
            delta_df,
            cmap=kwargs.get('cmap', 'RdBu_r'),
            center=0,
            vmin=-max_val,
            vmax=max_val,
            annot=True,
            fmt=".2f",
            linewidths=1,
            linecolor='white',
            cbar_kws={'label': f'Δ Score'},
            square=True,
            ax=ax,
        )
        
        ax.set_title(f"{target} vs {ref}", fontsize=12, pad=10)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=9)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=9)
        
        # Save individual if requested
        if save_prefix:
            fig_individual = plt.figure(figsize=figsize_per_plot)
            ax_individual = fig_individual.add_subplot(111)
            sns.heatmap(
                delta_df,
                cmap=kwargs.get('cmap', 'RdBu_r'),
                center=0,
                vmin=-max_val,
                vmax=max_val,
                annot=True,
                fmt=".2f",
                linewidths=1,
                linecolor='white',
                cbar_kws={'label': f'Δ Score ({target} - {ref})'},
                square=True,
                ax=ax_individual,
            )
            ax_individual.set_title(f"Functional Remodeling: {target} vs {ref}")
            plt.tight_layout()
            plt.savefig(f"{save_prefix}_{ref}_vs_{target}.pdf", dpi=300, bbox_inches='tight')
            plt.close(fig_individual)
        
        results.append((delta_df, {'ref': ref, 'target': target, 'n_groups': len(common_groups)}))
    
    # Hide unused subplots
    for idx in range(n_pairs, nrows * ncols):
        row = idx // ncols
        col = idx % ncols
        axes[row, col].axis('off')
    
    plt.tight_layout()
    
    if save_prefix:
        plt.savefig(f"{save_prefix}_combined.pdf", dpi=300, bbox_inches='tight')
        log.info(f"Saved combined figure to {save_prefix}_combined.pdf")
    
    plt.show()
    return results


def plot_score_violin_with_stats(
    adata: AnnData,
    score_key: str,
    groupby: str,
    group1: str,
    group2: str,
    subset_key: Optional[str] = None,
    subset_value: Optional[str] = None,
    test: Literal["wilcoxon", "ttest"] = "wilcoxon",
    palette: Optional[List[str]] = None,
    figsize: Tuple[float, float] = (6, 7),
    show_points: bool = True,
    sig_thresholds: Tuple[float, float, float] = (0.05, 0.01, 0.001),
    save: Optional[str] = None,
) -> Tuple[plt.Axes, Dict[str, float]]:
    """
    Plot violin plot with statistical test results.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    score_key : str
        Score column in adata.obs
    groupby : str
        Column to group by
    group1 : str
        First group name
    group2 : str
        Second group name
    subset_key : str, optional
        Additional column to subset (e.g., 'celltype')
    subset_value : str, optional
        Value to subset on
    test : str
        Statistical test ('wilcoxon' or 'ttest')
    palette : list, optional
        Color palette for groups
    figsize : tuple
        Figure size
    show_points : bool
        Overlay individual points
    sig_thresholds : tuple of floats
        Thresholds for *, **, *** significance markers (default: 0.05, 0.01, 0.001)
    save : str, optional
        Path to save figure

    Returns
    -------
    tuple
        (matplotlib axes, statistics dict)
    """
    # Subset data if requested
    plot_df = adata.obs.copy()
    if subset_key and subset_value:
        plot_df = plot_df[plot_df[subset_key] == subset_value]
        title_suffix = f" in {subset_value}"
    else:
        title_suffix = ""

    if plot_df.empty:
        raise ValueError("No data available after subsetting")

    # Extract group data
    group1_data = plot_df[plot_df[groupby] == group1][score_key].dropna()
    group2_data = plot_df[plot_df[groupby] == group2][score_key].dropna()

    # Perform statistical test
    if test == "wilcoxon":
        stat, pval = mannwhitneyu(group1_data, group2_data, alternative="two-sided")
        test_name = "Mann-Whitney U"
    elif test == "ttest":
        stat, pval = ttest_ind(group1_data, group2_data, equal_var=False)
        test_name = "Welch's t-test"
    else:
        raise ValueError("test must be 'wilcoxon' or 'ttest'")

    # Calculate effect size
    cohens_d = _cohens_d(group1_data.values, group2_data.values)

    # Store statistics
    stats = {
        "statistic": float(stat),
        "pvalue": float(pval),
        "n_group1": len(group1_data),
        "n_group2": len(group2_data),
        "mean_group1": float(group1_data.mean()),
        "mean_group2": float(group2_data.mean()),
        "cohens_d": float(cohens_d) if cohens_d is not None else None,
        "test": test_name,
    }

    # Print statistics
    log.info("\n=== Statistical Test Results ===")
    log.info(f"Test: {test_name}")
    log.info(
        f"Group 1 ({group1}): n={stats['n_group1']}, mean={stats['mean_group1']:.4f}"
    )
    log.info(
        f"Group 2 ({group2}): n={stats['n_group2']}, mean={stats['mean_group2']:.4f}"
    )
    log.info(f"Statistic: {stats['statistic']:.4f}")
    log.info(f"P-value: {stats['pvalue']:.6f}")
    if cohens_d is not None:
        log.info(f"Cohen's d: {stats['cohens_d']:.4f}")
    log.info("=" * 32)

    # Plot
    if palette is None:
        palette = ["#40CCB2", "#E04110"]

    fig, ax = plt.subplots(figsize=figsize)

    # Violin plot
    sns.violinplot(
        data=plot_df,
        x=groupby,
        y=score_key,
        order=[group1, group2],
        palette=palette,
        ax=ax,
    )

    # Overlay points if requested
    if show_points:
        sns.stripplot(
            data=plot_df,
            x=groupby,
            y=score_key,
            order=[group1, group2],
            color="black",
            jitter=0.1,
            size=2.5,
            alpha=0.5,
            ax=ax,
        )

    # Add statistical annotation
    y_max = plot_df[score_key].max()
    y_min = plot_df[score_key].min()
    y_range = y_max - y_min

    # Draw significance line
    line_height = y_max + 0.05 * y_range
    ax.plot([0, 1], [line_height, line_height], "k-", linewidth=1.5)

    # Add p-value text
    if pval < sig_thresholds[2]:
        sig_text = "***"
    elif pval < sig_thresholds[1]:
        sig_text = "**"
    elif pval < sig_thresholds[0]:
        sig_text = "*"
    else:
        sig_text = "ns"

    ax.text(
        0.5,
        line_height + 0.02 * y_range,
        f"{sig_text}\np={pval:.4f}",
        ha="center",
        va="bottom",
        fontsize=10,
    )

    # Set labels and title
    ax.set_title(f"{score_key}{title_suffix}", fontsize=16)
    ax.set_xlabel("Group", fontsize=12)
    ax.set_ylabel("Score", fontsize=12)

    plt.tight_layout()

    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")
        log.info(f"Saved violin plot to {save}")

    plt.show()
    return ax, stats


# ===================== Optional Functions =====================


def _calculate_group_stats(
    adata: AnnData,
    score_key: str,
    groupby: str,
    group1: str,
    group2: str,
    method: Literal["ttest", "wilcoxon"] = "wilcoxon",
) -> pd.DataFrame:
    """
    Calculate statistics comparing two groups for a given score.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    score_key : str
        Score column in adata.obs
    groupby : str
        Column to group by
    group1 : str
        First group name (or "rest" for one-vs-rest)
    group2 : str
        Second group name (or "rest" for one-vs-rest)
    method : str
        Statistical test method ('ttest' or 'wilcoxon')

    Returns
    -------
    pd.DataFrame
        Single-row DataFrame with comparison statistics

    Examples
    --------
    >>> stats = _calculate_group_stats(
    ...     adata, 'Cytotoxicity_score', 'tissue', 'Tumor', 'Normal', 'wilcoxon'
    ... )
    >>> print(stats[['group1', 'group2', 'pvalue', 'cohens_d']])
    """
    # Validate input
    if score_key not in adata.obs.columns:
        log.warning(f"Score '{score_key}' not found in adata.obs")
        return pd.DataFrame()

    if groupby not in adata.obs.columns:
        log.warning(f"Groupby column '{groupby}' not found in adata.obs")
        return pd.DataFrame()

    # Handle "rest" comparison (one-vs-rest)
    if group2 == "rest":
        mask1 = adata.obs[groupby] == group1
        mask2 = ~mask1
        g2_label = "rest"
    else:
        # Validate both groups exist
        available_groups = adata.obs[groupby].unique()
        if group1 not in available_groups:
            log.warning(f"Group '{group1}' not found in {groupby}")
            return pd.DataFrame()
        if group2 not in available_groups:
            log.warning(f"Group '{group2}' not found in {groupby}")
            return pd.DataFrame()

        mask1 = adata.obs[groupby] == group1
        mask2 = adata.obs[groupby] == group2
        g2_label = group2

    # Extract data
    data1 = adata.obs.loc[mask1, score_key].dropna()
    data2 = adata.obs.loc[mask2, score_key].dropna()

    # Check sample sizes
    if len(data1) < 2 or len(data2) < 2:
        log.warning(
            f"Insufficient data for '{score_key}': "
            f"{group1} (n={len(data1)}) vs {g2_label} (n={len(data2)}). "
            f"Minimum 2 samples required per group."
        )
        return pd.DataFrame()

    # Perform statistical test
    try:
        if method == "wilcoxon":
            stat, pval = mannwhitneyu(data1, data2, alternative="two-sided")
            test_name = "Mann-Whitney U"
        elif method == "ttest":
            stat, pval = ttest_ind(data1, data2, equal_var=False)
            test_name = "Welch's t-test"
        else:
            raise ValueError(f"Unknown method: {method}. Use 'ttest' or 'wilcoxon'.")
    except Exception as e:
        log.warning(
            f"Statistical test failed for '{score_key}' ({group1} vs {g2_label}): {e}"
        )
        return pd.DataFrame()

    # Calculate effect size
    effect_size = _cohens_d(data1.values, data2.values)

    # Calculate fold change (for log-transformed data, this is additive difference)
    mean_diff = float(data1.mean() - data2.mean())

    # Return results as single-row DataFrame
    return pd.DataFrame(
        [
            {
                "score": score_key,
                "group1": group1,
                "group2": g2_label,
                "n_group1": int(len(data1)),
                "n_group2": int(len(data2)),
                "mean_group1": float(data1.mean()),
                "mean_group2": float(data2.mean()),
                "std_group1": float(data1.std()),
                "std_group2": float(data2.std()),
                "mean_diff": mean_diff,
                "statistic": float(stat),
                "pvalue": float(pval),
                "cohens_d": float(effect_size) if effect_size is not None else np.nan,
                "method": test_name,
            }
        ]
    )


def batch_compare_scores(
    adata: AnnData,
    score_keys: List[str],
    groupby: str,
    group_pairs: Optional[List[Tuple[str, str]]] = None,
    method: Literal["ttest", "wilcoxon"] = "wilcoxon",
    adjust_pvalues: bool = True,
) -> pd.DataFrame:
    """
    Batch differential comparison for multiple scores and group pairs.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    score_keys : list of str
        List of score column names in adata.obs
    groupby : str
        Column name to group by
    group_pairs : list of tuples, optional
        List of (group1, group2) pairs to compare.
        If None, performs one-vs-rest for all groups.
        Use "rest" as group2 for explicit one-vs-rest.
    method : str
        Statistical test method ('ttest' or 'wilcoxon')
    adjust_pvalues : bool
        If True, apply FDR correction (requires statsmodels).
        If False, return uncorrected p-values.

    Returns
    -------
    pd.DataFrame
        Combined results with columns: score, group1, group2, pvalue, cohens_d, etc.

    Examples
    --------
    >>> # Compare all signatures between two conditions
    >>> results = batch_compare_scores(
    ...     adata,
    ...     score_keys=['Cytotoxicity_score', 'Exhausted_score', 'M1_score'],
    ...     groupby='tissue',
    ...     group_pairs=[('Tumor', 'Normal')],
    ...     method='wilcoxon'
    ... )
    >>>
    >>> # Filter significant results
    >>> sig = results[(results['pvalue'] < 0.01) & (results['cohens_d'].abs() > 0.5)]
    >>> print(sig[['score', 'pvalue', 'cohens_d']])
    """
    ns = _ensure_scoring_namespace(adata)

    # Validate score keys
    available_scores = [s for s in score_keys if s in adata.obs.columns]
    missing_scores = [s for s in score_keys if s not in adata.obs.columns]

    if missing_scores:
        log.warning(f"Missing score columns (skipped): {missing_scores}")

    if len(available_scores) == 0:
        log.warning("No valid score columns found. Returning empty DataFrame.")
        return pd.DataFrame()

    # Generate group pairs if not provided
    if group_pairs is None:
        groups = adata.obs[groupby].astype(str).unique()
        group_pairs = [(g, "rest") for g in groups]
        log.info(f"Auto-generated {len(group_pairs)} one-vs-rest comparisons")

    # Batch compare
    results = []
    total_comparisons = len(available_scores) * len(group_pairs)

    log.info(
        f"Running {total_comparisons} comparisons ({len(available_scores)} scores × {len(group_pairs)} pairs)..."
    )

    for score in available_scores:
        for g1, g2 in group_pairs:
            df = _calculate_group_stats(adata, score, groupby, g1, g2, method)
            if not df.empty:
                results.append(df)

    # Combine results
    if results:
        all_results = pd.concat(results, ignore_index=True)

        if adjust_pvalues:
            try:
                from statsmodels.stats.multitest import multipletests
                
                _, pvals_corrected, _, _ = multipletests(
                    all_results["pvalue"],
                    method="fdr_bh",
                )
                all_results["pvalue_adjusted"] = pvals_corrected
                all_results = all_results.sort_values("pvalue_adjusted")
                log.info("✅ Applied FDR correction (Benjamini-Hochberg)")
                
            except ImportError:
                log.warning(
                    "⚠️ statsmodels not installed. Skipping FDR correction. "
                    "Install with: pip install statsmodels"
                )
                all_results = all_results.sort_values("pvalue")
        else:
            all_results = all_results.sort_values("pvalue")
        
        # Save to namespace
        ns["batch_compare_results"] = sanitize_for_hdf5(
            {
                "data": all_results.to_dict(orient="records"),
                "n_comparisons": len(all_results),
                "method": method,
                "fdr_corrected": adjust_pvalues and "pvalue_adjusted" in all_results.columns,
            }
        )
        
        log.info(f"✅ Completed {len(all_results)} comparisons")
        return all_results
    else:
        log.warning("No valid comparisons produced.")
        return pd.DataFrame()
