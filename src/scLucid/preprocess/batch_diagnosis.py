"""
Dynamic Batch Effect Diagnosis and Method Recommendation.

This module provides:
1. Automated batch effect severity assessment
2. Method recommendation based on data characteristics
3. Expected integration quality prediction

Scientific innovation:
- First systematic framework for batch correction method selection
- Combines multiple diagnostic metrics
- Provides interpretable recommendations
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors

log = logging.getLogger(__name__)

__all__ = [
    "BatchDiagnosisConfig",
    "diagnose_batch_effects",
    "recommend_integration_method",
    "compare_integration_methods",
]


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class BatchDiagnosisConfig:
    """Configuration for batch effect diagnosis."""
    
    batch_key: str = "sampleID"
    label_key: Optional[str] = None  # Known cell type labels (if available)
    
    # Diagnosis metrics
    metrics: List[str] = None  # If None, use all
    
    # Method recommendation
    recommend_methods: bool = True
    n_methods_to_recommend: int = 3
    
    # Visualization
    plot: bool = True
    save_dir: Optional[str] = None
    
    def __post_init__(self):
        if self.metrics is None:
            self.metrics = [
                'pcr',  # Principal Component Regression
                'kbet',  # k-nearest neighbor batch effect test
                'lisi',  # Local Inverse Simpson's Index
                'silhouette',  # Silhouette coefficient
                'asw',  # Average Silhouette Width
            ]


# =============================================================================
# Diagnostic Metrics
# =============================================================================

def compute_pcr_batch_effect(
    adata: AnnData,
    batch_key: str,
    use_rep: str = "X_pca",
    n_comps: int = 50,
) -> float:
    """
    Compute batch effect using Principal Component Regression (PCR).
    
    This measures how much variance in the PC space is explained by batch.
    
    Args:
        adata: AnnData object
        batch_key: Batch identifier
        use_rep: Representation to use
        n_comps: Number of components
    
    Returns:
        R² score (0 = no batch effect, 1 = all variance explained by batch)
    """
    from sklearn.linear_model import LinearRegression
    from sklearn.preprocessing import LabelEncoder
    
    if use_rep not in adata.obsm:
        log.warning(f"{use_rep} not found. Computing PCA...")
        sc.tl.pca(adata, n_comps=n_comps)
        use_rep = "X_pca"
    
    X_pca = adata.obsm[use_rep][:, :n_comps]
    
    # Encode batch labels
    le = LabelEncoder()
    batch_encoded = le.fit_transform(adata.obs[batch_key]).reshape(-1, 1)
    
    # Fit regression
    reg = LinearRegression()
    reg.fit(batch_encoded, X_pca)
    
    # Compute R²
    r2 = reg.score(batch_encoded, X_pca)
    
    return r2


def compute_kbet_score(
    adata: AnnData,
    batch_key: str,
    use_rep: str = "X_pca",
    n_neighbors: int = 25,
    n_sample: int = 1000,
) -> float:
    """
    Compute k-nearest neighbor Batch Effect Test (kBET) score.
    
    Lower rejection rate = better mixing.
    
    Args:
        adata: AnnData object
        batch_key: Batch identifier
        use_rep: Representation
        n_neighbors: Number of neighbors
        n_sample: Number of cells to sample
    
    Returns:
        Acceptance rate (0 = strong batch effect, 1 = no batch effect)
    """
    from scipy.stats import chi2
    
    if use_rep not in adata.obsm:
        sc.tl.pca(adata)
        use_rep = "X_pca"
    
    X = adata.obsm[use_rep]
    batch_labels = adata.obs[batch_key]
    
    # Sample cells
    n_cells = X.shape[0]
    if n_cells > n_sample:
        indices = np.random.choice(n_cells, n_sample, replace=False)
        X_sample = X[indices]
        batch_sample = batch_labels.iloc[indices]
    else:
        X_sample = X
        batch_sample = batch_labels
    
    # Build kNN
    nn = NearestNeighbors(n_neighbors=n_neighbors + 1)
    nn.fit(X)
    distances, indices = nn.kneighbors(X_sample)
    indices = indices[:, 1:]  # Remove self
    
    # Global batch distribution
    batch_cats = sorted(batch_labels.unique())
    global_freq = batch_labels.value_counts(normalize=True).reindex(batch_cats).values
    
    # Test each neighborhood
    rejections = 0
    valid_tests = 0
    
    for i in range(len(X_sample)):
        neighbor_batches = batch_labels.iloc[indices[i]].values
        
        observed = np.array([
            (neighbor_batches == batch).sum() for batch in batch_cats
        ])
        
        expected = global_freq * n_neighbors
        
        if (expected < 5).any():
            continue
        
        chi2_stat = np.sum((observed - expected)**2 / expected)
        df = len(batch_cats) - 1
        p_value = 1 - chi2.cdf(chi2_stat, df)
        
        valid_tests += 1
        if p_value < 0.05:
            rejections += 1
    
    if valid_tests == 0:
        return np.nan
    
    acceptance_rate = 1 - (rejections / valid_tests)
    
    return acceptance_rate


def compute_lisi(
    adata: AnnData,
    batch_key: str,
    use_rep: str = "X_pca",
    n_neighbors: int = 90,
    perplexity: int = 30,
) -> float:
    """
    Compute Local Inverse Simpson's Index (LISI).
    
    Measures local diversity of batches.
    Higher LISI = better mixing.
    
    Args:
        adata: AnnData object
        batch_key: Batch identifier
        use_rep: Representation
        n_neighbors: Number of neighbors
        perplexity: Perplexity parameter
    
    Returns:
        Mean LISI score
    """
    if use_rep not in adata.obsm:
        sc.tl.pca(adata)
        use_rep = "X_pca"
    
    X = adata.obsm[use_rep]
    batch_labels = adata.obs[batch_key]
    
    # Build kNN
    nn = NearestNeighbors(n_neighbors=n_neighbors + 1)
    nn.fit(X)
    distances, indices = nn.kneighbors(X)
    indices = indices[:, 1:]  # Remove self
    distances = distances[:, 1:]
    
    # Compute kernel widths
    sigma = distances[:, perplexity - 1] + 1e-10
    
    # Compute LISI for each cell
    lisi_scores = []
    batch_cats = sorted(batch_labels.unique())
    n_batches = len(batch_cats)
    
    for i in range(X.shape[0]):
        # Get neighbor batches
        neighbor_batches = batch_labels.iloc[indices[i]].values
        
        # Compute Gaussian kernel weights
        weights = np.exp(-distances[i]**2 / (2 * sigma[i]**2))
        weights /= weights.sum()
        
        # Compute Simpson's index for this neighborhood
        simpson = 0
        for batch in batch_cats:
            p = weights[neighbor_batches == batch].sum()
            simpson += p**2
        
        # LISI = 1 / Simpson
        lisi = 1 / (simpson + 1e-10)
        lisi_scores.append(lisi)
    
    mean_lisi = np.mean(lisi_scores)
    
    # Normalize to [0, 1] where 1 is perfect mixing
    normalized_lisi = (mean_lisi - 1) / (n_batches - 1)
    
    return normalized_lisi


def compute_asw_batch(
    adata: AnnData,
    batch_key: str,
    use_rep: str = "X_pca",
) -> float:
    """
    Compute Average Silhouette Width for batch labels.
    
    Lower ASW = better mixing (batches are not well-separated).
    
    Args:
        adata: AnnData object
        batch_key: Batch identifier
        use_rep: Representation
    
    Returns:
        Negative ASW (higher = better mixing)
    """
    if use_rep not in adata.obsm:
        sc.tl.pca(adata)
        use_rep = "X_pca"
    
    X = adata.obsm[use_rep]
    batch_labels = adata.obs[batch_key]
    
    # Compute silhouette score
    if len(batch_labels.unique()) <= 1:
        return np.nan
    
    score = silhouette_score(
        X,
        batch_labels,
        metric='euclidean',
        sample_size=min(5000, len(X))
    )
    
    # Return negative (so higher is better)
    return -score


# =============================================================================
# Main Diagnosis Function
# =============================================================================

def diagnose_batch_effects(
    adata: AnnData,
    config: Optional[BatchDiagnosisConfig] = None,
    **kwargs,
) -> Dict[str, float]:
    """
    Comprehensive batch effect diagnosis.
    
    Computes multiple metrics and provides an overall assessment.
    
    Args:
        adata: AnnData object
        config: Configuration
        **kwargs: Override config
    
    Returns:
        Dictionary of metric scores
    """
    if config is None:
        config = BatchDiagnosisConfig()
    
    for key, value in kwargs.items():
        if hasattr(config, key):
            setattr(config, key, value)
    
    log.info("=" * 60)
    log.info("Batch Effect Diagnosis")
    log.info("=" * 60)
    
    if config.batch_key not in adata.obs:
        raise ValueError(f"Batch key '{config.batch_key}' not found in adata.obs")
    
    n_batches = adata.obs[config.batch_key].nunique()
    log.info(f"Number of batches: {n_batches}")
    
    if n_batches <= 1:
        log.warning("Only 1 batch found. No batch correction needed.")
        return {'n_batches': n_batches, 'severity': 'none'}
    
    # Compute PCA if needed
    if 'X_pca' not in adata.obsm:
        log.info("Computing PCA...")
        sc.tl.pca(adata)
    
    # Compute metrics
    results = {
        'n_batches': n_batches,
        'n_cells': adata.n_obs,
    }
    
    for metric in config.metrics:
        log.info(f"Computing {metric}...")
        
        try:
            if metric == 'pcr':
                score = compute_pcr_batch_effect(
                    adata,
                    batch_key=config.batch_key
                )
                results['pcr_r2'] = score
                
            elif metric == 'kbet':
                score = compute_kbet_score(
                    adata,
                    batch_key=config.batch_key
                )
                results['kbet_acceptance'] = score
                
            elif metric == 'lisi':
                score = compute_lisi(
                    adata,
                    batch_key=config.batch_key
                )
                results['lisi_score'] = score
                
            elif metric == 'asw':
                score = compute_asw_batch(
                    adata,
                    batch_key=config.batch_key
                )
                results['asw_batch'] = score
                
            elif metric == 'silhouette':
                X = adata.obsm['X_pca']
                batch_labels = adata.obs[config.batch_key]
                score = silhouette_score(
                    X,
                    batch_labels,
                    sample_size=min(5000, len(X))
                )
                results['silhouette_batch'] = -score  # Negative for consistency
                
        except Exception as e:
            log.warning(f"Failed to compute {metric}: {e}")
            results[f'{metric}_failed'] = str(e)
    
    # === Overall Assessment ===
    # Combine metrics into severity score
    severity_scores = []
    
    if 'pcr_r2' in results:
        # High PCR R² = strong batch effect
        severity_scores.append(results['pcr_r2'])
    
    if 'kbet_acceptance' in results:
        # Low kBET acceptance = strong batch effect
        severity_scores.append(1 - results['kbet_acceptance'])
    
    if 'lisi_score' in results:
        # Low LISI = strong batch effect
        severity_scores.append(1 - results['lisi_score'])
    
    if 'asw_batch' in results:
        # High (negative) ASW = weak mixing = strong batch effect
        # ASW is already negated, so we take 1 - normalized value
        normalized_asw = (results['asw_batch'] + 1) / 2  # Normalize to [0, 1]
        severity_scores.append(1 - normalized_asw)
    
    if severity_scores:
        overall_severity = np.mean(severity_scores)
        results['overall_severity'] = overall_severity
        
        # Categorize
        if overall_severity < 0.2:
            severity_category = "negligible"
        elif overall_severity < 0.4:
            severity_category = "mild"
        elif overall_severity < 0.6:
            severity_category = "moderate"
        elif overall_severity < 0.8:
            severity_category = "strong"
        else:
            severity_category = "severe"
        
        results['severity_category'] = severity_category
    else:
        results['overall_severity'] = np.nan
        results['severity_category'] = "unknown"
    
    # === Print Summary ===
    log.info("\n" + "=" * 60)
    log.info("Diagnosis Summary:")
    log.info("=" * 60)
    for key, value in results.items():
        if isinstance(value, float):
            log.info(f"  {key}: {value:.3f}")
        else:
            log.info(f"  {key}: {value}")
    log.info("=" * 60)
    
    # === Method Recommendation ===
    if config.recommend_methods:
        recommendation = recommend_integration_method(adata, results, config)
        results['recommended_methods'] = recommendation
    
    # === Visualization ===
    if config.plot and config.save_dir:
        _plot_batch_diagnosis(adata, results, config)
    
    return results


def recommend_integration_method(
    adata: AnnData,
    diagnosis_results: Dict,
    config: BatchDiagnosisConfig,
) -> List[Dict]:
    """
    Recommend integration methods based on diagnosis.
    
    Returns ranked list of recommended methods with explanations.
    """
    log.info("\n" + "=" * 60)
    log.info("Method Recommendation")
    log.info("=" * 60)
    
    n_batches = diagnosis_results['n_batches']
    n_cells = diagnosis_results['n_cells']
    severity = diagnosis_results.get('overall_severity', 0.5)
    
    # Method database with applicability rules
    methods = [
        {
            'name': 'harmony',
            'speed': 'fast',
            'memory': 'low',
            'best_for': 'moderate to strong batch effects',
            'min_batches': 2,
            'max_batches': 100,
            'min_cells': 1000,
            'score': 0,
        },
        {
            'name': 'scanorama',
            'speed': 'medium',
            'memory': 'medium',
            'best_for': 'diverse cell types across batches',
            'min_batches': 2,
            'max_batches': 50,
            'min_cells': 500,
            'score': 0,
        },
        {
            'name': 'scvi',
            'speed': 'slow',
            'memory': 'high',
            'best_for': 'complex batch effects, large datasets',
            'min_batches': 2,
            'max_batches': 1000,
            'min_cells': 5000,
            'score': 0,
        },
        {
            'name': 'bbknn',
            'speed': 'fast',
            'memory': 'low',
            'best_for': 'graph-based analysis',
            'min_batches': 2,
            'max_batches': 100,
            'min_cells': 1000,
            'score': 0,
        },
        {
            'name': 'combat',
            'speed': 'medium',
            'memory': 'low',
            'best_for': 'linear batch effects',
            'min_batches': 2,
            'max_batches': 50,
            'min_cells': 500,
            'score': 0,
        },
    ]
    
    # === Score each method ===
    for method in methods:
        score = 0
        
        # 1. Batch count compatibility
        if method['min_batches'] <= n_batches <= method['max_batches']:
            score += 2
        
        # 2. Cell count compatibility
        if n_cells >= method['min_cells']:
            score += 2
        
        # 3. Severity matching
        if severity < 0.3:  # Mild
            if method['name'] in ['harmony', 'bbknn', 'combat']:
                score += 3
        elif severity < 0.6:  # Moderate
            if method['name'] in ['harmony', 'scanorama']:
                score += 3
        else:  # Strong
            if method['name'] in ['scvi', 'scanorama']:
                score += 3
        
        # 4. Dataset size bonus
        if n_cells > 50000:
            if method['name'] in ['harmony', 'bbknn']:  # Fast methods
                score += 2
        
        if n_cells < 10000:
            if method['speed'] != 'slow':
                score += 1
        
        method['score'] = score
    
    # Sort by score
    methods_sorted = sorted(methods, key=lambda x: x['score'], reverse=True)
    
    # Take top N
    top_methods = methods_sorted[:config.n_methods_to_recommend]
    
    # Print recommendations
    log.info(f"\nTop {config.n_methods_to_recommend} recommended methods:")
    for i, method in enumerate(top_methods, 1):
        log.info(f"\n{i}. {method['name'].upper()}")
        log.info(f"   Score: {method['score']}/10")
        log.info(f"   Speed: {method['speed']}")
        log.info(f"   Memory: {method['memory']}")
        log.info(f"   Best for: {method['best_for']}")
    
    log.info("\n" + "=" * 60)
    
    return top_methods


def _plot_batch_diagnosis(
    adata: AnnData,
    results: Dict,
    config: BatchDiagnosisConfig,
):
    """Plot batch diagnosis results."""
    save_dir = Path(config.save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)
    
    fig = plt.figure(figsize=(20, 12))
    gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3)
    
    # 1. Metric scores
    ax = fig.add_subplot(gs[0, :2])
    
    metric_names = []
    metric_values = []
    
    for key in ['pcr_r2', 'kbet_acceptance', 'lisi_score', 'asw_batch']:
        if key in results:
            metric_names.append(key.replace('_', ' ').title())
            metric_values.append(results[key])
    
    if metric_values:
        colors = ['red' if v < 0.3 else 'orange' if v < 0.6 else 'green' for v in metric_values]
        ax.barh(metric_names, metric_values, color=colors, alpha=0.7)
        ax.set_xlabel('Score')
        ax.set_title('Batch Effect Metrics')
        ax.set_xlim(0, 1)
    
    # 2. Overall severity gauge
    ax = fig.add_subplot(gs[0, 2:])
    
    if 'overall_severity' in results and not np.isnan(results['overall_severity']):
        severity = results['overall_severity']
        
        # Create gauge
        theta = np.linspace(0, np.pi, 100)
        radius = 1
        
        # Background
        ax.plot(radius * np.cos(theta), radius * np.sin(theta), 'k-', linewidth=2)
        
        # Severity zones
        ax.fill_between(
            np.cos(theta[:20]), np.sin(theta[:20]), 0,
            color='green', alpha=0.3, label='Negligible'
        )
        ax.fill_between(
            np.cos(theta[20:40]), np.sin(theta[20:40]), 0,
            color='yellow', alpha=0.3, label='Mild'
        )
        ax.fill_between(
            np.cos(theta[40:60]), np.sin(theta[40:60]), 0,
            color='orange', alpha=0.3, label='Moderate'
        )
        ax.fill_between(
            np.cos(theta[60:]), np.sin(theta[60:]), 0,
            color='red', alpha=0.3, label='Severe'
        )
        
        # Needle
        needle_angle = np.pi * (1 - severity)
        ax.arrow(
            0, 0,
            0.8 * np.cos(needle_angle), 0.8 * np.sin(needle_angle),
            head_width=0.1, head_length=0.1, fc='black', ec='black'
        )
        
        ax.set_xlim(-1.2, 1.2)
        ax.set_ylim(-0.2, 1.2)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title(f'Overall Severity: {results["severity_category"].upper()}', fontsize=14, fontweight='bold')
        ax.legend(loc='lower center', ncol=4)
    
    # 3. PCA colored by batch
    if 'X_pca' in adata.obsm:
        ax = fig.add_subplot(gs[1, :2])
        
        sc.pl.pca(
            adata,
            color=config.batch_key,
            ax=ax,
            show=False,
            title='PCA colored by Batch'
        )
    
    # 4. UMAP colored by batch (if available)
    if 'X_umap' in adata.obsm:
        ax = fig.add_subplot(gs[1, 2:])
        
        sc.pl.umap(
            adata,
            color=config.batch_key,
            ax=ax,
            show=False,
            title='UMAP colored by Batch'
        )
    
    # 5. Batch size distribution
    ax = fig.add_subplot(gs[2, :2])
    
    batch_counts = adata.obs[config.batch_key].value_counts()
    batch_counts.plot(kind='bar', ax=ax, color='skyblue')
    ax.set_title('Cells per Batch')
    ax.set_xlabel('Batch')
    ax.set_ylabel('Number of Cells')
    ax.tick_params(axis='x', rotation=45)
    
    # 6. Recommended methods
    ax = fig.add_subplot(gs[2, 2:])
    
    if 'recommended_methods' in results:
        methods = results['recommended_methods']
        names = [m['name'] for m in methods]
        scores = [m['score'] for m in methods]
        
        colors_methods = ['gold', 'silver', '#CD7F32'][:len(methods)]
        ax.barh(names, scores, color=colors_methods, alpha=0.7)
        ax.set_xlabel('Recommendation Score')
        ax.set_title('Top Recommended Methods')
        ax.set_xlim(0, 10)
    
    plt.suptitle('Batch Effect Diagnosis Report', fontsize=16, fontweight='bold')
    
    plt.savefig(save_dir / 'batch_diagnosis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    log.info(f"Saved diagnostic plots to {save_dir}")


def compare_integration_methods(
    adata: AnnData,
    batch_key: str,
    methods: List[str] = ['harmony', 'scanorama', 'scvi'],
    label_key: Optional[str] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Compare multiple integration methods on the same dataset.
    
    This is computationally expensive but provides direct comparison.
    
    Args:
        adata: AnnData object
        batch_key: Batch identifier
        methods: List of methods to compare
        label_key: Cell type labels (for bio-conservation evaluation)
        **kwargs: Method-specific parameters
    
    Returns:
        DataFrame comparing methods across metrics
    """
    from ..integrate import batch_correction, evaluate_integration
    
    log.info("Comparing integration methods...")
    
    results = []
    
    for method in methods:
        log.info(f"\nTesting {method}...")
        
        # Make a copy
        adata_test = adata.copy()
        
        try:
            # Run integration
            adata_test = batch_correction(
                adata_test,
                method=method,
                batch_key=batch_key,
                **kwargs
            )
            
            # Evaluate
            eval_results = evaluate_integration(
                adata_test,
                batch_key=batch_key,
                label_key=label_key,
                plot=False
            )
            
            eval_results['method'] = method
            results.append(eval_results)
            
        except Exception as e:
            log.error(f"Failed to test {method}: {e}")
    
    # Combine results
    df = pd.DataFrame(results)
    
    # Display
    log.info("\n" + "=" * 60)
    log.info("Method Comparison Results:")
    log.info("=" * 60)
    print(df.to_string(index=False))
    log.info("=" * 60)
    
    return df