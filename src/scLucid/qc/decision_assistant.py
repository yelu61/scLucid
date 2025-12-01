# src/sclucid/qc/decision_assistant.py (新文件)
"""
Interactive QC decision assistant using reinforcement learning concepts.
"""

import logging
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from anndata import AnnData
from dataclasses import dataclass

log = logging.getLogger(__name__)


@dataclass
class QCDecision:
    """A single QC filtering decision with its consequences."""
    threshold_name: str
    threshold_value: float
    cells_removed: int
    genes_affected: int
    cluster_coherence_change: float  # Change in silhouette score
    marker_retention: float  # % of known marker genes retained


class QCDecisionAssistant:
    """
    Intelligent assistant for optimizing QC thresholds.
    
    Uses a simulation-based approach to predict the consequences
    of different QC decisions on downstream analysis.
    """
    
    def __init__(
        self,
        adata: AnnData,
        known_markers: Optional[Dict[str, List[str]]] = None
    ):
        """
        Initialize decision assistant.
        
        Args:
            adata: AnnData object (before filtering)
            known_markers: Dict mapping cell types to marker genes
        """
        self.adata = adata.copy()
        self.known_markers = known_markers or {}
        self.decision_history = []
        
    def simulate_threshold_impact(
        self,
        threshold_dict: Dict[str, float],
        n_simulations: int = 5
    ) -> Dict[str, float]:
        """
        Simulate the impact of proposed thresholds on clustering quality.
        
        This runs a quick clustering on the filtered data to assess quality.
        
        Returns:
            Dictionary with impact metrics
        """
        import scanpy as sc
        from sklearn.metrics import silhouette_score
        
        if self.adata.n_obs > 10000:
            if not hasattr(self, '_simulation_subset'):
                log.info("Subsampling data to 10k cells for rapid simulation...")
                sc.pp.subsample(self.adata, n_obs=10000, copy=True)
        
        # Apply thresholds to create filtered data
        adata_filtered = self.adata.copy()
        
        keep_mask = pd.Series(True, index=adata_filtered.obs_names)
        
        for metric, threshold_value in threshold_dict.items():
            if metric == 'min_genes':
                keep_mask &= adata_filtered.obs['n_genes_by_counts'] >= threshold_value
            elif metric == 'max_mt':
                keep_mask &= adata_filtered.obs['pct_counts_mt'] <= threshold_value
            # Add more metrics as needed
        
        adata_filtered = adata_filtered[keep_mask].copy()
        
        # Quick preprocessing and clustering
        sc.pp.normalize_total(adata_filtered, target_sum=1e4)
        sc.pp.log1p(adata_filtered)
        sc.pp.highly_variable_genes(adata_filtered, n_top_genes=2000)
        sc.pp.pca(adata_filtered, n_comps=30)
        sc.pp.neighbors(adata_filtered, n_neighbors=15)
        sc.tl.leiden(adata_filtered, resolution=0.5)
        
        # Calculate quality metrics
        metrics = {
            'cells_remaining': adata_filtered.n_obs,
            'cells_removed': self.adata.n_obs - adata_filtered.n_obs,
            'removal_rate': (self.adata.n_obs - adata_filtered.n_obs) / self.adata.n_obs,
            'n_clusters': len(adata_filtered.obs['leiden'].unique())
        }
        
        # Silhouette score (clustering coherence)
        try:
            X_pca = adata_filtered.obsm['X_pca']
            labels = adata_filtered.obs['leiden'].astype('category').cat.codes
            metrics['silhouette_score'] = silhouette_score(X_pca, labels)
        except:
            metrics['silhouette_score'] = np.nan
        
        # Marker gene retention
        if self.known_markers:
            all_markers = set(sum(self.known_markers.values(), []))
            markers_retained = all_markers & set(adata_filtered.var_names)
            metrics['marker_retention'] = len(markers_retained) / len(all_markers)
        else:
            metrics['marker_retention'] = np.nan
        
        return metrics
    
    def suggest_optimal_thresholds(
        self,
        metric_ranges: Dict[str, Tuple[float, float, float]],
        optimization_target: str = 'silhouette_score',
        constraints: Optional[Dict[str, float]] = None
    ) -> Dict[str, float]:
        """
        Use grid search to find optimal thresholds.
        
        Args:
            metric_ranges: Dict mapping metric names to (min, max, step) tuples
            optimization_target: Metric to optimize ('silhouette_score', 'marker_retention')
            constraints: Dict of constraints (e.g., {'removal_rate': 0.3})
            
        Returns:
            Dictionary of optimal thresholds
        """
        from itertools import product
        
        log.info("Starting threshold optimization...")
        log.info(f"Optimization target: {optimization_target}")
        
        # Generate all combinations
        metric_names = list(metric_ranges.keys())
        value_ranges = [
            np.arange(min_val, max_val, step)
            for min_val, max_val, step in metric_ranges.values()
        ]
        
        all_combinations = list(product(*value_ranges))
        log.info(f"Testing {len(all_combinations)} threshold combinations...")
        
        best_score = -np.inf
        best_thresholds = None
        results = []
        
        for i, values in enumerate(all_combinations):
            if i % 10 == 0:
                log.info(f"Progress: {i}/{len(all_combinations)}")
            
            threshold_dict = dict(zip(metric_names, values))
            
            try:
                metrics = self.simulate_threshold_impact(threshold_dict)
                
                # Check constraints
                if constraints:
                    violates_constraint = False
                    for constraint_name, constraint_value in constraints.items():
                        if metrics[constraint_name] > constraint_value:
                            violates_constraint = True
                            break
                    
                    if violates_constraint:
                        continue
                
                # Evaluate score
                score = metrics.get(optimization_target, -np.inf)
                
                results.append({
                    **threshold_dict,
                    **metrics,
                    'optimization_score': score
                })
                
                if score > best_score:
                    best_score = score
                    best_thresholds = threshold_dict.copy()
                
            except Exception as e:
                log.warning(f"Simulation failed for {threshold_dict}: {e}")
                continue
        
        results_df = pd.DataFrame(results)
        self.optimization_results = results_df
        
        log.info(f"Optimization complete!")
        log.info(f"Best thresholds: {best_thresholds}")
        log.info(f"Best {optimization_target}: {best_score:.4f}")
        
        return best_thresholds
    
    def plot_optimization_results(self, save_path: Optional[str] = None):
        """
        Visualize the relationship between thresholds and outcomes.
        """
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        if not hasattr(self, 'optimization_results'):
            log.warning("No optimization results to plot. Run suggest_optimal_thresholds first.")
            return
        
        df = self.optimization_results
        
        # Create a comprehensive visualization
        fig = plt.figure(figsize=(15, 10))
        gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
        
        # 1. Pareto front: removal_rate vs silhouette_score
        ax1 = fig.add_subplot(gs[0, :2])
        scatter = ax1.scatter(
            df['removal_rate'] * 100,
            df['silhouette_score'],
            c=df['marker_retention'] * 100,
            s=50,
            alpha=0.6,
            cmap='viridis'
        )
        ax1.set_xlabel('Cells Removed (%)')
        ax1.set_ylabel('Silhouette Score')
        ax1.set_title('QC Trade-off: Quality vs. Data Retention')
        plt.colorbar(scatter, ax=ax1, label='Marker Retention (%)')
        
        # 2. Distribution of silhouette scores
        ax2 = fig.add_subplot(gs[0, 2])
        ax2.hist(df['silhouette_score'].dropna(), bins=30, edgecolor='black')
        ax2.axvline(
            df['silhouette_score'].max(),
            color='red',
            linestyle='--',
            label='Best'
        )
        ax2.set_xlabel('Silhouette Score')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Score Distribution')
        ax2.legend()
        
        # 3-5. Threshold impact plots
        threshold_cols = [col for col in df.columns if col not in [
            'cells_remaining', 'cells_removed', 'removal_rate',
            'n_clusters', 'silhouette_score', 'marker_retention',
            'optimization_score'
        ]]
        
        for i, col in enumerate(threshold_cols[:3]):
            ax = fig.add_subplot(gs[1, i])
            ax.scatter(df[col], df['silhouette_score'], alpha=0.5, s=30)
            ax.set_xlabel(col.replace('_', ' ').title())
            ax.set_ylabel('Silhouette Score')
            ax.set_title(f'Impact of {col}')
        
        # 6. Cells removed vs clusters formed
        ax6 = fig.add_subplot(gs[2, 0])
        ax6.scatter(df['cells_removed'], df['n_clusters'], alpha=0.5, s=30)
        ax6.set_xlabel('Cells Removed')
        ax6.set_ylabel('Number of Clusters')
        ax6.set_title('Filtering Impact on Clustering')
        
        # 7. Marker retention across strategies
        ax7 = fig.add_subplot(gs[2, 1:])
        top_strategies = df.nlargest(20, 'optimization_score')
        ax7.barh(range(len(top_strategies)), top_strategies['marker_retention'] * 100)
        ax7.set_xlabel('Marker Retention (%)')
        ax7.set_ylabel('Strategy Rank')
        ax7.set_title('Top 20 Strategies by Marker Retention')
        
        plt.suptitle('QC Threshold Optimization Results', fontsize=16, fontweight='bold')
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            log.info(f"Saved optimization plot to {save_path}")
        
        return fig


# Usage example
def run_intelligent_qc(
    adata: AnnData,
    known_markers: Optional[Dict[str, List[str]]] = None,
    max_removal_rate: float = 0.3,
    save_dir: Optional[str] = None
) -> Tuple[AnnData, Dict[str, float]]:
    """
    Run QC with intelligent threshold optimization.
    
    Returns:
        Tuple of (filtered_adata, optimal_thresholds)
    """
    assistant = QCDecisionAssistant(adata, known_markers)
    
    # Define search space
    metric_ranges = {
        'min_genes': (200, 1000, 100),  # (min, max, step)
        'max_mt': (5, 25, 2.5)
    }
    
    # Find optimal thresholds
    optimal_thresholds = assistant.suggest_optimal_thresholds(
        metric_ranges=metric_ranges,
        optimization_target='silhouette_score',
        constraints={'removal_rate': max_removal_rate}
    )
    
    # Visualize results
    if save_dir:
        from pathlib import Path
        Path(save_dir).mkdir(parents=True, exist_ok=True)
        assistant.plot_optimization_results(
            save_path=Path(save_dir) / 'qc_optimization.png'
        )
    
    # Apply optimal thresholds
    # (Use existing filtering functions with optimal_thresholds)
    
    return adata, optimal_thresholds