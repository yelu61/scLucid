"""
CellChat: Python toolkit for inference and analysis of cell-cell communication
from single-cell and spatially resolved transcriptomics
"""

__version__ = "2.1.0"

from .core import CellChat
from .database import CellChatDB
from .preprocessing import preprocess_expression
from .visualization import *
from .analysis import *

__all__ = [
    'CellChat',
    'CellChatDB',
    'preprocess_expression'
]

# cellchat/core.py
"""
Core CellChat class implementation
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from typing import Optional, Dict, List, Union, Tuple
import warnings
from dataclasses import dataclass, field

@dataclass
class CellChatConfig:
    """Configuration for CellChat analysis"""
    species: str = "human"
    min_cells: int = 10
    thresh: float = 0.05
    population_size: bool = False
    distance_threshold: Optional[float] = None
    spatial_factors: Dict[str, float] = field(default_factory=dict)
    
class CellChat:
    """
    Main CellChat class for cell-cell communication analysis
    
    Parameters
    ----------
    data : Union[np.ndarray, sparse.spmatrix, pd.DataFrame]
        Gene expression matrix (genes x cells)
    meta : pd.DataFrame
        Cell metadata with cell labels
    group_by : str
        Column name in meta for cell grouping
    spatial_coords : Optional[pd.DataFrame]
        Spatial coordinates for spatial transcriptomics (cells x coordinates)
    """
    
    def __init__(
        self,
        data: Union[np.ndarray, sparse.spmatrix, pd.DataFrame],
        meta: pd.DataFrame,
        group_by: str,
        spatial_coords: Optional[pd.DataFrame] = None,
        config: Optional[CellChatConfig] = None
    ):
        self.config = config or CellChatConfig()
        
        # Process expression data
        if isinstance(data, pd.DataFrame):
            self.gene_names = data.index.tolist()
            self.cell_names = data.columns.tolist()
            self.data_expr = data.values
        else:
            self.gene_names = None
            self.cell_names = None
            self.data_expr = data
        
        # Store metadata
        self.meta = meta.copy()
        self.group_by = group_by
        self.cell_groups = meta[group_by].values
        self.unique_groups = np.unique(self.cell_groups)
        
        # Spatial information
        self.spatial_coords = spatial_coords
        self.is_spatial = spatial_coords is not None
        
        # Initialize containers for results
        self.db = None
        self.LR = None  # Ligand-receptor pairs
        self.net = {}  # Communication network
        self.netP = {}  # Pathway-level network
        self.idents = None
        self.data_signaling = None
        self.data_project = None
        self.distance_matrix = None
        
    def set_database(self, db: 'CellChatDB'):
        """Set CellChatDB database"""
        self.db = db
        self.LR = db.interaction
        
    def preprocess_data(
        self,
        subset_data: bool = True,
        do_sparse: bool = True
    ):
        """
        Preprocess expression data
        
        Parameters
        ----------
        subset_data : bool
            Whether to subset data to genes in database
        do_sparse : bool
            Whether to convert to sparse matrix
        """
        if self.db is None:
            raise ValueError("Please set database first using set_database()")
        
        # Get genes in database
        db_genes = set(self.db.get_all_genes())
        
        if self.gene_names is not None:
            available_genes = set(self.gene_names)
            common_genes = list(db_genes & available_genes)
            
            if subset_data:
                # Subset to common genes
                gene_idx = [self.gene_names.index(g) for g in common_genes]
                self.data_expr = self.data_expr[gene_idx, :]
                self.gene_names = common_genes
        
        # Convert to sparse if needed
        if do_sparse and not sparse.issparse(self.data_expr):
            self.data_expr = sparse.csr_matrix(self.data_expr)
        
        print(f"Preprocessed data: {self.data_expr.shape[0]} genes x {self.data_expr.shape[1]} cells")
        
    def identify_overexpressed_genes(self, thresh: float = 0.05):
        """
        Identify overexpressed genes for each cell group
        
        Parameters
        ----------
        thresh : float
            Threshold for detecting expression
        """
        n_genes = self.data_expr.shape[0]
        n_groups = len(self.unique_groups)
        
        # Convert to dense for computation
        if sparse.issparse(self.data_expr):
            data = self.data_expr.toarray()
        else:
            data = self.data_expr
        
        # Calculate mean expression for each group
        group_means = np.zeros((n_genes, n_groups))
        for i, group in enumerate(self.unique_groups):
            mask = self.cell_groups == group
            group_means[:, i] = data[:, mask].mean(axis=1)
        
        # Identify overexpressed genes
        self.data_signaling = group_means
        
        return group_means
    
    def compute_communication_prob(
        self,
        type: str = "truncatedMean",
        trim: float = 0.1,
        population_size: bool = False,
        distance_threshold: Optional[float] = None
    ):
        """
        Compute communication probability
        
        Parameters
        ----------
        type : str
            Method for computing average expression ("truncatedMean", "median", "mean")
        trim : float
            Fraction to trim for truncated mean
        population_size : bool
            Whether to consider population size
        distance_threshold : Optional[float]
            Distance threshold for spatial data
        """
        if self.data_signaling is None:
            self.identify_overexpressed_genes()
        
        n_groups = len(self.unique_groups)
        n_interactions = len(self.LR)
        
        # Initialize probability matrix
        prob = np.zeros((n_interactions, n_groups, n_groups))
        pval = np.ones((n_interactions, n_groups, n_groups))
        
        # Compute distance matrix if spatial
        if self.is_spatial and distance_threshold is not None:
            self.distance_matrix = self._compute_distance_matrix()
        
        # For each LR pair
        for idx, (lr_name, lr_info) in enumerate(self.LR.iterrows()):
            ligand = lr_info['ligand']
            receptor = lr_info['receptor']
            
            # Get expression levels
            ligand_expr = self._get_gene_expression(ligand)
            receptor_expr = self._get_gene_expression(receptor)
            
            if ligand_expr is None or receptor_expr is None:
                continue
            
            # Compute communication probability
            for i, group_i in enumerate(self.unique_groups):
                for j, group_j in enumerate(self.unique_groups):
                    # Ligand from group i, receptor in group j
                    L_expr = ligand_expr[i]
                    R_expr = receptor_expr[j]
                    
                    # Apply spatial constraint if needed
                    spatial_weight = 1.0
                    if self.is_spatial and distance_threshold is not None:
                        spatial_weight = self._compute_spatial_weight(
                            i, j, distance_threshold
                        )
                    
                    # Compute probability using mass action law
                    prob[idx, i, j] = L_expr * R_expr * spatial_weight
        
        # Store results
        self.net['prob'] = prob
        self.net['pval'] = pval
        
        # Compute pathway-level probability
        self._compute_pathway_prob()
        
        return prob
    
    def _get_gene_expression(self, gene_symbol: str) -> Optional[np.ndarray]:
        """Get expression of a gene across all groups"""
        if self.gene_names is None or gene_symbol not in self.gene_names:
            return None
        
        gene_idx = self.gene_names.index(gene_symbol)
        
        if sparse.issparse(self.data_expr):
            expr = self.data_expr[gene_idx, :].toarray().flatten()
        else:
            expr = self.data_expr[gene_idx, :]
        
        # Compute group averages
        group_expr = np.zeros(len(self.unique_groups))
        for i, group in enumerate(self.unique_groups):
            mask = self.cell_groups == group
            group_expr[i] = expr[mask].mean()
        
        return group_expr
    
    def _compute_distance_matrix(self) -> np.ndarray:
        """Compute pairwise distance matrix between cell groups"""
        from scipy.spatial.distance import cdist
        
        n_groups = len(self.unique_groups)
        centroids = np.zeros((n_groups, self.spatial_coords.shape[1]))
        
        # Compute centroids for each group
        for i, group in enumerate(self.unique_groups):
            mask = self.cell_groups == group
            centroids[i] = self.spatial_coords[mask].mean(axis=0)
        
        # Compute pairwise distances
        dist_matrix = cdist(centroids, centroids, metric='euclidean')
        
        return dist_matrix
    
    def _compute_spatial_weight(
        self,
        group_i: int,
        group_j: int,
        distance_threshold: float
    ) -> float:
        """Compute spatial weight based on distance"""
        if self.distance_matrix is None:
            return 1.0
        
        dist = self.distance_matrix[group_i, group_j]
        
        # Exponential decay
        weight = np.exp(-dist / distance_threshold)
        
        return weight
    
    def _compute_pathway_prob(self):
        """Compute pathway-level communication probability"""
        if 'prob' not in self.net:
            raise ValueError("Run compute_communication_prob first")
        
        # Group interactions by pathway
        pathways = self.LR['pathway_name'].unique()
        n_groups = len(self.unique_groups)
        
        pathway_prob = {}
        pathway_pval = {}
        
        for pathway in pathways:
            # Get interactions in this pathway
            pathway_mask = self.LR['pathway_name'] == pathway
            pathway_indices = np.where(pathway_mask)[0]
            
            # Aggregate probabilities
            prob_pathway = self.net['prob'][pathway_indices].sum(axis=0)
            pval_pathway = self.net['pval'][pathway_indices].min(axis=0)
            
            pathway_prob[pathway] = prob_pathway
            pathway_pval[pathway] = pval_pathway
        
        self.netP['prob'] = pathway_prob
        self.netP['pval'] = pathway_pval
    
    def filter_communication(
        self,
        min_cells: int = 10,
        thresh_p: float = 0.05
    ):
        """
        Filter communications based on significance
        
        Parameters
        ----------
        min_cells : int
            Minimum number of cells in a group
        thresh_p : float
            P-value threshold
        """
        # Filter by cell number
        group_sizes = pd.Series(self.cell_groups).value_counts()
        valid_groups = group_sizes[group_sizes >= min_cells].index.tolist()
        
        # Filter by p-value
        if 'pval' in self.net:
            prob = self.net['prob']
            pval = self.net['pval']
            
            # Set non-significant to 0
            prob[pval > thresh_p] = 0
            self.net['prob'] = prob
    
    def compute_network_centrality(self):
        """Compute network centrality measures"""
        from .analysis import compute_centrality
        
        if 'prob' not in self.net:
            raise ValueError("Run compute_communication_prob first")
        
        centrality_results = compute_centrality(self.net['prob'], self.unique_groups)
        self.net['centrality'] = centrality_results
        
        return centrality_results
    
    def identify_signaling_roles(self, pattern: str = "outgoing"):
        """
        Identify signaling roles of cell groups
        
        Parameters
        ----------
        pattern : str
            "outgoing", "incoming", or "overall"
        """
        from .analysis import identify_roles
        
        if 'prob' not in self.netP:
            raise ValueError("No pathway-level network computed")
        
        roles = identify_roles(self.netP['prob'], pattern=pattern)
        
        return roles
    
    def compare_interactions(
        self,
        other: 'CellChat',
        comparison_type: str = "functional"
    ):
        """
        Compare with another CellChat object
        
        Parameters
        ----------
        other : CellChat
            Another CellChat object to compare
        comparison_type : str
            Type of comparison ("functional", "structural")
        """
        from .comparison import compare_cellchat_objects
        
        comparison_results = compare_cellchat_objects(
            self, other, comparison_type=comparison_type
        )
        
        return comparison_results
    
    def update_cellchat(self):
        """Update CellChat object to latest version"""
        # Update spatial_factors if using old scale_factors
        if hasattr(self, 'scale_factors'):
            self.spatial_factors = self.scale_factors
            delattr(self, 'scale_factors')
        
        # Update slices to samples
        if 'slices' in self.meta.columns:
            self.meta['samples'] = self.meta['slices']
            self.meta.drop('slices', axis=1, inplace=True)
        
        print("CellChat object updated to version 2.1.0")
 

# cellchat/database.py
"""
CellChatDB database implementation
"""

import pandas as pd
import numpy as np
from typing import Optional, List, Dict
import pkg_resources
import json

class CellChatDB:
    """
    CellChatDB: Database of ligand-receptor interactions
    
    Parameters
    ----------
    species : str
        Species ("human" or "mouse")
    version : str
        Database version
    """
    
    def __init__(self, species: str = "human", version: str = "v2"):
        self.species = species
        self.version = version
        
        # Load database
        self.interaction = self._load_interaction_db()
        self.complex = self._load_complex_db()
        self.cofactor = self._load_cofactor_db()
        self.geneInfo = self._load_gene_info()
        
    def _load_interaction_db(self) -> pd.DataFrame:
        """Load ligand-receptor interaction database"""
        # This would load from package data
        # For demonstration, creating a sample structure
        
        interactions = pd.DataFrame({
            'interaction_name': [],
            'pathway_name': [],
            'ligand': [],
            'receptor': [],
            'evidence': [],
            'annotation': [],
            'interaction_type': [],  # Secreted, ECM-Receptor, Cell-Cell Contact
        })
        
        # In real implementation, load from CSV/JSON files
        # interactions = pd.read_csv('data/interaction_input_{}.csv'.format(self.species))
        
        return interactions
    
    def _load_complex_db(self) -> pd.DataFrame:
        """Load complex composition database"""
        complex_db = pd.DataFrame({
            'complex_name': [],
            'subunit': [],
            'subunit_type': [],  # core, accessory
        })
        
        return complex_db
    
    def _load_cofactor_db(self) -> pd.DataFrame:
        """Load cofactor database"""
        cofactor_db = pd.DataFrame({
            'cofactor': [],
            'interaction_name': [],
        })
        
        return cofactor_db
    
    def _load_gene_info(self) -> pd.DataFrame:
        """Load gene information"""
        gene_info = pd.DataFrame({
            'Symbol': [],
            'GeneName': [],
            'GeneID': [],
        })
        
        return gene_info
    
    def get_all_genes(self) -> List[str]:
        """Get all genes in database"""
        genes = set()
        
        # Add ligands
        genes.update(self.interaction['ligand'].dropna().unique())
        
        # Add receptors
        genes.update(self.interaction['receptor'].dropna().unique())
        
        # Add complex subunits
        if not self.complex.empty:
            genes.update(self.complex['subunit'].dropna().unique())
        
        # Add cofactors
        if not self.cofactor.empty:
            genes.update(self.cofactor['cofactor'].dropna().unique())
        
        return list(genes)
    
    def subset_db(
        self,
        interaction_types: Optional[List[str]] = None,
        pathways: Optional[List[str]] = None
    ) -> 'CellChatDB':
        """
        Subset database by interaction types or pathways
        
        Parameters
        ----------
        interaction_types : Optional[List[str]]
            List of interaction types to keep
        pathways : Optional[List[str]]
            List of pathways to keep
        """
        subset_db = CellChatDB(species=self.species, version=self.version)
        subset_db.interaction = self.interaction.copy()
        
        if interaction_types is not None:
            subset_db.interaction = subset_db.interaction[
                subset_db.interaction['interaction_type'].isin(interaction_types)
            ]
        
        if pathways is not None:
            subset_db.interaction = subset_db.interaction[
                subset_db.interaction['pathway_name'].isin(pathways)
            ]
        
        return subset_db
    
    def update_db(self, custom_interactions: pd.DataFrame):
        """
        Update database with custom interactions
        
        Parameters
        ----------
        custom_interactions : pd.DataFrame
            Custom interactions to add
        """
        required_cols = ['interaction_name', 'pathway_name', 'ligand', 'receptor']
        
        if not all(col in custom_interactions.columns for col in required_cols):
            raise ValueError(f"Custom interactions must contain: {required_cols}")
        
        self.interaction = pd.concat([
            self.interaction,
            custom_interactions
        ], ignore_index=True)
        
        print(f"Added {len(custom_interactions)} custom interactions")
    
    def search_interaction(
        self,
        ligand: Optional[str] = None,
        receptor: Optional[str] = None,
        pathway: Optional[str] = None
    ) -> pd.DataFrame:
        """Search interactions by ligand, receptor, or pathway"""
        results = self.interaction.copy()
        
        if ligand is not None:
            results = results[results['ligand'] == ligand]
        
        if receptor is not None:
            results = results[results['receptor'] == receptor]
        
        if pathway is not None:
            results = results[results['pathway_name'].str.contains(pathway, case=False)]
        
        return results

def update_cellchat_db(
    db: CellChatDB,
    new_interactions: pd.DataFrame
) -> CellChatDB:
    """
    Update CellChatDB with new interactions
    
    Parameters
    ----------
    db : CellChatDB
        Original database
    new_interactions : pd.DataFrame
        New interactions to add
    """
    updated_db = CellChatDB(species=db.species, version=db.version)
    updated_db.interaction = pd.concat([
        db.interaction,
        new_interactions
    ], ignore_index=True).drop_duplicates()
    
    return updated_db

# cellchat/analysis.py
"""
Analysis functions for CellChat
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import NMF
from sklearn.manifold import MDS
from typing import Dict, List, Tuple, Optional

def compute_centrality(
    prob_matrix: np.ndarray,
    group_names: List[str]
) -> Dict[str, pd.DataFrame]:
    """
    Compute network centrality measures
    
    Parameters
    ----------
    prob_matrix : np.ndarray
        Probability matrix (n_interactions x n_groups x n_groups)
    group_names : List[str]
        Names of cell groups
    
    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary containing centrality measures
    """
    n_groups = len(group_names)
    
    # Sum across all interactions to get overall network
    network = prob_matrix.sum(axis=0)
    
    # Compute out-degree (outgoing signals)
    out_degree = network.sum(axis=1)
    
    # Compute in-degree (incoming signals)
    in_degree = network.sum(axis=0)
    
    # Compute betweenness centrality (simplified version)
    betweenness = np.zeros(n_groups)
    for i in range(n_groups):
        # Simplified: sum of signals passing through node i
        betweenness[i] = (network[i, :].sum() + network[:, i].sum()) / 2
    
    # Normalize
    total = network.sum()
    if total > 0:
        out_degree = out_degree / total
        in_degree = in_degree / total
        betweenness = betweenness / total
    
    # Create results DataFrame
    results = pd.DataFrame({
        'group': group_names,
        'out_degree': out_degree,
        'in_degree': in_degree,
        'betweenness': betweenness,
        'total_strength': out_degree + in_degree
    })
    
    return {
        'centrality': results,
        'network_matrix': network
    }

def identify_roles(
    pathway_prob: Dict[str, np.ndarray],
    pattern: str = "outgoing",
    k: int = 5
) -> Dict:
    """
    Identify signaling roles using pattern recognition
    
    Parameters
    ----------
    pathway_prob : Dict[str, np.ndarray]
        Pathway-level probability matrices
    pattern : str
        "outgoing", "incoming", or "overall"
    k : int
        Number of patterns to identify
    """
    # Concatenate all pathway matrices
    pathways = list(pathway_prob.keys())
    n_pathways = len(pathways)
    n_groups = pathway_prob[pathways[0]].shape[0]
    
    # Create feature matrix
    if pattern == "outgoing":
        features = np.zeros((n_groups, n_pathways))
        for i, pathway in enumerate(pathways):
            features[:, i] = pathway_prob[pathway].sum(axis=1)
    elif pattern == "incoming":
        features = np.zeros((n_groups, n_pathways))
        for i, pathway in enumerate(pathways):
            features[:, i] = pathway_prob[pathway].sum(axis=0)
    else:  # overall
        features = np.zeros((n_groups, n_pathways))
        for i, pathway in enumerate(pathways):
            features[:, i] = pathway_prob[pathway].sum(axis=0) + pathway_prob[pathway].sum(axis=1)
    
    # Apply NMF for pattern recognition
    nmf = NMF(n_components=k, init='nndsvda', random_state=42)
    W = nmf.fit_transform(features)  # Cell group patterns
    H = nmf.components_  # Pathway patterns
    
    # Assign dominant pattern to each cell group
    dominant_pattern = W.argmax(axis=1)
    
    results = {
        'patterns': W,
        'pathway_patterns': H,
        'dominant_pattern': dominant_pattern,
        'feature_matrix': features,
        'pathways': pathways
    }
    
    return results

def identify_signaling_patterns(
    cellchat_obj,
    pattern: str = "outgoing",
    k: int = 5,
    height: float = 10
) -> Dict:
    """
    Identify and cluster signaling patterns
    
    Parameters
    ----------
    cellchat_obj : CellChat
        CellChat object
    pattern : str
        Pattern type
    k : int
        Number of patterns
    height : float
        Height for dendrogram cutting
    """
    from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
    from scipy.spatial.distance import pdist
    
    # Get roles
    roles = identify_roles(cellchat_obj.netP['prob'], pattern=pattern, k=k)
    
    # Hierarchical clustering
    distance_matrix = pdist(roles['patterns'], metric='euclidean')
    linkage_matrix = linkage(distance_matrix, method='ward')
    
    # Cut dendrogram
    clusters = fcluster(linkage_matrix, height, criterion='distance')
    
    results = {
        **roles,
        'linkage': linkage_matrix,
        'clusters': clusters
    }
    
    return results

def compute_network_similarity(
    cellchat1,
    cellchat2,
    type: str = "functional"
) -> Dict:
    """
    Compute similarity between two networks
    
    Parameters
    ----------
    cellchat1, cellchat2 : CellChat
        CellChat objects to compare
    type : str
        "functional" or "structural"
    """
    if type == "functional":
        # Compare pathway activities
        pathways1 = set(cellchat1.netP['prob'].keys())
        pathways2 = set(cellchat2.netP['prob'].keys())
        common_pathways = pathways1 & pathways2
        
        # Compute correlation for common pathways
        similarities = {}
        for pathway in common_pathways:
            prob1 = cellchat1.netP['prob'][pathway].flatten()
            prob2 = cellchat2.netP['prob'][pathway].flatten()
            
            if len(prob1) == len(prob2):
                corr = np.corrcoef(prob1, prob2)[0, 1]
                similarities[pathway] = corr
        
        return {
            'type': 'functional',
            'pathway_similarity': similarities,
            'mean_similarity': np.mean(list(similarities.values())) if similarities else 0
        }
    
    else:  # structural
        # Compare network topology
        net1 = cellchat1.net['prob'].sum(axis=0)
        net2 = cellchat2.net['prob'].sum(axis=0)
        
        # Flatten and correlate
        corr = np.corrcoef(net1.flatten(), net2.flatten())[0, 1]
        
        return {
            'type': 'structural',
            'correlation': corr
        }       

# cellchat/visualization.py
"""
Visualization functions for CellChat
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, List, Tuple, Dict
import warnings

def plot_circle_network(
    cellchat_obj,
    sources_use: Optional[List[str]] = None,
    targets_use: Optional[List[str]] = None,
    signaling: Optional[List[str]] = None,
    remove_isolate: bool = True,
    thresh: float = 0.05,
    figsize: Tuple[int, int] = (8, 8),
    **kwargs
):
    """
    Plot circle plot of cell-cell communication network
    
    Parameters
    ----------
    cellchat_obj : CellChat
        CellChat object
    sources_use : Optional[List[str]]
        Source cell groups to show
    targets_use : Optional[List[str]]
        Target cell groups to show
    signaling : Optional[List[str]]
        Specific signaling pathways to show
    remove_isolate : bool
        Remove isolated nodes
    thresh : float
        Threshold for edge display
    """
    try:
        import circlize
    except ImportError:
        print("circlize not available, using alternative visualization")
        return plot_chord_diagram(cellchat_obj, signaling, thresh, figsize)
    
    # Get network data
    if signaling is not None:
        # Aggregate specified pathways
        prob_matrix = np.zeros_like(list(cellchat_obj.netP['prob'].values())[0])
        for pathway in signaling:
            if pathway in cellchat_obj.netP['prob']:
                prob_matrix += cellchat_obj.netP['prob'][pathway]
    else:
        prob_matrix = cellchat_obj.net['prob'].sum(axis=0)
    
    # Filter by threshold
    prob_matrix[prob_matrix < thresh] = 0
    
    # Filter sources and targets
    groups = cellchat_obj.unique_groups
    if sources_use is not None:
        source_idx = [i for i, g in enumerate(groups) if g in sources_use]
        prob_matrix = prob_matrix[source_idx, :]
        groups = [groups[i] for i in source_idx]
    
    if targets_use is not None:
        target_idx = [i for i, g in enumerate(groups) if g in targets_use]
        prob_matrix = prob_matrix[:, target_idx]
    
    # Remove isolated nodes
    if remove_isolate:
        active_nodes = (prob_matrix.sum(axis=0) > 0) | (prob_matrix.sum(axis=1) > 0)
        prob_matrix = prob_matrix[active_nodes][:, active_nodes]
        groups = [g for g, a in zip(groups, active_nodes) if a]
    
    # Create circular layout (simplified implementation)
    fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(projection='polar'))
    
    n_groups = len(groups)
    theta = np.linspace(0, 2*np.pi, n_groups, endpoint=False)
    
    # Plot nodes
    ax.scatter(theta, np.ones(n_groups), s=1000, alpha=0.6)
    
    # Add labels
    for i, (t, g) in enumerate(zip(theta, groups)):
        ax.text(t, 1.2, g, ha='center', va='center')
    
    # Plot edges (simplified - would use bezier curves in full implementation)
    for i in range(n_groups):
        for j in range(n_groups):
            if prob_matrix[i, j] > thresh:
                # Draw connection
                ax.plot([theta[i], theta[j]], [1, 1], 
                       alpha=prob_matrix[i, j]/prob_matrix.max(),
                       linewidth=2)
    
    ax.set_ylim(0, 1.5)
    ax.axis('off')
    plt.tight_layout()
    
    return fig, ax

def plot_chord_diagram(
    cellchat_obj,
    signaling: Optional[List[str]] = None,
    thresh: float = 0.05,
    figsize: Tuple[int, int] = (10, 10)
):
    """
    Plot chord diagram of communication network
    """
    from matplotlib.patches import Arc, Wedge
    from matplotlib.collections import PatchCollection
    
    # Get network data
    if signaling is not None:
        prob_matrix = np.zeros_like(list(cellchat_obj.netP['prob'].values())[0])
        for pathway in signaling:
            if pathway in cellchat_obj.netP['prob']:
                prob_matrix += cellchat_obj.netP['prob'][pathway]
    else:
        prob_matrix = cellchat_obj.net['prob'].sum(axis=0)
    
    prob_matrix[prob_matrix < thresh] = 0
    groups = cellchat_obj.unique_groups
    n_groups = len(groups)
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Calculate angles for each group
    gap = 0.02
    total_angle = 2 * np.pi * (1 - gap * n_groups)
    angles = []
    current_angle = 0
    
    for i in range(n_groups):
        group_size = total_angle / n_groups
        angles.append((current_angle, current_angle + group_size))
        current_angle += group_size + 2 * np.pi * gap
    
    # Draw group arcs
    colors = plt.cm.tab20(np.linspace(0, 1, n_groups))
    
    for i, ((start, end), color, group) in enumerate(zip(angles, colors, groups)):
        wedge = Wedge((0, 0), 1, np.degrees(start), np.degrees(end),
                     width=0.1, facecolor=color, edgecolor='white', linewidth=2)
        ax.add_patch(wedge)
        
        # Add label
        mid_angle = (start + end) / 2
        x = 1.2 * np.cos(mid_angle)
        y = 1.2 * np.sin(mid_angle)
        ax.text(x, y, group, ha='center', va='center', fontsize=10)
    
    # Draw ribbons (connections)
    for i in range(n_groups):
        for j in range(n_groups):
            if prob_matrix[i, j] > thresh:
                start_angle = (angles[i][0] + angles[i][1]) / 2
                end_angle = (angles[j][0] + angles[j][1]) / 2
                
                # Create bezier curve for ribbon
                t = np.linspace(0, 1, 100)
                # Simplified ribbon - would use proper bezier in full implementation
                x = (1 - t) * 0.9 * np.cos(start_angle) + t * 0.9 * np.cos(end_angle)
                y = (1 - t) * 0.9 * np.sin(start_angle) + t * 0.9 * np.sin(end_angle)
                
                alpha = prob_matrix[i, j] / prob_matrix.max()
                ax.plot(x, y, color=colors[i], alpha=alpha*0.5, linewidth=2)
    
    plt.tight_layout()
    return fig, ax

def plot_heatmap(
    cellchat_obj,
    signaling: Optional[List[str]] = None,
    sources_use: Optional[List[str]] = None,
    targets_use: Optional[List[str]] = None,
    thresh: float = 0.05,
    figsize: Tuple[int, int] = (8, 6),
    **kwargs
):
    """
    Plot heatmap of communication probability
    """
    # Get network data
    if signaling is not None:
        prob_matrix = np.zeros_like(list(cellchat_obj.netP['prob'].values())[0])
        for pathway in signaling:
            if pathway in cellchat_obj.netP['prob']:
                prob_matrix += cellchat_obj.netP['prob'][pathway]
    else:
        prob_matrix = cellchat_obj.net['prob'].sum(axis=0)
    
    # Create DataFrame for easier plotting
    groups = cellchat_obj.unique_groups
    df = pd.DataFrame(prob_matrix, index=groups, columns=groups)
    
    # Filter
    if sources_use is not None:
        df = df.loc[sources_use, :]
    if targets_use is not None:
        df = df.loc[:, targets_use]
    
    # Apply threshold
    df[df < thresh] = 0
    
    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(df, annot=True, fmt='.2f', cmap='Reds', 
                square=True, linewidths=0.5, ax=ax,
                cbar_kws={'label': 'Communication Probability'})
    
    ax.set_xlabel('Target')
    ax.set_ylabel('Source')
    ax.set_title('Cell-Cell Communication Heatmap')
    
    plt.tight_layout()
    return fig, ax

def plot_bubble(
    cellchat_obj,
    sources_use: Optional[List[str]] = None,
    targets_use: Optional[List[str]] = None,
    signaling: Optional[List[str]] = None,
    pairLR_use: Optional[List[str]] = None,
    thresh: float = 0.05,
    figsize: Tuple[int, int] = (12, 8),
    **kwargs
):
    """
    Plot bubble plot of L-R pairs
    """
    # Prepare data
    if signaling is not None:
        lr_pairs = cellchat_obj.LR[cellchat_obj.LR['pathway_name'].isin(signaling)]
    else:
        lr_pairs = cellchat_obj.LR
    
    if pairLR_use is not None:
        lr_pairs = lr_pairs[lr_pairs['interaction_name'].isin(pairLR_use)]
    
    # Get probabilities
    prob_data = []
    pval_data = []
    
    for idx in lr_pairs.index:
        lr_name = lr_pairs.loc[idx, 'interaction_name']
        prob = cellchat_obj.net['prob'][idx]
        pval = cellchat_obj.net['pval'][idx]
        
        for i, source in enumerate(cellchat_obj.unique_groups):
            for j, target in enumerate(cellchat_obj.unique_groups):
                if sources_use is None or source in sources_use:
                    if targets_use is None or target in targets_use:
                        if prob[i, j] > thresh:
                            prob_data.append({
                                'source': source,
                                'target': target,
                                'lr_pair': lr_name,
                                'prob': prob[i, j],
                                'pval': pval[i, j]
                            })
    
    df = pd.DataFrame(prob_data)
    
    if df.empty:
        print("No significant interactions to plot")
        return None, None
    
    # Create bubble plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create interaction labels
    df['interaction'] = df['source'] + ' -> ' + df['target']
    
    # Plot
    scatter = ax.scatter(
        range(len(df)),
        df['lr_pair'].astype('category').cat.codes,
        s=df['prob'] * 1000,
        c=-np.log10(df['pval'] + 1e-10),
        cmap='Reds',
        alpha=0.6,
        edgecolors='black',
        linewidth=0.5
    )
    
    # Labels
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df['interaction'], rotation=90)
    ax.set_yticks(range(len(df['lr_pair'].unique())))
    ax.set_yticklabels(df['lr_pair'].unique())
    
    plt.colorbar(scatter, label='-log10(p-value)', ax=ax)
    ax.set_xlabel('Cell-Cell Interactions')
    ax.set_ylabel('L-R Pairs')
    ax.set_title('Communication Probability (bubble size)')
    
    plt.tight_layout()
    return fig, ax

def plot_river(
    cellchat_list: List,
    pathway: str,
    thresh: float = 0.05,
    figsize: Tuple[int, int] = (12, 6)
):
    """
    Plot river (Sankey) plot showing information flow changes
    """
    import matplotlib.patches as mpatches
    from matplotlib.sankey import Sankey
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Simplified implementation - full version would use proper Sankey diagram
    n_conditions = len(cellchat_list)
    
    for i, cellchat in enumerate(cellchat_list):
        if pathway in cellchat.netP['prob']:
            prob = cellchat.netP['prob'][pathway]
            
            # Calculate total flow for each group
            outgoing = prob.sum(axis=1)
            incoming = prob.sum(axis=0)
            
            x_pos = i / (n_conditions - 1) if n_conditions > 1 else 0.5
            
            # Plot as bars
            ax.barh(range(len(outgoing)), outgoing, left=x_pos-0.05, 
                   height=0.05, alpha=0.6, label=f'Condition {i+1}')
    
    ax.set_xlabel('Condition')
    ax.set_ylabel('Cell Groups')
    ax.set_title(f'Information Flow: {pathway}')
    ax.legend()
    
    plt.tight_layout()
    return fig, ax

def plot_contribution(
    cellchat_obj,
    signaling: str,
    thresh: float = 0.05,
    figsize: Tuple[int, int] = (10, 6)
):
    """
    Plot contribution of each L-R pair to a pathway
    """
    # Get L-R pairs in this pathway
    lr_pairs = cellchat_obj.LR[cellchat_obj.LR['pathway_name'] == signaling]
    
    # Calculate contribution
    contributions = []
    
    for idx in lr_pairs.index:
        lr_name = lr_pairs.loc[idx, 'interaction_name']
        prob = cellchat_obj.net['prob'][idx].sum()
        contributions.append({'lr_pair': lr_name, 'contribution': prob})
    
    df = pd.DataFrame(contributions).sort_values('contribution', ascending=False)
    
    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    ax.barh(range(len(df)), df['contribution'])
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df['lr_pair'])
    ax.set_xlabel('Contribution')
    ax.set_title(f'L-R Pair Contribution to {signaling}')
    
    plt.tight_layout()
    return fig, ax

def plot_signaling_gene_expression(
    cellchat_obj,
    signaling: str,
    enriched_only: bool = True,
    figsize: Tuple[int, int] = (12, 8)
):
    """
    Plot expression of signaling genes
    """
    # Get genes in this pathway
    lr_pairs = cellchat_obj.LR[cellchat_obj.LR['pathway_name'] == signaling]
    
    genes = set()
    for _, row in lr_pairs.iterrows():
        genes.add(row['ligand'])
        genes.add(row['receptor'])
    
    # Get expression
    gene_expr = []
    for gene in genes:
        expr = cellchat_obj._get_gene_expression(gene)
        if expr is not None:
            gene_expr.append(pd.Series(expr, name=gene, 
                                      index=cellchat_obj.unique_groups))
    
    if not gene_expr:
        print("No gene expression data available")
        return None, None
    
    df = pd.DataFrame(gene_expr).T
    
    # Plot heatmap
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(df.T, cmap='RdYlBu_r', center=0, annot=True, fmt='.2f',
                cbar_kws={'label': 'Expression Level'}, ax=ax)
    
    ax.set_xlabel('Genes')
    ax.set_ylabel('Cell Groups')
    ax.set_title(f'Signaling Gene Expression: {signaling}')
    
    plt.tight_layout()
    return fig, ax      

# cellchat/comparison.py
"""
Comparison analysis between CellChat objects
"""

import numpy as np
import pandas as pd
from typing import List, Dict
from scipy.stats import mannwhitneyu

def compare_cellchat_objects(
    cellchat1,
    cellchat2,
    comparison_type: str = "functional"
) -> Dict:
    """
    Compare two CellChat objects
    """
    from .analysis import compute_network_similarity
    
    # Compute similarity
    similarity = compute_network_similarity(cellchat1, cellchat2, type=comparison_type)
    
    # Identify differential pathways
    diff_pathways = identify_differential_pathways([cellchat1, cellchat2])
    
    # Identify conserved pathways
    conserved = identify_conserved_pathways([cellchat1, cellchat2])
    
    results = {
        'similarity': similarity,
        'differential_pathways': diff_pathways,
        'conserved_pathways': conserved
    }
    
    return results

def identify_differential_pathways(
    cellchat_list: List,
    thresh: float = 0.05
) -> pd.DataFrame:
    """
    Identify differential pathways across conditions
    """
    # Get common pathways
    all_pathways = set(cellchat_list[0].netP['prob'].keys())
    for cellchat in cellchat_list[1:]:
        all_pathways &= set(cellchat.netP['prob'].keys())
    
    results = []
    
    for pathway in all_pathways:
        # Get pathway strengths
        strengths = []
        for cellchat in cellchat_list:
            strength = cellchat.netP['prob'][pathway].sum()
            strengths.append(strength)
        
        # Test for difference
        if len(strengths) == 2:
            # Simple fold change for two conditions
            fc = strengths[1] / (strengths[0] + 1e-10)
            log_fc = np.log2(fc)
            
            results.append({
                'pathway': pathway,
                'condition1_strength': strengths[0],
                'condition2_strength': strengths[1],
                'log2_fc': log_fc,
                'significant': abs(log_fc) > 1
            })
    
    return pd.DataFrame(results).sort_values('log2_fc', ascending=False)

def identify_conserved_pathways(
    cellchat_list: List,
    correlation_thresh: float = 0.7
) -> List[str]:
    """
    Identify conserved pathways across conditions
    """
    all_pathways = set(cellchat_list[0].netP['prob'].keys())
    for cellchat in cellchat_list[1:]:
        all_pathways &= set(cellchat.netP['prob'].keys())
    
    conserved = []
    
    for pathway in all_pathways:
        # Get pathway matrices
        matrices = [cc.netP['prob'][pathway].flatten() for cc in cellchat_list]
        
        # Compute pairwise correlations
        correlations = []
        for i in range(len(matrices)):
            for j in range(i+1, len(matrices)):
                corr = np.corrcoef(matrices[i], matrices[j])[0, 1]
                correlations.append(corr)
        
        # Check if conserved
        if all(c > correlation_thresh for c in correlations):
            conserved.append(pathway)
    
    return conserved

def rank_net(
    cellchat_list: List,
    mode: str = "comparison",
    stacked: bool = True
) -> Dict:
    """
    Rank networks by information flow
    """
    results = {}
    
    for i, cellchat in enumerate(cellchat_list):
        # Compute total information flow for each pathway
        pathway_flows = {}
        for pathway, prob in cellchat.netP['prob'].items():
            pathway_flows[pathway] = prob.sum()
        
        # Rank pathways
        ranked = sorted(pathway_flows.items(), key=lambda x: x[1], reverse=True)
        results[f'condition_{i+1}'] = ranked
    
    return results

# cellchat/spatial.py
"""
Spatial analysis functions for spatially resolved transcriptomics
"""

import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix, KDTree
from typing import Optional, Dict, List

def compute_spatial_distance(
    coords: pd.DataFrame,
    metric: str = 'euclidean'
) -> np.ndarray:
    """
    Compute pairwise spatial distances
    """
    from scipy.spatial.distance import pdist, squareform
    
    dist = pdist(coords.values, metric=metric)
    dist_matrix = squareform(dist)
    
    return dist_matrix

def identify_spatial_neighbors(
    coords: pd.DataFrame,
    radius: Optional[float] = None,
    k: Optional[int] = None
) -> Dict[int, List[int]]:
    """
    Identify spatial neighbors
    
    Parameters
    ----------
    coords : pd.DataFrame
        Spatial coordinates
    radius : Optional[float]
        Radius for neighbor search
    k : Optional[int]
        Number of nearest neighbors
    """
    tree = KDTree(coords.values)
    neighbors = {}
    
    if radius is not None:
        # Radius-based neighbors
        for i in range(len(coords)):
            indices = tree.query_ball_point(coords.iloc[i].values, radius)
            neighbors[i] = [idx for idx in indices if idx != i]
    
    elif k is not None:
        # K-nearest neighbors
        for i in range(len(coords)):
            distances, indices = tree.query(coords.iloc[i].values, k=k+1)
            neighbors[i] = indices[1:].tolist()  # Exclude self
    
    return neighbors

def compute_spatial_autocorrelation(
    cellchat_obj,
    pathway: str,
    method: str = 'moran'
) -> float:
    """
    Compute spatial autocorrelation of pathway activity
    """
    if not cellchat_obj.is_spatial:
        raise ValueError("CellChat object does not contain spatial information")
    
    if pathway not in cellchat_obj.netP['prob']:
        raise ValueError(f"Pathway {pathway} not found")
    
    # Get pathway activity for each cell
    prob = cellchat_obj.netP['prob'][pathway]
    
    # Simplified Moran's I calculation
    if method == 'moran':
        # Compute spatial weights
        dist = compute_spatial_distance(cellchat_obj.spatial_coords)
        weights = 1 / (dist + 1)  # Inverse distance weighting
        np.fill_diagonal(weights, 0)
        
        # Normalize weights
        weights = weights / weights.sum(axis=1, keepdims=True)
        
        # Compute activity
        activity = prob.sum(axis=1)  # Outgoing activity
        
        # Moran's I
        n = len(activity)
        mean_activity = activity.mean()
        
        numerator = 0
        denominator = 0
        
        for i in range(n):
            for j in range(n):
                numerator += weights[i, j] * (activity[i] - mean_activity) * (activity[j] - mean_activity)
            denominator += (activity[i] - mean_activity) ** 2
        
        morans_i = (n / weights.sum()) * (numerator / denominator)
        
        return morans_i
    
    return 0.0

def identify_spatial_patterns(
    cellchat_obj,
    n_patterns: int = 5
) -> Dict:
    """
    Identify spatial patterns of cell-cell communication
    """
    if not cellchat_obj.is_spatial:
        raise ValueError("CellChat object does not contain spatial information")
    
    from sklearn.decomposition import NMF
    
    # Aggregate all pathway activities
    all_activities = []
    pathway_names = []
    
    for pathway, prob in cellchat_obj.netP['prob'].items():
        activity = prob.sum(axis=1)  # Outgoing activity
        all_activities.append(activity)
        pathway_names.append(pathway)
    
    activity_matrix = np.array(all_activities).T
    
    # Apply NMF to identify patterns
    nmf = NMF(n_components=n_patterns, random_state=42)
    patterns = nmf.fit_transform(activity_matrix)
    pattern_composition = nmf.components_
    
    results = {
        'patterns': patterns,
        'pattern_composition': pattern_composition,
        'pathway_names': pathway_names
    }
    
    return results

# cellchat/utils.py
"""
Utility functions
"""

import numpy as np
import pandas as pd
from typing import Optional, Union

def create_cellchat_from_seurat(
    seurat_obj,
    group_by: str,
    assay: str = "RNA",
    slot: str = "data"
):
    """
    Create CellChat object from Seurat object (using rpy2)
    """
    try:
        from rpy2 import robjects as ro
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()
        
        # This is a placeholder - actual implementation would use rpy2
        # to interface with R Seurat object
        
        raise NotImplementedError("Seurat interface requires rpy2")
    
    except ImportError:
        raise ImportError("rpy2 required for Seurat interface")

def create_cellchat_from_scanpy(
    adata,
    group_by: str,
    use_raw: bool = False,
    spatial_key: Optional[str] = None
):
    """
    Create CellChat object from Scanpy AnnData object
    
    Parameters
    ----------
    adata : anndata.AnnData
        Scanpy AnnData object
    group_by : str
        Column in adata.obs for cell grouping
    use_raw : bool
        Use raw counts
    spatial_key : Optional[str]
        Key for spatial coordinates in adata.obsm
    """
    # Get expression data
    if use_raw and adata.raw is not None:
        expr = adata.raw.X.T
        gene_names = adata.raw.var_names.tolist()
    else:
        expr = adata.X.T
        gene_names = adata.var_names.tolist()
    
    # Create expression DataFrame
    expr_df = pd.DataFrame(
        expr,
        index=gene_names,
        columns=adata.obs_names
    )
    
    # Get metadata
    meta = adata.obs[[group_by]].copy()
    
    # Get spatial coordinates if available
    spatial_coords = None
    if spatial_key is not None and spatial_key in adata.obsm:
        spatial_coords = pd.DataFrame(
            adata.obsm[spatial_key],
            index=adata.obs_names
        )
    
    # Create CellChat object
    from .core import CellChat
    
    cellchat = CellChat(
        data=expr_df,
        meta=meta,
        group_by=group_by,
        spatial_coords=spatial_coords
    )
    
    return cellchat

def merge_cellchat_objects(
    cellchat_list: list,
    add_names: Optional[list] = None
):
    """
    Merge multiple CellChat objects for comparison
    """
    if add_names is None:
        add_names = [f"Condition_{i+1}" for i in range(len(cellchat_list))]
    
    merged = {
        'objects': cellchat_list,
        'names': add_names
    }
    
    return merged

def export_to_cytoscape(
    cellchat_obj,
    pathway: str,
    filename: str,
    thresh: float = 0.05
):
    """
    Export network to Cytoscape format
    """
    if pathway not in cellchat_obj.netP['prob']:
        raise ValueError(f"Pathway {pathway} not found")
    
    prob = cellchat_obj.netP['prob'][pathway]
    groups = cellchat_obj.unique_groups
    
    # Create edge list
    edges = []
    for i, source in enumerate(groups):
        for j, target in enumerate(groups):
            if prob[i, j] > thresh:
                edges.append({
                    'source': source,
                    'target': target,
                    'weight': prob[i, j],
                    'interaction': 'pp'  # protein-protein
                })
    
    edge_df = pd.DataFrame(edges)
    edge_df.to_csv(filename, index=False)
    
    print(f"Exported {len(edges)} edges to {filename}")

def save_cellchat(cellchat_obj, filename: str):
    """Save CellChat object to file"""
    import pickle
    
    with open(filename, 'wb') as f:
        pickle.dump(cellchat_obj, f)
    
    print(f"Saved CellChat object to {filename}")

def load_cellchat(filename: str):
    """Load CellChat object from file"""
    import pickle
    
    with open(filename, 'rb') as f:
        cellchat_obj = pickle.load(f)
    
    print(f"Loaded CellChat object from {filename}")
    return cellchat_obj

# examples/basic_usage.py
"""
Basic usage example of CellChat Python implementation
"""

import numpy as np
import pandas as pd
from scLucid.tools.pycellchat import CellChat, CellChatDB
from cellchat.preprocessing import preprocess_expression
from cellchat.visualization import *

# Example 1: Basic workflow
def basic_workflow():
    # Load data (example with random data)
    n_genes = 2000
    n_cells = 500
    
    expr_data = np.random.rand(n_genes, n_cells)
    gene_names = [f"Gene_{i}" for i in range(n_genes)]
    cell_names = [f"Cell_{i}" for i in range(n_cells)]
    
    expr_df = pd.DataFrame(expr_data, index=gene_names, columns=cell_names)
    
    # Create metadata
    cell_types = np.random.choice(['TypeA', 'TypeB', 'TypeC'], n_cells)
    meta = pd.DataFrame({'cell_type': cell_types}, index=cell_names)
    
    # Create CellChat object
    cellchat = CellChat(
        data=expr_df,
        meta=meta,
        group_by='cell_type'
    )
    
    # Set database
    db = CellChatDB(species='human')
    cellchat.set_database(db)
    
    # Preprocess
    cellchat.preprocess_data()
    
    # Identify overexpressed genes
    cellchat.identify_overexpressed_genes()
    
    # Compute communication probability
    cellchat.compute_communication_prob()
    
    # Filter communications
    cellchat.filter_communication()
    
    # Visualize
    plot_circle_network(cellchat)
    plot_heatmap(cellchat)
    
    return cellchat

# Example 2: Spatial transcriptomics
def spatial_workflow():
    # Load spatial data
    n_cells = 300
    
    # Generate random spatial coordinates
    coords = pd.DataFrame(
        np.random.rand(n_cells, 2) * 100,
        columns=['x', 'y']
    )
    
    # Generate expression data
    expr_data = np.random.rand(2000, n_cells)
    
    # Create metadata
    meta = pd.DataFrame({
        'cell_type': np.random.choice(['A', 'B', 'C'], n_cells)
    })
    
    # Create CellChat object with spatial info
    cellchat = CellChat(
        data=expr_data,
        meta=meta,
        group_by='cell_type',
        spatial_coords=coords
    )
    
    # Set database
    db = CellChatDB(species='human')
    cellchat.set_database(db)
    
    # Preprocess
    cellchat.preprocess_data()
    
    # Compute communication with spatial constraint
    cellchat.compute_communication_prob(distance_threshold=10.0)
    
    return cellchat

# Example 3: Comparison analysis
def comparison_workflow():
    # Create two CellChat objects
    cellchat1 = basic_workflow()
    cellchat2 = basic_workflow()
    
    # Compare
    from cellchat.comparison import compare_cellchat_objects
    
    results = compare_cellchat_objects(cellchat1, cellchat2)
    
    print("Comparison results:")
    print(results['similarity'])
    print(results['differential_pathways'])
    
    return results

if __name__ == "__main__":
    # Run basic workflow
    cellchat = basic_workflow()
    
    # Run spatial workflow
    # cellchat_spatial = spatial_workflow()
    
    # Run comparison
    # comparison_results = comparison_workflow()