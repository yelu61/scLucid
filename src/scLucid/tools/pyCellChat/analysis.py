"""
Analysis functions for CellChat (R-free)
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import NMF
from typing import Dict, List, Tuple, Optional
import logging

log = logging.getLogger(__name__)


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
    if n_groups < k:
        k = n_groups
        log.warning(f"Reduced k to {k} due to small number of groups")

    if k == 0:
        return {
            'patterns': np.zeros((n_groups, 1)),
            'pathway_patterns': np.zeros((1, n_pathways)),
            'dominant_pattern': np.zeros(n_groups, dtype=int),
            'feature_matrix': features,
            'pathways': pathways
        }

    nmf = NMF(n_components=k, init='nndsvda', random_state=42, max_iter=500)
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
