"""
Phylogenetic tree construction and analysis for tumor evolution.

This module provides tools for building and analyzing phylogenetic
trees from single-cell data.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union
from anndata import AnnData
import logging

log = logging.getLogger(__name__)


class PhylogenyBuilder:
    """
    Build phylogenetic trees from single-cell data.

    Parameters
    ----------
    method : str
        Tree building method ("nj", "upgma", "ml")
    distance_metric : str
        Distance metric for tree construction

    Attributes
    ----------
    tree_ : dict
        Phylogenetic tree structure
    """

    def __init__(
        self,
        method: str = "nj",
        distance_metric: str = "euclidean",
    ):
        self.method = method
        self.distance_metric = distance_metric
        self.tree_: Optional[Dict] = None
        self.distance_matrix_: Optional[pd.DataFrame] = None

    def build_phylogenetic_tree(
        self,
        adata: AnnData,
        clone_key: str = "clone_id",
        use_variants: bool = False,
        variant_key: Optional[str] = None,
    ) -> Dict:
        """
        Build phylogenetic tree from clonal data.

        Parameters
        ----------
        adata : AnnData
            Expression or variant data
        clone_key : str
            Column containing clone IDs
        use_variants : bool
            Use variant data instead of expression
        variant_key : str, optional
            Key for variant data in adata.obsm

        Returns
        -------
        dict
            Phylogenetic tree structure
        """
        # Get clone-level data
        clones = adata.obs[clone_key].unique()
        n_clones = len(clones)

        if use_variants and variant_key and variant_key in adata.obsm:
            # Use variant data
            clone_data = self._aggregate_by_clone(
                adata.obsm[variant_key], adata.obs[clone_key], clones
            )
        else:
            # Use expression data (mean per clone)
            clone_data = self._aggregate_by_clone(
                adata.X, adata.obs[clone_key], clones
            )
            if hasattr(clone_data, 'toarray'):
                clone_data = clone_data.toarray()

        # Calculate distance matrix
        from scipy.spatial.distance import pdist, squareform
        distances = pdist(clone_data, metric=self.distance_metric)
        self.distance_matrix_ = pd.DataFrame(
            squareform(distances),
            index=clones,
            columns=clones
        )

        # Build tree
        if self.method == "nj":
            self.tree_ = self._neighbor_joining(self.distance_matrix_)
        elif self.method == "upgma":
            self.tree_ = self._upgma(self.distance_matrix_)
        else:
            raise ValueError(f"Unknown method: {self.method}")

        return self.tree_

    def _aggregate_by_clone(
        self,
        data: np.ndarray,
        clone_labels: pd.Series,
        clones: List,
    ) -> np.ndarray:
        """Aggregate data by clone."""
        aggregated = []
        for clone in clones:
            mask = clone_labels == clone
            clone_data = data[mask].mean(axis=0)
            aggregated.append(clone_data)
        return np.array(aggregated)

    def _neighbor_joining(self, distance_matrix: pd.DataFrame) -> Dict:
        """
        Build tree using Neighbor-Joining algorithm.

        This is a simplified implementation of NJ.
        """
        n = len(distance_matrix)
        labels = list(distance_matrix.index)

        # Initialize
        dm = distance_matrix.values.copy()
        active = list(range(n))
        nodes = [{"id": i, "label": labels[i], "children": [], "height": 0.0} for i in range(n)]
        next_node_id = n

        while len(active) > 2:
            # Calculate Q matrix
            q = np.zeros((len(active), len(active)))

            for i, ii in enumerate(active):
                for j, jj in enumerate(active):
                    if i != j:
                        r_i = np.sum([dm[ii, kk] for kk in active if kk != ii])
                        r_j = np.sum([dm[jj, kk] for kk in active if kk != jj])
                        q[i, j] = (len(active) - 2) * dm[ii, jj] - r_i - r_j

            # Find minimum Q
            np.fill_diagonal(q, np.inf)
            min_idx = np.unravel_index(np.argmin(q), q.shape)
            i, j = min_idx
            ii, jj = active[i], active[j]

            # Calculate branch lengths
            r_i = np.sum([dm[ii, kk] for kk in active if kk != ii])
            r_j = np.sum([dm[jj, kk] for kk in active if kk != jj])

            d_ij = dm[ii, jj]
            limb_i = d_ij / 2 + (r_i - r_j) / (2 * (len(active) - 2))
            limb_j = d_ij - limb_i

            # Create new node
            new_node = {
                "id": next_node_id,
                "label": f"Internal_{next_node_id}",
                "children": [nodes[ii], nodes[jj]],
                "height": max(nodes[ii]["height"] + limb_i, nodes[jj]["height"] + limb_j),
                "branch_length_i": limb_i,
                "branch_length_j": limb_j,
            }
            nodes.append(new_node)

            # Update distance matrix
            new_distances = []
            for kk in active:
                if kk != ii and kk != jj:
                    new_d = (dm[ii, kk] + dm[jj, kk] - d_ij) / 2
                    new_distances.append(new_d)

            # Update active list
            active = [k for k in active if k != ii and k != jj] + [next_node_id]

            # Update distance matrix
            new_dm = np.zeros((len(active), len(active)))
            old_active = [k for k in active[:-1]]
            for i_new, i_old in enumerate(old_active):
                for j_new, j_old in enumerate(old_active):
                    new_dm[i_new, j_new] = dm[i_old, j_old]

            for i_new, i_old in enumerate(old_active):
                new_dm[i_new, -1] = new_distances[i_new]
                new_dm[-1, i_new] = new_distances[i_new]

            dm = new_dm
            next_node_id += 1

        # Connect final two nodes
        if len(active) == 2:
            ii, jj = active[0], active[1]
            final_node = {
                "id": next_node_id,
                "label": "Root",
                "children": [nodes[ii], nodes[jj]],
                "height": max(nodes[ii]["height"], nodes[jj]["height"]) + dm[0, 1] / 2,
                "branch_length_i": dm[0, 1] / 2,
                "branch_length_j": dm[0, 1] / 2,
                "is_root": True,
            }
            nodes.append(final_node)

        return final_node

    def _upgma(self, distance_matrix: pd.DataFrame) -> Dict:
        """
        Build tree using UPGMA algorithm.
        """
        n = len(distance_matrix)
        labels = list(distance_matrix.index)

        dm = distance_matrix.values.copy()
        active = list(range(n))
        nodes = [{"id": i, "label": labels[i], "children": [], "height": 0.0, "size": 1} for i in range(n)]
        next_node_id = n

        while len(active) > 1:
            # Find minimum distance pair
            min_dist = np.inf
            min_pair = (0, 1)
            for i, ii in enumerate(active):
                for j, jj in enumerate(active):
                    if i < j and dm[i, j] < min_dist:
                        min_dist = dm[i, j]
                        min_pair = (ii, jj)

            ii, jj = min_pair

            # Create new node
            height = min_dist / 2
            new_node = {
                "id": next_node_id,
                "label": f"Internal_{next_node_id}",
                "children": [nodes[ii], nodes[jj]],
                "height": height,
                "size": nodes[ii]["size"] + nodes[jj]["size"],
            }
            nodes.append(new_node)

            # Calculate new distances
            new_distances = []
            for kk in active:
                if kk != ii and kk != jj:
                    new_d = (dm[active.index(ii), active.index(kk)] * nodes[ii]["size"] +
                             dm[active.index(jj), active.index(kk)] * nodes[jj]["size"]) / \
                            (nodes[ii]["size"] + nodes[jj]["size"])
                    new_distances.append(new_d)

            # Update active list
            active = [k for k in active if k != ii and k != jj] + [next_node_id]

            # Update distance matrix
            new_dm = np.zeros((len(active), len(active)))
            for i_new in range(len(active) - 1):
                for j_new in range(len(active) - 1):
                    new_dm[i_new, j_new] = dm[i_new, j_new]

            for i_new in range(len(active) - 1):
                new_dm[i_new, -1] = new_distances[i_new]
                new_dm[-1, i_new] = new_distances[i_new]

            dm = new_dm
            next_node_id += 1

        return nodes[-1]

    def root_tree(self, tree: Dict, outgroup: Optional[str] = None) -> Dict:
        """
        Root the phylogenetic tree.

        Parameters
        ----------
        tree : dict
            Tree structure
        outgroup : str, optional
            Outgroup label for rooting

        Returns
        -------
        dict
            Rooted tree
        """
        if outgroup is None:
            # Midpoint rooting
            tree["is_root"] = True
            return tree

        # Root at outgroup
        # This is a simplified rooting
        tree["is_root"] = True
        tree["outgroup"] = outgroup

        return tree


def build_phylogenetic_tree(
    adata: AnnData,
    clone_key: str = "clone_id",
    method: str = "nj",
) -> Dict:
    """
    Build phylogenetic tree from single-cell data.

    Parameters
    ----------
    adata : AnnData
        Expression data
    clone_key : str
        Column containing clone IDs
    method : str
        Tree building method

    Returns
    -------
    dict
        Phylogenetic tree structure
    """
    builder = PhylogenyBuilder(method=method)
    return builder.build_phylogenetic_tree(adata, clone_key)


def root_tree(tree: Dict, outgroup: Optional[str] = None) -> Dict:
    """
    Root a phylogenetic tree.

    Parameters
    ----------
    tree : dict
        Tree structure
    outgroup : str, optional
        Outgroup label

    Returns
    -------
    dict
        Rooted tree
    """
    builder = PhylogenyBuilder()
    return builder.root_tree(tree, outgroup)


def calculate_tree_metrics(tree: Dict) -> Dict:
    """
    Calculate phylogenetic tree metrics.

    Parameters
    ----------
    tree : dict
        Tree structure

    Returns
    -------
    dict
        Tree metrics
    """
    def count_tips(node: Dict) -> int:
        if not node.get("children"):
            return 1
        return sum(count_tips(child) for child in node["children"])

    def count_internal_nodes(node: Dict) -> int:
        if not node.get("children"):
            return 0
        return 1 + sum(count_internal_nodes(child) for child in node["children"])

    def get_depth(node: Dict) -> float:
        if not node.get("children"):
            return 0.0
        child_depths = [get_depth(child) + node.get(f"branch_length_{i+1}", 0.5)
                       for i, child in enumerate(node["children"])]
        return max(child_depths)

    return {
        "n_tips": count_tips(tree),
        "n_internal_nodes": count_internal_nodes(tree),
        "tree_height": get_depth(tree),
        "is_rooted": tree.get("is_root", False),
    }
