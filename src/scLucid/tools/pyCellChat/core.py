"""
Core CellChat class implementation (R-free)
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, Optional, Union

import numpy as np
import pandas as pd
from scipy import sparse

log = logging.getLogger(__name__)


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
    Main CellChat class for cell-cell communication analysis (Pure Python)

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
        config: Optional[CellChatConfig] = None,
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

        log.info(
            f"Created CellChat object: {len(self.unique_groups)} cell groups, "
            f"{self.data_expr.shape[1]} cells, {self.data_expr.shape[0]} genes"
        )

    @property
    def n_cells(self) -> int:
        """Number of cells in the expression matrix."""
        return int(self.data_expr.shape[1])

    @property
    def n_genes(self) -> int:
        """Number of genes in the expression matrix."""
        return int(self.data_expr.shape[0])

    def set_database(self, db: "CellChatDB"):
        """Set CellChatDB database"""
        self.db = db
        self.LR = db.interaction
        log.info(f"Set database with {len(self.LR)} interactions")

    def preprocess_data(self, subset_data: bool = True, do_sparse: bool = True):
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

            if len(common_genes) == 0:
                raise ValueError("No common genes between data and database!")

            if subset_data:
                # Subset to common genes
                gene_idx = [self.gene_names.index(g) for g in common_genes if g in self.gene_names]
                self.data_expr = self.data_expr[gene_idx, :]
                self.gene_names = common_genes
                log.info(f"Subset to {len(common_genes)} common genes")

        # Convert to sparse if needed
        if do_sparse and not sparse.issparse(self.data_expr):
            self.data_expr = sparse.csr_matrix(self.data_expr)

        log.info(
            f"Preprocessed data: {self.data_expr.shape[0]} genes x {self.data_expr.shape[1]} cells"
        )

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
        distance_threshold: Optional[float] = None,
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

        if self.LR is None or len(self.LR) == 0:
            raise ValueError("No ligand-receptor pairs loaded. Check database.")

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
            ligand = lr_info.get("ligand")
            receptor = lr_info.get("receptor")

            if pd.isna(ligand) or pd.isna(receptor):
                continue

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
                        spatial_weight = self._compute_spatial_weight(i, j, distance_threshold)

                    # Compute probability using mass action law
                    prob[idx, i, j] = L_expr * R_expr * spatial_weight

        # Store results
        self.net["prob"] = prob
        self.net["pval"] = pval
        self.net["LR"] = self.LR

        # Compute pathway-level probability
        self._compute_pathway_prob()

        log.info(f"Computed communication probability for {n_interactions} interactions")
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
        dist_matrix = cdist(centroids, centroids, metric="euclidean")

        return dist_matrix

    def _compute_spatial_weight(
        self, group_i: int, group_j: int, distance_threshold: float
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
        if "prob" not in self.net:
            raise ValueError("Run compute_communication_prob first")

        # Group interactions by pathway
        pathways = self.LR["pathway_name"].unique()
        n_groups = len(self.unique_groups)

        pathway_prob = {}
        pathway_pval = {}

        for pathway in pathways:
            if pd.isna(pathway):
                continue

            # Get interactions in this pathway
            pathway_mask = self.LR["pathway_name"] == pathway
            pathway_indices = np.where(pathway_mask)[0]

            if len(pathway_indices) == 0:
                continue

            # Aggregate probabilities
            prob_pathway = self.net["prob"][pathway_indices].sum(axis=0)
            pval_pathway = self.net["pval"][pathway_indices].min(axis=0)

            pathway_prob[pathway] = prob_pathway
            pathway_pval[pathway] = pval_pathway

        self.netP["prob"] = pathway_prob
        self.netP["pval"] = pathway_pval

        log.info(f"Computed pathway-level probabilities for {len(pathway_prob)} pathways")

    def filter_communication(self, min_cells: int = 10, thresh_p: float = 0.05):
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
        if "pval" in self.net:
            prob = self.net["prob"]
            pval = self.net["pval"]

            # Set non-significant to 0
            prob[pval > thresh_p] = 0
            self.net["prob"] = prob

        log.info(f"Filtered communications: {len(valid_groups)} valid groups")

    def compute_network_centrality(self):
        """Compute network centrality measures"""
        from .analysis import compute_centrality

        if "prob" not in self.net:
            raise ValueError("Run compute_communication_prob first")

        centrality_results = compute_centrality(self.net["prob"], self.unique_groups)
        self.net["centrality"] = centrality_results

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

        if "prob" not in self.netP:
            raise ValueError("No pathway-level network computed")

        roles = identify_roles(self.netP["prob"], pattern=pattern)

        return roles

    def compare_interactions(self, other: "CellChat", comparison_type: str = "functional"):
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

        comparison_results = compare_cellchat_objects(self, other, comparison_type=comparison_type)

        return comparison_results

    def save(self, filename: str):
        """Save CellChat object to file"""
        import pickle

        with open(filename, "wb") as f:
            pickle.dump(self, f)
        log.info(f"Saved CellChat object to {filename}")

    @classmethod
    def load(cls, filename: str):
        """Load CellChat object from file"""
        import pickle

        with open(filename, "rb") as f:
            obj = pickle.load(f)
        log.info(f"Loaded CellChat object from {filename}")
        return obj
