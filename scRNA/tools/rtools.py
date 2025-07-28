"""
A bridge for running R-based single-cell analysis tools.

This module uses rpy2 to create a robust interface between Python/AnnData and
popular R/Bioconductor packages like CopyKAT, Monocle3, and CellChat.
"""

import logging
from typing import Literal, Optional

import anndata
import numpy as np
import pandas as pd
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr

log = logging.getLogger(__name__)

# Activate R-to-Python data conversions
pandas2ri.activate()
numpy2ri.activate()


class RTools:
    """A class to manage and run R-based analysis tools via rpy2."""

    def __init__(self):
        """Initializes the R environment and checks for required packages."""
        log.info("Initializing R environment via rpy2...")
        try:
            self.R = importr("base")
            self.SCE = importr("SingleCellExperiment")
            log.info("R environment successfully initialized.")
        except ImportError as e:
            raise ImportError(
                f"rpy2 initialization failed. Is R installed and in your PATH? Error: {e}"
            )
        # Add CellChat to the list of required packages
        self._check_r_packages(
            ["Seurat", "SingleCellExperiment", "CopyKAT", "monocle3", "CellChat"]
        )

    def _check_r_packages(self, packages: list):
        """Checks if a list of R packages are installed."""
        for pkg in packages:
            try:
                importr(pkg)
                log.info(f"R package '{pkg}' found.")
            except ImportError:
                log.warning(f"R package '{pkg}' not found.")
                raise ImportError(
                    f"Required R package '{pkg}' is not installed. "
                    f"Please install it in R (e.g., using BiocManager::install('{pkg}'))."
                )

    def _anndata_to_sce(
        self,
        adata: anndata.AnnData,
        use_raw: bool = True,
        use_rep: Optional[str] = None,
    ):
        """Converts an AnnData object to an R SingleCellExperiment object."""
        # ... (this helper function remains the same as before) ...
        log.info("Converting AnnData to R SingleCellExperiment...")
        from scipy.sparse import csc_matrix, issparse

        # Prepare assays
        assays = {}
        if use_raw and adata.raw is not None:
            counts_matrix = adata.raw.X.T
            gene_names = adata.raw.var_names
            assays["counts"] = (
                csc_matrix(counts_matrix) if issparse(counts_matrix) else counts_matrix
            )
        else:  # Use adata.X as counts if no raw
            counts_matrix = adata.X.T
            gene_names = adata.var_names
            assays["counts"] = (
                csc_matrix(counts_matrix) if issparse(counts_matrix) else counts_matrix
            )

        # Add logcounts if available
        if "log1p_norm" in adata.layers:
            logcounts_matrix = adata.layers["log1p_norm"].T
            assays["logcounts"] = (
                csc_matrix(logcounts_matrix)
                if issparse(logcounts_matrix)
                else logcounts_matrix
            )

        sce = self.SCE.SingleCellExperiment(
            assays=self.R.list(**assays),
            rowData=pandas2ri.py2rpy(pd.DataFrame(index=gene_names)),
            colData=pandas2ri.py2rpy(adata.obs),
        )

        # Add embeddings
        if use_rep and use_rep in adata.obsm:
            self.R.reducedDim(sce, use_rep.upper(), adata.obsm[use_rep])

        return sce

    def _sce_to_anndata(self, sce: "rpy2.robjects.RObject") -> anndata.AnnData:
        """Converts an R SingleCellExperiment object back to an AnnData object."""
        log.info("Converting R SingleCellExperiment to AnnData...")

        # Extract assays
        assays = self.R("assays")(sce)
        assay_dict = {}
        for name in self.R("names")(assays):
            matrix = self.R("assay")(sce, name)
            if "dgCMatrix" in self.R("class")(matrix):  # Check for sparse matrix
                matrix = matrix.T
            assay_dict[name] = matrix

        # Extract obs and var
        obs_df = pandas2ri.rpy2py(self.R("colData")(sce))
        var_df = pandas2ri.rpy2py(self.R("rowData")(sce))

        # Create AnnData object
        adata = anndata.AnnData(
            X=assay_dict.pop("counts", None), obs=obs_df, var=var_df, layers=assay_dict
        )

        # Extract embeddings
        reduced_dims = self.R("reducedDims")(sce)
        for name in self.R("names")(reduced_dims):
            adata.obsm[f"X_{name.lower()}"] = self.R("reducedDim")(sce, name)

        return adata

    def run_r_script(self, script: str, args: dict):
        """Executes an R script with specified arguments."""
        # ... (this helper function remains the same as before) ...
        try:
            r_func = self.R(script)
            result = r_func(**args)
            return result
        except Exception as e:
            r_error_message = str(e)
            log.error(
                f"An error occurred while executing the R script: {r_error_message}"
            )
            raise RuntimeError(f"R execution failed: {r_error_message}")

    def anndata_to_seurat(
        self, adata: anndata.AnnData, use_raw: bool = True
    ) -> "rpy2.robjects.RObject":
        """
        Converts an AnnData object to an R Seurat object.

        This is achieved by first converting to SingleCellExperiment, then to Seurat.

        Args:
            adata: The AnnData object to convert.
            use_raw: Whether to use adata.raw for the counts matrix.

        Returns:
            An rpy2 object representing the Seurat object in the R environment.
        """
        log.info("Converting AnnData to Seurat object...")
        sce = self._anndata_to_sce(adata, use_raw=use_raw)

        r_script = """
        function(sce) {
            library(Seurat)
            seurat_obj <- as.Seurat(sce, counts = "counts", data = "logcounts")
            return(seurat_obj)
        }
        """
        seurat_obj = self.run_r_script(r_script, {"sce": sce})
        log.info("Conversion to Seurat object complete.")
        return seurat_obj

    def seurat_to_anndata(self, seurat_obj: "rpy2.robjects.RObject") -> anndata.AnnData:
        """
        Converts an R Seurat object back to an AnnData object.

        This is achieved by first converting to SingleCellExperiment, then to AnnData.

        Args:
            seurat_obj: An rpy2 object representing a Seurat object.

        Returns:
            A new AnnData object.
        """
        log.info("Converting Seurat object back to AnnData...")
        r_script = """
        function(seurat_obj) {
            library(Seurat)
            sce <- as.SingleCellExperiment(seurat_obj)
            return(sce)
        }
        """
        sce = self.run_r_script(r_script, {"seurat_obj": seurat_obj})

        # Now convert SCE to AnnData
        adata = self._sce_to_anndata(sce)
        log.info("Conversion to AnnData complete.")
        return adata

    def run_copykat(
        self, adata: anndata.AnnData, key_added: str = "copykat_prediction", **kwargs
    ):
        """Runs the CopyKAT algorithm to classify tumor/normal cells."""
        # ... (this function remains the same as before) ...
        log.info("Running CopyKAT analysis via R...")
        sce = self._anndata_to_sce(adata, use_raw=True)

        r_script = """
        function(sce) {
            library(CopyKAT)
            counts <- as.matrix(assay(sce, "counts"))
            copykat.test <- copykat(rawmat=counts, ngene.chr=5, sam.name="test")
            pred.test <- data.frame(copykat.test$prediction)
            return(pred.test)
        }
        """
        predictions_df = self.run_r_script(r_script, {"sce": sce})
        adata.obs[key_added] = predictions_df.loc[adata.obs_names, "copykat.pred"]
        log.info(f"CopyKAT results added to adata.obs['{key_added}']")
        return adata

    def run_monocle3(
        self,
        adata: anndata.AnnData,
        root_group_key: str,
        root_group_name: str,
        key_added: str = "monocle3_pseudotime",
        **kwargs,
    ):
        """Runs the Monocle3 workflow to calculate pseudotime."""
        # ... (this function remains the same as before) ...
        log.info("Running Monocle3 workflow via R...")
        sce = self._anndata_to_sce(adata, use_raw=False, use_rep="X_umap")

        r_script = """
        function(sce, root_group_key, root_group_name) {
            library(monocle3)
            cds <- as.cell_data_set(sce)
            cds <- cluster_cells(cds, reduction_method = "UMAP")
            cds <- learn_graph(cds)
            cell_ids <- which(colData(cds)[, root_group_key] == root_group_name)
            closest_vertex <- cds@principal_graph_aux[['UMAP']]$pr_graph_cell_proj_closest_vertex
            root_pr_nodes <- unique(closest_vertex[colnames(cds)[cell_ids], 1])
            if (length(root_pr_nodes) == 0) {
                stop("Could not find a valid root node for the specified root group.")
            }
            cds <- order_cells(cds, root_pr_nodes=root_pr_nodes)
            return(pseudotime(cds))
        }
        """
        pseudotime_values = self.run_r_script(
            r_script,
            {
                "sce": sce,
                "root_group_key": root_group_key,
                "root_group_name": root_group_name,
            },
        )

        pseudotime_series = pd.Series(
            np.array(pseudotime_values), index=adata.obs_names
        )
        adata.obs[key_added] = pseudotime_series
        log.info(f"Monocle3 pseudotime added to adata.obs['{key_added}']")
        return adata

    def run_cellchat(
        self,
        adata: anndata.AnnData,
        groupby: str,
        species: Literal["human", "mouse"],
        key_added: str = "cellchat",
    ):
        """
        Run CellChat cell-cell communication analysis via R.

        Args:
            adata: AnnData object with normalized data in .X or .raw.X.
            groupby: Column in adata.obs for cell type annotation.
            species: Species of the data ('human' or 'mouse').
            key_added: Key in adata.uns to store CellChat results.

        Returns:
            AnnData object with CellChat results in adata.uns[key_added].
        """
        log.info(f"Running CellChat analysis for {species} via R...")

        # CellChat needs normalized data, prefer raw if available
        sce = self._anndata_to_sce(adata, use_raw=True)

        r_script = """
        function(sce, groupby, species) {
            library(CellChat)
            
            # Extract data and metadata
            data.input <- assay(sce, "counts")
            meta <- as.data.frame(colData(sce))
            
            # Create CellChat object
            cellchat <- createCellChat(object = data.input, meta = meta, group.by = groupby)
            
            # Set the species-specific ligand-receptor database
            if (species == "human") {
                CellChatDB <- CellChatDB.human
            } else if (species == "mouse") {
                CellChatDB <- CellChatDB.mouse
            } else {
                stop("Species must be 'human' or 'mouse'")
            }
            cellchat@DB <- CellChatDB
            
            # Preprocessing
            cellchat <- subsetData(cellchat)
            cellchat <- identifyOverExpressedGenes(cellchat)
            cellchat <- identifyOverExpressedInteractions(cellchat)
            
            # Infer communication network
            cellchat <- computeCommunProb(cellchat)
            cellchat <- filterCommunication(cellchat, min.cells = 10)
            cellchat <- computeCommunProbPathway(cellchat)
            cellchat <- aggregateNet(cellchat)
            
            # Extract the final communication dataframe
            df.net <- subsetCommunication(cellchat)
            
            return(df.net)
        }
        """

        interactions_df = self.run_r_script(
            r_script, {"sce": sce, "groupby": groupby, "species": species}
        )

        # Store results
        adata.uns[key_added] = {
            "interactions": interactions_df,
            "params": {"groupby": groupby, "species": species},
        }
        log.info(
            f"CellChat analysis complete. Found {len(interactions_df)} significant interactions."
        )
        log.info(f"Results stored in adata.uns['{key_added}']['interactions']")

        return adata
