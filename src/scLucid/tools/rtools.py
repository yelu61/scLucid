"""
RTools: Python ↔ R bridge for single-cell advanced analysis.

Provides convenient, robust access to R-based packages such as
Monocle3, Slingshot, CopyKAT, CellChat, BayesPrism, DWLS, DESeq2, etc.

All conversions and method calls are handled automatically.
"""

import logging
from typing import Optional, Literal

import anndata
import numpy as np
import pandas as pd
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr

# Activate automatic conversion between Python and R objects
pandas2ri.activate()
numpy2ri.activate()

log = logging.getLogger(__name__)

class RTools:
    """
    Main interface for running R-based single-cell analysis from Python.
    """

    def __init__(self, verbose: bool = True):
        """
        Initialize the R environment and check for required packages.
        """
        log.info("Initializing RTools R environment")
        try:
            self.R = importr("base")
            self.SCE = importr("SingleCellExperiment")
        except Exception as e:
            raise ImportError(f"Failed to import core R packages: {e}")

        if verbose:
            self._check_r_packages([
                "Seurat", "SingleCellExperiment", "monocle3", "slingshot",
                "CopyKAT", "CellChat", "BayesPrism", "DWLS", "DESeq2"
            ])

    def _check_r_packages(self, pkgs):
        """
        Check if a list of R packages are installed.
        """
        for pkg in pkgs:
            try:
                importr(pkg)
                log.info(f"R package '{pkg}' found")
            except Exception:
                log.warning(f"R package '{pkg}' not found. Please install in R.")

    def _anndata_to_sce(
        self,
        adata: anndata.AnnData,
        use_raw: bool = True,
        use_rep: Optional[str] = None,
    ):
        """
        Convert AnnData (Python) to SingleCellExperiment (R).
        Optionally, add reducedDims if available.
        """
        from scipy.sparse import csc_matrix, issparse

        log.info("Converting AnnData to SingleCellExperiment")
        assays = {}
        # Use raw or X for counts
        if use_raw and adata.raw is not None:
            counts_matrix = adata.raw.X.T
            gene_names = adata.raw.var_names
            assays["counts"] = csc_matrix(counts_matrix) if issparse(counts_matrix) else counts_matrix
        else:
            counts_matrix = adata.X.T
            gene_names = adata.var_names
            assays["counts"] = csc_matrix(counts_matrix) if issparse(counts_matrix) else counts_matrix

        # Optionally add logcounts layer
        if "log1p_norm" in adata.layers:
            logcounts_matrix = adata.layers["log1p_norm"].T
            assays["logcounts"] = csc_matrix(logcounts_matrix) if issparse(logcounts_matrix) else logcounts_matrix

        # Create the SingleCellExperiment object in R
        sce = self.SCE.SingleCellExperiment(
            assays=self.R.list(**assays),
            rowData=pandas2ri.py2rpy(pd.DataFrame(index=gene_names)),
            colData=pandas2ri.py2rpy(adata.obs),
        )

        # Optionally add reducedDims (e.g. UMAP)
        if use_rep and use_rep in adata.obsm:
            try:
                self.R.reducedDim(sce, use_rep.upper(), adata.obsm[use_rep])
            except Exception:
                log.warning(f"Failed to add reducedDim({use_rep}). Skipped.")

        return sce

    def _sce_to_anndata(self, sce):
        """
        Convert SingleCellExperiment (R) back to AnnData (Python).
        """
        log.info("Converting SingleCellExperiment to AnnData")
        assays = self.R("assays")(sce)
        assay_dict = {}
        for name in self.R("names")(assays):
            matrix = self.R("assay")(sce, name)
            assay_dict[name] = matrix
        obs_df = pandas2ri.rpy2py(self.R("colData")(sce))
        var_df = pandas2ri.rpy2py(self.R("rowData")(sce))
        adata = anndata.AnnData(
            X=assay_dict.pop("counts", None), obs=obs_df, var=var_df, layers=assay_dict
        )
        # Try to convert reducedDims
        try:
            reduced_dims = self.R("reducedDims")(sce)
            for name in self.R("names")(reduced_dims):
                adata.obsm[f"X_{name.lower()}"] = self.R("reducedDim")(sce, name)
        except Exception:
            pass
        return adata

    def run_r_script(self, script: str, args: dict):
        """
        Run a raw R function using the given script and arguments.
        Returns the result as an R object.
        """
        try:
            r_func = self.R(script)
            result = r_func(**args)
            return result
        except Exception as e:
            msg = str(e)
            log.error(f"Error executing R script: {msg}")
            raise RuntimeError(f"R execution error: {msg}")

    def run_monocle3(
        self,
        adata: anndata.AnnData,
        root_group_key: str,
        root_group_name: str,
        key_added: str = "monocle3_pseudotime",
        **kwargs,
    ):
        """
        Run Monocle3 trajectory inference in R, return pseudotime to adata.obs.

        Args:
            adata: AnnData object (must have UMAP).
            root_group_key: Column in adata.obs for root group.
            root_group_name: Name of group to use as origin.
            key_added: obs column to add pseudotime.

        Returns:
            AnnData with pseudotime in obs.
        """
        log.info("Running Monocle3 trajectory inference in R")
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
        log.info(f"Monocle3 pseudotime written to obs['{key_added}']")
        return adata

    def run_slingshot(
        self,
        adata: anndata.AnnData,
        groupby: str,
        start: Optional[str] = None,
        key_added: str = "slingshot_pseudotime",
    ) -> anndata.AnnData:
        """
        Run Slingshot trajectory inference in R, return pseudotime(s) to adata.obs.

        Args:
            adata: AnnData object (should have UMAP in obsm).
            groupby: Column in adata.obs (clusters).
            start: Cluster name to use as root (optional).
            key_added: obs column prefix for pseudotime.

        Returns:
            AnnData with pseudotime(s) in obs.
        """
        log.info(f"Running Slingshot: groupby={groupby}, start={start}")
        sce = self._anndata_to_sce(adata, use_raw=False, use_rep="X_umap")
        r_script = """
        function(sce, groupby, start) {
            library(slingshot)
            cluster <- as.factor(colData(sce)[[groupby]])
            reduced <- reducedDims(sce)[["UMAP"]]
            if (!is.null(start) && start != "") {
                ss <- slingshot(reduced, clusterLabels=cluster, start.clus=start)
            } else {
                ss <- slingshot(reduced, clusterLabels=cluster)
            }
            pt <- as.data.frame(slingPseudotime(ss))
            rownames(pt) <- colnames(sce)
            return(pt)
        }
        """
        pt_df = self.run_r_script(
            r_script, {"sce": sce, "groupby": groupby, "start": start if start else ""}
        )
        pt_df = pandas2ri.rpy2py(pt_df)
        # Write all lineages to obs (e.g. slingshot_pseudotime_1, _2, ...)
        for idx, col in enumerate(pt_df.columns):
            adata.obs[f"{key_added}_{idx+1}"] = pt_df[col].values
        # Default: first lineage as main pseudotime
        adata.obs[key_added] = pt_df.iloc[:, 0].values
        log.info(f"Slingshot pseudotime(s) written to obs['{key_added}']")
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
    
    def run_dwls_deconvolution(
        self,
        adata_ref: anndata.AnnData,
        bulk_data: pd.DataFrame,
        cell_type_key: str,
    ) -> pd.DataFrame:
        """Performs bulk data deconvolution using the DWLS R package."""
        log.info("Running DWLS deconvolution via R...")

        sce_ref = self._anndata_to_sce(adata_ref, use_raw=True)
        
        r_script_bulk = """
        function(bulk_df) {
            return(as.matrix(bulk_df))
        }
        """
        bulk_matrix_r = self.run_r_script(r_script_bulk, {"bulk_df": bulk_data})
        
        r_script_dwls = """
        function(sce, bulk.matrix, cell_type_key) {
            library(DWLS)
            library(SingleCellExperiment)
            
            # Build signature matrix from reference
            signature <- buildSignatureMatrixMAST(
                scdata = counts(sce),
                sc_cluster_info = colData(sce)[[cell_type_key]],
                subject_IDs = NULL # Assuming no sample grouping within reference
            )
            
            # Run deconvolution
            res <- trimm(x = signature, y = bulk.matrix)
            res_proportions <- res$relative.proportion
            
            return(as.data.frame(res_proportions))
        }
        """
        proportions_df = self.run_r_script(
            r_script_dwls,
            {
                "sce": sce_ref,
                "bulk.matrix": bulk_matrix_r,
                "cell_type_key": cell_type_key,
            },
        )
        log.info("DWLS deconvolution complete.")
        return proportions_df

    def run_bayesprism_deconvolution(
        self,
        adata_ref: anndata.AnnData,
        bulk_data: pd.DataFrame,
        cell_type_key: str,
        sample_key: str,
    ) -> pd.DataFrame:
        """Performs bulk data deconvolution using the BayesPrism R package."""
        log.info("Running BayesPrism deconvolution via R...")
        
        # BayesPrism has specific data format requirements
        # 1. Bulk data (matrix)
        r_script_bulk = "function(df) { as.matrix(df) }"
        bulk_matrix_r = self.run_r_script(r_script_bulk, {"df": bulk_data})
        
        # 2. scRNA reference data (matrix)
        sc_matrix = adata_ref.raw.X.T if adata_ref.raw is not None else adata_ref.X.T
        if scipy.sparse.issparse(sc_matrix):
            sc_matrix = sc_matrix.toarray()
        sc_matrix_r = numpy2ri.py2rpy(sc_matrix)
        self.R.rownames(sc_matrix_r, adata_ref.var_names)
        self.R.colnames(sc_matrix_r, adata_ref.obs_names)

        # 3. Cell type and sample labels
        cell_type_labels_r = pandas2ri.py2rpy(adata_ref.obs[cell_type_key])
        sample_labels_r = pandas2ri.py2rpy(adata_ref.obs[sample_key])
        
        r_script_bayesprism = """
        function(bk.dat, sc.dat, cell.type.labels, cell.state.labels) {
            library(BayesPrism)

            # Clean up gene names for compatibility
            sc.dat <- cleanup.genes(sc.dat, species="hs")
            bk.dat <- cleanup.genes(bk.dat, species="hs")
            
            # Align data
            prism.data <- new.prism(
                reference=sc.dat,
                mixture=bk.dat,
                input.type="count",
                cell.type.labels = cell.type.labels, 
                cell.state.labels = cell.state.labels,
                key=NULL # No tumor-specific analysis for now
            )
            
            # Run the deconvolution
            res <- run.prism(prism=prism.data, n.cores=4) # Use 4 cores for speed
            
            # Return the estimated proportions
            theta <- get.prism(res, "theta")
            return(as.data.frame(theta))
        }
        """
        proportions_df = self.run_r_script(
            r_script_bayesprism,
            {
                "bk.dat": bulk_matrix_r,
                "sc.dat": sc_matrix_r,
                "cell.type.labels": cell_type_labels_r,
                "cell.state.labels": sample_labels_r
            }
        )
        log.info("BayesPrism deconvolution complete.")
        return proportions_df
