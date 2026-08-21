#!/usr/bin/env Rscript

# Export the official sc_mixology SingleCellExperiment objects into a small,
# dependency-light interchange bundle. The Python companion writes the H5AD.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: export_mixology_rdata.R <sincell_with_class.RData> <output_dir>")
}

source_path <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("Matrix", quietly = TRUE) ||
    !requireNamespace("SingleCellExperiment", quietly = TRUE) ||
    !requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  stop("Matrix, SingleCellExperiment, and SummarizedExperiment are required")
}

source_env <- new.env(parent = emptyenv())
load(source_path, envir = source_env)
object_names <- c("sce_sc_10x_qc", "sce_sc_CELseq2_qc", "sce_sc_Dropseq_qc")
missing_objects <- object_names[!vapply(
  object_names,
  exists,
  logical(1),
  envir = source_env,
  inherits = FALSE
)]
if (length(missing_objects) > 0L) {
  stop(sprintf("Missing expected objects: %s", paste(missing_objects, collapse = ", ")))
}

objects <- list(
  `10x` = get("sce_sc_10x_qc", envir = source_env, inherits = FALSE),
  CELseq2 = get("sce_sc_CELseq2_qc", envir = source_env, inherits = FALSE),
  Dropseq = get("sce_sc_Dropseq_qc", envir = source_env, inherits = FALSE)
)
genes <- rownames(objects[[1]])
genes <- genes[vapply(genes, function(gene) {
  all(vapply(objects[-1], function(sce) gene %in% rownames(sce), logical(1)))
}, logical(1))]
if (length(genes) == 0L || anyDuplicated(genes)) {
  stop("Common gene intersection must be non-empty and unique")
}

count_matrices <- list()
metadata_rows <- list()
for (protocol in names(objects)) {
  sce <- objects[[protocol]]
  counts <- SummarizedExperiment::assay(sce[genes, ], "counts")
  counts <- methods::as(counts, "dgCMatrix")
  cell_id <- paste(protocol, colnames(sce), sep = "::")
  colnames(counts) <- cell_id
  count_matrices[[protocol]] <- counts

  metadata <- as.data.frame(SummarizedExperiment::colData(sce), stringsAsFactors = FALSE)
  required <- c("cell_line_demuxlet", "cell_line", "demuxlet_cls")
  missing_metadata <- setdiff(required, colnames(metadata))
  if (length(missing_metadata) > 0L) {
    stop(sprintf("%s lacks metadata: %s", protocol, paste(missing_metadata, collapse = ", ")))
  }
  optional <- intersect(
    c("outliers", "mapped_to_MT", "number_of_genes", "total_count_per_cell", "non_mt_percent"),
    colnames(metadata)
  )
  metadata <- metadata[, c(required, optional), drop = FALSE]
  metadata$cell_id <- cell_id
  metadata$protocol <- protocol
  metadata_rows[[protocol]] <- metadata
}

combined_counts <- do.call(cbind, count_matrices)
combined_metadata <- do.call(rbind, metadata_rows)
rownames(combined_metadata) <- NULL
if (!identical(colnames(combined_counts), combined_metadata$cell_id)) {
  stop("Count columns and metadata rows are not aligned")
}

invisible(Matrix::writeMM(combined_counts, file.path(output_dir, "counts.mtx")))
write.table(
  data.frame(gene = genes),
  file.path(output_dir, "genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  combined_metadata,
  file.path(output_dir, "cells.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

cat(sprintf(
  "Exported %d genes x %d cells across %d protocols\n",
  nrow(combined_counts), ncol(combined_counts), length(objects)
))
