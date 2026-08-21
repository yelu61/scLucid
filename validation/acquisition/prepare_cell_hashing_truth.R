#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: prepare_cell_hashing_truth.R HTO_COUNTS RNA_COUNTS OUTPUT_TSV")
}

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

hto_path <- args[[1]]
rna_path <- args[[2]]
output_path <- args[[3]]

hto_table <- read.csv(gzfile(hto_path), row.names = 1, check.names = FALSE)
hto_counts <- as.matrix(hto_table[, seq_len(12), drop = FALSE])
rna_barcodes <- strsplit(readLines(gzfile(rna_path), n = 1), "\t", fixed = TRUE)[[1]]
joint_barcodes <- intersect(rna_barcodes, rownames(hto_counts))
if (length(joint_barcodes) < 2) {
  stop("fewer than two RNA/HTO barcodes overlap")
}

dummy <- sparseMatrix(
  i = rep(1L, length(joint_barcodes)),
  j = seq_along(joint_barcodes),
  x = rep(1, length(joint_barcodes)),
  dims = c(1L, length(joint_barcodes)),
  dimnames = list("barcode_placeholder", joint_barcodes)
)
object <- CreateSeuratObject(counts = dummy, min.cells = 0, min.features = 0)
object[["HTO"]] <- CreateAssayObject(counts = t(hto_counts[joint_barcodes, , drop = FALSE]))
object <- NormalizeData(
  object,
  assay = "HTO",
  normalization.method = "CLR",
  verbose = FALSE
)
object <- HTODemux(
  object,
  assay = "HTO",
  positive.quantile = 0.99,
  seed = 42,
  verbose = FALSE
)

metadata <- object[[]]
result <- data.frame(
  barcode = rownames(metadata),
  hto_max_id = metadata$HTO_maxID,
  hto_second_id = metadata$HTO_secondID,
  hto_margin = metadata$HTO_margin,
  hto_classification = metadata$HTO_classification,
  hto_classification_global = metadata$HTO_classification.global,
  hash_id = metadata$hash.ID,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.table(result, output_path, sep = "\t", quote = FALSE, row.names = FALSE)

