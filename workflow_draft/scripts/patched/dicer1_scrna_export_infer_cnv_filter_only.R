# dicer1_scrna_export_infer_cnv_filter_only.R
# Memory-optimized version: writes ONLY filter_combined matrix + annotation TSV
# (what run_infer_cnv.R actually needs - skips unfilter matrices that OOM)

library(Seurat)
library(dplyr)
library(Matrix)

out_dir <- "../output/infer_cnv/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading control (unfilter) object...\n")
seuat_unfilter_normal <- readRDS(file = "../output/scrna/dicer1_scrna_unfilter_control.rds")
cat("Control cells:", ncol(seuat_unfilter_normal), "\n")

cat("Loading filter_tumor object...\n")
seurat_filter_tumor <- readRDS(file = "../output/scrna/dicer1_scrna_filter_tumor.rds")
cat("Tumor cells:", ncol(seurat_filter_tumor), "\n")

# Feature map for Ensembl ID rename
feature_map <- read.table(
  file = "../data/PX3085_TTAATACGCG-ACCCGAGGTG/PX3085_TTAATACGCG-ACCCGAGGTG_features.tsv.gz",
  sep = "\t") %>% rename("ensembl_id" = "V1", "gene_name" = "V2", "type" = "V3")

# Extract sparse counts (no dense conversion)
get_counts_sparse <- function(obj) {
  tryCatch(
    GetAssayData(object = obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(object = obj, assay = "RNA", slot = "counts")
  )
}

cat("Extracting sparse count matrices...\n")
counts_normal <- get_counts_sparse(seuat_unfilter_normal)
counts_filter_tumor <- get_counts_sparse(seurat_filter_tumor)

# Rename rownames: gene symbols → Ensembl IDs
gene_order_idx <- match(make.unique(feature_map$gene_name), rownames(counts_normal))
valid_rows <- !is.na(gene_order_idx)
cat("Valid gene mappings:", sum(valid_rows), "\n")

rownames(counts_normal)[valid_rows] <- feature_map$ensembl_id[valid_rows]
rownames(counts_filter_tumor)[valid_rows] <- feature_map$ensembl_id[valid_rows]

# Cell name formatting (dashes → dots for inferCNV)
colnames(counts_normal) <- gsub("-", ".", colnames(counts_normal))
colnames(counts_filter_tumor) <- gsub("-", ".", colnames(counts_filter_tumor))

# Combine: normal (control) + filter_tumor
cat("Combining control + filter_tumor matrices (sparse)...\n")
counts_combined_filter <- cbind(counts_normal, counts_filter_tumor)
cat("Combined matrix dims:", nrow(counts_combined_filter), "x", ncol(counts_combined_filter), "\n")

# Write combined filter matrix (tab-delimited dense text, as inferCNV requires)
# For memory efficiency, write one chunk at a time using R connections
cat("Writing filter_combined counts matrix (dense text format)...\n")
cat("Expected file size: ~", round(nrow(counts_combined_filter) * ncol(counts_combined_filter) * 5 / 1e9, 2), "GB\n")

# Write using file connection in chunks to avoid huge single allocation
out_file <- paste0(out_dir, "dicer1_scrna_filter_combined_counts.matrix")
n_genes <- nrow(counts_combined_filter)
chunk_size <- 1000  # genes per chunk

# Write header (cell names)
con <- file(out_file, open = "wt")
writeLines(paste(colnames(counts_combined_filter), collapse = "\t"), con)
close(con)

con <- file(out_file, open = "at")
for (start in seq(1, n_genes, by = chunk_size)) {
  end <- min(start + chunk_size - 1, n_genes)
  chunk <- as.matrix(counts_combined_filter[start:end, ])
  chunk <- round(chunk, digits = 3)
  write.table(chunk, con, sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
  if (start %% 5000 == 1) cat("Progress:", start, "/", n_genes, "genes\n")
}
close(con)
rm(counts_combined_filter); gc()

cat("Filter combined matrix written.\n")

# Cell annotation files
infer_cnv_cells_normal <- data.frame(
  cell_name = colnames(counts_normal),
  cell_type = paste0("Control ", seuat_unfilter_normal$cell_type)
)
infer_cnv_cells_tumor_filter <- data.frame(
  cell_name = colnames(counts_filter_tumor),
  cell_type = paste0("Tumor ", seurat_filter_tumor$cell_type)
)
infer_cnv_cells_combine_filter <- rbind(infer_cnv_cells_normal, infer_cnv_cells_tumor_filter)

write.table(infer_cnv_cells_combine_filter,
            file = paste0(out_dir, "dicer1_scrna_filter_combined_infer_cnv_cells.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

cat("Cell annotation table written.\n")
cat("inferCNV filter-combined export complete.\n")
