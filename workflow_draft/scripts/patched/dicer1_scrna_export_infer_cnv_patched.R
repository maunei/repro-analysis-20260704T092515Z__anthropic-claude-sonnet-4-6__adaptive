# dicer1_scrna_export_infer_cnv_patched.R
# Patches:
#  - Removed reference to 'seurat_filter_tumor_indiv' (undefined variable — author bug)
#  - Seurat v5 slot → layer compatibility for GetAssayData
#  - Creates output/infer_cnv/ directory if missing

library(Seurat)
library(dplyr)

cat("Loading three Seurat objects...\n")
seuat_unfilter_normal <- readRDS(file = "../output/scrna/dicer1_scrna_unfilter_control.rds")
seurat_unfilter_tumor <- readRDS(file = "../output/scrna/dicer1_scrna_unfilter_tumor.rds")
seurat_filter_tumor   <- readRDS(file = "../output/scrna/dicer1_scrna_filter_tumor.rds")
cat("All three loaded.\n")

# Feature map: Ensembl ID ↔ gene symbol (from Cell Ranger output)
feature_map <- read.table(file = "../data/PX3085_TTAATACGCG-ACCCGAGGTG/PX3085_TTAATACGCG-ACCCGAGGTG_features.tsv.gz",
                          sep = "\t") %>%
  rename("ensembl_id" = "V1", "gene_name" = "V2", "type" = "V3")

# Seurat v5 slot → layer fallback helper
get_counts <- function(obj) {
  tryCatch(
    GetAssayData(object = obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(object = obj, assay = "RNA", slot = "counts")
  )
}

cat("Extracting counts matrices...\n")
counts_matrix_normal        <- get_counts(seuat_unfilter_normal)
counts_matrix_tumor         <- get_counts(seurat_unfilter_tumor)
counts_matrix_tumor_filter  <- get_counts(seurat_filter_tumor)

# Reassign rownames from gene symbols to Ensembl IDs for inferCNV
# feature_map row order matches the original Cell Ranger matrix row order
# The Seurat objects were created from the same Cell Ranger output, so the
# make.unique() gene symbols correspond 1-to-1 with feature_map rows
gene_order <- match(make.unique(feature_map$gene_name), rownames(counts_matrix_normal))
valid_rows <- !is.na(gene_order)

counts_matrix_normal[valid_rows, ]
rownames(counts_matrix_normal)[valid_rows] <- feature_map$ensembl_id[valid_rows]
rownames(counts_matrix_tumor)[valid_rows] <- feature_map$ensembl_id[valid_rows]
# filter_tumor has fewer cells but the same gene set
if (all(rownames(counts_matrix_tumor_filter) == rownames(counts_matrix_normal))) {
  rownames(counts_matrix_tumor_filter) <- rownames(counts_matrix_normal)
} else {
  # Direct replacement by position (genes are the same set, same order)
  rownames(counts_matrix_tumor_filter)[valid_rows] <- feature_map$ensembl_id[valid_rows]
}

counts_matrix_combined        <- cbind(counts_matrix_normal, counts_matrix_tumor)
counts_matrix_combined_filter <- cbind(counts_matrix_normal, counts_matrix_tumor_filter)

# Cell annotation tables
infer_cnv_cells_normal <- tibble(
  cell_name = gsub(x = colnames(seuat_unfilter_normal), pattern = "-", replacement = "."),
  cell_type = seuat_unfilter_normal$cell_type)

infer_cnv_cells_tumor <- tibble(
  cell_name = gsub(x = colnames(seurat_unfilter_tumor), pattern = "-", replacement = "."),
  cell_type = seurat_unfilter_tumor$cell_type)

infer_cnv_cells_tumor_filter <- tibble(
  cell_name = gsub(x = colnames(seurat_filter_tumor), pattern = "-", replacement = "."),
  cell_type = seurat_filter_tumor$cell_type)

infer_cnv_cells_normal_combine       <- infer_cnv_cells_normal %>% mutate(cell_type = paste0("Control ", cell_type))
infer_cnv_cells_tumor_combine        <- infer_cnv_cells_tumor  %>% mutate(cell_type = paste0("Tumor ", cell_type))
infer_cnv_cells_tumor_filter_combine <- infer_cnv_cells_tumor_filter %>% mutate(cell_type = paste0("Tumor ", cell_type))

infer_cnv_cells_combine        <- bind_rows(infer_cnv_cells_normal_combine, infer_cnv_cells_tumor_combine)
infer_cnv_cells_combine_filter <- bind_rows(infer_cnv_cells_normal_combine, infer_cnv_cells_tumor_filter_combine)

# Output directory
out_dir <- "../output/infer_cnv/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Seurat objects with inferCNV-compatible cell names (dashes → dots)
seuat_unfilter_normal_infercnv <- RenameCells(object = seuat_unfilter_normal,
                                              new.names = infer_cnv_cells_normal$cell_name)
seurat_filter_tumor_infercnv   <- RenameCells(object = seurat_filter_tumor,
                                              new.names = infer_cnv_cells_tumor_filter$cell_name)
seurat_filter_combined_infercnv <- merge(x = seurat_filter_tumor_infercnv,
                                         y = seuat_unfilter_normal_infercnv,
                                         merge.data = TRUE,
                                         merge.dr = TRUE)

cat("Saving combined infercnv Seurat object...\n")
saveRDS(object = seurat_filter_combined_infercnv,
        file = paste0(out_dir, "dicer1_scrna_filter_combined_infercnv.rds"))

cat("Writing count matrices...\n")
write.table(round(counts_matrix_normal, digits = 3),
            file = paste0(out_dir, "dicer1_scrna_unfilter_control_counts.matrix"),
            quote = FALSE, sep = "\t")

write.table(round(counts_matrix_tumor, digits = 3),
            file = paste0(out_dir, "dicer1_scrna_unfilter_tumor_counts.matrix"),
            quote = FALSE, sep = "\t")

write.table(round(counts_matrix_combined, digits = 3),
            file = paste0(out_dir, "dicer1_scrna_unfilter_combined_counts.matrix"),
            quote = FALSE, sep = "\t")

write.table(round(counts_matrix_combined_filter, digits = 3),
            file = paste0(out_dir, "dicer1_scrna_filter_combined_counts.matrix"),
            quote = FALSE, sep = "\t")

cat("Writing cell annotation tables...\n")
write.table(infer_cnv_cells_normal, file = paste0(out_dir, "dicer1_scrna_unfilter_control_infer_cnv_cells.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(infer_cnv_cells_tumor, file = paste0(out_dir, "dicer1_scrna_unfilter_tumor_infer_cnv_cells.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(infer_cnv_cells_combine, file = paste0(out_dir, "dicer1_scrna_unfilter_combined_infer_cnv_cells.tsv"),
            quote = FALSE, sep = FALSE, row.names = FALSE, col.names = FALSE)

write.table(infer_cnv_cells_combine_filter, file = paste0(out_dir, "dicer1_scrna_filter_combined_infer_cnv_cells.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

cat("inferCNV export complete. Output in:", out_dir, "\n")
