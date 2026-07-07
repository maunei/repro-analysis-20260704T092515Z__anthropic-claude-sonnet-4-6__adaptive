# Script to export the required inputs for inferCNV

# Required packages
library(Seurat)
library(dplyr)

# Read in the data and feature annotations
seuat_unfilter_normal <- readRDS(file = "../output/scrna/dicer1_scrna_unfilter_control.rds")
seurat_unfilter_tumor <- readRDS(file = "../output/scrna/dicer1_scrna_unfilter_tumor.rds")
seurat_filter_tumor <- readRDS(file = "../output/scrna/dicer1_scrna_filter_tumor.rds")

# The features (genes) are the same for each sample; therefore we could use any for this purpose
feature_map <- read.table(file = "../data/PX3085_TTAATACGCG-ACCCGAGGTG/PX3085_TTAATACGCG-ACCCGAGGTG_features.tsv.gz",
                          sep = "\t") %>%
  rename("ensembl_id" = "V1",
         "gene_name" = "V2",
         "type" = "V3")

# Extract required fields for inferCNV
counts_matrix_normal = GetAssayData(object = seuat_unfilter_normal,
                                    slot = "counts",
                                    assay = "RNA")

counts_matrix_tumor = GetAssayData(object = seurat_unfilter_tumor,
                                   slot = "counts",
                                   assay = "RNA")

counts_matrix_tumor_filter = GetAssayData(object = seurat_filter_tumor,
                                          slot = "counts",
                                          assay = "RNA")

counts_matrix_tumor_filter_indiv <- GetAssayData(object = seurat_filter_tumor_indiv,
                                                 slot = "counts",
                                                 assay = "RNA")

rownames(counts_matrix_normal) <- feature_map$ensembl_id
rownames(counts_matrix_tumor) <- feature_map$ensembl_id
rownames(counts_matrix_tumor_filter) <- feature_map$ensembl_id

counts_matrix_combined <- cbind(counts_matrix_normal, counts_matrix_tumor)
counts_matrix_combined_filter <- cbind(counts_matrix_normal, counts_matrix_tumor_filter)

infer_cnv_cells_normal <- tibble(cell_name = gsub(x = colnames(seuat_unfilter_normal),
                                                  pattern = "-",
                                                  replacement = "."),
                                 cell_type = seuat_unfilter_normal$cell_type)

infer_cnv_cells_tumor <- tibble(cell_name = gsub(x = colnames(seurat_unfilter_tumor),
                                                 pattern = "-",
                                                 replacement = "."),
                                cell_type = seurat_unfilter_tumor$cell_type)

infer_cnv_cells_tumor_filter <- tibble(cell_name = gsub(x = colnames(seurat_filter_tumor),
                                                        pattern = "-",
                                                        replacement = "."),
                                       cell_type = seurat_filter_tumor$cell_type)

infer_cnv_cells_normal_combine <- infer_cnv_cells_normal %>%
  mutate(cell_type = paste0("Control ", cell_type))

infer_cnv_cells_tumor_combine <- infer_cnv_cells_tumor %>%
  mutate(cell_type = paste0("Tumor ", cell_type))

infer_cnv_cells_tumor_filter_combine <- infer_cnv_cells_tumor_filter %>%
  mutate(cell_type = paste0("Tumor ", cell_type))

infer_cnv_cells_combine <- bind_rows(infer_cnv_cells_normal_combine,
                                     infer_cnv_cells_tumor_combine)

infer_cnv_cells_combine_filter <- bind_rows(infer_cnv_cells_normal_combine,
                                            infer_cnv_cells_tumor_filter_combine)

# Write to file
out_dir <- "../output/infer_cnv/"

# RDS objects with the same cell names as inferCNV requires, so that the results can be added back easily after running inferCNV and visualized
seuat_unfilter_normal_infercnv <- RenameCells(object = seuat_unfilter_normal,
                                              new.names = infer_cnv_cells_normal$cell_name)
seurat_filter_tumor_infercnv <- RenameCells(object = seurat_filter_tumor,
                                            new.names = infer_cnv_cells_tumor_filter$cell_name)
seurat_filter_combined_infercnv <- merge(x = seurat_filter_tumor_infercnv,
                                         y = seuat_unfilter_normal_infercnv,
                                         merge.data = TRUE,
                                         merge.dr = TRUE)

saveRDS(object = seurat_filter_combined_infercnv,
        file = paste0(out_dir, "dicer1_scrna_filter_combined_infercnv.rds"))

# Tables and matrices for input to inferCNV; various options to run inferCNV in various ways
write.table(round(counts_matrix_normal,
                  digits = 3),
            file = paste0(out_dir, "dicer1_scrna_unfilter_control_counts.matrix"),
            quote = FALSE,
            sep = "\t")

write.table(round(counts_matrix_tumor,
                  digits = 3),
            file = paste0(out_dir, "dicer1_scrna_unfilter_tumor_counts.matrix"),
            quote = FALSE,
            sep = "\t")

write.table(round(counts_matrix_combined,
                  digits = 3),
            file = paste0(out_dir, "dicer1_scrna_unfilter_combined_counts.matrix"),
            quote = FALSE,
            sep = "\t")

write.table(round(counts_matrix_combined_filter,
                  digits = 3),
            file = paste0(out_dir, "dicer1_scrna_filter_combined_counts.matrix"),
            quote = FALSE,
            sep = "\t")

write.table(infer_cnv_cells_normal,
            file = paste0(out_dir, "dicer1_scrna_unfilter_control_infer_cnv_cells.tsv"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

write.table(infer_cnv_cells_tumor,
            file = paste0(out_dir, "dicer1_scrna_unfilter_tumor_infer_cnv_cells.tsv"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

write.table(infer_cnv_cells_combine,
            file = paste0(out_dir, "dicer1_scrna_unfilter_combined_infer_cnv_cells.tsv"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

write.table(infer_cnv_cells_combine_filter,
            file = paste0(out_dir, "dicer1_scrna_filter_combined_infer_cnv_cells.tsv"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
