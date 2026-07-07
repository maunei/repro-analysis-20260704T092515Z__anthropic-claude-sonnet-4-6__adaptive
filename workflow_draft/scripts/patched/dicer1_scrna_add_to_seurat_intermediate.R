# Intermediate script: runs infercnv::add_to_seurat() to create
# dicer1_scrna_filter_combined_infercnv_metadata.rds
# This file is required by dicer1_scrna_add_infer_cnv_seurat.R
# Run from: workflow_draft/repos/dicer1-cell-of-origin/code/

library(Seurat)
library(infercnv)

infercnv_output_dir <- "../output/infer_cnv/combined_filter/"
cat("InferCNV output dir:", infercnv_output_dir, "\n")
cat("Files:\n"); print(list.files(infercnv_output_dir))

# Load the Seurat objects for inferCNV cells
cat("Loading Seurat objects...\n")
seurat_unfilter_normal <- readRDS("../output/scrna/dicer1_scrna_unfilter_control.rds")
seurat_filter_tumor    <- readRDS("../output/scrna/dicer1_scrna_filter_tumor.rds")
gc()

# Rename cells: inferCNV replaces "-" with "."
seurat_unfilter_normal <- RenameCells(seurat_unfilter_normal,
  new.names = gsub("-", ".", colnames(seurat_unfilter_normal)))
seurat_filter_tumor    <- RenameCells(seurat_filter_tumor,
  new.names = gsub("-", ".", colnames(seurat_filter_tumor)))

# Merge into combined object (all cells used in inferCNV)
cat("Merging Seurat objects...\n")
seurat_combined <- merge(seurat_unfilter_normal, seurat_filter_tumor)
rm(seurat_unfilter_normal, seurat_filter_tumor); gc()

# Add inferCNV HMM predictions to metadata
cat("Running add_to_seurat...\n")
seurat_combined <- infercnv::add_to_seurat(
  seurat_obj          = seurat_combined,
  infercnv_output_path = infercnv_output_dir,
  top_n               = 10  # top 10 genes per segment for marker-gene features
)

cat("Metadata columns added:\n")
print(names(seurat_combined@meta.data))

# Save as the expected input for dicer1_scrna_add_infer_cnv_seurat.R
saveRDS(seurat_combined, "../output/infer_cnv/dicer1_scrna_filter_combined_infercnv_metadata.rds")
cat("Saved dicer1_scrna_filter_combined_infercnv_metadata.rds\n")
