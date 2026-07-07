# Missing link in inferCNV pipeline.
# Reads the completed inferCNV object from combined_filter/, adds per-cell
# HMM chromosome duplication metadata to the merged Seurat object, and
# saves as dicer1_scrna_filter_combined_infercnv_metadata.rds — which is
# what dicer1_scrna_add_infer_cnv_seurat.R expects.
#
# Run order:
#   1. dicer1_scrna_export_infer_cnv.R          (done)
#   2. dicer1_scrna_run_infer_cnv_v2_patched.R  (resume from step 14)
#   3. THIS SCRIPT
#   4. dicer1_scrna_add_infer_cnv_seurat.R      (produces final SVG plots)

library(infercnv)
library(Seurat)

out_dir <- "../output/infer_cnv/"

cat("Reading merged Seurat object...\n")
seurat_combined <- readRDS(file = paste0(out_dir, "dicer1_scrna_filter_combined_infercnv.rds"))

cat("Adding inferCNV HMM predictions to Seurat metadata...\n")
seurat_combined <- infercnv::add_to_seurat(
    seurat_obj           = seurat_combined,
    infercnv_output_path = paste0(out_dir, "combined_filter")
)

cat("Columns added:\n")
new_cols <- grep("proportion_dupli|cnv", colnames(seurat_combined@meta.data), value = TRUE)
print(new_cols)

cat("Saving metadata-enriched Seurat object...\n")
saveRDS(seurat_combined,
        file = paste0(out_dir, "dicer1_scrna_filter_combined_infercnv_metadata.rds"))

cat("Done. Saved: dicer1_scrna_filter_combined_infercnv_metadata.rds\n")
