# Script to add the results from inferCNV back to the seurat metadata for visualization

# Required packages
library(Seurat)

# Read in the data
seurat_unfilter_normal <- readRDS(file = "../output/scrna/dicer1_scrna_unfilter_control.rds")
seurat_filter_tumor <- readRDS(file = "../output/scrna/dicer1_scrna_filter_tumor.rds")
seurat_combined_filter_infercnv <- readRDS(file = "../output/infer_cnv/dicer1_scrna_filter_combined_infercnv_metadata.rds")

# Rename cells as named for inferCNV
seurat_unfilter_normal <- RenameCells(object = seurat_unfilter_normal,
                                      new.names = gsub(x = colnames(seurat_unfilter_normal),
                                                       pattern = "-",
                                                       replacement = "."))
seurat_filter_tumor <- RenameCells(object = seurat_filter_tumor,
                                   new.names = gsub(x = colnames(seurat_filter_tumor),
                                                    pattern = "-",
                                                    replacement = "."))

# Map the results per-cell back to our tumor-only object for plotting
cell_index <- match(colnames(seurat_filter_tumor), colnames(seurat_combined_filter_infercnv))
diff_meta <- setdiff(colnames(seurat_combined_filter_infercnv@meta.data),
                     colnames(seurat_filter_tumor@meta.data))
metadata <- seurat_combined_filter_infercnv@meta.data[cell_index, diff_meta]

stopifnot(all(rownames(metadata) == colnames(seurat_filter_tumor))) # All cells must match
seurat_filter_tumor_infercnv <- AddMetaData(object = seurat_filter_tumor,
                                            metadata = metadata)

# Visualize key features of interest
global_theme <- theme(plot.title = element_blank(),
                      axis.line = element_line(size = 0.3))

chr_plot <- FeaturePlot(object = seurat_filter_tumor_infercnv,
                        features = c("proportion_dupli_chr1",
                                     "proportion_dupli_chr6",
                                     "proportion_dupli_chr15",
                                     "proportion_dupli_chr2"),
                        order = TRUE,
                        keep.scale = "all",
                        ncol = 4) +
  global_theme

chr_plot_split_sample <- FeaturePlot(object = seurat_filter_tumor_infercnv,
                                     features = c("proportion_dupli_chr1",
                                                  "proportion_dupli_chr6",
                                                  "proportion_dupli_chr15",
                                                  "proportion_dupli_chr2"),
                                     order = TRUE,
                                     keep.scale = "all",
                                     ncol = 4,
                                     split.by = "Name",
                                     by.col = FALSE)

# Save plots
out_dir <- "../output/infer_cnv/"

ggsave(filename = "dicer1_scrna_filter_tumor_infercnv.svg", plot = chr_plot, path = out_dir, width = 24, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_infercnv_split_sample.svg", plot = chr_plot_split_sample, path = out_dir, width = 24, height = 16, device = "svg")
