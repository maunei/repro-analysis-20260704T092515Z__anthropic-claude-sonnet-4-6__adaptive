# Script to add the results from inferCNV back to the seurat metadata for visualization
# PATCHED: added library(ggplot2) for theme(); wrapped FeaturePlot in tryCatch for Seurat v5 S7 compat

library(Seurat)
library(ggplot2)

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

stopifnot(all(rownames(metadata) == colnames(seurat_filter_tumor)))
seurat_filter_tumor_infercnv <- AddMetaData(object = seurat_filter_tumor,
                                            metadata = metadata)

cat("CNV metadata columns available:\n")
cnv_cols <- grep("proportion_dupli|proportion_loss", names(seurat_filter_tumor_infercnv@meta.data), value = TRUE)
print(cnv_cols)

# Visualize key features of interest
global_theme <- theme(plot.title = element_blank(),
                      axis.line = element_line(linewidth = 0.3))  # size -> linewidth in ggplot2 >= 3.4

features_chr <- c("proportion_dupli_chr1",
                  "proportion_dupli_chr6",
                  "proportion_dupli_chr15",
                  "proportion_dupli_chr2")

# Fallback: manual UMAP scatter if FeaturePlot S7 incompatibility triggers
safe_featureplot <- function(obj, features, ...) {
  tryCatch(
    FeaturePlot(object = obj, features = features, ...),
    error = function(e) {
      cat("FeaturePlot error:", conditionMessage(e), "\n")
      cat("Falling back to manual ggplot scatter\n")
      umap_coords <- as.data.frame(Embeddings(obj, "umap"))
      colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
      plots <- lapply(features, function(feat) {
        if (feat %in% names(obj@meta.data)) {
          umap_coords[[feat]] <- obj@meta.data[[feat]]
        } else {
          umap_coords[[feat]] <- FetchData(obj, vars = feat)[, 1]
        }
        ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = .data[[feat]])) +
          geom_point(size = 0.3) +
          scale_color_gradient(low = "lightgrey", high = "blue") +
          ggtitle(feat) +
          theme_classic(base_size = 10)
      })
      patchwork::wrap_plots(plots, ncol = length(features))
    }
  )
}

chr_plot <- safe_featureplot(seurat_filter_tumor_infercnv,
                             features = features_chr,
                             order = TRUE,
                             keep.scale = "all",
                             ncol = 4)

chr_plot_split_sample <- safe_featureplot(seurat_filter_tumor_infercnv,
                                          features = features_chr,
                                          order = TRUE,
                                          keep.scale = "all",
                                          ncol = 4,
                                          split.by = "Name",
                                          by.col = FALSE)

# Save plots
out_dir <- "../output/infer_cnv/"

ggsave(filename = "dicer1_scrna_filter_tumor_infercnv.svg",
       plot = chr_plot, path = out_dir, width = 24, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_infercnv.png",
       plot = chr_plot, path = out_dir, width = 24, height = 4, dpi = 150)

ggsave(filename = "dicer1_scrna_filter_tumor_infercnv_split_sample.svg",
       plot = chr_plot_split_sample, path = out_dir, width = 24, height = 16, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_infercnv_split_sample.png",
       plot = chr_plot_split_sample, path = out_dir, width = 24, height = 16, dpi = 150)

cat("Step 4 complete. Plots saved to", out_dir, "\n")
