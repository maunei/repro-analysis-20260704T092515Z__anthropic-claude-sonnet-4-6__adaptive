# dicer1_xenium_proseg_mus_musculus_patched.R
# Patches applied:
#   - DESeq2 + SeuratDisk library() wrapped in tryCatch (not installed in 653_renv2; not used in script)
#   - UpdateSeuratObject() added after readRDS (FOV class S4 slot migration for older RDS)
#   - Two-pass memory management: process WT first (load+subset+rm+gc), reload for MUT
#     avoids OOM on 15GB system with inferCNV running (~5GB already used)
#   - QC filter + subset combined in one step (colnames-based index) to minimize peak RAM
#   - DimHeatmap/ElbowPlot wrapped in tryCatch (display-device calls fail in batch mode)
#   - Path bug: ./output/xenium/ -> ../output/xenium/ (misplaced saveRDS/readRDS blocks)
#   - Removed dead xenium.obj.mut save before it is created (author copy-paste artifact)
#   - SCT_snn_res.0.5 replaced with seurat_clusters (object reclustered at 0.6 by that point)
#   - dir.create() for output dirs added
#   - future.globals.maxSize set globally

library(Seurat)
library(tidyverse)
library(BPCells)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
tryCatch(library(DESeq2),    error = function(e) cat("DESeq2 not installed (not used):", conditionMessage(e), "\n"))
library(SeuratWrappers)
library(monocle3)
library(reshape2)
library(MetBrewer)
library(RColorBrewer)
library(ggplot2)
library(viridis)
library(glmGamPoi)
library(future)
library(reshape2)
library(patchwork)
library(cowplot)
library(ggpubr)
library(yaml)
library(arrow)
tryCatch(library(SeuratDisk), error = function(e) cat("SeuratDisk not installed (not used):", conditionMessage(e), "\n"))

set.seed(404)
options(future.globals.maxSize = 3 * 1024^3)

dir.create("../output/plots/",  recursive = TRUE, showWarnings = FALSE)
dir.create("../output/tables/", recursive = TRUE, showWarnings = FALSE)
dir.create("../output/xenium/", recursive = TRUE, showWarnings = FALSE)

output_dir <- "../output/plots/"

#################################################################################
## Shared color definitions
colors <- c( "#C77CFF", "#00BB4E", "#00B0F6", "#00C1A3", "#FF61CC",
             "#A3A500", "#D89000", "#FD7083", "#CE9500", "#39B606",
             "#00B0E1", "#E28A00", "#D973FC", "#B2A000", "#00C091",
             "#F27D53", "#FF6A98", "#B186FF", "#00BAE0", "#F265E8",
             "#00B5ED", "#FF62BC", "#CD9600", "#00B92A", "#00C0B5",
             "#EA8331", "#00A9FF", "#7099FF", "#00BE67", "#35A2FF",
             "#E76BF3", "#91AA00", "#00BDD4", "#00C0B6", "#FF65AA",
             "#FD6F88", "#61B200", "#C09B06", "#39B607", "#F8766D")
colors_secondary <- c(  "#FC8D62", "#66C2A5", "#8DA0CB", "#E78AC3",
                        "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
                        "#FDB462", "#BC80BD", "#1B9E77", "#D95F02",
                        "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
                        "#A6761D", "#666666", "#A6CEE3", "#FB9A99",
                        "#B2DF8A", "#CAB2D6", "#FFFF99", "#FDBF6F",
                        "#FF7F00")
cluster_colors <- c(
  "Prox Tub" = "#91AA00", "Dist Tub" = "#66A61E", "Loop Henle" = "#00C091",
  "Endo" = "#F8766D", "Immune" = "#B186FF", "Fibro" = "#CE9500",
  "Collect Duct" = "#F27D53", "Trans Epi" = "#FF6A98", "Uni Fibro" = "#00BAE0",
  "Podo" = "#D973FC", "Ground" = "#00BB4E", "Progen" = "#C77CFF",
  "Diff Myo" = "#00C1A3", "Prolif" = "#FF61CC", "Diff Fibro" = "#A3A500",
  "TR Diff Myo" = "#7099FF", "TR Diff Chondro" = "#00A9FF", "Mural" = "#39B606",
  "Uni Fibro (Pi16)" = "#B186FF", "Uni Fibro (Col15a1)" = "#FD7083",
  "Cycle" = "#00BAE0", "Mac Densa" = "#39B607", "Adipo" = "#C09B06",
  "Prolif 1" = "#FF61CC", "Prolif 2" = "#FF65AA"
)

prefixes <- c("TMA1", "TMA2", "K1", "K2", "K3", "K4")

xenium_rds_path <- "../output/xenium/dicer1_xenium_mus_musculus_all.rds"

cell_type_markers_mus_musculus_kidney <-
  read_yaml("../reference/cell_type_markers_mus_musculus_kidney.yaml")

#################################################################################
#################################################################################
# PASS 1: WT ANALYSIS
# Load full object, immediately subset to WT+QC, delete full object
#################################################################################
cat("=== PASS 1: Loading Xenium object for WT analysis ===\n")
xenium.obj <- readRDS(xenium_rds_path)
cat("Loaded. Cells:", ncol(xenium.obj), "\n")
xenium.obj <- UpdateSeuratObject(xenium.obj)
cat("UpdateSeuratObject done.\n")
gc()

# Print global cell counts (use full object before subsetting)
cat("Sample_type counts:\n"); print(table(xenium.obj@meta.data$Sample_type))
cat("Gender counts:\n"); print(table(xenium.obj@meta.data$Gender))

# PATCH: combine QC filter + WT subset in ONE step to minimize peak RAM
# This avoids creating a temporary QC-filtered object of 1.5M cells
wt_qc_idx <- which(
  xenium.obj$nCount_Xenium > 5 &
  xenium.obj$nFeature_Xenium > 5 &
  xenium.obj$cell_area > 20 &
  xenium.obj$Inclusion_control != "no"
)
cat("WT+QC cells:", length(wt_qc_idx), "\n")
xenium.obj.wt <- xenium.obj[, wt_qc_idx]
rm(xenium.obj, wt_qc_idx)
gc()
cat("Full object removed. WT cells:", ncol(xenium.obj.wt), "\n")

# Print WT cell counts
cat("WT Lesion_ID counts:\n"); print(table(xenium.obj.wt@meta.data$Lesion_ID))
cat("WT core counts:\n"); print(table(xenium.obj.wt@meta.data$core))

# SCTransform
xenium.obj.wt <- SCTransform(xenium.obj.wt, assay = "Xenium")
# RunPCA
xenium.obj.wt <- RunPCA(xenium.obj.wt, verbose = FALSE)

tryCatch({
  png(paste0(output_dir, "dicer1_xenium_wt_dimheatmap.png"), width=1600, height=1200, res=150)
  DimHeatmap(xenium.obj.wt, dims = 1:20, cells = 500, balanced = TRUE)
  dev.off()
}, error = function(e) { try(dev.off(), silent=TRUE); cat("DimHeatmap wt skipped:", conditionMessage(e), "\n") })

tryCatch({
  p_elbow_wt <- ElbowPlot(xenium.obj.wt, ndims = 50, reduction = "pca")
  ggsave("dicer1_xenium_wt_elbowplot.png", plot=p_elbow_wt, path=output_dir, width=6, height=4)
}, error = function(e) cat("ElbowPlot wt skipped:", conditionMessage(e), "\n"))

xenium.obj.wt <- RunUMAP(xenium.obj.wt, reduction = "pca", dims = 1:30)
xenium.obj.wt <- FindNeighbors(xenium.obj.wt, reduction = "pca", dims = 1:30)
xenium.obj.wt <- FindClusters(xenium.obj.wt, resolution = 0.3)

## DimPlots
(xenium.obj.wt_umap_clusters <-
   DimPlot(xenium.obj.wt, label = T, repel = TRUE, reduction = "umap",
           cols = colors, shuffle = T, group.by = "seurat_clusters") +
   theme(plot.title = element_blank()))
ggsave("dicer1_xenium_wt_umap_clusters.svg", plot = xenium.obj.wt_umap_clusters,
       width = 5, height = 4, device = "svg", path = output_dir)

(xenium.obj.wt_umap_sample_ID <-
    DimPlot(xenium.obj.wt, label = F, repel = TRUE, reduction = "umap",
            cols = colors_secondary, shuffle = T, group.by = "Sample_ID") +
    theme(plot.title = element_blank()))
ggsave("dicer1_xenium_wt_umap_sample_id.svg", plot = xenium.obj.wt_umap_sample_ID,
       width = 5.6, height = 4, device = "svg", path = output_dir)

cat("WT cells:", ncol(xenium.obj.wt), "\n")

## Extract cluster annotations
clusters_combined_wt <- xenium.obj.wt$SCT_snn_res.0.3
all_cell_ids_wt <- names(clusters_combined_wt)
for (prefix in prefixes) {
  matching_cells_wt <- grep(paste0("^", prefix, "_"), all_cell_ids_wt, value = TRUE)
  cluster_subset_wt <- clusters_combined_wt[matching_cells_wt]
  cleaned_cell_ids_wt <- sub(paste0("^", prefix, "_"), "", names(cluster_subset_wt))
  df_export_wt <- data.frame(cell_id = cleaned_cell_ids_wt, group = cluster_subset_wt)
  file_name <- paste0("../output/tables/dicer1_xenium_wt_", prefix, "_annotation_clusters.csv")
  write.csv(df_export_wt, file = file_name, row.names = FALSE)
  cat("Exported:", file_name, "\n")
}

## FindAllMarkers
if ("SCT" %in% names(xenium.obj.wt@assays)) {
  for (i in seq_along(xenium.obj.wt@assays$SCT@SCTModel.list)) {
    slot(xenium.obj.wt@assays$SCT@SCTModel.list[[i]], name = "umi.assay") <- "Xenium"
  }
}
xenium.obj.wt <- PrepSCTFindMarkers(xenium.obj.wt)
Idents(xenium.obj.wt) <- xenium.obj.wt@meta.data$seurat_clusters
markers_wt_clusters <- FindAllMarkers(xenium.obj.wt, only.pos = TRUE)
markers_wt_clusters %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)
markers_wt_clusters %>% group_by(cluster) %>% filter(avg_log2FC > 1) %>%
  slice_head(n = 3) %>% ungroup() -> top3_wt_cluster

maxcells.wt = min(table(Idents(xenium.obj.wt)))
xenium.obj.wt_top3_genes_clusters <-
  DoHeatmap(subset(xenium.obj.wt, downsample = maxcells.wt),
            features = top3_wt_cluster$gene, raster = TRUE, angle = 45,
            draw.line = F, label = T) +
  viridis::scale_fill_viridis(option = "G")
ggsave("dicer1_xenium_wt_top3_genes_clusters.svg", plot = xenium.obj.wt_top3_genes_clusters,
       width = 15, height = 7, device = "svg", path = output_dir)
ggsave("dicer1_xenium_wt_top3_genes_clusters.png", plot = xenium.obj.wt_top3_genes_clusters,
       width = 10, height = 5, dpi = 300, path = output_dir)
write_csv(x = markers_wt_clusters, file = "../output/tables/dicer1_xenium_wt_markers_clusters.csv")

## DotPlot
(xenium.obj.wt_dotplot_clusters <-
    DotPlot(xenium.obj.wt, features = cell_type_markers_mus_musculus_kidney)
  + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  + theme(strip.text.x = element_text(angle = 90, hjust = 0)))
ggsave("dicer1_xenium_wt_dotplot_clusters.svg", plot = xenium.obj.wt_dotplot_clusters,
       width = 10, height = 8, device = "svg", path = output_dir)

## Cell type annotation
xenium.obj.wt_cell_type <- list(
  "0"="Endo", "1"="Prox Tub", "2"="Loop Henle", "3"="Prox Tub", "4"="Fibro",
  "5"="Loop Henle", "6"="Prox Tub", "7"="Collect Duct", "8"="Immune", "9"="Dist Tub",
  "10"="Prox Tub", "11"="Prox Tub", "12"="Mural", "13"="Loop Henle", "14"="Prox Tub",
  "15"="Uni Fibro", "16"="Prox Tub", "17"="Loop Henle", "18"="Collect Duct",
  "19"="Adipo", "20"="Podo", "21"="Collect Duct", "22"="Collect Duct",
  "23"="Trans Epi", "24"="Endo", "25"="Endo", "26"="Loop Henle"
)
xenium.obj.wt@meta.data$cell_type <-
  factor(x = xenium.obj.wt_cell_type[as.character(xenium.obj.wt$seurat_clusters)],
         levels = c("Uni Fibro","Fibro","Mesang","Mural","Endo","Podo","Prox Tub",
                    "Dist Tub","Loop Henle","Collect Duct","Trans Epi","Immune","Adipo"))

(xenium.obj.wt_umap_cell_type <- DimPlot(
  xenium.obj.wt, label = T, repel = TRUE, reduction = "umap",
  shuffle = T, cols = cluster_colors, group.by = "cell_type") +
   theme(plot.title = element_blank()))
ggsave("dicer1_xenium_wt_umap_cell_type.svg", plot = xenium.obj.wt_umap_cell_type,
       width = 5, height = 4, device = "svg", path = output_dir)

xenium.obj.wt_subset <- subset(xenium.obj.wt,
  subset = Sample_ID %in% c("HDT967_kidney_2", "HDT294_kidney_2"))
(xenium.obj.wt_umap_cell_type_split <- DimPlot(
  xenium.obj.wt_subset, label = T, repel = TRUE, reduction = "umap",
  shuffle = T, cols = cluster_colors, group.by = "cell_type", split.by = "Sample_age") +
    theme(plot.title = element_blank()))
ggsave("dicer1_xenium_wt_umap_cell_type_split.svg", plot = xenium.obj.wt_umap_cell_type_split,
       width = 6.6, height = 4, device = "svg", path = output_dir)
cell_counts <- xenium.obj.wt_subset@meta.data %>%
  group_by(Sample_age, cell_type) %>% summarise(n_cells = n(), .groups = "drop")
print(cell_counts)
write_csv(x = cell_counts, file = "../output/tables/dicer1_xenium_wt_subset_time_split_counts.csv")
rm(xenium.obj.wt_subset)

## Cell type annotation export
cell_types_wt_annotation <- xenium.obj.wt$cell_type
all_cell_ids_wt <- names(cell_types_wt_annotation)
for (prefix in prefixes) {
  matching_cells_wt <- grep(paste0("^", prefix, "_"), all_cell_ids_wt, value = TRUE)
  annotation_subset_wt <- cell_types_wt_annotation[matching_cells_wt]
  cleaned_cell_ids_wt <- sub(paste0("^", prefix, "_"), "", names(annotation_subset_wt))
  df_export_wt <- data.frame(cell_id = cleaned_cell_ids_wt, group = annotation_subset_wt)
  file_name <- paste0("../output/tables/dicer1_xenium_wt_", prefix, "_annotation_cell_type.csv")
  write.csv(df_export_wt, file = file_name, row.names = FALSE)
  cat("Exported:", file_name, "\n")
}

## FindAllMarkers by cell type
Idents(xenium.obj.wt) <- xenium.obj.wt@meta.data$cell_type
markers_wt_cell_type <- FindAllMarkers(xenium.obj.wt, only.pos = TRUE)
markers_wt_cell_type %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)
markers_wt_cell_type %>% group_by(cluster) %>% filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>% ungroup() -> top5_wt_cell_type

maxcells.wt <- min(table(Idents(xenium.obj.wt)))
xenium.obj.wt_top5_genes_cell_type <-
  DoHeatmap(subset(xenium.obj.wt, downsample = maxcells.wt),
            features = top5_wt_cell_type$gene, raster = T, angle = 45, draw.line = F, label = T,
            group.colors = unlist(cluster_colors[levels(xenium.obj.wt@active.ident)])) +
  viridis::scale_fill_viridis(option = "G")
ggsave("dicer1_xenium_wt_top5_genes_cell_type.svg", plot = xenium.obj.wt_top5_genes_cell_type,
       width = 15, height = 7, device = "svg", path = output_dir)
ggsave("dicer1_xenium_wt_top5_genes_cell_type.png", plot = xenium.obj.wt_top5_genes_cell_type,
       width = 15, height = 7, device = "png", path = output_dir)
write_csv(x = markers_wt_cell_type, file = "../output/tables/dicer1_xenium_wt_markers_cell_type.csv")

(xenium.obj.wt_dotplot_cell_types <-
    DotPlot(xenium.obj.wt, features = cell_type_markers_mus_musculus_kidney)
  + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  + theme(strip.text.x = element_text(angle = 90, hjust = 0)))
ggsave("dicer1_xenium_wt_dotplot_cell_types.svg", plot = xenium.obj.wt_dotplot_cell_types,
       width = 11, height = 6, device = "svg", path = output_dir)

(xenium.obj.wt_feature_plot = FeaturePlot(xenium.obj.wt, order = T,
   features = c("Bgn", "Dpt", "Mfap4", "Itga8", "Cfh", "Myh11"), keep.scale = "all"))
ggsave("dicer1_xenium_wt_feature_plot.svg", plot = xenium.obj.wt_feature_plot,
       width = 6.66, height = 8, device = "svg", path = output_dir)

(xenium.obj.wt_Vln_plot = VlnPlot(xenium.obj.wt, pt.size = 0, cols = colors_secondary,
    stack = TRUE, sort = F, flip = TRUE,
    features = c("Bgn", "Dpt", "Pi16", "Mfap4", "Itga8", "Cfh", "Myh11")))
ggsave("dicer1_xenium_wt_vln_plot.svg", plot = xenium.obj.wt_Vln_plot,
       width = 9, height = 6, device = "svg", path = output_dir)

# Save WT object
saveRDS(xenium.obj.wt, "../output/xenium/dicer1_xenium_wt.rds")
cat("Saved dicer1_xenium_wt.rds\n")

# Free memory before loading MUT
rm(xenium.obj.wt, xenium.obj.wt_umap_clusters, xenium.obj.wt_umap_sample_ID,
   xenium.obj.wt_umap_cell_type, xenium.obj.wt_umap_cell_type_split,
   xenium.obj.wt_top3_genes_clusters, xenium.obj.wt_top5_genes_cell_type,
   xenium.obj.wt_dotplot_clusters, xenium.obj.wt_dotplot_cell_types,
   xenium.obj.wt_feature_plot, xenium.obj.wt_Vln_plot,
   markers_wt_clusters, markers_wt_cell_type, top3_wt_cluster, top5_wt_cell_type)
gc()
cat("WT analysis complete. Starting MUT analysis...\n")


#################################################################################
#################################################################################
# PASS 2: MUT ANALYSIS
# Re-load full object, immediately subset to MUT+QC, delete full object
#################################################################################
cat("=== PASS 2: Loading Xenium object for MUT analysis ===\n")
xenium.obj <- readRDS(xenium_rds_path)
xenium.obj <- UpdateSeuratObject(xenium.obj)
gc()

# PATCH: combine QC filter + MUT subset in ONE step
mut_qc_idx <- which(
  xenium.obj$nCount_Xenium > 5 &
  xenium.obj$nFeature_Xenium > 5 &
  xenium.obj$cell_area > 20 &
  xenium.obj$Inclusion_mutant_primary != "no"
)
cat("MUT+QC cells:", length(mut_qc_idx), "\n")
xenium.obj.mut <- xenium.obj[, mut_qc_idx]
rm(xenium.obj, mut_qc_idx)
gc()
cat("Full object removed. MUT cells:", ncol(xenium.obj.mut), "\n")

# SCTransform
xenium.obj.mut <- SCTransform(xenium.obj.mut, assay = "Xenium")
xenium.obj.mut <- RunPCA(xenium.obj.mut, verbose = FALSE)

tryCatch({
  png(paste0(output_dir, "dicer1_xenium_mut_dimheatmap.png"), width=1600, height=1200, res=150)
  DimHeatmap(xenium.obj.mut, dims = 1:20, cells = 500, balanced = TRUE)
  dev.off()
}, error = function(e) { try(dev.off(), silent=TRUE); cat("DimHeatmap mut skipped:", conditionMessage(e), "\n") })

tryCatch({
  p_elbow_mut <- ElbowPlot(xenium.obj.mut, ndims = 50, reduction = "pca")
  ggsave("dicer1_xenium_mut_elbowplot.png", plot=p_elbow_mut, path=output_dir, width=6, height=4)
}, error = function(e) cat("ElbowPlot mut skipped:", conditionMessage(e), "\n"))

xenium.obj.mut <- RunUMAP(xenium.obj.mut, reduction = "pca", dims = 1:30)
xenium.obj.mut <- FindNeighbors(xenium.obj.mut, reduction = "pca", dims = 1:30)
# First round FindClusters at 0.5 (for cluster 13 exclusion)
xenium.obj.mut <- FindClusters(xenium.obj.mut, resolution = 0.5)

## Exclude cluster 13 (low feature count / absent positive markers)
xenium.obj.mut <- xenium.obj.mut [, xenium.obj.mut$seurat_clusters != 13]
# Re-run SCT + PCA + UMAP + clustering at 0.6
xenium.obj.mut <- SCTransform(xenium.obj.mut, assay = "Xenium")
xenium.obj.mut <- RunPCA(xenium.obj.mut, verbose = FALSE)

tryCatch({
  png(paste0(output_dir, "dicer1_xenium_mut_dimheatmap_v2.png"), width=1600, height=1200, res=150)
  DimHeatmap(xenium.obj.mut, dims = 1:20, cells = 500, balanced = TRUE)
  dev.off()
}, error = function(e) { try(dev.off(), silent=TRUE); cat("DimHeatmap mut v2 skipped:", conditionMessage(e), "\n") })

tryCatch({
  p_elbow_mut2 <- ElbowPlot(xenium.obj.mut, ndims = 50, reduction = "pca")
  ggsave("dicer1_xenium_mut_elbowplot_v2.png", plot=p_elbow_mut2, path=output_dir, width=6, height=4)
}, error = function(e) cat("ElbowPlot mut v2 skipped:", conditionMessage(e), "\n"))

xenium.obj.mut <- RunUMAP(xenium.obj.mut, reduction = "pca", dims = 1:30)
xenium.obj.mut <- FindNeighbors(xenium.obj.mut, reduction = "pca", dims = 1:30)
xenium.obj.mut <- FindClusters(xenium.obj.mut, resolution = 0.6)

## DimPlots
(xenium.obj.mut_umap_clusters <- DimPlot(
  xenium.obj.mut, label = T, repel = TRUE, reduction = "umap",
  shuffle = T, group.by = "seurat_clusters") +
   theme(plot.title = element_blank()))
ggsave("dicer1_xenium_mut_umap_clusters.svg", plot = xenium.obj.mut_umap_clusters,
       width = 5.5, height = 4, device = "svg", path = output_dir)

(xenium.obj.mut_umap_lesion_ID <- DimPlot(
  xenium.obj.mut, label = F, repel = TRUE, reduction = "umap",
  shuffle = T, group.by = "Sample_ID") +
    guides(colour = guide_legend(ncol = 2)) +
    theme(plot.title = element_blank()))
ggsave("dicer1_xenium_mut_umap_sample_id.svg", plot = xenium.obj.mut_umap_lesion_ID,
       width = 9, height = 4, device = "svg", path = output_dir)

(xenium.obj.mut_umap_sample_type <- DimPlot(
  xenium.obj.mut, label = F, repel = TRUE, reduction = "umap",
  shuffle = T, cols = colors_secondary, group.by = "Sample_type") +
    theme(plot.title = element_blank()))
ggsave("dicer1_xenium_mut_umap_sample_type.svg", plot = xenium.obj.mut_umap_sample_type,
       width = 6.5, height = 4, device = "svg", path = output_dir)

(xenium.obj.mut_umap_histotype <- DimPlot(
  xenium.obj.mut, label = F, repel = TRUE, reduction = "umap",
  shuffle = T, cols = colors_secondary, group.by = "Histotype") +
    theme(plot.title = element_blank()))
ggsave("dicer1_xenium_mut_umap_histotype.svg", plot = xenium.obj.mut_umap_histotype,
       width = 7.5, height = 4, device = "svg", path = output_dir)

cat("MUT cells:", ncol(xenium.obj.mut), "\n")

## Extract cluster annotations - use seurat_clusters (active after 0.6 re-clustering)
# PATCH: original used SCT_snn_res.0.5 but object was reclustered at 0.6
clusters_combined_mut <- xenium.obj.mut$seurat_clusters
all_cell_ids_mut <- names(clusters_combined_mut)
for (prefix in prefixes) {
  matching_cells_mut <- grep(paste0("^", prefix, "_"), all_cell_ids_mut, value = TRUE)
  cluster_subset_mut <- clusters_combined_mut[matching_cells_mut]
  cleaned_cell_ids_mut <- sub(paste0("^", prefix, "_"), "", names(cluster_subset_mut))
  df_export_mut <- data.frame(cell_id = cleaned_cell_ids_mut, group = cluster_subset_mut)
  file_name <- paste0("../output/tables/dicer1_xenium_mut_", prefix, "_annotation_clusters.csv")
  write.csv(df_export_mut, file = file_name, row.names = FALSE)
  cat("Exported:", file_name, "\n")
}

## FindAllMarkers for clusters
if ("SCT" %in% names(xenium.obj.mut@assays)) {
  for (i in seq_along(xenium.obj.mut@assays$SCT@SCTModel.list)) {
    slot(xenium.obj.mut@assays$SCT@SCTModel.list[[i]], name = "umi.assay") <- "Xenium"
  }
}
xenium.obj.mut <- PrepSCTFindMarkers(xenium.obj.mut)
Idents(xenium.obj.mut) <- xenium.obj.mut$seurat_clusters
markers_mut_clusters <- FindAllMarkers(xenium.obj.mut, only.pos = TRUE)
markers_mut_clusters %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)
markers_mut_clusters %>% group_by(cluster) %>% filter(avg_log2FC > 1) %>%
  slice_head(n = 3) %>% ungroup() -> top3_mut_clusters

maxcells.mut <- min(table(Idents(xenium.obj.mut)))
xenium.obj.mut_top3_genes_clusters <-
  DoHeatmap(subset(xenium.obj.mut, downsample = maxcells.mut),
            features = top3_mut_clusters$gene, raster = TRUE, angle = 45,
            draw.line = F, label = T, group.colors = colors) +
  viridis::scale_fill_viridis(option = "G")
ggsave("dicer1_xenium_mut_top3_genes_clusters.svg", plot = xenium.obj.mut_top3_genes_clusters,
       width = 15, height = 7, device = "svg", path = output_dir)
ggsave("dicer1_xenium_mut_top3_genes_clusters.png", plot = xenium.obj.mut_top3_genes_clusters,
       width = 10, height = 5, dpi = 300, path = output_dir)
write_csv(x = markers_mut_clusters, file = "../output/tables/dicer1_xenium_mut_markers_clusters.csv")

## DotPlot clusters
Idents(xenium.obj.mut) <- xenium.obj.mut$seurat_clusters
(xenium.obj.mut_dotplot_all <-
    DotPlot(xenium.obj.mut, features = cell_type_markers_mus_musculus_kidney)
  + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  + theme(strip.text.x = element_text(angle = 90, hjust = 0)))
ggsave("dicer1_xenium_mut_dotplot_all.svg", plot = xenium.obj.mut_dotplot_all,
       width = 11, height = 10, device = "svg", path = output_dir)
ggsave("dicer1_xenium_mut_dotplot_all.png", plot = xenium.obj.mut_dotplot_all,
       width = 10, height = 10, dpi = 300, path = output_dir)

## Cell type annotation
xenium.obj.mut_cell_type <- list(
  `0`="Ground", `1`="Progen", `2`="Endo", `3`="Immune", `4`="Fibro",
  `5`="Prox Tub", `6`="Loop Henle", `7`="Immune", `8`="Prox Tub", `9`="Diff Myo",
  `10`="Progen", `11`="Endo", `12`="Prox Tub", `13`="TR Diff Myo", `14`="Prox Tub",
  `15`="Loop Henle", `16`="Collect Duct", `17`="Trans Epi", `18`="Loop Henle",
  `19`="Prox Tub", `20`="Prox Tub", `21`="Ground", `22`="Prox Tub", `23`="Immune",
  `24`="Prolif", `25`="Loop Henle", `26`="Mural", `27`="Prox Tub",
  `28`="Collect Duct", `29`="Collect Duct", `30`="Endo", `31`="Dist Tub",
  `32`="Collect Duct", `33`="Podo", `34`="Immune"
)
xenium.obj.mut@meta.data$cell_type <-
  factor(x = xenium.obj.mut_cell_type[as.character(xenium.obj.mut$seurat_clusters)],
         levels = c("Progen","Ground","Prolif","TR Diff Myo","Diff Myo",
                    "Fibro","Mural","Endo","Podo","Prox Tub","Dist Tub",
                    "Loop Henle","Collect Duct","Trans Epi","Immune"))

xenium.obj.mut_cell_count_cell_type <- table(xenium.obj.mut@meta.data$cell_type)
print(xenium.obj.mut_cell_count_cell_type)

Idents(xenium.obj.mut) <- xenium.obj.mut$cell_type
(xenium.obj.mut_umap_cell_type <- DimPlot(
  xenium.obj.mut, label = T, repel = TRUE, reduction = "umap",
  shuffle = T, group.by = "cell_type",
  cols = unlist(cluster_colors[levels(xenium.obj.mut@active.ident)])) +
    theme(plot.title = element_blank()))
ggsave("dicer1_xenium_mut_umap_cell_type.svg", plot = xenium.obj.mut_umap_cell_type,
       width = 5.7, height = 4, device = "svg", path = output_dir)

## Cell type annotation export
cell_type_mut_annotation <- xenium.obj.mut$cell_type
all_cell_ids_mut <- names(cell_type_mut_annotation)
for (prefix in prefixes) {
  matching_cells_mut <- grep(paste0("^", prefix, "_"), all_cell_ids_mut, value = TRUE)
  if (length(matching_cells_mut) == 0) {
    cat("No matching cells for prefix:", prefix, "\n"); next
  }
  cell_type_subset <- cell_type_mut_annotation[matching_cells_mut]
  cleaned_cell_ids_mut <- sub(paste0("^", prefix, "_"), "", names(cell_type_subset))
  df_export_mut <- data.frame(cell_id = cleaned_cell_ids_mut, group = cell_type_subset)
  file_name <- paste0("../output/tables/dicer1_xenium_mut_", prefix, "_annotation_cell_type.csv")
  write.csv(df_export_mut, file = file_name, row.names = FALSE)
  cat("Exported:", file_name, "\n")
}

## FindAllMarkers by cell type
xenium.obj.mut <- PrepSCTFindMarkers(xenium.obj.mut)
Idents(xenium.obj.mut) <- xenium.obj.mut$cell_type
markers_mut_cell_type <- FindAllMarkers(xenium.obj.mut, only.pos = TRUE)
markers_mut_cell_type %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)
markers_mut_cell_type %>% group_by(cluster) %>% filter(avg_log2FC > 1) %>%
  slice_head(n = 3) %>% ungroup() -> top3_mut_cell_type

maxcells.mut <- min(table(Idents(xenium.obj.mut)))
xenium.obj.mut_top3_genes_cell_type <-
  DoHeatmap(subset(xenium.obj.mut, downsample = maxcells.mut),
            features = top3_mut_cell_type$gene, raster = TRUE, angle = 45,
            draw.line = F, label = T,
            group.colors = unlist(cluster_colors[levels(xenium.obj.mut@active.ident)])) +
  viridis::scale_fill_viridis(option = "G")
ggsave("dicer1_xenium_mut_top3_genes_cell_type.svg", plot = xenium.obj.mut_top3_genes_cell_type,
       width = 7, height = 7, device = "svg", path = output_dir)
ggsave("dicer1_xenium_mut_top3_genes_cell_type.png", plot = xenium.obj.mut_top3_genes_cell_type,
       width = 7, height = 7, device = "png", path = output_dir)
write_csv(x = markers_mut_cell_type, file = "../output/tables/dicer1_xenium_mut_markers_cell_types.csv")

(xenium.obj.mut_dotplot_cell_type <-
    DotPlot(xenium.obj.mut, features = cell_type_markers_mus_musculus_kidney)
  + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  + theme(strip.text.x = element_text(angle = 90, hjust = 0)))
ggsave("dicer1_xenium_mut_dotplot_cell_type.svg", plot = xenium.obj.mut_dotplot_cell_type,
       width = 11, height = 7, device = "svg", path = output_dir)

# Save MUT object
saveRDS(xenium.obj.mut, "../output/xenium/dicer1_xenium_mut.rds")
cat("Saved dicer1_xenium_mut.rds\n")

## MSC lineage subset
xenium.obj.mut.subset <-
  subset(xenium.obj.mut, idents = c("Progen","Ground","Prolif","TR Diff Myo",
                                    "Diff Myo","Fibro","Mural"))

(xenium.obj.mut.subset_feature_plot = FeaturePlot(xenium.obj.mut.subset, order = T,
   features = c("Mfap4","Itga8","Myoz2","Cenpf","Cfh","Myh11"), keep.scale = "all"))
ggsave("dicer1_xenium_mut_subset_feature_plot.svg", plot = xenium.obj.mut.subset_feature_plot,
       width = 6.66, height = 8, device = "svg", path = output_dir)

(xenium.obj.mut.subset_Vln_plot = VlnPlot(xenium.obj.mut.subset, pt.size = 0,
    cols = colors_secondary, stack = TRUE, sort = F, flip = TRUE,
    features = c("Dpt","Pi16","Mfap4","Itga8","Cfh","Myh11","Des","Car3","Myoz2","Cenpf","Sox9")))
ggsave("dicer1_xenium_mut_subset_vln_plot.svg", plot = xenium.obj.mut.subset_Vln_plot,
       width = 9, height = 6, device = "svg", path = output_dir)

cat("Xenium mus musculus analysis complete.\n")
