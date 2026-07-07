---
title: "Dicer1_mus_musculus"
author: "Felix Kommoss"
date: "2025-11-09"
---
  
  #################################################################################
## load libraries
library(Seurat)
library(tidyverse)
library(BPCells)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(DESeq2)
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
library(SeuratDisk)

#################################################################################
## set seed
set.seed(404)

#################################################################################
## load data 
# readRDS
xenium.obj <- readRDS("../output/xenium/dicer1_xenium_mus_musculus_all.rds")

#################################################################################
## set themes and output directory
no_legend_theme <- theme(legend.position = "none")

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
  "Prox Tub" = "#91AA00",
  "Dist Tub" = "#66A61E",
  "Loop Henle" = "#00C091",
  "Endo" = "#F8766D",
  "Immune" = "#B186FF",
  "Fibro" = "#CE9500",
  "Collect Duct" = "#F27D53",
  "Trans Epi" = "#FF6A98",
  "Uni Fibro" = "#00BAE0",
  "Podo" = "#D973FC",
  "Ground" = "#00BB4E",
  "Progen" = "#C77CFF",
  "Diff Myo" = "#00C1A3",
  "Prolif" = "#FF61CC",
  "Diff Fibro" = "#A3A500",
  "TR Diff Myo" = "#7099FF",
  "TR Diff Chondro" = "#00A9FF",
  "Mural" = "#39B606",
  "Uni Fibro (Pi16)" = "#B186FF",
  "Uni Fibro (Col15a1)" = "#FD7083",
  "Cycle" = "#00BAE0",
  "Mac Densa" = "#39B607",
  "Adipo" = "#C09B06",
  "Prolif 1" = "#FF61CC",
  "Prolif 2" = "#FF65AA"
)

output_dir <- "../output/plots/"
#################################################################################
## QC
# remove cells with less than  5 counts, and 5 features, and area smaller than 20 
xenium.obj <- subset(
  xenium.obj, 
  subset = nCount_Xenium > 5 & 
    nFeature_Xenium > 5 & 
    cell_area > 20
)                                         

# cell counts for xenium.obj
cell_count_core <- table(xenium.obj@meta.data$core)
print(cell_count_core)
cell_count_Lesion_ID <- table(xenium.obj@meta.data$Lesion_ID)
print(cell_count_Lesion_ID)
cell_count_Sample_type <- table(xenium.obj@meta.data$Sample_type)
print(cell_count_Sample_type)
cell_count_Gender <- table(xenium.obj@meta.data$Gender)
print(cell_count_Gender)


#################################################################################
#################################################################################
# CONTROL
#################################################################################
## subset xenium.obj to samples with wild type genotype
# subset
xenium.obj.wt <- xenium.obj [, xenium.obj$Inclusion_control != "no"]

# settings
options(future.globals.maxSize = 3 * 1024^3)

# SCTransform
xenium.obj.wt <- SCTransform(xenium.obj.wt, assay = "Xenium")
# RunPCA
xenium.obj.wt <- RunPCA(xenium.obj.wt, verbose = FALSE)
DimHeatmap(xenium.obj.wt, dims = 1:20, cells = 500, balanced = TRUE)
# determine dims with ElbowPlot
ElbowPlot(xenium.obj.wt, ndims = 50, reduction = "pca")
# RunUMAP
xenium.obj.wt <- RunUMAP(xenium.obj.wt, reduction = "pca", dims = 1:30)
xenium.obj.wt <- FindNeighbors(xenium.obj.wt, reduction = "pca", dims = 1:30)
# FindClusters using louvain algorithm
xenium.obj.wt <- FindClusters(xenium.obj.wt, resolution = 0.3)

#################################################################################
## plot clusters and lesion_ID
# DimPLot clusters
(xenium.obj.wt_umap_clusters <- 
   DimPlot(xenium.obj.wt, label = T, repel = TRUE,
           reduction = "umap", cols = colors, shuffle = T, 
           group.by = "seurat_clusters") +
   theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_wt_umap_clusters.svg", 
       plot = xenium.obj.wt_umap_clusters, width = 5, height = 4, 
       device = "svg", path = output_dir)
# DimPLot Sample_ID
(xenium.obj.wt_umap_sample_ID <- 
    DimPlot(xenium.obj.wt, label = F, repel = TRUE, 
            reduction = "umap", cols = colors_secondary,
            shuffle = T, group.by = "Sample_ID") +
    theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_wt_umap_sample_id.svg", 
       plot = xenium.obj.wt_umap_sample_ID, width = 5.6, height = 4, 
       device = "svg", path = output_dir)

## count cells in DimPlot 
# n = 622209
total_cells_wt <- ncol(xenium.obj.wt)
print(total_cells_wt)

#################################################################################
## extract annotation for cells visualization in Xenium explorer
# write CSV for individual assays

# Define prefixes to filter
prefixes <- c("TMA1", "TMA2", "K1", "K2", "K3", "K4")

# Extract cluster annotations from the combined Seurat object
clusters_combined_wt <- xenium.obj.wt$SCT_snn_res.0.3

# Get all cell IDs from cluster vector
all_cell_ids_wt <- names(clusters_combined_wt)

# Loop over each prefix and export its subset
for (prefix in prefixes) {
  # Identify cells that match the prefix
  matching_cells_wt <- grep(paste0("^", prefix, "_"), all_cell_ids_wt, value = TRUE)
  # Extract cluster values for those cells
  cluster_subset_wt <- clusters_combined_wt[matching_cells_wt]
  # Remove prefix from cell IDs
  cleaned_cell_ids_wt <- sub(paste0("^", prefix, "_"), "", names(cluster_subset_wt))
  # Build dataframe in required format
  df_export_wt <- data.frame(cell_id = cleaned_cell_ids_wt, group = cluster_subset_wt)
  # Create export filename
  file_name <- paste0("../output/tables/dicer1_xenium_wt_", prefix, "_annotation_clusters.csv")
  # Export
  write.csv(df_export_wt, file = file_name, row.names = FALSE)
  cat("Exported:", file_name, "\n")
}


#################################################################################
## find markers for all cluster

# Correct SCT models to Xenium
if ("SCT" %in% names(xenium.obj.wt@assays)) {
  for (i in seq_along(xenium.obj.wt@assays$SCT@SCTModel.list)) {
    slot(xenium.obj.wt@assays$SCT@SCTModel.list[[i]], name = "umi.assay") <- "Xenium"
  }
}

# FindAllMarkers
xenium.obj.wt <- PrepSCTFindMarkers(xenium.obj.wt)
Idents(xenium.obj.wt) <- xenium.obj.wt@meta.data$seurat_clusters
markers_wt_clusters <- FindAllMarkers(xenium.obj.wt, only.pos = TRUE)
markers_wt_clusters %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
# top 3 markers
markers_wt_clusters %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 3) %>%
  ungroup() -> top3_wt_cluster

## plot top 3 markers
# DoHeatmap (downsampled)
maxcells.wt = min(table(Idents(xenium.obj.wt)))
xenium.obj.wt_top3_genes_clusters <- 
  DoHeatmap(subset(xenium.obj.wt, downsample = maxcells.wt), 
            features = top3_wt_cluster$gene, raster = TRUE, angle = 45, 
            draw.line = F, label = T) + 
  viridis::scale_fill_viridis(option = "G") #, group.colors = colors
ggsave(filename = "dicer1_xenium_wt_top3_genes_clusters.svg", 
       plot = xenium.obj.wt_top3_genes_clusters, width = 15, height = 7, 
       device = "svg", path = output_dir)
ggsave(
  filename = "dicer1_xenium_wt_top3_genes_clusters.png",
  plot = xenium.obj.wt_top3_genes_clusters,
  width = 10,
  height = 5,
  dpi = 300,
  path = output_dir
)

## print markers for all clusters
# write CSV
write_csv(x = markers_wt_clusters, file = "../output/tables/dicer1_xenium_wt_markers_clusters.csv")

#################################################################################
## annotate cluster using broad cell type markers
# Path to YAML file for cell type markers
cell_type_markers_mus_musculus_kidney <- 
  read_yaml("../reference/cell_type_markers_mus_musculus_kidney.yaml")

## plot cell type markers
# DotPlot
(xenium.obj.wt_dotplot_clusters <- 
    DotPlot(xenium.obj.wt, features = cell_type_markers_mus_musculus_kidney) 
  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  + theme(strip.text.x = element_text(angle = 90, hjust = 0)))
ggsave(filename = "dicer1_xenium_wt_dotplot_clusters.svg", 
       plot = xenium.obj.wt_dotplot_clusters, width = 10, height = 8, 
       device = "svg", path = output_dir)

#################################################################################
## batch and rename clusters based on broad cell type markers
# cluster names

xenium.obj.wt_cell_type <- list(
  "0"  = "Endo",         # Endothelial cells
  "1"  = "Prox Tub",     # Proximal tubule
  "2"  = "Loop Henle",   # Loop of Henle
  "3"  = "Prox Tub",     # Proximal tubule
  "4"  = "Fibro",        # Fibroblasts
  "5"  = "Loop Henle",   # Loop of Henle
  "6"  = "Prox Tub",     # Proximal tubule
  "7"  = "Collect Duct", # Collecting duct
  "8"  = "Immune",       # Immune cells
  "9"  = "Dist Tub",     # Distal tubule
  "10" = "Prox Tub",     # Proximal tubule
  "11" = "Prox Tub",     # Proximal tubule
  "12" = "Mural",        # Mural cells
  "13" = "Loop Henle",   # Loop of Henle
  "14" = "Prox Tub",     # Proximal tubule
  "15" = "Uni Fibro",    # Universal fibroblasts
  "16" = "Prox Tub",     # Proximal tubule
  "17" = "Loop Henle",   # Loop of Henle
  "18" = "Collect Duct", # Collecting duct
  "19" = "Adipo",        # Adipocytes
  "20" = "Podo",         # Podocytes
  "21" = "Collect Duct", # Collecting duct
  "22" = "Collect Duct", # Collecting duct
  "23" = "Trans Epi",    # Transitional epithelium
  "24" = "Endo",         # Endothelial cells
  "25" = "Endo",         # Endothelial cells
  "26" = "Loop Henle"    # Loop of Henle
)

xenium.obj.wt@meta.data$cell_type <- 
  factor(x = 
           xenium.obj.wt_cell_type[as.character(xenium.obj.wt$seurat_clusters)],
         levels = c("Uni Fibro", "Fibro", "Mesang", "Mural", "Endo", "Podo", "Prox Tub", 
                    "Dist Tub", "Loop Henle", "Collect Duct", "Trans Epi", "Immune", "Adipo"))

#################################################################################
## plot cell type
# DimPlot cell types
(xenium.obj.wt_umap_cell_type <- DimPlot(
  xenium.obj.wt, label = T, repel = TRUE,reduction = "umap",
  shuffle = T, cols = cluster_colors, group.by = "cell_type") +
   theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_wt_umap_cell_type.svg", 
       plot = xenium.obj.wt_umap_cell_type, width = 5, height = 4, 
       device = "svg", path = output_dir)

# DimPlot cell types split by age for kidney 1,2 and 3
# Subset the Seurat object to include only the three desired sample_IDs
xenium.obj.wt_subset <- subset(xenium.obj.wt, subset = Sample_ID %in% c("HDT967_kidney_2", "HDT294_kidney_2"))
# DimPlot cell types
(xenium.obj.wt_umap_cell_type_split <- DimPlot(
  xenium.obj.wt_subset, label = T, repel = TRUE,reduction = "umap",
  shuffle = T, cols = cluster_colors, group.by = "cell_type",
  split.by = "Sample_age" ) +
    theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_wt_umap_cell_type_split.svg", 
       plot = xenium.obj.wt_umap_cell_type_split, width = 6.6, height = 4, 
       device = "svg", path = output_dir)
# Count cells by sample_age and cell_type in the subset object
cell_counts <- xenium.obj.wt_subset@meta.data %>%
  group_by(Sample_age, cell_type) %>%
  summarise(n_cells = n(), .groups = "drop")

# View the result
print(cell_counts)
# Save as Excel file in the output directory
write_csv(x = cell_counts, file = "../output/tables/dicer1_xenium_wt_subset_time_split_counts.csv")
#################################################################################
## extract annotation for cells visualization in Xenium explorer
# write CSV

# Define prefixes to filter
prefixes <- c("TMA1", "TMA2", "K1", "K2", "K3", "K4")

# Extract cell type annotations
cell_types_wt_annotation <- xenium.obj.wt$cell_type

# Get all cell IDs
all_cell_ids_wt <- names(cell_types_wt_annotation)

# Loop over each prefix and export its subset
for (prefix in prefixes) {
  # Identify cells that match the prefix
  matching_cells_wt <- grep(paste0("^", prefix, "_"), all_cell_ids_wt, value = TRUE)
  # Extract cell type values for those cells
  annotation_subset_wt <- cell_types_wt_annotation[matching_cells_wt]
  # Remove prefix from cell IDs (to match Xenium format)
  cleaned_cell_ids_wt <- sub(paste0("^", prefix, "_"), "", names(annotation_subset_wt))
  # Build dataframe
  df_export_wt <- data.frame(cell_id = cleaned_cell_ids_wt, group = annotation_subset_wt)
  # Create export filename
  file_name <- paste0("../output/tables/dicer1_xenium_wt_", prefix, "_annotation_cell_type.csv")
  # Export to CSV
  write.csv(df_export_wt, file = file_name, row.names = FALSE)
  cat("Exported:", file_name, "\n")
}

#################################################################################
## find markers for batched clusters
# FindAllMarkers
Idents(xenium.obj.wt) <- xenium.obj.wt@meta.data$cell_type
markers_wt_cell_type <- FindAllMarkers(xenium.obj.wt, only.pos = TRUE)
markers_wt_cell_type %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
# top 5 markers
markers_wt_cell_type %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5_wt_cell_type

## plot top 5 markers
# DoHeatmap (downsampled)
maxcells.wt <- min(table(Idents(xenium.obj.wt)))
xenium.obj.wt_top5_genes_cell_type <- 
  DoHeatmap(subset(xenium.obj.wt, downsample = maxcells.wt), 
            features = top5_wt_cell_type$gene, raster = T, angle = 45, 
            draw.line = F, label = T, group.colors = 
              unlist(cluster_colors[levels(xenium.obj.wt@active.ident)])) + 
  viridis::scale_fill_viridis(option = "G") 
ggsave(filename = "dicer1_xenium_wt_top5_genes_cell_type.svg", 
       plot = xenium.obj.wt_top5_genes_cell_type, width = 15, height = 7, 
       device = "svg", path = output_dir)
ggsave(filename = "dicer1_xenium_wt_top5_genes_cell_type.png", 
       plot = xenium.obj.wt_top5_genes_cell_type, width = 15, height = 7, 
       device = "png", path = output_dir)

## print markers for 10 batched clusters
# write CSV
write_csv(x = markers_wt_cell_type, file = "../output/tables/dicer1_xenium_wt_markers_cell_type.csv")

#################################################################################
## annotate cluster using broad cell type markers
# Path to YAML file for cell type markers
cell_type_markers_mus_musculus_kidney <- 
  read_yaml("../reference/cell_type_markers_mus_musculus_kidney.yaml")

## plot cell type markers
# DotPlot
(xenium.obj.wt_dotplot_cell_types <- 
    DotPlot(xenium.obj.wt, features = cell_type_markers_mus_musculus_kidney) 
  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  + theme(strip.text.x = element_text(angle = 90, hjust = 0)))
ggsave(filename = "dicer1_xenium_wt_dotplot_cell_types.svg", 
       plot = xenium.obj.wt_dotplot_cell_types, width = 11, height = 6, 
       device = "svg", path = output_dir)

#################################################################################
## plot key genes
# FeaturePlot
(xenium.obj.wt_feature_plot = 
   FeaturePlot(xenium.obj.wt, 
               order = T, 
               features = c("Bgn", "Dpt", "Mfap4", "Itga8", "Cfh", "Myh11"), 
               keep.scale = "all"))
ggsave(filename = "dicer1_xenium_wt_feature_plot.svg", 
       plot = xenium.obj.wt_feature_plot, width = 6.66, 
       height = 8, device = "svg", path = output_dir)

# VlnPlot
(xenium.obj.wt_Vln_plot = 
    VlnPlot(xenium.obj.wt, 
            pt.size = 0, cols = colors_secondary, stack = TRUE, sort = F, 
            flip = TRUE, features = 
              c("Bgn", "Dpt", "Pi16", "Mfap4", "Itga8", "Cfh", "Myh11")))
ggsave(filename = "dicer1_xenium_wt_vln_plot.svg", 
       plot = xenium.obj.wt_Vln_plot, width = 9, 
       height = 6, device = "svg", path = output_dir)

#################################################################################
## save xenium.obj.wt on disk
# SaveRDS
saveRDS(xenium.obj.wt, "../output/xenium/dicer1_xenium_wt.rds")
xenium.obj.wt <- readRDS("../output/xenium/dicer1_xenium_wt.rds")



#################################################################################
#################################################################################
#HDT
#################################################################################
saveRDS(xenium.obj.mut, "./output/xenium/dicer1_xenium_mut.rds")
xenium.obj.mut <- readRDS("./output/xenium/dicer1_xenium_mut.rds")

## subset xenium.obj to primary samples with mutant genotype
# subset
xenium.obj.mut <- xenium.obj [, xenium.obj$Inclusion_mutant_primary != "no"]
# SCTransform
xenium.obj.mut <- SCTransform(xenium.obj.mut, assay = "Xenium")
# RunPCA
xenium.obj.mut <- RunPCA(xenium.obj.mut, verbose = FALSE)
DimHeatmap(xenium.obj.mut, dims = 1:20, cells = 500, balanced = TRUE)
# determine dims with ElbowPlot
ElbowPlot(xenium.obj.mut, ndims = 50, reduction = "pca")
# RunUMAP
xenium.obj.mut <- RunUMAP(xenium.obj.mut, reduction = "pca", dims = 1:30)
xenium.obj.mut <- FindNeighbors(xenium.obj.mut, reduction = "pca", dims = 1:30)
# FindClusters using louvain algorithm
xenium.obj.mut <- FindClusters(xenium.obj.mut, resolution = 0.5) # 0.6

(FeaturePlot(xenium.obj.mut, features = c("nFeature_Xenium")))

#################################################################################
## exclude cluster 13
## low feature count detected and absence of positive markers above background expression
# subset
xenium.obj.mut <- xenium.obj.mut [, xenium.obj.mut$seurat_clusters != 13]
# SCTransform
xenium.obj.mut <- SCTransform(xenium.obj.mut, assay = "Xenium")
# RunPCA
xenium.obj.mut <- RunPCA(xenium.obj.mut, verbose = FALSE)
DimHeatmap(xenium.obj.mut, dims = 1:20, cells = 500, balanced = TRUE)
# determine dims with ElbowPlot
ElbowPlot(xenium.obj.mut, ndims = 50, reduction = "pca")
# RunUMAP
xenium.obj.mut <- RunUMAP(xenium.obj.mut, reduction = "pca", dims = 1:30)
xenium.obj.mut <- FindNeighbors(xenium.obj.mut, reduction = "pca", dims = 1:30)
# FindClusters using louvain algorithm
xenium.obj.mut <- FindClusters(xenium.obj.mut, resolution = 0.6)

#################################################################################
## plot clusters, lesion_ID and sample type
# DimPlot clusters
(xenium.obj.mut_umap_clusters <- DimPlot(
  xenium.obj.mut, label = T, repel = TRUE, reduction = "umap",
  shuffle = T, group.by = "seurat_clusters") +
   theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_mut_umap_clusters.svg", 
       plot = xenium.obj.mut_umap_clusters, width = 5.5, height = 4, 
       device = "svg", path = output_dir)
# DimPlot sample_ID
(xenium.obj.mut_umap_lesion_ID <- DimPlot(
  xenium.obj.mut, label = F, repel = TRUE, reduction = "umap",
  shuffle = T, group.by = "Sample_ID") + 
    guides(colour = guide_legend(ncol = 2)) +
    theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_mut_umap_sample_id.svg", 
       plot = xenium.obj.mut_umap_lesion_ID, width = 9, height = 4, 
       device = "svg", path = output_dir)
# DimPlot sample type
(xenium.obj.mut_umap_sample_type <- DimPlot(
  xenium.obj.mut, label = F, repel = TRUE, reduction = "umap",
  shuffle = T, cols = colors_secondary, group.by = "Sample_type") +
    theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_mut_umap_sample_type.svg", 
       plot = xenium.obj.mut_umap_sample_type, width = 6.5, height = 4, 
       device = "svg", path = output_dir)
# DimPlot histotype
(xenium.obj.mut_umap_sample_type <- DimPlot(
  xenium.obj.mut, label = F, repel = TRUE, reduction = "umap",
  shuffle = T, cols = colors_secondary, group.by = "Histotype") +
    theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_mut_umap_histotype.svg", 
       plot = xenium.obj.mut_umap_sample_type, width = 7.5, height = 4, 
       device = "svg", path = output_dir)

## count cells in DimPlot
# n = 692374
total_cells_mut <- ncol(xenium.obj.mut)
print(total_cells_mut)

#################################################################################
## extract annotation for cells visualization in Xenium explorer
# write CSV

# Define prefixes to filter
prefixes <- c("TMA1", "TMA2", "K1", "K2", "K3", "K4")

# Extract cluster annotations from the combined Seurat object
clusters_combined_mut <- xenium.obj.mut$SCT_snn_res.0.5

# Get all cell IDs from cluster vector
all_cell_ids_mut <- names(clusters_combined_mut)

# Loop over each prefix and export its subset
for (prefix in prefixes) {
  # Identify cells that match the prefix
  matching_cells_mut <- grep(paste0("^", prefix, "_"), all_cell_ids_mut, value = TRUE)
  # Extract cluster values for those cells
  cluster_subset_mut <- clusters_combined_mut[matching_cells_mut]
  # Remove prefix from cell IDs
  cleaned_cell_ids_mut <- sub(paste0("^", prefix, "_"), "", names(cluster_subset_mut))
  # Build dataframe in required format‚
  df_export_mut <- data.frame(cell_id = cleaned_cell_ids_mut, group = cluster_subset_mut)
  # Create export filename
  file_name <- paste0("../output/tables/dicer1_xenium_mut_", prefix, "_annotation_clusters.csv")
  # Export
  write.csv(df_export_mut, file = file_name, row.names = FALSE)
  cat("Exported:", file_name, "\n")
}

#################################################################################
## find markers for all cluster

# Correct SCT models to Xenium
if ("SCT" %in% names(xenium.obj.mut@assays)) {
  for (i in seq_along(xenium.obj.mut@assays$SCT@SCTModel.list)) {
    slot(xenium.obj.mut@assays$SCT@SCTModel.list[[i]], name = "umi.assay") <- "Xenium"
  }
}


# FindAllMarkers
xenium.obj.mut <- PrepSCTFindMarkers(xenium.obj.mut)
Idents(xenium.obj.mut) <- xenium.obj.mut$seurat_clusters
markers_mut_clusters <- FindAllMarkers(xenium.obj.mut, only.pos = TRUE)
markers_mut_clusters %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
# top 3 markers
markers_mut_clusters %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 3) %>%
  ungroup() -> top3_mut_clusters

## plot top 3 markers per cluster 
# DoHeatmap (downsampled)
maxcells.mut <- min(table(Idents(xenium.obj.mut)))
xenium.obj.mut_top3_genes_clusters <- 
  DoHeatmap(subset(xenium.obj.mut, downsample = maxcells.mut), 
            features = top3_mut_clusters$gene, raster = TRUE, angle = 45, 
            draw.line = F, label = T, group.colors = colors) + 
  viridis::scale_fill_viridis(option = "G" )
ggsave(filename = "dicer1_xenium_mut_top3_genes_clusters.svg", 
       plot = xenium.obj.mut_top3_genes_clusters, width = 15, height = 7, 
       device = "svg", path = output_dir)
ggsave(
  filename = "dicer1_xenium_mut_top3_genes_clusters.png",
  plot = xenium.obj.mut_top3_genes_clusters,
  width = 10,
  height = 5,
  dpi = 300,
  path = output_dir
)
## print markers for all 33 clusters
# write CSV
write_csv(x = markers_mut_clusters, file = "../output/tables/dicer1_xenium_mut_markers_clusters.csv")

#################################################################################
## annotate cluster using broad cell type markers
# Path to YAML file for cell type markers
cell_type_markers_mus_musculus_kidney <- 
  read_yaml("../reference/cell_type_markers_mus_musculus_kidney.yaml")

## plot cell type markers‚
# DotPlot
Idents(xenium.obj.mut) <- xenium.obj.mut$seurat_clusters
(xenium.obj.mut_dotplot_all <- 
    DotPlot(xenium.obj.mut, features = cell_type_markers_mus_musculus_kidney) 
  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  + theme(strip.text.x = element_text(angle = 90, hjust = 0)))
ggsave(filename = "dicer1_xenium_mut_dotplot_all.svg", 
       plot = xenium.obj.mut_dotplot_all, width = 11, height = 10, 
       device = "svg", path = output_dir)
ggsave(
  filename = "dicer1_xenium_mut_dotplot_all.png",
  plot = xenium.obj.mut_dotplot_all,
  width = 10,
  height = 10,
  dpi = 300,
  path = output_dir
)
#################################################################################
## batch and rename clusters based on broad cell type markers
# cluster names

xenium.obj.mut_cell_type <- list(
  `0`  = "Ground",        # Ground
  `1`  = "Progen",        # Progenitor
  `2`  = "Endo",          # Endothelium
  `3`  = "Immune",        # Immune
  `4`  = "Fibro",         # Fibroblast
  `5`  = "Prox Tub",      # Proximal tubule
  `6`  = "Loop Henle",    # Loop of Henle
  `7`  = "Immune",        # Immune
  `8`  = "Prox Tub",      # Proximal tubule
  `9`  = "Diff Myo",      # Differentiated myogenic
  `10` = "Progen",        # Progenitor
  `11` = "Endo",          # Endothelium
  `12` = "Prox Tub",      # Proximal tubule
  `13` = "TR Diff Myo",   # Transiting differentiated myogenic
  `14` = "Prox Tub",      # Proximal tubule
  `15` = "Loop Henle",    # Loop of Henle
  `16` = "Collect Duct",  # Collecting duct
  `17` = "Trans Epi",     # Transitional epithelium
  `18` = "Loop Henle",    # Loop of Henle
  `19` = "Prox Tub",      # Proximal tubule
  `20` = "Prox Tub",      # Proximal tubule
  `21` = "Ground",        # Ground
  `22` = "Prox Tub",      # Proximal tubule
  `23` = "Immune",        # Immune
  `24` = "Prolif",        # Proliferative
  `25` = "Loop Henle",    # Loop of Henle
  `26` = "Mural",         # Mural cells
  `27` = "Prox Tub",      # Proximal tubule
  `28` = "Collect Duct",  # Collecting duct
  `29` = "Collect Duct",  # Collecting duct
  `30` = "Endo",          # Endothelium
  `31` = "Dist Tub",      # Distal tubule
  `32` = "Collect Duct",  # Collecting duct
  `33` = "Podo",          # Podocyte
  `34` = "Immune"         # Immune
)

xenium.obj.mut@meta.data$cell_type <- 
  factor(x = 
           xenium.obj.mut_cell_type[as.character(xenium.obj.mut$seurat_clusters)],
         levels = c("Progen", "Ground", "Prolif", "TR Diff Myo", "Diff Myo", 
                    "Fibro", "Mural", "Endo", "Podo", "Prox Tub", "Dist Tub",
                    "Loop Henle", "Collect Duct", "Trans Epi", "Immune"))

# count cells in cell types
xenium.obj.mut_cell_count_cell_type <- table(xenium.obj.mut@meta.data$cell_type)
print(xenium.obj.mut_cell_count_cell_type)

#################################################################################
## plot cell types
# DimPLot
Idents(xenium.obj.mut) <- xenium.obj.mut$cell_type
(xenium.obj.mut_umap_cell_type <- DimPlot(
  xenium.obj.mut, label = T, repel = TRUE, reduction = "umap",
  shuffle = T, group.by = "cell_type", cols = 
    unlist(cluster_colors[levels(xenium.obj.mut@active.ident)])) +
    theme(plot.title = element_blank()))
ggsave(filename = "dicer1_xenium_mut_umap_cell_type.svg", 
       plot = xenium.obj.mut_umap_cell_type, width = 5.7, height = 4, 
       device = "svg", path = output_dir)

#################################################################################
## extract annotation for cells visualization in Xenium explorer
# write CSV

# Define prefixes to filter
prefixes <- c("TMA1", "TMA2", "K1", "K2", "K3", "K4")

# Extract cell type annotations
cell_type_mut_annotation <- xenium.obj.mut$cell_type

# Get all cell IDs from annotation
all_cell_ids_mut <- names(cell_type_mut_annotation)

# Loop over each prefix and export its subset
for (prefix in prefixes) {
  # Identify cells that match the prefix (e.g., "K1_", "TMA2_", etc.)
  matching_cells_mut <- grep(paste0("^", prefix, "_"), all_cell_ids_mut, value = TRUE)
  
  # Skip if no matching cells found
  if (length(matching_cells_mut) == 0) {
    cat("No matching cells found for prefix:", prefix, "\n")
    next
  }
  # Extract cell type values for matching cells
  cell_type_subset <- cell_type_mut_annotation[matching_cells_mut]
  # Remove prefix from cell IDs (e.g., "K1_123" → "123")
  cleaned_cell_ids_mut <- sub(paste0("^", prefix, "_"), "", names(cell_type_subset))
  # Build dataframe
  df_export_mut <- data.frame(cell_id = cleaned_cell_ids_mut, group = cell_type_subset)
  # Create export filename
  file_name <- paste0("../output/tables/dicer1_xenium_mut_", prefix, "_annotation_cell_type.csv")
  # Export CSV
  write.csv(df_export_mut, file = file_name, row.names = FALSE)
  cat("Exported:", file_name, "\n")
}


#################################################################################
## find markers for cell types

# Correct SCT models to Xenium
if ("SCT" %in% names(xenium.obj.mut@assays)) {
  for (i in seq_along(xenium.obj.mut@assays$SCT@SCTModel.list)) {
    slot(xenium.obj.mut@assays$SCT@SCTModel.list[[i]], name = "umi.assay") <- "Xenium"
  }
}

# FindAllMarkers
xenium.obj.mut <- PrepSCTFindMarkers(xenium.obj.mut)
Idents(xenium.obj.mut) <- xenium.obj.mut$cell_type
markers_mut_cell_type <- FindAllMarkers(xenium.obj.mut, only.pos = TRUE)
markers_mut_cell_type %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
# top 3 markers
markers_mut_cell_type %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 3) %>%
  ungroup() -> top3_mut_cell_type

## plot top 3 markers
# DoHeatmap (downsampled)
maxcells.mut <- min(table(Idents(xenium.obj.mut)))
xenium.obj.mut_top3_genes_cell_type <- 
  DoHeatmap(subset(xenium.obj.mut, downsample = maxcells.mut), 
            features = top3_mut_cell_type$gene, raster = TRUE, angle = 45, 
            draw.line = F, label = T, group.colors = 
              unlist(cluster_colors[levels(xenium.obj.mut@active.ident)])) + 
  viridis::scale_fill_viridis(option = "G") 
ggsave(filename = "dicer1_xenium_mut_top3_genes_cell_type.svg", 
       plot = xenium.obj.mut_top3_genes_cell_type, width = 7, height = 7, 
       device = "svg", path = output_dir)
ggsave(filename = "dicer1_xenium_mut_top3_genes_cell_type.png", 
       plot = xenium.obj.mut_top3_genes_cell_type, width = 7, height = 7, 
       device = "png", path = output_dir)

## print markers for cell types
# write CSV
write_csv(x = markers_mut_cell_type, file = "../output/tables/dicer1_xenium_mut_markers_cell_types.csv")

#################################################################################
## annotate cluster using broad cell type markers
# Path to YAML file for cell type markers
cell_type_markers_mus_musculus_kidney <- 
  read_yaml("../reference/cell_type_markers_mus_musculus_kidney.yaml")

## plot cell type markers
# DotPlot
(xenium.obj.mut_dotplot_cell_type <- 
    DotPlot(xenium.obj.mut, features = cell_type_markers_mus_musculus_kidney) 
  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  + theme(strip.text.x = element_text(angle = 90, hjust = 0)))
ggsave(filename = "dicer1_xenium_mut_dotplot_cell_type.svg", 
       plot = xenium.obj.mut_dotplot_cell_type, width = 11, height = 7, 
       device = "svg", path = output_dir)

#################################################################################
## save xenium.obj on disk 
saveRDS(xenium.obj.mut, "./output/xenium/dicer1_xenium_mut.rds")
xenium.obj.mut <- readRDS("./output/xenium/dicer1_xenium_mut.rds")

#################################################################################
## subset to MSC lineage
# subset
xenium.obj.mut.subset <- 
  subset(xenium.obj.mut, idents = c("Progen", "Ground", "Prolif", "TR Diff Myo",
                                    "Diff Myo", "Fibro", "Mural"))

#################################################################################
## plot key genes
# FeaturePlot
(xenium.obj.mut.subset_feature_plot = 
   FeaturePlot(xenium.obj.mut.subset, 
               order = T, 
               features = c("Mfap4", "Itga8", "Myoz2", "Cenpf", "Cfh", "Myh11"), 
               keep.scale = "all"))
ggsave(filename = "dicer1_xenium_mut_subset_feature_plot.svg", 
       plot = xenium.obj.mut.subset_feature_plot, width = 6.66, 
       height = 8, device = "svg", path = output_dir)

# VlnPlot
(xenium.obj.mut.subset_Vln_plot = 
    VlnPlot(xenium.obj.mut.subset, 
            pt.size = 0, cols = colors_secondary, stack = TRUE, sort = F, 
            flip = TRUE, features = 
              c("Dpt", "Pi16", "Mfap4", "Itga8", "Cfh", "Myh11", "Des", "Car3", "Myoz2", 
                "Cenpf", "Sox9")))
ggsave(filename = "dicer1_xenium_mut_subset_vln_plot.svg", 
       plot = xenium.obj.mut.subset_Vln_plot, width = 9, 
       height = 6, device = "svg", path = output_dir)