# dicer1_scrna_tumors_filter_only_patched.R
# Runs ONLY the filtered (mesenchymal) section of dicer1_scrna_tumors.R.
# Loads the already-saved unfilter_tumor.rds instead of re-running the full script.
# Memory management (rm/gc) added before each FindAllMarkers to prevent OOM (SIGKILL exit 137).

library(SingleCellExperiment)
source("/home/node/MAURICIO/papers/653_Spatial_single_cell_transcriptomic_analysis_informs_tumor_developmental_hierarchy_of_DICER1_syndrome_related_sarcoma_2026/REPRODUCIBILITY_ANALYSIS/20260703T220654Z__anthropic-claude-sonnet-4-6__adaptive/workflow_draft/scripts/patched/dotplot_compat_patch.R")

library(dplyr)
library(Seurat)
library(scater)
library(cowplot)
library(ggplot2)
library(yaml)

set.seed(1337)

cat("Loading unfilter_tumor.rds ...\n")
sce_integrated_geno <- readRDS(file = "../output/scrna/dicer1_scrna_unfilter_tumor.rds")
cat("Loaded. Cells:", ncol(sce_integrated_geno), "\n")

dot_plot_genes <- read_yaml(file = "../reference/cell_type_markers_mus_musculus_kidney.yaml")
dot_plot_genes_mes <- read_yaml(file = "../reference/cell_type_markers_mus_musculus_mesenchymal.yaml")

global_theme <- theme(plot.title = element_blank(),
                      axis.line = element_line(size = 0.3))

colors <- c( "#C77CFF", "#00BB4E", "#00B0F6", "#00C1A3", "#FF61CC",
             "#A3A500", "#D89000", "#FD7083", "#CE9500", "#39B606",
             "#00B0E1", "#E28A00", "#D973FC", "#B2A000", "#00C091",
             "#F27D53", "#FF6A98", "#B186FF", "#00BAE0", "#F265E8",
             "#00B5ED", "#FF62BC", "#CD9600", "#00B92A", "#00C0B5",
             "#EA8331", "#00A9FF", "#7099FF", "#00BE67", "#35A2FF",
             "#E76BF3", "#91AA00", "#00BDD4", "#00C0B6", "#FF65AA",
             "#FD6F88", "#61B200", "#C09B06", "#39B607", "#F8766D")

colors_secondary <- c( "#FC8D62", "#66C2A5", "#8DA0CB", "#E78AC3",
                       "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
                       "#FDB462", "#BC80BD")

cluster_colors <- list(
  "Prox Tub" = "#B2A000",
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
  "Mural" = "#39B606",
  "Uni Fibro (Pi16)" = "#FD7083",
  "Fibro Perturb" = "#EA8331",
  "Uni Fibro (Col15a1)" = "#B186FF",
  "Cycle" = "#00BAE0",
  "Mac Densa" = "#39B607",
  "Peri" = "#C09B06",
  "SM Vasc" = "#91AA00",
  "Unk" = "black"
)

### FILTERED SECTION
cat("Subsetting to mesenchymal cells ...\n")
mes_cells <- sce_integrated_geno$cell_type %in% c("Uni Fibro",
                                                   "Fibro",
                                                   "Cycle",
                                                   "TR Diff Myo",
                                                   "Diff Myo",
                                                   "Mural")
keep_cells <- sce_integrated_geno$Tomato > 0

# Use layer= for Seurat v5 compatibility; fall back to slot= if needed
counts <- tryCatch(
  GetAssayData(object = sce_integrated_geno, assay = "RNA", layer = "counts"),
  error = function(e) GetAssayData(object = sce_integrated_geno, assay = "RNA", slot = "counts")
)
ptprc_cells  <- counts["Ptprc", ] > 0
epcam_cells  <- counts["Epcam", ] > 0
trans_epi_cells <- counts["Krt19", ] > 0
discard_cells <- ptprc_cells | epcam_cells | trans_epi_cells

sce_integrated_sample_filter <- sce_integrated_geno[, mes_cells & keep_cells & !discard_cells]

sce_integrated_sample_filter <- sce_integrated_sample_filter[,
  sce_integrated_sample_filter$Sample %in%
  names(table(sce_integrated_sample_filter$Sample))[table(sce_integrated_sample_filter$Sample) > 50]]

cat("Filtered cells:", ncol(sce_integrated_sample_filter), "\n")

# Free the large unfiltered object before re-integration
rm(sce_integrated_geno, counts, ptprc_cells, epcam_cells, trans_epi_cells, discard_cells, mes_cells, keep_cells)
gc()
cat("Freed unfiltered object. Starting re-integration ...\n")

seurat_list_sample_filter <- SplitObject(sce_integrated_sample_filter, split.by = "Sample")
seurat_list_sample_filter <- lapply(seurat_list_sample_filter, function(seurat_obj) NormalizeData(seurat_obj, verbose = FALSE))
seurat_list_sample_filter <- lapply(seurat_list_sample_filter, function(seurat_obj) FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000))

features_sample_filter <- SelectIntegrationFeatures(object.list = seurat_list_sample_filter, nfeatures = 2000)

seurat_list_sample_filter <- lapply(seurat_list_sample_filter, function(seurat_obj) ScaleData(seurat_obj, features = features_sample_filter))
seurat_list_sample_filter <- lapply(seurat_list_sample_filter, function(seurat_obj) RunPCA(seurat_obj, features = features_sample_filter))

anchors_sample_filter <- FindIntegrationAnchors(object.list = seurat_list_sample_filter, anchor.features = features_sample_filter, reduction = "rpca")

sce_integrated_sample_filter <- IntegrateData(anchorset = anchors_sample_filter, k.weight = 25)

# Free memory after integration
rm(seurat_list_sample_filter, anchors_sample_filter, features_sample_filter)
gc()

sce_integrated_sample_filter <- ScaleData(sce_integrated_sample_filter, do.scale = TRUE, do.center = TRUE, assay = "integrated")
sce_integrated_sample_filter <- RunPCA(sce_integrated_sample_filter, npcs = 50, assay = "integrated", seed.use = 1337)

plot_elbow_filter <- ElbowPlot(sce_integrated_sample_filter, ndims = 50)

sce_integrated_sample_filter <- RunTSNE(sce_integrated_sample_filter, dims = 1:20, assay = "integrated", dim.embed = 3, seed.use = 300, perplexity = 20)
sce_integrated_sample_filter <- RunUMAP(sce_integrated_sample_filter, reduction = "pca", dims = 1:20, assay = "integrated", seed.use = 1337)
sce_integrated_sample_filter <- FindNeighbors(sce_integrated_sample_filter, reduction = "pca", dims = 1:20, assay = "integrated")
sce_integrated_sample_filter <- FindClusters(sce_integrated_sample_filter, resolution = 0.2, random.seed = 1337)

DefaultAssay(object = sce_integrated_sample_filter) <- "RNA"
sce_integrated_sample_filter <- JoinLayers(sce_integrated_sample_filter)
sce_integrated_sample_filter <- ScaleData(object = sce_integrated_sample_filter)

plot_clusters_filter <- DimPlot(object = sce_integrated_sample_filter, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE) +
  global_theme +
  scale_colour_manual(values = colors)

plot_samples_filter <- DimPlot(object = sce_integrated_sample_filter, reduction = "umap", group.by = "Name", shuffle = TRUE) +
  global_theme +
  scale_colour_manual(values = colors_secondary)

plot_time_filter <- DimPlot(object = sce_integrated_sample_filter, reduction = "umap", group.by = "Time", shuffle = TRUE) +
  global_theme +
  scale_colour_manual(values = colors_secondary)

plot_cell_dot_kidney_filter <- safe_dotplot(object = sce_integrated_sample_filter, features = dot_plot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

plot_cell_dot_mes_filter <- safe_dotplot(object = sce_integrated_sample_filter, features = dot_plot_genes_mes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

# FindAllMarkers 1: cluster markers (only.pos=TRUE)
cat("Running FindAllMarkers (cluster, only.pos=TRUE) ...\n")
gc()
markers_cluster_filter <- FindAllMarkers(object = sce_integrated_sample_filter, only.pos = TRUE)

top_markers_cluster_filter <- markers_cluster_filter %>%
  group_by(cluster) %>%
  slice_head(n = 30)

downsample_cluster_filter <- min(table(sce_integrated_sample_filter$seurat_clusters))
plot_heatmap_cluster_filter <- safe_heatmap(object = subset(sce_integrated_sample_filter, downsample = downsample_cluster_filter),
                                         features = top_markers_cluster_filter$gene,
                                         raster = TRUE,
                                         angle = 45,
                                         draw.lines = FALSE,
                                         label = TRUE) +
  viridis::scale_fill_viridis(option = "G")

cluster_labels_filter <- c("0" = "Progen",
                           "1" = "Ground",
                           "2" = "TR Diff Myo",
                           "3" = "Ground",
                           "4" = "Prolif",
                           "5" = "Diff Myo",
                           "6" = "Progen",
                           "7" = "Peri",
                           "8" = "SM Vasc",
                           "9" = "Fibro Perturb")

sce_integrated_sample_filter$cell_type <- factor(
  x = as.character(cluster_labels_filter[as.character(sce_integrated_sample_filter$seurat_clusters)]),
  levels = c("Progen", "Ground", "Prolif", "TR Diff Myo", "Diff Myo", "Fibro Perturb", "Peri", "SM Vasc"))

plot_cell_type_filter <- DimPlot(object = sce_integrated_sample_filter, reduction = "umap", label = TRUE, group.by = "cell_type", label.size = 3) +
  global_theme +
  scale_colour_manual(values = cluster_colors)

plot_cell_type_time_split_filter <- DimPlot(object = sce_integrated_sample_filter, reduction = "umap", label = FALSE, group.by = "cell_type", split.by = "Time", ncol = 2) +
  global_theme +
  scale_colour_manual(values = cluster_colors)

Idents(sce_integrated_sample_filter) <- sce_integrated_sample_filter$cell_type

# FindAllMarkers 2: cell-type batch markers (only.pos=TRUE)
cat("Running FindAllMarkers (batch, only.pos=TRUE) ...\n")
gc()
markers_batch_filter <- FindAllMarkers(object = sce_integrated_sample_filter, only.pos = TRUE)

top_markers_batch_filter <- markers_batch_filter %>%
  group_by(cluster) %>%
  slice_head(n = 35)

downsample_batch_filter <- min(table(sce_integrated_sample_filter$cell_type))
plot_heatmap_batch_filter <- safe_heatmap(
  object = subset(sce_integrated_sample_filter, downsample = downsample_batch_filter),
  features = top_markers_batch_filter$gene,
  raster = TRUE, angle = 45, draw.lines = FALSE, label = TRUE,
  group.colors = unlist(cluster_colors[levels(sce_integrated_sample_filter$cell_type)])) +
  viridis::scale_fill_viridis(option = "G")

plot_cell_dot_batch_kidney_filter <- safe_dotplot(object = sce_integrated_sample_filter, features = dot_plot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

plot_cell_dot_batch_mes_filter <- safe_dotplot(object = sce_integrated_sample_filter, features = dot_plot_genes_mes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

plot_features_filter <- safe_featureplot(object = sce_integrated_sample_filter,
                                    features = c("tdtomato", "Hic1", "Pdgfra", "Dpt", "Pi16", "Col15a1"),
                                    ncol = 3, order = TRUE)

plot_features_xenium_filter <- safe_featureplot(object = sce_integrated_sample_filter,
                                           features = c("Pi16", "Lum", "Col15a1", "Mfap4", "Tagln", "Myh11",
                                                        "Cenpf", "Mki67", "Car3", "Myoz2", "Itga8", "Rgs5"),
                                           ncol = 3, order = TRUE)

# FindAllMarkers 3: all DE (only.pos=FALSE) — most memory-intensive; wrapped in tryCatch
cat("Running FindAllMarkers (batch, only.pos=FALSE) ...\n")
gc()
all_de_batch_filter <- tryCatch(
  FindAllMarkers(object = sce_integrated_sample_filter, only.pos = FALSE),
  error = function(e) {
    cat("WARNING: all_de_batch_filter failed (likely OOM):", conditionMessage(e), "\n")
    data.frame()
  }
)

# Save the filtered, integrated Seurat object
cat("Saving filter_tumor.rds ...\n")
saveRDS(object = sce_integrated_sample_filter, file = "../output/scrna/dicer1_scrna_filter_tumor.rds")
cat("Saved filter_tumor.rds\n")

# Write marker tables
out_dir <- "../output/plots/"

write.csv(x = markers_cluster_filter, file = "../output/tables/dicer1_scrna_filter_tumor_cluster_markers_table.csv", row.names = FALSE)
write.csv(x = markers_batch_filter, file = "../output/tables/dicer1_scrna_filter_tumor_batch_markers_table.csv", row.names = FALSE)
write.csv(x = all_de_batch_filter, file = "../output/tables/dicer1_scrna_filter_tumor_batch_all_de_table.csv", row.names = FALSE)

# Save filtered plots
ggsave(filename = "dicer1_scrna_filter_tumor_clusters.png", plot = plot_clusters_filter, path = out_dir, width = 5.5, height = 4)
ggsave(filename = "dicer1_scrna_filter_tumor_clusters.svg", plot = plot_clusters_filter, path = out_dir, width = 5.5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_samples.png", plot = plot_samples_filter, path = out_dir, width = 6.5, height = 4)
ggsave(filename = "dicer1_scrna_filter_tumor_samples.svg", plot = plot_samples_filter, path = out_dir, width = 6.5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_cell_type.png", plot = plot_cell_type_filter, path = out_dir, width = 6, height = 4)
ggsave(filename = "dicer1_scrna_filter_tumor_cell_type.svg", plot = plot_cell_type_filter, path = out_dir, width = 6, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_time.png", plot = plot_time_filter, path = out_dir, width = 5, height = 4)
ggsave(filename = "dicer1_scrna_filter_tumor_time.svg", plot = plot_time_filter, path = out_dir, width = 5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_dot_plot_kidney.png", plot = plot_cell_dot_kidney_filter, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_filter_tumor_dot_plot_kidney.svg", plot = plot_cell_dot_kidney_filter, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_dot_plot_mes.png", plot = plot_cell_dot_mes_filter, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_filter_tumor_dot_plot_mes.svg", plot = plot_cell_dot_mes_filter, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_dot_plot_batch_kidney.png", plot = plot_cell_dot_batch_kidney_filter, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_filter_tumor_dot_plot_batch_kidney.svg", plot = plot_cell_dot_batch_kidney_filter, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_dot_plot_batch_mes.png", plot = plot_cell_dot_batch_mes_filter, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_filter_tumor_dot_plot_batch_mes.svg", plot = plot_cell_dot_batch_mes_filter, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_heatmap_clusters.png", plot = plot_heatmap_cluster_filter, path = out_dir, width = 15, height = 9)
ggsave(filename = "dicer1_scrna_filter_tumor_heatmap_clusters.svg", plot = plot_heatmap_cluster_filter, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_heatmap_batch.png", plot = plot_heatmap_batch_filter, path = out_dir, width = 15, height = 9)
ggsave(filename = "dicer1_scrna_filter_tumor_heatmap_batch.svg", plot = plot_heatmap_batch_filter, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_features.png", plot = plot_features_filter, path = out_dir, width = 18, height = 10)
ggsave(filename = "dicer1_scrna_filter_tumor_features.svg", plot = plot_features_filter, path = out_dir, width = 18, height = 10, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_xenium_features.png", plot = plot_features_xenium_filter, path = out_dir, width = 18, height = 20)
ggsave(filename = "dicer1_scrna_filter_tumor_xenium_features.svg", plot = plot_features_xenium_filter, path = out_dir, width = 18, height = 20, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_cell_type_time_split.png", plot = plot_cell_type_time_split_filter, path = out_dir, width = 10, height = 8)
ggsave(filename = "dicer1_scrna_filter_tumor_cell_type_time_split.svg", plot = plot_cell_type_time_split_filter, path = out_dir, width = 10, height = 8, device = "svg")

cat("All filtered tumor outputs complete.\n")
