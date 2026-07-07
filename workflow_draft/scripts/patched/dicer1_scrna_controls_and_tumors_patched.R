# dicer1_scrna_controls_and_tumors_patched.R
# Patches: JoinLayers (Seurat v5), safe_heatmap/safe_dotplot/safe_featureplot wrappers,
#          source dotplot_compat_patch.R for ggplot2 ≥3.4 compatibility

source("/home/node/MAURICIO/papers/653_Spatial_single_cell_transcriptomic_analysis_informs_tumor_developmental_hierarchy_of_DICER1_syndrome_related_sarcoma_2026/REPRODUCIBILITY_ANALYSIS/20260703T220654Z__anthropic-claude-sonnet-4-6__adaptive/workflow_draft/scripts/patched/dotplot_compat_patch.R")

### CONTROLS + MUTANTS
library(dplyr)
library(Seurat)
library(scater)
library(cowplot)
library(ggplot2)
library(yaml)

set.seed(1337)

### FILTERED
seurat_control_filter <- readRDS(file = "../output/scrna/dicer1_scrna_filter_control.rds")
seurat_mutant_filter  <- readRDS(file = "../output/scrna/dicer1_scrna_filter_tumor.rds")

seurat_list_control_filter <- SplitObject(seurat_control_filter, split.by = "Sample")
seurat_list_mutant_filter  <- SplitObject(seurat_mutant_filter,  split.by = "Sample")

seurat_list_sample_filter <- c(seurat_list_control_filter, seurat_list_mutant_filter)
seurat_list_sample_filter <- lapply(seurat_list_sample_filter, function(seurat_obj) NormalizeData(seurat_obj, verbose = FALSE))
seurat_list_sample_filter <- lapply(seurat_list_sample_filter, function(seurat_obj) FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000))

features_sample_filter <- SelectIntegrationFeatures(object.list = seurat_list_sample_filter, nfeatures = 2000)

seurat_list_sample_filter <- lapply(seurat_list_sample_filter, function(seurat_obj) ScaleData(seurat_obj, features = features_sample_filter))
seurat_list_sample_filter <- lapply(seurat_list_sample_filter, function(seurat_obj) RunPCA(seurat_obj, features = features_sample_filter))

anchors_sample_filter <- FindIntegrationAnchors(object.list = seurat_list_sample_filter, anchor.features = features_sample_filter, reduction = "rpca")

sce_integrated_sample_filter <- IntegrateData(anchorset = anchors_sample_filter, k.weight = 50)

rm(seurat_list_sample_filter, anchors_sample_filter, features_sample_filter, seurat_control_filter, seurat_mutant_filter,
   seurat_list_control_filter, seurat_list_mutant_filter)
gc()

sce_integrated_sample_filter <- ScaleData(sce_integrated_sample_filter, do.scale = TRUE, do.center = TRUE, assay = "integrated")
sce_integrated_sample_filter <- RunPCA(sce_integrated_sample_filter, npcs = 50, assay = "integrated", seed.use = 1337)

plot_elbow_filter <- ElbowPlot(sce_integrated_sample_filter, ndims = 50)

sce_integrated_sample_filter <- RunTSNE(sce_integrated_sample_filter, dims = 1:30, assay = "integrated", dim.embed = 3, seed.use = 300, perplexity = 20)
sce_integrated_sample_filter <- RunUMAP(sce_integrated_sample_filter, reduction = "pca", dims = 1:30, assay = "integrated", seed.use = 1337)
sce_integrated_sample_filter <- FindNeighbors(sce_integrated_sample_filter, reduction = "pca", dims = 1:30, assay = "integrated")
sce_integrated_sample_filter <- FindClusters(sce_integrated_sample_filter, resolution = 0.2, random.seed = 1337)

# PATCH: JoinLayers required in Seurat v5 before FindAllMarkers/ScaleData on RNA assay
DefaultAssay(object = sce_integrated_sample_filter) <- "RNA"
sce_integrated_sample_filter <- JoinLayers(sce_integrated_sample_filter)
sce_integrated_sample_filter <- ScaleData(object = sce_integrated_sample_filter)

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
  "Prox Tub" = "#B2A000", "Loop Henle" = "#00C091", "Endo" = "#F8766D",
  "Immune" = "#B186FF", "Fibro" = "#CE9500", "Collect Duct" = "#F27D53",
  "Trans Epi" = "#FF6A98", "Uni Fibro" = "#00BAE0", "Podo" = "#D973FC",
  "Ground" = "#00BB4E", "Progen" = "#C77CFF", "Diff Myo" = "#00C1A3",
  "Prolif" = "#FF61CC", "Diff Fibro" = "#A3A500", "TR Diff Myo" = "#7099FF",
  "Mural" = "#39B606", "Uni Fibro (Pi16)" = "#FD7083", "Fibro Perturb" = "#EA8331",
  "Uni Fibro (Col15a1)" = "#B186FF", "Cycle" = "#00BAE0", "Mac Densa" = "#39B607",
  "Peri" = "#C09B06", "SM Vasc" = "#91AA00", "Unk" = "black"
)

plot_clusters_filter <- DimPlot(object = sce_integrated_sample_filter, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE) +
  global_theme + scale_colour_manual(values = colors)

plot_samples_filter <- DimPlot(object = sce_integrated_sample_filter, reduction = "umap", group.by = "Name", shuffle = TRUE) +
  global_theme + scale_colour_manual(values = colors_secondary)

plot_time_filter <- DimPlot(object = sce_integrated_sample_filter, reduction = "umap", group.by = "Time", shuffle = TRUE) +
  global_theme + scale_colour_manual(values = colors_secondary)

dot_plot_genes     <- read_yaml(file = "../reference/cell_type_markers_mus_musculus_kidney.yaml")
dot_plot_genes_mes <- read_yaml(file = "../reference/cell_type_markers_mus_musculus_mesenchymal.yaml")

plot_cell_dot_kidney_filter <- safe_dotplot(object = sce_integrated_sample_filter, features = dot_plot_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

plot_cell_dot_mes_filter <- safe_dotplot(object = sce_integrated_sample_filter, features = dot_plot_genes_mes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

cat("Running FindAllMarkers (cluster) ...\n")
gc()
markers_cluster_filter <- FindAllMarkers(object = sce_integrated_sample_filter, only.pos = TRUE)

top_markers_cluster_filter <- markers_cluster_filter %>%
  group_by(cluster) %>%
  slice_head(n = 10)

downsample_cluster_filter <- min(table(sce_integrated_sample_filter$seurat_clusters))
plot_heatmap_cluster_filter <- safe_heatmap(
  object = subset(sce_integrated_sample_filter, downsample = downsample_cluster_filter),
  features = top_markers_cluster_filter$gene,
  raster = TRUE, angle = 45, draw.lines = FALSE, label = TRUE) +
  viridis::scale_fill_viridis(option = "G")

plot_features_filter <- safe_featureplot(object = sce_integrated_sample_filter,
                                    features = c("tdtomato", "Hic1", "Pdgfra", "Dpt", "Pi16", "Col15a1"),
                                    ncol = 3, order = TRUE, min.cutoff = 0)

plot_features_xenium_filter <- safe_featureplot(object = sce_integrated_sample_filter,
                                           features = c("Pi16", "Lum", "Col15a1", "Mfap4", "Tagln", "Myh11",
                                                        "Cenpf", "Mki67", "Car3", "Myoz2", "Itga8", "Rgs5"),
                                           ncol = 3, order = TRUE, min.cutoff = 0)

# Plot old cell type annotations (per-dataset, before combined integration)
sce_integrated_sample_filter$old_cell_type <- factor(sce_integrated_sample_filter$cell_type,
  levels = c("Uni Fibro (Col15a1)", "Progen", "Fibro", "Ground", "Cycle", "Prolif",
             "TR Diff Myo", "Diff Myo", "Uni Fibro (Pi16)", "Fibro Perturb", "Mural", "Peri", "SM Vasc"))

plot_cell_type_filter <- DimPlot(object = sce_integrated_sample_filter, reduction = "umap", label = TRUE,
                                 group.by = "old_cell_type", label.size = 3, shuffle = TRUE) +
  global_theme + scale_colour_manual(values = colors)

plot_genotype_filter <- DimPlot(object = sce_integrated_sample_filter, reduction = "umap",
                                group.by = "Genotype", shuffle = TRUE) +
  global_theme + scale_colour_manual(values = colors_secondary)

saveRDS(object = sce_integrated_sample_filter, file = "../output/scrna/dicer1_scrna_filter_control_and_tumor.rds")

write.csv(x = markers_cluster_filter, file = "../output/tables/dicer1_scrna_filter_control_and_tumor_cluster_markers_table.csv", row.names = FALSE)

out_dir <- "../output/plots/"

ggsave(filename = "dicer1_scrna_filter_control_and_tumor_clusters.png", plot = plot_clusters_filter, path = out_dir, width = 5.5, height = 4)
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_clusters.svg", plot = plot_clusters_filter, path = out_dir, width = 5.5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_samples.png", plot = plot_samples_filter, path = out_dir, width = 6.5, height = 4)
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_samples.svg", plot = plot_samples_filter, path = out_dir, width = 6.5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_cell_type.png", plot = plot_cell_type_filter, path = out_dir, width = 6, height = 4)
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_cell_type.svg", plot = plot_cell_type_filter, path = out_dir, width = 6, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_time.png", plot = plot_time_filter, path = out_dir, width = 5, height = 4)
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_time.svg", plot = plot_time_filter, path = out_dir, width = 5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_genotype.png", plot = plot_genotype_filter, path = out_dir, width = 5, height = 4)
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_dot_plot_kidney.png", plot = plot_cell_dot_kidney_filter, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_dot_plot_kidney.svg", plot = plot_cell_dot_kidney_filter, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_dot_plot_mes.png", plot = plot_cell_dot_mes_filter, path = out_dir, width = 15, height = 9, bg = "white")
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_dot_plot_mes.svg", plot = plot_cell_dot_mes_filter, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_heatmap_clusters.png", plot = plot_heatmap_cluster_filter, path = out_dir, width = 15, height = 9)
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_heatmap_clusters.svg", plot = plot_heatmap_cluster_filter, path = out_dir, width = 15, height = 9, device = "svg")
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_features.png", plot = plot_features_filter, path = out_dir, width = 18, height = 10)
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_features.svg", plot = plot_features_filter, path = out_dir, width = 18, height = 10, device = "svg")
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_xenium_features.png", plot = plot_features_xenium_filter, path = out_dir, width = 18, height = 20)
ggsave(filename = "dicer1_scrna_filter_control_and_tumor_xenium_features.svg", plot = plot_features_xenium_filter, path = out_dir, width = 18, height = 20, device = "svg")

cat("All controls+tumors outputs complete.\n")
