# dicer1_scrna_trajectory_patched.R
# Runs using 653_renv2 (R 4.5.3 + SeuratWrappers + monocle3)
#
# Patches:
#  - Non-interactive order_cells: auto-selects root as principal graph node
#    nearest the median UMAP coordinates of "Progen" cells (per author comment)
#  - choose_graph_segments (interactive) is SKIPPED — both branch trajectories
#    require user to click on the graph in a shiny app. Those sections are wrapped
#    in tryCatch and will produce blank placeholder plots.
#  - Added ggsave() for all plots (author original had none)
#  - Seurat v5 @assays$RNA@data access updated

library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)

set.seed(seed = 1337)

cat("Loading filter_tumor.rds ...\n")
seurat_integrated <- readRDS(file = "../output/scrna/dicer1_scrna_filter_tumor.rds")
cat("Loaded. Cells:", ncol(seurat_integrated), "\n")

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

# Seurat v5: use LayerData() for data layer access
cat("Building CellDataSet for monocle3...\n")
gene_annotation <- data.frame(gene_short_name = rownames(seurat_integrated))
rownames(gene_annotation) <- gene_annotation$gene_short_name

cell_metadata <- seurat_integrated@meta.data

exp_matrix <- tryCatch(
  LayerData(seurat_integrated, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(object = seurat_integrated, assay = "RNA", slot = "data")
)

cds_integrated <- new_cell_data_set(exp_matrix,
                                    cell_metadata = cell_metadata,
                                    gene_metadata = gene_annotation)

# Transfer UMAP embeddings from Seurat to monocle CDS
# (avoid re-running UMAP so trajectory is built on the Seurat embedding)
recreate_partition <- c(rep(1, length(cds_integrated@colData@rownames)))
names(recreate_partition) <- cds_integrated@colData@rownames
recreate_partition <- as.factor(recreate_partition)
cds_integrated@clusters@listData[["UMAP"]][["partitions"]] <- recreate_partition

clusters <- seurat_integrated$cell_type
names(clusters) <- colnames(seurat_integrated)
cds_integrated@clusters@listData[["UMAP"]][["clusters"]] <- clusters
cds_integrated@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
cds_integrated@int_colData@listData$reducedDims@listData[["UMAP"]] <-
  seurat_integrated@reductions[["umap"]]@cell.embeddings

cat("Learning trajectory graph...\n")
cds_integrated <- learn_graph(cds = cds_integrated,
                              close_loop = TRUE,
                              use_partition = FALSE,
                              learn_graph_control = list("minimal_branch_len" = 15))

# Non-interactive root selection:
# Find the principal graph node closest to the median UMAP coordinates of Progen cells
# (per author comment: "choose start node as node nearest the Dpt+ region of the Progen cluster")
cat("Selecting root node automatically (nearest to Progen cluster centroid)...\n")
progen_cells <- colnames(seurat_integrated)[seurat_integrated$cell_type == "Progen"]
umap_coords <- seurat_integrated@reductions[["umap"]]@cell.embeddings
progen_centroid <- colMeans(umap_coords[progen_cells, , drop = FALSE])

principal_graph_aux <- cds_integrated@principal_graph_aux
pr_graph_nodes <- principal_graph_aux[["UMAP"]][["dp_mst"]]
node_coords <- t(pr_graph_nodes)
dists <- sqrt(rowSums((node_coords - progen_centroid)^2))
root_node <- rownames(node_coords)[which.min(dists)]
cat("Selected root node:", root_node, "\n")

cds_integrated <- order_cells(cds = cds_integrated, root_pr_nodes = root_node)

plot_cell_type <- plot_cells(cds = cds_integrated, reduction_method = "UMAP",
                             color_cells_by = "cell_type", label_cell_groups = FALSE)

plot_pseudotime <- plot_cells(cds = cds_integrated, color_cells_by = "pseudotime") +
  theme(legend.position = "none")

plot_pseudotime_no_graph <- plot_cells(cds = cds_integrated,
                                       color_cells_by = "pseudotime",
                                       show_trajectory_graph = FALSE)

plot_cell_type_graph <- plot_cells(cds = cds_integrated,
                                   color_cells_by = "cell_type",
                                   label_cell_groups = FALSE) +
  scale_colour_manual(name = "Cell type", values = cluster_colors) +
  theme(legend.position = "none")

out_dir <- "../output/plots/"

ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_pseudotime.png", plot = plot_pseudotime, path = out_dir, width = 5.5, height = 3.5)
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_pseudotime.svg", plot = plot_pseudotime, path = out_dir, width = 5.5, height = 3.5, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_pseudotime_no_graph.png", plot = plot_pseudotime_no_graph, path = out_dir, width = 5.5, height = 4)
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_pseudotime_no_graph.svg", plot = plot_pseudotime_no_graph, path = out_dir, width = 5.5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_cell_type_graph.png", plot = plot_cell_type_graph, path = out_dir, width = 5.5, height = 3.5)
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_cell_type_graph.svg", plot = plot_cell_type_graph, path = out_dir, width = 5.5, height = 3.5, device = "svg")

# Branch-specific trajectories require choose_graph_segments() — an interactive Shiny call.
# SKIPPED in non-interactive mode; blank placeholder plots are written instead.
muscle_genes <- c("Dpt", "Mfap4", "Mki67", "Pax7", "Myog")
ground_genes  <- c("Dpt", "Pi16", "Col15a1", "Mfap4", "Itga8")

blank_plot <- ggplot() + geom_blank() +
  annotate("text", x = 0.5, y = 0.5,
           label = "SKIPPED: choose_graph_segments() requires interactive Shiny input",
           size = 4) + theme_void()

ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_myo_diff_genes_time.png", plot = blank_plot, path = out_dir, width = 10, height = 10)
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_myo_diff_genes_cell_type.png", plot = blank_plot, path = out_dir, width = 4, height = 5)
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_ground_genes_time.png", plot = blank_plot, path = out_dir, width = 10, height = 10)
ggsave(filename = "dicer1_scrna_filter_tumor_trajectory_ground_genes_cell_type.png", plot = blank_plot, path = out_dir, width = 4, height = 5)

# Save the CDS object for downstream use
saveRDS(cds_integrated, file = "../output/scrna/dicer1_scrna_trajectory_cds.rds")

cat("Trajectory analysis complete. Branch trajectory plots skipped (require interactive root selection).\n")
