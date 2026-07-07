# Script to pathway activity inference per cluster

# Required packages
library(Seurat)
library(decoupleR)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(patchwork)

# Read in the Seurat object and set cell identities to the desired cluster naming convention
rds_file <- "../output/scrna/dicer1_scrna_filter_tumor.rds"
seurat_integrated <- readRDS(file = rds_file)

ident_field <- "cell_type"
Idents(seurat_integrated) <- seurat_integrated[[ident_field]]

# Get the pathway information
pathways <- get_progeny(organism = "mouse", top = 500)

# Extract the normalized log-transformed counts
log_counts <- as.matrix(seurat_integrated@assays$RNA@data)

# Model the gene expression as a function of the pathway weights using multivariate linear modeling
pathway_activity <- run_mlm(mat = log_counts, network = pathways, .source = "source", .target = "target", .mor = "weight", minsize = 5)

# Store pathway activity as an assay in the Seurat object
seurat_integrated[['pathwaysmlm']] <- pathway_activity %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

DefaultAssay(object = seurat_integrated) <- "pathwaysmlm"
seurat_integrated <- ScaleData(object = seurat_integrated)
seurat_integrated@assays$pathwaysmlm@data <- seurat_integrated@assays$pathwaysmlm@scale.data

# Plot some activity data across the clusters
cluster_plot <- DimPlot(object = seurat_integrated, reduction = "umap", label = TRUE) +
  NoLegend() +
  ggtitle("Cell type")

tgfb_activity_plot <- FeaturePlot(object = seurat_integrated, features = c("TGFb")) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
  ggtitle("TGFb activity")

tnfa_activity_plot <- FeaturePlot(object = seurat_integrated, features = c("TNFa")) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
  ggtitle("TNFa activity")

nfkb_activity_plot <- FeaturePlot(object = seurat_integrated, features = c("NFkB")) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
  ggtitle("NFkB activity")

cluster_activity_plot <- cluster_plot | activity_plot

# Plot as a heat map
cluster_acitivity_df <- t(as.matrix(seurat_integrated@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(seurat_integrated)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

cluster_acitivity_df_wide <- cluster_acitivity_df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

palette_length = 100
my_colours = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

pheatmap(cluster_acitivity_df_wide, border_color = NA, color= my_colours, breaks = my_breaks)
