# dicer1_scrna_run_pathway_analysis_patched.R
# Patches:
#  - Fixed undefined 'activity_plot' variable (author bug)
#  - Added ggsave() and pdf() calls for all plots (original sent all to display device)
#  - Seurat v5 assay data access via LayerData() with slot= fallback
#  - CreateAssayObject replaced with CreateAssay5Object for v5 compatibility
#  - BYPASS OmnipathR: load PROGENy data from local TSV downloaded from omnipathdb.org API
#    (decoupleR 2.4.0 get_progeny() requires OmnipathR which is not installed)
#  - Explicit Idents assignment via $ rather than [[]] to avoid data.frame coercion error

library(Seurat)
library(decoupleR)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(patchwork)
library(stringr)

rds_file <- "../output/scrna/dicer1_scrna_filter_tumor.rds"
seurat_integrated <- readRDS(file = rds_file)
cat("Loaded filter_tumor.rds. Cells:", ncol(seurat_integrated), "\n")

Idents(seurat_integrated) <- seurat_integrated$cell_type  # PATCH: [[]] returns df, $ returns vector

# PATCH: local get_progeny() replacement — reads from pre-downloaded OmniPath TSV
# DOC: OmniPath REST API https://omnipathdb.org/annotations?resources=PROGENy&format=tsv
# Downloaded to ../data/progeny_omnipath.tsv  (record_id groups weight+p_value+pathway per gene)
local_get_progeny <- function(organism = "mouse", top = 500) {
  tsv_path <- "../data/progeny_omnipath.tsv"
  raw <- read.table(tsv_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  # pivot: one row per (record_id, genesymbol); spread label -> value
  wide <- raw %>%
    dplyr::select(record_id, genesymbol, label, value) %>%
    tidyr::pivot_wider(names_from = label, values_from = value) %>%
    dplyr::mutate(
      weight   = as.double(weight),
      p_value  = as.double(p_value)
    ) %>%
    dplyr::rename(source = pathway) %>%
    dplyr::distinct(source, genesymbol, .keep_all = TRUE)
  # Apply mouse sentence case (same as get_progeny)
  if (tolower(organism) == "mouse") {
    wide$genesymbol <- stringr::str_to_sentence(wide$genesymbol)
  }
  # Take top N genes per pathway by lowest p_value
  result <- wide %>%
    dplyr::group_by(source) %>%
    dplyr::arrange(p_value) %>%
    dplyr::slice_head(n = top) %>%
    dplyr::ungroup() %>%
    dplyr::select(source, genesymbol, weight, p_value)
  colnames(result) <- c("source", "target", "weight", "p_value")
  return(result)
}

cat("Loading PROGENy mouse pathways from local TSV...\n")
pathways <- local_get_progeny(organism = "mouse", top = 500)
cat("PROGENy pathways loaded. Rows:", nrow(pathways), "Pathways:", length(unique(pathways$source)), "\n")

# Seurat v5 compatible count data extraction
cat("Extracting normalized log counts...\n")
log_counts <- tryCatch(
  as.matrix(LayerData(seurat_integrated, assay = "RNA", layer = "data")),
  error = function(e) as.matrix(GetAssayData(object = seurat_integrated, assay = "RNA", slot = "data"))
)

cat("Running multivariate linear model (MLM)...\n")
pathway_activity <- run_mlm(mat = log_counts, network = pathways,
                             .source = "source", .target = "target",
                             .mor = "weight", minsize = 5)

rm(log_counts)
gc()

# Build wide pathway activity matrix
pathway_mat <- pathway_activity %>%
  pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  column_to_rownames('source')

# Store as Seurat assay — use CreateAssay5Object if available (Seurat v5), else CreateAssayObject
# PATCH v3: Seurat v5 Assay5 requires data layer set explicitly before ScaleData
cat("Storing pathway activity as Seurat assay...\n")
pathway_assay_obj <- tryCatch(
  CreateAssay5Object(pathway_mat),
  error = function(e) CreateAssayObject(counts = pathway_mat, data = pathway_mat)
)
seurat_integrated[['pathwaysmlm']] <- pathway_assay_obj

# For Seurat v5 Assay5: set data layer = pathway_mat (activity scores, no normalization needed)
tryCatch({
  seurat_integrated <- SetAssayData(seurat_integrated, assay = "pathwaysmlm",
                                    layer = "data", new.data = pathway_mat)
  cat("Data layer set via SetAssayData (Seurat v5 path)\n")
}, error = function(e) {
  cat("SetAssayData skipped (Seurat v4 path):", conditionMessage(e), "\n")
})

DefaultAssay(object = seurat_integrated) <- "pathwaysmlm"
seurat_integrated <- ScaleData(object = seurat_integrated)

# Store scaled data back into @data slot (safe for both v4 and v5)
tryCatch({
  seurat_integrated@assays$pathwaysmlm@data <- seurat_integrated@assays$pathwaysmlm@scale.data
  cat("scale.data copied to data slot (Seurat v4 path)\n")
}, error = function(e) {
  # Seurat v5 Assay5: use SetAssayData to set the data layer from scale.data
  scaled <- tryCatch(
    LayerData(seurat_integrated, assay = "pathwaysmlm", layer = "scale.data"),
    error = function(e2) LayerData(seurat_integrated, assay = "pathwaysmlm", layer = "data")
  )
  seurat_integrated <<- SetAssayData(seurat_integrated, assay = "pathwaysmlm",
                                     layer = "data", new.data = scaled)
  cat("scale.data set via SetAssayData (Seurat v5 path)\n")
})

# UMAP cluster plot
cluster_plot <- DimPlot(object = seurat_integrated, reduction = "umap", label = TRUE) +
  NoLegend() + ggtitle("Cell type")

# Feature plots for individual pathways
# PATCH: FeaturePlot with S7 Seurat objects may fail on & operator internally
# Wrap in tryCatch and fall back to manual ggplot overlay if needed
safe_featureplot_pathway <- function(seurat_obj, feature, title) {
  # Get UMAP coords and pathway score into a data frame
  umap_df <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  score_vec <- tryCatch(
    FetchData(seurat_obj, vars = feature)[, 1],
    error = function(e) rep(0, ncol(seurat_obj))
  )
  umap_df$score <- score_vec
  ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = score)) +
    geom_point(size = 0.3) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",
                           name = feature) +
    ggtitle(title) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
}

tgfb_activity_plot <- tryCatch(
  FeaturePlot(object = seurat_integrated, features = c("TGFb")) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
    ggtitle("TGFb activity"),
  error = function(e) {
    cat("FeaturePlot fallback for TGFb:", conditionMessage(e), "\n")
    safe_featureplot_pathway(seurat_integrated, "TGFb", "TGFb activity")
  }
)

tnfa_activity_plot <- tryCatch(
  FeaturePlot(object = seurat_integrated, features = c("TNFa")) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
    ggtitle("TNFa activity"),
  error = function(e) {
    cat("FeaturePlot fallback for TNFa:", conditionMessage(e), "\n")
    safe_featureplot_pathway(seurat_integrated, "TNFa", "TNFa activity")
  }
)

nfkb_activity_plot <- tryCatch(
  FeaturePlot(object = seurat_integrated, features = c("NFkB")) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
    ggtitle("NFkB activity"),
  error = function(e) {
    cat("FeaturePlot fallback for NFkB:", conditionMessage(e), "\n")
    safe_featureplot_pathway(seurat_integrated, "NFkB", "NFkB activity")
  }
)

# PATCH: author bug — 'activity_plot' was never defined; use tgfb as default
cluster_tgfb_plot <- cluster_plot | tgfb_activity_plot

# Heatmap of mean pathway activity per cluster
# PATCH: use LayerData for Seurat v5 Assay5 compatibility
pathway_data_mat <- tryCatch(
  as.matrix(LayerData(seurat_integrated, assay = "pathwaysmlm", layer = "data")),
  error = function(e) as.matrix(seurat_integrated@assays$pathwaysmlm@data)
)
cluster_acitivity_df <- t(pathway_data_mat) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(seurat_integrated)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score), .groups = "drop")

cluster_acitivity_df_wide <- cluster_acitivity_df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source', values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

palette_length <- 100
my_colours <- colorRampPalette(c("Darkblue", "white", "red"))(palette_length)
my_breaks <- c(seq(-3, 0, length.out = ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out = floor(palette_length/2)))

out_dir <- "../output/plots/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(filename = "dicer1_scrna_filter_tumor_pathway_tgfb.png", plot = tgfb_activity_plot, path = out_dir, width = 5, height = 4)
ggsave(filename = "dicer1_scrna_filter_tumor_pathway_tgfb.svg", plot = tgfb_activity_plot, path = out_dir, width = 5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_pathway_tnfa.png", plot = tnfa_activity_plot, path = out_dir, width = 5, height = 4)
ggsave(filename = "dicer1_scrna_filter_tumor_pathway_tnfa.svg", plot = tnfa_activity_plot, path = out_dir, width = 5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_pathway_nfkb.png", plot = nfkb_activity_plot, path = out_dir, width = 5, height = 4)
ggsave(filename = "dicer1_scrna_filter_tumor_pathway_nfkb.svg", plot = nfkb_activity_plot, path = out_dir, width = 5, height = 4, device = "svg")
ggsave(filename = "dicer1_scrna_filter_tumor_pathway_cluster_tgfb.png", plot = cluster_tgfb_plot, path = out_dir, width = 10, height = 4)
ggsave(filename = "dicer1_scrna_filter_tumor_pathway_cluster_tgfb.svg", plot = cluster_tgfb_plot, path = out_dir, width = 10, height = 4, device = "svg")

# Save pheatmap as PNG (author original only sent to display device)
png(filename = paste0(out_dir, "dicer1_scrna_filter_tumor_pathway_heatmap.png"), width = 2400, height = 1800, res = 200)
pheatmap(cluster_acitivity_df_wide, border_color = NA, color = my_colours, breaks = my_breaks)
dev.off()

svg(filename = paste0(out_dir, "dicer1_scrna_filter_tumor_pathway_heatmap.svg"), width = 12, height = 9)
pheatmap(cluster_acitivity_df_wide, border_color = NA, color = my_colours, breaks = my_breaks)
dev.off()

# Save pathway activity table
dir.create("../output/tables/", recursive = TRUE, showWarnings = FALSE)
write.csv(cluster_acitivity_df_wide, file = "../output/tables/dicer1_scrna_filter_tumor_pathway_activity.csv", row.names = TRUE)

cat("Pathway analysis complete.\n")
