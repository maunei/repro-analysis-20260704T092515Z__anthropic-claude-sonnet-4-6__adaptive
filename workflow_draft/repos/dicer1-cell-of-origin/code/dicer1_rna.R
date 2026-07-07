# Script to do DE with DESeq2, followed by enrichment

# Required packages
library(tximport)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ComplexHeatmap)

# Read in metadata and data
metadata <- read.delim(file = "../metadata/dicer1_rna_metadata.csv",
                       sep = ",") %>%
  mutate(genotype = factor(genotype),
         media = factor(media),
         geno_media = factor(paste0(genotype,
                                    "_",
                                    media),
                             levels = c("control_mesencult",
                                        "control_dmem",
                                        "mutant_mesencult",
                                        "mutant_dmem")))
  
tx2gene <- read.delim(file = "../data/star_salmon/tx2gene.tsv", 
                      sep = "\t",
                      header = FALSE)

count_files <- paste0("../data/star_salmon/", 
                      metadata$sample, 
                      "/quant.sf")

counts <- tximport(files = count_files, 
                   type = "salmon", 
                   txIn = TRUE, 
                   countsFromAbundance = "no", 
                   tx2gene = tx2gene)

# Build our DESeq objects
dds <- DESeqDataSetFromTximport(txi = counts, 
                                colData = metadata, 
                                design = ~ geno_media)

keep_genes <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep_genes, ]
dds <- DESeq(object = dds)

result_mesencult <- results(object = dds,
                            contrast = c("geno_media", 
                                         "mutant_mesencult", 
                                         "control_mesencult"))

result_dmem <- results(object = dds,
                       contrast = c("geno_media",
                                    "mutant_dmem",
                                    "control_dmem"))

result_mesencult_gene_name <- result_mesencult
result_dmem_gene_name <- result_dmem

rownames(result_mesencult_gene_name) <- tx2gene$V3[match(rownames(result_mesencult_gene_name), tx2gene$V2)]
rownames(result_dmem_gene_name) <- tx2gene$V3[match(rownames(result_dmem_gene_name), tx2gene$V2)]

# Normalised counts for export
norm_counts <- counts(dds,
                      normalized = TRUE) %>%
  as.data.frame()

names(norm_counts) <- colData(dds)$sample
norm_counts <- norm_counts %>%
  tibble::rownames_to_column(var = "ensembl_id") %>%
  left_join(y = tx2gene,
            by = join_by(ensembl_id == V2),
            multiple = "first") %>%
  dplyr::select(!V1) %>%
  rename(gene_name = V3) %>%
  relocate(gene_name)
  
# Quick PCA plot to visualize variance w/r/t genotype and media type
vst_norm <- vst(object = dds)
pca_plot <- plotPCA(object = vst_norm,
                    intgroup = "geno_media") +
  scale_colour_manual(name = "Genotype & condition",
                    labels = c("control_mesencult" = "Control, Mesencult",
                               "control_dmem" = "Control, DMEM",
                               "mutant_mesencult" = "Mutant, Mesencult",
                               "mutant_dmem" = "Mutant, DMEM"),
                    values = c("red",
                    "pink",
                    "blue",
                    "lightblue")) +
  theme_classic(base_size = 18)

# Select marker genes from the spatial transcriptomics and scRNA-seq and visualize
control_cluster_markers <- read.csv( "../output/tables/seurat_filter_control_batch_markers_table.csv")

top_markers <- control_cluster_markers %>%
  filter(gene %in% rownames(result_mesencult_gene_name)) %>% # Only keep genes that are actually expressed in the bulk data
  group_by(cluster) %>%
  slice_head(n = 50)

vst_counts <- vst(dds)
vst_mat <- assay(vst_counts)
rownames(vst_mat) <- tx2gene$V3[match(rownames(vst_mat), tx2gene$V2)]
vst_mat <- vst_mat[top_markers$gene, ]
vst_mat_scale <- t(scale(t(vst_mat)))


col_annot_df <- data.frame(Genotype = vst_counts$genotype,
                    Media = vst_counts$media)
col_annot <- HeatmapAnnotation(df = col_annot_df,
                           which = "col")

colours <- list(
  Cluster = c('Cycle' = 'blue', 'Fibro' = 'red', 'Mural' = 'green3', 'Uni Fibro (Col15a1)' = 'yellow', 'Uni Fibro (Pi16)' = 'orange'))

row_annot_df <- data.frame(Cluster = top_markers$cluster)
row_annot <- HeatmapAnnotation(df = row_annot_df,
                               which = "row",
                               col = colours)

heatmap <- Heatmap(matrix = vst_mat_scale,
                   top_annotation = col_annot,
                   right_annotation = row_annot,
                   cluster_rows = FALSE,
                   show_row_names = TRUE)

# Dot plot of genes of interest
goi <- factor(x = c("Dpt",
         "Pi16",
         "Dcn",
         "Ly6a",
         "Ly6c1",
         "Cfh",
         "Col15a1",
         "Myh11",
         "Acta2",
         "Tagln"),
         levels = c("Dpt",
                    "Dcn",
                    "Pi16",
                    "Ly6a",
                    "Ly6c1",
                    "Cfh",
                    "Col15a1",
                    "Itga8",
                    "Myh11",
                    "Acta2",
                    "Tagln"))

goi_filter <- goi[goi %in% rownames(result_mesencult_gene_name)] # Keep genes that are expressed in the bulk data

result_mesencult_goi_df <- as.data.frame(result_mesencult_gene_name) %>%
  tibble::rownames_to_column(var = "gene_name") %>%
  filter(gene_name %in% goi_filter) %>%
  mutate(condition = "Mesencult")

result_dmem_goi_df <- as.data.frame(result_dmem_gene_name) %>%
  tibble::rownames_to_column(var = "gene_name") %>%
  filter(gene_name %in% goi_filter) %>%
  mutate(condition = "DMEM")

result_goi_df <- bind_rows(result_mesencult_goi_df,
                           result_dmem_goi_df) %>%
  mutate(sig = if_else(padj < 0.05,
                       TRUE,
                       FALSE),
         gene_name = factor(x = gene_name,
                            levels = levels(goi)),
         group = case_when(gene_name %in% c("Dpt", "Dcn") ~ "Universal",
                           gene_name %in% c("Pi16", "Ly6a", "Ly6c1") ~ "Universal (Pi16)",
                           gene_name %in% c("Cfh", "Col15a1") ~ "Universal (Col15a1)",
                           gene_name %in% c("Itga8") ~ "Fibro",
                           gene_name %in% c("Myh11", "Acta2", "Tagln") ~ "Peri"),
         group = factor(x = group,
                        levels = c("Universal",
                                   "Universal (Pi16)",
                                   "Universal (Col15a1)",
                                   "Fibro",
                                   "Peri")))

dot_plot_all <- ggplot(data = as.data.frame(result_goi_df),
                             aes(x = log2FoldChange,
                                 y = gene_name,
                                 size = -log10(padj),
                                 colour = condition,
                                 shape = sig)) +
  geom_point() +
  scale_shape_manual(values = c("TRUE" = 16,
                                "FALSE" = 1)) +
  theme_classic(base_size = 18) +
  labs(title = "Dicer1 mutant vs. control",
       x = "log2(fold-change)",
       y = "Gene",
       col = "Condition",
       shape = "Significant (padj)") +
  xlim(-9, 9) +
  coord_flip() +
  facet_grid(~ group,
             scales = "free_x",
             space = "free",
             switch = "x") +
  theme(strip.placement = "outside",
    strip.background = element_blank(),
    panel.spacing.x = unit(0, "cm"),
    axis.ticks.x = element_blank())

# Generate volcano plots for DE
volcano_mesencult <- EnhancedVolcano(result_mesencult_gene_name, 
                                     lab = rownames(result_mesencult_gene_name), 
                                     x = "log2FoldChange", 
                                     y = "padj", 
                                     pCutoff = 0.05,
                                     legendPosition = "right") +
  labs(title = "Dicer1 -/D1693N vs. Dicer1 -/+",
       subtitle = "Mesencult")

volcano_dmem <- EnhancedVolcano(result_dmem_gene_name, 
                                     lab = rownames(result_dmem_gene_name), 
                                     x = "log2FoldChange", 
                                     y = "padj", 
                                     pCutoff = 0.05,
                                     legendPosition = "right") +
  labs(title = "Dicer1 -/D1693N vs. Dicer1 -/+",
       subtitle = "Dmem")

# GSEA and some related functions
doGsea <- function(df, gene_sets, max_sets = 25, alpha = 0.05, ...) {
  gene_rank <- df$log2FoldChange
  names(gene_rank) <- rownames(df)
  
  gsea_df <- fgsea(pathways = gene_sets, stats = gene_rank, ...) %>%
    arrange(padj)
  
  sig_pathways <- gene_sets[gsea_df$pathway[gsea_df$padj < alpha]][1:max_sets]
  gsea_plot <- plotGseaTable(pathways = sig_pathways, stats = gene_rank, fgseaRes = gsea_df)
  
  return(list("table" = gsea_df, "plot" = gsea_plot))
}

# MSigDB gene sets
msigdb_hallmarks <- msigdbr(db_species = "MM",
                            species = "Mus musculus", 
                            category = "MH")
msigdb_hallmarks_split = split(x = msigdb_hallmarks$gene_symbol, 
                               f = msigdb_hallmarks$gs_name)
msigdb_oncogenic <- msigdbr(species = "Mus musculus", 
                            category = "C6")
msigdb_oncogenic_split <- split(x = msigdb_oncogenic$gene_symbol, 
                                f = msigdb_oncogenic$gs_name)
msigdb_tf <- msigdbr(species = "Mus musculus", 
                     category = "C3", 
                     subcategory = "TFT:GTRD")
msigdb_tf_split <- split(x = msigdb_tf$gene_symbol, 
                         f = msigdb_tf$gs_name)
msigdb_reactome <- msigdbr(species = "Mus musculus", 
                           category = "C2", 
                           subcategory = "CP:REACTOME")
msigdb_reactome_split <- split(x = msigdb_reactome$gene_symbol, 
                               f = msigdb_reactome$gs_name)
msigdb_mir <- msigdbr(db_species = "MM",
                      species = "Mus musculus", 
                      category = "M3")
msigdb_mir_split <- split(x = msigdb_mir$gene_symbol, 
                          f = msigdb_mir$gs_name)
msigdb_ont <- msigdbr(db_species = "MM",
                      species = "Mus musculus", 
                      category = "M5",
                      subcollection = "GO:BP")
msigdb_ont_split <- split(x = msigdb_ont$gene_symbol, 
                          f = msigdb_ont$gs_name)
msigdb_ont_muscle_split <- msigdb_ont_split[grepl(pattern = "MUSCLE", x = names(msigdb_ont_split))]

# Results & plotting, GSEA
hallmarks_mesencult <- doGsea(df = as.data.frame(result_mesencult_gene_name),
                              gene_sets = msigdb_hallmarks_split,
                              alpha = 0.05,
                              nPermSimple = 10000)
onocgenic_mesencult <- doGsea(df = as.data.frame(result_mesencult_gene_name),
                              gene_sets = msigdb_oncogenic_split,
                              alpha = 0.05,
                              nPermSimple = 10000)
tf_mesencult <- doGsea(df = as.data.frame(result_mesencult_gene_name),
                              gene_sets = msigdb_tf_split,
                              alpha = 0.05,
                              nPermSimple = 10000)
reactome_mesencult <- doGsea(df = as.data.frame(result_mesencult_gene_name),
                       gene_sets = msigdb_reactome_split,
                       alpha = 0.05,
                       nPermSimple = 10000)
mir_mesencult <- doGsea(df = as.data.frame(result_mesencult_gene_name),
                       gene_sets = msigdb_mir_split,
                       alpha = 0.05,
                       nPermSimple = 10000)
ont_mesencult <- doGsea(df = as.data.frame(result_mesencult_gene_name),
                        gene_sets = msigdb_ont_split,
                        alpha = 0.05,
                        nPermSimple = 10000)
muscle_mesencult <- doGsea(df = as.data.frame(result_mesencult_gene_name),
                           gene_sets = msigdb_ont_muscle_split,
                           alpha = 0.05,
                           nPermSimple = 10000)

hallmarks_dmem <- doGsea(df = as.data.frame(result_dmem_gene_name),
                         gene_sets = msigdb_hallmarks_split,
                         alpha = 0.05,
                         nPermSimple = 10000)
onocgenic_dmem <- doGsea(df = as.data.frame(result_dmem_gene_name),
                         gene_sets = msigdb_oncogenic_split,
                         alpha = 0.05,
                         nPermSimple = 10000)
tf_dmem <- doGsea(df = as.data.frame(result_dmem_gene_name),
                  gene_sets = msigdb_tf_split,
                  alpha = 0.05,
                  nPermSimple = 10000)
reactome_dmem <- doGsea(df = as.data.frame(result_dmem_gene_name),
                        gene_sets = msigdb_reactome_split,
                        alpha = 0.05,
                        nPermSimple = 10000)
mir_dmem <- doGsea(df = as.data.frame(result_dmem_gene_name),
                        gene_sets = msigdb_mir_split,
                        alpha = 0.05,
                        nPermSimple = 10000)
ont_dmem <- doGsea(df = as.data.frame(result_dmem_gene_name),
                        gene_sets = msigdb_ont_split,
                        alpha = 0.05,
                        nPermSimple = 10000)
muscle_dmem <- doGsea(df = as.data.frame(result_dmem_gene_name),
                      gene_sets = msigdb_ont_muscle_split,
                      alpha = 0.05,
                      nPermSimple = 10000)

# Results and plotting, ORA
mesencult_up <- result_mesencult %>%
  as.data.frame() %>%
  dplyr::filter(log2FoldChange > 1,
         padj < 0.05)

mesencult_down <- result_mesencult %>%
  as.data.frame() %>%
  dplyr::filter(log2FoldChange < -1,
                padj < 0.05)

ora_go_mesencult_up <- enrichGO(gene = rownames(mesencult_up),
                             OrgDb = org.Mm.eg.db,
                             keyType = "ENSEMBL",
                             ont = "BP",
                             qvalueCutoff = 0.05,
                             pAdjustMethod = "BH",
                             universe = names(keep_genes))
ora_go_mesencult_up <- ora_go_mesencult_up@result

ora_go_mesencult_down <- enrichGO(gene = rownames(mesencult_down),
                                OrgDb = org.Mm.eg.db,
                                keyType = "ENSEMBL",
                                ont = "BP",
                                qvalueCutoff = 0.05,
                                pAdjustMethod = "BH",
                                universe = names(keep_genes))
ora_go_mesencult_down <- ora_go_mesencult_down@result

# Save some plots and tables
out_path <- "../output/"

write.csv(x = norm_counts, file = paste0(out_path, "tables/dicer1_rna_norm_counts_no_log.csv"), row.names = FALSE)
write.csv(x = result_mesencult_gene_name, file = paste0(out_path, "tables/dicer1_rna_mesencult_de_mut_vs_wt.csv"), row.names = TRUE)
write.csv(x = result_dmem_gene_name, file = paste0(out_path, "tables/dicer1_rna_dmem_de_mut_vs_wt.csv"), row.names = TRUE)
write.csv(x = result_dmem_gene_name, file = paste0(out_path, "tables/dicer1_rna_mesencult_up_ora_go_bp.csv"), row.names = FALSE)
write.csv(x = ora_go_mesencult_up, file = paste0(out_path, "tables/dicer1_rna_mesencult_up_ora_go_bp.csv"), row.names = FALSE)
write.csv(x = ora_go_mesencult_down, file = paste0(out_path, "tables/dicer1_rna_mesencult_dn_ora_go_bp.csv"), row.names = FALSE)

out_path_plots <- "../output/plots/"
ggsave(filename = "dicer1_rna_mesencult_gsea_hallmarks.png", path = out_path_plots, plot = hallmarks_mesencult$plot, bg = "white")
ggsave(filename = "dicer1_rna_mesencult_gsea_oncogenic.png", path = out_path_plots, plot = onocgenic_mesencult$plot, bg = "white")
ggsave(filename = "dicer1_rna_mesencult_gsea_tf.png", path = out_path_plots, plot = tf_mesencult$plot, bg = "white")
ggsave(filename = "dicer1_rna_mesencult_gsea_reactome.png", path = out_path_plots, plot = reactome_mesencult$plot, bg = "white")
ggsave(filename = "dicer1_rna_mesencult_gsea_mir.png", path = out_path_plots, plot = mir_mesencult$plot, bg = "white")
ggsave(filename = "dicer1_rna_mesencult_gsea_ont.png", path = out_path_plots, plot = ont_mesencult$plot, bg = "white")

ggsave(filename = "dicer1_rna_dmem_gsea_hallmarks.png", path = out_path_plots, plot = hallmarks_dmem$plot, bg = "white")
ggsave(filename = "dicer1_rna_dmem_gsea_oncogenic.png", path = out_path_plots, plot = onocgenic_dmem$plot, bg = "white")
ggsave(filename = "dicer1_rna_dmem_gsea_tf.png", path = out_path_plots, plot = tf_dmem$plot, bg = "white")
ggsave(filename = "dicer1_rna_dmem_gsea_reactome.png", path = out_path_plots, plot = reactome_dmem$plot, bg = "white")
ggsave(filename = "dicer1_rna_dmem_gsea_mir.png", path = out_path_plots, plot = mir_dmem$plot, bg = "white")
ggsave(filename = "dicer1_rna_dmem_gsea_ont.png", path = out_path_plots, plot = ont_dmem$plot, bg = "white")

ggsave(filename = "dicer1_rna_mut_vs_control.svg", path = out_path_plots, plot = dot_plot_all, width = 12, height = 8)
ggsave(filename = "dicer1_bulk_rna_pca.svg", path = out_path_plots, plot = pca_plot, width = 12, height = 8)
