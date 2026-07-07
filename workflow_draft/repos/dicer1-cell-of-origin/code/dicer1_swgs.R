# Required packages
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(QDNAseq)
library(utanos)
library(exact2x2)

# Restructure metadata and filter for samples of interest; metadata as in table S2
master <- read.csv("../metadata/dicer1_swgs_sample_metadata.csv") %>%
  mutate(across(everything(), ~ ifelse(. == "", NA, as.character(.)))) %>%
  mutate(across(everything(), ~ ifelse(. == "n.a.", NA, as.character(.)))) %>%
  dplyr::rename("sample_id" = "Sample_ID",
                "core_id" = "Core_ID") %>%
  mutate(Genotype = str_trim(Genotype),
         Kras = case_when(is.na(Kras) ~ "Unknown",
                          Kras == "wt" ~ "Wt",
                          .default = "Missense"),
         Trp53 = case_when(is.na(Trp53) ~ "Unknown",
                           Trp53 == "wt" ~ "Wt",
                          .default = "Missense"),
         Hras = case_when(is.na(Hras) ~ "Unknown",
                          Hras == "wt" ~ "Wt",
                           .default = "Missense"),
         Nras = case_when(is.na(Nras) ~ "Unknown",
                          Nras == "wt" ~ "Wt",
                           .default = "Missense"),
         Dicer1 = case_when(is.na(Dicer1) ~ "Unknown",
                            Dicer1 == "wt" ~ "Wt",
                           .default = "Missense"),
         Kras.FISH = case_when(is.na(Kras.FISH) ~ "Unknown",
                               Kras.FISH == "amplifcation" ~ "Amplification",
                               .default = "No amplification"),
         Histotype = factor(Histotype,
                               levels = c("LGMT-like",
                                          "SARC-like"))) %>%
  filter(Trp53 != "Unknown") %>%
  mutate(sample_id = factor(sample_id, levels = sample_id)) %>%
  arrange(Histotype)

# Read in the CN data and extract the required fields, reformat
rcn_obj <- readRDS(file = "../data/dicer1_swgs_50kb_noXY_copyNumbersSegmented_filtered_called.rds")

rcn_segs_df <- ExportBinsQDNAObj(object = rcn_obj, type = "segments") %>%
  pivot_longer(cols = !any_of(c("feature", "chromosome", "start", "end")), names_to = "sample_id", values_to = "segmented") %>%
  arrange(sample_id) %>%
  dplyr::select(!feature) %>%
  relocate(sample_id) %>%
  mutate(segmented = log2(segmented))

rcn_calls_df <- ExportBinsQDNAObj(object = rcn_obj, type = "calls") %>%
  pivot_longer(cols = !any_of(c("feature", "chromosome", "start", "end")), names_to = "sample_id", values_to = "call") %>%
  arrange(sample_id) %>%
  dplyr::select(!feature) %>%
  relocate(sample_id)

segs_calls_df_collapse <- CopyNumberSegments(rcn_segs_df) %>%
  left_join(y = rcn_calls_df, by = join_by(sample_id, chromosome, start <= start, end >= end)) %>%
  select(!c(start.y, end.y)) %>%
  rename("start" = "start.x",
         "end" = "end.x") %>%
  distinct() %>%
  mutate(cn_call = case_when(call == -1 ~ "loss", call == 1 ~ "gain", call == -2 ~ "deletion", call == 2 ~ "amplification", .default = "neutral"))

cn_summary <- segs_calls_df_collapse %>% # Combining loss/del and gain/amp for summary -- we keep them separate when assessing focal changes (ex. Kras)
  left_join(y = master, by = join_by(sample_id == core_id)) %>% # Map the original sample core IDs to the paper IDs
  mutate(cn_call = case_when(cn_call ==  "deletion" ~ "loss", cn_call == "amplification" ~ "gain", .default = cn_call)) %>%
  group_by(sample_id, chromosome) %>%
  mutate(total_bins = sum(bin_count)) %>%
  ungroup() %>%
  pivot_wider(names_from = cn_call, values_from = bin_count,
              values_fill = 0) %>%
  group_by(sample_id, chromosome, total_bins) %>%
  summarise(gain = sum(gain),
            neutral = sum(neutral),
            loss = sum(loss)) %>%
  ungroup() %>%
         mutate(status = case_when(gain/total_bins > 0.9 ~ "Gain", # Require 90% of the chromosome to be altered for changes to be considered whole-chromosome
                            loss/total_bins > 0.9 ~ "Loss",
                            .default = "Neutral")) %>%
  right_join(y = master, by = join_by(sample_id == core_id)) %>% # Here we add back HDT395_kidney_1_tumor for plotting, which failed QC after sequencing
  rename("paper_id" = "sample_id.y")

chr1_data <- cn_summary %>%
  filter(chromosome == "1" | is.na(chromosome)) %>%
  mutate(paper_id = factor(x = paper_id, levels = levels(master$sample_id)), # Order the samples in the same way as the master table
         status = if_else(is.na(status), "Unknown", status)) # Handle the QC-failed missing case here

chr6_data <- cn_summary %>%
  filter(chromosome == "6" | is.na(chromosome)) %>%
  mutate(paper_id = factor(x = paper_id, levels = levels(master$sample_id)), # Order the samples in the same way as the master table
         status = if_else(is.na(status), "Unknown", status)) # Handle the missing QC-failed case here

# Panel-style plot of genotype, histology, Tp53 mutations, Kras mutations, CN gains
chr1_plot <- ggplot(data = chr1_data,
                    mapping = aes(x = paper_id, y = 1, fill = as.factor(status))) +
  geom_tile() +
  scale_fill_manual(values = c("Gain" = "#be3324",
                               "Neutral" = "#dedede",
                               "Unknown" = "darkgray"),
                    name = "Chr1 CN alterations") +
  theme_void() +
  theme(plot.title = element_text(size = 7)) +
  ggtitle("Chr1 CN alterations") +
  coord_flip()

chr6_plot <- ggplot(data = chr6_data,
                    mapping = aes(x = paper_id, y = 1, fill = as.factor(status))) +
  geom_tile() +
  scale_fill_manual(values = c("Gain" = "#be3324",
                               "Neutral" = "#dedede",
                               "Unknown" = "darkgray"),
                    name = "Chr6 CN alterations") +
  theme_void() +
  theme(plot.title = element_text(size = 7)) +
  ggtitle("Chr6 CN alterations") +
  coord_flip()

gender_plot <- ggplot(data = master,
                      mapping = aes(x = sample_id, y = 1, fill = as.factor(Gender))) +
  geom_tile() +
  scale_fill_manual(values = c("M" = "#98632e",
                               "F" = "#3f8370"),
                    name = "Gender") +
  theme_void() +
  theme(plot.title = element_text(size = 7),
        axis.text.y = element_text(size = 7)) +
  ggtitle("Gender") +
  coord_flip()

histotype_plot <- ggplot(data = master,
                         mapping = aes(x = sample_id, y = 1, fill = as.factor(Histotype))) +
  geom_tile() +
  scale_fill_manual(values = c("SARC-like" = "#e47249",
                               "LGMT-like" = "#f0b76d"),
                    name = "Histotype") +
  theme_void() +
  theme(plot.title = element_text(size = 7)) +
  ggtitle("Histotype") +
  coord_flip()

tp53_plot <- ggplot(data = master,
                    mapping = aes(x = sample_id, y = 1, fill = as.factor(Trp53))) +
  geom_tile() +
  scale_fill_manual(values = c("Complex" = "#e47249",
                               "Missense" = "#f0b76d",
                               "Wt" = "#dedede",
                               "Unknown" = "darkgray"),
                    name = "Trp53 alterations") +
  theme_void() +
  theme(plot.title = element_text(size = 7)) +
  ggtitle("Trp53 alterations") +
  coord_flip()

kras_plot <- ggplot(data = master,
                    mapping = aes(x = sample_id, y = 1, fill = as.factor(Kras))) +
  geom_tile() +
  scale_fill_manual(values = c("Missense" = "#f0b76d",
                               "Wt" = "#dedede",
                               "Unknown" = "darkgray"),
                    name = "Kras alterations") +
  theme_void() +
  theme(plot.title = element_text(size = 7)) +
  ggtitle("Kras alterations") +
  coord_flip()

hras_plot <- ggplot(data = master,
                    mapping = aes(x = sample_id, y = 1, fill = as.factor(Hras))) +
  geom_tile() +
  scale_fill_manual(values = c("Missense" = "#f0b76d",
                               "Wt" = "#dedede",
                               "Unknown" = "darkgray"),
                    name = "Hras alterations") +
  theme_void() +
  theme(plot.title = element_text(size = 7)) +
  ggtitle("Hras alterations") +
  coord_flip()

nras_plot <- ggplot(data = master,
                    mapping = aes(x = sample_id, y = 1, fill = as.factor(Nras))) +
  geom_tile() +
  scale_fill_manual(values = c("Missense" = "#f0b76d",
                               "Wt" = "#dedede",
                               "Unknown" = "darkgray"),
                    name = "Nras alterations") +
  theme_void() +
  theme(plot.title = element_text(size = 7)) +
  ggtitle("Nras alterations") +
  coord_flip()

dicer1_plot <- ggplot(data = master,
                    mapping = aes(x = sample_id, y = 1, fill = as.factor(Dicer1))) +
  geom_tile() +
  scale_fill_manual(values = c("Missense" = "#f0b76d",
                               "Wt" = "#dedede",
                               "Unknown" = "darkgray"),
                    name = "Dicer1 alterations") +
  theme_void() +
  theme(plot.title = element_text(size = 7)) +
  ggtitle("Dicer1 alterations") +
  coord_flip()

kras_amp_plot <- ggplot(data = master,
                        mapping = aes(x = sample_id, y = 1, fill = as.factor(Kras.FISH))) +
  geom_tile() +
  scale_fill_manual(values = c("Amplification" = "#be3324",
                               "No amplifcation" = "#dedede",
                               "Unknown" = "darkgray"),
                    name = "Kras CN alterations") +
  theme_void() +
  theme(plot.title = element_text(size = 7)) +
  ggtitle("Kras CN alterations") +
  coord_flip()

panel_plot <- gender_plot +
  histotype_plot +
  dicer1_plot +
  tp53_plot +
  kras_plot +
  kras_amp_plot +
  chr1_plot +
  chr6_plot +
  plot_layout(nrow = 1,
              ncol = 8,
              guides = "collect",
              heights = c(100, 100, 100, 100, 100, 100, 100, 100))

# CN frequency plots by histotype
sarc_samples <- master %>%
  filter(Histotype == "SARC-like",
        core_id %in% rcn_obj@phenoData@data$name)

lgmt_samples <- master %>%
  filter(Histotype == "LGMT-like",
         core_id %in% rcn_obj@phenoData@data$name)

rcn_sarc <- rcn_obj[, sarc_samples$core_id]
rcn_lgmt <- rcn_obj[, lgmt_samples$core_id]

SummaryCNPlot(x = rcn_sarc, main = glue::glue('Relative copy-number frequency (SARC-like)'),
                      summarytype = 'frequency',
                      maskprob = 0, maskaberr = 0,
                      gaincol='red', losscol='blue', misscol=NA,
                      build='mm10', plotXY = FALSE)

SummaryCNPlot(x = rcn_lgmt, main = glue::glue('Relative copy-number frequency (LGMT-like)'),
              summarytype = 'frequency',
              maskprob = 0, maskaberr = 0,
              gaincol='red', losscol='blue', misscol=NA,
              build='mm10', plotXY = FALSE)

# Statistics on chromosome 1 gains
chr1_data_filter <- chr1_data %>%
  filter(status != "Unknown")

chr1_table <- table(chr1_data_filter$Histotype,
                    chr1_data_filter$status)
chr1_table_reorder <- rbind(chr1_table[2, ],
                            chr1_table[1, ])

chr1_test <- fisher.exact(x = chr1_table_reorder,
                          alternative = "two.sided",
                          conf.level = 0.95)

# Export plots
ggsave(filename = "dicer1_figure_1_panel_plot.svg", plot = panel_plot, path = "../output/plots/", width = 10, height = 11.5, device = "svg")
ggsave(filename = "dicer1_figure_1_panel_plot.png", plot = panel_plot, path = "../output/plots/", width = 10, height = 11.5, device = "png")
