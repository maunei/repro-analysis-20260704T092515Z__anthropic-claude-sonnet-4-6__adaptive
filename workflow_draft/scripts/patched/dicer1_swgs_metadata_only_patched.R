# dicer1_swgs_metadata_only_patched.R
# Independent reimplementation — metadata/clinical plots only
# BLOCKED: dicer1_swgs_50kb_noXY_copyNumbersSegmented_filtered_called.rds is not publicly
# available on GEO. Raw data (PRJNA1221971) requires QDNAseq preprocessing which cannot
# run on this system. Only the clinical metadata plots from 'master' are produced here.

library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)

master <- read.csv("../metadata/dicer1_swgs_sample_metadata.csv") %>%
  mutate(across(everything(), ~ ifelse(. == "", NA, as.character(.)))) %>%
  mutate(across(everything(), ~ ifelse(. == "n.a.", NA, as.character(.)))) %>%
  dplyr::rename("sample_id" = "Sample_ID",
                "core_id"   = "Core_ID") %>%
  mutate(Genotype   = str_trim(Genotype),
         Kras = case_when(is.na(Kras) ~ "Unknown", Kras == "wt" ~ "Wt", .default = "Missense"),
         Trp53 = case_when(is.na(Trp53) ~ "Unknown", Trp53 == "wt" ~ "Wt", .default = "Missense"),
         Hras = case_when(is.na(Hras) ~ "Unknown", Hras == "wt" ~ "Wt", .default = "Missense"),
         Nras = case_when(is.na(Nras) ~ "Unknown", Nras == "wt" ~ "Wt", .default = "Missense"),
         Dicer1 = case_when(is.na(Dicer1) ~ "Unknown", Dicer1 == "wt" ~ "Wt", .default = "Missense"),
         Kras.FISH = case_when(is.na(`Kras.FISH`) ~ "Unknown",
                               `Kras.FISH` == "amplifcation" ~ "Amplification",
                               .default = "No amplification"),
         Histotype = factor(Histotype, levels = c("LGMT-like", "SARC-like"))) %>%
  filter(Trp53 != "Unknown") %>%
  mutate(sample_id = factor(sample_id, levels = sample_id))

cat("Master metadata loaded:", nrow(master), "samples\n")
print(head(master))

out_dir <- "../output/plots/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Clinical metadata plots (independent of CN data)
gender_plot <- ggplot(data = master, mapping = aes(y = sample_id, x = Gender)) +
  geom_tile(mapping = aes(fill = Gender)) +
  labs(title = "Gender") +
  scale_x_discrete(guide = guide_axis(angle = 0)) +
  theme_classic(base_size = 9) +
  theme(axis.title = element_blank(), legend.position = "none", axis.text.x = element_blank())

histotype_plot <- ggplot(data = master, mapping = aes(y = sample_id, x = Histotype)) +
  geom_tile(mapping = aes(fill = Histotype)) +
  labs(title = "Histotype") +
  scale_x_discrete(guide = guide_axis(angle = 0)) +
  theme_classic(base_size = 9) +
  theme(axis.title = element_blank(), legend.position = "none", axis.text.x = element_blank())

tp53_plot <- ggplot(data = master, mapping = aes(y = sample_id, x = Trp53)) +
  geom_tile(mapping = aes(fill = Trp53)) +
  labs(title = "Trp53") +
  scale_x_discrete(guide = guide_axis(angle = 0)) +
  theme_classic(base_size = 9) +
  theme(axis.title = element_blank(), legend.position = "none", axis.text.x = element_blank())

kras_plot <- ggplot(data = master, mapping = aes(y = sample_id, x = Kras)) +
  geom_tile(mapping = aes(fill = Kras)) +
  labs(title = "Kras") +
  scale_x_discrete(guide = guide_axis(angle = 0)) +
  theme_classic(base_size = 9) +
  theme(axis.title = element_blank(), legend.position = "none", axis.text.x = element_blank())

dicer1_plot <- ggplot(data = master, mapping = aes(y = sample_id, x = Dicer1)) +
  geom_tile(mapping = aes(fill = Dicer1)) +
  labs(title = "Dicer1") +
  scale_x_discrete(guide = guide_axis(angle = 0)) +
  theme_classic(base_size = 9) +
  theme(axis.title = element_blank(), legend.position = "none", axis.text.x = element_blank())

# Combine metadata panels
metadata_panel <- gender_plot | histotype_plot | tp53_plot | kras_plot | dicer1_plot
ggsave("dicer1_swgs_metadata_panel.svg", plot=metadata_panel,
       path=out_dir, width=10, height=7, device="svg")
ggsave("dicer1_swgs_metadata_panel.png", plot=metadata_panel,
       path=out_dir, width=10, height=7, dpi=200)

write.csv(master, "../output/tables/dicer1_swgs_sample_metadata.csv", row.names=FALSE)

cat("Metadata plots saved.\n")
cat("NOTE: CN plots (chr1/chr6) BLOCKED — dicer1_swgs_50kb_noXY_copyNumbersSegmented_filtered_called.rds\n")
cat("      is not publicly available on GEO. Raw data (PRJNA1221971) requires QDNAseq preprocessing.\n")
