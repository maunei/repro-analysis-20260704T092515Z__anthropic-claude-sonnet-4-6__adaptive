# PATCHED v2: Changed tumor_subcluster_partition_method from leiden (default) to random_trees
# to avoid Seurat v5 incompatibility inside infercnv::define_signif_tumor_subclusters()
# (ScaleData fails with "No layer matching pattern 'data' found" in Seurat v5 because
#  leiden subclustering omits NormalizeData before ScaleData).
# resume_mode=TRUE (default) picks up from checkpoint step 14 - only step 15 re-runs.
# num_threads=2: limits parallel workers during random_trees subclustering to avoid OOM
# (default was 4 threads x ~2.5 GB each = ~10 GB, exhausting 15 GB system RAM + all swap)

library(infercnv)

infercnv_obj_combined = CreateInfercnvObject(
  raw_counts_matrix = "../output/infer_cnv/dicer1_scrna_filter_combined_counts.matrix",
  annotations_file  = "../output/infer_cnv/dicer1_scrna_filter_combined_infer_cnv_cells.tsv",
  delim             = "\t",
  gene_order_file   = "../reference/gene_order.tsv",
  ref_group_names   = c("Control Immune",
                        "Control Fibro",
                        "Control Loop Henle",
                        "Control Collect Duct",
                        "Control Prox Tub",
                        "Control Endo",
                        "Control Mac Densa",
                        "Control Cycle",
                        "Control Trans Epi",
                        "Control Uni Fibro",
                        "Control Podo",
                        "Control Mural"))

infercnv_obj_result = infercnv::run(
  infercnv_obj_combined,
  cutoff                            = 0.1,
  out_dir                           = "../output/infer_cnv/combined_filter",
  cluster_by_groups                 = TRUE,
  denoise                           = TRUE,
  HMM                               = TRUE,
  tumor_subcluster_partition_method = "random_trees",  # PATCHED: avoids Seurat v5 leiden bug
  num_threads                       = 2,               # PATCHED: limit workers to avoid OOM (was 4)
  resume_mode                       = TRUE              # pick up from existing checkpoints
)

saveRDS(infercnv_obj_result, "../output/infer_cnv/combined_filter/run_final_infercnv_obj.rds")
cat("InferCNV run complete. Final object saved.\n")
