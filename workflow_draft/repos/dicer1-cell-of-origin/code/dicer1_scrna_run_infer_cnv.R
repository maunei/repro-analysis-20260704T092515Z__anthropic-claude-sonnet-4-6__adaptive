# Quick script to run inferCNV

# Required packages
library(infercnv)

# Read in the data; we set the reference to all the control cells
# If alternative comparisons are desired (ex. running all the unfiltered cells combined), specify the appropriate matrix and annotations here
infercnv_obj_combined = CreateInfercnvObject(raw_counts_matrix = "../output/infer_cnv/dicer1_scrna_filter_combined_counts.matrix",
                                             annotations_file = "../output/infer_cnv/dicer1_scrna_filter_combined_infer_cnv_cells.tsv",
                                             delim = "\t",
                                             gene_order_file = "../reference/gene_order.tsv",
                                             ref_group_names = c("Control Immune",
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

# Run inferCNV with standard parameters for 10X data
infercnv_obj_result = infercnv::run(infercnv_obj_combined,
                                    cutoff = 0.1,
                                    out_dir = "../output/infer_cnv/combined_filter",
                                    cluster_by_groups = T,
                                    denoise = T,
                                    HMM = T)
