# RUN ORDER — DICER1 sarcoma (paper 653)

Run folder: `20260704T030158Z__anthropic-claude-sonnet-4-6__adaptive`  
Resumed from: `20260703T220654Z__anthropic-claude-sonnet-4-6__adaptive`

---

## Step 1 — scRNA-seq QC
**Script:** `code/dicer1_scrna_qc.R`  
**Runtime:** 653_renv (R 4.2.2)  
**Status:** ✅ COMPLETED  
**Inputs:** 8 Cell Ranger output dirs (GSE288990) + sample metadata CSV  
**Outputs:** `output/scrna/dicer1_scrna_sce_qc.rds`  

---

## Step 2a — Control sample integration
**Script:** `workflow_draft/scripts/patched/dicer1_scrna_controls_patched.R`  
**Runtime:** 653_renv (R 4.2.2)  
**Status:** ✅ COMPLETED (k.weight patched from 80→50 due to anchor count)  
**Inputs:** `dicer1_scrna_sce_qc.rds`  
**Outputs:** `dicer1_scrna_unfilter_control.rds`, `dicer1_scrna_filter_control.rds`  

---

## Step 2b — Tumor sample integration
**Scripts:** `dicer1_scrna_tumors.R` → `dicer1_scrna_tumors_filter_only_patched.R`  
**Runtime:** 653_renv (R 4.2.2)  
**Status:** ✅ COMPLETED  
**Inputs:** `dicer1_scrna_sce_qc.rds`  
**Outputs:** `dicer1_scrna_unfilter_tumor.rds`, `dicer1_scrna_filter_tumor.rds`  

---

## Step 3a — Bulk RNA-seq DE + GSEA
**Script:** `workflow_draft/scripts/patched/dicer1_rna_patched.R`  
**Runtime:** 653_renv (R 4.2.2)  
**Status:** ✅ COMPLETED (msigdbr 7.5.1 category codes corrected: MH→H, M3→C3, M5→C5)  
**Inputs:** `data/star_salmon/*/quant.sf`; `tx2gene.tsv`  
**Outputs:** `output/tables/dicer1_rna_*.csv`; `output/plots/dicer1_rna_*.png`  

---

## Step 3b — Joint control+tumor integration
**Script:** `workflow_draft/scripts/patched/dicer1_scrna_controls_and_tumors_patched.R`  
**Runtime:** 653_renv (R 4.2.2)  
**Status:** ✅ COMPLETED  
**Inputs:** `dicer1_scrna_filter_control.rds` + `dicer1_scrna_filter_tumor.rds`  
**Outputs:** `dicer1_scrna_filter_control_and_tumor.rds`  

---

## Step 3c — Pseudotime trajectory (Monocle3)
**Script:** `workflow_draft/scripts/patched/dicer1_scrna_trajectory_patched.R`  
**Runtime:** 653_renv2 (R 4.5.3)  
**Status:** ✅ COMPLETED (branch endpoints automated; root node approximated from Progen centroid)  
**Inputs:** `dicer1_scrna_filter_tumor.rds`  
**Outputs:** `dicer1_scrna_trajectory_cds.rds`; `output/plots/dicer1_scrna_filter_tumor_trajectory_*.png`  

---

## Step 3d — PROGENy pathway activity (decoupleR MLM)
**Script:** `workflow_draft/scripts/patched/dicer1_scrna_run_pathway_analysis_patched.R`  
**Runtime:** 653_renv (R 4.2.2)  
**Status:** ✅ COMPLETED (OmnipathR bypassed; local TSV from OmniPath REST API; S7 FeaturePlot fallback)  
**Inputs:** `dicer1_scrna_filter_tumor.rds`; `data/progeny_omnipath.tsv`  
**Outputs:** `dicer1_scrna_filter_tumor_pathway_activity.csv`; `output/plots/dicer1_scrna_filter_tumor_pathway_*.{png,svg}`  

---

## Step 4a — inferCNV matrix export
**Script:** `workflow_draft/scripts/patched/dicer1_scrna_export_infer_cnv_filter_only.R`  
**Runtime:** 653_renv (R 4.2.2)  
**Status:** ✅ COMPLETED — produced `filter_combined_counts.matrix` (unfiltered control + filtered tumor cells; 624MB), which is the only matrix consumed by `dicer1_scrna_run_infer_cnv.R`. Author's script also exports three unfiltered-only matrices (`unfilter_control`, `unfilter_tumor`, `unfilter_combined`) — these were not written (OOM) and are not needed by any downstream script.  
**Inputs:** `dicer1_scrna_unfilter_control.rds` + `dicer1_scrna_filter_tumor.rds`  
**Outputs:** `output/infer_cnv/dicer1_scrna_filter_combined_counts.matrix` (624MB); `dicer1_scrna_filter_combined_infer_cnv_cells.tsv`  

---

## Step 4b — inferCNV run
**Script:** `code/dicer1_scrna_run_infer_cnv_v2_patched.R` (patched: random_trees subclustering)  
**Runtime:** 653_renv (R 4.2.2)  
**Status:** ✅ COMPLETED (completed July 2026, post-initial-run, with 2-thread constraint to reduce memory pressure)  
**Inputs:** `dicer1_scrna_filter_combined_counts.matrix`; `dicer1_scrna_filter_combined_infer_cnv_cells.tsv`; `reference/gene_order.tsv`  
**Outputs:** `output/infer_cnv/combined_filter/run.final.infercnv_obj`; `run_final_infercnv_obj.rds`; HMM CNV heatmaps (`infercnv.png`, `infercnv.17_HMM_pred...png`, `infercnv.19_HMM_pred...png`, `infercnv.20_HMM_pred...png`)  
**Note:** Seurat v5 leiden incompatibility at step 15 (`define_signif_tumor_subclusters`) required `tumor_subcluster_partition_method = "random_trees"`; cache invalidated forcing full re-run from step 1.

---

## Step 4c — Add inferCNV to Seurat
**Script:** `workflow_draft/scripts/patched/dicer1_scrna_add_to_seurat_intermediate.R` → `code/dicer1_scrna_add_infer_cnv_seurat.R`  
**Runtime:** 653_renv (R 4.2.2)  
**Status:** ✅ COMPLETED (completed July 2026, post-initial-run)  
**Outputs:** `output/infer_cnv/dicer1_scrna_filter_combined_infercnv_metadata.rds`; `dicer1_scrna_filter_tumor_infercnv.{png,svg}`; `dicer1_scrna_filter_tumor_infercnv_split_sample.{png,svg}`  
**Note:** `infercnv::add_to_seurat()` step was absent from author's `dicer1_scrna_run_infer_cnv.R`; intermediate script created to produce the required `infercnv_metadata.rds` before the author's add_infer_cnv_seurat script ran.

---

## Step X1 — Xenium mouse spatial (mus musculus)
**Script:** `workflow_draft/scripts/patched/dicer1_xenium_proseg_mus_musculus_patched.R`  
**Runtime:** 653_renv2 (R 4.5.3)  
**Status:** ❌ BLOCKED — OOM; requires ≥32 GB RAM  
**Inputs:** `output/xenium/dicer1_xenium_mus_musculus_all.rds` (2.2GB, GSE289001) — downloaded and present  
**Blocker:** After QC filtering the WT subset has 622K cells. SCTransform on 622K cells requires ~30–50 GB RAM for the dense residual matrix; this 15 GB system cannot run it. Three attempts failed (v1: FOV class error; v2: OOM on 1.5M full object; v3: OOM on 622K WT subset).  
**Ready to run:** patched script (v3) is ready; no code or data gap. Run on ≥32 GB RAM system.  

---

## Step X2 — Xenium human spatial (homo sapiens)
**Status:** ❌ BLOCKED — Human Xenium RDS not deposited in any public repository  

---

## Step X3 — Shallow WGS copy number
**Status:** ❌ BLOCKED — Processed sWGS RDS not deposited; raw WGS (PRJNA1221971) excluded per directive  

