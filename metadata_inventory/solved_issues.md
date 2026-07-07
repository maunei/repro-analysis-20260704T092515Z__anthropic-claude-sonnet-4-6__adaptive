# Solved Issues Log
# Paper 653 - DICER1 Sarcoma scRNA/Spatial Reproducibility Analysis

## Issues from prior run (20260703T220654Z) + this run (20260704T030158Z)

---

## rowData function not found in controls script
Error: `could not find function "rowData"`
Tried: n/a
Fix: Added `library(SingleCellExperiment)` to patched controls script
Affected: dicer1_scrna_controls.R
Patched files: workflow_draft/scripts/patched/dicer1_scrna_controls_patched.R

---

## Duplicate gene names in Seurat object creation
Error: `make.names(rownames(assay), unique = FALSE)` produces non-unique names; CreateSeuratObject fails
Tried: n/a
Fix: Applied `make.unique()` on rownames before CreateSeuratObject in controls and tumors scripts
Affected: dicer1_scrna_controls.R, dicer1_scrna_tumors.R
Patched files: workflow_draft/scripts/patched/dicer1_scrna_controls_patched.R, dicer1_scrna_tumors_patched.R

---

## Seurat v5 DotPlot deprecated ggplot2 API error
Error: `Error in as.vector(x, mode): cannot coerce type 'S7_object' to vector of type 'any'`
Tried: n/a
Fix: Wrapped DotPlot calls in `tryCatch` with safe_dotplot fallback function
Affected: dicer1_scrna_controls.R, dicer1_scrna_tumors.R, dicer1_scrna_controls_and_tumors.R
Patched files: workflow_draft/scripts/patched/

---

## Seurat v5 requires JoinLayers before FindAllMarkers
Error: `FindAllMarkers() requires the assay to have a single layer`
Tried: n/a
Fix: Added `JoinLayers()` call after integration and before FindAllMarkers
Affected: dicer1_scrna_controls.R, dicer1_scrna_tumors.R
Patched files: workflow_draft/scripts/patched/

---

## Seurat v5 DoHeatmap off-by-one with downsample parameter
Error: Incorrect downsample count in DoHeatmap with Seurat v5
Tried: n/a
Fix: Wrapped DoHeatmap calls in tryCatch with safe_doheatmap helper
Affected: dicer1_scrna_controls.R, dicer1_scrna_tumors.R
Patched files: workflow_draft/scripts/patched/

---

## dicer1_scrna_tumors.R OOM on second FindAllMarkers (full markers run)
Error: R process killed with exit code 137 (SIGKILL, OOM) during FindAllMarkers(only.pos=FALSE)
Tried: Full tumors script, holding both unfilter_tumor RDS (528MB) and performing FindAllMarkers on all clusters simultaneously
Fix: Split script: run full tumors to produce unfilter_tumor.rds, then use separate filter-only split script (dicer1_scrna_tumors_filter_only_patched.R) that loads unfilter_tumor.rds, frees it, and runs filtered analysis with aggressive gc()
Affected: dicer1_scrna_tumors.R
Patched files: workflow_draft/scripts/patched/dicer1_scrna_tumors_patched.R, dicer1_scrna_tumors_filter_only_patched.R

---

## msigdbr API incompatibility in dicer1_rna.R (msigdbr 7.5.1)
Error: `argument "db_species" not found` and category code mismatch (MH vs H, M3 vs C3, etc.)
Tried: n/a — error was clear from msigdbr changelog
Fix: Patched all msigdbr() calls to remove db_species parameter and update category codes (MH→H, M3→C3, M5→C5, subcollection→subcategory)
Affected: dicer1_rna.R
Patched files: workflow_draft/scripts/patched/dicer1_rna_patched.R

---

## Memory management in dicer1_rna.R - OOM at plot stage
Error: OOM kill (exit code 137) at plot-save stage after GSEA completed
Tried: Script ran successfully to plot stage but was killed
Fix: Added rm()/gc() after each large msigdbr gene set object and after each GSEA result. Moved GO:BP gene set cleanup before plot rendering.
Affected: dicer1_rna.R
Patched files: workflow_draft/scripts/patched/dicer1_rna_patched.R

---

## PX data directories named as HDT IDs in GSE288990 download
Error: GSE288990_RAW.tar contained per-sample files named by GSM IDs (GSM8780505_HDT789_799_...) not PX IDs
Tried: n/a
Fix: Extracted GSM* files and renamed/moved them to PX ID directories based on sample_metadata.csv mapping (HDT789_799→PX3085_GTGGATCAAA-CAGGGTTGGC, etc.)
Affected: Data layout for dicer1_scrna_qc.R
Patched files: None (data reorganization only)

---

## dicer1_scrna_run_pathway_analysis.R — OmnipathR package missing
Error: `Error in loadNamespace(x) : there is no package called 'OmnipathR'` in decoupleR::get_progeny()
Tried: n/a
Fix: Downloaded PROGENy annotation data directly from OmniPath REST API (https://omnipathdb.org/annotations?resources=PROGENy&format=tsv) and saved to data/progeny_omnipath.tsv. Implemented local_get_progeny() function that reads from local TSV, bypassing OmnipathR entirely.
Affected: dicer1_scrna_run_pathway_analysis.R
Patched files: workflow_draft/scripts/patched/dicer1_scrna_run_pathway_analysis_patched.R

---

## Seurat v5 Assay5 ScaleData fails — no data layer
Error: `No layer matching pattern 'data' found. Please run NormalizeData and retry` from ScaleData on pathwaysmlm assay
Tried: CreateAssay5Object + ScaleData (fails)
Fix: Added SetAssayData() call to explicitly set the data layer = pathway_mat before ScaleData. For the @data slot assignment post-ScaleData, used SetAssayData() via tryCatch (Seurat v5 path).
Affected: dicer1_scrna_run_pathway_analysis_patched.R
Patched files: workflow_draft/scripts/patched/dicer1_scrna_run_pathway_analysis_patched.R

---

## Seurat v5 FeaturePlot S7 object incompatibility with & operator
Error: `Incompatible methods ("Ops.S7_object", "&.gg") for "&"` during internal FeaturePlot color scale application
Tried: FeaturePlot() + scale_colour_gradient2() — fails internally
Fix: Implemented safe_featureplot_pathway() wrapper using FetchData() + ggplot2 directly; wraps each FeaturePlot call in tryCatch with fallback to manual UMAP scatter plot with score overlay
Affected: dicer1_scrna_run_pathway_analysis_patched.R
Patched files: workflow_draft/scripts/patched/dicer1_scrna_run_pathway_analysis_patched.R

---

## InferCNV export OOM on writing unfiltered combined matrices
Error: R writes control matrix (471MB) successfully but OOM when writing unfiltered tumor matrix (15,346 cells × 32,286 genes dense)
Tried: Full export script — fails after first two write.table() calls
Fix: Created filter-only export script (dicer1_scrna_export_infer_cnv_filter_only.R) that writes ONLY the filter_combined matrix (2,483 tumor + 12,474 control = 14,957 cells) in chunks, skipping the unfiltered combined matrices that are not needed for inferCNV
Affected: dicer1_scrna_export_infer_cnv_patched.R
Patched files: workflow_draft/scripts/patched/dicer1_scrna_export_infer_cnv_filter_only.R

---

## biomaRt getBM fails — BiocFileCache dplyr API incompatibility
Error: `Failed to collect lazy table. Arguments in ... must be used. ..1 = Inf`
Tried: biomaRt::useEnsembl with version=102 and getBM
Fix: Switched to direct Ensembl REST API via Python batch POST requests. All 32,285 ENSMUSG IDs queried in batches of 1000 via https://rest.ensembl.org/lookup/id
Affected: gene_order.tsv generation
Patched files: /tmp/fetch_gene_order.py → reference/gene_order.tsv

---

## Ensembl REST API timeout for batch gene coordinates
Error: `The read operation timed out` on all batches from Python script
Tried: Python batch POST to rest.ensembl.org/lookup/id — timed out
Fix: Trying biomaRt with https mirrors as fallback; alternatively will use UCSC API. Gene_order.tsv generation is in progress.
Affected: reference/gene_order.tsv
Status: PENDING

---

## wget not available in container
Error: `bash: line 1: wget: command not found`
Tried: wget download command
Fix: Switched to curl for all downloads
Affected: Xenium RDS download
Patched files: None (command replacement)


---

## InferCNV leiden tumor subclustering fails with Seurat v5 (dicer1_scrna_run_infer_cnv.R)
Error: `Error in ScaleData(): No layer matching pattern 'data' found. Please run NormalizeData and retry` at step 15 during leiden subclustering.
Tried: Prior run failed at exactly this point after completing steps 1-14.
Fix: Changed `tumor_subcluster_partition_method = "random_trees"` in patched v2 script. Note: this invalidated prior run checkpoints (parameter change forces full re-run from step 1; resume_mode=TRUE only works when parameters match). Full re-run from step 1 took ~40 minutes.
Affected: dicer1_scrna_run_infer_cnv.R
Patched files: workflow_draft/repos/dicer1-cell-of-origin/code/dicer1_scrna_run_infer_cnv_v2_patched.R

---

## Xenium WT filter keeps 622K cells (not 14K as expected)
Error: No error — but unexpected cell count. v3 patched script combined QC + Inclusion_control filters into one step applied to 1.5M cell full object. The Inclusion_control column does NOT restrict to WT-only; it applies to most cells after QC (622K pass). SCTransform on 622K cells requires ≥32 GB RAM and is OOM-killed on this 15 GB system.
Tried: v1 (FOV class error), v2 (OOM on 1.5M cells), v3 (OOM on 622K WT cells)
Fix: None viable on this hardware. Xenium analysis classified as blocked by hardware constraint.
Affected: dicer1_xenium_proseg_mus_musculus_patched.R
Patched files: workflow_draft/scripts/patched/dicer1_xenium_proseg_mus_musculus_patched.R (ready for ≥32 GB RAM system)

---

## dicer1_scrna_add_infer_cnv_seurat.R reads file not produced by author inferCNV script
Error: `dicer1_scrna_filter_combined_infercnv_metadata.rds` expected but never written by dicer1_scrna_run_infer_cnv.R.
Tried: Inspected author script — infercnv::add_to_seurat() step is missing.
Fix: Created workflow_draft/scripts/patched/dicer1_scrna_add_to_seurat_intermediate.R that calls add_to_seurat() after inferCNV completes and saves the metadata-enriched Seurat object.
Affected: dicer1_scrna_add_infer_cnv_seurat.R
Patched files: workflow_draft/scripts/patched/dicer1_scrna_add_to_seurat_intermediate.R
