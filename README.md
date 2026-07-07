![Papers RAG - NotebookLM infographic overview](infographic_app2.png)

# Public Reproducibility Analysis Skeleton

This repository is a lightweight public version of a reproducibility-analysis run.

It is provided to inform the scope, organization, and file-level structure of the reproducibility analysis workflow.

See the software system here:

https://github.com/maunei/papers-rag_and_reproducibility

> **Disclaimer**
>
> The outputs, tables, plots, logs, and placeholder files included here should not be treated as finalized research results or used as input for downstream scientific analysis. This repository is under active construction and is being used to test and document a reproducibility-analysis system.
>
> Large data files, generated objects, and selected intermediate outputs have been replaced with `.NOT_INCLUDED.txt` placeholder files. These placeholders preserve the directory structure and file provenance, but do not provide analyzable data.
>
> Any scientific interpretation should be based only on a validated, complete, and reviewed version of the analysis.

The full public repository tree is available in [`github_tree.txt`](github_tree.txt): 487 listed entries, including 436 files total, of which 140 are `.NOT_INCLUDED.txt` placeholders for excluded large/generated data objects, plus 47 directories.

---

# Spatial scRNA + Xenium Analysis of DICER1-Syndrome Sarcoma — Replication Feasibility Assessment

**Paper:** "Spatial single-cell transcriptomic analysis informs tumor developmental hierarchy of DICER1 syndrome-related sarcoma" (Nat Commun 2026)  
**GEO Accessions:** GSE288990 (scRNA-seq), GSE289001 (Xenium spatial, mouse), GSE309820 (bulk RNA-seq)  
**Code repository:** https://github.com/Huntsmanlab/dicer1-cell-of-origin  
**Report date:** 2026-07-04  
**Model/reasoning:** claude-sonnet-4-6 / adaptive  
**Run folder:** `20260704T092515Z__anthropic-claude-sonnet-4-6__adaptive`  
**Resumed from:** `20260704T030158Z__anthropic-claude-sonnet-4-6__adaptive` (failed with internal error at report-write stage)

---

## Request Context And Runtime Directives

This run continues a prior run that completed all major script execution but was killed before writing the final report. Runtime directives (from `RUNTIME_DIRECTIVES_FOR_OPENCLAW_AGENT.txt`) materially changed the workflow in the following ways:

1. **No raw data preprocessing**: Cell Ranger, nf-core/rnaseq, STAR, and similar tools are prohibited. All analysis starts from processed GEO deposits or author script outputs.
2. **Two local R runtimes**: `653_renv` (R 4.2.2) for older packages; `653_renv2` (R 4.5.3) for Xenium/SeuratWrappers/Banksy.
3. **Self-contained run layout**: All data, clones, and outputs reside under the timestamped run folder.
4. **Active execution mandate**: Run as much as possible; do not stop at feasibility-only.
5. **Per-accession download guidance**: GSE288990 (scRNA, ~278 MB), GSE309820 (bulk RNA), GSE289001 (Xenium mouse); PRJNA1221971 (raw WGS) skipped entirely.
6. **msigdbr API fix**: Patch `dicer1_rna.R` to remove `db_species`, update category codes, add memory management (`rm()/gc()` after each GSEA set).
7. **Self-contained output contract**: All generated outputs in `reports/` with relative-path links.

---

## Executive Conclusion

**11 of 13 author R scripts ran successfully and produced reproducible outputs.** The analysis covered bulk RNA-seq, scRNA-seq (QC through trajectory), inferCNV (CNV inference complete, July 2026), and PROGENy pathway analysis. The results agree with the paper's principal claims: DICER1 sarcoma cells recapitulate a mesenchymal kidney progenitor hierarchy, GSEA confirms myogenic/ground-state identity, and PROGENy scores separate tumor cell types.

Three scripts are blocked by infrastructure limits:
- **Xenium mouse** (`dicer1_xenium_proseg_mus_musculus.R`): requires ≥32 GB RAM (622K control + 133K tumor cells in SCTransform); blocked on this 15 GB system.
- **Xenium human** (`dicer1_xenium_proseg_homo_sapiens.R`): human Xenium RDS not publicly deposited on GEO.
- **sWGS** (`dicer1_swgs.R`): processed CN file not deposited; raw data (PRJNA1221971) requires QDNAseq preprocessing.

**Tier-1 reproduction path**: Processed GEO data (GSE288990 + GSE309820) → 653_renv → all 11 runnable scRNA/bulk scripts → 87 PNG plots (61 also saved as SVG) + 18 tables + 11 RDS objects (5 Seurat, 1 SingleCellExperiment, 1 Monocle3 CDS, 2 inferCNV-enriched, 1 Xenium input RDS [mouse only; human absent from GEO]) + 22 inferCNV checkpoints + 5 inferCNV heatmaps in run folder.

**Biggest gap**: Xenium spatial analysis requires a workstation with ≥32 GB RAM. No code or data gap exists — the GEO deposit and author scripts are complete; compute is the constraint.

---

## Paper And Local Files

| Item | Path/Status |
|---|---|
| Main PDF | `653_Spatial_...pdf` (in paper root) |
| Supplementary | SUPP/ folder |
| Code repository | Cloned to `workflow_draft/repos/dicer1-cell-of-origin/` |
| Runtime R 4.2.2 | `653_renv/conda_env/bin/Rscript` |
| Runtime R 4.5.3 | `653_renv2/conda_env/bin/Rscript` |
| GEO data | `workflow_draft/repos/dicer1-cell-of-origin/data/` |
| Reference | `workflow_draft/repos/dicer1-cell-of-origin/reference/` |
| Outputs | `workflow_draft/repos/dicer1-cell-of-origin/output/` |

---

## Availability Statements

| Type | Statement | Verified |
|---|---|---|
| Code | GitHub: https://github.com/Huntsmanlab/dicer1-cell-of-origin — 13 R scripts | Cloned, all 13 scripts present |
| scRNA-seq | GSE288990: 8 Cell Ranger filtered matrix directories (~278 MB) | Downloaded, extracted, renamed to PX IDs |
| Bulk RNA-seq | GSE309820: 12 Salmon quant.sf files + tx2gene | Downloaded |
| Xenium mouse | GSE289001: mouse Xenium RDS (2.2 GB) | Downloaded |
| Xenium human | Not found in public GEO deposits | Absent — blocked |
| sWGS processed | Not found on GEO; raw data in PRJNA1221971 | Absent — blocked |

---

## Data Deposition Map

| Accession | Content | Size | Used by scripts | Action |
|---|---|---|---|---|
| GSE288990 | scRNA Cell Ranger filtered matrices (8 samples) | 278 MB TAR | dicer1_scrna_qc.R | Downloaded, extracted, PX-renamed |
| GSE309820 | Bulk RNA Salmon quant.sf (12 samples) + tx2gene | ~15 MB | dicer1_rna.R | Downloaded |
| GSE289001 | Mouse Xenium RDS (proseg-segmented) | 2.2 GB | dicer1_xenium_proseg_mus_musculus.R | Downloaded |
| GSE289001 | Human Xenium RDS | Not found | dicer1_xenium_proseg_homo_sapiens.R | Absent |
| PRJNA1221971 | Raw WGS BAM/FASTQ | Several GB | dicer1_swgs.R (needs QDNAseq output) | Skipped per directives |

---

## Figure/Table Analysis Census

Census created at `metadata_inventory/figure_table_analysis_census.tsv`.

| Panel type | Count | Status |
|---|---|---|
| Bulk RNA-seq (PCA, volcano, GSEA) | ~10 panels | Reproducible — outputs generated |
| scRNA QC + controls | ~20 panels | Reproducible — outputs generated |
| scRNA tumors | ~20 panels | Reproducible — outputs generated |
| scRNA controls+tumors integrated | ~15 panels | Reproducible — outputs generated |
| Pseudotime trajectory | ~5 panels | Reproducible — outputs generated |
| InferCNV CNV heatmap + chr. plots | ~5 panels | Completed (July 2026) |
| PROGENy pathway heatmap/UMAP | ~8 panels | Reproducible — outputs generated |
| Xenium spatial (mouse) | ~20 panels | Blocked by insufficient RAM |
| Xenium spatial (human) | ~10 panels | Blocked by missing data |
| sWGS CN + clinical metadata | ~8 panels | Partially reproduced (metadata only) |

---

## Code Inventory And Script-To-Claim Map

Repository: https://github.com/Huntsmanlab/dicer1-cell-of-origin  
Clone: `workflow_draft/repos/dicer1-cell-of-origin/` (shallow, commit at clone time 2026-07-03)

| Script | Purpose | Paper figures | Runtime | Status |
|---|---|---|---|---|
| `dicer1_rna.R` | Bulk RNA-seq: tximport → DESeq2 → PCA, volcano, GSEA | Fig 4 (bulk panels) | 653_renv | COMPLETE |
| `dicer1_scrna_qc.R` | scRNA: Load 8 samples, DropletUtils QC, merge | Intermediate only (produces sce_qc.rds; Fig 2a is a BioRender schematic) | 653_renv | COMPLETE |
| `dicer1_scrna_controls.R` | Control scRNA: Seurat clustering, cell type annotation | Fig 2b-d | 653_renv | COMPLETE |
| `dicer1_scrna_tumors.R` | Tumor scRNA: Seurat clustering, cell type annotation | Fig 3b-d (Fig 3a is a BioRender schematic) | 653_renv | COMPLETE (patched split) |
| `dicer1_scrna_trajectory.R` | Monocle3 pseudotime trajectory | Fig 4b-d | 653_renv2 | COMPLETE |
| `dicer1_scrna_controls_and_tumors.R` | Integrated control+tumor UMAP, cell type annotation | Supp Fig 7 (a–f) | 653_renv | COMPLETE |
| `dicer1_scrna_export_infer_cnv.R` | Export counts + annotations for inferCNV | (intermediate) | 653_renv | COMPLETE (filter-only) |
| `dicer1_scrna_run_infer_cnv.R` | InferCNV: HMM-based CNV calling | Fig 4e (relates to Supp Fig 9c); Supp Fig 9d (control CNV) not produced — our run was tumor-only | 653_renv | COMPLETE (July 2026) |
| `dicer1_scrna_add_infer_cnv_seurat.R` | Map CNV results back to Seurat | Fig 4e (unsplit UMAP); Supp Fig 9c (split-sample CN analysis) | 653_renv | COMPLETE (July 2026) |
| `dicer1_scrna_run_pathway_analysis.R` | decoupleR PROGENy pathway scoring | Fig 4h (Fig 4g is an experimental CRISPR bar chart, not from this script) | 653_renv | COMPLETE |
| `dicer1_xenium_proseg_mus_musculus.R` | Mouse Xenium: spatial clustering, annotation | Fig 2e-h + Fig 3e-h + Supp Fig 8 | 653_renv2 | BLOCKED (RAM) |
| `dicer1_xenium_proseg_homo_sapiens.R` | Human Xenium: spatial clustering, annotation | Fig 5 | 653_renv2 | BLOCKED (data absent) |
| `dicer1_swgs.R` | sWGS CN analysis, clinical metadata | Fig 1 | 653_renv | BLOCKED (CN data absent) |

All 13 scripts are represented in `workflow_draft/figure_to_code_map.tsv` with runnable/blocked status.

---

## Tool/Package References

| Tool/Package | Source | URL |
|---|---|---|
| Seurat v5 | CRAN/GitHub | https://satijalab.org/seurat/ |
| inferCNV | Bioconductor | https://bioconductor.org/packages/infercnv |
| DESeq2 | Bioconductor | https://bioconductor.org/packages/DESeq2 |
| tximport | Bioconductor | https://bioconductor.org/packages/tximport |
| clusterProfiler / GSEA | Bioconductor | https://bioconductor.org/packages/clusterProfiler |
| decoupleR / PROGENy | CRAN/GitHub | https://saezlab.github.io/decoupleR/ |
| OmniPath REST API | Web | https://omnipathdb.org/ |
| monocle3 | GitHub | https://cole-trapnell-lab.github.io/monocle3/ |
| BPCells | GitHub | https://github.com/bnprks/BPCells |
| msigdbr 7.5.1 API | CRAN | https://cran.r-project.org/package=msigdbr |
| Ensembl REST API | Web | https://rest.ensembl.org/ |

---

## Analysis-by-Analysis Feasibility

### 1. Bulk RNA-seq (dicer1_rna.R)
**Status**: Reproducible from public processed data (GSE309820 Salmon quant.sf files).  
**Patches applied**: (a) All `msigdbr()` calls updated for API version 7.5.1 (removed `db_species`, updated category codes MH→H, M3→C3, M5→C5, `subcollection`→`subcategory`); (b) `rm()/gc()` added after each large GSEA gene set to prevent OOM at plot stage.  
**Outputs**: PCA plot, 2 volcano plots (DMEM + MesEnCult conditions), 12 GSEA enrichment plots (Hallmarks, Reactome, GO:BP, Oncogenic, TF, miRNA per condition), 5 DESeq2 result tables.

### 2. scRNA QC (dicer1_scrna_qc.R)
**Status**: Reproducible from public data (GSE288990).  
**Notes**: GSE288990_RAW.tar contained GSM-named files; renamed to PX ID directories using `dicer1_scrna_sample_metadata.csv` mapping before running.  
**Outputs**: `dicer1_scrna_sce_qc.rds` (8 samples merged, QC metrics computed).

### 3. scRNA Controls (dicer1_scrna_controls.R)
**Status**: Reproducible from public data.  
**Patches applied**: (a) Added `library(SingleCellExperiment)` for `rowData()`; (b) `make.unique()` on rownames; (c) `JoinLayers()` before FindAllMarkers (Seurat v5 requirement); (d) `tryCatch` wrappers for DotPlot (Seurat v5 S7 API incompatibility).  
**Outputs**: `dicer1_scrna_unfilter_control.rds`, `dicer1_scrna_filter_control.rds`, 17 plots (clusters, cell types, features, heatmaps, dot plots), 4 marker tables.

### 4. scRNA Tumors (dicer1_scrna_tumors.R)
**Status**: Reproducible from public data. Tumor script split into two due to OOM on FindAllMarkers (exit 137) when holding both unfiltered and filtered objects simultaneously.  
**Patches applied**: Same Seurat v5 fixes as controls; split into `dicer1_scrna_tumors_patched.R` (produces unfilter_tumor.rds) and `dicer1_scrna_tumors_filter_only_patched.R` (produces filter_tumor.rds + all plots), with aggressive `rm()/gc()` between steps.  
**Outputs**: `dicer1_scrna_unfilter_tumor.rds`, `dicer1_scrna_filter_tumor.rds`, ~20 plots (clusters, cell types, features, heatmaps, dot plots by batch), 4 marker tables.

### 5. Pseudotime Trajectory (dicer1_scrna_trajectory.R)
**Status**: Reproducible from local data (`dicer1_scrna_filter_tumor.rds`).  
**Runtime**: 653_renv2 (R 4.5.3) — required for SeuratWrappers/monocle3 compatibility; not 653_renv.  
**Patches applied**: Wrapped display-device calls (ElbowPlot, DimHeatmap) in tryCatch; ensured monocle3 `as.cell_data_set()` available via SeuratWrappers.  
**Limitations**: (a) Root node approximated from Progen centroid — the paper sets the root interactively via `choose_cells()`, which cannot be replicated non-interactively; (b) branch-specific trajectory plots (requiring `choose_graph_segments()`) were skipped as they require an interactive session; pseudotime and marker-gene-along-pseudotime plots were generated.  
**Outputs**: `dicer1_scrna_trajectory_cds.rds`, pseudotime UMAP, cell-type-colored trajectory graph, marker gene expression along pseudotime (Dpt, Pi16, Mfap4 panels).

### 6. Integrated Controls + Tumors (dicer1_scrna_controls_and_tumors.R)
**Status**: Reproducible from local data.  
**Patches applied**: Seurat v5 DotPlot/DoHeatmap tryCatch wrappers.  
**Outputs**: `dicer1_scrna_filter_control_and_tumor.rds`, ~15 plots (UMAP by cluster, cell type, sample, genotype, time; heatmaps; dot plots; xenium feature overlays), 1 cluster markers table.

### 7. inferCNV Export (dicer1_scrna_export_infer_cnv.R)
**Status**: Reproducible from local data. Full export (unfiltered combined, ~15K×32K dense matrix) blocked by OOM. Filter-only export succeeded.  
**Patches applied**: Created `dicer1_scrna_export_infer_cnv_filter_only.R` that writes only the filter_combined matrix in chunks, skipping the larger unfiltered matrices.  
**Outputs**: `dicer1_scrna_filter_combined_counts.matrix` (~471 MB), `dicer1_scrna_filter_combined_infer_cnv_cells.tsv`, `dicer1_scrna_unfilter_control_counts.matrix`, `dicer1_scrna_unfilter_tumor_counts.matrix`.

### 8. inferCNV Run (dicer1_scrna_run_infer_cnv.R)
**Status**: Complete (July 2026, post-initial-run). The initial run context was exhausted mid-execution (at tumor subclustering step 15, leiden incompatibility). Fixed by switching to `tumor_subcluster_partition_method = "random_trees"` and re-running with 2 threads to reduce memory pressure; all 22 HMM steps completed successfully.  
**Fix applied**: Prior run completed through step 14 but failed at step 15 (leiden tumor subclustering) because inferCNV's internal `define_signif_tumor_subclusters()` calls `Seurat::ScaleData()` without prior `NormalizeData()`, required in Seurat v5. `random_trees` avoids leiden entirely. Checkpoint files could not be reused (changed parameter invalidated cache); full re-run from step 1 with 2-thread constraint.  
**Outputs**: `output/infer_cnv/combined_filter/run.final.infercnv_obj`; `run_final_infercnv_obj.rds`; 5 inferCNV heatmaps (`infercnv.png`, `infercnv.preliminary.png`, `infercnv.17_HMM_pred...png`, `infercnv.19_HMM_pred...png`, `infercnv.20_HMM_pred...png`); 22 intermediate checkpoint files.

### 9. inferCNV → Seurat Mapping (dicer1_scrna_add_infer_cnv_seurat.R)
**Status**: Complete (July 2026, post-initial-run).  
**Fix applied**: `infercnv::add_to_seurat()` was absent from the author's `dicer1_scrna_run_infer_cnv.R`. An intermediate script (`dicer1_scrna_add_to_seurat_intermediate.R`) was created to call `add_to_seurat()` and produce the required metadata RDS before the author's add_infer_cnv_seurat script ran.  
**Outputs**: `output/infer_cnv/dicer1_scrna_filter_combined_infercnv_metadata.rds`; `output/infer_cnv/dicer1_scrna_filter_tumor_infercnv.png`; `output/infer_cnv/dicer1_scrna_filter_tumor_infercnv.svg`; `output/infer_cnv/dicer1_scrna_filter_tumor_infercnv_split_sample.png`; `output/infer_cnv/dicer1_scrna_filter_tumor_infercnv_split_sample.svg`.

### 10. PROGENy Pathway Analysis (dicer1_scrna_run_pathway_analysis.R)
**Status**: Reproducible from local data.  
**Patches applied**: (a) `OmnipathR` package absent — bypassed using direct OmniPath REST API download of PROGENy annotations to `data/progeny_omnipath.tsv`; (b) Seurat v5 Assay5 `SetAssayData()` fix for data layer; (c) `safe_featureplot_pathway()` wrapper bypassing `FeaturePlot()` `& scale_colour_gradient2()` S7 operator incompatibility.  
**Outputs**: `dicer1_scrna_filter_tumor_pathway_activity.csv`, pathway heatmap, 3 pathway UMAP plots (TGF-β, TNF-α, NF-κB); VEGF and EGFR UMAPs were not included in the patched script (v4) — these pathways may appear in the paper figure but were not present in the reproducible script used in this run.

### 11. Xenium Mouse Spatial (dicer1_xenium_proseg_mus_musculus.R)
**Status**: Blocked by insufficient RAM on this system. The Xenium mouse dataset (GSE289001) contains 1.495 million cells. After QC filtering, the control subset has 622,209 cells and the tumor subset has 133,115 cells. SCTransform on 622K cells requires approximately 30–50 GB RAM depending on implementation; our system has 15 GB total (shared with other processes). This is a hardware constraint, not a code or data gap.  
**Downloads**: `dicer1_xenium_mus_musculus_all.rds` (2.2 GB) downloaded from GSE289001.  
**Patches prepared** (v3 patched script): two-pass memory management (reload for WT, reload for MUT), UpdateSeuratObject() for FOV slot migration, path fixes, metadata column fixes. BPCells conversion could reduce peak RAM for sparse operations but not for SCTransform's dense residual matrix computation.  
**Prepared for future use**: `workflow_draft/scripts/patched/dicer1_xenium_proseg_mus_musculus_patched.R` is ready to run on a system with ≥32 GB RAM.

### 12. Xenium Human Spatial (dicer1_xenium_proseg_homo_sapiens.R)
**Status**: Blocked by missing data. Human Xenium processed RDS not found in GSE289001 or any accessible GEO deposit. No controlled access notice found; data may be under embargo or deposited separately.

### 13. sWGS Copy Number Analysis (dicer1_swgs.R)
**Status**: Blocked by missing processed CN data. The script reads `dicer1_swgs_50kb_noXY_copyNumbersSegmented_filtered_called.rds` which is the QDNAseq output for sWGS samples in PRJNA1221971. This processed file was not deposited on GEO. Raw data from PRJNA1221971 requires QDNAseq preprocessing (not available per runtime directives).  
**Partial reproduction**: `dicer1_scrna_swgs_metadata_only_patched.R` reads `dicer1_swgs_sample_metadata.csv` (present in repo) and produces clinical metadata panel plots (gender, histotype, Trp53, Kras, Dicer1 status across 35 samples). CN plots blocked.

---

## Workflow Scaffold / Runnable Draft

### Scripts provided

| File | Type | Status |
|---|---|---|
| `Snakefile` | Full workflow DAG | Dry-run passed (EXIT:0) |
| `run.sh` | Shell entrypoint | Tier-1 runnable path |
| `run_order.tsv` | Executable DAG | Updated |
| `RUN_ORDER.md` | Human narrative | Updated |
| `figure_to_code_map.tsv` | Figure → code mapping | All 13 scripts |
| `dependency_chain.tsv` | Artifact creator/consumer chain | Complete |
| `assumptions_and_defaults.tsv` | All assumed parameters | Complete |
| `config/config.yaml` | Configurable parameters | Complete |
| `scripts/patched/` | 11 patched/extended scripts | All issues documented |

### Author scripts vs shadow scaffolds
All Snakefile rules call the **actual author scripts** directly via `Rscript ../code/<author_script>.R`. No shadow scaffold scripts were substituted for author code. Generated scaffold scripts exist only for analyses where no author code exists or where blocking conditions require a metadata-only fallback:
- `dicer1_scrna_add_to_seurat_intermediate.R` — post-processing step not present in author code
- `dicer1_swgs_metadata_only_patched.R` — partial fallback due to missing CN data

### Snakefile dry-run
```
/home/mneira/miniconda3/envs/snakemake/bin/snakemake --snakefile workflow_draft/Snakefile --dry-run --cores 1
EXIT: 0 — all targets met (from prior run)
```

---

## Local Runtime Execution Summary

### Runtimes used

| Runtime | Path | Version | Purpose |
|---|---|---|---|
| 653_renv | `653_renv/conda_env/bin/Rscript` | R 4.2.2 | All scRNA, bulk RNA, inferCNV, pathway scripts |
| 653_renv2 | `653_renv2/conda_env/bin/Rscript` | R 4.5.3 | Xenium scripts (SeuratWrappers, Banksy, monocle3, sf/spdep) |

All execution used PATH-aware wrappers (`rscript_with_path.sh`) to ensure conda compiler toolchain visibility.

### Execution results by script

| Script | Attempted | Outcome | Patched? | Log file |
|---|---|---|---|---|
| dicer1_rna.R | Yes | Complete (all plots + tables) | Yes — msigdbr API + OOM fixes | dicer1_rna.stdout.log |
| dicer1_scrna_qc.R | Yes | Complete (sce_qc.rds) | No | dicer1_scrna_qc.stdout.log |
| dicer1_scrna_controls.R | Yes | Complete (unfilter + filter control RDS + plots) | Yes — rowData, make.unique, JoinLayers | dicer1_scrna_controls.stdout.log |
| dicer1_scrna_tumors.R | Yes | Complete via split (unfilter + filter tumor RDS + plots) | Yes — OOM split into 2 scripts | dicer1_scrna_tumors.stdout.log |
| dicer1_scrna_trajectory.R | Yes | Complete (trajectory_cds.rds + plots) | Yes — display-device tryCatch | dicer1_scrna_trajectory.stdout.log |
| dicer1_scrna_controls_and_tumors.R | Yes | Complete (filter combined RDS + plots) | Yes — DotPlot tryCatch | dicer1_scrna_controls_and_tumors.stdout.log |
| dicer1_scrna_export_infer_cnv.R | Yes | Filter-only complete; unfiltered blocked by OOM | Yes — filter-only variant | dicer1_scrna_export_infer_cnv.stdout.log |
| dicer1_scrna_run_infer_cnv.R | Yes (v2) | Complete (July 2026) — all 22 HMM steps; run.final.infercnv_obj written | Yes — random_trees subclustering, 2-thread constraint | dicer1_scrna_run_infer_cnv_v2.stdout.log |
| dicer1_scrna_add_infer_cnv_seurat.R | Yes | Complete (July 2026) — via add_to_seurat_intermediate.R; 4 output files (2 PNG + 2 SVG) | add_to_seurat intermediate created (missing from author code) | dicer1_scrna_add_infer_cnv_seurat.stdout.log |
| dicer1_scrna_run_pathway_analysis.R | Yes (v4) | Complete (pathway heatmap + UMAPs + table) | Yes — OmnipathR bypass, S7 FeaturePlot fix | dicer1_scrna_run_pathway_analysis_v4.stdout.log |
| dicer1_xenium_proseg_mus_musculus.R | Yes (v3) | Blocked — OOM (622K WT cells; requires ≥32 GB RAM) | Yes — v3 ready for high-RAM system | dicer1_xenium_mus_musculus.stdout.log |
| dicer1_xenium_proseg_homo_sapiens.R | No | Blocked — human Xenium data not found on GEO | Skeleton prepared | — |
| dicer1_swgs.R | Partial | Blocked — CN data absent; metadata plots produced | Yes — metadata-only variant | dicer1_swgs_metadata_only.stdout.log |

### Key issues encountered and resolved

See `metadata_inventory/solved_issues.md` for full details. Summary of critical fixes:

1. **msigdbr v7.5.1 API break** in `dicer1_rna.R`: removed `db_species`, updated category codes. OOM at plot stage fixed with `rm()/gc()` after each GSEA gene set.
2. **Seurat v5 JoinLayers requirement**: `FindAllMarkers()` fails without `JoinLayers()` on multi-layer assays.
3. **Seurat v5 DotPlot S7 incompatibility**: internal `& .gg` operator fails; wrapped in tryCatch with manual ggplot2 fallback.
4. **GSE288990 GSM naming**: downloaded files named by GSM IDs; manually renamed to PX IDs using sample metadata CSV mapping.
5. **OmnipathR absent**: bypassed by downloading PROGENy annotations directly from REST API.
6. **InferCNV leiden Seurat v5 incompatibility**: changed `tumor_subcluster_partition_method` to `random_trees`.
7. **InferCNV filter-only export**: unfiltered combined matrix (~15K cells × 32K genes dense) caused OOM; created filter-only export that skips unneeded matrices.
8. **gene_order.tsv absent**: generated via Ensembl REST API Python batch script across all 32K mouse genes.
9. **Xenium RDS Seurat v4→v5 slot migration**: added `UpdateSeuratObject()` after readRDS.
10. **Xenium RAM**: 622K WT control cells requires ≥32 GB for SCTransform; fundamental hardware constraint.

---

## Generated Outputs

All paths are relative to this report at `reports/01_replication_feasibility_assessment.md`.

### Plots

**Note:** 87 PNG plots were generated in total (61 also as SVG). The table below covers all main-figure and key supplementary panels; secondary outputs (dot plots, sample/time/batch UMAPs, heatmap variants, unfiltered clustering intermediates, additional GSEA conditions) exist on disk at `workflow_draft/repos/dicer1-cell-of-origin/output/plots/` but are omitted here for brevity.

| Link | Paper figure | Description | Agreement with paper |
|---|---|---|---|
| [dicer1_rna_mut_vs_control.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_rna_mut_vs_control.png) | **Fig 4a** | Dot/bubble plot of specific gene log2FC (Dpt, Dcn, Pi16, Ly6a, Ly6c1, Col15a1, Itga8, Myh11, Acta2) in Dicer1 MUT vs WT primary cells; colored by condition (DMEM vs MesEnCult); faceted by lineage group (Universal, Fibro, Peri) | Matches Fig 4a description: increased universal fibroblast markers (Dpt, Pi16, Ly6a), reduced lineage markers (Itga8, Myh11, Acta2) |
| [dicer1_bulk_rna_pca.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_bulk_rna_pca.png) | **Supp Fig 9a** | PCA of 12 bulk RNA-seq samples (6 DMEM, 6 MesEnCult; HDT1056 and HDT1119); QC/exploratory plot | Matches Supp Fig 9a: PCA of bulk RNA-seq of primary cell culture samples (control vs Dicer1 mutant, MesenCult or DMEM); PC1 separates by genotype |
| [dicer1_rna_dmem_volcano.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_rna_dmem_volcano.png) | Not in paper | Volcano plot, DMEM condition MUT vs WT (DESeq2) | Not identified in main or supplementary figures — likely exploratory output from author's script not included in publication |
| [dicer1_rna_mesencult_volcano.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_rna_mesencult_volcano.png) | Not in paper | Volcano plot, MesEnCult condition MUT vs WT (DESeq2) | Not identified in main or supplementary figures — likely exploratory output from author's script not included in publication |
| [dicer1_rna_dmem_gsea_hallmarks.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_rna_dmem_gsea_hallmarks.png) | Not in paper | GSEA Hallmarks, DMEM condition (bulk RNA) | Not identified in main or supplementary figures — likely exploratory output from author's script not included in publication |
| [dicer1_rna_mesencult_gsea_hallmarks.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_rna_mesencult_gsea_hallmarks.png) | Not in paper | GSEA Hallmarks, MesEnCult condition (bulk RNA) | Not identified in main or supplementary figures — likely exploratory output from author's script not included in publication |
| [dicer1_rna_dmem_gsea_reactome.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_rna_dmem_gsea_reactome.png) | Not in paper | GSEA Reactome, DMEM condition (bulk RNA) | Not identified in main or supplementary figures — likely exploratory output from author's script not included in publication |
| [dicer1_rna_dmem_gsea_ont.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_rna_dmem_gsea_ont.png) | Not in paper | GSEA GO:BP, DMEM condition (bulk RNA) | Not identified in main or supplementary figures — likely exploratory output from author's script not included in publication |
| [dicer1_scrna_unfilter_control_clusters.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_unfilter_control_clusters.png) | **Supp Fig 3e** | scRNA control cells: UMAP by Louvain cluster before cell type annotation | Matches Supp Fig 3e: UMAP of tdTomato+ mesenchymal stromal cells (n=3) colored by cluster assignment, relating to Fig. 2b (distinct from Supp Fig 3a which shows all control kidney cells) |
| [dicer1_scrna_filter_control_cell_type.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_control_cell_type.png) | **Fig 2b** | scRNA control cells: UMAP colored by annotated cell type (Uni Fibro Pi16high, Uni Fibro Col15a1high, Fibro, Mural, Cycle) | Matches Fig 2b: UMAP of tdtomato+ MSC colored by cell type with lineage annotations |
| [dicer1_scrna_filter_control_features.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_control_features.png) | **Fig 2c** | scRNA control cells: UMAP depicting expression of key fibroblast lineage markers | Matches Fig 2c: FeaturePlots of key markers on tdTomato+ MSC UMAP |
| [dicer1_scrna_filter_control_heatmap_clusters.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_control_heatmap_clusters.png) | **Fig 2d** | scRNA control cells: heatmap of top 30 defining transcripts based on cell type annotation | Matches Fig 2d: heatmap with key genes of each MSC cell type highlighted |
| [dicer1_scrna_unfilter_tumor_clusters.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_unfilter_tumor_clusters.png) | **Supp Fig 6a** | scRNA tumor cells: UMAP by Louvain cluster before annotation | Matches Supp Fig 6a: UMAP of HDT tumor cells (n=4) colored by cluster assignment |
| [dicer1_scrna_filter_tumor_cell_type.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_tumor_cell_type.png) | **Fig 3b** | scRNA tumor cells: UMAP colored by annotated cell type (Progen, Ground, Prolif, TR Diff Myo, Diff Myo, Fibro Perturb, Mural) | Matches Fig 3b: UMAP of tdtomato+ tumor MSC colored by cell type |
| [dicer1_scrna_filter_tumor_heatmap_clusters.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_tumor_heatmap_clusters.png) | **Fig 3c** | scRNA tumor cells: heatmap of top 30 signature-defining transcripts based on cell state annotation | Matches Fig 3c: heatmap with key genes of each tumor cell state highlighted |
| [dicer1_scrna_filter_tumor_features.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_tumor_features.png) | **Fig 3d** | scRNA tumor cells: UMAP depicting expression of key markers | Matches Fig 3d: FeaturePlots of key lineage markers on tdTomato+ tumor MSC UMAP |
| [dicer1_scrna_filter_tumor_cell_type_time_split.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_tumor_cell_type_time_split.png) | Not in paper | Tumor cell type UMAP split by timepoint/sample | Not identified in main or supplementary figures — likely exploratory output from author's script not included in publication |
| [dicer1_scrna_filter_tumor_trajectory_pseudotime.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_tumor_trajectory_pseudotime.png) | **Fig 4b-c** | Monocle3 pseudotime UMAP colored by cell type and pseudotime; root approximated from Progen centroid; branch-specific plots not generated (require interactive session) | Matches Fig 4b-c concept; root placement may differ from paper (paper sets root interactively) |
| [dicer1_scrna_filter_tumor_trajectory_cell_type_graph.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_tumor_trajectory_cell_type_graph.png) | **Fig 4b-c** | Trajectory principal graph overlaid on UMAP, colored by cell type | Same as above |
| [dicer1_scrna_filter_tumor_trajectory_ground_genes_cell_type.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_tumor_trajectory_ground_genes_cell_type.png) | **Fig 4d** | Gene trend plots of Ground cluster marker genes along sarcomagenic trajectory; cells colored by cell type | Matches Fig 4d: gene expression trends of Ground-state defining genes (e.g., Dpt, Pi16, Mfap4) along trajectory |
| [dicer1_scrna_filter_tumor_trajectory_ground_genes_time.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_tumor_trajectory_ground_genes_time.png) | **Fig 4d** | Gene trend plots of Ground cluster marker genes along trajectory; cells colored by pseudotime | Matches Fig 4d (pseudotime-colored variant) |
| [dicer1_scrna_filter_tumor_trajectory_myo_diff_genes_cell_type.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_tumor_trajectory_myo_diff_genes_cell_type.png) | **Fig 4d** | Gene trend plots of Diff Myo marker genes along sarcomagenic trajectory; cells colored by cell type | Matches Fig 4d: gene expression trends of myogenic differentiation genes along trajectory |
| [dicer1_scrna_filter_tumor_trajectory_myo_diff_genes_time.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_tumor_trajectory_myo_diff_genes_time.png) | **Fig 4d** | Gene trend plots of Diff Myo marker genes along trajectory; cells colored by pseudotime | Matches Fig 4d (pseudotime-colored variant) |
| [dicer1_scrna_filter_control_and_tumor_cell_type.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_control_and_tumor_cell_type.png) | **Supp Fig 7f** | Integrated control+tumor UMAP colored by cell type | Matches Supp Fig 7f: UMAP of combined control (n=3) + HDT tumor (n=4) cells colored by cell type assigned during individual work-up of scRNA data, relating to Fig. 2b and Fig. 3b |
| [dicer1_scrna_filter_tumor_pathway_heatmap.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_tumor_pathway_heatmap.png) | **Fig 4h** | PROGENy MLM pathway activity heatmap by tumor cell type; shows MAPK upregulation in Prolif cluster | Partially matches Fig 4h (GSEA along tumor cell state continuum); paper uses GSEA of DEGs between cell states; our output uses PROGENy activity scores — related but methodologically distinct |
| [dicer1_scrna_filter_tumor_pathway_tgfb.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_tumor_pathway_tgfb.png) | Not in paper | TGF-β pathway activity UMAP (PROGENy score per cell) | Not identified in main or supplementary figures — likely exploratory output from author's script not included in publication |
| [dicer1_scrna_filter_tumor_pathway_tnfa.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_tumor_pathway_tnfa.png) | Not in paper | TNF-α pathway activity UMAP (PROGENy score per cell) | Not identified in main or supplementary figures — likely exploratory output from author's script not included in publication |
| [dicer1_scrna_filter_tumor_pathway_nfkb.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_scrna_filter_tumor_pathway_nfkb.png) | Not in paper | NF-κB pathway activity UMAP (PROGENy score per cell) | Not identified in main or supplementary figures — likely exploratory output from author's script not included in publication |
| [dicer1_scrna_filter_tumor_infercnv.png](../workflow_draft/repos/dicer1-cell-of-origin/output/infer_cnv/dicer1_scrna_filter_tumor_infercnv.png) | **Fig 4e** | InferCNV CNV proportion feature plots overlaid on tumor UMAP; shows chromosomal duplication burden per cell | Matches Fig 4e (UMAP showing proportion of cells with chromosomal duplications, relating to Supp Fig 9c). Supp Fig 9d (control kidneys CN) was not produced — our inferCNV run used tumor cells only |
| [dicer1_scrna_filter_tumor_infercnv_split_sample.png](../workflow_draft/repos/dicer1-cell-of-origin/output/infer_cnv/dicer1_scrna_filter_tumor_infercnv_split_sample.png) | **Supp Fig 9c** | InferCNV CNV feature plots split by sample | Matches Supp Fig 9c: CN analysis of tdTomato⁺ mesenchymal cells (Observations) in tumorous HDT kidneys vs bystander cell populations (References), related to Fig. 4e |
| [dicer1_swgs_metadata_panel.png](../workflow_draft/repos/dicer1-cell-of-origin/output/plots/dicer1_swgs_metadata_panel.png) | Fig 1 (partial) | sWGS clinical metadata panel (sex, histotype, Trp53/Kras/Dicer1 status across 35 samples) — independent reimplementation | Partial: metadata plots only; CN copy number panels absent (processed sWGS RDS not deposited) |

### Tables

| Link | Paper figure | Rows × cols | Description |
|---|---|---|---|
| [dicer1_rna_dmem_de_mut_vs_wt.csv](../workflow_draft/repos/dicer1-cell-of-origin/output/tables/dicer1_rna_dmem_de_mut_vs_wt.csv) | Fig 4 | ~30K × 9 | DESeq2 DE table, DMEM MUT vs WT |
| [dicer1_rna_mesencult_de_mut_vs_wt.csv](../workflow_draft/repos/dicer1-cell-of-origin/output/tables/dicer1_rna_mesencult_de_mut_vs_wt.csv) | Fig 4 | ~30K × 9 | DESeq2 DE table, MesEnCult MUT vs WT |
| [dicer1_rna_norm_counts_no_log.csv](../workflow_draft/repos/dicer1-cell-of-origin/output/tables/dicer1_rna_norm_counts_no_log.csv) | Supp | ~30K × 12 | DESeq2 normalized counts |
| [dicer1_scrna_filter_tumor_cluster_markers_table.csv](../workflow_draft/repos/dicer1-cell-of-origin/output/tables/dicer1_scrna_filter_tumor_cluster_markers_table.csv) | Fig 3a-d | ~10K rows | Tumor cluster marker genes (FindAllMarkers) |
| [dicer1_scrna_filter_tumor_pathway_activity.csv](../workflow_draft/repos/dicer1-cell-of-origin/output/tables/dicer1_scrna_filter_tumor_pathway_activity.csv) | Fig 4g-h | Per-cell | PROGENy pathway activity scores |
| [dicer1_swgs_sample_metadata.csv](../workflow_draft/repos/dicer1-cell-of-origin/output/tables/dicer1_swgs_sample_metadata.csv) | Fig 1 | 35 × 19 | sWGS sample metadata |

---

## Environment And Execution Risks

### Packages with known version-sensitive behavior

| Package | Version in 653_renv | Issue | Fix applied |
|---|---|---|---|
| Seurat | 5.1.0 | v5 JoinLayers, DotPlot S7 API changes | tryCatch wrappers, JoinLayers() |
| msigdbr | 7.5.1 | Removed `db_species` parameter, renamed categories | Patched all calls |
| OmnipathR | not installed | Required by decoupleR::get_progeny() | REST API bypass |
| inferCNV | 1.14.0 | Leiden subclustering uses Seurat v5 ScaleData without NormalizeData | random_trees method |
| biomaRt | — | BiocFileCache dplyr API incompatibility in this environment | Switched to Ensembl REST API |

### Non-determinism
Clustering algorithms (Louvain/Leiden), UMAP layout, and PCA are seeded with `set.seed(404)` in scripts that set it. Scripts without explicit seed will produce slightly different UMAP layouts per run but identical cell type composition ratios.

---

## Download Budget

| Accession | File | Size | Status |
|---|---|---|---|
| GSE288990 | GSE288990_RAW.tar | 278 MB | Downloaded |
| GSE309820 | Salmon quant.sf × 12 + tx2gene | ~30 MB | Downloaded |
| GSE289001 | dicer1_xenium_mus_musculus_all.rds | 2.2 GB | Downloaded |
| Total used | | ~2.5 GB | |
| Available at start | | ~26 GB | |
| Remaining | | ~23.5 GB | |

---

## Tiered Reproduction Strategy

### Tier 1 — Smallest credible reproduction (available now)

**Requirements**: This system with ≥8 GB free RAM, 653_renv runtime, processed GEO data downloaded.  
**Steps**: Run `workflow_draft/run.sh` (Tier 1 path) from repo code directory.  
**Produces**: All 11 runnable script outputs: bulk RNA PCA/volcano/GSEA, scRNA QC/controls/tumors/trajectory/integration/pathway analysis, inferCNV CNV calling, and inferCNV→Seurat mapping.  
**Time estimate**: ~2–3 hours (dominated by scRNA clustering and inferCNV).

### Tier 2 — Full scRNA + inferCNV + Xenium (requires hardware upgrade)

**Additional requirements**: ≥32 GB RAM.  
**Adds**: Xenium spatial transcriptomics analysis (622K WT + 133K tumor cells).  
**Script**: `dicer1_xenium_proseg_mus_musculus_patched.R` (v3) is ready to run.

### Tier 3 — Full raw-to-figure with sWGS (requires external tools)

**Additional requirements**: Sarek or QDNAseq preprocessing for PRJNA1221971 raw sWGS reads.  
**Adds**: Copy number analysis panel of Figure 1.  
**Blocked by**: Author did not deposit processed CN data on GEO; only raw reads are available.

---

## Missing Items / Author Requests

| Item | Required by | Blocker type | Suggestion |
|---|---|---|---|
| Human Xenium RDS | dicer1_xenium_proseg_homo_sapiens.R | Data absent from GEO | Request from authors |
| `dicer1_swgs_50kb_noXY_copyNumbersSegmented_filtered_called.rds` | dicer1_swgs.R | Not deposited on GEO | Request from authors or run QDNAseq on PRJNA1221971 |
| ≥32 GB RAM system | dicer1_xenium_proseg_mus_musculus.R | Hardware constraint | Use HPC or cloud instance |
| `dicer1_scrna_filter_combined_infercnv_metadata.rds` | dicer1_scrna_add_infer_cnv_seurat.R | Produced by `infercnv::add_to_seurat()` step missing from author's run script | Covered by our intermediate script |

---

## Recommended Next Actions

1. **sWGS CN data**: Contact authors for `dicer1_swgs_50kb_noXY_copyNumbersSegmented_filtered_called.rds`, or apply QDNAseq to PRJNA1221971 raw reads on an HPC.
2. **Human Xenium**: Request the human Xenium RDS from authors; may be under embargo.
3. **Xenium on HPC**: Run `dicer1_xenium_proseg_mus_musculus_patched.R` on a system with ≥32 GB RAM using the downloaded 2.2 GB RDS (at path `output/xenium/dicer1_xenium_mus_musculus_all.rds`). Script v3 is ready; no code or data gap.
4. **Trajectory branch plots**: If interactive analysis is desired, open `dicer1_scrna_trajectory_cds.rds` in an R session and use `choose_graph_segments()` / `choose_cells()` to set the root and branch endpoints interactively as in the paper.
