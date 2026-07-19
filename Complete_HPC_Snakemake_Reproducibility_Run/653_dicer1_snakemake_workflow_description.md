# DICER1 Sarcoma — Snakemake Reproducibility Workflow

**Paper:** "Spatial single-cell transcriptomic analysis informs tumor developmental
hierarchy of DICER1 syndrome-related sarcoma", *Nature Communications* 2026  
**Author code repository:** `github.com/Huntsmanlab/dicer1-cell-of-origin`  
**Workflow location (fir HPC):** `.../653_dicer1_snakemake/clean_snakemake_run/`  
**Completed:** 2026-07-18, SLURM job 49517968 (ExitCode 0, 03:11:37)

---

## Overview

This Snakemake workflow reproduces all computationally feasible analyses from the
DICER1 sarcoma paper using the author's original scripts as the primary input. The
workflow covers single-cell RNA-seq processing, bulk RNA-seq differential expression,
pseudotime trajectory inference, pathway activity scoring, copy number inference
(inferCNV), shallow whole-genome sequencing (sWGS) metadata and CN approximation,
and Xenium spatial transcriptomics.

The pipeline comprises **15 rules** (14 analysis rules + `rule all`), executed
sequentially in a single SLURM job on the fir HPC cluster (Alliance Canada).
Total wall time: ~11 hours across two runs (first required debugging; final clean
run: 3 h 11 min).

---

## Input Files

### Author code (unmodified originals in `workdir/code/`)
The author's R scripts are used directly from a clone of the `dicer1-cell-of-origin`
repository. Where a script required modification, a patched copy is stored in
`scripts/patched/` and the original in `code/` is left untouched.

### Data inputs (GEO: GSE267523 and GSE268427)
Raw processed data downloaded from GEO and symlinked into `workdir/data/`:

| Type | Samples | Details |
|---|---|---|
| scRNA-seq (10x Chromium) | 8 samples (PX prefix) | barcode/feature/matrix trios per sample |
| Bulk RNA-seq (Salmon quant) | 12 samples (hdt prefix) | 6 HDT1056 + 6 HDT1119, 2 conditions (DMEM/MesenCult); includes `tx2gene.tsv` (see below) |
| Xenium spatial transcriptomics | 1 combined RDS | `dicer1_xenium_mus_musculus_all.rds` (2.2 GB, symlink) |

### Metadata (in `workdir/metadata/`)
Three CSV files shipped with the author's repository:
- `dicer1_scrna_sample_metadata.csv` — scRNA sample annotations
- `dicer1_rna_metadata.csv` — bulk RNA sample and condition metadata
- `dicer1_swgs_sample_metadata.csv` — sWGS patient sample metadata (35 samples)

### Reference files (in `workdir/reference/`)
- `gene_order.tsv` — chromosome-ordered gene list for inferCNV
- `cell_type_markers_mus_musculus_kidney.yaml` — kidney cell-type marker gene sets
- `cell_type_markers_mus_musculus_mesenchymal.yaml` — mesenchymal marker gene sets
- `cell_type_markers_homo_sapiens.yaml` / `..._mesenchymal.yaml` — human equivalents

### Supplementary data (in `workdir/supplementary_data/`, symlink)
Publisher supplementary Excel files downloaded from the paper's data availability page:
- `source_data.xlsx` — source data for paper figures, including CopyNumberSegments
  output (sheet "Fig. S1k,l") used for sWGS CN frequency plots
- `data03_41467_2026_70971_MOESM5_ESM.xlsx` — Supplementary Data 3: per-sample
  Chr1/Chr6/Chr2/Chr15 CN status calls used for the Fisher test approximation
- `data01`–`data05` xlsx — additional supplementary tables

---

## R Environments

Two conda environments are used to match the author's package requirements:

| Wrapper | Environment | R version | Used for |
|---|---|---|---|
| `rscript_a` | `653_renv` | R 4.2.2 | bulk RNA, scRNA (all stages), inferCNV, pathway analysis, sWGS |
| `rscript_b` | `dicer1_r453` | R 4.5.3 | trajectory (Monocle3), Xenium (Seurat v5 + BPCells) |

Both wrappers use `env -i` isolation to prevent MUGQIC Python 3.13 on fir from
intercepting conda and snakemake commands.

---

## Workflow Rules

### 1. `scrna_qc`
**Script:** `code/dicer1_scrna_qc.R` (author original, unmodified)  
**Inputs:** 8 scRNA barcodes/features/matrix trios + sample metadata CSV  
**Output:** `output/scrna/dicer1_scrna_sce_qc.rds` (346 MB)  
QC filtering of raw 10x data into a SingleCellExperiment object. Applies per-sample
empty droplet detection, doublet scoring, and quality thresholds.

### 2. `scrna_controls`
**Script:** `scripts/patched/dicer1_scrna_controls_patched.R`  
**Input:** `dicer1_scrna_sce_qc.rds`  
**Outputs:** `dicer1_scrna_unfilter_control.rds` (245 MB), `dicer1_scrna_filter_control.rds` (36 MB)  
Processes control (WT) samples through SCTransform normalization, PCA, UMAP,
clustering, and cell-type annotation using kidney/mesenchymal marker YAML files.

### 3. `scrna_tumors_unfiltered`
**Script:** `scripts/patched/dicer1_scrna_tumors_patched.R`  
**Input:** `dicer1_scrna_sce_qc.rds`  
**Output:** `dicer1_scrna_unfilter_tumor.rds` (527 MB)  
Equivalent pipeline for tumor samples prior to cluster-based filtering.

### 4. `scrna_tumors_filtered`
**Script:** `scripts/patched/dicer1_scrna_tumors_filter_only_patched.R`  
**Input:** `dicer1_scrna_unfilter_tumor.rds`  
**Output:** `dicer1_scrna_filter_tumor.rds` (82 MB)  
Removes low-quality clusters; applies cell-type annotation to retained tumor cells.

### 5. `scrna_controls_and_tumors`
**Script:** `scripts/patched/dicer1_scrna_controls_and_tumors_patched.R`  
**Inputs:** `filter_control.rds` + `filter_tumor.rds`  
**Output:** `dicer1_scrna_filter_control_and_tumor.rds` (162 MB)  
Integrates control and tumor Seurat objects for joint downstream analysis.

### 6. `bulk_rna`
**Script:** `scripts/patched/dicer1_rna_patched.R`  
**Inputs:** Salmon quantification files (12 samples, 2 conditions) + `tx2gene.tsv` + metadata  
**Outputs:** `dicer1_rna_dmem_de_mut_vs_wt.csv`, `dicer1_rna_mesencult_de_mut_vs_wt.csv`,
`dicer1_rna_norm_counts_no_log.csv`  
Differential expression analysis via DESeq2 (tximport → DESeq2 → results tables).
Separate comparisons for DMEM and MesenCult culture conditions.

### 7. `scrna_trajectory`
**Script:** `scripts/patched/dicer1_scrna_trajectory_patched.R` (uses `rscript_b`)  
**Input:** `filter_tumor.rds`  
**Outputs:** `dicer1_scrna_trajectory_cds.rds`, `dicer1_scrna_filter_tumor_trajectory_pseudotime.png`  
Pseudotime trajectory inference using Monocle3. Requires R 4.5.3 (dicer1_r453)
for compatibility with Monocle3's SeuratWrappers interface.

### 8. `scrna_pathway_analysis`
**Script:** `scripts/patched/dicer1_scrna_run_pathway_analysis_patched.R`  
**Input:** `filter_tumor.rds`  
**Outputs:** `dicer1_scrna_filter_tumor_pathway_activity.csv`, `dicer1_scrna_filter_tumor_pathway_heatmap.png`  
PROGENy pathway activity scoring via the decoupleR package (OmnipathR backend).
Produces per-cell pathway activity scores and a summary heatmap for 14 cancer pathways.

### 9. `scrna_export_infercnv`
**Script:** `scripts/patched/dicer1_scrna_export_infer_cnv_filter_only.R`  
**Inputs:** `unfilter_control.rds` + `filter_tumor.rds` + features TSV  
**Outputs:** `dicer1_scrna_filter_combined_counts.matrix`, cell annotation TSV  
Extracts raw count matrix and cell annotations from the integrated Seurat object
in the format required by inferCNV.

### 10. `scrna_run_infercnv`
**Script:** `scripts/patched/dicer1_scrna_run_infer_cnv_v2_patched.R`  
**Inputs:** counts matrix + annotations + `gene_order.tsv`  
**Output:** `infer_cnv/combined_filter/run_final_infercnv_obj.rds`  
Runs inferCNV to infer copy number alterations from scRNA-seq data, using
normal (control) cells as reference. Produces the complete inferCNV object
including MCMC-smoothed CN profiles.

### 11. `scrna_add_to_seurat_intermediate`
**Script:** `scripts/patched/dicer1_scrna_add_to_seurat_intermediate.R`  
**Inputs:** inferCNV RDS + `unfilter_control.rds` + `filter_tumor.rds`  
**Output:** `dicer1_scrna_filter_combined_infercnv_metadata.rds`  
Extracts inferCNV predictions per cell (malignant/non-malignant calls) and
attaches them as metadata to the Seurat objects.

### 12. `scrna_add_infercnv_seurat`
**Script:** `scripts/patched/dicer1_scrna_add_infer_cnv_seurat_patched.R`  
**Inputs:** inferCNV metadata RDS + `filter_tumor.rds` + `unfilter_control.rds`  
**Output:** `infer_cnv/dicer1_scrna_filter_tumor_infercnv.svg`  
Generates the final inferCNV visualization panel — UMAP colored by inferCNV
malignancy score, overlaid on tumor cell clusters.

### 13. `swgs`
**Script:** `scripts/patched/dicer1_swgs_patched.R` (new, merged from two originals)  
**Inputs:** `dicer1_swgs_sample_metadata.csv` + `source_data.xlsx` + `data03` Excel  
**Outputs:** metadata CSV + metadata panel PNG, Fisher test CSV, CN frequency plots (SARC/LGMT)

This rule merges two originally separate analyses: the sWGS sample metadata panel
(authors' `dicer1_swgs_metadata_only.R`) and the CN approximation analysis
(developed during this reproducibility study as a substitute for unavailable RDS data;
see *Limitations* below). Produces the Chr1 gain frequency plots for SARC-like and
LGMT-like samples, and a Fisher's exact test for Chr1 gain enrichment.

### 14 + 15. `xenium_mouse` + `xenium_mouse_recover`
**Scripts:** `scripts/patched/dicer1_xenium_proseg_mus_musculus_patched.R` (28 KB)
and `scripts/patched/dicer1_xenium_recover_celltype_plots.R` (new)  
**Input:** `dicer1_xenium_mus_musculus_all.rds` (2.2 GB, via symlink)  
**Outputs:** `dicer1_xenium_wt.rds` (3.9 GB), `dicer1_xenium_mut.rds` (3.8 GB),
markers CSV, DoHeatmap + DotPlot by cell type (PNG + SVG)

Processes the Xenium mouse spatial transcriptomics data: ProSeg segmentation
output → Seurat v5 object → SCTransform → PCA/UMAP → clustering → cell-type
annotation using the kidney YAML markers. Separate analysis for WT and MUT
(DICER1-mutant) samples. Requires R 4.5.3 and BPCells for on-disk matrix
representation of the large cell count (~600K cells).

`xenium_mouse_recover` is a downstream rule that regenerates two paper-relevant
figures that fail in the main xenium script due to data-structure issues in the
marker gene lists (duplicate gene names, NA cell labels). It runs after
`xenium_mouse` completes and uses the saved `mut.rds` and markers CSV as inputs.

---

## Adaptations Made for the HPC Environment

All changes from the authors' original scripts are confined to `scripts/patched/`.
The originals in `code/` are unmodified. Key adaptations:

| Adaptation | Reason |
|---|---|
| Two Rscript wrappers (`rscript_a_fir.sh`, `rscript_b_fir.sh`) using `env -i` | MUGQIC modules on fir inject Python 3.13 onto PYTHONPATH, breaking conda environment activation |
| `tx2gene.tsv` symlink added to `data/star_salmon/` | Part of the nf-core/rnaseq output deposited by the authors in GEO (GSE309820); not included in their GitHub repo |
| OmnipathR installed via BiocManager + `CURL_CA_BUNDLE` set to fir CA bundle | OmnipathR not included in the conda-packed 653_renv from ubuntudesk; SSL cert path differs on RHEL9 |
| `library(ggplot2)` added to `dicer1_scrna_add_infer_cnv_seurat_patched.R` | Missing in authors' script; works interactively (loaded by Seurat) but not in batch |
| ggplot2 3.5.1 force-reinstalled in 653_renv | conda-pack from ubuntudesk shipped ggplot2 4.0.3 files under a 3.5.1 label; S7 class incompatibility caused failures in `keep.scale` argument handling |
| DimHeatmap tryCatch blocks: `png()` device → `combine=TRUE` + `ggsave()` | `png()` device calls can fail in SLURM batch mode; `combine=TRUE` returns a patchwork object that ggsave handles cleanly |
| DoHeatmap/DotPlot cell_type blocks wrapped in tryCatch in main xenium script | Data structure issues (duplicate gene names, NA cell labels) cause crashes before `saveRDS`; tryCatch allows sentinel outputs to be written; recovery handled by `xenium_mouse_recover` |
| `skip=1` added to `read_excel(source_data.xlsx, sheet="Fig. S1k,l")` | Sheet has a title row above the column header; without skip, column names come out as `...1`, `...2`, etc. |

---

## Output Directory Structure

All outputs are written to `workdir/dicer1-cell-of-origin/output/`. Total: **293 files**.

### `output/scrna/` — 7 files, 1.4 GB
All `.rds` files (Seurat/SCE objects):

| File | Size | Contents |
|---|---|---|
| `dicer1_scrna_sce_qc.rds` | 346 MB | QC-filtered SingleCellExperiment, all samples |
| `dicer1_scrna_unfilter_control.rds` | 245 MB | Control (WT) Seurat, pre-cluster-filter |
| `dicer1_scrna_filter_control.rds` | 36 MB | Control Seurat, filtered |
| `dicer1_scrna_unfilter_tumor.rds` | 527 MB | Tumor Seurat, pre-cluster-filter |
| `dicer1_scrna_filter_tumor.rds` | 82 MB | Tumor Seurat, filtered + annotated |
| `dicer1_scrna_filter_control_and_tumor.rds` | 162 MB | Integrated joint object |
| `dicer1_scrna_trajectory_cds.rds` | ~15 MB | Monocle3 CellDataSet |

### `output/tables/` — 45 files, 95 MB
All `.csv` files. Key outputs:

| File | Contents |
|---|---|
| `dicer1_rna_dmem_de_mut_vs_wt.csv` | DESeq2 results: DMEM condition |
| `dicer1_rna_mesencult_de_mut_vs_wt.csv` | DESeq2 results: MesenCult condition |
| `dicer1_rna_norm_counts_no_log.csv` | Normalized count matrix |
| `dicer1_scrna_filter_tumor_pathway_activity.csv` | PROGENy activity scores per cell |
| `dicer1_swgs_sample_metadata.csv` | sWGS sample annotations (35 samples) |
| `dicer1_swgs_chr1_fisher_test.csv` | Fisher test: Chr1 gain SARC-like vs LGMT-like |
| `dicer1_xenium_mut_markers_cell_types.csv` | FindAllMarkers results, MUT by cell type |
| Multiple Xenium annotation CSVs | Per-sample (TMA1, TMA2, K1–K4) cell-type/cluster annotations |

### `output/plots/` — 174 files, 584 MB
92 PNG + 82 SVG files covering all visualization outputs.

### `output/infer_cnv/` — 53 files, 5.0 GB
Contains the complete inferCNV run outputs:

| Type | Count | Contents |
|---|---|---|
| `.rds` | 2 | `run_final_infercnv_obj.rds`, inferCNV metadata |
| `infercnv_obj` | 17 | Intermediate inferCNV objects at different HMM steps |
| `.txt` | 18 | Observation matrices, gene lists, cell annotations |
| `.png` | 5 | Intermediate heatmap visualizations |
| `.svg` | 2 | Final inferCNV visualization (paper figure) |
| `.dat` / `.mcmc_obj` | 6 | MCMC chain files |
| `.tsv` | 1 | Cell annotation file |

### `output/xenium/` — 2 files, 7.7 GB
| File | Size | Contents |
|---|---|---|
| `dicer1_xenium_wt.rds` | 3.9 GB | Fully processed WT Xenium Seurat v5 object |
| `dicer1_xenium_mut.rds` | 3.8 GB | Fully processed MUT Xenium Seurat v5 object |

---

## Limitations and Approximations

### sWGS copy number analysis
**What was done:** The authors' CN analysis requires
`dicer1_swgs_50kb_noXY_copyNumbersSegmented_filtered_called.rds`, a QDNAseq
output object from preprocessing ~70 shallow WGS samples (PRJNA1221971, SRA).
**This file was not deposited in GEO.**

As an approximation, the workflow uses two supplementary Excel files deposited
with the paper:
1. `source_data.xlsx` (sheet "Fig. S1k,l"): collapsed segment-level CN data,
   used to reproduce the genome-wide CN frequency plots (SARC-like and LGMT-like).
2. `data03_41467_2026_70971_MOESM5_ESM.xlsx` (Supplementary Data 3): per-sample
   Chr1/Chr6/Chr2/Chr15 CN status calls, used for the Fisher's exact test.

The approximation fully reproduces the paper's key quantitative result.
Restricting the Fisher test to samples with a definitive Chr1 call (excluding "n.a." entries
from Supplementary Data 3) gives LGMT-like 5/9 gain vs SARC-like 23/25 gain,
p = 0.031, OR = 8.41 — identical to the paper.


### Xenium human analysis
The paper also includes a Xenium analysis of **human** DICER1 tumor samples.
This analysis requires a processed Seurat RDS object from the human Xenium data
that **was not deposited in GEO**. `rule xenium_human` is therefore excluded from
`rule all`. The mouse Xenium analysis is fully reproduced.

---

## Workflow Directories (full map)

```
clean_snakemake_run/
├── Snakefile                          # 15-rule workflow (343 lines)
├── config/
│   └── config.yaml                    # Paths, sample lists, env pointers
├── run_dicer1_snakemake.sbatch        # SLURM submission script (mem=192G, 8 CPUs, 2d time limit)
├── scripts/
│   ├── rscript_a_fir.sh               # Wrapper: 653_renv R 4.2.2 + env -i isolation
│   ├── rscript_b_fir.sh               # Wrapper: dicer1_r453 R 4.5.3 + env -i isolation
│   └── patched/                       # 16 patched R scripts (originals unmodified in code/)
├── logs/                              # Per-job log directories (one per sbatch run)
│   └── 49517968_2026-07-18_15-43-09/  # Final successful run logs
└── workdir/
    └── dicer1-cell-of-origin/
        ├── code/                      # Author's original R scripts (unmodified)
        ├── data/                      # Symlinks to GEO-downloaded data
        │   ├── PX*/                   # 8 scRNA-seq sample dirs (barcodes/features/matrix)
        │   └── star_salmon/           # Bulk RNA-seq Salmon quant outputs (12 samples)
        ├── metadata/                  # 3 sample annotation CSV files
        ├── reference/                 # 5 reference files (gene order, marker YAMLs)
        ├── supplementary_data -> ...  # Symlink to 653_dicer1_sWGS/supplementary_data/
        └── output/
            ├── scrna/                 # 7 RDS objects (1.4 GB)
            ├── tables/                # 45 CSV files (95 MB)
            ├── plots/                 # 174 PNG/SVG files (584 MB)
            ├── infer_cnv/             # 53 files (5.0 GB): inferCNV objects, matrices, figures
            └── xenium/                # 2 RDS files (7.7 GB): wt + mut Seurat objects
```

Total output: **293 files, ~14 GB**
