# DICER1 Sarcoma Reproducibility Assessment — v8

**Paper:** Kommoss FKF, Zhang JYH, Lynch BJ, Chen SY, Senz J, Moscovitz Y, et al.  
"Spatial single-cell transcriptomic analysis informs tumor developmental hierarchy of DICER1 syndrome-related sarcoma."  
*Nature Communications* 2026. PMID: 41912509 · PMC: PMC13199431 · DOI: 10.1038/s41467-026-70971-6

**Author repository:** github.com/Huntsmanlab/dicer1-cell-of-origin  
**Assessment maintained at:** `…/claude-science/dicer1_docs/`  
**Analysis root:** `…/claude-science/`  
**Last updated:** 2026-08-28

---

## Part 1: Understanding the paper

### 1.1 DICER1 syndrome and the biological question

DICER1 syndrome (DICER1-related tumor predisposition, DRTP) is a rare hereditary cancer syndrome caused by germline pathogenic variants in the *DICER1* gene. DICER1 encodes a ribonuclease (RNase III) central to microRNA biogenesis: it cleaves pre-microRNAs into mature ~22 nt miRNAs that silence target messenger RNAs post-transcriptionally. Affected individuals carry one defective germline copy; tumor development requires a somatic "second hit," almost exclusively a hotspot missense mutation in the RNase IIIb domain. The resulting compound heterozygous state causes a global collapse of the 5p-strand miRNA pool.

Clinically, DICER1 syndrome produces pleiotropic tumors — pleuropulmonary blastoma (PPB), cystic nephroma, embryonal rhabdomyosarcoma, Sertoli-Leydig cell tumors, and others — all occurring predominantly in children and young adults and all sharing similar molecular features.

The paper asks: What is the cell of origin of DICER1-associated sarcomas, and **what developmental trajectory do these tumors follow?** Identifying the cell of origin is critical for understanding why sarcomas driven by the same molecular defect arise at diverse anatomical sites and for identifying developmental states and signaling pathways that might be therapeutically targeted. Using a tamoxifen-inducible, lineage-traceable mouse model that generates a compound-heterozygous *Dicer1* state in Hic1+ mesenchymal stromal cells—one RNase IIIb-mutant D1693N allele and one disrupted allele—the authors show that this alteration produces renal tumors closely resembling the histological and molecular continuum of human DICER1-associated mesenchymal neoplasms, from low-grade lesions to overt sarcoma. Their single-cell and spatial analyses identify Hic1+ Pdgfra+ Dpt+ Pi16+ universal fibroblast-like progenitors—particularly cells associated with the subepithelial basement membrane of transitional epithelium—as the putative cell of origin in the mouse model. DICER1 disruption appears to perturb the differentiation of these progenitors and generate divergent trajectories: one branch undergoes aberrant, stepwise rhabdomyoblastic differentiation through transitional myogenic states, whereas another progresses through a fibroblast-like “Ground” state toward highly proliferative sarcoma, accompanied by copy-number alterations, reduced p53-associated activity, and increased MAPK signaling. Analysis of 16 human DICER1-associated tumors from several anatomical sites identifies analogous progenitor, myogenic, Ground, and proliferative states with similar spatial relationships, suggesting that this developmental hierarchy is conserved in human disease. The distribution of universal fibroblasts across subepithelial, capsular, and perivascular niches also offers a possible explanation for why related tumors can arise in both epithelial-lined and non-epithelial organs. However, the paper establishes the fibroblastic origin directly only in its renal mouse model; in human tumors and at other anatomical sites, the proposed common origin remains an inference from shared transcriptional and spatial cell states rather than definitive lineage tracing.

### 1.2 The mouse model — a triple transgenic system

To answer this question, the authors created a triple-transgenic genetically engineered mouse model (GEMM). Understanding each component is necessary to interpret what the data means:

**Component 1 — Hic1CreERT2 (cell specificity):** HIC1 (Hypermethylated in Cancer 1) is a transcription factor that marks quiescent mesenchymal stromal/progenitor populations in many vascularized adult tissues, including fibroblastic progenitors, perivascular stromal cells, and pericytes. The CreERT2 fusion protein (Cre recombinase fused to a mutated estrogen receptor ligand-binding domain) remains cytoplasmic and inactive until tamoxifen is administered. Tamoxifen binds the mutated ER domain, causing nuclear translocation and activation. Because CreERT2 is driven by the Hic1 promoter, only Hic1-expressing cells respond to tamoxifen. This means the entire model is restricted to the Hic1+ mesenchymal compartment — no other cell type is affected.



**Component 2 — Dicer1^fl/D1693N (the mutation):** This is a carefully engineered conditional knock-in allele. Before tamoxifen, the D1693N-bearing exon is inverted in the genome and not positioned for normal transcription — the D1693N mutation is therefore silent in all cells until Cre acts. A second *Dicer1* allele carries a loxP-flanked exon (the floxed allele). When tamoxifen activates CreERT2 in Hic1+ cells, Cre-mediated recombination simultaneously (a) flips the D1693N-bearing exon into the correct orientation so it is now expressed, and (b) deletes the floxed exon in the other allele, creating a loss-of-function. After tamoxifen, Hic1+ cells therefore become compound heterozygous — one allele expressing the D1693N hotspot mutation, the other disrupted — exactly the molecular state of human DICER1 sarcomas.

**Component 3 — Rosa26^LSL-tdTomato (lineage tracing):** The Rosa26 locus is ubiquitously expressed. A loxP-Stop-loxP cassette prevents tdTomato expression until Cre is activated. After tamoxifen, Hic1+ cells permanently express tdTomato — and all their descendants do too, regardless of whether they still express Hic1. This is the lineage tracer: any cell derived from a Hic1+ precursor is tdTomato+ for life.

**How it works together:** Three weeks after birth, mice are injected with tamoxifen. In Hic1+ cells, Cre recombination simultaneously activates D1693N expression, disables the floxed Dicer1 allele, and switches on tdTomato. Over the following weeks, these Hic1+ compound-heterozygous cells can be tracked by FACS (tdTomato+) and, in kidneys where sarcomas arise, their progression from normal MSC to tumor cell can be studied by single-cell and spatial transcriptomics.

### 1.3 The bulk RNA-seq cell lines

For the bulk RNA-seq experiment, the authors cultured primary mesenchymal cells derived from the GEMM model. Tumor and control kidney tissues from DICER1-mutant (HDT) and control mice were dissociated, and tdTomato+ cells were enriched by FACS, then established as primary cultures in two media conditions — MesenCult and DMEM/F12 + 10% FBS. These are therefore not immortalized cell lines in the classical sense, but primary cultures derived directly from the same GEMM mice used for all other experiments. Additional cell lines were created by CRISPR-mediated *Trp53* knockout using lentiCRISPRv2GFP and gRNAs targeting *Trp53* exons 6 and 10, modelling *Trp53*-driven progression seen in aggressive DICER1 sarcomas. Both mutant and control primary cultures appear in the accession GSE309820.

### 1.4 What the experiments aim to show

The paper is organized around a central claim: DICER1 sarcoma originates from a subset of Hic1+ mesenchymal progenitors resembling Dpt+ Pi16+ universal fibroblasts, whose differentiation is disrupted by DICER1 mutation and subsequently diverges toward either aberrant rhabdomyoblastic differentiation or, through a fibroblast-like “Ground” state, proliferative sarcoma. The evidence comes from four complementary deposited murine datasets, plus human spatial transcriptomic data whose raw files were **not publicly deposited** but are available from the authors upon request.

| Data type | GEO/SRA accession | Platform | Key purpose | Deposited |
|-----------|------------------|-----------|-------------|-----------|
| scRNA-seq | GSE288990 | 10x Genomics Chromium | Cell type atlas of Hic1+ tumor and control cells | ✅ |
| Xenium spatial (mouse) | GSE289001 | 10x Genomics Xenium In Situ (379 genes) | Spatial organisation of cell states within tumor tissue | ✅ |
| Bulk RNA-seq | GSE309820 | DNBSEQ-G400RS (Complete Genomics) | Transcriptional reprogramming of mutant vs control MSC cultures | ✅ |
| sWGS | PRJNA1221971 | Illumina short-read | Somatic copy-number landscape across 35 tumors | ✅ |
| Xenium spatial (human TMA) | not deposited | 10x Genomics Xenium In Situ (Human Multi-Tissue and Cancer Panel PN1000626, 379 genes) | Validation of mouse cell-type markers in human DICER1 sarcoma tissue | ❌ |

The human Xenium TMA experiment covered a fourth Xenium slide containing a section of a human tumor tissue microarray. It was used to confirm that the cell-type markers identified in the mouse model are also detectable in human DICER1 sarcoma tissue. Because this data was not deposited on GEO, we did not attempt to reproduce this part of the analysis.

In addition, a targeted sequencing panel for variant calling (KRAS/NRAS/HRAS/TP53/DICER1 hotspots) and a qPCR/FISH assay for copy number validation at two loci are described in the paper. These are addressed in the independent analyses section (§7.1).

---

## Part 2: What the authors deposited for reproducibility

### 2.1 Repository and scripts

The GitHub repository (github.com/Huntsmanlab/dicer1-cell-of-origin) contains the R analysis scripts used in the study, together with a small set of supporting metadata/reference files. The repository does not contain the primary analysis data or generated analysis outputs; its data/ directory is essentially empty, and the output/ subdirectories contain only .gitkeep placeholders.

The key scripts present in the repository at the time of our analysis were:

**Scripts from the author repository:**
- `dicer1_rna.R` — bulk RNA-seq import with tximport + DESeq2 differential expression + downstream enrichment analysis
- `dicer1_scrna_controls.R` — scRNA-seq processing + integration + clustering + annotation (controls)
- `dicer1_scrna_controls_and_tumors.R` — joint integration + comparison (controls + tumors)
- `dicer1_scrna_export_infer_cnv.R` — InferCNV input file generation/export
- `dicer1_scrna_qc.R` — initial scRNA-seq QC + filtering + normalization + exploratory visualization
- `dicer1_scrna_add_infer_cnv_seurat.R` — addition of InferCNV-derived results back into Seurat metadata for visualization
- `dicer1_scrna_run_infer_cnv.R` — InferCNV invocation
- `dicer1_scrna_run_pathway_analysis.R` — PROGENy pathway scoring via decoupleR (`get_progeny` + `run_mlm`)
- `dicer1_scrna_trajectory.R` — Monocle3 trajectory inference
- `dicer1_scrna_tumors.R` — tumor scRNA-seq processing + integration + clustering + cell type/state annotation
- `dicer1_swgs.R` — sWGS QDNAseq copy-number pipeline
- `dicer1_xenium_proseg_mus_musculus.R` — Xenium mouse spatial analysis with ProSeg
- `dicer1_xenium_proseg_homo_sapiens.R` — Xenium human spatial analysis with ProSeg

**Scripts created entirely by us (not in the GitHub repository):**
- `dicer1_swgs_from_rds.R` — re-run of the sWGS analysis starting from our own QDNAseq-segmented RDS, necessary after an OST filesystem failure on the HPC (see §5.1 rule `swgs`)
- `dicer1_xenium_banksy.R` — Banksy spatial domain analysis (not in the author repo; implemented by us based on the paper's description)
- `dicer1_xenium_recover_celltype_plots.R` — a recovery script to regenerate cell-type-level dotplots from a previously saved `mut.rds` after a job failure
- `dicer1_scrna_export_infer_cnv_filter_only.R` — a variant of the InferCNV export for the filtered-only pipeline branch
- `dicer1_scrna_process_infer_cnv.R` — InferCNV post-processing script created by us; not present in the author GitHub repository
- `dicer1_scrna_add_to_seurat_intermediate.R` — an intermediate step script for InferCNV annotation insertion into the Seurat object
- `dotplot_compat_patch.R` — a shared helper that patches DotPlot color scale for Seurat 5.x compatibility


### 2.2 R environment

The repository does not include a `renv.lock` file or any description of the R environment. The R environments used for this analysis were built by us, initially on a local workstation (ubuntudesk, Ubuntu 20.04) and then reconstructed on the HPC. This distinction matters: any package version difference between the author's original environment and ours is our deviation, not theirs, because we have no record of exactly which versions the authors used.

Two conda-backed R environments were used on the HPC:
- **`653_renv` / `dicer1_r453`** — the primary environment, containing Seurat 5.x, Monocle3, InferCNV, DESeq2, tximport, decoupleR, and all other packages needed by the author scripts.
- One package was added during the analysis: **OmnipathR 3.5.25**, required at runtime by `decoupleR::get_progeny()` for PROGENy pathway scoring. This package is not listed in the authors' GitHub repository. It was installed separately because the R environment was initially built on the local workstation — which had internet access — and OmnipathR was installed interactively at the time it was found to be missing from the initial build. On the HPC it was installed on the login node via `rscript_a_fir.sh`.

### 2.3 Deposited data — per accession

**GSE288990 (scRNA-seq):** Raw FASTQ files via SRA (consent=public); processed Cell Ranger matrix trios (barcodes, features, matrix) from the authors' own Cell Ranger run deposited on GEO; two intermediate RDS files:
- `dicer1_scrna_filter_tumor.rds` — 2,483 Hic1+ tumor cells with cell type labels (authors' final analysis object, 103 MB)
- `dicer1_scrna_unfilter_tumor.rds` — 15,346 unfiltered tumor cells (528 MB)
These intermediate files are critical — see §3.1 for why.

**GSE309820 (bulk RNA-seq):** Raw FASTQ files via SRA (24 SRR runs, 12 samples × 2 runs each); per-sample `quant.sf.gz` files (Salmon output) deposited on GEO; the transcript-to-gene mapping `GSE309820_tx2gene.tsv.gz`.

**GSE289001 (murine Xenium spatial):** Authors' deposited processed Seurat RDS (`dicer1_xenium_mus_musculus_all.rds`, ~2.2 GB), together with deposited Xenium Ranger-derived files for the mouse samples, including transcript-level spatial data such as `tma_2_transcripts.zarr.zip`; cell × gene count-matrix components `tma_2_matrix.mtx.gz`, `tma_2_barcodes.tsv.gz`, and `tma_2_features.tsv.gz`; and morphology images `tma_2_morphology.ome.tif` and `tma_2_morphology_focus.ome.tif`. The deposition also includes `dicer1_xenium_mus_musculus_slide_map.xlsx`. Note: the corresponding Xenium Ranger-derived files for the murine K3 sample are absent from GEO (§6.2).



**PRJNA1221971 (sWGS):** Raw FASTQ files via SRA; Excel supplementary data files containing copy-number segment results from the authors' QDNAseq analysis.

**Panel variant calling:** FASTQs via SRA; no processed VCF files deposited.

**What was NOT deposited:** The authors' run-specific ProSeg-derived cell metadata, transcript-to-cell assignments, and segmentation geometry were not deposited as intermediates; only the final processed Seurat RDS contains the `cell_area` field used by the analysis (§6.3). Depending on invocation, ProSeg can write tabular cell metadata as CSV or Parquet, while polygon geometry is normally represented separately (for example as GeoJSON or in a spatial-data Zarr), so these missing outputs should not be described collectively as “polygon `.parquet` files.” The exact tdTomato sequence appended to the Cell Ranger reference used in the scRNA-seq pipeline was also not deposited (§3.1). For the Xenium Ranger transcript Zarr archives, our attempts did not recover a comprehensive mapping from the compact UUID fields to the accompanying matrix barcodes (§6.4). Consequently, the deposited pre-ProSeg Xenium files and the authors' ProSeg-derived RDS are not connected by a demonstrated, end-to-end reconstruction in this project.

---

## Part 3: The central challenge — the Hic1+ cell identification problem

### 3.1 The tdTomato story

The scRNA-seq pipeline depends critically on identifying which cells in the sequencing data are Hic1+ (tdTomato+). In the mouse experiment, these cells are sorted out by FACS before sequencing — but the sequencing also records tdTomato transcript counts as an additional feature in the Cell Ranger output, allowing the bioinformatic pipeline to confirm the isolation by requiring ≥5% tdTomato+ cells per sample before proceeding.

When we ran Cell Ranger on the public SRA FASTQ data, tdTomato counts were nearly absent: only **263 / 22,897 cells (1.1%)** had any tdTomato counts. This is far below the authors' ≥5% threshold. The pipeline consequently falls through to a Level 4 fallback (using all mesenchymal cells regardless of tdTomato expression), which does not replicate the biologically specific Hic1+ selection the authors performed.

The explanation for the 1.1% detection rate: the Cell Ranger reference genome used to generate the SRA data included a custom tdTomato sequence appended to the standard GRCm38 reference. The exact sequence of that custom tdTomato construct — which specific fluorescent protein variant, which codon optimisation, which cloning vector backbone — was not deposited with the data. We constructed our Cell Ranger reference using a candidate sequence (KJ081243.1 from GenBank), which is apparently not the exact cassette used by the authors. As a result, reads that should align to tdTomato fail to do so, producing the near-zero detection rate. A full description of our Cell Ranger pipeline and reference construction is given in `02_nfcore_scrnaseq_assessment.md` in the nextflow analysis directory.

This is a genuine reproducibility gap. Without the exact tdTomato FASTA sequence, we cannot replicate the Hic1+ isolation step from raw data alone.

**Note:** Our initial tdTomato reference choice (`KJ081243.1`) was likely inappropriate; a more faithful reconstruction should instead explore the publicly documented Ai9 reporter sequence, including the tdTomato/WPRE cassette. We did not rerun Cell Ranger with that alternative at this stage; instead, we proceeded by matching the Hic1+/tdTomato+ cell IDs available in the authors' intermediate files as described below:


### 3.2 How we addressed it — the barcode bridge

The authors deposited their processed Hic1+ tumor object (`dicer1_scrna_filter_tumor.rds`, 103 MB) on GEO. This object contains the 2,483 cells that passed their tdTomato filter, labelled with their seven published cell types.

Our approach was to use this intermediate file as a barcode bridge: extract the cell barcodes from the authors' filtered object, compare them with our independently produced Seurat object (`dicer1_scrna_seurat_ours.rds`, 22,860 cells from our Cell Ranger → Seurat pipeline; §7.4), and use the overlap to define an Hic1+-enriched query subset. Because the objects use incompatible sample prefixes and separators, the implemented match discarded those prefixes and retained only the 16-character ACGT barcode. On the author side, 2,396 of 2,483 filtered cells had a 16-mer present in our object (96.5%). Subsetting our object by membership in that 16-mer set selected 2,410 full cell identifiers.

The two totals measure different sides of a non-unique lookup. The 2,410 query cells contain 2,395 unique 16-mers: 15 barcode sequences occur in two different query samples, creating 15 additional query-cell occurrences. These cross-sample collisions are ambiguous because sample identity was not included in the key. Likewise, the 87 author cells without a matching 16-mer were absent from our final object; differences in cell calling and QC are plausible causes, but the diagnostic does not establish the cause for each cell.

The authors' RDS supplied cell identities and published cell-type labels, not expression values, to the independent bridge analysis. The resulting `dicer1_hic1_bridge.rds` contains our RNA data, SCTransform result, PCA, Harmony integration, clustering, and UMAP, together with an author-derived `authors_cell_type` metadata field. Most assignments are unambiguous, but assignments associated with the 15 duplicated query 16-mers should not be described as exact one-to-one matches. The Phase 2 scRNA analysis therefore uses independently processed expression data but author-assisted cell selection and comparison labels.

The full barcode bridge analysis and its results are described in §7.6 and in the dedicated report `ANALYSIS_REPORT.md` at `653_dicer1_nextflow/scRNA/barcode_bridge/`.

### 3.3 The Xenium cell_area trick — and its limits

A parallel challenge arose for the Xenium spatial analysis. ProSeg generated revised transcript-to-cell assignments and cell-segmentation information; the resulting processed object contains `cell_area`, which the authors' R script subsequently filters at `cell_area > 20 µm²`. The run-specific ProSeg intermediates were not deposited, so this field was absent from our matrix-derived reconstruction.

We tested a spatial nearest-neighbour bridge: for every reconstructed cell with valid coordinates, find the nearest cell in the authors' RDS within 5 µm and transfer its `cell_area` and selected metadata. TMA1 had explicit `cell_id`, `x_location`, and `y_location` fields in the deposited transcript CSV, allowing coordinates to be reconstructed for essentially all cells. The bridge found an author cell within 5 µm for 279,041 of 376,264 reconstructed TMA1 cells (74.2%; mean distance 1.37 µm). The algorithm did not enforce one-to-one matching, so this number is a count of reconstructed query cells with a nearby author cell—not a count of distinct author cells recovered.

Across TMA2, K1, K2, and K4, only 512 additional query cells were assigned `cell_area` in the final v3 attempt because the compact UUID fields in the deposited transcript Zarr archives could not be decoded into comprehensive matrix-cell centroids (§6.4). The v3 bridge output contains 279,553 matches across all five sections, of which 279,041 are TMA1. Importantly, the v3 workflow did **not** apply `cell_area > 20`: its file named `dicer1_xenium_bridged_area_filtered_v3.rds` was filtered only by `nCount_Xenium > 5` and `nFeature_Xenium > 5`, and all 1,504,937 cells passed. The bridge is therefore a partial metadata-transfer experiment, not a reproduction of the authors' area-filtered cell population.

---

## Part 4: Overall project structure and phases

Our reproducibility project proceeded in two broad phases.

**Phase 1** was an exploratory period in which we ran individual analyses on a local workstation (ubuntudesk) and on the HPC, developing code, identifying bugs, and building up understanding of the data. Phase 1 used the authors' patched scripts together with some of the authors' own intermediate files as starting points, and included an early Snakemake workflow started on the local workstation and then adapted for the HPC. This was a development and exploration phase, not a clean reproducibility demonstration.

**Phase 2** is the main subject of this report. It is a completed HPC workflow combining independent reconstruction with explicitly identified uses of deposited processed objects. The scRNA-seq, bulk RNA-seq, and sWGS branches largely begin from public FASTQs, although the scRNA Hic1+ subset is defined through the author-assisted barcode bridge (§3.2). The principal Snakemake Xenium branch starts from the authors' deposited ProSeg-derived Seurat RDS. Separately, we reconstructed and analyzed a five-section Xenium object from deposited matrices; that independent branch omits K3 and the `cell_area > 20` filter (§6). Thus, computational completion should not be equated with uniform end-to-end reproduction from raw data.

---

## Part 5: Phase 2 — completed Snakemake workflow

The Phase 2 Snakemake workflow lives at:  
`…/claude-science/653_dicer1_snakemake_complete/`

Orchestrator sbatch: `sbatch/11_run_snakemake.sbatch`  
Snakefile: `Snakefile` (project root)  
All scripts: `workdir/dicer1-cell-of-origin/code/` (author originals) and `workdir/dicer1-cell-of-origin/code/scripts/patched/` (our modified versions)  
All outputs: `workdir/dicer1-cell-of-origin/output/`  
SLURM logs: `logs/smk_slurm/`

### 5.1 Rule-by-rule description

For each rule, **inputs are distinguished by origin**: files produced by our pipeline versus files from the authors' GEO deposit. Dependencies on authors' processed objects are identified for each rule rather than reduced to a count of “exceptions.” When a rule's output is consumed by a later rule, that dependency is noted. Our reconstructed outputs reproduce the authors' file naming conventions, so the same filename may refer to either our version or the authors' deposited version; context makes clear which is which.

---

**Rule: `scrna_qc`**

*Input (ours):* Cell Ranger output matrix trios (barcodes, features, matrix) produced by our own nf-core/scrnaseq Cell Ranger run on the 7 tumor and 4 control SRA FASTQ files from GSE288990. The authors also deposited their own Cell Ranger matrix trios on GEO, but we used our independently generated matrices. Our Cell Ranger run is documented in `02_nfcore_scrnaseq_assessment.md`.  
*Script:* `dicer1_scrna_qc.R` (author script, lightly modified for path references).  
*What it does:* Reads each sample into SingleCellExperiment (SCE) objects; computes per-cell QC metrics; applies adaptive outlier-based filtering (scran/scuttle); creates a merged SCE object.  
*Authors' deposited equivalent for this step:* None — the authors deposited their Cell Ranger matrix trios but not a merged QC'd SCE object.  
*Output (ours, feeds into `scrna_controls`, `scrna_tumors_unfiltered`, `scrna_run_infercnv`):* `dicer1_scrna_sce_qc.rds` (346 MB).  
*Status:* ✅ COMPLETE

---

**Rule: `scrna_controls`**

*Input (ours):* `dicer1_scrna_sce_qc.rds` from rule `scrna_qc` above; control sample IDs.  
*Script:* `dicer1_scrna_controls_patched.R` (author's `dicer1_scrna_controls.R` with 8 patches; see §5.2.1).  
*What it does:* Subsets control samples; runs Seurat normalization, PCA, and integration across samples (FindIntegrationAnchors + IntegrateData); clusters integrated object; identifies marker genes.  
*Authors' deposited equivalent:* `dicer1_scrna_filter_control.rds` (34 MB) and `dicer1_scrna_unfilter_control.rds` (245 MB) on GEO. We reproduce these under the same filenames, but they are our independently computed outputs.  
*Output (ours, feeds into `scrna_controls_and_tumors`):* `dicer1_scrna_filter_control.rds` (31 MB) + `dicer1_scrna_unfilter_control.rds` (245 MB).  
*Status:* ✅ COMPLETE (SLURM job 56926133, ~3 hours)

---

**Rule: `scrna_tumors_unfiltered`**

*Input (mixed):* `dicer1_scrna_sce_qc.rds` from rule `scrna_qc` (ours); tumor sample IDs; plus the **approved exception**: authors' `dicer1_scrna_filter_tumor.rds` from GEO, used as the barcode source for 16-mer Hic1+ identification (§3.2).  
*Script:* `dicer1_scrna_tumors_patched.R` (author's `dicer1_scrna_tumors.R` with 4 patches + bridge patch; see §5.2.2).  
*What it does:* Subsets tumor samples; applies the 16-mer barcode bridge to identify Hic1+ cells — replacing the tdTomato filter that fails on SRA data; runs Seurat normalization and PCA on the full unfiltered tumor object.  
*Authors' deposited equivalent:* `dicer1_scrna_unfilter_tumor.rds` (528 MB) on GEO.  
*Output (ours, feeds into `scrna_tumors_filtered`):* `dicer1_scrna_unfilter_tumor.rds` (535 MB).  
*Status:* ✅ COMPLETE (SLURM job 56926134, ~1 h 27 min)

---

**Rule: `scrna_tumors_filtered`**

*Input (mixed):* `dicer1_scrna_unfilter_tumor.rds` from rule `scrna_tumors_unfiltered` (ours); plus the **approved exception**: authors' `dicer1_scrna_filter_tumor.rds` from GEO, used as barcode source for Hic1+ subsetting (§3.2). Expression values and all analyses use only our data.  
*Script:* `dicer1_scrna_tumors_filter_only_patched.R` (our script implementing the filter-only branch; see §5.2.3).  
*What it does:* Applies the 16-mer bridge to subset to Hic1+ cells; runs SCTransform v2, PCA, integration; clusters at multiple resolutions; assigns cluster labels; runs Monocle3 trajectory inference.  
*Authors' deposited equivalent:* `dicer1_scrna_filter_tumor.rds` (103 MB) on GEO. We produce a file of the same name, but it is our independently processed object.  
*Output (ours, feeds into `scrna_controls_and_tumors`, `scrna_trajectory`, `scrna_add_infercnv_seurat`):* `dicer1_scrna_filter_tumor.rds` (103 MB) + `dicer1_scrna_trajectory_cds.rds` (20 MB) + 56 plots.  
*Status:* ✅ COMPLETE (SLURM job 56953186, ~25 min)

---

**Rule: `scrna_controls_and_tumors`**

*Input (ours):* `dicer1_scrna_filter_control.rds` from rule `scrna_controls` + `dicer1_scrna_filter_tumor.rds` from rule `scrna_tumors_filtered` — both our reconstructed versions, not the GEO files.  
*Script:* `dicer1_scrna_controls_and_tumors_patched.R` (author's script with 5 patches; see §5.2.4).  
*What it does:* Merges control and tumor filtered objects; runs joint integration across all samples; produces the final combined Seurat object.  
*Authors' deposited equivalent:* None. The authors did not deposit a `dicer1_scrna_filter_control_and_tumor.rds` on GEO.  
*Output (ours, feeds into `scrna_pathway_analysis`):* `dicer1_scrna_filter_control_and_tumor.rds` (181 MB).  
*Status:* ✅ COMPLETE (SLURM job 57009208 child, completed 22:08 PDT Aug 26 2026)

---

**Rule: `scrna_trajectory`**

*Input (ours):* `dicer1_scrna_filter_tumor.rds` + `dicer1_scrna_trajectory_cds.rds`, both from rule `scrna_tumors_filtered`.  
*Script:* `dicer1_scrna_trajectory_patched.R` (author's script with minor path patches).  
*What it does:* Runs Monocle3 trajectory analysis on the Hic1+ tumor cells; infers pseudotime and developmental ordering; produces trajectory UMAP plots.  
*Authors' deposited equivalent:* `dicer1_scrna_trajectory_cds.rds` (20 MB) on GEO. We produce it under the same name.  
*Output (ours):* Trajectory UMAP plots.  
*Status:* ✅ COMPLETE

---

**Rule: `scrna_pathway_analysis`**

*Input (ours):* `dicer1_scrna_filter_control_and_tumor.rds` from rule `scrna_controls_and_tumors`.  
*Script:* `dicer1_scrna_run_pathway_analysis_patched.R` (author's script extended for OmnipathR).  
*What it does:* Computes PROGENy pathway activity scores per cell using decoupleR. The authors' script calls `get_progeny()`, which requires OmnipathR (not listed in the authors' GitHub) to fetch the PROGENy gene-weight matrix from the OmniPath database at runtime. OmnipathR 3.5.25 was installed to satisfy this dependency (see §2.2). PROGENy scores 14 cancer-relevant signalling pathways (EGFR, WNT, MAPK, JAK-STAT, TGFb, p53, and others) per cell.  
*Authors' deposited equivalent:* None deposited.  
*Output (ours):* PROGENy pathway score plots.  
*Status:* ✅ COMPLETE (child of SLURM orchestrator 56962679, ~4 min)

---

**Rule: `scrna_run_infercnv`**

*Input (ours):* `dicer1_scrna_sce_qc.rds` from rule `scrna_qc`; our cell annotation file; gene position file.  
*Script:* `dicer1_scrna_run_infer_cnv_v2_patched.R` (with dynamic reference group detection).  
*What it does:* Runs InferCNV to infer somatic copy-number variations from scRNA-seq read depth patterns; uses control cell clusters as reference. The patch replaced a hardcoded list of 12 reference group labels with auto-detection from the annotation file (2 of the 12 hardcoded labels — "Control Trans Epi", "Control Mural" — did not exist in our 15-group annotation, causing failure).  
*Authors' deposited equivalent:* InferCNV outputs deposited on GEO. We did not use those.  
*Output (ours, feeds into `scrna_export_infercnv`):* InferCNV HMM states, CNV score matrices.  
*Status:* ✅ COMPLETE (child of SLURM orchestrator 56962679, ~3 hours)

---

**Rule: `scrna_export_infercnv`**

*Input (ours):* InferCNV output matrices from rule `scrna_run_infercnv`; `dicer1_scrna_unfilter_tumor.rds` from rule `scrna_tumors_unfiltered`.  
*Script:* `dicer1_scrna_export_infer_cnv_patched.R`.  
*What it does:* Exports InferCNV CNV scores per cell back into the Seurat object metadata.  
*Output (ours, feeds into `scrna_add_to_seurat_intermediate`):* CNV-annotated metadata tables.  
*Status:* ✅ COMPLETE

---

**Rule: `scrna_add_to_seurat_intermediate`**

*Input (ours):* CNV metadata from rule `scrna_export_infercnv`.  
*Script:* `dicer1_scrna_add_to_seurat_intermediate.R` (our script; no author original for this step).  
*What it does:* Inserts InferCNV annotations into an intermediate Seurat object, bridging between the InferCNV export step and the final annotation step.  
*Output (ours, feeds into `scrna_add_infercnv_seurat`):* InferCNV-annotated intermediate Seurat object.  
*Status:* ✅ COMPLETE

---

**Rule: `scrna_add_infercnv_seurat`**

*Input (ours):* Outputs from `scrna_add_to_seurat_intermediate` + `dicer1_scrna_filter_tumor.rds` from rule `scrna_tumors_filtered`.  
*Script:* `dicer1_scrna_add_infer_cnv_seurat_patched.R`.  
*What it does:* Adds InferCNV-derived CNV labels to the filtered tumor Seurat object, producing the final CNV-annotated version used for figures.  
*Output (ours):* InferCNV-annotated Seurat object.  
*Status:* ✅ COMPLETE

---

**Rule: `bulk_rna`**

*Input (ours):* Salmon `quant.sf` files produced by our own nf-core/rnaseq run on the 24 SRR FASTQ files from GSE309820 (see §7.3); `GSE309820_tx2gene.tsv.gz` downloaded directly from GEO.  
*Script:* `dicer1_rna_patched.R` (author's `dicer1_rna.R` with minor path corrections).  
*What it does:* Imports per-sample Salmon quantification via tximport; runs DESeq2 differential expression separately for DMEM and MesenCult conditions (mutant HDT1056 vs control HDT1119, 3 replicates per group); generates volcano plots, heatmaps, GSEA results.  
*Authors' deposited equivalent:* Per-sample `quant.sf.gz` files on GEO. We used our independently generated files. The high correlation (r ≥ 0.999) between the two confirms they are effectively equivalent (§7.3).  
*Output (ours, reproduces authors' output filenames):* DESeq2 result tables, volcano plots, GSEA plots.  
*Status:* ✅ COMPLETE

---

**Rule: `swgs`**

*Input (ours):* Aligned BAMs from our own BWA-MEM alignment of downloaded SRA FASTQ files (PRJNA1221971); supplementary Excel files for metadata and comparison.  
*Script:* `dicer1_swgs_from_rds.R` (our script, based on author's `dicer1_swgs.R`). The script runs from our QDNAseq-segmented RDS (`dicer1_swgs_50kb_noXY_copyNumbersSegmented_filtered_called.rds`) — a file the authors used internally but **never deposited**. We reconstructed this intermediate by running QDNAseq on our own aligned BAMs. The `_from_rds` variant was necessary after an OST filesystem failure on the HPC destroyed earlier intermediate files.  
*What it does:* Runs QDNAseq copy-number analysis; produces CN frequency plots showing the copy-number landscape across 35 HDT tumor and 7 non-neoplastic kidney samples; applies Fisher's exact test to assess enrichment of recurrent chromosomal gains — notably Chr1q — in tumor vs. non-neoplastic samples (result: p=0.031, OR=8.41, consistent with the paper's reported enrichment); generates the metadata panel figure.  
*Authors' deposited equivalent:* Excel supplementary files contain the final copy-number segments. The authors' intermediate QDNAseq RDS was not deposited.  
*Output (ours):* CN frequency plots, Fisher CSV (p=0.031, OR=8.41), metadata panel (PDF).  
*Status:* ✅ COMPLETE

---

**Rule: `xenium_mouse`**

*Input (authors' processed object):* Authors' deposited Xenium RDS (`dicer1_xenium_mus_musculus_all.rds`, 2.2 GB, from GEO GSE289001). The main Snakemake Xenium rule is a downstream reanalysis of this complete ProSeg-derived expression object, not merely an identity lookup. Although the deposited pre-ProSeg transcript coordinates appear sufficient in principle to perform a new ProSeg run, we did not perform that reconstruction, and the authors' exact run-specific segmentation cannot be recovered directly from the deposited intermediates. Separately, we reconstructed and analyzed a five-section Xenium object from deposited matrices (§6.3 and §7.2); that independent branch omitted K3 and the `cell_area` filter.  
*Script:* `dicer1_xenium_proseg_mus_musculus_patched.R` (author's script with DimHeatmap and DotPlot patches; see §5.2.5).  
*What it does:* Normalizes and clusters the full Xenium object; splits into wild-type (normal kidney) and mutant (tumor) subsets; runs UMAP and marker analysis for each; saves `wt.rds` and `mut.rds`.  
*Output (ours, feeds into `xenium_mouse_recover` and `xenium_banksy`):* `dicer1_xenium_wt.rds` (3.9 GB) + `dicer1_xenium_mut.rds` (3.8 GB) + WT/MUT plots.  
*Status:* ✅ COMPLETE

---

**Rule: `xenium_mouse_recover`**

*Input (ours):* `dicer1_xenium_mut.rds` from rule `xenium_mouse`; marker CSVs from the same rule.  
*Script:* `dicer1_xenium_recover_celltype_plots.R` (our script).  
*What it does:* Regenerates cell-type-level dotplots and heatmaps that crashed during `xenium_mouse` due to a Xenium subset with fewer than 3 expressed panel genes. Loading the saved `mut.rds` allows these plots to be produced in a short job without re-running the expensive clustering.  
*Output (ours):* Cell-type dotplot, top-3-gene-per-cell-type plots.  
*Status:* ✅ COMPLETE

---

**Rule: `xenium_banksy`**

*Input (ours):* `dicer1_xenium_wt.rds` + `dicer1_xenium_mut.rds` from rule `xenium_mouse`.  
*Script:* `dicer1_xenium_banksy.R` (our script; Banksy spatial domain analysis not in the author repo).  
*What it does:* Applies Banksy (version 1.9.1 from prabhakarlab/Banksy@devel — the Bioconductor 3.22 version 1.6.0 is incompatible with our R 4.5 environment, and version 1.8.1 requires R 4.6) to define spatial domains within the Xenium sections. A patch was required: `computeHarmonics()` returns complex-valued matrices (using `j=sqrt(as.complex(-1))`) even when M=0, triggering a `zgeMatrix` dispatch removed in Matrix 1.7+. The fix wraps the result with `Re()` and `as.matrix()` before matrix multiplication.  
*Output (ours):* `dicer1_xenium_banksy_wt.rds` (5.6 GB) + `dicer1_xenium_banksy_mut.rds` (5.4 GB) + spatial domain plots.  
*Status:* ✅ COMPLETE

---

### 5.2 Patch summaries by script

#### 5.2.1 `dicer1_scrna_controls_patched.R` — 8 patches

All patches address the same root cause: small samples (some with < 50 cells after QC) cause default Seurat integration parameters to request more neighbours than cells exist. Patches cap PCA dimensions, k.filter, k.score, and k.weight dynamically based on per-sample cell counts.

| Patch | Fix |
|-------|-----|
| 7 | `npcs = min(50, max_pcs)` — cap PCA dims per sample |
| 9 | `dims = 1:min_pcs_sample` in FindIntegrationAnchors |
| 10 | Physically truncate PCA embeddings to `min_pcs_filter` columns before integration |
| 11 | `k.filter = min(200, min_cells_filter - 1)` |
| 12 | `k.filter = NA` when `min_cells_filter < 50` |
| 13 | `k.score = min(30, min_cells_filter - 1)` |
| 14 | `k.weight = min(50, n_anchors - 1)` dynamically |
| 15 | Additional `k.weight = min(50, min_cells_filter - 1)` guard |

#### 5.2.2 `dicer1_scrna_tumors_patched.R` — 4 patches + bridge

| Patch | Fix |
|-------|-----|
| 8 | Fallback: `mes_cells & !discard_cells` when tdTomato zero count |
| 10 | `GetAssayData(slot=)` deprecated → `tryCatch(layer=, slot=)` |
| 11 | Empty mes_cells → 4-level cascading fallback |
| 12 | 26 clusters, cluster_labels_filter maps only 0-9 → map extras to "Unk" + `safe_subset_downsample` |
| Bridge (14) | Replace entire tdTomato/fallback block with 16-mer barcode bridge (§3.2) |

#### 5.2.3 `dicer1_scrna_tumors_filter_only_patched.R` — 3 patches

| Patch | Fix |
|-------|-----|
| 16 | Replace `Tomato > 0` filter with 16-mer bridge |
| 16b | Remove `mes_cells &` from subset intersection |
| 16c | Fix 16-mer regex: `.*[._]([ACGT]{16})-[0-9]+$` — achieved 2,491/15,327 matched (16.3%) |

#### 5.2.4 `dicer1_scrna_controls_and_tumors_patched.R` — 5 patches (same root cause as controls)

| Patch | Fix |
|-------|-----|
| 1 | Per-sample `npcs_s = min(50, ncell-1, nfeature-1)` in RunPCA lapply |
| 2 | `min_pcs_s` → `dims_s` for FindIntegrationAnchors |
| 3 | Truncate all PCA embeddings to `min_pcs_s` cols before FindIntegrationAnchors (fixes rbind column mismatch) |
| 4 | `k.filter = NA` when min_cells < 50; else `min(200, ncell-1)` |
| 5 | `k.weight = min(50, ncell-1)` in IntegrateData |

#### 5.2.5 `dicer1_xenium_proseg_mus_musculus_patched.R` — 2 patches

| Patch | Fix |
|-------|-----|
| DimHeatmap | Wrapped `DimHeatmap()` output with `ggsave()` (Seurat 5.x returns a ggplot object; `png()`/`dev.off()` cannot capture it) |
| DotPlot subsets | `tryCatch()` wrappers around DotPlot calls that fail when a cell-type subset has fewer than 3 expressed panel genes |

---

## Part 6: Xenium spatial analysis — the segmentation gap

### 6.1 Data structure and what the authors deposited

The Xenium dataset (GSE289001) covers six tissue sections: four whole kidney sections (K1, K2, K3, K4) and two TMAs (TMA1 = 32 tumor cores, TMA2 = 36 tumor cores). The authors deposited a pre-processed 2.2 GB Seurat RDS as well as selected Xenium Ranger-derived matrix, transcript-coordinate, and morphology files. A fourth slide with human TMA data was run but not deposited (§1.4).

### 6.2 K3 — completely absent from GEO

K3 (whole kidney section 3) raw data was not deposited on GEO. Its transcripts are absent from the combined GEO files, and no per-section K3 files exist. The authors' RDS contains 113,730 K3 cells, which were processed internally from a copy never uploaded. K3 cannot be reconstructed from public data — this is a genuine reproducibility limitation.

### 6.3 Our independent reconstruction — the ProSeg gap

Starting from the deposited pre-ProSeg Xenium Ranger matrix and transcript-coordinate files, we independently reconstructed the multi-section Xenium Seurat object using `LoadXenium()` and metadata from Supplementary Data 5. The reconstruction (job 54544591, 2 h 04 min) produced 1,504,937 cells across five sections (K1, K2, K4, TMA1, TMA2). The downstream analysis was subsequently completed on this object in job 54691173, with the `cell_area > 20` condition removed and K3 absent. This independent branch is distinct from the principal Snakemake Xenium rule, which used the authors' RDS.

The main deviation from the authors' object arises from the ProSeg re-segmentation step. GEO states that the deposited Xenium Ranger-derived files reflect the original Xenium Ranger 1.7.1.1 segmentation using a 5-µm nuclear expansion, whereas the deposited Seurat RDS contains data after re-segmentation with ProSeg. The authors' analysis script subsequently filters this ProSeg-derived object using `nCount_Xenium > 5`, `nFeature_Xenium > 5`, and `cell_area > 20`.

The intermediate ProSeg outputs used by the authors—including run-specific per-cell metadata, transcript assignments, and segmentation geometry—were not deposited in GEO. The deposited pre-ProSeg transcript coordinates appear sufficient in principle to perform a new ProSeg segmentation. Such a run would be an independent re-segmentation, however, rather than recovery of the authors' exact result; ProSeg is sampling-based, and the complete invocation and run-specific outputs are unavailable.

We did not rerun ProSeg from the deposited Xenium data. Consequently, our reconstruction does not reproduce the authors' ProSeg-derived segmentation or `cell_area` metadata and cannot apply `cell_area > 20` in the same way. This contributes to differences in the retained cell population, particularly near TMA core boundaries. The limitation is therefore an unperformed independent re-segmentation combined with absence of the authors' exact intermediates—not proof that public inputs are intrinsically insufficient to run ProSeg.


The cell count comparison (GEO raw vs Authors' RDS vs Our reconstruction):

| Section | GEO raw barcodes | Authors' RDS | Our reconstruction | Delta (ours vs authors) |
|---------|------------------|--------------|--------------------|------------------------|
| K1 | 184,925 | 186,042 | 184,782 | −0.7% ✅ |
| K2 | 334,792 | 331,230 | 333,779 | +0.8% ✅ |
| K3 | absent | 113,730 | absent | — |
| K4 | 180,576 | 162,445 | 175,936 | +8.3% ⚠️ |
| TMA1 | absent as labeled file; data in combined files | 281,638 | 376,264 | +33.6% ⚠️ |
| TMA2 | 436,397 | 419,931 | 434,176 | +3.4% ✅ |

The TMA1 discrepancy (+33.6%, approximately 94,600 extra cells) is the most visible difference. Small or boundary-associated cells are a plausible contributor, and the authors' downstream `cell_area > 20` criterion would remove some of them, but the current bridge did not apply that threshold or establish that it accounts for the entire difference (§3.3). For K1/K2/K4/TMA2, see §6.4.

### 6.4 The zarr bridge problem — TMA2/K1/K2/K4

TMA1 could be bridged because its deposited transcript CSV contains explicit `cell_id`, `x_location`, and `y_location` fields, allowing centroids to be reconstructed for essentially every matrix cell. The other four public sections use transcript Zarr archives with compact `uuid_w0` and `uuid_w1` fields. Two decoding attempts were made. The v2 attempt incorrectly treated `uuid_w0` as a cell index and produced misleading broad matches. The v3 attempt used `uuid_w1 - 1`, but recovered only 295 unique identifiers for TMA2, 74 for K1, 133 for K2, and 73 for K4 despite hundreds of thousands of matrix barcodes per section. The diagnostic pattern—few unique `uuid_w1` values across tens of millions of transcripts—shows that this field does not provide the required comprehensive matrix-cell index under that interpretation.

Accordingly, the immediate demonstrated limitation is that our decoding attempts did not connect the Zarr transcript assignments comprehensively to the matrix cell identifiers; it is too strong to state simply that the Zarr fields are codewords or that spatial centroids are impossible in principle. In the final v3 nearest-neighbour bridge, only 73 K1, 122 K2, 62 K4, and 255 TMA2 query cells were matched within 5 µm, compared with 279,041 TMA1 query cells. These 512 non-TMA1 transfers are too sparse to support section-wide area filtering.

The independent five-section analysis therefore ran without `cell_area > 20`, using the available transcript-count criteria instead. K3 remains absent from public data entirely. The main Snakemake Xenium analysis, by contrast, starts from the authors' processed RDS and therefore includes the authors' ProSeg-derived data and K3.

---

## Part 7: Independent analyses

These analyses were conducted outside the Snakemake workflow, as separate investigations to assess specific reproducibility claims in the paper.

### 7.1 Analysis 1 — Panel variant calling (KRAS concordance)

The paper reports targeted sequencing of a cancer panel across 35 DICER1 sarcoma samples to identify somatic variants. We downloaded the panel FASTQs and ran our own variant calling pipeline (BWA-MEM → GATK → SnpEff annotation), then compared our KRAS calls against the Kras mutation status column in the paper's supplementary metadata.

**Result:** 35/35 concordant. Every KRAS-mutant sample in the metadata was called as mutant in our pipeline, and every wild-type sample as wild-type. This is strong evidence that the variant calling is straightforward and reproducible from the raw data.

**Assessment directory:** `…/claude-science/653_dicer1_sWGS/`

### 7.2 Analysis 2 — Xenium cell bridge (TMA1 cell_area transfer)

Described in §3.3 and §6.3–6.4. In TMA1, 279,041 of 376,264 reconstructed query cells (74.2%) had an author-RDS cell within 5 µm. Because the algorithm independently assigns the nearest author cell to each query cell and does not enforce one-to-one matching, this is not a count of distinct author cells recovered. The bridge is a partial, many-to-one metadata transfer; it did not apply `cell_area > 20`. Coverage for TMA2/K1/K2/K4 was negligible because comprehensive matrix-cell centroids could not be decoded from the deposited Zarr structures.

**Assessment directory:** `…/claude-science/653_dicer1_nextflow/Xenium/`  
**Detail report:** `04_xenium_assessment.md` in that directory.

### 7.3 Analysis 3 — Bulk RNA-seq quantification comparison

We ran the complete nf-core/rnaseq pipeline (v3.21.0, matching the authors' v3.14.0 as closely as possible) on the 24 raw SRR FASTQ files from GSE309820. The pipeline used STAR alignment + Salmon quantification with `--seqBias --gcBias` flags as specified in the paper methods. The resulting `quant.sf` files were compared to the authors' deposited `quant.sf.gz` files on GEO.

**Transcript-level comparison:** Our Salmon-estimated read counts per transcript correlate with the deposited GEO counts at r ≥ 0.999 for all 12 samples — confirming that nf-core/rnaseq v3.21.0 reproduces the authors' Salmon output.

**DESeq2 log2 fold change correlation:**  
After running the author's `dicer1_rna.R` with our reproduced `quant.sf` files:  
- DMEM condition: Pearson r = 0.9794  
- MesenCult condition: Pearson r = 0.9747

These high correlations confirm that the differential expression signal underlying the paper's biological conclusions is reproduced. The minor version difference (3.14.0 vs 3.21.0) has no material impact.

**Assessment directory:** `…/claude-science/653_dicer1_nextflow/RNASeq/`  
**Detail report:** `01_nfcore_rnaseq_bulk_assessment.md` in that directory.

### 7.4 Analysis 4 — Independent scRNA-seq Seurat reconstruction

We ran Cell Ranger on the raw scRNA-seq FASTQs from GSE288990 and built an independent Seurat object from scratch using SCTransform v2 normalization, PCA, and Harmony batch correction (by sample). This produced `dicer1_scrna_seurat_ours.rds` (371 MB, 22,860 cells after QC). Clustering at multiple resolutions identified 23 clusters. FindAllMarkers yielded 42,265 differentially expressed genes. This object is the source of cells used in the barcode bridge (§3.2 and §7.6).

**Assessment directory:** `…/claude-science/653_dicer1_nextflow/scRNA/seurat_ours/`

### 7.5 Analysis 5 — Banksy spatial domain analysis

Banksy (Bioconductor) was applied to the authors' deposited Xenium RDS to define spatial domains. All 15 mutant cell types described in the paper were recovered. A technical patch was required (see §5.1 rule `xenium_banksy` for details).

**Assessment directory:** `…/claude-science/653_banksy/`

### 7.6 Analysis 6 — Hic1+ barcode bridge clustering (ARI = 0.338)

This was the primary test of whether independent expression processing and clustering reproduce the authors' seven-cell-type classification after author-assisted selection of the Hic1+-enriched query subset. The barcode bridge approach and its collision limitation are described in §3.2.

The full pipeline selected 2,410 query cells through the 16-mer bridge, then ran SCTransform v2 → PCA (30 dimensions) → Harmony batch correction by sample → FindNeighbors/FindClusters/UMAP at resolution 0.6. The cluster assignments were compared with author-derived cell-type labels using the adjusted Rand index (ARI). Expression values and clustering were ours; cell selection and comparison labels depended on the authors' filtered object. The subset contains 2,395 unique 16-mers, including 15 sequences duplicated across query samples.

**ARI = 0.338.** An ARI of 1.0 indicates perfect agreement and 0 indicates agreement expected by chance. A score of 0.338 means our clustering captures the same broad populations but does not reproduce the authors' exact partition. The confusion matrix shows over-partitioning of Ground and Progenitor states and splitting of TR Diff Myo across clusters; all seven author label categories are represented. Likely contributors include incomplete and non-one-to-one barcode recovery, cell-calling and QC differences, feature selection, normalization, Harmony rather than the authors' integration strategy, clustering parameters, and software versions. FACS enrichment occurred experimentally before sequencing: failure to quantify tdTomato computationally prevented direct reconstruction of the authors' filtering rule but did not alter the original FACS-enriched population.

**Output files in** `full_analysis/`:  
- `dicer1_hic1_bridge.rds` (135 MB)  
- `umap_our_clusters.pdf`, `umap_authors_labels.pdf`, `umap_sidebyside.pdf`  
- `confusion_matrix.pdf` / `confusion_matrix.csv`  
- `cluster_markers.csv` (14,799 marker genes across 10 clusters)  
- `bridge_summary.txt`

**Assessment directory:** `…/claude-science/653_dicer1_nextflow/scRNA/barcode_bridge/`  
**Full report:** `ANALYSIS_REPORT.md` in that directory.

---

## Part 8: MapQuery investigation — raising ARI to 0.699

Analysis 6 raised a natural question: how close could we get to the authors' result if we used a more direct label-transfer approach? This investigation was run as a separate, exploratory analysis outside the Snakemake workflow — it uses the authors' intermediate files as input rather than our independently processed data, so it is a methodological investigation, not a true reproducibility measure. Detail reports and reasoning chain are documented in the `mapquery_investigation/` subdirectory.

### 8.1 What is MapQuery?

MapQuery is a Seurat function (since Seurat v4; confirmed present in our environment with Seurat 5.5.1) that performs label transfer from a reference object to a query object via `FindTransferAnchors` (CCA-based anchor identification) and `TransferData` (label propagation). An optional `ProjectUMAP` step projects the query into the reference's UMAP geometry, but requires that the reference was saved with `return.model=TRUE` in `RunUMAP`.

### 8.2 What the diagnostic revealed

We ran a diagnostic job (57059461) to inspect the authors' reference objects. Neither `dicer1_scrna_filter_tumor.rds` (2,483 cells) nor `dicer1_scrna_unfilter_tumor.rds` (15,346 cells) were saved with `return.model=TRUE`. ProjectUMAP is therefore not available. We proceeded with `FindTransferAnchors + TransferData` (label transfer only).

### 8.3 Results

Reference: authors' `dicer1_scrna_filter_tumor.rds` (2,483 Hic1+ cells with cell type labels).  
Query: our bridge-subset 2,410 cells.

| Metric | Value |
|--------|-------|
| Anchors found | 5,470 |
| **ARI (predicted vs authors' labels)** | **0.699** |
| Exact label accuracy | 2,125 / 2,410 = **88.2%** |
| Mean prediction score | 0.883 |

The improvement from 0.338 to 0.699 is substantial, but the two numbers answer different questions. The 0.338 analysis independently recalculates the expression representation and clusters after author-assisted cell selection. The 0.699 result is a reference-assisted consistency analysis: the authors' filtered object supplies the MapQuery reference and, indirectly through the barcode bridge, the comparison labels. It shows that our expression profiles are broadly compatible with the authors' classification, but it is not an independent validation that every query cell is correctly identified. The main disagreements occur among biologically adjacent Ground, Progenitor, and TR Diff Myo states.

**Detail reports** in `…/claude-science/653_dicer1_nextflow/scRNA/mapquery_investigation/`:  
- `PLAN.md` — rationale and decision logic  
- `DIAGNOSTIC_ANALYSIS.md` — full diagnostic analysis  
- `MAPQUERY_CONCLUSION.md` — conclusion report

---

## Part 9: Output inventory (Phase 2 Snakemake)

All outputs at: `…/claude-science/653_dicer1_snakemake_complete/workdir/dicer1-cell-of-origin/output/`

### Reconstructed RDS files

| Our file | Size | Authors' deposited equivalent | Author size |
|----------|------|-------------------------------|-------------|
| `dicer1_scrna_sce_qc.rds` | 346 MB | none deposited | — |
| `dicer1_scrna_filter_control.rds` | 31 MB | `dicer1_scrna_filter_control.rds` (GEO) | 34 MB |
| `dicer1_scrna_unfilter_control.rds` | 245 MB | `dicer1_scrna_unfilter_control.rds` (GEO) | 245 MB |
| `dicer1_scrna_unfilter_tumor.rds` | 535 MB | `dicer1_scrna_unfilter_tumor.rds` (GEO) | 528 MB |
| `dicer1_scrna_filter_tumor.rds` | 103 MB | `dicer1_scrna_filter_tumor.rds` (GEO) | 103 MB |
| `dicer1_scrna_trajectory_cds.rds` | 20 MB | `dicer1_scrna_trajectory_cds.rds` (GEO) | 20 MB |
| `dicer1_scrna_filter_control_and_tumor.rds` | 181 MB | none deposited | — |
| `dicer1_xenium_wt.rds` | 3.9 GB | none equivalent deposited | — |
| `dicer1_xenium_mut.rds` | 3.8 GB | none equivalent deposited | — |
| `dicer1_xenium_banksy_wt.rds` | 5.6 GB | none (analysis not in GitHub) | — |
| `dicer1_xenium_banksy_mut.rds` | 5.4 GB | none (analysis not in GitHub) | — |

### Plots and tables

178 total output files in `output/plots/`, including:
- Xenium: UMAP (clusters, cell types, sample, histotype), dotplots, DimHeatmaps, MSC violin/feature plots
- scRNA: UMAP plots, trajectory plots, PROGENy pathway heatmaps, InferCNV CNV score plots
- Bulk RNA: Volcano plots, GSEA enrichment plots, PCA
- sWGS: CN frequency plots, metadata panel, Fisher's exact test results (p=0.031, OR=8.41)

---

## Part 10: Overall assessment

### 10.1 Summary of results

| Analysis | Result | Verdict |
|----------|--------|---------|
| Panel variant calling (35 KRAS calls) | 35/35 concordant | ✅ Fully reproduced |
| Xenium bridge (TMA1) | 279,041/376,264 query cells had an author cell within 5 µm; many-to-one | ⚠️ Partial metadata transfer |
| Bulk RNA DESeq2 (log2FC Pearson r) | r=0.9794 DMEM, r=0.9747 MesenCult | ✅ Strongly reproduced |
| Banksy spatial domains | 15/15 MUT cell types recovered | ✅ Reproduced |
| Independent clustering of the author-assisted Hic1+ bridge subset | ARI = 0.338 | ⚠️ Partial |
| MapQuery label transfer (exploratory) | ARI = 0.699, accuracy 88.2% | ℹ️ Reference-assisted, non-independent |

### 10.2 Where tricks were required

1. **Hic1+ cell identification (scRNA):** The tdTomato sequence appended to the Cell Ranger reference was not deposited. Without it, our Cell Ranger run recovers only 1.1% tdTomato+ cells (vs ≥5% needed). We used the authors' deposited `dicer1_scrna_filter_tumor.rds` as a barcode and label source. Author-side recovery was 2,396/2,483 (96.5%); membership-based subsetting produced 2,410 query cells because 15 16-mers occur in two query samples. Downstream expression processing and clustering used our data, but cell selection and comparison labels were author-assisted. See §3.1–3.2.

2. **Xenium `cell_area` metadata:** The authors' run-specific ProSeg intermediates were not deposited. A 5-µm nearest-neighbour bridge transferred author-RDS metadata to 279,041 TMA1 query cells, but matching was many-to-one and did not measure unique author-cell recovery. Comprehensive centroids could not be decoded from the other sections' transcript Zarr structures. The v3 bridge did not apply `cell_area > 20`. See §3.3 and §6.3–6.4.

3. **Xenium Snakemake workflow input:** Rule `xenium_mouse` starts from the authors' complete deposited ProSeg-derived RDS and is therefore a downstream reanalysis of an author-processed object. A separate five-section reconstruction and downstream analysis were completed from deposited matrices without K3 or the area filter. A new ProSeg run appears possible in principle from deposited transcript coordinates, but was not performed, and would not recover the authors' exact run-specific result automatically.

4. **MapQuery ARI boost** (exploratory, §8): Using the authors' processed Hic1+ object as a computational reference raised agreement from ARI=0.338 to ARI=0.699. Because the reference and the barcode-derived comparison labels ultimately come from the same author object, this is a reference-assisted consistency test, not independent reproducibility.

### 10.3 Key reproducibility limitations

- **K3 Xenium section:** Raw Xenium data for the murine K3 section are absent from the GEO deposit; this section therefore cannot be reconstructed from the publicly deposited raw data.

- **tdTomato reference reconstruction:** The exact custom tdTomato representation used in the authors' Cell Ranger reference was not deposited or explicitly documented. Our initial choice (`KJ081243.1`) was likely inappropriate; the publicly documented Ai9 reporter sequence, including the tdTomato/WPRE cassette, provides a more appropriate alternative that should be explored. We did not rerun Cell Ranger with this alternative at this stage and instead recovered Hic1+/tdTomato+ cell identities from the authors' intermediate files.

- **ProSeg re-segmentation:** The authors' run-specific ProSeg cell metadata, transcript assignments, and segmentation geometry were not deposited. The pre-ProSeg transcript coordinates appear sufficient in principle for a new ProSeg run, but we did not perform one; such a run would be an independent re-segmentation, not automatic recovery of the authors' exact sampling result. The authors' `cell_area > 20 µm²` filtering was therefore not reproduced in our reconstructed branch.

- **Human Xenium TMA data:** The human Xenium raw data were not publicly deposited; this part of the analysis therefore cannot be reproduced from publicly available data alone.

- **Xenium Zarr decoding and RDS linkage:** TMA1's transcript CSV supplied explicit cell IDs and coordinates, enabling broad spatial matching. For TMA2/K1/K2/K4, neither tested interpretation of the compact Zarr UUID fields produced comprehensive matrix-cell centroids; the v3 attempt yielded only a few hundred decoded identifiers per section. The public pre-ProSeg representation was therefore not linked comprehensively to the authors' ProSeg-derived RDS in this project.

- **scRNA barcode collisions:** Matching solely on the 16-mer ignores library identity. Fifteen barcode sequences occurred in two query samples, so the 2,410-cell subset contains 2,395 unique 16-mers and includes ambiguous cross-sample assignments. A rigorous one-to-one bridge would require a verified PX-to-HDT sample crosswalk or explicit exclusion of collision-associated cells.

- **nf-core/rnaseq version difference:** Our reconstruction used nf-core/rnaseq 3.21.0 rather than the authors' 3.14.0. This version difference had a negligible effect on the relevant outputs, as confirmed by comparison of the resulting `quant.sf` files.


### 10.4 Infrastructure notes

All Phase 2 computation was performed on the HPC cluster (fir, Alliance Canada / DRAC) using SLURM. Several challenges were encountered:
- Lustre OST filesystem failures (recovered via file recovery and job resubmission)
- Python 3.13 `re._constants` ModuleNotFoundError preventing `conda activate` on compute nodes (fixed by direct `Rscript` path: `/home/maunei/miniconda3/envs/dicer1_r453/bin/Rscript`)
- Banksy Matrix incompatibility (fixed by `computeHarmonics` patch)
- Complex PCA integration failures for small samples (fixed by 8-patch series in controls and controls_and_tumors scripts)
