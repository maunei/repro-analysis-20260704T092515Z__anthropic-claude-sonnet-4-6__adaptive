# Pending runs — 20260704T092515Z run

Scripts must be run from:
`workflow_draft/repos/dicer1-cell-of-origin/code/`

R binary:
`/home/node/MAURICIO/papers/653_Spatial_single_cell_transcriptomic_analysis_informs_tumor_developmental_hierarchy_of_DICER1_syndrome_related_sarcoma_2026/653_renv/conda_env/bin/Rscript`

---

## Block A — inferCNV (Fig 4e + Supp Fig 9 c/d heatmaps)

### Status
- Step 1 `dicer1_scrna_export_infer_cnv.R` — ✅ DONE
- Step 2 `dicer1_scrna_run_infer_cnv_v2_patched.R` — ⏳ stopped at checkpoint step 14 (resume_mode=TRUE, will pick up automatically)
- Step 3 `dicer1_scrna_process_infer_cnv.R` — ✅ WRITTEN, not yet run (new script, was missing from pipeline)
- Step 4 `dicer1_scrna_add_infer_cnv_seurat.R` — ⏳ not yet run

### What step 3 does
`infercnv::add_to_seurat()` reads the HMM prediction files from `output/infer_cnv/combined_filter/`
and adds per-cell chromosome duplication proportions (`proportion_dupli_chrN`) to the merged
Seurat object, saving it as `dicer1_scrna_filter_combined_infercnv_metadata.rds`.
This file is what step 4 reads — it did not exist before and would have caused a crash.

### Memory + execution notes

**Why the previous run was killed:** It appeared stuck at step 15 (`define_signif_tumor_subclusters`).
This step is the longest (can take hours with no visible output) and uses Leiden algorithm internally —
which has a known hang bug in Seurat v5. The patched script uses `random_trees` instead, which fixes
that. So the previous kill was likely premature; the step may complete fine now.

**This is NOT an OOM risk** (unlike Xenium). InferCNV's HMM operates on sparse matrices and fits
comfortably in 15 GB. The concern is just runtime and potential hangs.

**zram active (as of 2026-07-05):** `/dev/zram0` — 4 GB zstd compressed swap at priority 100.
Effectively adds ~8–12 GB fast swap before the system touches the slower swapfile. This provides
a comfortable buffer for any transient memory spikes during inferCNV steps.

### How to run (background, survives session disconnect)

All paths below are relative to the **run folder root**:
`/home/node/MAURICIO/papers/653_.../REPRODUCIBILITY_ANALYSIS/20260704T092515Z__anthropic-claude-sonnet-4-6__adaptive/`

```bash
RUN=/home/node/MAURICIO/papers/653_Spatial_single_cell_transcriptomic_analysis_informs_tumor_developmental_hierarchy_of_DICER1_syndrome_related_sarcoma_2026/REPRODUCIBILITY_ANALYSIS/20260704T092515Z__anthropic-claude-sonnet-4-6__adaptive
RSCRIPT=/home/node/MAURICIO/papers/653_Spatial_single_cell_transcriptomic_analysis_informs_tumor_developmental_hierarchy_of_DICER1_syndrome_related_sarcoma_2026/653_renv/conda_env/bin/Rscript
CODE=$RUN/workflow_draft/repos/dicer1-cell-of-origin/code
OUT=$RUN/workflow_draft/repos/dicer1-cell-of-origin/output/infer_cnv

# Step 2 — resume from step 14 checkpoint (3–8 hours; runs in background)
cd $CODE
nohup $RSCRIPT dicer1_scrna_run_infer_cnv_v2_patched.R > $OUT/run_v2_patched.log 2>&1 &
echo $! > $OUT/step2.pid
echo "Step 2 PID: $(cat $OUT/step2.pid)"

# Steps 3 and 4 — run ONLY after step 2 finishes (~10 min each)
# Check step 2 is done first (see monitoring section below), then:
cd $CODE
$RSCRIPT dicer1_scrna_process_infer_cnv.R > $OUT/process_infer_cnv.log 2>&1
$RSCRIPT dicer1_scrna_add_infer_cnv_seurat.R > $OUT/add_seurat.log 2>&1
```

### How to monitor step 2

```bash
OUT=/home/node/MAURICIO/papers/653_Spatial_single_cell_transcriptomic_analysis_informs_tumor_developmental_hierarchy_of_DICER1_syndrome_related_sarcoma_2026/REPRODUCIBILITY_ANALYSIS/20260704T092515Z__anthropic-claude-sonnet-4-6__adaptive/workflow_draft/repos/dicer1-cell-of-origin/output/infer_cnv

# Is it still running?
ps aux | grep "dicer1_scrna_run_infer_cnv"

# What is it doing right now? (last 20 lines of log)
tail -20 $OUT/run_v2_patched.log

# What checkpoints have been written? (each file = one completed inferCNV step)
ls -lht $OUT/combined_filter/*.infercnv_obj 2>/dev/null | head -5

# RAM usage right now
free -h
```

**Step 2 is done when** the log ends with:
`"Wrote run_final_infercnv_obj.rds"` or similar, AND
`$OUT/combined_filter/run_final_infercnv_obj.rds` exists.

### How to kill step 2 if needed

```bash
OUT=/home/node/MAURICIO/papers/653_Spatial_single_cell_transcriptomic_analysis_informs_tumor_developmental_hierarchy_of_DICER1_syndrome_related_sarcoma_2026/REPRODUCIBILITY_ANALYSIS/20260704T092515Z__anthropic-claude-sonnet-4-6__adaptive/workflow_draft/repos/dicer1-cell-of-origin/output/infer_cnv

# Option 1 — kill using saved PID (preferred)
kill $(cat $OUT/step2.pid)

# Option 2 — kill by name if PID file is gone
pkill -f "dicer1_scrna_run_infer_cnv"

# Verify it's dead
ps aux | grep "dicer1_scrna_run_infer_cnv"
```

**After killing:** `resume_mode=TRUE` is set. Re-running step 2 with the same nohup command above
will automatically pick up from the last saved checkpoint — no data is lost.

### How to check where it will resume from

```bash
OUT=/home/node/MAURICIO/papers/653_Spatial_single_cell_transcriptomic_analysis_informs_tumor_developmental_hierarchy_of_DICER1_syndrome_related_sarcoma_2026/REPRODUCIBILITY_ANALYSIS/20260704T092515Z__anthropic-claude-sonnet-4-6__adaptive/workflow_draft/repos/dicer1-cell-of-origin/output/infer_cnv

ls -lht $OUT/combined_filter/*.infercnv_obj 2>/dev/null | head -3
# The highest-numbered file is the last checkpoint. Step 15 is the next to run.
```

### Output plots produced by step 4
- `output/infer_cnv/dicer1_scrna_filter_tumor_infercnv.svg` → **Fig 4e**
- `output/infer_cnv/dicer1_scrna_filter_tumor_infercnv_split_sample.svg` → Fig 4e (split by sample)

### Note on Supp Fig 9 c/d heatmaps
InferCNV itself generates large heatmap PNGs during the run (step 2). These appear in
`output/infer_cnv/combined_filter/` as `infercnv.png` / `infercnv.preliminary.png` once
the run completes. These are the Supp Fig 9 c/d reference figures.

---

## Block B — Xenium mouse (Fig 2e–h, Fig 3e–h, Supp Fig 8)

### Status — ❌ HARDWARE-BLOCKED (requires ≥32 GB RAM)

**Cannot run on this machine.** The agent attempted this 3 times during the initial run:
- v1: FOV class error (Seurat v4→v5 migration issue)
- v2: SIGKILL during subset (OOM)
- v3: two-pass memory management — also SIGKILL (OOM)

Root cause: After QC filtering, the WT control subset has **622,209 cells**. SCTransform on 622K cells requires ~30–50 GB RAM for the dense residual matrix computation. This machine has **15 GB total**. No software workaround exists for this step — BPCells cannot reduce peak RAM during SCTransform.

- Input RDS exists: `output/xenium/dicer1_xenium_mus_musculus_all.rds` (2.2 GB) ✅
- Script: `workflow_draft/scripts/patched/dicer1_xenium_proseg_mus_musculus_patched.R` (v3) ✅ ready
- **To run**: needs HPC or cloud instance with ≥32 GB RAM; no code or data gap

### Commands (for future high-RAM system)

```bash
RSCRIPT=...  # 653_renv2/conda_env/bin/Rscript (R 4.5.3 with SeuratWrappers/Banksy)
CODE=workflow_draft/repos/dicer1-cell-of-origin/code
LOG=workflow_draft/repos/dicer1-cell-of-origin/output

$RSCRIPT $CODE/dicer1_xenium_proseg_mus_musculus_patched.R \
  > $LOG/xenium_mus_musculus.log 2>&1
```

### Output plots expected (once unblocked)
- `output/plots/dicer1_xenium_mut_umap_cell_type.svg` → Fig 3e–h / Supp Fig 8
- `output/plots/dicer1_xenium_mut_dotplot_cell_type.svg`
- `output/plots/dicer1_xenium_mut_top3_genes_clusters.svg`
- `output/plots/dicer1_xenium_mut_subset_feature_plot.svg` → Fig 3e–h markers
- `output/plots/dicer1_xenium_mut_subset_vln_plot.svg`
- (and corresponding PNGs)

---

## Permanently blocked (no data)

- **Fig 5 + Supp Fig 10** — human Xenium RDS not deposited on GEO. Contact authors.

---

## After Block A completes — update the feasibility report

Once step 4 finishes and `output/infer_cnv/dicer1_scrna_filter_tumor_infercnv.svg` exists:

1. Update `reports/01_replication_feasibility_assessment.md`:
   - Mark inferCNV row as ✅ complete in execution results table
   - Add Fig 4e SVG to Generated Outputs section
   - Add HMM heatmap PNGs (`combined_filter/infercnv.png`) as Supp Fig 9 c/d outputs
