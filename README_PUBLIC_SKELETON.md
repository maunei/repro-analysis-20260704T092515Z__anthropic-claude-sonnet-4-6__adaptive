# Public Reproducibility Analysis Skeleton

This directory is a lightweight public version of the original reproducibility analysis.

It is provided to inform the scope, organization, and file-level structure of the reproducibility analysis workflow.

The original analysis included raw data, generated data objects, compressed matrices, large RDS objects, inferCNV intermediate objects, quantification outputs, Snakemake runtime state, and other large outputs. Those files are not included here.

Instead, excluded files are represented by small placeholder files ending in:

`.NOT_INCLUDED.txt`

Each placeholder records the original path, filename, size in bytes, and reason for exclusion.

Private notes, including `my_notes_post_analysis`, are intentionally excluded. Workflow runtime/cache folders such as `.snakemake` and notebook checkpoint folders are also intentionally excluded.

## Disclaimer

The outputs, tables, plots, logs, and placeholder files included here should not be treated as finalized research results or used as input for downstream scientific analysis. This repository is under active construction and is being used to test and document a reproducibility-analysis system.

Large data files, generated objects, and selected intermediate outputs have been replaced with `.NOT_INCLUDED.txt` placeholder files. These placeholders are intended to preserve the directory structure and file provenance, not to provide analyzable data.

Any scientific interpretation should be based only on a validated, complete, and reviewed version of the analysis.
