# Bulk RNA-seq: SRP043694

This workflow processes the 14 samples in `metadata/sample_sheet.tsv` (7 LGD and 7 SIM). It retrieves paired FASTQs when needed, trims reads, aligns them, counts exon-overlapping fragments, and performs DESeq2 and enrichment analyses.

## Inputs and reference

FASTQs are stored in `raw_fastq/` as `SAMPLE_1.fastq.gz` and `SAMPLE_2.fastq.gz`. Their ENA paths are in the sample sheet. `scripts/00_download_fastq.sh` downloads missing files and checks gzip integrity; valid local copies are skipped.

Build or obtain a GRCh38 HISAT2 index and a matching Ensembl release 110 GTF (the release recorded in the original featureCounts command). Set `HISAT2_INDEX` to the index prefix, `REF_GTF` to the annotation path, and optionally `THREADS` (default 2).

## Environment

Use Linux or WSL with Conda/Mamba:

```bash
conda env create -f environment.yml
conda activate rnaseq-srp043694
```

## Pipeline

Run from the repository root. First copy the example configuration, set the local reference paths in it, and source it in the same shell:

```bash
cp scripts/config.example.sh scripts/config.local.sh
# Edit scripts/config.local.sh with the local reference paths.
source scripts/config.local.sh
bash scripts/00_download_fastq.sh
bash scripts/01_fastp.sh
bash scripts/02_hisat2_sort.sh
bash scripts/03_featurecounts_matrix.sh
Rscript scripts/run_deseq2.R
Rscript scripts/run_enrichment.R
Rscript scripts/plot_deseq2.R
```

The featureCounts step requires both mates to map, excludes chimeric fragments, and counts read pairs. Its matrix excludes summary rows beginning with `__`. Enrichment separates genes higher and lower in SIM and uses all tested genes as background. Confirm the sample sheet and study design before interpreting the `~ condition` contrast as independent biological samples.

The tracked `results/deseq2/summary.txt` reports 23,551 genes after filtering and 136 significant genes at `padj < 0.05` and `|log2FC| >= 1` (52 higher and 84 lower in SIM). These are historical outputs and have not been regenerated during this audit.

Code is licensed under the MIT License; see [LICENSE](LICENSE).

## Configuration and reruns

`scripts/config.local.sh` is a machine-specific local file and is ignored by Git. `PROJECT_ROOT` selects the project directory; the older `PROJECT` variable is also accepted. Run the R steps from the repository root.

Trimming publishes outputs only after fastp succeeds and gzip validation passes. A completion marker is required to reuse trimming outputs. Existing outputs without a marker must be inspected and moved aside or processed into a separate output directory. Alignment writes a temporary BAM and checks it before publishing. Reuse checks detect incomplete files, but do not establish that the reference, inputs, or parameters are unchanged; use fresh output directories when changing them.

`STRANDNESS` accepts 0, 1, or 2. Count filenames include that value. For a value other than 0, export `COUNTS_FILE` to the resulting matrix before DESeq2. Fragment counting uses `-p --countReadPairs`; the historical command lacked `--countReadPairs`, so old outputs do not demonstrate the corrected counting workflow. Plotting uses the fitted `results/deseq2/dds.rds` produced by the same analysis.

The study design and sample independence still require confirmation from the study metadata. Rerunning does not validate those assumptions.

Author: Narges Shayesteh.
