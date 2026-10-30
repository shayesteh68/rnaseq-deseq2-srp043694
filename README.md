# Bulk RNA-seq: SRP043694

This workflow processes the 14 samples in `metadata/sample_sheet.tsv` (7 LGD and 7 SIM). It retrieves paired FASTQs when needed, trims reads, aligns them, counts exon-overlapping fragments, and performs DESeq2 and enrichment analyses.

## Inputs and reference

FASTQs are stored in `raw_fastq/` as `SAMPLE_1.fastq.gz` and `SAMPLE_2.fastq.gz`. Their ENA paths are in the sample sheet; `metadata/fastq_manifest.tsv` records the expected byte counts and MD5 checksums from ENA. `scripts/00_download_fastq.sh` checks size, MD5, and gzip integrity before reusing or publishing a file. Existing files that fail validation are preserved; use a clean `RAW_FASTQ_DIR` to retrieve verified replacements.

Build or obtain a GRCh38 HISAT2 index and a matching Ensembl release 110 GTF (the release recorded in the original featureCounts command). Set `HISAT2_INDEX` to the index prefix, `REF_GTF` to the annotation path, and optionally `THREADS` (default 2).

## Source samples and study design

The source series, GEO GSE58963 (SRA study SRP043694), contains 21 samples: 7 high-grade dysplasia (HGD), 7 low-grade dysplasia (LGD), and 7 specialized intestinal metaplasia (SIM). This workflow analyzes the 14 LGD and SIM runs. `metadata/source_accessions.tsv` maps each run to its GEO sample and BioSample accessions and submitted sample title.

The public sample records do not provide an explicit patient-level pairing field for the analyzed runs. The current DESeq2 model uses an unpaired design (`~ condition`), which assumes independent biological samples; that assumption is not confirmed by the available metadata. Do not interpret the contrast as a within-patient comparison, and confirm patient-level independence with the source investigators before making clinical claims.

## Environment

Use Linux or WSL with Conda/Mamba:

```bash
conda env create -f environment.yml
conda activate rnaseq-srp043694
```

The exact Linux x86_64 environments used for this audit are pinned in `environment-analysis-linux-64.lock.txt` (R/DESeq2/enrichment) and `environment-tools-linux-64.lock.txt` (fastp/HISAT2/samtools/featureCounts). Recreate either environment with `conda create -n ENV_NAME --file LOCK_FILE`.

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

The featureCounts step requires both mates to map, excludes chimeric fragments, and counts read pairs. Its matrix excludes summary rows beginning with `__`. Enrichment separates genes higher and lower in SIM and uses all tested genes mapped to Entrez IDs as the background.

## Audited run and interpretation

The current audit counted the 14 existing BAM files with featureCounts using paired-fragment counting (`-p --countReadPairs -B -C -s 0`) and the Ensembl release 110 GTF. This validates counting and downstream analysis from those BAMs; it does not establish that the FASTQ trimming and alignment steps were rerun from raw reads in this audit. The generated matrix has 62,754 gene rows and 14 sample columns.

The DESeq2 run retained 20,842 genes after low-count filtering. The standard Wald test with `padj < 0.05` and `|MLE log2FC| >= 1` identified 147 candidates (57 higher and 90 lower in SIM). After `ashr` shrinkage, 44 exploratory candidates met `padj < 0.05` and `|shrunken log2FC| >= 1` (16 higher and 28 lower in SIM). This combined filter is a screening rule; it is not a formal test that the true effect exceeds one log2 unit. The separate DESeq2 threshold test (`lfcThreshold = 1`, `altHypothesis = "greaterAbs"`) found one gene at BH-adjusted `p < 0.05`. Its counts are nonzero in only two SIM samples and zero in all seven LGD samples, so it needs independent validation and should not be treated as a robust cohort-wide finding.

The sample PCA and per-sample metrics are in `results/deseq2/pca_plot.png` and `results/deseq2/sample_qc.tsv`. `SRR1502800` is a strong PCA outlier with a squamous/keratin expression program; it was retained rather than removed automatically. The sample set is small, patient independence is unresolved, and the outlier may reflect sample or tissue heterogeneity. Treat gene and pathway findings as exploratory.

Enrichment uses the tested, Entrez-mapped genes as background. In the current run, the 28 lower-in-SIM candidates produced 227 GO biological-process terms at the configured cutoffs; the 16 higher-in-SIM candidates produced none, and neither direction produced a significant KEGG term. About 15.7% of tested Ensembl IDs did not map to Entrez IDs. KEGG results depend on the live KEGG service and may change as its database changes. The tables and plots are under `results/enrichment/` and `results/plots/`.

The current summary, candidate tables, PCA, volcano plot, sample QC table, and enrichment outputs are regenerated by the scripts above. `results/deseq2/summary.txt` records the thresholds and counts; `results/deseq2/top10_genes.txt` labels its ranked lists as exploratory and reports the formal threshold-test result separately.

Code is licensed under the MIT License; see [LICENSE](LICENSE).

## Configuration and reruns

`scripts/config.local.sh` is a machine-specific local file and is ignored by Git. `PROJECT_ROOT` selects the project directory; the older `PROJECT` variable is also accepted. Run the R steps from the repository root.

Trimming publishes outputs only after fastp succeeds and gzip validation passes. A completion marker is required to reuse trimming outputs. Existing outputs without a marker must be inspected and moved aside or processed into a separate output directory. Alignment writes a temporary BAM and checks it before publishing. Reuse checks detect incomplete files, but do not establish that the reference, inputs, or parameters are unchanged; use fresh output directories when changing them.

`STRANDNESS` accepts 0, 1, or 2. Count filenames include that value. For a value other than 0, export `COUNTS_FILE` to the resulting matrix before DESeq2. Fragment counting uses `-p --countReadPairs`; the historical command lacked `--countReadPairs`, so old outputs do not demonstrate the corrected counting workflow. Plotting uses the fitted `results/deseq2/dds.rds` produced by the same analysis.

The study design and sample independence still require confirmation from the study metadata. Rerunning does not validate those assumptions.

Author: Narges Shayesteh.
