# 🧬 End-to-End Bulk RNA-Seq Analysis Pipeline (SRP043694)

![Pipeline Status](https://img.shields.io/badge/Pipeline-Production--Ready-success?style=for-the-badge&logo=gnubash)
![R Version](https://img.shields.io/badge/R-4.3%2B-blue?style=for-the-badge&logo=r)
![Bioconductor](https://img.shields.io/badge/Bioconductor-DESeq2-brightgreen?style=for-the-badge)
![License](https://img.shields.io/badge/License-MIT-orange?style=for-the-badge)

A robust, reproducible, and automated bulk RNA-Seq pipeline for Differential Expression Analysis (DEA) and functional pathway enrichment, profiling disease progression in **Barrett's Esophagus** ($n = 14$).

## Key Results

| Finding | Result |
| --- | --- |
| Total Features Filtered | 24,115 genes |
| Statistically Significant DEGs | 142 genes |
| Upregulated in SIM | 90 genes |
| Downregulated in SIM | 52 genes |
| Top Enriched GO Terms | Keratinization, Epithelial Development |
| Top KEGG Pathways | ECM-receptor interaction, Focal adhesion |

## Repository Structure

```text
.
├── scripts/
│   ├── 01_fastp.sh
│   ├── 02_hisat2_sort.sh
│   ├── config.example.sh
│   ├── run_deseq2.R
│   ├── run_enrichment.R
│   └── plot_deseq2.R
├── plots/
│   ├── pca.png
│   ├── volcano.png
│   ├── ma.png
│   └── heatmap_top_genes.png
├── tables/
└── README.md
```

## Requirements

- Bash
- R 4.3 or later
- Bioconductor and DESeq2
- fastp
- HISAT2 and samtools
- R packages required by the analysis and enrichment scripts

## Execution

Run the commands from the repository root. Configure paths, sample information, reference genome, and annotation in a local copy of `scripts/config.example.sh` before running the upstream processing scripts.

### 1. Configure and run read quality control

```bash
cp scripts/config.example.sh scripts/config.sh
# Edit scripts/config.sh for your input files, output directories, and reference settings.
bash scripts/01_fastp.sh
```

### 2. Align reads and sort alignments

```bash
bash scripts/02_hisat2_sort.sh
```

### 3. Run differential expression, enrichment, and visualization

```bash
Rscript scripts/run_deseq2.R
Rscript scripts/run_enrichment.R
Rscript scripts/plot_deseq2.R
```

## Visualizations

### PCA of normalized expression

![PCA of normalized sample expression](plots/pca.png)

### Differential-expression overview

![Volcano plot of DESeq2 results](plots/volcano.png)

### Mean expression and fold change

![MA plot of DESeq2 results](plots/ma.png)

### Expression pattern of selected genes

![Heatmap of selected differentially expressed genes](plots/heatmap_top_genes.png)

## Author

**Narges Shayesteh**  
Biophysical & Bioinformatics Scientist | Specializing in NGS Data Pipelines & Drug Discovery

- GitHub: [@shayesteh68](https://github.com/shayesteh68)
- LinkedIn: [Narges Shayesteh](https://www.linkedin.com/in/narges-shayesteh)
