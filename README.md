# RNA-Seq Differential Expression Analysis: SIM vs LGD in Barrett's Esophagus

## Project Overview
This repository provides a reproducible end-to-end differential expression analysis comparing **Specialized Intestinal Metaplasia (SIM)** against **Low-Grade Dysplasia (LGD)** in Barrett's Esophagus patients using public RNA-Seq data (**NCBI SRA: SRP043694**).

---

## Dataset Description
- **Accession:** NCBI SRA SRP043694
- **Cohort Size:** 14 mucosal biopsy samples (7 SIM vs 7 LGD)
- **Sequencing:** Illumina paired-end RNA-Seq
- **Reference Genome:** Human GRCh38 / Ensembl

---

## Project Structure
```text
RNAseq_SRP043694/
├── metadata/
│   ├── sample_sheet.tsv
│   ├── ena_fastq_manifest.tsv
│   └── reference_paths.txt
├── results/
│   ├── counts/
│   │   ├── SRP043694_count_matrix_s0.tsv
│   │   ├── SRP043694_featureCounts_s0.txt
│   │   └── SRP043694_featureCounts_s0.txt.summary
│   └── deseq2/
│       ├── deseq2_results_all.csv
│       ├── deg_SIM_vs_LGD.csv
│       ├── deg_upregulated_SIM_vs_LGD.csv
│       ├── deg_downregulated_SIM_vs_LGD.csv
│       ├── pca_plot.png
│       ├── volcano_plot.png
│       ├── summary.txt
│       └── top10_genes.txt
├── scripts/
│   ├── 01_fastp.sh
│   ├── 02_hisat2_sort.sh
│   ├── run_deseq2.R
│   └── plot_deseq2.R
└── README.md

Statistical Methodology
Differential expression modeling was conducted in R using DESeq2:

Low-count filtering: Genes with read count >= 10 in at least 2 samples were retained (23,551 genes).
Design formula: ~ condition (Reference level: LGD).
Contrast: SIM vs LGD (positive log2FoldChange indicates higher expression in SIM).
Thresholds: Benjamini-Hochberg adjusted p-value (padj) < 0.05 and |log2FoldChange| >= 1.0.
Results & Quality Control
Summary

Metric	Value
Genes tested (post-filter)	23,551
Total Significant DEGs	136
Upregulated in SIM	52
Downregulated in SIM	84
Visualizations
<p align=“center”>

<img src=“results/deseq2/pca_plot.png” width=“48%” />

<img src=“results/deseq2/volcano_plot.png” width=“48%” />

</p>

Key Biological Findings
Top Upregulated Genes in Metaplasia (SIM)
PGA4 (Pepsinogen A4): log2FC = 25.28, padj = 4.62e-13. Major gastric secretory marker, indicating preserved specialized secretory phenotype.
SSTR1 (Somatostatin Receptor 1): Endocrine regulator governing cellular proliferation and secretion.
IGHV1-46 & LTK: Markers associated with immune infiltration and mucosal inflammation.
Top Downregulated Genes in SIM (Elevated in LGD)
ZBTB16 (PLZF): log2FC = -1.97, padj = 2.20e-03. Key tumor suppressor and transcription factor; loss indicates progression toward dysplasia.
PRRX1 & TBX18: Mesenchymal transcription factors involved in developmental pathways and epithelial remodeling.
PYGO1: Component of the canonical Wnt/beta-catenin signaling cascade linked to dysplastic transformation.
How to Reproduce
1. Environment Setup
bash
conda create -n r_deseq_env -c conda-forge -c bioconda r-base=4.3 r-deseq2
conda activate r_deseq_env
2. Execute Pipeline
bash
# Run statistical testing
Rscript scripts/run_deseq2.R

# Generate plots
Rscript scripts/plot_deseq2.R
Author
Narges Shayesteh

Bioinformatics & Computational Biology

