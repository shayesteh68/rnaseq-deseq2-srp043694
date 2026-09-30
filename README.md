# 🧬 End-to-End Bulk RNA-Seq Analysis Pipeline (SRP043694)

![Pipeline Status](https://img.shields.io/badge/Pipeline-Production--Ready-success?style=for-the-badge&logo=snakemake)
![R Version](https://img.shields.io/badge/R-4.3%2B-blue?style=for-the-badge&logo=r)
![Bioconductor](https://img.shields.io/badge/Bioconductor-DESeq2-brightgreen?style=for-the-badge)
![License](https://img.shields.io/badge/License-MIT-orange?style=for-the-badge)

A robust, reproducible, and automated bulk RNA-Seq pipeline for Differential Expression Analysis (DEA) and functional pathway enrichment, profiling disease progression in **Barrett's Esophagus** ($n = 14$).

## 🔬 Dataset & Experimental Context

* **NCBI SRA Accession:** `SRP043694`
* **Study Cohort:** 14 paired/grouped human endoscopic biopsy profiles representing Barrett's Esophagus histological stages.
* **Primary Contrast:** Specialized Intestinal Metaplasia (**SIM**) vs. Low-Grade Dysplasia (**LGD**).
* **Reference Genome:** GRCh38 / Ensembl Release 109.

---

## 📈 Key Quantitative Discoveries

| Metric / Analysis | Value / Finding | Interpretation |
| :--- | :--- | :--- |
| **Total Features Filtered** | 24,115 genes | Pre-filtered ($>10$ read counts across $\ge 3$ samples) |
| **Statistically Significant DEGs** | **142 genes** | Threshold: $|\log_2\text{FC}| \ge 1.0$, $p_{\text{adj}} < 0.05$ (BH FDR) |
| **Upregulated in SIM** | **90 genes** | Epithelial differentiation & barrier markers |
| **Downregulated in SIM** | **52 genes** | Dysplasia-associated & cell-cycle deregulation markers |
| **Top Enriched GO Terms** | Keratinization, Epithelial Development | Confirms phenotypic mucosal remodeling |
| **Top KEGG Pathways** | ECM-receptor interaction, Focal adhesion | Structural/cytoskeletal reorganization |

---

## 📁 Repository Architecture

```text
RNAseq_SRP043694/
├── metadata/
│   └── samples.csv                 # Sample annotations & clinical phenotype metadata
├── scripts/
│   ├── 01_deseq2_analysis.R        # Normalization, dispersion estimation & DEA
│   ├── 02_pathway_enrichment.R     # ClusterProfiler (GO BP, MF, CC & KEGG)
│   └── 03_generate_publication_plots.R # Volcano, PCA, MA & Heatmap generator
├── results/
│   ├── tables/                     # Full & filtered DEG tables (TSV/CSV)
│   └── plots/                      # 300+ DPI publication figures
├── environment.yml                # Conda / Mamba reproducible environment
├── LICENSE
└── README.md
```

---

## ⚙️ Installation & Usage

### 1. Clone & Set Up Environment
```bash
git clone https://github.com/shayesteh68/rnaseq-deseq2-srp043694.git
cd rnaseq-deseq2-srp043694

# Create conda environment
conda env create -f environment.yml
conda activate bulk-rnaseq-env
```

### 2. Execute Full Workflow
```bash
# Run Differential Expression Analysis
Rscript scripts/01_deseq2_analysis.R

# Run Functional Enrichment Analysis
Rscript scripts/02_pathway_enrichment.R

# Generate Figures
Rscript scripts/03_generate_publication_plots.R
```

---

## Visualizations & Results

<p align="center">
  <img src="results/plots/pca_plot.png" width="48%" alt="PCA Plot" />
  <img src="results/plots/volcano_plot.png" width="48%" alt="Volcano Plot" />
</p>
<p align="center">
  <img src="results/plots/enrichment_GO_BP_dotplot.png" width="48%" alt="GO Enrichment" />
  <img src="results/plots/enrichment_KEGG_dotplot.png" width="48%" alt="KEGG Pathway" />
</p>

## 👩‍💻 Author & Contact

**Narges Shayesteh**  
*Biophysical & Bioinformatics Scientist | Specializing in NGS Data Pipelines & Drug Discovery*  
* **GitHub:** [@shayesteh68](https://github.com/shayesteh68)  
* **LinkedIn:** [Narges Shayesteh](https://www.linkedin.com/in/shayesteh)
