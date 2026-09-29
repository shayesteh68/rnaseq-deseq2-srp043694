# End-to-End Bulk RNA-Seq Pipeline: Barrett's Esophagus Progression (SRP043694)

![RNA-Seq](https://img.shields.io/badge/Bioinformatics-RNA--Seq-blue.svg)
![R](https://img.shields.io/badge/Language-R%20%7C%20Bash-green.svg)
![Bioconductor](https://img.shields.io/badge/Platform-Bioconductor-teal.svg)
![License](https://img.shields.io/badge/License-MIT-orange.svg)

## 🔬 Project Overview
This repository contains a reproducible, end-to-end bulk RNA-Seq pipeline analyzing transcriptomic alterations in **Barrett's Esophagus (BE)** progression, focusing on the transition from **Specialized Intestinal Metaplasia (SIM)** to **Low-Grade Dysplasia (LGD)** using dataset **SRP043694** (GSE58640).

The objective is to delineate molecular dysregulations driving disease progression and characterize early diagnostic biomarkers prior to malignant transformation into Esophageal Adenocarcinoma (EAC).

---

## 📊 Workflow & Methodology
1. **Raw Read Quality Assessment:** `FastQC` (v0.12) & `MultiQC` for read quality, duplication rates, and base composition.
2. **Quality Trimming:** `Trimmomatic` (adapter clipping, SLIDINGWINDOW:4:20, MINLEN:36).
3. **Splice-Aware Alignment:** `HISAT2` indexed against Ensembl GRCh38 / GENCODE reference genome.
4. **Quantification:** `featureCounts` (Subread package) summarized at the gene level.
5. **Differential Expression Analysis:** `DESeq2` implementing the Wald test, independent filtering, and Benjamini-Hochberg FDR correction ($p_{\text{adj}} < 0.05$, $|\log_2\text{FC}| \ge 1$).
6. **Functional Annotation & Pathway Enrichment:** `clusterProfiler` & `org.Hs.eg.db` for Gene Ontology (GO - Biological Process) and KEGG pathway over-representation analysis.

---

## 📈 Key Findings & Biological Insights

### 1. Differential Expression Summary
From **23,551** quantified genes passing independent pre-filtering:
- **136 Significant Differentially Expressed Genes (DEGs)** ($p_{\text{adj}} < 0.05$, $|\log_2\text{FC}| \ge 1$)
  - **52 Upregulated** in LGD compared to SIM
  - **84 Downregulated** in LGD compared to SIM

| Feature | Gene Symbol | Log2FC | Adjusted p-value | Biological Relevance |
| :--- | :--- | :--- | :--- | :--- |
| **Top Upregulated** | `PGA4` | +6.24 | $3.12 \times 10^{-6}$ | Gastric mucosal phenotypic shift |
| **Top Downregulated** | `ZBTB16` | -3.85 | $7.84 \times 10^{-5}$ | PLZF transcription factor; tumor suppressor loss |

### 2. Functional & Pathway Enrichment (GO & KEGG)
- **BMP Signaling Disruption (`GO:0030514`, Fold Enrichment > 10):** Direct dysregulation of the Bone Morphogenetic Protein pathway, an established mechanism driving columnar differentiation and dysplastic epithelial reprogramming in Barrett's esophagus.
- **Digestive System Remodeling:** Significant enrichment in specialized gastrointestinal tract functional modules reflecting morphological glandular restructuring.

---

## 🖼 Visualizations

<div align="center">
  <table>
    <tr>
      <td align="center"><b>Sample Distances (PCA)</b></td>
      <td align="center"><b>Volcano Plot (DEGs)</b></td>
    </tr>
    <tr>
      <td><img src="results/plots/pca_plot.png" width="400"/></td>
      <td><img src="results/plots/volcano_plot.png" width="400"/></td>
    </tr>
    <tr>
      <td align="center"><b>GO Biological Process</b></td>
      <td align="center"><b>KEGG Pathway Enrichment</b></td>
    </tr>
    <tr>
      <td><img src="results/plots/enrichment_GO_BP_dotplot.png" width="400"/></td>
      <td><img src="results/plots/enrichment_KEGG_dotplot.png" width="400"/></td>
    </tr>
  </table>
</div>

---

## 📁 Repository Structure
```text
.
├── metadata/                  # Sample sheet and experimental metadata
│   └── SraRunTable.txt
├── scripts/                   # Reproducible pipeline scripts
│   ├── run_deseq2.R           # Statistical DEG testing via DESeq2
│   ├── plot_deseq2.R          # PCA and Volcano generation
│   └── run_enrichment.R       # GO & KEGG enrichment analysis
├── results/
│   ├── deseq2/                # Differential expression output tables
│   │   ├── deg_SIM_vs_LGD.csv
│   │   └── deseq2_results_all.csv
│   ├── enrichment/            # GO & KEGG functional output tables
│   │   ├── GO_BP_results.csv
│   │   └── KEGG_results.csv
│   └── plots/                 # High-resolution (300 DPI) publication figures
│       ├── pca_plot.png
│       ├── volcano_plot.png
│       ├── enrichment_GO_BP_dotplot.png
│       └── enrichment_KEGG_dotplot.png
└── README.md
```
## 💻 Environment & Reproducibility
- **R Version:** >= 4.3.0
- **Key Packages:** `DESeq2`, `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`, `ggplot2`
- Managed via `Miniforge/Mamba` environment.
```
# Clone the repository
git clone https://github.com/shayesteh68/rnaseq-deseq2-srp043694.git
cd rnaseq-deseq2-srp043694

# Run Enrichment Analysis
Rscript scripts/run_enrichment.R
```
**Author:** Narges Shayesteh
