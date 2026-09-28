#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
})

# 0) Paths
counts_file   <- "results/counts/SRP043694_count_matrix_s0.tsv"
metadata_file <- "metadata/sample_sheet.tsv"
output_dir    <- "results/deseq2"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1) Reading Metadata
if (!file.exists(counts_file)) stop("Counts file not found: ", counts_file)
if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)

meta <- read.delim(metadata_file, stringsAsFactors = FALSE)
if (!all(c("sample_id", "condition") %in% colnames(meta))) {
  stop("Metadata file is missing required columns: sample_id or condition")
}

meta <- meta[!duplicated(meta$sample_id), ]
meta$sample_id <- as.character(meta$sample_id)
meta$condition <- as.factor(meta$condition)

# 2) Reading Counts Matrix
counts_df <- read.delim(counts_file, row.names = 1, check.names = FALSE)

# 3) Align Samples
count_samples <- colnames(counts_df)
meta_samples  <- meta$sample_id

if (!all(count_samples %in% meta_samples)) stop("Counts contain samples not in metadata!")
if (!all(meta_samples %in% count_samples)) stop("Metadata contains samples not in counts!")

meta <- meta[match(count_samples, meta$sample_id), , drop = FALSE]
rownames(meta) <- meta$sample_id

# 4) Reference Level
if (!("LGD" %in% levels(meta$condition))) stop("Reference level 'LGD' not found.")
meta$condition <- relevel(meta$condition, ref = "LGD")

# 5) Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_df)),
  colData   = meta,
  design    = ~ condition
)

# Filter low counts: at least 2 samples with count >= 10
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
dds <- DESeq(dds)

# 6) Results: SIM vs LGD
res <- results(dds, contrast = c("condition", "SIM", "LGD"))

# 7) DEG Identification (padj < 0.05 & |log2FC| >= 1)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

deg <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1, ]
deg_up   <- deg[deg$log2FoldChange > 0, ]
deg_down <- deg[deg$log2FoldChange < 0, ]

# 8) Save Results
write.csv(res_df, file.path(output_dir, "deseq2_results_all.csv"), row.names = FALSE)
write.csv(deg, file.path(output_dir, "deg_SIM_vs_LGD.csv"), row.names = FALSE)
write.csv(deg_up, file.path(output_dir, "deg_upregulated_SIM_vs_LGD.csv"), row.names = FALSE)
write.csv(deg_down, file.path(output_dir, "deg_downregulated_SIM_vs_LGD.csv"), row.names = FALSE)

# 9) Summary
summary_txt <- c(
  "=== DESeq2 Analysis Summary ===",
  paste0("Total genes after low-count filtering: ", nrow(dds)),
  paste0("Significant DEGs (padj < 0.05, |log2FC| >= 1): ", nrow(deg)),
  paste0("Upregulated in SIM: ", nrow(deg_up)),
  paste0("Downregulated in SIM: ", nrow(deg_down)),
  "Design: ~ condition",
  "Contrast: SIM vs LGD (Ref: LGD)"
)

writeLines(summary_txt, file.path(output_dir, "summary.txt"))
cat("\n", paste(summary_txt, collapse = "\n"), "\n\n")
