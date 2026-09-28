#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
})

output_dir <- "results/deseq2"
counts_file <- "results/counts/SRP043694_count_matrix_s0.tsv"
metadata_file <- "metadata/sample_sheet.tsv"
res_file <- file.path(output_dir, "deseq2_results_all.csv")

# 1) Load data & run VST for PCA
meta <- read.delim(metadata_file, stringsAsFactors = FALSE)
meta <- meta[!duplicated(meta$sample_id), ]
rownames(meta) <- meta$sample_id

counts_df <- read.delim(counts_file, row.names = 1, check.names = FALSE)
meta <- meta[match(colnames(counts_df), meta$sample_id), , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_df)),
  colData   = meta,
  design    = ~ condition
)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Save PCA Plot
png(file.path(output_dir, "pca_plot.png"), width = 2400, height = 2000, res = 300)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

colors <- ifelse(pca_data$condition == "SIM", "#D95F02", "#1B9E77")
plot(pca_data$PC1, pca_data$PC2,
     xlab = paste0("PC1: ", percent_var[1], "% variance"),
     ylab = paste0("PC2: ", percent_var[2], "% variance"),
     main = "PCA: SIM vs LGD (VST normalized)",
     col = colors, pch = 19, cex = 1.8,
     xlim = range(pca_data$PC1) * 1.2, ylim = range(pca_data$PC2) * 1.2)
text(pca_data$PC1, pca_data$PC2, labels = pca_data$name, pos = 3, cex = 0.7)
legend("topright", legend = levels(factor(pca_data$condition)),
       col = c("#1B9E77", "#D95F02"), pch = 19, bty = "n")
grid()
dev.off()

# 2) Volcano Plot
res <- read.csv(res_file, stringsAsFactors = FALSE)
res <- res[!is.na(res$padj), ]

res$log_padj <- -log10(res$padj)

col_vec <- rep("grey70", nrow(res))
col_vec[res$padj < 0.05 & res$log2FoldChange >= 1]  <- "#D95F02" # Up
col_vec[res$padj < 0.05 & res$log2FoldChange <= -1] <- "#1B9E77" # Down

png(file.path(output_dir, "volcano_plot.png"), width = 2400, height = 2000, res = 300)
plot(res$log2FoldChange, res$log_padj,
     col = col_vec, pch = 20, cex = 0.8,
     xlab = expression(log[2] ~ "Fold Change (SIM / LGD)"),
     ylab = expression(-log[10] ~ "(Adjusted P-value)"),
     main = "Volcano Plot: SIM vs LGD",
     ylim = c(0, max(res$log_padj, na.rm = TRUE) * 1.05))
abline(h = -log10(0.05), col = "black", lty = 2, lwd = 1)
abline(v = c(-1, 1), col = "black", lty = 2, lwd = 1)
legend("topleft",
       legend = c("Up in SIM (52)", "Down in SIM (84)", "Not Significant"),
       col = c("#D95F02", "#1B9E77", "grey70"), pch = 19, bty = "n")
grid()
dev.off()

cat("Plots generated successfully in", output_dir, "\n")
