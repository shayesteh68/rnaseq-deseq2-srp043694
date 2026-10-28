#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DESeq2))
output_dir <- "results/deseq2"
fitted_file <- file.path(output_dir, "dds.rds")
if (!file.exists(fitted_file)) stop("Run scripts/run_deseq2.R first; fitted object missing: ", fitted_file)
dds <- readRDS(fitted_file)
if (!inherits(dds, "DESeqDataSet") || !"condition" %in% names(colData(dds))) {
  stop("The fitted object is not a compatible DESeq2 analysis.")
}
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
pad_axis <- function(x, fraction = 0.08) {
  limits <- range(x, finite = TRUE)
  span <- diff(limits)
  if (!is.finite(span) || span == 0) span <- max(abs(limits), 1)
  limits + c(-1, 1) * span * fraction
}

png(file.path(output_dir, "pca_plot.png"), width = 2400, height = 2000, res = 300)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
colors <- ifelse(pca_data$condition == "SIM", "#D95F02", "#1B9E77")
plot(pca_data$PC1, pca_data$PC2,
     xlab = paste0("PC1: ", percent_var[1], "% variance"),
     ylab = paste0("PC2: ", percent_var[2], "% variance"),
     main = "PCA: SIM vs LGD (VST normalized)",
     col = colors, pch = 19, cex = 1.8,
     xlim = pad_axis(pca_data$PC1), ylim = pad_axis(pca_data$PC2))
text(pca_data$PC1, pca_data$PC2, labels = pca_data$name, pos = 3, cex = 0.7)
legend("topright", legend = levels(factor(pca_data$condition)),
       col = c("#1B9E77", "#D95F02"), pch = 19, bty = "n")
grid()
dev.off()

res <- as.data.frame(results(dds, contrast = c("condition", "SIM", "LGD"), alpha = 0.05))
res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), , drop = FALSE]
if (!nrow(res)) stop("No adjusted DESeq2 results available for the volcano plot.")
res$log_padj <- -log10(pmax(res$padj, .Machine$double.xmin))
col_vec <- rep("grey70", nrow(res))
up <- !is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange >= 1
down <- !is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange <= -1
col_vec[up] <- "#D95F02"
col_vec[down] <- "#1B9E77"

png(file.path(output_dir, "volcano_plot.png"), width = 2400, height = 2000, res = 300)
max_log_padj <- max(res$log_padj, na.rm = TRUE)
if (!is.finite(max_log_padj) || max_log_padj <= 0) max_log_padj <- 1
plot(res$log2FoldChange, res$log_padj,
     col = col_vec, pch = 20, cex = 0.8,
     xlab = expression(log[2] ~ "Fold Change (SIM / LGD)"),
     ylab = expression(-log[10] ~ "(Adjusted P-value)"),
     main = "Volcano Plot: SIM vs LGD",
     ylim = c(0, max_log_padj * 1.05))
abline(h = -log10(0.05), col = "black", lty = 2, lwd = 1)
abline(v = c(-1, 1), col = "black", lty = 2, lwd = 1)
legend("topleft",
       legend = c(paste0("Up in SIM (", sum(up), ")"), paste0("Down in SIM (", sum(down), ")"), "Not significant"),
       col = c("#D95F02", "#1B9E77", "grey70"), pch = 19, bty = "n")
grid()
dev.off()
message("Plots generated in ", output_dir)
