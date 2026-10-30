#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
})
output_dir <- "results/deseq2"
plots_dir <- "results/plots"
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
fitted_file <- file.path(output_dir, "dds.rds")
if (!file.exists(fitted_file)) stop("Run scripts/run_deseq2.R first; fitted object missing: ", fitted_file)
dds <- readRDS(fitted_file)
if (!inherits(dds, "DESeqDataSet") || !"condition" %in% names(colData(dds))) {
  stop("The fitted object is not a compatible DESeq2 analysis.")
}
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pca_data$condition <- factor(pca_data$condition, levels = c("LGD", "SIM"))
percent_var <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text_repel(
    size = 3.2,
    box.padding = 0.45,
    point.padding = 0.35,
    min.segment.length = 0,
    max.overlaps = Inf,
    max.time = 5,
    seed = 42
  ) +
  scale_color_manual(values = c(LGD = "#1B9E77", SIM = "#D95F02")) +
  labs(
    title = "PCA: SIM vs LGD (VST normalized)",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Condition"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top",
    plot.margin = margin(10, 24, 10, 10)
  )
for (plot_dir in c(output_dir, plots_dir)) {
  ggsave(file.path(plot_dir, "pca_plot.png"), p_pca, width = 10, height = 7, dpi = 300, bg = "white")
}

sample_index <- match(pca_data$name, colnames(dds))
sample_qc <- data.frame(
  sample_id = pca_data$name,
  condition = as.character(pca_data$condition),
  filtered_count_sum = as.numeric(colSums(counts(dds, normalized = FALSE))[sample_index]),
  size_factor = as.numeric(sizeFactors(dds)[sample_index]),
  PC1 = pca_data$PC1,
  PC2 = pca_data$PC2,
  stringsAsFactors = FALSE
)
write.table(sample_qc, file.path(output_dir, "sample_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

results_file <- file.path(output_dir, "deseq2_results_all.csv")
if (!file.exists(results_file)) stop("Run scripts/run_deseq2.R first; results missing: ", results_file)
res <- read.csv(results_file, stringsAsFactors = FALSE)
if (!all(c("gene_id", "padj", "log2FoldChange") %in% names(res))) {
  stop("DESeq2 result table must include gene_id, padj, and ashr-shrunken log2FoldChange.")
}
res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), , drop = FALSE]
if (!nrow(res)) stop("No adjusted DESeq2 results available for the volcano plot.")
res$log_padj <- -log10(pmax(res$padj, .Machine$double.xmin))
col_vec <- rep("grey70", nrow(res))
up <- !is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange >= 1
down <- !is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange <= -1
col_vec[up] <- "#D95F02"
col_vec[down] <- "#1B9E77"

png(file.path(output_dir, "volcano_plot.png"), width = 2400, height = 2000, res = 300, bg = "white")
max_log_padj <- max(res$log_padj, na.rm = TRUE)
if (!is.finite(max_log_padj) || max_log_padj <= 0) max_log_padj <- 1
plot(res$log2FoldChange, res$log_padj,
     col = col_vec, pch = 20, cex = 0.8,
     xlab = expression(log[2] ~ "Fold Change (ashr-shrunken SIM / LGD)"),
     ylab = expression(-log[10] ~ "(Adjusted P-value)"),
     main = "Volcano Plot: SIM vs LGD (ashr-shrunken LFC)",
     ylim = c(0, max_log_padj * 1.05))
formal <- !is.na(res$padj_lfc_gt_1) & res$padj_lfc_gt_1 < 0.05
if (any(formal)) {
  text(res$log2FoldChange[formal], res$log_padj[formal], labels = res$gene_id[formal],
       pos = 2, cex = 0.75, col = "#A64B00")
}
abline(h = -log10(0.05), col = "black", lty = 2, lwd = 1)
abline(v = c(-1, 1), col = "black", lty = 2, lwd = 1)
legend("topleft",
       legend = c(paste0("Up in SIM (", sum(up), ")"), paste0("Down in SIM (", sum(down), ")"), "Not significant"),
       col = c("#D95F02", "#1B9E77", "grey70"), pch = 19, bty = "n", cex = 0.85)
grid()
dev.off()
if (!file.copy(file.path(output_dir, "volcano_plot.png"),
               file.path(plots_dir, "volcano_plot.png"), overwrite = TRUE)) {
  stop("Could not synchronize volcano_plot.png to ", plots_dir)
}
message("Plots and per-sample QC generated in ", output_dir)
