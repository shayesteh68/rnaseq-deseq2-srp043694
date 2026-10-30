#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DESeq2))

counts_file <- Sys.getenv("COUNTS_FILE", "results/counts/SRP043694_count_matrix_s0.tsv")
metadata_file <- "metadata/sample_sheet.tsv"
output_dir <- "results/deseq2"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(counts_file)) stop("Counts file not found: ", counts_file)
if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)
meta <- read.delim(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("sample_id", "condition") %in% names(meta))) stop("Metadata needs sample_id and condition columns.")
if (anyNA(meta$sample_id) || any(!nzchar(meta$sample_id)) || anyDuplicated(meta$sample_id)) {
  stop("Metadata sample_id values must be present and unique.")
}
if (anyNA(meta$condition) || !all(meta$condition %in% c("LGD", "SIM"))) {
  stop("Conditions must be non-missing LGD or SIM values.")
}

counts_df <- read.delim(counts_file, row.names = 1, check.names = FALSE)
if (anyDuplicated(rownames(counts_df))) stop("Count matrix gene IDs must be unique.")
count_samples <- colnames(counts_df)
if (anyDuplicated(count_samples)) stop("Count matrix sample columns must be unique.")
if (!setequal(count_samples, meta$sample_id)) {
  stop("Count matrix columns and metadata sample IDs must match exactly.")
}
meta <- meta[match(count_samples, meta$sample_id), , drop = FALSE]
rownames(meta) <- meta$sample_id
count_matrix <- as.matrix(counts_df)
storage.mode(count_matrix) <- "numeric"
if (anyNA(count_matrix) || any(!is.finite(count_matrix)) ||
    any(count_matrix < 0) || any(count_matrix > .Machine$integer.max) ||
    any(count_matrix != floor(count_matrix))) {
  stop("Raw counts must be finite non-negative integers.")
}
storage.mode(count_matrix) <- "integer"
condition_counts <- table(factor(meta$condition, levels = c("LGD", "SIM")))
if (any(condition_counts < 2L)) {
  stop("At least two independent biological samples are required in each condition.")
}
meta$condition <- relevel(factor(meta$condition), ref = "LGD")

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = meta, design = ~ condition)
keep <- rowSums(counts(dds) >= 10) >= 2L
dds <- dds[keep, ]
if (nrow(dds) < 2L) stop("Too few genes remain after low-count filtering.")
dds <- DESeq(dds)
saveRDS(dds, file.path(output_dir, "dds.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

contrast <- c("condition", "SIM", "LGD")
res_mle <- results(dds, contrast = contrast, alpha = 0.05)
coef_name <- "condition_SIM_vs_LGD"
if (!coef_name %in% resultsNames(dds)) {
  stop("Expected coefficient missing; inspect resultsNames(dds): ", paste(resultsNames(dds), collapse = ", "))
}
res_shrunk <- lfcShrink(dds, coef = coef_name, type = "ashr")
res_lfc1 <- results(dds, contrast = contrast, alpha = 0.05,
                    lfcThreshold = 1, altHypothesis = "greaterAbs")
res_df <- as.data.frame(res_mle)
res_df$gene_id <- rownames(res_df)
idx_shrunk <- match(res_df$gene_id, rownames(res_shrunk))
idx_lfc1 <- match(res_df$gene_id, rownames(res_lfc1))
if (anyNA(idx_shrunk) || anyNA(idx_lfc1)) stop("DESeq2 result rows could not be aligned by gene ID.")
res_df$log2FoldChange_MLE <- res_df$log2FoldChange
res_df$lfcSE_MLE <- res_df$lfcSE
res_df$log2FoldChange <- res_shrunk$log2FoldChange[idx_shrunk]
res_df$lfcSE <- res_shrunk$lfcSE[idx_shrunk]
res_df$pvalue_lfc_gt_1 <- res_lfc1$pvalue[idx_lfc1]
res_df$padj_lfc_gt_1 <- res_lfc1$padj[idx_lfc1]

candidate <- !is.na(res_df$padj) & res_df$padj < 0.05 &
  !is.na(res_df$log2FoldChange) & abs(res_df$log2FoldChange) >= 1
deg <- res_df[candidate, , drop = FALSE]
deg_up <- deg[deg$log2FoldChange > 0, , drop = FALSE]
deg_down <- deg[deg$log2FoldChange < 0, , drop = FALSE]
formal_threshold <- !is.na(res_df$padj_lfc_gt_1) & res_df$padj_lfc_gt_1 < 0.05
deg_formal <- res_df[formal_threshold, , drop = FALSE]

write.csv(res_df, file.path(output_dir, "deseq2_results_all.csv"), row.names = FALSE)
write.csv(deg, file.path(output_dir, "deg_SIM_vs_LGD.csv"), row.names = FALSE)
write.csv(deg_up, file.path(output_dir, "deg_upregulated_SIM_vs_LGD.csv"), row.names = FALSE)
write.csv(deg_down, file.path(output_dir, "deg_downregulated_SIM_vs_LGD.csv"), row.names = FALSE)
write.csv(deg_formal, file.path(output_dir, "deg_formal_lfc_gt_1_SIM_vs_LGD.csv"), row.names = FALSE)
top_lines <- c(
  "Exploratory candidates ranked by Wald-test adjusted p-value and ashr-shrunken effect size.",
  "Candidate filter: padj < 0.05 and |ashr-shrunken log2FC| >= 1.",
  paste0("Formal Wald test of |log2FC| > 1 at BH padj < 0.05: ", nrow(deg_formal), " gene(s)."),
  ""
)
for (direction in c("upregulated_SIM", "downregulated_SIM")) {
  subset <- if (direction == "upregulated_SIM") deg_up else deg_down
  subset <- subset[order(subset$padj, -abs(subset$log2FoldChange)), , drop = FALSE]
  subset <- head(subset, 10L)
  top_lines <- c(top_lines, paste0("=== TOP 10 ", toupper(direction), " ==="))
  if (!nrow(subset)) {
    top_lines <- c(top_lines, "No genes passed the exploratory candidate filter.", "")
    next
  }
  top_table <- data.frame(
    rank = seq_len(nrow(subset)),
    gene_id = subset$gene_id,
    baseMean = subset$baseMean,
    log2FoldChange_ashr = subset$log2FoldChange,
    log2FoldChange_MLE = subset$log2FoldChange_MLE,
    padj = subset$padj,
    padj_lfc_gt_1 = subset$padj_lfc_gt_1,
    check.names = FALSE
  )
  top_lines <- c(top_lines, capture.output(write.table(
    top_table, sep = "\t", quote = FALSE, row.names = FALSE
  )), "")
}
while (length(top_lines) && !nzchar(tail(top_lines, 1L))) {
  top_lines <- head(top_lines, -1L)
}
writeLines(top_lines, file.path(output_dir, "top10_genes.txt"))
writeLines(c(
  "=== DESeq2 Analysis Summary ===",
  paste0("Samples analyzed: ", ncol(dds), " (LGD=", condition_counts[["LGD"]], ", SIM=", condition_counts[["SIM"]], ")"),
  paste0("Genes after low-count filtering: ", nrow(dds)),
  "Design: ~ condition",
  "Contrast: SIM vs LGD (reference: LGD)",
  "Wald-test p-values and padj test a non-zero log2 fold change.",
  paste0("MLE candidate filter (padj < 0.05, |MLE log2FC| >= 1): ",
         sum(!is.na(res_mle$padj) & res_mle$padj < 0.05 & abs(res_mle$log2FoldChange) >= 1)),
  paste0("ashr-shrunken candidate filter (padj < 0.05, |shrunken log2FC| >= 1): ",
         nrow(deg), " (up in SIM: ", nrow(deg_up), "; down in SIM: ", nrow(deg_down), ")"),
  paste0("Formal Wald test of |log2FC| > 1 at BH padj < 0.05: ", nrow(deg_formal)),
  "Reported log2FoldChange is ashr-shrunken; log2FoldChange_MLE and lfcSE_MLE preserve the unshrunken estimates.",
  "The formal threshold-test results are in padj_lfc_gt_1 and deg_formal_lfc_gt_1_SIM_vs_LGD.csv."
), file.path(output_dir, "summary.txt"))
message("DESeq2 analysis completed for ", ncol(dds), " samples.")
