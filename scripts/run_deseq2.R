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
if (any(table(factor(meta$condition, levels = c("LGD", "SIM"))) < 2L)) {
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
res <- results(dds, contrast = c("condition", "SIM", "LGD"), alpha = 0.05)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
deg <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1, , drop = FALSE]
deg_up <- deg[deg$log2FoldChange > 0, , drop = FALSE]
deg_down <- deg[deg$log2FoldChange < 0, , drop = FALSE]

write.csv(res_df, file.path(output_dir, "deseq2_results_all.csv"), row.names = FALSE)
write.csv(deg, file.path(output_dir, "deg_SIM_vs_LGD.csv"), row.names = FALSE)
write.csv(deg_up, file.path(output_dir, "deg_upregulated_SIM_vs_LGD.csv"), row.names = FALSE)
write.csv(deg_down, file.path(output_dir, "deg_downregulated_SIM_vs_LGD.csv"), row.names = FALSE)
writeLines(c(
  "=== DESeq2 Analysis Summary ===",
  paste0("Total genes after low-count filtering: ", nrow(dds)),
  paste0("Significant DEGs (padj < 0.05, |log2FC| >= 1): ", nrow(deg)),
  paste0("Upregulated in SIM: ", nrow(deg_up)),
  paste0("Downregulated in SIM: ", nrow(deg_down)),
  "Design: ~ condition",
  "Contrast: SIM vs LGD (reference: LGD)"
), file.path(output_dir, "summary.txt"))
message("DESeq2 analysis completed for ", ncol(dds), " samples.")