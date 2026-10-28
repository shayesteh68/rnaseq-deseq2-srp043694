#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

results_file <- "results/deseq2/deseq2_results_all.csv"
out_dir <- "results/enrichment"
plots_dir <- "results/plots"
if (!file.exists(results_file)) stop("DESeq2 results not found: ", results_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

res <- read.csv(results_file, stringsAsFactors = FALSE)
required <- c("gene_id", "padj", "log2FoldChange")
if (!all(required %in% names(res))) stop("DE results must contain: ", paste(required, collapse = ", "))
res <- res[!is.na(res$gene_id) & !is.na(res$log2FoldChange), , drop = FALSE]
res$ensembl <- sub("\\..*$", "", as.character(res$gene_id))
res <- res[nzchar(res$ensembl), , drop = FALSE]
if (!nrow(res)) stop("No usable Ensembl IDs in DE results.")

all_map <- bitr(unique(res$ensembl), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
universe <- unique(na.omit(all_map$ENTREZID))
if (!length(universe)) stop("No tested genes mapped to Entrez IDs.")

for (direction in c("upregulated_SIM", "downregulated_SIM")) {
  if (direction == "upregulated_SIM") {
    selected <- unique(res$ensembl[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange >= 1])
  } else {
    selected <- unique(res$ensembl[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange <= -1])
  }
  if (!length(selected)) {
    message(direction, ": no significant genes; skipping.")
    next
  }
  mapped <- bitr(selected, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  genes <- unique(na.omit(mapped$ENTREZID))
  if (length(genes) < 2L) {
    message(direction, ": fewer than two mapped genes; skipping.")
    next
  }

  go <- enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  if (!is.null(go) && nrow(as.data.frame(go))) {
    write.csv(as.data.frame(go), file.path(out_dir, paste0("GO_BP_", direction, ".csv")), row.names = FALSE)
    ggsave(file.path(plots_dir, paste0("enrichment_GO_BP_", direction, ".png")),
           dotplot(go, showCategory = 15, title = paste("GO BP:", direction)), width = 8, height = 7, dpi = 300)
  }

  kegg <- tryCatch(
    enrichKEGG(gene = genes, universe = universe, organism = "hsa", pAdjustMethod = "BH",
               pvalueCutoff = 0.05, qvalueCutoff = 0.2),
    error = function(e) {
      message("KEGG ", direction, " failed: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(kegg) && nrow(as.data.frame(kegg))) {
    readable <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    write.csv(as.data.frame(readable), file.path(out_dir, paste0("KEGG_", direction, ".csv")), row.names = FALSE)
    ggsave(file.path(plots_dir, paste0("enrichment_KEGG_", direction, ".png")),
           dotplot(readable, showCategory = 15, title = paste("KEGG:", direction)), width = 8, height = 6, dpi = 300)
  }
}
message("Direction-specific enrichment finished.")