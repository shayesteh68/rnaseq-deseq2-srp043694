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

# Remove only this script's prior outputs so an empty result cannot leave a stale table or plot.
for (direction in c("upregulated_SIM", "downregulated_SIM")) {
  unlink(file.path(out_dir, paste0("GO_BP_", direction, ".csv")))
  unlink(file.path(out_dir, paste0("KEGG_", direction, ".csv")))
  unlink(file.path(plots_dir, paste0("enrichment_GO_BP_", direction, ".png")))
  unlink(file.path(plots_dir, paste0("enrichment_KEGG_", direction, ".png")))
}
unlink(file.path(out_dir, c("GO_BP_results.csv", "KEGG_results.csv", "summary.txt")))
unlink(file.path(plots_dir, c("enrichment_GO_BP_dotplot.png", "enrichment_KEGG_dotplot.png")))

res <- read.csv(results_file, stringsAsFactors = FALSE)
required <- c("gene_id", "padj", "log2FoldChange", "log2FoldChange_MLE", "padj_lfc_gt_1")
if (!all(required %in% names(res))) stop("DE results must contain: ", paste(required, collapse = ", "))
res <- res[!is.na(res$gene_id) & !is.na(res$log2FoldChange), , drop = FALSE]
res$ensembl <- sub("\\..*$", "", as.character(res$gene_id))
res <- res[nzchar(res$ensembl), , drop = FALSE]
if (!nrow(res)) stop("No usable Ensembl IDs in DE results.")

all_map <- bitr(unique(res$ensembl), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
universe <- unique(na.omit(all_map$ENTREZID))
if (!length(universe)) stop("No tested genes mapped to Entrez IDs.")
summary_lines <- c(
  "Directional enrichment using ashr-shrunken log2 fold changes.",
  "Candidate threshold: DESeq2 Wald padj < 0.05 and |shrunken log2FC| >= 1.",
  "A separate formal |log2FC| > 1 test is reported in results/deseq2/summary.txt.",
  paste0("DESeq2 result rows supplied: ", nrow(res)),
  paste0("Unique tested Entrez background IDs: ", length(universe)),
  paste0("Ensembl IDs that failed Entrez mapping: ", length(setdiff(unique(res$ensembl), unique(all_map$ENSEMBL))))
)
go_tables <- list()
kegg_tables <- list()
ratio_value <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(parts, function(p) as.numeric(p[1]) / as.numeric(p[2]), numeric(1))
}
aggregate_plot <- function(tbl, path, title) {
  if (!nrow(tbl)) {
    png(path, width = 2100, height = 1200, res = 300, bg = "white")
    plot.new()
    text(0.5, 0.5, "No terms passed the enrichment thresholds", cex = 1.2)
    title(main = title)
    dev.off()
    return(invisible(NULL))
  }
  tbl$ratio_numeric <- ratio_value(tbl$GeneRatio)
  tbl$score <- -log10(pmax(as.numeric(tbl$p.adjust), .Machine$double.xmin))
  plot_data <- do.call(rbind, lapply(split(tbl, tbl$direction), function(x) {
    x <- x[order(x$p.adjust, -x$Count), , drop = FALSE]
    head(x, 15L)
  }))
  p <- ggplot(plot_data, aes(x = ratio_numeric, y = reorder(Description, ratio_numeric),
                             color = direction, size = Count)) +
    geom_point(alpha = 0.85) +
    facet_wrap(~ direction, scales = "free_y") +
    theme_bw(base_size = 11) +
    labs(title = title, x = "Gene ratio", y = NULL, color = "Direction", size = "Genes")
  ggsave(path, p, width = 11, height = 8, dpi = 300, bg = "white")
}

for (direction in c("upregulated_SIM", "downregulated_SIM")) {
  is_up <- direction == "upregulated_SIM"
  selected <- unique(res$ensembl[!is.na(res$padj) & res$padj < 0.05 &
    if (is_up) res$log2FoldChange >= 1 else res$log2FoldChange <= -1])
  if (!length(selected)) {
    summary_lines <- c(summary_lines, paste0(direction, ": no candidate genes."))
    next
  }
  mapped <- bitr(selected, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  genes <- unique(na.omit(mapped$ENTREZID))
  if (length(genes) < 2L) {
    summary_lines <- c(summary_lines, paste0(direction, ": ", length(selected),
      " candidate Ensembl genes; fewer than two mapped Entrez genes; enrichment skipped."))
    next
  }
  summary_lines <- c(summary_lines, paste0(direction, ": ", length(selected),
    " candidate Ensembl genes; ", length(genes), " unique Entrez genes."))

  go <- enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  go_df <- if (is.null(go)) data.frame() else as.data.frame(go)
  if (nrow(go_df)) {
    go_df$direction <- direction
    go_tables[[direction]] <- go_df
    write.csv(go_df, file.path(out_dir, paste0("GO_BP_", direction, ".csv")), row.names = FALSE)
    ggsave(file.path(plots_dir, paste0("enrichment_GO_BP_", direction, ".png")),
           dotplot(go, showCategory = 15, title = paste("GO BP:", direction)),
           width = 8, height = 7, dpi = 300, bg = "white")
  }
  summary_lines <- c(summary_lines, paste0(direction, ": GO BP terms passing cutoffs = ", nrow(go_df), "."))

  kegg <- tryCatch(
    enrichKEGG(gene = genes, universe = universe, organism = "hsa", pAdjustMethod = "BH",
               pvalueCutoff = 0.05, qvalueCutoff = 0.2),
    error = function(e) {
      summary_lines <<- c(summary_lines, paste0(direction, ": KEGG request failed: ", conditionMessage(e)))
      NULL
    }
  )
  kegg_df <- if (is.null(kegg)) data.frame() else as.data.frame(kegg)
  if (nrow(kegg_df)) {
    readable <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    kegg_df <- as.data.frame(readable)
    kegg_df$direction <- direction
    kegg_tables[[direction]] <- kegg_df
    write.csv(kegg_df, file.path(out_dir, paste0("KEGG_", direction, ".csv")), row.names = FALSE)
    ggsave(file.path(plots_dir, paste0("enrichment_KEGG_", direction, ".png")),
           dotplot(readable, showCategory = 15, title = paste("KEGG:", direction)),
           width = 8, height = 6, dpi = 300, bg = "white")
  }
  summary_lines <- c(summary_lines, paste0(direction, ": KEGG terms passing cutoffs = ", nrow(kegg_df), "."))
}

combine_tables <- function(x) {
  x <- Filter(function(z) nrow(z) > 0L, x)
  if (!length(x)) return(data.frame(direction = character(), Description = character()))
  do.call(rbind, x)
}
go_all <- combine_tables(go_tables)
kegg_all <- combine_tables(kegg_tables)
write.csv(go_all, file.path(out_dir, "GO_BP_results.csv"), row.names = FALSE)
write.csv(kegg_all, file.path(out_dir, "KEGG_results.csv"), row.names = FALSE)
aggregate_plot(go_all, file.path(plots_dir, "enrichment_GO_BP_dotplot.png"), "GO BP enrichment by direction")
aggregate_plot(kegg_all, file.path(plots_dir, "enrichment_KEGG_dotplot.png"), "KEGG enrichment by direction")
summary_lines <- c(summary_lines, "KEGG uses the live KEGG REST service; results can change when the pathway database changes.")
writeLines(summary_lines, file.path(out_dir, "summary.txt"))
message("Direction-specific and combined enrichment outputs saved.")