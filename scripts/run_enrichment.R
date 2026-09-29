#!/usr/bin/env Rscript
# ==============================================================================
# Functional Enrichment Analysis (GO & KEGG) using clusterProfiler
# Project: SRP043694 (Barrett's Esophagus: SIM vs LGD)
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

# 1. Define paths
input_file <- "results/deseq2/deg_SIM_vs_LGD.csv"
out_dir    <- "results/enrichment"
plots_dir  <- "results/plots"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

message(">>> Loading significant DEGs from: ", input_file)
if (!file.exists(input_file)) {
  stop("Error: File not found: ", input_file)
}

deg_data <- read.csv(input_file, stringsAsFactors = FALSE)

# Clean Ensembl ID
deg_data$ensembl <- sub("\\..*$", "", as.character(deg_data$gene_id))
all_sig_ensembl  <- unique(na.omit(deg_data$ensembl))
all_sig_ensembl  <- all_sig_ensembl[all_sig_ensembl != ""]

message(sprintf(">>> Loaded %d unique Ensembl gene IDs.", length(all_sig_ensembl)))

# 2. Convert Ensembl IDs to Entrez IDs & Symbols
message(">>> Mapping Ensembl IDs to Entrez IDs via org.Hs.eg.db...")
id_map <- bitr(all_sig_ensembl,
               fromType = "ENSEMBL",
               toType   = c("ENTREZID", "SYMBOL"),
               OrgDb    = org.Hs.eg.db)

entrez_genes <- unique(na.omit(id_map$ENTREZID))
message(sprintf(">>> Successfully mapped %d genes to Entrez IDs.", length(entrez_genes)))

# 3. GO Biological Process (BP)
message(">>> Running GO Biological Process enrichment...")
ego_bp <- enrichGO(gene          = entrez_genes,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)

if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
  write.csv(as.data.frame(ego_bp), file.path(out_dir, "GO_BP_results.csv"), row.names = FALSE)
  p_bp <- dotplot(ego_bp, showCategory = 15, title = "GO: Biological Process (SIM vs LGD)") +
          theme_minimal(base_size = 11)
  ggsave(file.path(plots_dir, "enrichment_GO_BP_dotplot.png"), plot = p_bp, width = 8, height = 7, dpi = 300)
  message("  [+] GO BP results and plot saved.")
} else {
  message("  [-] No significant GO BP terms found at pvalueCutoff = 0.05.")
}

# 4. KEGG Pathway Enrichment
message(">>> Running KEGG Pathway enrichment...")
ekegg <- tryCatch({
  enrichKEGG(gene         = entrez_genes,
             organism     = 'hsa',
             pvalueCutoff = 0.05)
}, error = function(e) {
  message("  [!] KEGG query note: ", e$message)
  return(NULL)
})

if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  ekegg_readable <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.csv(as.data.frame(ekegg_readable), file.path(out_dir, "KEGG_results.csv"), row.names = FALSE)
  p_kegg <- dotplot(ekegg, showCategory = 15, title = "KEGG Pathway Enrichment (SIM vs LGD)") +
            theme_minimal(base_size = 11)
  ggsave(file.path(plots_dir, "enrichment_KEGG_dotplot.png"), plot = p_kegg, width = 8, height = 6, dpi = 300)
  message("  [+] KEGG results and plot saved.")
} else {
  message("  [-] No significant KEGG pathways found at threshold 0.05.")
}

message(">>> Enrichment analysis completed successfully!")
