###############################################################################
# GO/KEGG enrichment analysis for eGenes (clusterProfiler)
#
# Input:
#   eGenes.txt  (must contain a column named: gene_id; ENSEMBL IDs for human)
#
# Output:
#   eGenes_GO.txt
#   GO_enrichment_barplot.pdf
#   eGenes_KEGG_BH.txt
#   KEGG_enrichment_barplot.pdf
#
# Notes:
# - This script does not include raw data.
# - Set the working directory to the folder containing eGenes.txt.
###############################################################################

# ---- Working directory ----
setwd("D:/file/KEGG_GO")

# ---- Packages ----
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(tools)

# ---- Input ----
gene_list <- fread("eGenes.txt")   # column: gene_id (ENSEMBL)

###############################################################################
# Part 1: GO enrichment (ENSEMBL -> GO)
###############################################################################
ego <- enrichGO(
  gene          = gene_list$gene_id,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

write.table(as.data.frame(ego), "eGenes_GO.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# Select top 8 terms per ontology (BP/CC/MF) by -log10(P)
go_tbl <- as.data.table(as.data.frame(ego))
go_tbl$Description <- sapply(go_tbl$Description, toTitleCase)
go_tbl[, `-log10P` := -log10(pvalue)]

go_top <- go_tbl[, .SD[order(-`-log10P`)][1:min(.N, 8)], by = ONTOLOGY]
go_top$Description <- factor(go_top$Description, levels = rev(unique(go_top$Description)))

go_bar <- ggplot(go_top, aes(x = Description, y = `-log10P`, fill = ONTOLOGY)) +
  geom_bar(stat = "identity", width = 0.8, colour = "black", size = 0.15) +
  coord_flip() +
  theme_bw() +
  labs(x = "GO term", y = "-log10(P)", title = "GO enrichment (top terms per ontology)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 13),
    plot.title = element_text(size = 13, hjust = 0.5, face = "bold")
  )

ggsave("GO_enrichment_barplot.pdf", plot = go_bar, width = 8.5, height = 8)

###############################################################################
# Part 2: KEGG enrichment (ENSEMBL -> ENTREZ -> KEGG)
###############################################################################
# Convert ENSEMBL to ENTREZID
gi <- bitr(gene_list$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Keep only mapped genes for KEGG enrichment
kegg <- enrichKEGG(
  gene          = gi$ENTREZID,
  keyType       = "kegg",
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  qvalueCutoff  = 1
)

write.table(as.data.frame(kegg), "eGenes_KEGG_BH.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# Plot top 20 pathways by p-value
kegg_tbl <- as.data.table(as.data.frame(kegg))
kegg_tbl <- kegg_tbl[order(pvalue)][1:min(.N, 20)]
kegg_tbl$Description <- sapply(kegg_tbl$Description, toTitleCase)
kegg_tbl[, `-log10P` := -log10(pvalue)]
kegg_tbl$Description <- factor(kegg_tbl$Description, levels = rev(kegg_tbl$Description))

kegg_bar <- ggplot(kegg_tbl, aes(x = Description, y = `-log10P`)) +
  geom_bar(stat = "identity", width = 0.8, colour = "black", size = 0.15, fill = "#d0e7ed") +
  coord_flip() +
  theme_bw() +
  labs(x = "Pathway", y = "-log10(P)", title = "KEGG enrichment (top 20 pathways)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 13),
    plot.title = element_text(size = 13, hjust = 0.5, face = "bold")
  )

ggsave("KEGG_enrichment_barplot.pdf", plot = kegg_bar, width = 7.2, height = 8)
