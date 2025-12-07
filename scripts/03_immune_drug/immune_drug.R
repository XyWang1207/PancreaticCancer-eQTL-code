###############################################################################
# eGene annotation: immune infiltration / drug association / mutation 
# (simple copy-and-run version for GitHub)
#
# What this script does (high-level):
# 1) Correlate eGene expression with immune infiltration scores (Spearman)
# 2) Merge correlation results with eGenes list; apply multiple testing correction (BH/FDR)
# 3) Prepare random subset for dot-plot visualization (optional)
# 4) Drug-gene correlation and pathway mapping (table joins + summary)
# 5) MAF/TMB summary and oncoplot (maftools)
#
# Notes:
# - Please change setwd() paths to your local folders.
# - This repo does NOT distribute raw TCGA data or proprietary portal outputs.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

###############################################################################
# Part A: Immune infiltration correlation (Spearman)
###############################################################################
setwd("D:/file/immune")

# Inputs (rename your local files accordingly)
# - immune_scores.txt: immune infiltration matrix (samples x infiltration features)
# - expression_matrix.txt: expression matrix (genes x samples), first column = gene name
immune_scores <- fread("immune_scores.txt", header = TRUE)
expr_raw      <- fread("expression_matrix.txt", header = TRUE)

# Assume first column is gene name; remaining columns are samples
gene_name <- expr_raw[[1]]
gene_expr <- expr_raw[, -1, with = FALSE]   # genes x samples (numeric)

# Immune scores: keep numeric columns only (exclude sample ID column if present)
# If the first column is sample ID, drop it:
immu_mat <- immune_scores[, -1, with = FALSE]
feature_names <- colnames(immu_mat)

# Compute Spearman correlation for each (feature, gene)
outTab <- data.frame()
for (k in seq_along(feature_names)) {
  x <- as.numeric(immu_mat[[k]])
  for (j in seq_along(gene_name)) {
    t <- as.numeric(gene_expr[j, ])
    n <- cor.test(t, x, method = "spearman")
    outTab <- rbind(outTab, data.frame(
      feature = feature_names[k],
      gene_name = gene_name[j],
      p = n$p.value,
      r = as.numeric(n$estimate),
      stringsAsFactors = FALSE
    ))
  }
}

write.table(outTab, "immune_gene_spearman.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

###############################################################################
# Part B: Merge with eGenes and apply FDR (BH) per immune feature
###############################################################################
setwd("D:/file/immune")

eGenes <- fread("eGenes.txt", header = TRUE)              # must contain a column named "gene_name"
immu_corr <- fread("immune_gene_spearman.txt", header = TRUE)

# Ensure column name exists
if (!("gene_name" %in% names(immu_corr))) {
  stop("immune_gene_spearman.txt must contain a column named 'gene_name'.")
}
if (!("gene_name" %in% names(eGenes))) {
  # If your eGenes file uses a different column name, change it here:
  # names(eGenes)[1] <- "gene_name"
  stop("eGenes.txt must contain a column named 'gene_name'. Please rename accordingly.")
}

gene_imm <- merge(eGenes, immu_corr, by = "gene_name")

# P<0.05 subset (optional output)
gene_imm_p <- gene_imm[gene_imm$p < 0.05, ]

# Apply BH FDR within each immune feature (more robust than hard-coded slicing)
gene_imm$FDR <- ave(gene_imm$p, gene_imm$feature, FUN = function(v) p.adjust(v, method = "BH"))

gene_imm_fdr <- gene_imm[gene_imm$FDR < 0.05, ]

fwrite(gene_imm_fdr, "gene_immune_corr_FDR.txt", sep = "\t", quote = FALSE)
fwrite(gene_imm_p,   "gene_immune_corr_P005.txt", sep = "\t", quote = FALSE)

# Summary counts per immune feature
num <- as.data.frame(table(gene_imm_fdr$feature))
colnames(num) <- c("feature", "n_genes")
fwrite(num, "gene_immune_corr_FDR_counts.txt", sep = "\t", quote = FALSE)

###############################################################################
# Part C (optional): sample 30 genes for dot-plot data export
###############################################################################
set.seed(1)
if (nrow(gene_imm_fdr) >= 30) {
  sel30 <- data.frame(gene_name = sample(unique(gene_imm_fdr$gene_name), size = 30, replace = FALSE))
  sel30 <- merge(sel30, gene_imm_fdr, by = "gene_name")
  write.table(sel30, file = "immune_correlation_plot_data.txt",
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  # If you need a fixed order of immune features, set factor levels here (optional).
  # Otherwise, comment this out.
  sel30$feature <- factor(sel30$feature, levels = unique(sel30$feature))
  
  p1 <- ggplot(sel30, aes(x = gene_name, y = feature, size = -log10(p), color = r)) +
    geom_point() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p1)
}

###############################################################################
# Part D: Drug-gene correlation / pathway mapping (table joins)
###############################################################################
setwd("D:/file/drug")

drug_pair <- fread("drug_pair.txt", header = TRUE)  # must contain Compound/RNA mapping columns
eGenes2   <- fread("eGenes.txt", header = TRUE)

if (ncol(eGenes2) >= 2 && !("RNA.molecule" %in% names(eGenes2))) {
  names(eGenes2)[2] <- "RNA.molecule"
}

# Merge with drug-gene pairs
gene_drug <- merge(eGenes2, drug_pair, by = "RNA.molecule")

# Portal output (ANOVA) file: rename to a generic name locally
path <- fread("PAAD_ANOVA.csv", header = TRUE)
path <- path[, c(1,4,5,6), with = FALSE]   # keep minimal columns as in your original
path <- path[!duplicated(path[[1]]), ]     # drop duplicated drug names

# Make sure drug name column is "Compound" for joining
names(path)[1] <- "Compound"

gene_drug2 <- left_join(gene_drug, path, by = "Compound")
gene_drug3 <- na.omit(gene_drug2)

# Correlation direction (1=positive, 0=non-positive) based on Pearson.est.
if (!("Pearson.est." %in% names(gene_drug3))) {
  # If your column name differs, change it here.
  stop("Expected column 'Pearson.est.' not found in merged drug table.")
}
gene_drug3$correlation <- ifelse(gene_drug3$Pearson.est. > 0, "1", "0")

# CrossTable summary (requires gmodels)
if (!requireNamespace("gmodels", quietly = TRUE)) {
  message("Package 'gmodels' not installed; skipping CrossTable output.")
} else {
  library(gmodels)
  tab <- data.frame(CrossTable(gene_drug3$Compound, gene_drug3$correlation))
  fwrite(tab, "drug_correlation_crosstable.txt", sep = "\t", quote = FALSE)
}

fwrite(gene_drug3, "drug_correlation_all.txt", sep = "\t", quote = FALSE)

###############################################################################
# Part E (optional): select drugs of interest and prepare plot table
###############################################################################

gene_drug_wide <- gene_drug3

drug_plot <- c(
  "GSK690693","PIK-93","YM201636","OSI-027","AKT inhibitor VIII","Idelalisib",
  "ZSTK474","AZD8055","Omipalisib","AS605240","Temsirolimus","PF-4708671",
  "AZD6482","Dactolisib","MK-2206","Pictilisib","Motesanib"
)

plot_df <- gene_drug_wide %>% filter(Compound %in% drug_plot)

# Sample 30 genes for plotting (optional)
set.seed(1)
if (nrow(plot_df) >= 30) {
  sample30 <- data.frame(RNA.molecule = sample(unique(plot_df$RNA.molecule), size = 30, replace = FALSE))
  plot_data <- merge(sample30, plot_df, by = "RNA.molecule")
  
  # The following columns must exist; adjust if your portal output differs:
  # - Pearson.fdr, Pearson.est.
  if (!all(c("Pearson.fdr","Pearson.est.") %in% names(plot_data))) {
    message("Pearson.fdr / Pearson.est. not found; skip ggplot.")
  } else {
    p2 <- ggplot(plot_data, aes(x = Compound, y = RNA.molecule,
                                size = -log10(Pearson.fdr), color = Pearson.est.)) +
      geom_point() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(p2)
  }
  
  # Export minimal columns for downstream figure generation
  keep_cols <- intersect(c("RNA.molecule","Compound","Pearson.est.","Pearson.fdr"), names(plot_data))
  write.table(plot_data[, keep_cols, drop = FALSE],
              "drug_gene_correlation_plot_data.txt",
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

###############################################################################
# Part F: MAF/TMB summary and oncoplot (maftools)
###############################################################################
setwd("D:/file/tmb_maf")

# Inputs:
# - PAAD_maftools.maf : MAF file
# - top20_genes.txt   : one gene symbol per line (no header)
if (!requireNamespace("maftools", quietly = TRUE)) {
  stop("Package 'maftools' is required for the MAF section.")
}
library(maftools)

top20 <- fread("top20_genes.txt", header = FALSE)

acc <- read.maf(maf = "PAAD_maftools.maf")
getSampleSummary(acc)
getGeneSummary(acc)

write.mafSummary(maf = acc, basename = "PAAD")
plotmafSummary(maf = acc, rmOutlier = TRUE, addStat = "median",
               dashboard = TRUE, titvRaw = FALSE)

# Variant classification colors (optional; can be removed if not needed)
vc_cols <- c("#3288BD","#DE5A64","#FEE08B","#437C69","#E6F598",
             "#F46D43","#5E4FA2","#66C2A5","#FFFFbF")
names(vc_cols) <- c(
  "Nonsense_Mutation","Missense_Mutation","Frame_Shift_Del","Multi_Hit",
  "Frame_Shift_Ins","In_Frame_Ins","Splice_Site","In_Frame_Del"
)

oncoplot(maf = acc, colors = vc_cols, genes = top20$V1)
