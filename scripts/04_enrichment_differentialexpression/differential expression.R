###############################################################################
# GEO differential expression (non-paired and paired) for selected target genes
#
# This script follows the GEO2R-style workflow:
#   1) Download GEO expression matrices via GEOquery
#   2) Apply log2 transform if needed (GEO2R logic)
#   3) Map probes to gene IDs and collapse probes to gene-level expression (limma::avereps)
#   4) Extract expression of predefined target genes
#   5) Differential expression with limma:
#        - GSE62165: non-paired tumor vs normal
#        - GSE15471: paired tumor vs normal (patient blocking)
#
# Outputs include:
#   - probe-level and gene-level expression matrices
#   - GraphPad-friendly tables for the selected genes
#   - limma P values and BH-adjusted P values for the selected genes
###############################################################################

###############################################################################
# Part 1: GSE62165 (non-paired tumor vs normal)
###############################################################################
setwd("D:/file/GSE62165")
rm(list = ls())

library(GEOquery)
library(limma)

## Download expression matrix from GEO
gset <- getGEO("GSE62165", GSEMatrix = TRUE, AnnotGPL = FALSE)[[1]]
ex <- exprs(gset)

## Check whether log2 transformation is needed (GEO2R logic)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

## Probe -> Ensembl (ENSG) mapping from platform annotation
fd <- gset@featureData@data
map <- data.frame(PROBEID = rownames(fd), ENSG = fd$Ensembl, stringsAsFactors = FALSE)
map <- map[!is.na(map$ENSG) & map$ENSG != "" & map$ENSG != "---", ]

## Remove ambiguous probes annotated to multiple genes ("///")
map <- map[!grepl("///", map$ENSG), ]
map <- map[grepl("^ENSG", map$ENSG), ]

## Collapse probe-level expression to gene-level expression (average across probes)
keep <- intersect(map$PROBEID, rownames(ex))
ex2 <- ex[keep, , drop = FALSE]
map2 <- map[match(rownames(ex2), map$PROBEID), ]
gene_ex <- avereps(ex2, ID = map2$ENSG)

## Export probe-level and gene-level expression matrices
write.table(
  data.frame(PROBEID = rownames(ex), ex, check.names = FALSE),
  "GSE62165_ex_probe.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(ENSG = rownames(gene_ex), gene_ex, check.names = FALSE),
  "GSE62165_ex_gene.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

## Extract expression values for target genes (ENSG IDs)
genes15 <- c(
  "ENSG00000127463","ENSG00000007341","ENSG00000143971","ENSG00000122641",
  "ENSG00000139515","ENSG00000279314","ENSG00000185024","ENSG00000184887",
  "ENSG00000211896","ENSG00000090863","ENSG00000240338","ENSG00000153774",
  "ENSG00000261783","ENSG00000161395","ENSG00000073605"
)

expr15 <- gene_ex[intersect(genes15, rownames(gene_ex)), , drop = FALSE]
missing <- setdiff(genes15, rownames(expr15))
if (length(missing)) message("These ENSG IDs were not found in the platform annotation/expression matrix: ",
                             paste(missing, collapse = ", "))

## Sample grouping (GEO2R group membership string)
gsms <- paste0("10111111111101111111111101111111111111101111110111",
               "11111111111111111111011111111111111101111111111111",
               "1111011011101111101111011110111")
sml <- strsplit(gsms, split = "")[[1]]
gs <- factor(sml)
levels(gs) <- make.names(c("normal", "tumor"))
group <- factor(gs, levels = c("normal", "tumor"))

## Export GraphPad-friendly table (one row per sample)
out_wide <- data.frame(Sample = colnames(expr15), Group = group, t(expr15), check.names = FALSE)
write.csv(out_wide, "GSE62165_15ENSG_forGraphPad.csv", row.names = FALSE)

## Differential expression testing using limma (tumor vs normal)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr15, design)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(tumor - normal, levels = design)))

res <- topTable(fit2, number = Inf, sort.by = "none")
res$ENSG <- rownames(res)

## Export limma results (P values and BH-adjusted P values)
write.csv(res[, c("ENSG","logFC","AveExpr","t","P.Value","adj.P.Val")],
          "GSE62165_15ENSG_limma_Pvalues.csv", row.names = FALSE)

###############################################################################
# Part 2: GSE15471 (paired tumor vs normal)
###############################################################################
setwd("D:/file/GSE15471")
rm(list = ls())

library(GEOquery)
library(limma)

## Download expression matrix from GEO (choose the matrix with the most samples if multiple exist)
gset_list <- getGEO("GSE15471", GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gset_list) > 1) {
  ns <- sapply(gset_list, function(x) ncol(exprs(x)))
  gset <- gset_list[[which.max(ns)]]
} else {
  gset <- gset_list[[1]]
}

## Make feature variable names syntactically valid (as in GEO2R)
fvarLabels(gset) <- make.names(fvarLabels(gset))

ex <- exprs(gset)

## Check whether log2 transformation is needed (GEO2R logic)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
  exprs(gset) <- ex
}

## Probe -> Gene mapping (auto-detect Ensembl Gene ID or Entrez Gene ID)
fd <- gset@featureData@data

ens_col <- grep("ensembl", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
gid_col <- grep("^Gene\\.ID$|^Gene ID$|entrez", colnames(fd), ignore.case = TRUE, value = TRUE)[1]

if (!is.na(ens_col)) {
  map <- data.frame(PROBEID = rownames(fd), GENEID = fd[[ens_col]], stringsAsFactors = FALSE)
  id_label <- "ENSG"
} else if (!is.na(gid_col)) {
  map <- data.frame(PROBEID = rownames(fd), GENEID = fd[[gid_col]], stringsAsFactors = FALSE)
  id_label <- "GeneID"
} else {
  stop("No Ensembl or Gene ID column found in featureData(gset). Please provide the platform annotation table to map probes.")
}

map <- map[!is.na(map$GENEID) & map$GENEID != "" & map$GENEID != "---", ]
map$GENEID <- trimws(map$GENEID)

## Remove ambiguous probes annotated to multiple genes ("///")
map <- map[!grepl("///", map$GENEID), ]

## Keep only valid IDs
if (id_label == "ENSG") {
  map <- map[grepl("^ENSG", map$GENEID), ]
} else {
  map <- map[grepl("^\\d+$", map$GENEID), ]
}
map <- unique(map)

## Collapse probe-level expression to gene-level expression (average across probes)
keep <- intersect(map$PROBEID, rownames(ex))
ex2 <- ex[keep, , drop = FALSE]
map2 <- map[match(rownames(ex2), map$PROBEID), ]
gene_ex <- avereps(ex2, ID = map2$GENEID)

## Export probe-level and gene-level expression matrices
write.table(
  data.frame(PROBEID = rownames(ex), ex, check.names = FALSE),
  "GSE15471_ex_probe.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(ID = rownames(gene_ex), gene_ex, check.names = FALSE),
  "GSE15471_ex_gene.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

## Extract expression values for target genes and report missing IDs
genes15 <- c(
  "54879","55876","2734","3624","23065","3651",
  "54465","10428","93210","90135","2972","3500"
)

expr15 <- gene_ex[intersect(genes15, rownames(gene_ex)), , drop = FALSE]
missing <- setdiff(genes15, rownames(expr15))
if (length(missing)) message("The following target IDs were not found in the platform annotation/expression matrix: ",
                             paste(missing, collapse = ", "))

## Sample grouping (GEO2R group membership string)
gsms <- paste0("00000000000000000000000000000000000000011111111111",
               "1111111111111111111111111111")
sml <- strsplit(gsms, split = "")[[1]]

if (length(sml) != ncol(ex)) stop("Length of gsms does not match the number of samples in the expression matrix.")

gs <- factor(sml)
levels(gs) <- make.names(c("normal", "tumor"))
group <- factor(gs, levels = c("normal", "tumor"))

## Define pairing (patient blocking factor) for paired analysis
## Assumes normals are listed first and tumors second, and pairs are matched by within-group order.
idx_n <- which(group == "normal")
idx_t <- which(group == "tumor")

if (length(idx_n) != length(idx_t)) stop("Unequal numbers of normal and tumor samples; cannot form complete pairs.")

npair <- length(idx_n)
patient <- rep(NA, length(group))
patient[idx_n] <- seq_len(npair)
patient[idx_t] <- seq_len(npair)
patient <- factor(patient)

## Sanity check: each patient should have both groups
pair_ok <- all(tapply(as.character(group), patient, function(x) length(unique(x)) == 2))
if (!pair_ok) warning("Pairing check failed: please verify that sample order correctly matches normal/tumor pairs.")

## Export GraphPad-friendly table (includes patient ID for paired plotting)
out_wide <- data.frame(Sample = colnames(expr15), Patient = patient, Group = group, t(expr15), check.names = FALSE)
write.csv(out_wide, "GSE15471_15genes_forGraphPad.csv", row.names = FALSE)

## Paired differential expression testing using limma (patient blocking)
design <- model.matrix(~ patient + group)
fit <- lmFit(expr15, design)
fit2 <- eBayes(fit)

coef_name <- grep("^group", colnames(design), value = TRUE)
coef_name <- coef_name[grep("tumor", coef_name, ignore.case = TRUE)][1]

res <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "none")
res$ID <- rownames(res)

## Export limma results (P values and BH-adjusted P values) for the selected genes
write.csv(res[, c("ID","logFC","AveExpr","t","P.Value","adj.P.Val")],
          "GSE15471_15genes_limma_Pvalues.csv", row.names = FALSE)
