###############################################################################
# TCGA eQTL pipeline (simple, copy-and-run version)
# Part 1: Inverse normal transformation (INT) of expression matrix
# Part 2: Adjust expression for methylation and CNV (linear regression residuals)
# Part 3: eQTL mapping using MatrixEQTL
#
# NOTE:
# - Please change the setwd() paths to where YOUR local files are stored.
# - Genotype PCs are assumed to be already included in Covariates.txt.
###############################################################################

require(data.table)
options(scipen = 200)

###############################################################################
# Part 1: Inverse normal transformation (INT) of expression matrix
###############################################################################
setwd("D:/file")

# Expression matrix: rows = genes, columns = samples
tpm <- read.csv("tpm.csv", header = TRUE, row.names = 1, check.names = FALSE)

# INT transformation gene-by-gene
y <- qnorm((rank(tpm[1, ], na.last = "keep") - 0.5) / sum(!is.na(tpm[1, ])))
for (i in 2:nrow(tpm)) {
  ys <- qnorm((rank(tpm[i, ], na.last = "keep") - 0.5) / sum(!is.na(tpm[i, ])))
  y <- rbind(y, ys)
}
write.csv(y, "tpm_INT.csv", row.names = TRUE)

# Shapiro test on INT-transformed matrix (P-value per gene; may be slow)
data_int <- read.csv("tpm_INT.csv", row.names = 1, check.names = FALSE)

p <- shapiro.test(as.numeric(data_int[1, ]))
for (i in 2:nrow(data_int)) {
  ps <- shapiro.test(as.numeric(data_int[i, ]))
  p <- rbind(p, ps)
}
write.csv(p, "tpm_INT_shapiroP.csv", row.names = TRUE)

###############################################################################
# Part 2: Adjust expression for methylation and CNV
###############################################################################
# Transpose to samples x genes, and export intermediate files
# Input matrices are expected to have rows = genes, columns = samples.

mRNA  <- read.table("fpkm_INT.txt", row.names = 1, check.names = FALSE, header = TRUE,
                    sep = "\t", quote = "")
methy <- read.table("methylation_matrix.txt", row.names = 1, check.names = FALSE, header = TRUE,
                    sep = "\t", quote = "")
cnv   <- read.table("cnv_matrix.txt", row.names = 1, check.names = FALSE, header = TRUE,
                    sep = "\t", quote = "")

mRNA[is.na(mRNA)]   <- 0
methy[is.na(methy)] <- 0
cnv[is.na(cnv)]     <- 0

mRNA  <- data.frame(t(mRNA),  check.names = FALSE)  # samples x genes
methy <- data.frame(t(methy), check.names = FALSE)  # samples x genes
cnv   <- data.frame(t(cnv),   check.names = FALSE)  # samples x genes

write.csv(mRNA,  "fpkm_INT_samples_x_genes.csv", row.names = TRUE)
write.csv(methy, "methylation_samples_x_genes.csv", row.names = TRUE)
write.csv(cnv,   "cnv_samples_x_genes.csv", row.names = TRUE)

# Linear regression residuals: expression ~ methylation + CNV
# Output residuals matrix is the adjusted expression used for eQTL.
mRNA  <- read.csv("fpkm_INT_samples_x_genes.csv", row.names = 1, check.names = FALSE)
methy <- read.csv("methylation_samples_x_genes.csv", row.names = 1, check.names = FALSE)
cnv   <- read.csv("cnv_samples_x_genes.csv", row.names = 1, check.names = FALSE)

lm.re <- lm(as.vector(as.matrix(mRNA[, 1])) ~ as.vector(as.matrix(methy[, 1])) + as.vector(as.matrix(cnv[, 1])))
e <- data.frame(resid(lm.re))

for (i in 2:ncol(mRNA)) {
  lm.re <- lm(as.vector(as.matrix(mRNA[, i])) ~ as.vector(as.matrix(methy[, i])) + as.vector(as.matrix(cnv[, i])))
  es <- data.frame(resid(lm.re))
  e <- cbind(e, es)
}

colnames(e) <- colnames(mRNA)
rownames(e) <- rownames(mRNA)

# Save residuals (samples x genes) and transpose to genes x samples
write.csv(e, "residual_exp_matrix.csv", row.names = TRUE)

e <- read.csv("residual_exp_matrix.csv", row.names = 1, check.names = FALSE)
e_t <- data.frame(t(e), check.names = FALSE)  # genes x samples

# Ensure sample IDs match the original expression columns (from fpkm_INT.txt)
col_ref <- read.table("fpkm_INT.txt", row.names = 1, check.names = FALSE, header = TRUE,
                      sep = "\t", quote = "")
colnames(e_t) <- colnames(col_ref)

write.csv(e_t,  "residual_exp_matrix_genes_x_samples.csv", row.names = TRUE)
write.table(e_t, "residual_exp_matrix_genes_x_samples.txt", row.names = TRUE,
            sep = "\t", quote = FALSE)

# Shapiro test for residual expression (P-value per gene; may be slow)
p <- shapiro.test(as.numeric(e_t[1, ]))
for (i in 2:nrow(e_t)) {
  ps <- shapiro.test(as.numeric(e_t[i, ]))
  p <- rbind(p, ps)
}
write.csv(p, "residual_shapiroP.csv", row.names = TRUE)

###############################################################################
# Part 3: eQTL mapping using MatrixEQTL
###############################################################################
library(MatrixEQTL)
require(data.table)

useModel <- modelLINEAR

SNP_file_name <- "SNP.txt"
snps_location_file_name <- "snpsloc.txt"

expression_file_name <- "GE.txt"
gene_location_file_name <- "geneloc.txt"

covariates_file_name <- "Covariates.txt"

output_file_name_cis <- tempfile()
output_file_name_tra <- tempfile()

pvOutputThreshold_cis <- 5e-2
pvOutputThreshold_tra <- 1e-20

errorCovariance <- numeric()
cisDist <- 1e6

snps <- SlicedData$new()
snps$fileDelimiter <- "\t"
snps$fileOmitCharacters <- "NA"
snps$fileSkipRows <- 1
snps$fileSkipColumns <- 1
snps$fileSliceSize <- 200000
snps$LoadFile(SNP_file_name)

gene <- SlicedData$new()
gene$fileDelimiter <- "\t"
gene$fileOmitCharacters <- "NA"
gene$fileSkipRows <- 1
gene$fileSkipColumns <- 1
gene$fileSliceSize <- 2000
gene$LoadFile(expression_file_name)

cvrt <- SlicedData$new()
cvrt$fileDelimiter <- "\t"
cvrt$fileOmitCharacters <- "NA"
cvrt$fileSkipRows <- 1
cvrt$fileSkipColumns <- 1
if (length(covariates_file_name) > 0) {
  cvrt$LoadFile(covariates_file_name)
}

snpspos <- read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos <- read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)

me <- Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)

unlink(output_file_name_tra)
unlink(output_file_name_cis)

cat("Analysis done in:", me$time.in.sec, "seconds", "\n")
cat("Detected local (cis) eQTLs:", "\n")
show(me$cis$eqtls)
cat("Detected distant (trans) eQTLs:", "\n")
show(me$trans$eqtls)

plot(me)
write.csv(me$cis$eqtls,  "cis.csv",  row.names = FALSE)
write.csv(me$trans$eqtls, "trans.csv", row.names = FALSE)
