# TCGA eQTL

## Input files (place in the corresponding folders and/or update `setwd()` paths in the R script)
**Part 1 (INT transform; optional)**  
- `tpm.csv` : expression matrix (rows = genes, columns = samples; header = sample IDs; first column = gene IDs)

**Part 2 (adjust for methylation + CNV)**  
- `fpkm_INT.txt` : INT-transformed expression (rows = genes, columns = samples; tab-delimited)  
- `methylation_matrix.txt` : methylation matrix aligned to the same genes/samples  
- `cnv_matrix.txt` : CNV matrix aligned to the same genes/samples  

**Part 3 (MatrixEQTL)**  
- `SNP.txt` : genotype dosage matrix (rows = SNPs, columns = samples; first column = SNP IDs; tab-delimited)  
- `snpsloc.txt` : SNP genomic positions (required by MatrixEQTL)  
- `GE.txt` : expression matrix used for eQTL (genes × samples; tab-delimited)  
- `geneloc.txt` : gene genomic positions (required by MatrixEQTL)  
- `Covariates.txt` : covariates matrix (covariates × samples; includes Age, Gender, 5 genotype PCs, 15 PEER factors, CNAs, methylation)

## Run
Open R (or RStudio), set working directories inside the script (the three `setwd("D:/file/...")` lines), then run:
- `source("run_tcga_eqtl.R")`

## Output
The script writes intermediate files (INT matrix, residual expression matrix) and MatrixEQTL results:
- `cis.csv` (cis-eQTLs)  
- `trans.csv` (trans-eQTLs)  
plus optional Shapiro test P-value tables and diagnostic plots.
