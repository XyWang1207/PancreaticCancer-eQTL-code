# GWAS association testing (PLINK)
# Note: Replace paths with your local paths and provide your own genotype files.

# 1) QC (example)
plink --bfile <input_genotypes> \
  --geno 0.05 --mind 0.02 --maf 0.01 --hwe 1e-5 \
  --make-bed --out <qc_out>

# 2) PCA (example; adapt to your workflow)
plink --bfile <qc_out> --pca 20 --out <pca_out>

# 3) Logistic regression with covariates (example)
plink --bfile <qc_out> \
  --logistic hide-covar \
  --covar <covariates.txt> \
  --covar-name Age,Sex,PC1,PC2,PC3,PC4,PC5 \
  --pheno <pheno.txt> --pheno-name case_control \
  --out <gwas_out>
