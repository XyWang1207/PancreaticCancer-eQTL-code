# eQTL meta-analysis (GTEx + TCGA)

This script performs a simple pairwise meta-analysis for eQTL effects using `meta::metagen`, assuming each eQTL appears as two consecutive rows (e.g., GTEx and TCGA).

---

## Input
Place the input file in the working directory (or edit `setwd()` in the script):

- `meta_data.csv` (CSV with header)

Required columns:
- `id` : identifier/label for each row (e.g., `rsID_gene_GTEx`, `rsID_gene_TCGA`)
- `beta` : effect size (TE)
- `se` : standard error of beta (seTE)
- `pval` : single-study p-value (used only as input metadata; meta p-values are computed by `metagen`)

**Important assumption:** rows are ordered in pairs:
- row 1–2 = eQTL 1 (study A, study B)
- row 3–4 = eQTL 2 (study A, study B)
- ...

If your file is not ordered in pairs, please sort/group by eQTL before running.

---

## Run
1) Edit the path in the script if needed:
```r
setwd("D:/file/meta_eqtl")
2) Run in R/RStudio:
source("eqtl_meta.R")
