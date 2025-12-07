# eGene annotation (immune infiltration / drug association / MAF summary)

This script performs downstream annotation for prioritized eGenes, including:
1) Spearman correlation between eGene expression and immune infiltration scores
2) BH FDR correction (within each immune feature) and summary outputs
3) Drug–gene association table joins and optional dot-plot table export
4) MAF summary and oncoplot (maftools)

---

## Input files (prepare locally; no raw TCGA data are distributed)

### A) Immune correlation (`D:/file/immune`)
- `immune_scores.txt`  
  Immune infiltration matrix **(samples × features)**.  
  The **first column** is sample ID; remaining columns are immune features.
- `expression_matrix.txt`  
  Expression matrix **(genes × samples)**.  
  The **first column** is `gene_name`; remaining columns are sample IDs (same samples as above).
- `eGenes.txt`  
  eGene list/table. Must contain a column named `gene_name`.

### B) Drug association (`D:/file/drug`)
- `drug_pair.txt`  
  Mapping table used to link genes to drugs. Must include column `RNA.molecule` and a drug column that matches `Compound`.
- `eGenes.txt`  
  Same as above. If the 2nd column is the gene identifier used for drug pairing, it will be renamed to `RNA.molecule` by the script.
- `PAAD_ANOVA.csv`  
  Drug response/statistics table exported from your portal/tool.  
  The script keeps columns 1,4,5,6 and renames column 1 to `Compound`.

### C) MAF / mutation summary (`D:/file/tmb_maf`)
- `PAAD_maftools.maf`  
  TCGA PAAD MAF file (or a MAF formatted for maftools).
- `top20_genes.txt`  
  A plain text file with one gene symbol per line (no header).

---

## Run
1) Edit the three `setwd("D:/file/...")` lines in the script to match your local folders.
2) In R/RStudio:
```r
source("run_egene_annotation.R")
