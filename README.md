# Code availability

This repository provides the custom scripts and command-line records used in the main computational analyses of our manuscript, including:
1) TCGA pancreatic eQTL mapping
2) eQTL meta-analysis
3) eGene downstream analyses (immune infiltration and drug sensitivity)
4) GO/KEGG enrichment and differential expression analyses
5) Circos plot preparation/visualization
6) GWAS association testing and GWAS meta-analysis (command-line workflows)

## Important note on data
No raw or controlled-access data are redistributed in this repository.
Please obtain TCGA/GTEx/public resource data through authorized channels as described in the manuscript Methods.

## Folder structure
- `scripts/01_tcga_eqtl/` : TCGA eQTL mapping scripts
- `scripts/02_eqtl_meta/` : eQTL meta-analysis scripts
- `scripts/03_immune_drug/` : immune infiltration & drug sensitivity analyses
- `scripts/04_enrichment_differentialexpression/` : GO/KEGG enrichment and differential expression
- `scripts/05_circos/` : Circos input preparation and plotting
- `scripts/06_gwas/` : GWAS association testing (PLINK) and GWAS meta-analysis (META v1.7) command-line workflows

## How this maps to figures/tables
- Fig 1 (eQTL results): `01_tcga_eqtl`, `02_eqtl_meta`
- Fig 3 (immune/drug analyses): `03_immune_drug`
- Fig S1 (GO/KEGG): `04_enrichment_differentialexpression`
- Fig 4 (differential expression): `04_enrichment_differentialexpression`
- Fig 4 (Circos): `05_circos`
- GWAS / GWAS meta (results used in the manuscript): `06_gwas`

## Software
- R scripts: R version and package details are documented in `sessionInfo_R.txt`.
- GWAS association testing: PLINK (commands provided under `scripts/06_gwas/`).
- GWAS meta-analysis: META v1.7 (commands provided under `scripts/06_gwas/`).
```
