# 04_enrichment_differentialexpression

This folder contains two analysis scripts used for (i) functional enrichment of prioritized eGenes and (ii) differential expression analyses in public GEO datasets (paired and non-paired designs).

> Note: This repository does **not** redistribute raw TCGA/GEO-derived data.  
> GEO expression data are downloaded directly from NCBI GEO by `GEOquery` when the GEO scripts are run.

---

## Contents

- `run_enrichment_GO_KEGG.R`  
  GO and KEGG enrichment analysis for the eGene list (clusterProfiler).

- `run_GEO_DE_limma.R`  
  Differential expression of selected genes in:
  - **GSE62165** (non-paired tumor vs normal)
  - **GSE15471** (paired tumor vs normal with patient blocking)

---

## A) GO / KEGG enrichment

### Input
Place the following file in the working directory (or edit `setwd()` inside the script):

- `eGenes.txt`  
  Must contain a column named `gene_id` with **human ENSEMBL IDs** (e.g., `ENSG...`).

### Run
In R/RStudio:
```r
source("run_enrichment_GO_KEGG.R")

````

### Output

* `eGenes_GO.txt` : GO enrichment results (ALL ontologies; BH correction)
* `GO_enrichment_barplot.pdf` : GO barplot (top terms per ontology)
* `eGenes_KEGG_BH.txt` : KEGG enrichment results (BH correction; human `hsa`)
* `KEGG_enrichment_barplot.pdf` : KEGG barplot (top 20 pathways)

### Dependencies

```r
install.packages(c("data.table","dplyr","ggplot2","tools"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler","org.Hs.eg.db"))
```

---

## B) GEO differential expression (limma)

This script follows a GEO2R-style workflow:
(i) download expression matrices via `GEOquery`, (ii) log2-transform if needed, (iii) map probes to gene IDs and collapse to gene-level expression (`limma::avereps`), and (iv) test differential expression using `limma`.

### Datasets

* **GSE62165**: non-paired tumor vs normal
* **GSE15471**: paired tumor vs normal (patient blocking)

### Run

Edit `setwd()` to point to your local folders (recommended: one folder per dataset), then run:

```r
source("run_GEO_DE_limma.R")
```

### Output (per dataset)

* `*_ex_probe.txt` : probe-level expression matrix
* `*_ex_gene.txt` : gene-level expression matrix (after collapsing probes)
* `*_forGraphPad.csv` : GraphPad-friendly table for the selected genes
* `*_limma_Pvalues.csv` : limma results (P values and BH-adjusted P values) for the selected genes

### Dependencies

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery","limma"))
```

```
```
