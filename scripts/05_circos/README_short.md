# Circos integration plot (RCircos)

This script generates a Circos-style integration plot using **RCircos**, combining:
- GWAS association signals (scatter track)
- eQTL signals (scatter track)
- locus-to-gene links (link/ribbon track)
- candidate gene labels (connector + gene name track)

The script uses the built-in **UCSC hg19 cytoband ideogram** provided by RCircos.

---

## Files

### Script
- `circos.R`

### Required input files
Place the following input files in the working directory (or edit `setwd()` in the script):

- `candidate_genes.txt`  
  Candidate gene positions and labels for plotting (RCircos gene plot format).

- `circos_eQTL.txt`  
  eQTL points for the scatter track (RCircos scatter plot format).

- `circos_GWAS.txt`  
  GWAS points for the scatter track (RCircos scatter plot format).

- `circos_link.txt`  
  Links/ribbons connecting loci and genes (RCircos link/ribbon plot format).

> Note: The exact column structure must follow RCircos input specifications for each plot type.

---

## Run

1) Edit the working directory inside the script if needed:
```r
setwd("D:/file/circos_processing")
````

2) Run in R/RStudio:
```r
source("circos.R")
```

---

## Output

* `Circos_plot_allGWAS_black.tif`
  High-resolution TIFF (600 dpi; LZW compression), suitable for manuscript figures.

---

## Dependencies

R packages:

* `RCircos`
* `data.table`

Install if needed:

```r
install.packages("data.table")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("RCircos")
```

---

## Notes

* Sex chromosomes (`chrX`, `chrY`) are excluded by default.
* The plot is generated using the **hg19** ideogram (`UCSC.HG19.Human.CytoBandIdeogram`) bundled with RCircos.
* Track positions/heights are set in the script and can be adjusted via `RCircos.Get.Plot.Parameters()` and `RCircos.Reset.Plot.Parameters()`.

```
```

