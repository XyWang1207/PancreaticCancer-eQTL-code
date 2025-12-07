###############################################################################
# Circos plot (RCircos): integrate GWAS, eQTL, links, and candidate genes
#
# Input files (tab-delimited recommended):
#   - candidate_genes_15.txt     : candidate gene positions and labels for plotting
#   - circos_eQTL_FDRmin.txt     : eQTL points for the scatter track
#   - circos_GWAS.txt            : GWAS points for the scatter track
#   - circos_link.txt            : link/ribbon connections between loci and genes
#
# Output:
#   - Circos_plot_allGWAS_black.tif (600 dpi)
#
# Notes:
# - Input tables must follow RCircos-required formats for each plotting function.
# - This script uses the built-in UCSC HG19 ideogram in RCircos.
###############################################################################

rm(list = ls())

# Set your local working directory (contains the input files)
setwd("D:/file/circos_processing")

library(RCircos)
library(data.table)

# ---- Load input tables ----
gene <- fread("candidate_genes.txt")
eQTL <- fread("circos_eQTL.txt")
link <- fread("circos_link.txt")
GWAS <- fread("circos_GWAS.txt")

# ---- Load human ideogram (hg19) ----
data(UCSC.HG19.Human.CytoBandIdeogram)
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram

# Exclude sex chromosomes
chr.exclude <- c("chrX", "chrY")

# ---- Initialize RCircos core components ----
tracks.inside <- 8
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

# ---- Output device ----
tiff(
  file = "Circos_plot_allGWAS_black.tif",
  height = 120, width = 120, units = "mm",
  res = 600, compression = "lzw", pointsize = 5
)

# ---- Plot ideogram ----
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

# ---- Global track style ----
params <- RCircos.Get.Plot.Parameters()
params$track.height <- 0.25
params$track.background <- NULL
RCircos.Reset.Plot.Parameters(params)

# ---- Track 1: eQTL scatter ----
data.col <- 4
track.num <- 1
side <- "in"
RCircos.Scatter.Plot(eQTL, data.col, track.num, side)

# ---- Track 2: GWAS scatter ----
data.col <- 4
track.num <- 2
side <- "in"
RCircos.Scatter.Plot(GWAS, data.col, track.num, side)

# ---- Links / ribbons ----
params <- RCircos.Get.Plot.Parameters()
params$track.height <- 0.22
params$track.background <- NULL
RCircos.Reset.Plot.Parameters(params)

track.num <- 5
RCircos.Link.Plot(link, track.num, TRUE)
RCircos.Ribbon.Plot(ribbon.data = link, track.num = 5, by.chromosome = FALSE, twist = FALSE)

# ---- Gene connectors and labels ----
params <- RCircos.Get.Plot.Parameters()
params$text.size <- 0.6
params$track.height <- 0.11
RCircos.Reset.Plot.Parameters(params)

side <- "in"
track.num <- 5
RCircos.Gene.Connector.Plot(gene, track.num, side)

params <- RCircos.Get.Plot.Parameters()
params$track.height <- 0.2
params$track.background <- NULL
RCircos.Reset.Plot.Parameters(params)

name.col <- 4
track.num <- 4
RCircos.Gene.Name.Plot(gene, name.col, track.num, side)

dev.off()
rm(list = ls())

