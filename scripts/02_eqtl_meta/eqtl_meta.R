###############################################################################
# eQTL meta-analysis (simple, pairwise fixed/random effects using meta::metagen)
#
# Input:
#   meta_data.csv (tabular file with at least: id, beta, se, pval)
#   Each eQTL appears in TWO rows (e.g., GTEx and TCGA) in consecutive order:
#     row 1 = study A
#     row 2 = study B
#     row 3 = study A
#     row 4 = study B
#     ...
#
# Output:
#   meta_result.csv
#   Columns: id, p_heterogeneity_Q, p_fixed, p_random
#
# Notes:
# - This script assumes consecutive paired rows per id (minimal-change version).
# - If your data are not ordered in pairs, sort/group by id before running.
###############################################################################

require(data.table)
library(meta)

setwd("D:/file/meta_eqtl")

# Read input (CSV). If yours is TSV, change fread() arguments accordingly.
data <- fread("meta_data.csv")

# Basic checks
required_cols <- c("id", "beta", "se", "pval")
missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols) > 0) {
  stop("Missing required columns in meta_data.csv: ", paste(missing_cols, collapse = ", "))
}

if (nrow(data) < 2) stop("meta_data.csv has <2 rows; cannot run pairwise meta-analysis.")

if (nrow(data) %% 2 == 1) {
  warning("Number of rows is odd. The last row will be ignored to keep pairwise meta-analysis.")
}

result <- data.frame(
  id = character(),
  pval_Q = numeric(),
  pval_fixed = numeric(),
  pval_random = numeric(),
  stringsAsFactors = FALSE
)

# Pairwise meta-analysis: rows (1,2), (3,4), ...
for (i in seq(1, nrow(data) - 1, by = 2)) {
  
  dat2 <- data[i:(i + 1), ]
  
  # Use the first row's id as the identifier (assumes the pair shares the same id)
  pair_id <- as.character(dat2$id[1])
  
  # Run meta-analysis (beta +/- SE)
  m <- metagen(
    TE = dat2$beta,
    seTE = dat2$se,
    studlab = dat2$id,      # study label; if you have a separate study column, use that instead
    data = dat2,
    sm = "SMD"              # effect scale label; not used for inference here, but required by metagen
  )
  
  # Extract P-values:
  # - m$pval.Q: Cochran's Q test for heterogeneity
  # - m$pval.fixed: fixed-effect meta-analysis P
  # - m$pval.random: random-effects meta-analysis P
  result <- rbind(result, data.frame(
    id = pair_id,
    I2 = as.numeric(m$I2), # I-squared (%)
    pval_Q = m$pval.Q,
    pval_fixed = m$pval.fixed,
    pval_random = m$pval.random,
    stringsAsFactors = FALSE
  ))
}

write.csv(result, "meta_result.csv", row.names = FALSE)
