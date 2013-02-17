## ----------------------------------------------------------------
## Convert HITChip profiling output to QIIME format
## ----------------------------------------------------------------

## Instructions
# 1. Specify the input and output files and parameters below
# 2. Run source('HITChip2QIIME.R')

# Parameters
level <- "species"

# Add NMF and 006, frpa with affinities etc.

for (method in c("sum", "rpa", "frpa", "ave")) {
  source("generate.qiime.files.R")
}

# source("compare.hitchip.ngs.R")