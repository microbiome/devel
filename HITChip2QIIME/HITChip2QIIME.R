## ----------------------------------------------------------------
## Convert HITChip profiling output to QIIME format
## ----------------------------------------------------------------

## Instructions
# 1. Specify the input and output files and parameters below
# 2. Run source('HITChip2QIIME.R')

# Parameters
level <- "species"

for (method in c("sum", "rpa", "frpa", "ave", "nmf")) {
  source("generate.qiime.files.R")
}

# source("compare.hitchip.ngs.R")