## ----------------------------------------------------------------
## Convert HITChip profiling output to QIIME format
## ----------------------------------------------------------------

## Instructions
# 1. Specify the input and output files and parameters below
# 2. Run source('HITChip2QIIME.R')

# Parameters
level <- "species"

# Add NMF and 006, frpa with affinities etc., direct RPA etc.

outputdirs <- c()
for (method in c("sum", "rpa", "frpa", "ave")) {
  source("generate.qiime.files.R")
  outputdirs[[method]] <- outputdir
}

# source("compare.hitchip.ngs.R")