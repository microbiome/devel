## ----------------------------------------------------------------
## Convert HITChip profiling output to QIIME format
## ----------------------------------------------------------------

## Instructions
# 1. Specify the input and output files and parameters below
# 2. Run source('HITChip2QIIME.R')

# Parameters
level <- "species"


# Add(?) frpa.with.affinities, frpa.with.affinities.direct, 006, atlas3268, atlas5000
outputdirs <- c()
summarization.methods = c("sum", "ave", "rpa", "frpa", "sum.through.species", "ave.through.species", "rpa.direct", "rpa.with.affinities", "rpa.with.affinities.direct") 
for (method in summarization.methods) {
  source("generate.qiime.files.R")
  outputdirs[[method]] <- outputdir
}



# source("compare.hitchip.ngs.R")



