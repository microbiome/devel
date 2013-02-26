## ----------------------------------------------------------------
## Convert HITChip profiling output to QIIME format
## ----------------------------------------------------------------

## Instructions
# 1. Specify the input and output files and parameters below
# 2. Run source('HITChip2QIIME.R')

# Parameters
level <- "species"

outputdirs <- c()

summarization.methods = c("sum", "ave", "rpa", "frpa") 

hitchip.data.dir <- "hitchip.files" # HITChip profiling output dir
for (method in summarization.methods) {
  outputdir <- paste("qiime", level, method, sep = "-") # QIIME results
  # HITChip in QIIME format output file
  hitchip.output.file <- paste("hitchip_", level, "_", method, ".txt", sep = "")

  source("generate.qiime.files.R")
  outputdirs[[method]] <- outputdir
}

hitchip.data.dir <- "hitchip.files.quantiles" # HITChip profiling output dir
summarization.methods = c("sum", "ave", "rpa", "frpa") 
for (method in summarization.methods) {
  outputdir <- paste("qiime-quantiles", level, method, sep = "-") # QIIME results
  # HITChip in QIIME format output file
  hitchip.output.file <- paste("hitchip_", level, "_", method, ".txt", sep = "")
  source("generate.qiime.files.R")
  outputdirs[[paste(method, "quantiles", sep = "-")]] <- outputdir
}

# -------------------------------

# For 006 script:
#hitchip.data.dir <- "FECAL006" # HITChip profiling output dir
#for (method in c("ave", "sum")) {
#  outputdir <- paste("qiime-006", level, method, sep = "-") # QIIME results
#  # HITChip in QIIME format output file
#  hitchip.output.file <- paste("hitchip_", level, "_", method, ".txt", sep = "")
#  source("generate.qiime.files.R")
#  outputdirs[[paste(method, "006", sep = "-")]] <- outputdir
#}


# --------------------------------


#hitchip.data.dir <- "hitchip.files.specificoligos" # HITChip profiling output dir
#summarization.methods = c("sum", "ave", "rpa", "frpa") 
#for (method in summarization.methods) {
#  outputdir <- paste("qiime-specificoligos", level, method, sep = "-") # QIIME results
#  source("generate.qiime.files.R")
#  outputdirs[[paste(method, "specificoligos", sep = "")]] <- outputdir
#}

# source("compare.hitchip.ngs.R")



