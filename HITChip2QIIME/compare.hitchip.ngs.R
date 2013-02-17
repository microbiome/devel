# COMPARE TO NGS DATA

source("funcs.R")

# Read NGS data for the same samples
genus.ngs <- get.qiime.matrix("~/projects/hitchip-benchmarking13/data", file.id = "otu_table") 
# Rename columns
colnames(genus.ngs) <- gsub("\\.", "", gsub(".l3.2", "", gsub(".l3.1", "", colnames(genus.ngs))))

# ----------------------------------------------

outputdirs <- c(ave = "qiime-species-ave", sum = "qiime-species-sum", rpa = "qiime-species-rpa")

par(mfrow = c(2, 2))
stats <- c()
for (nam in names(outputdirs)) {
  # Get Genus-level abundance matrix in QIIME format
  outputdir <- outputdirs[[nam]]
  genus.hitchip <- get.qiime.matrix(outputdir, file.id = "abundance_with_taxa") 
  source("benchmarks.R")		  		
  stats[[nam]] <- nonzero.corr.log10
}







