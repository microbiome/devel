# COMPARE TO NGS DATA

source("funcs.R")

# Read NGS data for the same samples
genus.ngs <- get.qiime.matrix("~/projects/hitchip-benchmarking13/data", file.id = "otu_table") 
# Rename columns
colnames(genus.ngs) <- gsub("\\.", "", gsub(".l3.2", "", gsub(".l3.1", "", colnames(genus.ngs))))

# ----------------------------------------------

par(mfrow = c(3, 3))
stats <- c()
for (method in setdiff(summarization.methods, "ave")) {
  
  # Get Genus-level abundance matrix in QIIME format
  outputdir <- paste("qiime-species", method, sep = "-")
  genus.hitchip <- get.qiime.matrix(outputdir, file.id = "abundance_with_taxa") 
  source("benchmarks.R")		  		
  stats[[method]] <- nonzero.corr.log10

}
dev.off()
pdf("OverallCorrelations.pdf")
par(mar = c(1, 10, 1, 1)); barplot(sort(stats), horiz = T, las = 1, main = "NGS-HITChip correlation")
dev.off()

# --------------------------------------------------

# Problematic due to compositionality effects
#source("taxa.correlations.R")

# --------------------------------------------------------

# Check RPA variants and ave?
source("sample.correlations.R")
