# COMPARE TO NGS DATA

summarization.methods = c("sum", "ave", "rpa", "frpa") 
outputdirs <- c("qiime-species-sum", "qiime-species-ave",
                "qiime-species-rpa", "qiime-species-frpa",
  		"qiime-quantiles-species-sum", "qiime-quantiles-species-ave", 
 		"qiime-quantiles-species-rpa", "qiime-quantiles-species-frpa")

source("funcs.R")

# Read NGS data for the same samples
genus.ngs <- get.qiime.matrix("~/projects/hitchip-benchmarking13/data", file.id = "otu_table") 
# Rename columns
colnames(genus.ngs) <- gsub("\\.", "", gsub(".l3.2", "", gsub(".l3.1", "", colnames(genus.ngs))))

# ----------------------------------------------

source("overall.correlations.R")
pdf("OverallCorrelations.pdf")
par(mar = c(3, 14, 1, 1)); barplot(sort(stats), horiz = T, las = 1, main = "NGS-HITChip correlation")
dev.off()

# --------------------------------------------------

# Check RPA variants and ave?
source("sample.correlations.R")

# --------------------------------------------------------

# Problematic due to compositionality effects
#source("taxa.correlations.R")



