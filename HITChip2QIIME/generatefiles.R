## ----------------------------------------------------------------
## Convert HITChip profiling output to QIIME format
## ----------------------------------------------------------------

## Instructions
# 1. Specify the input and output files and parameters below
# 2. Run source('HITChip2QIIME.R')

# Parameters
level <- "species"
method <- "sum"

# Output files / dirs 
outputdir <- paste("qiime", level, method, sep = "-")
hitchip.output.file <- paste("hitchip_", level, "_", method, ".txt", sep = "")

# Input files / directories
sample.mapping.file <- "qiime_sample_mapping_file.txt"
phylogeny.file <- "HITChipPhylogeny.csv"
hitchip.taxonomy.file <- "hitchip_taxonomic_assignments/HITChip16S_tax_assignments.txt"
hitchip.data.dir <- "hitchip.files"

# -----------------------------------------------------------

# Load functions
source("funcs.R")

# Convert HITChip matrix into QIIME format and write to text file
sample.names <- HITChip2QIIME(data.dir = hitchip.data.dir, level = level, method = method, phylogeny.file = phylogeny.file, hitchip.output.file = hitchip.output.file)

# Write sample mapping file for QIIME
tmp <- write.QIIME.sample.mapping.file(sample.names, output.file = sample.mapping.file) 

