print(" Specify output files / dirs ")

  # Sample mapping file for QIIME
  sample.mapping.file <- "qiime_sample_mapping_file.txt"

print(" Input files / directories")
  # Precalculated HITChip taxonomy file
  hitchip.taxonomy.file <- "hitchip_taxonomic_assignments/HITChip16S_tax_assignments.txt"
  # HITChip probe-taxa phylogenies CSV file 
  phylogeny.file <- "HITChipPhylogeny.csv"

# -----------------------------------------------------------

print("Load functions")
source("funcs.R")

print("Convert HITChip matrix into QIIME format and write into text file")
sample.names <- HITChip2QIIME(data.dir = hitchip.data.dir, level = level, method = method, phylogeny.file = phylogeny.file, hitchip.output.file = hitchip.output.file)

# Write sample mapping file for QIIME
tmp <- write.QIIME.sample.mapping.file(sample.names, output.file = sample.mapping.file) 

print("Create BIOM abundance table")
system(paste("trflp_file_to_otu_table.py -i", hitchip.output.file, "-o abundance_otu.table.biom"))

print("Link observations and metadata")
system(paste("add_metadata.py -i abundance_otu.table.biom -o abundance_with_taxa.biom --sample_mapping_fp", sample.mapping.file, "--observation_mapping_fp", hitchip.taxonomy.file, "--observation_header OTUID,taxonomy,confidence --sc_separated taxonomy"))

print("Get the summaries from QIIME for HITChip")
system(paste("rm -rf", outputdir))
system(paste("summarize_taxa_through_plots.py -m", sample.mapping.file, "-i abundance_with_taxa.biom -o", outputdir))

# ------------------------------------------------------------
