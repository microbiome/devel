
  tab <- read.table(f, header = TRUE, sep = "\t", fill=T, stringsAsFactors=F)

  # Ensure the sample names contain only alphanumeric characters ('.' also OK), 
  # as required by QIIME
  colnames(tab) <- gsub("_", ".", colnames(tab))
  colnames(tab) <- gsub("\\.", "", colnames(tab))

  # Read mapping between accessions and names and replace species names with accessions
  phylogeny <- read.csv(phylogeny.file, header=T, fill=T, colClasses='character')

  # Level correspondence in HITChip vs. Phylogeny table
  if (level == "L1") {conversion.level <- "Phylum"}
  if (level == "L2") {conversion.level <- "Genus"}
  if (level == "species") {conversion.level <- "Phylotype"}

  # Get sequence accession for each HITChip taxa
  # Replace HITChip taxa by accessions
  tab[,1] <- phylogeny[match(tab[,1], phylogeny[,conversion.level]), "Accession"]

  # Round the observations in HITChip matrix into integers
  tab[,-1] <- round(tab[,-1])
  
  # Add column names directly to the table
  tab <- rbind(colnames(tab), tab)
      
  # Make samples x accessions
  tab <- t(tab)
  rownames(tab) <- NULL
  colnames(tab) <- NULL
  tab[1,1] <- ""


  # Write the matrix in output file in QIIME format 
  write.table(tab, sep='\t', file = hitchip.output.file, col.names=F, row.names = F, quote=F)

