get.qiime.matrix <- function (outputdir, file.id) {

  genus.file <- paste(outputdir, "/", file.id, "_L6.txt", sep = "")
  tab <- read.table(genus.file, header = TRUE, row.names = 1)

  # Extract genus names
  nams <- unname(sapply(rownames(tab), function (x) {levs <- strsplit(x, ";")[[1]]; ind <- grep("g__", levs); if (length(ind) > 0) {gsub("g__", "", levs[[ind]])} else {NA}}))
  nams[is.na(nams)] <- ""

  # Remove entries with no genus-level classification
  keep <- which(!nams == "")
  tab <- tab[keep, ]
  nams <- nams[keep]

  # Remove duplicated occurrences of a genus
  keep <- which(!duplicated(nams))
  tab <- tab[keep, ]
  nams <- nams[keep]   

  rownames(tab) <- nams

  tab

}



#' @param hitchip.matrix HITChip data matrix taxa x samples
HITChip2QIIME <- function (data.dir, level = "species", method = "rpa", phylogeny.file, hitchip.output.file = "hitchip_data.txt") {

  # level <- "oligo"; method = "sum"; data.dir = "test/"; log10 = TRUE
  if (level %in% c("L0", "L1", "L2", "species")) {
      f <- paste(data.dir, "/", level, "-", method, ".tab", sep = "")
  } else if (level == "oligo") {
      f <- paste(data.dir, "/oligoprofile.tab", sep = "")
  }
  
  #f <- "hitchip.files/speciesSumModified.tab"
  message(paste("Reading", f))

  # table from HITChip output: sample names by species accessions
  hit.data <- read.table(f, header=T, sep='\t', fill=T, stringsAsFactors=F)
  colnames(hit.data) <- gsub("_", ".", colnames(hit.data))

  # read mapping between accessions and names and replace species names with accessions
  phylogeny <- read.csv(phylogeny.file, header=T, fill=T, colClasses='character')

  # Level correspondence in HITChip vs. Phylogeny table
  if (level == "L1") {conversion.level <- "Phylum"}
  if (level == "L2") {conversion.level <- "Genus"}
  if (level == "species") {conversion.level <- "Phylotype"}

  # Get sequence accession for each HITChip taxa
  # Replace HITChip taxa by accessions
  hit.data[,1] <- phylogeny[match(hit.data[,1], phylogeny[,conversion.level]), "Accession"]

  hit.abund <- t(cbind(hit.data[,1], round((hit.data[,2:length(hit.data[1,])]))))
  sample.names <- colnames(hit.data)
  sample.names[1] <- ''

  # Ensure the sample names contain only alphanumeric characters ('.' also OK), 
  # as required by QIIME
  sample.names <- gsub("_", ".", sample.names)
  sample.names <- gsub("\\.", "", sample.names)
  rownames(hit.abund) <- sample.names

  write.table(hit.abund, sep='\t', file = hitchip.output.file, col.names = F, quote=F)

  sample.names[-1]

}

#' @param hitchip.matrix HITChip data matrix taxa x samples
HITChip2QIIME006 <- function (data.dir, level = "species", method = "sum", phylogeny.file, hitchip.output.file = "hitchip_data.txt") {

  # data.dir = hitchip.data.dir; 

  # level <- "oligo"; method = "sum"; data.dir = "test/"; log10 = TRUE
  if (level %in% c("L0", "L1", "L2", "species")) {
      f <- paste(data.dir, "/", level, "_", method, "_006.tab", sep = "")
  } else if (level == "oligo") {
      f <- paste(data.dir, "/oligoprofile.tab", sep = "")
  }
  
  #f <- "hitchip.files/speciesSumModified.tab"
  message(paste("Reading", f))

  # table from HITChip output: sample names by species accessions
  hit.data <- t(read.table(f, header=T, sep='\t', fill=T, stringsAsFactors=F))
  colnames(hit.data) <- gsub("_", ".", colnames(hit.data))

  # read mapping between accessions and names and replace species names with accessions
  phylogeny <- read.csv(phylogeny.file, header=T, fill=T, colClasses='character')

  # Level correspondence in HITChip vs. Phylogeny table
  if (level == "level 2") {conversion.level <- "Genus"}
  if (level == "species") {conversion.level <- "Phylotype"}

  # Get sequence accession for each HITChip taxa
  # Replace HITChip taxa by accessions
  hit.data <- cbind(rownames(hit.data), hit.data) #phylogeny[match(hit.data[,1], phylogeny[,conversion.level]), "Accession"]
  rownames(hit.data) <- NULL

  sample.names <- hit.data[1,]
  sample.names[1] <- ''
  hit.abund <- t(cbind(hit.data[-1,1], hit.data[-1,-1]))

  # Ensure the sample names contain only alphanumeric characters ('.' also OK), 
  # as required by QIIME
  sample.names <- gsub("_", ".", sample.names)
  sample.names <- gsub("\\.", "", sample.names)
  rownames(hit.abund) <- sample.names

  write.table(hit.abund, sep='\t', file = hitchip.output.file, col.names = F, quote=F)

  sample.names[-1]

}


# Generate sample mapping file for QIIME
write.QIIME.sample.mapping.file <- function (sample.names, output.file = "grouping_for_qiime.txt", Nbarcodes = 200) {

  # generate pseudo-barcodes
  bcodes <- unique(sapply(1:Nbarcodes,FUN=function(x) paste(sample(c('A','C','T','G'),5,replace=T),collapse='')))

  # Write sample mapping to output file
  message(paste("Writing sample mappings to", output.file))
  cat('#SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'Treatment', 'Description', sep='\t', file = output.file)
  cat('\n', sep='', file = output.file, append = T)
  for(nn in 1:length(sample.names)) {
    cat(sample.names[nn], bcodes[nn], 'GGT', 'unknown', sample.names[nn],
      sep='\t', file = output.file, append=T)
    cat('\n', file = output.file, append=T)
  }

  output.file

}


