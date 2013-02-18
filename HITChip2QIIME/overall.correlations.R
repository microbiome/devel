stats <- c()
for (output.dir.id in c("", "-quantiles")) {
 for (method in summarization.methods) {

  # Get Genus-level abundance matrix in QIIME format
  outputdir <- paste("qiime", output.dir.id, "-species-", method, sep = "")
  genus.hitchip <- get.qiime.matrix(outputdir, file.id = "abundance_with_taxa") 

  # Pick taxa that are available in both HITChip and NGS data
  common.taxa <- intersect(rownames(genus.hitchip), rownames(genus.ngs))
  common.samples <- intersect(colnames(genus.hitchip), colnames(genus.ngs))
  hit <- genus.hitchip[common.taxa, common.samples]; hit[hit < 1e-4] <- 0
  ngs <- genus.ngs[common.taxa, common.samples]

  # Ignore zero abundances and plot in log scale
  # Correlation with zeroes excluded and at log10
  nv <- log10(unlist(ngs)) 
  hv <- log10(unlist(hit))
  keep <- rowSums(is.infinite(cbind(nv, hv))) == 0
  stats[[outputdir]] <- cor(nv[keep], hv[keep])
 }
}

