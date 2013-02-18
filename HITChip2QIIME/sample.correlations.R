

samplecorrelations <- matrix(NA, nrow = length(common.samples), ncol = length(outputdirs))
rownames(samplecorrelations) <- common.samples
colnames(samplecorrelations) <- outputdirs
for (outputdir in outputdirs) {
  
  # Get Genus-level abundance matrix in QIIME format
  genus.hitchip <- get.qiime.matrix(outputdir, file.id = "abundance_with_taxa") 

  # Pick taxa that are available in both HITChip and NGS data
  common.taxa <- intersect(rownames(genus.hitchip), rownames(genus.ngs))
  common.samples <- intersect(colnames(genus.hitchip), colnames(genus.ngs))
  hit <- genus.hitchip[common.taxa, common.samples]; hit[hit < 1e-4] <- 0
  ngs <- genus.ngs[common.taxa, common.samples]

  samplecor <- c()
  for (pt in common.samples) {
    a <- log10(unlist(hit[, pt]))
    b <- log10(unlist(ngs[,pt]))
    keep <- which(rowSums(is.infinite(cbind(a, b))) == 0)
    samplecor[[pt]] <- cor(a[keep], b[keep], method = "pearson")
  }

  samplecorrelations[common.samples, outputdir] <- samplecor[common.samples]

}

pdf("samplecorrelations-histograms.pdf")
par(mfrow = c(3,3))
for (method in colnames(samplecorrelations)) {
  hist(samplecorrelations[common.samples, method], main = method, xlim = c(-1,1))
}
dev.off()

pdf("samplecorrelations-boxplot.pdf")
o <- order(apply(samplecorrelations, 2, median))
par(mar = c(14, 3, 3,1)); boxplot(samplecorrelations[, o], las = 2, main = "NGS-HITChip sample correlations")
dev.off()

#####################################################
#plot(samplecorrelations[, "frpa"], samplecorrelations[, "sum"]); abline(0,1)
#par(mar = c(2, 10, 1, 1)); barplot(sort(apply(samplecorrelations, 2, function (x) {mean(na.omit(x))})), las = 1, horiz = T)
#par(mar = c(2, 10, 1, 1)); barplot(sort(apply(samplecorrelations, 2, function (x) {median(na.omit(x))})), las = 1, horiz = T, main = "Median sample correlations NGS-HITChip")