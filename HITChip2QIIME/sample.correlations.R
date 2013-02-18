samplecorrelations <- matrix(NA, nrow = length(common.samples), ncol = length(summarization.methods))
rownames(samplecorrelations) <- common.samples
colnames(samplecorrelations) <- summarization.methods
for (method in setdiff(summarization.methods, "ave")) {
  
  # Get Genus-level abundance matrix in QIIME format
  outputdir <- paste("qiime-species", method, sep = "-")
  genus.hitchip <- get.qiime.matrix(outputdir, file.id = "abundance_with_taxa") 
  source("benchmarks.R")		  		
  samplecorrelations[common.samples, method] <- samplecor[common.samples]

}
dev.off()
pdf("samplecorrelations-histograms.pdf")
par(mfrow = c(3,3))
for (method in setdiff(summarization.methods, "ave")) {
  hist(samplecorrelations[common.samples, method], main = method, xlim = c(-1,1))
}
dev.off()

pdf("samplecorrelations-boxplot.pdf")
o <- order(apply(samplecorrelations, 2, median))
par(mar = c(10, 3, 1,1)); boxplot(samplecorrelations[, o], las = 2)
dev.off()
#plot(samplecorrelations[, "frpa"], samplecorrelations[, "sum"]); abline(0,1)
#par(mar = c(2, 10, 1, 1)); barplot(sort(apply(samplecorrelations, 2, function (x) {mean(na.omit(x))})), las = 1, horiz = T)
#par(mar = c(2, 10, 1, 1)); barplot(sort(apply(samplecorrelations, 2, function (x) {median(na.omit(x))})), las = 1, horiz = T, main = "Median sample correlations NGS-HITChip")