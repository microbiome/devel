taxa.correlations <- matrix(NA, nrow = length(common.taxa), ncol = length(summarization.methods))
rownames(taxa.correlations) <- common.taxa
colnames(taxa.correlations) <- summarization.methods
for (method in setdiff(summarization.methods, "ave")) {
  
  # Get Genus-level abundance matrix in QIIME format
  outputdir <- paste("qiime-species", method, sep = "-")
  genus.hitchip <- get.qiime.matrix(outputdir, file.id = "abundance_with_taxa") 
  source("benchmarks.R")		  		
  taxa.correlations[common.taxa, method] <- taxcor[common.taxa]

}

par(mfrow = c(3,3))
for (method in setdiff(summarization.methods, "ave")) {
  hist(taxa.correlations[common.taxa, method], main = method, xlim = c(-1,1))
}
par(mar = c(10, 3, 1,1)); boxplot(taxa.correlations, las = 2)
plot(taxa.correlations[, "frpa"], taxa.correlations[, "sum"]); abline(0,1)
par(mar = c(2, 10, 1, 1)); barplot(sort(apply(taxa.correlations, 2, function (x) {mean(na.omit(x))})), las = 1, horiz = T)
par(mar = c(2, 10, 1, 1)); barplot(sort(apply(taxa.correlations, 2, function (x) {median(na.omit(x))})), las = 1, horiz = T)

