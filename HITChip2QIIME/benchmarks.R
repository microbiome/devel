
# Pick taxa that are available in both HITChip and NGS data
common.taxa <- intersect(rownames(genus.hitchip), rownames(genus.ngs))
common.samples <- intersect(colnames(genus.hitchip), colnames(genus.ngs))
hit <- genus.hitchip[common.taxa, common.samples]; hit[hit < 1e-4] <- 0
ngs <- genus.ngs[common.taxa, common.samples]

taxcor <- c()
for (pt in common.taxa) {
  #taxcor[[pt]] <- cor(log10(unlist(hit[pt, ])), log10(unlist(ngs[pt,])), method = "pearson")
  taxcor[[pt]] <- cor(log10(unlist(hit[pt, ])), log10(unlist(ngs[pt,])), method = "spearman")
}

samplecor <- c()
for (pt in common.samples) {
  a <- log10(unlist(hit[, pt]))
  b <- log10(unlist(ngs[,pt]))
  keep <- which(rowSums(is.infinite(cbind(a, b))) == 0)
  samplecor[[pt]] <- cor(a[keep], b[keep], method = "pearson")
}

#plot(unlist(ngs), unlist(hit), main = paste("Correlation:", round(cor(log10(1 + unlist(ngs)), log10(1 + unlist(hit)), use = "pairwise.complete.obs"), 3))); abline(0,1)
nv <- log10(unlist(ngs)) 
hv <- log10(unlist(hit))

# Ignore zero abundances and plot in log scale
keep <- rowSums(is.infinite(cbind(nv, hv))) == 0
plot(nv[keep], hv[keep], main = paste(outputdir, "corr:", round(cor(nv[keep], hv[keep]), 3))); abline(0,1)

# Correlation with zeroes excluded and at log10
nonzero.corr.log10 <- cor(nv[keep], hv[keep])



# Fraction of samples where given taxon is seen only on one of the platforms
#par(mfrow = c(1,2))
#v <- sort(rowMeans(hit == 0 & ngs > 0)); v <- v[v > 0]
#par(mar = c(4, 10, 3, 1)); barplot(v, horiz = T, las = 1, main = "Only NGS", cex.names = 1, xlab = "Fraction of samples where seen", xlim = c(0,1))
#v <- sort(rowMeans(hit > 0 & ngs == 0)); v <- v[v > 0]
#par(mar = c(4, 10, 3, 1)); barplot(v, horiz = T, las = 1, main = "Only HITChip", cex.names = 1, xlab = "Fraction of samples", xlim = c(0,1))

# Fraction of samples where given taxon is seen only on one of the platforms
# in considerable amounts 
#par(mfrow = c(1,2))
#v <- sort(rowMeans(hit == 0 & ngs > 0.01)); v <- v[v > 0]
#par(mar = c(4, 10, 3, 1)); barplot(v, horiz = T, las = 1, main = "Only NGS", cex.names = 1, xlab = "Fraction of samples where seen", xlim = c(0,1))
#v <- sort(rowMeans(hit > 0.01 & ngs == 0)); v <- v[v > 0]
#par(mar = c(4, 10, 3, 1)); barplot(v, horiz = T, las = 1, main = "Only HITChip", cex.names = 1, xlab = "Fraction of samples", xlim = c(0,1))


# Compare individual taxa
# k <- k+1; plot(unlist(ngs[k, ]), unlist(hit[k, ])); abline(0,1)