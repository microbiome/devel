# COMPARE TO NGS DATA

# Get Genus-level abundance matrix in QIIME format
genus.hitchip <- get.qiime.matrix(outputdir, file.id = "abundance_with_taxa") 

# Read NGS data for the same samples
genus.ngs <- get.qiime.matrix("~/projects/hitchip-benchmarking13/data", file.id = "otu_table") 
# Rename columns
colnames(genus.ngs) <- gsub("\\.", "", gsub(".l3.2", "", gsub(".l3.1", "", colnames(genus.ngs))))

# Pick taxa that are available in both HITChip and NGS data
common.taxa <- intersect(rownames(genus.hitchip), rownames(genus.ngs))
common.samples <- intersect(colnames(genus.hitchip), colnames(genus.ngs))
hit <- genus.hitchip[common.taxa, common.samples]; hit[hit < 1e-4] <- 0
ngs <- genus.ngs[common.taxa, common.samples]

#plot(unlist(ngs), unlist(hit), main = paste("Correlation:", round(cor(log10(1 + unlist(ngs)), log10(1 + unlist(hit)), use = "pairwise.complete.obs"), 3))); abline(0,1)
nv <- log10(unlist(ngs)) 
hv <- log10(unlist(hit))

# Ignore zero abundances and plot in log scale
keep <- rowSums(is.infinite(cbind(nv, hv))) == 0
plot(nv[keep], hv[keep], main = paste("Correlation:", round(cor(nv[keep], hv[keep]), 3))); abline(0,1)


# Fraction of samples where given taxon is seen only on one of the platforms
par(mfrow = c(1,2))
v <- sort(rowMeans(hit == 0 & ngs > 0)); v <- v[v > 0]
par(mar = c(4, 10, 3, 1)); barplot(v, horiz = T, las = 1, main = "Only NGS", cex.names = 1, xlab = "Fraction of samples where seen", xlim = c(0,1))
v <- sort(rowMeans(hit > 0 & ngs == 0)); v <- v[v > 0]
par(mar = c(4, 10, 3, 1)); barplot(v, horiz = T, las = 1, main = "Only HITChip", cex.names = 1, xlab = "Fraction of samples", xlim = c(0,1))

# Fraction of samples where given taxon is seen only on one of the platforms
# in considerable amounts 
par(mfrow = c(1,2))
v <- sort(rowMeans(hit == 0 & ngs > 0.01)); v <- v[v > 0]
par(mar = c(4, 10, 3, 1)); barplot(v, horiz = T, las = 1, main = "Only NGS", cex.names = 1, xlab = "Fraction of samples where seen", xlim = c(0,1))
v <- sort(rowMeans(hit > 0.01 & ngs == 0)); v <- v[v > 0]
par(mar = c(4, 10, 3, 1)); barplot(v, horiz = T, las = 1, main = "Only HITChip", cex.names = 1, xlab = "Fraction of samples", xlim = c(0,1))


# Compare individual taxa
# k <- k+1; plot(unlist(ngs[k, ]), unlist(hit[k, ])); abline(0,1)

