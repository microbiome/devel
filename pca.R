trans <- "clr"
distmeth <- "euclidean"

# Transformation on Genus abundances
o <- abundances(transform(phfinrisk_genus, trans))

# Sample dissimilarities
d <- vegdist(t(o), distmeth)

# PCA
pca <- princomp(t(o))
