# Read data
hit <- read.csv("~/projects/flint11/20130429/FLN\ L2\ for\ hm.txt", sep = "\t")
colnames(hit) <- gsub("X", "", colnames(hit))
colnames(hit) <- gsub("\\.", "", colnames(hit))
colnames(hit) <- gsub("ID", "L2", colnames(hit))

# Subject 12 is missing from M diet, remove from others too
hit <- hit[, setdiff(colnames(hit), c("12N", "12R", "12W"))]
subjects <- sapply(colnames(hit[, -c(1,2)]), function (s) {substr(s, 1, nchar(s)-1)})
diets <- sapply(colnames(hit[, -c(1,2)]), function (s) {substr(s, nchar(s), nchar(s))})
# Matched matrices for the diets
dietmatrices <- list()
for (diet in unique(diets)) {
  ind <- which(diets == diet); 
  dm <- hit[, -c(1,2)][, ind]; 
  colnames(dm) <- subjects[ind]
  dietmatrices[[diet]] <- dm
}
# Contrast other diets to M subject-wise and remove M diet
dietmatrices <- lapply(dietmatrices, function (mat) {mat - dietmatrices[["M"]]})
dietdata <- do.call("cbind", dietmatrices)
dietdata <- cbind(hit[, 1:2], dietdata[, -grep("M", colnames(dietdata))])
hit <- dietdata

print("Read qvalues")
qv <- read.csv("~/projects/flint11/20130429/qvals\ for\ hm.txt", sep = "\t")

print("Filter out non-significants")
keep <- which(apply(qv, 1, min) <= 0.054)
qv <- qv[keep, ]
hit <- hit[keep, ]

# Add sign to qvalues
qvalue.signs <- sign(sapply(split(3:ncol(hit), unname(sapply(names(hit[, -c(1,2)]), function (s) {substr(s, 1, 1)}))), function (inds) {rowMeans(hit[, inds])}))



library(reshape2)
data.table <- melt(hit)
data.table$diet <- factor(sapply(as.character(data.table$variable), function (x) {substr(x, nchar(x), nchar(x))}))
data.table$subject <- factor(sapply(as.character(data.table$variable), function (x) {substr(x, 1, nchar(x)-1)}))
# Order the rows
data.table$L2 <- factor(data.table$L2, levels = rev(as.character(hit$L2)))


# Load visualization libraries
library(ggplot2)

# Set theme
theme_set(theme_bw(15))

# Pick only the correlations with q<0.054 Note: this will leave other cells
# empty

print("Arrange the figure")
p <- ggplot(data.table, aes(x = variable, y = L2, fill = value))
p <- p + geom_tile()
#m <- max(abs(c(min(data.table$value), max(data.table$value))))
lims <- c(-1, 1)
p <- p + scale_fill_gradientn("value", breaks = seq(from = min(lims), to = max(lims), by = 0.5), colours = c("darkblue", "blue", "white", "red", "darkred"), limits = lims, name = "Fold Change (Log10)")
p <- p + theme(axis.text.x = theme_text(angle = 90, size = 7), axis.text.y = theme_text(size = 7)) + xlab("") + ylab("")

pdf("~/FigA.pdf", width = 8, height = 8)
print(p)
dev.off()

print("Mark the most significant cells with stars")
qv.bu <- qv
qvs <- qv[, 3:5]
qvs[(qv[, 3:5] <= 0.054) & (sign(qvalue.signs) == -1)] <- -1
qvs[(qv[, 3:5] <= 0.054) & (sign(qvalue.signs) == 1)] <- 1
qvs[(abs(qv[, 3:5]) > 0.054)] <- 0

qv[, 3:5] <- qvs
names(qv) <- gsub("N.M.qval", "N", names(qv))
names(qv) <- gsub("R.M.qval", "R", names(qv))
names(qv) <- gsub("W.M.qval", "W", names(qv))
qtable <- melt(qv)
qtable$Significance <- factor(qtable$value, levels = c(-1, 0, 1))
qtable$L1 <- factor(qtable$L1, levels = rev(unique(as.character(qtable$L1))))
qtable$L2 <- factor(qtable$L2, levels = rev(unique(as.character(qtable$L2))))

print("Significances")
p2 <- ggplot(qtable, aes(x = variable, y = L2, fill = Significance))
p2 <- p2 + geom_tile()
p2 <- p2 + theme(axis.text.x = element_text(angle = 0, vjust = .5), axis.text.y = element_text(size = 7))
p2 <- p2 + xlab("") + ylab("")
p2 <- p2 + scale_fill_manual(values = c("-1" = "blue", "0" = "white", "1" = "red"), guide = "none")


pdf("~/FigB.pdf", width = 2.8, height = 8)
print(p2)
dev.off()
