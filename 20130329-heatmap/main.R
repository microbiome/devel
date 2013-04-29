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


# Read qvalues
qv <- read.csv("~/projects/flint11/20130429/qvals\ for\ hm.txt", sep = "\t")

# Filter out non-significants
keep <- which(apply(qv, 1, min) < 0.05)
qv <- qv[keep, ]
hit <- hit[keep, ]



# Visualize each bacteria against its abundance in M diet same subject
#hit[, -c(1, 2)] <- hit[, -c(1,2)] - matrix(rowMeans(hit[, grep("M", colnames(hit))]))
#hit[, -c(1, 2)] <- hit[, -c(1,2)] - matrix(rowMeans(hit[, grep("M", colnames(hit))]))

data.table <- melt(hit)
data.table$diet <- factor(sapply(as.character(data.table$variable), function (x) {substr(x, nchar(x), nchar(x))}))
data.table$subject <- factor(sapply(as.character(data.table$variable), function (x) {substr(x, 1, nchar(x)-1)}))
# Order the rows
data.table$L2 <- factor(data.table$L2, levels = rev(as.character(hit$L2)))


# Load visualization libraries
library(ggplot2)

# Set theme
theme_set(theme_bw(15))

# Pick only the correlations with q<0.05 Note: this will leave other cells
# empty

# Arrange the figure
p <- ggplot(data.table, aes(x = variable, y = L2, fill = value))
p <- p + geom_tile()
#m <- max(abs(c(min(data.table$value), max(data.table$value))))
lims <- c(-1.25, 1.25)
p <- p + scale_fill_gradientn("value", breaks = seq(from = min(lims), to = max(lims), by = 0.5), colours = c("darkblue", "blue", "white", "red", "darkred"), limits = lims)
p <- p + theme(axis.text.x = theme_text(angle = 90, size = 7), axis.text.y = theme_text(size = 7)) + xlab("") + ylab("")

pdf("~/FigA.pdf")
print(p)
dev.off()


# Mark the most significant cells with stars
qv[, 3:5] <- qv[, 3:5] < 0.05
qtable <- melt(qv)
qtable$L2 <- factor(qtable$L2, levels = rev(as.character(hit$L2)))

#p <- p + geom_text(data = subset(qtable, qvalue < 0.05), aes(x = X1, y = X2, label = "+"), col = "white", size = 5)

# Arrange the figure
p2 <- ggplot(qtable, aes(x = variable, y = L2, fill = value))
p2 <- p2 + geom_tile()
#lims <- c(0,1)
#p <- p + scale_fill_gradientn("value", breaks = seq(from = min(lims), to = max(lims), by = 0.5), colours = c("darkblue", "blue", "white", "red", "darkred"), limits = lims)
#p <- p + scale_fill_gradientn("value", breaks = seq(from = min(lims), to = max(lims), by = 0.5), colours = c("darkblue", "blue", "white", "red", "darkred"), limits = lims)
p2 <- p2 + theme(axis.text.x = theme_text(angle = 90), axis.text.y = theme_text(size = 7)) + xlab("") + ylab("")

pdf("~/FigB.pdf")
print(p2)
dev.off()