# Plot a Motion Chart using googleVis package - a brief example
# Leo Lahti 2012

# Load the library
# Note: requires flash and internet connection
# use install.packages("googleVis") to install the library if needed

# Load the data
load("~/Downloads/forleo180912.RData")

# Load libraries
library(microbiome)
library(googleVis)

# Form a motion chart
# NOTE: the data set must be given as data.frame
# which can contain ONLY NUMERIC and CHARACTER fields 
# (NO FACTORS, NO LOGICAL variables!). 
# The FIRST FOUR FIELDS must be provided in the following order: 
# idvar, timevar, two numeric fields, 
# then any number of numeric and character fields
# See help(gvisMotionChart) for details

df <- list()
# Ensure first four fields as ID, time, numeric variable 1, numeric variable 2
df$sampleID <- paste("Sample-", 1:nrow(d.i_tmp), sep = "")
df$time <- rep(1, nrow(d.i_tmp))
df$variable1 <- d.i_tmp$d.i.anm
# Then add the other variables
df <- cbind(df, d.i_tmp)
# Convert all factors to characters
i <- sapply(df, is.factor); 
df[i] <- lapply(df[i], as.character)


# Plot a Motion Chart using googleVis -package
library(googleVis)

#mchart <- gvisMotionChart(df, idvar="sampleID", timevar="constant")
mchart <- gvisMotionChart(df, idvar="sampleID", timevar="time")

# Plot immediately (opens in browser)
plot(mchart)

# Save as html (needs javascript to open!) 
print(mchart, file="~/roihu/MotionChart.html")
