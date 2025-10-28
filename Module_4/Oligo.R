# Oligo

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing "Oligo"
BiocManager::install(c("oligo"), force = TRUE)

# 3. Troubleshooting installing packages
BiocManager::valid()

BiocManager::install(c(
  "BiocParallel", "oligo", "promises", "reticulate", "Rsamtools"
), update = TRUE, ask = FALSE, force = TRUE)

# 4. Loading Packages
library(oligo)
library(GEOquery)

# 5. Get supplementary files from GEO
getGEOSuppFiles("GSE38792")

# 6. Displaying the TAR file
list.files("GSE38792")

# 7. Un-tar the files
untar("GSE38792/GSE38792_RAW.tar", exdir = "GSE38792/CEL")

# 8. list the untarred file
list.files("GSE38792/CEL")

# 9. List of file names here with the full path to it
celfiles <- list.files("GSE38792/CEL", full = TRUE)
celfiles

# 10. Convenience functions that read this data into a container
rawData <- read.celfiles(celfiles)

# 11. Display raw data
rawData

# 12. Displaying the gene class
Class("GeneFeatureSet")

# 13. Displaying the gene class
getClass("GeneFeatureSet")

# 14. Look at raw data expression access
exprs(rawData)[1:4,1:3]

# 15. Maximum value inside expression data
max(exprs(rawData))

# 16. Getting the filename from the raw data
filename <- sampleNames(rawData)
filename

# 17. Storing the raw data in pData
pData(rawData)$filename <- filename

# 18. Clean up sample names of the data set
sampleNames <- sub(".*_", "", filename)
sampleNames
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames
sampleNames(rawData) <- sampleNames
sampleNames

# 19. Use the sample names to create a group variable which will identify the 
# different groups
pData(rawData)$group <- ifelse(grepl("^OSA", sampleNames(rawData)),
                                                 "OSA", "control")
pData(rawData)

# 20. Create a Box plot of raw data (probeset intensities across samples)
boxplot(rawData, target = "core")

# 21. Use RMA to normalize micro array 
normData <- rma(rawData)

# 22. Displaying normData
normData

# 23. Displaying feature names of normData
featureNames(normData)[1:10]

# 24. Create a Box plot of normalized data (probeset intensities across samples)
boxplot(normData, target = "core")

# 25. Confirm micro array measurement are log transformed
exprs(normData)[1:3,1:3]









