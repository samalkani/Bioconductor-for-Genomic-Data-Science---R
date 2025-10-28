# GEO query

# 1. Installing Bioconductor
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing "GEOquery"
BiocManager::install(c("GEOquery"), force = TRUE)

# 3. Troubleshooting installing packages
BiocManager::valid()

BiocManager::install(c(
  "BiocParallel", "promises", "Rsamtools"
), update = TRUE, ask = FALSE, force = TRUE)

# 4. Load Libraries
library(GenomicRanges)
library(GEOquery)
library(Biobase)
library(GenomicFiles)

# 5. Downloading data from GEO
elist <- getGEO("GSE11675")

# 6. Class
class(elist)

# 7. Length of list
length(elist)

# 8. Name of elist
names(elist)

# 9. Get the data
eData = elist[[1]]
eData

# 10. The names of the phenotype data
names(pData(eData))

# 11. Accessing the raw (supplementary) data
eList2 <- getGEOSuppFiles("GSE11675")

# 12. Look at the downloaded tar files
eList2




