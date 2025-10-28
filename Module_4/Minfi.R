# Minfi

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing "minfi"
BiocManager::install(c("minfi"), force = TRUE)

# 3. Troubleshooting installing packages
BiocManager::valid()

BiocManager::install(c(
  "BiocParallel", "oligo", "promises", "Rsamtools"
), update = TRUE, ask = FALSE, force = TRUE)

# 4. Load Libraries
library(minfi)
library(GEOquery)

# 5. Download the supplementary files
getGEOSuppFiles("GSE68777")

# 6. Untar the files
untar("GSE68777/GSE68777_RAW.tar", exdir = "GSE68777/idat")

# 7. Display idat files
head(list.files("GSE68777/idat", pattern = "idat"))

# 8. The files are compressed in the "GSE68777/idat" directory
idatFiles <- list.files("GSE68777/idat", pattern = "idat.gz$", full = TRUE)
# sapply(idatFiles, gunzip, overwrite = TRUE)
library(R.utils)

for (x in idatFiles) {
  gunzip(x, remove=FALSE)
}

# 9. Read in uncompressed files (read.450k.exp() function defunct now)
# rgSet <- read.450k.exp("GSE68777/idat")
rgSet <- read.metharray.exp("GSE68777/idat")

# 10. Alternative substitute function (reading in uncompressed IDAT data file)
rgSet <- read.metharray.exp(base = "GSE68777/idat", recursive = TRUE)
rgSet

# 11. Accessing phenotype data
pData(rgSet)

# 12. Accessing the sample name from the rgSet (no phenotype information)
head(sampleNames(rgSet))

# 13. Downloading the original process data
geoMat <- getGEO("GSE68777")

# 14. Getting phenotype data from original process data
pD.all <- pData(geoMat[[1]])
pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1", "characteristics_ch1.2")]
head(pD)

# 15. Clean up the data rename the columns
names(pD)[c(3,4)] <- c("group", "sex")
pD$group <- sub("^diagnosis: ", "", pD$group)
pD$sex <- sub("^Sex: ", "", pD$sex)
head(pD)

# 16. Continue clean up on the rgSet
sampleNames(rgSet) <- sub(".*_5", "5", sampleNames(rgSet))
rownames(pD) <- pD$title
pD <- pD[sampleNames(rgSet),]

# 17. Display Sample Names
head(sampleNames(rgSet))

# 18. Display the first 10 records on the pD object
head(pD)

# 19. Ready to assign it
# pData(rgSet) <- pD
# rgSet
pD = pData(rgSet)
rgSet

# 20. Pre-processing for the 450K micro-array (normalization)
grSet <- preprocessQuantile(rgSet)

# 21. Summarised experiment set (genomic ratio set)
grSet

# 22. Getting genomic ranges on the CpG's (methylated parts of genome)
granges(grSet)

# 23. Getting information on the status of a CpG  
head(getIslandStatus(grSet))

# 24. Getting Beta values for the (normalized) rgSet
getBeta(grSet)[1:3,1:3]
























