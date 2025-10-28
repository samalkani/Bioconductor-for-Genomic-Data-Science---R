# Summarized Experiment

# 1. Install Summarized Experiment and airway packages

# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
BiocManager::install("airway")
BiocManager::install("GenomicFiles")

# 2. Load the Libraries
library(GenomicRanges)
library(airway)
library(Biobase)
library(GenomicFiles)
library(SummarizedExperiment)

# 3. Accessing data from the airway library
data(airway)
airway

# 4. Use colData() to access the same information from the summarized experiment as using
# pData() in the expression set
colData(airway)

# 5. Accessing a specific column
airway$cell

# 6. Look at Experiment Data using exptData()
# exptData(airway)
# exptData(), this function is now defunct


# 7. Alternative function metadata()
metadata(airway)

# 8. Obtaining names for the different samples
colnames(airway)

# 9. Obtaining feature names
head(rownames(airway))

# 10. Look at the airway data itself
airway

# 11. List of all the assays
assayNames(airway)

# 12. Obtaining the assay
assay(airway, "counts")[1:4,1:4]

# 13. Obtaining the length (number) of rows
length(rowRanges(airway))

# 14. Using rowRanges() function to obtain GRanges list
rowRanges(airway)

# 15. elementLengths() function is defunct use elementNROWS
head(elementNROWS(rowRanges(airway)))

# 16. How many exons do we have per gene?
sum(elementNROWS(rowRanges(airway)))

# 17. Start co-ordinates of each of the exons
start(airway)

# 18. Alternative use subset by overlaps function GRanges
gr <- GRanges("1", ranges = IRanges(start = 1, end = 10^7))
gr

# 19. Subset by overlaps
subsetByOverlaps(airway, gr)










