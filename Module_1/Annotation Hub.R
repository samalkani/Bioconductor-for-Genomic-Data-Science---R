# Annotation Hub

# 1. Installing the AnnotationHub package
BiocManager::install("AnnotationHub")

# 2. Load AnnotationHub
library(AnnotationHub)

# 3. Creating a local annotation hub
ah = AnnotationHub()
ah

# 4. Look at first element
ah[1]

# 5. Installing the rtracklayer package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rtracklayer")

# 6. Load rtracklayer package
library(rtracklayer)

# 7. Retrieving objects from annotation hub
ah[[1]]

# 8. Use dataprovider to display all the data providers
unique(ah$dataprovider)

# 9. Information on Species in the ah database
unique(ah$species)

# 10. Information Specific to Homo sapiens
ah = subset(ah, species == "Homo sapiens")
ah

# 11. Information specific to Histomodification
query(ah, "H3K4me3")

# 12. Information on a specific cell line and Histomodification
query(ah, c("H3K4me3", "Gm12878"))

# 13. Install BiocHubsShiny package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocHubsShiny")

# 14. Loading BiocHubsShiny package
library(BiocHubsShiny)

# 15. Using the spreadsheet-like interface
BiocHubsShiny()
write.table(ah_meta, file = 'ah_meta.txt')
print(ah_meta)
