# BSGenome

# 1. Installing Bioconductor (BiocManager) packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

# 2. Installing Bioconductor packages
BiocManager::install(c("BSgenome", "DCchoice", "IRanges", "RcppArmadillo",
                       "rtracklayer", "BSgenome.Scerevisiae.UCSC.sacCer2" ), 
                     force = TRUE)

# 3. Troubleshooting installing packages
BiocManager::valid()

# 4. Loading Libraries
library(rtracklayer)
library(RcppArmadillo)
library(DCchoice)
library(BSgenome)
library(BSgenome.Scerevisiae.UCSC.sacCer2)

# 5. Available genomes from the Bioconductor website
available.genomes()

# 6. Loading a specific yeast genome
library(BSgenome.Scerevisiae.UCSC.sacCer2)

# 7. Call the short name of the genome object
Scerevisiae

# 8. Get the names with seqnames and lengths with seqlengths
seqnames(Scerevisiae)
seqlengths(Scerevisiae)

# 9. Load a specific chromosome into memory
Scerevisiae$chrI

# 10. Now you can call standard functions on this chromosome e.g. GC content
letterFrequency(Scerevisiae$chrI, "GC")

# 11. GC content of the chromosome 1 as a percentage
letterFrequency(Scerevisiae$chrI, "GC", as.prob = TRUE)

# 12. Calculating GC content for the entire genome
param = new("BSParams", X = Scerevisiae, FUN = letterFrequency)
unlist(bsapply(param, "GC"))

# 13. Calculate the GC content over the entire genome
sum(unlist(bsapply(param, "GC"))) / sum(seqlengths(Scerevisiae))

# 14. Summary of GC content across all chromosomes
unlist(bsapply(param, "GC", as.prob = TRUE))
