# R Samtools

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing "Rsamtools"
BiocManager::install(c("Rsamtools"), force = TRUE)

# 3. Troubleshooting installing packages
BiocManager::valid()

BiocManager::install(c(
  "BiocParallel", "promises", "reticulate", "Rsamtools"
), update = TRUE, ask = FALSE, force = TRUE)

# 4. Load Library "RSamtools"
library(Rsamtools)

# 5. Pointer down to a BAM file
bamPath <- system.file("extdata", "ex1.bam", package="Rsamtools")
bamFile <- BamFile(bamPath)
bamFile

# 6. Which sequences of chromosomes the reads align to?
seqinfo(bamFile)

# 7. Use SCANBAM to read in data into BAM file
aln <- scanBam(bamFile)
length(aln)

# 8. Class
class(aln)

# 9. Names (this is the entire file)
names(aln)

# 10. Getting the first element of the list and then print out the names
aln <- aln[[1]]
names(aln)

# 11. Look at the first element of each of these components
lapply(aln, function(xx) xx[1])

# 12. Setting the yield size on the BAM file so you only read parts of the entire read
yieldSize(bamFile) <- 1
bamFile

# 13. Open the bamFile before calling SCANBAM
open(bamFile)

# 14. Calling SCANBAM to output only 1 read
scanBam(bamFile)[[1]]$seq

# 15. Calling SCANBAM to output only 1 read
scanBam(bamFile)[[1]]$seq

# 16. Calling SCANBAM to output only 1 read
scanBam(bamFile)[[1]]$seq

# 17. Cleanup
close(bamFile)
yieldSize(bamFile) <- NA

# 18. Setting up a GRanges object
gr <- GRanges(seqnames = "seq2",
              ranges = IRanges(start = c(100, 1000), end = c(1500,2000)))
gr

# 19. Setting up SCANBAMPARAM to set up which genomic region were going to read in
params <- ScanBamParam(which = gr, what = scanBamWhat())
scanBamWhat()

# 20. Displaying the genomic ranges inputted earlier in the SCANBAMPARAM
aln <- scanBam(bamFile, param = params)
names(aln)

# 21. Displaying the positions between 100 and 1500
aln[[1]]$pos
head(aln[[1]]$pos)
tail(aln[[1]]$pos)

# aln[[2]]$pos
# head(aln[[2]]$pos)
# tail(aln[[2]]$pos)

# 22. Displaying the first few positions
head(aln[[1]]$pos)

# 23. Setting up a bamview with a single file in it
bamView <- BamViews(bamPath)
bamView

# 24. Reading in data with SCANBAM
aln <- scanBam(bamView)

# 25. Gives the Outer
names(aln)

# 26. Gives a list (first level)
names(aln[[1]])

# 27. Gives a list (second level), content
names(aln[[1]][[1]])

# 28. Setting ranges for the BAM view we get the outer as before
bamRanges(bamView) <- gr
aln <- scanBam(bamView)
names(aln)

# 29. Select the file name and we get the two sequences
names(aln[[1]])

# 30. Getting a quick summary of what is in the file
quickBamFlagSummary(bamFile)







