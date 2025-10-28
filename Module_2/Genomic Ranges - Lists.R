# Genomic Ranges - Lists

# 1. Load library
library(GenomicRanges)

# 2. Constructing a GRanges list
gr1 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:4, width = 3))
gr2 <- GRanges(seqnames = "chr2", ranges = IRanges(start = 1:4, width = 3))
gL <- GRangesList(gr1 = gr1, gr2 = gr2)
gL

# 3. Indexing the 1st element in the GRanges list using [[]]
gL[[1]]

# 4. Indexing the 1st element in the GRanges list using the $ operator
gL$gr1

# 5. Using the start() access function
start(gL)

# 6. Using seqnames on GRanges list
seqnames(gL)

# 7. How long are each of the elements in GRanges (function out of date)

# First do this
elementLengths=elementNROWS

# Then run the function
elementLengths(gL)

# 8. Using sapply() function to find the length of each element in GRanges (slower)
sapply(gL, length)

# 9. Indoor Applies
shift(gL, 10)

# 10. FindOverlaps() function between GRanges List (gL) and Grange (gr2)
findOverlaps(gL, gr2)

# 11. FindOverlaps() function between GRanges List (gL) and GRanges List (gL)
findOverlaps(gL, gL)








