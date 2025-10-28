# Genomic Ranges - Seqinfo

# 1. Installing the GenomicRanges package

# BiocManager::install("GenomicRanges")

# 2. Load GenomicRanges
library(GenomicRanges)

# 3. Construct GRanges with two chromosomes
gr = GRanges(seqnames = c("chr1","chr2"), ranges = IRanges(start = 1:2, end = 4:5))
gr

# 4. Selecting a specific chromosome using seqlevels
# seqlevels(gr, force = TRUE) = "chr1"
# gr
# Error in `seqlevels<-`(`*tmp*`, force = TRUE, value = "chr1") : 
#  unused argument (force = TRUE)

# 5. Selecting a specific chromosome using seqlevels updated code
dropSeqlevels(gr, "chr1", pruning.mode = "coarse")
seqlevels(gr, pruning.mode="coarse") = "chr1"
gr

# 6. Recreate the object
gr = GRanges(seqnames = c("chr1","chr2"), ranges = IRanges(start = 1:2, end = 4:5))
gr

# 7. Drop Sequence levels drop chr2 same as point 5.
dropSeqlevels(gr, "chr1", pruning.mode = "coarse")
seqlevels(gr, pruning.mode="coarse") = "chr1"
gr

# 8. Recreate the object
gr = GRanges(seqnames = c("chr1","chr2"), ranges = IRanges(start = 1:2, end = 4:5))
gr

# 9. Keep Sequence levels
keepSeqlevels(gr, "chr1", pruning.mode = "coarse")
seqlevels(gr, pruning.mode="coarse") = "chr1"
gr

# 10. Recreate the object with a weird chromosome name
gr = GRanges(seqnames = c("chr1","chrU345"), ranges = IRanges(start = 1:2, end = 4:5))
gr
keepStandardChromosomes(gr, pruning.mode="coarse")

# 11. Recreate the object 
gr = GRanges(seqnames = c("chr1","chr2"), ranges = IRanges(start = 1:2, end = 4:5))
gr

# 12. Creating NCBI style for chromosome names (Step 1)
newStyle = mapSeqlevels(seqlevels(gr), "NCBI")
newStyle

# 13. Renaming the GRanges object (Step 2)
gr = renameSeqlevels(gr, newStyle)
gr





