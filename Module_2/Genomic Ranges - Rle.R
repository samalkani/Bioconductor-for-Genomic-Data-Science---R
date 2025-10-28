# Genomic Ranges - Rle

# 1. Load library
library(GenomicRanges)

# 2. Construct an RLE (Run Length Encoded vectors)
rl = Rle(c(1,1,1,1,1,1,2,2,2,2,2,4,4,2))
rl

# 3. Access lengths using runLength()
runLength(rl)

# 4. Access values using runValue()
runValue(rl)

# 5. Converting RLE vector to a normal vector
as.numeric(rl)

# 6. Simulation
# A. IRanges
ir = IRanges(start = c(2,8), width = 4)

# B. Aggregate
aggregate(rl, ir, FUN = mean)

# C. IRanges
ir

# 7. Original vector
vec = as.numeric(rl)
vec
mean(vec[2:5])
mean(vec[8:11])

# 8. Construct new IRanges
ir = IRanges(start = 1:5, width = 3)
ir

# 9. Coverage - how often does an integer appear in an IRange?
coverage(ir)

# 10. Original RLE
rl

# 11. Slice RLE at a value of 2
slice(rl, 2)

# 12. Slice RLE at a value of 3
slice(rl, 3)

# 13. Instantiate our own views
vi = Views(rl, IRanges(2,8))
vi

# 14. Instantiate a set of views, (original vector: 1 1 1 1 1 1 2 2 2 2 2 4 4 2)
vi = Views(rl, IRanges(c(2,8), width = 2))
vi

# 15. Take the mean of the view
mean(vi)

# 16. Construct RLE with Genomic Ranges
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:10, width = 3))
rl <- coverage(gr)
rl

# 17. Instantiate a view on a genomic range RLE list
vi = Views(rl, GRanges("chr1", ranges = IRanges(3,7)))
vi

# 18. Instantiate a view on a genomic range RLE list (Ranges list)
# vi = Views(rl, as(GRanges("chr1", ranges = IRanges(3,7)), "RangesList"))
# vi
# Error in as(GRanges("chr1", ranges = IRanges(3, 7)), "RangesList") : 
# no method or default for coercing “GRanges” to “RangesList”

# 19. Look at a GRange
gr





















