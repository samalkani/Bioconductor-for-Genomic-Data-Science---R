# GenomicRanges - Basic GRanges Usage

# 1. Installing the GenomicRanges package

# BiocManager::install("GenomicRanges")

# 2. Load GenomicRanges
library(GenomicRanges)

# 3. Construct IRanges and a data frame
ir = IRanges(start = 1:3, width = 2)
df = DataFrame(ir = ir, score = rnorm(3))
df

# 4. Subset IRange of length = 1
df[1,1]

# 5. Use the $ operator to print data frame
df$ir

# 6. Classic Data Frame in R with IRanges
df2 = data.frame(ir = ir)
df2

# 7. Constructing GRanges
gr <- GRanges(seqnames = "chr1", strand = c("+","-","+"), ranges = IRanges(start = c(1,3,5), width = 3))
gr

# 8. GRanges have something called values
values(gr) = DataFrame(score = rnorm(3))
gr

# 9. Accessing the values column in GRanges (gr)
values(gr)

# 10. Accessing the values column with mcols
mcols(gr)

# 11. Access a different column using dollar operator
gr$score

# 12. Assigning a new column
gr$score2 = gr$score / 3
gr

# 13. New GRanges
gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr1"), strand = "*", ranges = IRanges(start = c(1,3,5), width = 3))
gr2

# 14. Old GRanges
gr

# 15. Find Overlaps between gr and gr2
findOverlaps(gr, gr2)

# 16. Ignoring Strand in findOverlaps() function
findOverlaps(gr, gr2, ignore.strand = TRUE)

# 17. Subset Overlaps
subsetByOverlaps(gr, gr2)

# 18. Subset Overlaps (reverse)
subsetByOverlaps(gr2, gr)

# 19. Construct classic R data frame
df = data.frame(chr = "chr1", start = 1:3, end = 4:6, score = rnorm(3))
df

# 20. Converting classic R data frame to a GRanges
makeGRangesFromDataFrame(df)

# 21. Converting classic R data frame to GRanges keeping extra columns
makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)













