# Genomic Ranges - GRanges

# 1. Installing the IRanges package

BiocManager::install("GenomicRanges")

# 2. load IRanges Package
library(GenomicRanges)

# 3. Example of GRanges
gr = GRanges(seqnames = c("chr1"), strand = c("+","-","+"), ranges = IRanges(start = c(1,3,5), width = 3))
gr

# 4. Flanking sequence on GRange is relative to the direction of transcription
flank(gr, 5)

# 5. Promoters
promoters(gr)

# 6. Sequence Information
seqinfo(gr)

# 7. Give the Chromosome a length and then request for Sequence information
seqlengths(gr) = c("chr1" = 10)
seqinfo(gr)

# 8. Chromosome name given by seqlevels() function
seqlevels(gr)

# 9. Gaps() function covers everything not covered by IRanges and GRanges
gaps(gr)

# 10. Adding more Chromosomes
seqnames(gr) = c("chr1","chr2","chr1")

# 11. Use Seqlevels() to create two levels then assign seq values to it
seqlevels(gr) = c("chr1","chr2")
seqnames(gr) = c("chr1", "chr2", "chr1")
gr

# 12. Re-organize the seqlevels so "chr1" comes before "chr2" then sort
  # present order of chromosomes
  seqlevels(gr)
  # reorganise chromosomes so chr2 comes before chr1
  seqlevels(gr) = c("chr2","chr1")
  # Then sort order of chromosomes
  sort(gr)

# 13. Assigning a genome name
genome(gr) = "hg19"
gr

# 14. Displaying genome information of Seqinfo()
seqinfo(gr)

# 15. Make a copy of the genome
gr2 = gr
genome(gr2) = "hg18"

# 16. Find Overlaps between genome "hg19" and "hg18"
findOverlaps(gr, gr2)







