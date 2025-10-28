# Genomic Features

# 1. Installing Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

# 2. Installing Genomic features and TxDb.Hsapiens.UCSC.hg19.knownGene
BiocManager::install(c("GenomicFeatures", "TxDb.Hsapiens.UCSC.hg19.knownGene"), force = TRUE)

# 3. Troubleshooting installing packages
BiocManager::valid()

# 4. Load GenomicFeatures and TxDb.Hsapiens.UCSC.hg19.knownGene library
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# 5. Rename the DB object for convenience
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# 6. What is the DB object (TxDB)?
txdb

# 7. Select GRanges to examine on the DB object (TxDB)
gr <- GRanges(seqnames = "chr1", strand = "+", ranges = IRanges(start = 11874, end = 14409))

# 8. Call genes() on our DB object (txdB)
genes(txdb)

# 9. Let's look at what genes overlap with this Genomic range (gr)
subsetByOverlaps(genes(txdb), gr)

# 10. Let's look at what genes overlap with this Genomic range (gr) ignore strand
subsetByOverlaps(genes(txdb), gr, ignore.strand = TRUE)

# 11. Let's look at what transcripts overlap with this Genomic range (gr)
subsetByOverlaps(transcripts(txdb), gr)

# 12. Let's look at what exons overlap with this Genomic range (gr)
subsetByOverlaps(exons(txdb), gr)

# 13. How figure out how the exons are combined together to form transcripts?
subsetByOverlaps(exonsBy(txdb, by = "tx"), gr)

# 14. Let's look at what coding sequences (cds) overlap with this Genomic range (gr)
subsetByOverlaps(cds(txdb), gr)

# 15. Find out which transcript is the coding sequence (cds)
subsetByOverlaps(cdsBy(txdb, by = "tx"), gr)

# 16. Find out which exon is the coding sequence (cds)
subsetByOverlaps(exonsBy(txdb, by = "tx"), gr)["2"]

# 17. Find transcript lengths for the spliced RNA, not pre-mRNA
subset(transcriptLengths(txdb, with.cds_len = TRUE), gene_id == "100287102")

# 18. Checking the length of the coding sequence
sum(width(subsetByOverlaps(cdsBy(txdb, by = "tx"), gr)[["2"]]))















