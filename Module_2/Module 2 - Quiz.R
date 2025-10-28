# Module 2 - Quiz

# 1.	What is the GC content of “chr22” in the 
# “hg19” build of the human genome?
  
# Tip: The reference genome includes “N” bases; 
# you will need to exclude those.

# 1. A. Load Libraries
library(AnnotationHub)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(GenomicFeatures)

# 1. B. Load the hg19 library
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

# 1. C. Load the Library
library(BSgenome.Hsapiens.UCSC.hg19)

# 1. D. Count total bases on chr22
alphabet_frequency <- alphabetFrequency(Hsapiens$chr22)
alphabet_frequency
total_bases <- sum(alphabet_frequency[c('A','G','C','T')])
total_bases

# 1. E. Count GC bases on chr22 
letter_freq <- letterFrequency(Hsapiens$chr22, "GC")
letter_freq

# 1. F. Calculate GC ratio
GC_content <- letter_freq/total_bases
GC_content

# 2.	Background: In the previous assessment we studied H3K27me3 “narrowPeak” regions from the H1 cell line 
# (recall that the Roadmap ID for this cell line is “E003”). We want to examine whether the GC content of 
# the regions influence the signal; in other words whether the reported results appear biased by GC content.

# What is mean GC content of H3K27me3 “narrowPeak” regions from Epigenomics Roadmap from the H1 stem cell 
# line on chr 22.

# Clarification: Compute the GC content for each peak region as a percentage and then average those percentages 
# to compute a number between 0 and 1.

# 2. A. Retrieve record
ah <- AnnotationHub()
ah
H3K27me3_qh <- query(ah, c("H3K27me3", "E003", "narrowPeak"))
H3K27me3_qh
H3K27me3_record <- H3K27me3_qh[["AH29892"]]
H3K27me3_record

# 2. B. Extract chr 22 and sequences
H3K27me3_chr22 <- subset(H3K27me3_record, seqnames == "chr22")
H3K27me3_chr22
H3K27me3_chr22_views <- Views(Hsapiens, H3K27me3_chr22)
H3K27me3_chr22_views

# 2. C. Calculate mean GC content
GC_contents <- letterFrequency(H3K27me3_chr22_views, "GC", as.prob = TRUE)
head(GC_contents)
tail(GC_contents)
mean_GC <- mean(GC_contents)
mean_GC

# 3.	Background - The “narrowPeak” regions includes information on a value 
# they call “signalValue”.

# What is the correlation between GC content and “signalValue” of these 
# regions (on chr22)?
  
signal_value <- mcols(H3K27me3_chr22_views)$signalValue
head(signal_value)
tail(signal_value)
cor(signal_value, GC_contents)

# 4.	Background: - The “narrowPeak” regions are presumably reflective of a ChIP
# signal in these regions. To confirm this, we want to obtain the “fc.signal” 
# data from AnnotationHub package on the same cell line and histone modification. 
# This data represents a vector of fold-change enrichment of ChIP signal over input.

# What is the correlation between the “signalValue” of the “narrowPeak” regions 
# and the average “fc.signal” across the same regions?
  
# Clarification: First compute the average “fc.signal” for across each region, 
# for example using “Views”; this yields a single number of each region. Next 
# correlate these numbers with the “signalValue” of the “narrowPeaks”.

# 4. A. Retrieve record
H3K27me3_fc <- query(ah, c("H3K27me3", "E003", "fc.signal"))
H3K27me3_fc
H3K27me3_fc_record <- H3K27me3_fc[["AH32033"]]
H3K27me3_fc_record

# 4. B. Get subset data on chr22
gr22 <- GRanges(seqnames = "chr22", ranges = IRanges(start = start(Hsapiens$chr22), end = end(Hsapiens$chr22)))
gr22
H3K27me3_fc_gr <- import(H3K27me3_fc_record, which = gr22, as = "Rle")
H3K27me3_fc_gr
H3K27me3_fc_gr22 <- H3K27me3_fc_gr$chr22
H3K27me3_fc_gr22

# 4. C. View fc.signal data
fc.signal <- Views(H3K27me3_fc_gr22, start = start(H3K27me3_chr22), end = end(H3K27me3_chr22))
fc.signal

# 4. D. Calculate the correlation between the average of fc.signal and signalValue
fc.signal_mean <- mean(fc.signal)
head(fc.signal_mean)
tail(fc.signal_mean)
cor(fc.signal_mean, signal_value)

# 5.	Referring to the objects made and defined in the previous question.

# How many bases on chr22 have an fc.signal greater than or equal to 1?

sum(H3K27me3_fc_gr22 >= 1)

# 6.	Background: The H1 stem cell line is an embryonic stem cell line, a 
# so-called pluripotent cell. Many epigenetic marks change upon differentiation. 
# We will examine this. We choose the cell type with Roadmap ID “E055” which is 
# foreskin fibroblast primary cells.

# We will use the “fc.signal” for this cell type for the H3K27me3 mark, on chr22.
# We now have a signal track for E003 and a signal track for E055. We want to 
# identify regions of the genome which gain H3K27me3 upon differentiation. 
# These are regions which have a higher signal in E055 than in E003. To do this 
# properly, we would need to standardize (normalize) the signal across the two 
# samples; we will ignore this for now.

# Identify the regions of the genome where the signal in E003 is 0.5 or lower 
# and the signal in E055 is 2 or higher.

# Tip: If you end up with having to intersect two different Views, note that you
# will need to convert the Views to IRanges or GRanges first with \verb|ir <- as
# (vi, "IRanges")|ir <- as(vi, "IRanges").

# 6. A. Find the specific record
H3K27me3_E055 <- query(ah, c("H3K27me3", "E055"))
H3K27me3_E055
H3K27me3_E055_record <- H3K27me3_E055[["AH32470"]]
H3K27me3_E055_record

# 6. B. Get subset data on chr22
gr_chr22 <- GRanges(seqnames = "chr22", ranges = IRanges(start = start(Hsapiens$chr22), end = end(Hsapiens$chr22)))
gr_chr22
H3K27me3_fc_gr_E055 <- import(H3K27me3_E055_record, which = gr_chr22, as = "Rle")
H3K27me3_fc_gr_E055
H3K27me3_fc_gr22_E055 <- H3K27me3_fc_gr_E055$chr22
H3K27me3_fc_gr22_E055

# 6. C. Identify region
region_E003 <- as(slice(H3K27me3_fc_gr22, upper = 0.5), "IRanges")
region_E003
region_E055 <- as(slice(H3K27me3_fc_gr22_E055, lower = 2), "IRanges")
region_E055
inter_region <- intersect(region_E003, region_E055)
inter_region

sum(width(inter_region))

# 7.	Background: - CpG Islands are dense clusters of CpGs. The classic definition
# of a CpG Island compares the observed to the expected frequencies of CpG dinucleotides
# as well as the GC content.

# Specifically, the observed CpG frequency is just the number of “CG” dinucleotides
# in a region. The expected CpG frequency is defined as the frequency of C 
# multiplied by the frequency of G divided by the length of the region.

# What is the average observed-to-expected ratio of CpG dinucleotides for CpG 
# Islands on chromosome 22?
  
# 7. A. Retrieve the cpg
ah_human <- subset(ah, species == "Homo sapiens")
ah_human
ah_human_cpg <- query(ah_human, "CpG Islands")
ah_human_cpg
ah_human_cpg_record <- ah_human_cpg[["AH5086"]]
ah_human_cpg_record

# 7. B. Get subset data on chr22
ah_human_cpg_chr22 <- subset(ah_human_cpg_record, seqnames == "chr22")
ah_human_cpg_chr22
ah_human_cpg_chr22_views <- Views(Hsapiens, ah_human_cpg_chr22)
ah_human_cpg_chr22_views

# 7. C. Calculate observed GC bases
observed_GC <- dinucleotideFrequency(ah_human_cpg_chr22_views)[,7]/width(ah_human_cpg_chr22_views)
head(observed_GC)
tail(observed_GC)

# 7. D. Calculate expected GC bases
freq_C <- letterFrequency(ah_human_cpg_chr22_views, "C")
head(freq_C)
tail(freq_C)
freq_G <- letterFrequency(ah_human_cpg_chr22_views, "G")
head(freq_G)
tail(freq_G)
expected_GC <- (freq_C/width(ah_human_cpg_chr22_views))*(freq_G/width(ah_human_cpg_chr22_views))
head(expected_GC)
tail(expected_GC)

# 7. E. Calculate the average observed-to-expected ratio of CpG dinucleotides
mean(observed_GC/expected_GC)

# 8.	Background: A TATA box is a DNA element of the form “TATAAA”. Around 25% 
# of genes should have a TATA box in their promoter. We will examine this statement.

# How many TATA boxes are there on chr 22 of build hg19 of the human genome?
  
# Clarification: You need to remember to search both forward and reverse strands.

TATA_boxes <- countPattern("TATAAA", Hsapiens$chr22) + countPattern("TATAAA", reverseComplement(Hsapiens$chr22))
TATA_boxes

# 9.	How many promoters of transcripts on chromosome 22 containing a coding sequence,
# contains a TATA box on the same strand as the transcript?
  
# Clarification: Use the TxDb.Hsapiens.UCSC.hg19.knownGene package to define transcripts
# and coding sequence. Here, we defined a promoter to be 900bp upstream and 100bp downstream
# of the transcription start site.

# 9. A. Load the TxDb.Hsapiens.UCSC.hg19 library
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

# 9. B. Load the Library
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# 9. C. Identify the record and create the GRanges
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
gr <- GRanges(seqnames = "chr22", ranges = IRanges(start = start(Hsapiens$chr22), end = end(Hsapiens$chr22)))
gr

# 9. D. find promoters of transcripts on chr 22
gr_trans_chr22 <- subsetByOverlaps(transcripts(txdb), gr, ignore.strand = TRUE)
gr_trans_chr22
proms <- promoters(gr_trans_chr22, upstream = 900, downstream = 100)
proms

# 9.  E. find coding sequences on chr 22
gr_cds_chr22 <- subsetByOverlaps(cds(txdb), gr, ignore.strand = TRUE)
gr_cds_chr22

# 9. F. find overlaps between promoters of transcripts and coding sequences
proms_cds <- findOverlaps(proms, gr_cds_chr22)
proms_cds

# 9. G. calculate TATA box on overlaps
count = 0
for (i in unique(queryHits(proms_cds))){
  proms_cds_view <- Views(Hsapiens, proms[i])
  count = count + vcountPattern("TATAAA", DNAStringSet(proms_cds_view))
}
count

# 10.	Background: It is possible for two promoters from different transcripts to 
# overlap, in which case the regulatory features inside the overlap might affect
# both transcripts. This happens frequently in bacteria.

# How many bases on chr22 are part of more than one promoter of a coding sequence?
  
# Clarification: Use the TxDb.Hsapiens.UCSC.hg19.knownGene package to define 
# transcripts and coding sequence. Here, we define a promoter to be 900bp 
# upstream and 100bp downstream of the transcription start site. In this case, 
# ignore strand in the analysis.

# 10. A. Calculate transcript lengths
trans_len_chr22 <- transcriptLengths(txdb, with.cds_len = TRUE)
head(trans_len_chr22)
tail(trans_len_chr22)
trans_len_chr22 <- trans_len_chr22[trans_len_chr22$cds_len > 0,]
head(trans_len_chr22)
tail(trans_len_chr22)

# 10. B. Find promoters from different transcripts to overlap
trans_eval <- proms[mcols(proms)$tx_id %in% trans_len_chr22$tx_id]
trans_eval
result = sum(coverage(trans_eval) > 1)
result
result["chr22"]

  