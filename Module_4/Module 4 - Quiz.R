# Module 4 - Quiz

# Question 1

# 1.	The yeastRNASeq experiment data package contains FASTQ files from an RNA 
# seq experiment in yeast. When the package is installed, you can access one of 
# the FASTQ files by the path given by

# library(yeastRNASeq)
# fastqFilePath <- system.file("reads", "wt_1_f.fastq.gz", 
# package = "yeastRNASeq")

# What fraction of reads in this file has an A nucleotide in the 5th base of the
# read?

# 1. A. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 1. B. Installing "yeastRNASeq" and "leeBamViews"
BiocManager::install(c("yeastRNASeq", "leeBamViews", "zebrafishRNASeq"), force = TRUE)

# 1. C. Troubleshooting installing packages
BiocManager::valid()

BiocManager::install(c(
  "BiocParallel", "DESeq2", "oligo", "promises", "Rsamtools", "xml2"
), update = TRUE, ask = FALSE, force = TRUE)

# 1. D. Load "yeastRNASeq"
library(yeastRNASeq)

# 1. E. Pointer down to FASTQ file
fastqFilePath <- system.file("reads", "wt_1_f.fastq.gz", package = "yeastRNASeq")

# 1. F. Load Libraries
library(ShortRead)
library(Biostrings)

# 1. G. Read fastq file
reads <- readFastq(fastqFilePath)
reads
DNAStringSet <- sread(reads)
DNAStringSet

# 1. H. Fraction of reads has an A in the 5th base
cm <- consensusMatrix(DNAStringSet, as.prob=TRUE, baseOnly=TRUE)
cm
cm['A', 5]


# Question 2

# This is a continuation of Question 1.

# What is the average numeric quality value of the 5th base of these reads?

# 2. A. Average numeric quality value of base 5
mean(as(quality(reads), "matrix")[,5])

# Question 3

# The leeBamViews experiment data package contains aligned BAM files from an 
# RNA seq experiment in yeast (the same experiment as in Questions 1 and 2, but 
# that is not pertinent to the question). You can access one of the BAM files by 
# the path given by

# library(leeBamViews)
# bamFilePath <- system.file("bam", "isowt5_13e.bam", package="leeBamViews")

# These reads are short reads (36bp) and have been aligned to the genome using 
# a standard aligner, ie. potential junctions have been ignored (this makes some
# sense as yeast has very few junctions and the reads are very short).

# A read duplicated by position is a read where at least one more read shares the same position.

# We will focus on the interval from 800,000 to 801,000 on yeast chromosome 13.

# In this interval, how many reads are duplicated by position?

# 3. A. Load "leeBamViews"
library(leeBamViews)

# 3. B. Pointer down to BAM file
bamFilePath <- system.file("bam", "isowt5_13e.bam", package="leeBamViews")

# 3. C. Load “Rsamtools”
library(Rsamtools)

bamFile <- BamFile(bamFilePath)
bamFile

# 3. D. Focus on Scchr13, interval from 800,000 to 801,000
gr <- GRanges(seqnames = "Scchr13", ranges = IRanges(start = c(800000), end = c(801000)))
gr
params <- ScanBamParam(which = gr, what = scanBamWhat())
params
aln <- scanBam(bamFile, param = params)
aln

# Question 4

# This is a continuation of Question 3.

# The package contains 8 BAM files in total, representing 8 different samples 
# from 4 groups. A full list of file paths can be had as

# bpaths <- list.files(system.file("bam", package="leeBamViews"), 
#                      pattern = "bam$", full=TRUE)

# An objective of the original paper was the discovery of novel transcribed 
# regions in yeast. One such region is Scchr13:807762-808068.

# What is the average number of reads across the 8 samples falling in this 
# interval?
  
# 4. A. Pointer down to the BAM file 
bpaths <- list.files(system.file("bam", package="leeBamViews"), pattern = "bam$", full=TRUE)

# 4. B. Focus on the novel transcribed regions
bamView <- BamViews(bpaths)
bamView
gr_nt <- GRanges(seqnames="Scchr13", ranges=IRanges(start = c(807762), end = c(808068)))
gr_nt
bamRanges(bamView) <- gr_nt
aln_nt <- scanBam(bamView)
aln_nt

# 4. C. Get sequences for each sample
alns <- lapply(aln_nt, function(xx) xx[[1]]$seq)
alns

# 4. D. Calculate the average number of reads across 8 the samples
alns_len_sum = 0
for (i in 1:length(alns)){
  alns_len_sum = alns_len_sum + length(alns[i][[1]])
}
alns_len_sum / length(alns)

# Question 5

# 5.	In the lecture on the oligo package an ExpressionSet with 18 samples is 
# constructed, representing normalized data from an Affymetrix gene expression 
# microarray. The samples are divided into two groups given by the \verb|group|
# group variable.

# What is the average expression across samples in the control group for the 
# “8149273” probeset (this is a character identifier, not a row number).

# 5. A. Load Libraries
library(oligo)
library(GEOquery)

# 5. B. Get data
getGEOSuppFiles("GSE38792")
untar("GSE38792/GSE38792_RAW.tar", exdir = "GSE38792/CEL")

# 5. C. Read data
celfiles <- list.files("GSE38792/CEL", full = TRUE)
rawData <- read.celfiles(celfiles)

# 5. D. Parse pData
filename <- sampleNames(rawData)
filename
pData(rawData)$filename <- filename
filename
sampleNames <- sub(".*_", "", filename)
sampleNames
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames
sampleNames(rawData) <- sampleNames
sampleNames
pData(rawData)$group <- ifelse(grepl("^OSA", sampleNames(rawData)), "OSA", "Control")
rawData

# 5. E. Find "8149273" probeset
normData <- rma(rawData)
loc <- match("8149273", rownames(normData))
loc

# 5. F. Average expression in control group
mean(exprs(normData[loc,])[1:8])

# Question 6

# 6.	This is a continuation of Question 5.

# Use the limma package to fit a two group comparison between the control group 
# and the OSA group, and borrow strength across the genes using \verb|eBayes()|eBayes().
# Include all 18 samples in the model fit.

# What is the absolute value of the log foldchange (\verb|logFC|logFC) of the 
# gene with the lowest \verb|P.value|P.value.

# 6. A. Load Library
library(limma)

# 6. B. Use limma to fit between control group and OSA group
normData$group <- factor(normData$group)
normData$group
design <- model.matrix(~normData$group)
design
fit <- lmFit(normData, design)
fit
fit <- eBayes(fit)
fit

# 6. C. Absolute value of logFC which has lowest P.value
abs(topTable(fit)$logFC[1])

# Question 7

# 7.	This is a continuation of Question 6.

# How many genes are differentially expressed between the two groups at an 
# \verb|adj.P.value|adj.P.value cutoff of 0.05?
  
fit_toptable <- topTable(fit)
fit_toptable
de <- subset(fit_toptable, adj.P.Val < 0.05)
de

# Question 8

# 8.	An example 450k dataset is contained in the minfiData package. This dataset
# contains 6 samples; 3 cancer and 3 normals. Cancer has been shown to be globally 
# hypo-methylated (less methylated) compared to normal tissue of the same kind.

# Take the RGsetEx dataset in this package and preprocess it with the 
# preprocessFunnorm function. For each sample, compute the average Beta value 
# (percent methylation) across so-called OpenSea loci.

# What is the mean difference in beta values between the 3 normal samples and 
# the 3 cancer samples, across OpenSea CpGs?
  
# 8. A. Load Libraries
library(minfi)
library(minfiData)

# 8. B. Get OpenSea loci in RGsetEx with preprocess
rgSet <- preprocessFunnorm(RGsetEx)
rgSet
rg_opensea <- rgSet[c(getIslandStatus(rgSet) == "OpenSea")]
rg_opensea

# 8. C. Get Beta value in both group
rg_beta <- getBeta(rg_opensea)
head(rg_beta)
normal <- mean(rg_beta[, c(1,2,5)])
normal
cancer <- mean(rg_beta[, c(3,4,6)])
cancer

# 8. D. Mean difference between normal and cancer group
normal - cancer

# Question 9
  
# 9.	This is a continuation of Question 8.

# The Caco2 cell line is a colon cancer cell line profiled by ENCODE. Obtain the
# narrowPeak DNase hyper sensitive sites computed by the analysis working group 
# (AWG).

# How many of these DNase hypersensitive sites contain one or more CpGs on the 
# 450k array?

# 9. A. Load Library
library(AnnotationHub)

# 9. B. Get Caco2 data
ah <- AnnotationHub()
ah
ah <- subset(ah, species=="Homo sapiens")
ah
ah_Caco2 <- query(ah, c("Caco2", "AWG"))
ah_Caco2
ah_Caco2 <- ah_Caco2[["AH22442"]]
ah_Caco2

CpG_450K <- granges(rgSet)
CpG_450K

unique(findOverlaps(CpG_450K, ah_Caco2, type="within"))

# Question 10  

# 10.	The zebrafishRNASeq package contains summarized data from an RNA-seq 
# experiment in zebrafish in the form of a data.frame called \verb|zfGenes|zfGenes. 
# The experiment compared 3 control samples to 3 treatment samples.

# Each row is a transcript; the data.frame contains 92 rows with spikein 
# transcripts; these have a rowname starting with “ERCC”. Exclude these rows 
# from the analysis.

# Use DESeq2 to perform a differential expression analysis between control and 
# treatment. Do not discard (filter) genes and use the \verb|padj|padj results 
# output as the p-value.

# How many features are differentially expressed between control and treatment 
# (ie. \verb|padj <= 0.05|padj <= 0.05)?

# 10. A. Load Libraries
library(DESeq2)
library(zebrafishRNASeq)

# 10. B. Get and parse data
data("zfGenes")
data
zf <- zfGenes[grep("^ERCC", rownames(zfGenes), invert = T), ]
head(zf)
zf <- as.matrix(zf)
head(zf)
colData <- DataFrame(sampleID = colnames(zf), group = as.factor(c("control", "control", "control", "treatment", "treatment", "treatment")))
colData

# 10. C. Perform DESeq2
dds <- DESeqDataSetFromMatrix(zf, colData, design = ~ group)
dds
dds <- DESeq(dds)
dds

# 10. D. Find differentially expressed features
res <- results(dds)
res
sigRes <- subset(res, padj <= 0.05)
sigRes
dim(sigRes)[1]






  



