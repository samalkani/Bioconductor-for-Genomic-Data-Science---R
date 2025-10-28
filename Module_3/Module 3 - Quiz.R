# Module 3 - Quiz

# 1. What is the mean expression across all features for sample 5 in the ALL 
# dataset (from the ALL package)?

# 1. A. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 1. B. Installing ALL Package"
# BiocManager::install("ALL", force = TRUE)

# 1. C. Load the Libraries
library(ALL)

# 1. D. Create the expression set
data(ALL)
head(exprs(ALL[,5]))
tail(exprs(ALL[,5]))

# 1. E. Select sample 5 (column 5 from the expression matrix) find the mean value
mean(exprs(ALL[,5]))

# 2. We will use the biomaRt package to annotate an Affymetrix microarray. We 
# want our results in the hg19 build of the human genome and we therefore need 
# to connect to Ensembl 75 which is the latest release on this genome version. 
# How to connect to older versions of Ensembl is described in the biomaRt package 
# vignette; it can be archived with the command 

# mart <- useMart(host='feb2014.archive.ensembl.org', biomart = 
#                   "ENSEMBL_MART_ENSEMBL")|mart <- useMart(host=’feb2014.archive.ensembl.org’, 
#                                                         biomart = "ENSEMBL_MART_ENSEMBL")

# Using this version of Ensembl, annotate each feature of the ALL dataset with 
# the Ensembl gene id. How many probesets (features) are annotated with more than 
# one Ensembl gene id?

# 2. A. Installing "biomaRt", "hgu95av2.db" Packages
BiocManager::install(c("biomaRt", "hgu95av2.db"), force = TRUE)

# 2. B. Load the libraries
library(biomaRt)
library(hgu95av2.db)

# 2. C. list Marts
mart <- useMart(host='feb2014.archive.ensembl.org', biomart = "ENSEMBL_MART_ENSEMBL")
mart
ensembl <- useDataset("hsapiens_gene_ensembl", mart)
ensembl

# 2. D. annotate each feature
feature_name <- featureNames(ALL)
head(feature_name)
annotation_ALL <- getBM(attributes=c("ensembl_gene_id","affy_hg_u95av2"), 
                        filters="affy_hg_u95av2", values=feature_name, mart=ensembl)
head(annotation_ALL)
sum(table(annotation_ALL[,2])>1)

# 3. How many probesets (Affymetrix IDs) are annotated with one or more genes on 
# the autosomes (chromosomes 1 to 22).

# 3. A. list Attributes
attributes <- listAttributes(ensembl)
head(attributes)
filters <- listFilters(ensembl)
head(filters)

# 3. B. annotate autosomes
chrom <- c(1:22)
chrom
annotation_ALL_chr <- getBM(attributes=c("ensembl_gene_id", "affy_hg_u95av2", "chromosome_name"),
                            filters=c("affy_hg_u95av2","chromosome_name"), 
                            values=list(feature_name, chrom), mart=ensembl)
head(annotation_ALL_chr,)
head(table(annotation_ALL_chr[,2]))
head(table(table(annotation_ALL_chr[,2])), 17)
sum(table(table(annotation_ALL_chr[,2])))

# 4. Use the MsetEx dataset from the minfiData package. Part of this question is
# to use the help system to figure out how to address the question.

# What is the mean value of the Methylation channel across the features for sample 
# “5723646052_R04C01”?
  
# 4. A. Installing "minfiData", "minfi"
BiocManager::install(c("minfiData", "minfi"), force = TRUE)

# 4. B. Troubleshooting installing packages
BiocManager::valid()

# 4. C. Load minfiData library    
library(minfiData)
library(minfi)

# 4. D. Select column 2 of expression set (sample “5723646052_R04C01”)
# Get mean value of the sample (Mehylation channel) across all features
MsetEx
MsetEx[,2]
head(getMeth(MsetEx)[,2])
mean(getMeth(MsetEx)[,2])    

# 5. Access the processed data from NCBI GEO Accession number GSE788. 
# What is the mean expression level of sample GSM9024?

# 5. A. Load library Gene Expression Omnibus (GEO)    
library(GEOquery)

# 5. B. Download data from GEO using accession number "GSE788"
eList <- getGEO("GSE788")
eList
eData <- eList[[1]]
eData

# 5. C. Select column 2, sample "GSM9024" and find mean
head(exprs(eData)[1:5,1:5])
mean(exprs(eData)[,2])

# 6. We are using the airway dataset from the airway package.

# What is the average of the average length across the samples in the experiment?
  
# 6. A. Loading libraries    
library(airway)
library(GenomicRanges)

# 6. B. Create expression set
data(airway)
airway

# 6 C. Select column avgLength take the mean
airway$avgLength
mean(airway$avgLength)

# 7. We are using the airway dataset from the airway package. The features in this 
# dataset are Ensembl genes.

# What is the number of Ensembl genes which have a count of 1 read or more in 
# sample SRR1039512?

# 7. A. Sum the no. of genes that have a read count of >1 for sample SRR1039512 (column 3)
head(assay(airway))
sum(assay(airway)[,3]>=1)

# 8. The airway dataset contains more than 64k features. How many of these features 
# overlaps with transcripts on the autosomes (chromosomes 1-22) as represented by the 
# TxDb.Hsapiens.UCSC.hg19.knownGene package?
  
# Clarification: A feature has to overlap the actual transcript, not the intron 
# of a transcript.  So you will need to make sure that the transcript representation 
# does not contain introns.

# 8. A. Load libraries  
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# 8. B. exon data of txdb
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
txdb_exons <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb_exons

# 8. C. transcripts on the autosome
autosome <- paste0("chr", c(1:22))
autosome
txdb_exons_autosome <- keepSeqlevels(txdb_exons, autosome, pruning.mode = "coarse")
txdb_exons_autosome

# 8. D. rename in NCBI format
txdb_ncbi <- mapSeqlevels(seqlevels(txdb_exons), "NCBI")
txdb_ncbi
txdb_exons_ncbi <- renameSeqlevels(txdb_exons_autosome, txdb_ncbi)
txdb_exons_ncbi

dim(subsetByOverlaps(airway, txdb_exons_ncbi))
dim(subsetByOverlaps(airway, txdb_exons_ncbi))[1]  

# 9. The expression measures of the airway dataset are the number of reads mapping 
# to each feature. In the previous question we have established that many of these 
# features do not overlap autosomal transcripts from the TxDb.Hsapiens.UCSC.hg19.knownGene. 
# But how many reads map to features which overlaps these transcripts?
  
# For sample SRR1039508, how big a percentage (expressed as a number between 0 and 1) 
# of the total reads in the airway dataset for that sample, are part of a feature which 
# overlaps an autosomal TxDb.Hsapiens.UCSC.hg19.knownGene transcript?

# 9. A. How many reads map to features which overlaps these transcripts?
sample_SRR1039508 <- airway[, 1]
sample_SRR1039508
sample_SRR1039508_autosome <- subsetByOverlaps(sample_SRR1039508, txdb_exons_ncbi)
sample_SRR1039508_autosome

autosome_reads <- sum(assay(sample_SRR1039508_autosome, "counts"))
head(assay(sample_SRR1039508_autosome, "counts"))
autosome_reads
total_reads <- sum(assay(sample_SRR1039508, "counts"))
head(assay(sample_SRR1039508, "counts"))
total_reads

# 9. B. Percentage of the total reads in airway dataset for SRR1039508 which overlaps
# autosome of txdb
autosome_reads/total_reads

# 10. Consider sample SRR1039508 and only consider features which overlap autosomal 
# transcripts from TxDb.Hsapiens.UCSC.hg19.knownGene. We should be able to very roughly 
# divide these transcripts into expressed and non expressed transcript. Expressed 
# transcripts should be marked by H3K4me3 at their promoter. The airway dataset have 
# assayed “airway smooth muscle cells”. In the Roadmap Epigenomics data set, the E096 
# is supposed to be “lung”. Obtain the H3K4me3 narrowPeaks from the E096 sample using 
# the AnnotationHub package.

# What is the median number of counts per feature (for sample SRR1039508) containing 
# a H3K4me narrowPeak in their promoter (only features which overlap autosomal 
# transcripts from TxDb.Hsapiens.UCSC.hg19.knownGene are considered)?
  
# Clarification: We are using the standard 2.2kb default Bioconductor promotor setting.

# Conclusion Compare this to the median number of counts for features without a 
# H3K4me3 peak. Note that this short analysis has not taken transcript lengths 
# into account and it compares different genomic regions to each other; this is
# highly susceptible to bias such as sequence bias.

# 10. A. Selecting the "E096", "H3K4me3", "narrowPeak" records
library(AnnotationHub)
ah <- AnnotationHub()
ah
ah_E096 <- query(ah, c("E096", "H3K4me3", "narrowPeak"))
ah_E096
ah_record <- ah_E096[["AH30596"]]
ah_record

# 10. B. Selecting Autosomal transcripts (genes) for "E096", "H3K4me3", "narrowPeak" 
# and renaming by NCBI system
ah_record_autosome <- keepSeqlevels(ah_record, autosome, pruning.mode = "coarse")
ah_record_autosome
ah_record_ncbi <- renameSeqlevels(ah_record_autosome, txdb_ncbi)
ah_record_ncbi

# 10. C. Selecting Autosomes for sample SRR1039508
ncbi_group <- extractSeqlevelsByGroup(species = "Homo sapiens", style = "NCBI", group = "auto")
ncbi_group
sample_ncbi <- keepSeqlevels(range(rowRanges(sample_SRR1039508_autosome)), ncbi_group)
sample_ncbi

# 10. D. Finding the median no of counts per feature (gene) for sample SRR1039508 
# (using subset overlaps)
ov <- subsetByOverlaps(promoters(sample_ncbi), ah_record_ncbi)
ov
ov <- subsetByOverlaps(sample_SRR1039508, ov)
ov

median(assay(ov, "counts"))


# 10. E. Selecting the "E096", "narrowPeak" records
library(AnnotationHub)
ah <- AnnotationHub()
ah
ah_E096 <- query(ah, c("E096", "narrowPeak"))
ah_E096
ah_record <- ah_E096[["AH30596"]]
ah_record

# 10. F. Selecting Autosomal transcripts (genes) for "E096", "narrowPeak" 
# and renaming by NCBI system
ah_record_autosome <- keepSeqlevels(ah_record, autosome, pruning.mode = "coarse")
ah_record_autosome
ah_record_ncbi <- renameSeqlevels(ah_record_autosome, txdb_ncbi)
ah_record_ncbi

# 10. G. Selecting Autosomes for sample SRR1039508
ncbi_group <- extractSeqlevelsByGroup(species = "Homo sapiens", style = "NCBI", group = "auto")
ncbi_group
sample_ncbi <- keepSeqlevels(range(rowRanges(sample_SRR1039508_autosome)), ncbi_group)
sample_ncbi

# 10. H. Finding the median no of counts per feature (gene) for sample SRR1039508 
# (using subset overlaps)
ov <- subsetByOverlaps(promoters(sample_ncbi), ah_record_ncbi)
ov
ov <- subsetByOverlaps(sample_SRR1039508, ov)
ov

median(assay(ov, "counts"))





















