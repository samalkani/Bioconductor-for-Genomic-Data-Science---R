# Biostrings - Matching

# 1. Installing Biostrings Package
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::valid()

# BiocManager::install(c("BiocParallel", "Rsamtools"), update = TRUE, ask = FALSE, force = TRUE)

# BiocManager::install("Biostrings")

# 2. Loading Biostrings Package
library(BiocGenerics)
library(S4Vectors)
library(IRanges)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Scerevisiae.UCSC.sacCer2)

# 3. Getting a Genome and a string
dnaseq <- DNAString("ACGTACGT")
dnaseq

# 4. Matching (mapping) a single string to a single string
matchPattern(dnaseq, Scerevisiae$chrI)

# 5. Count Matches
countPattern(dnaseq, Scerevisiae$chrI)

# 6. Matching a single string against many sequences (chromosomes)
vmatchPattern(dnaseq, Scerevisiae)

# 7. Is dnaseq the reverse compliment of itself?
dnaseq == reverseComplement(dnaseq)

# 8. MatchPVM() function
# MatchPVM - precision weight matrix is also known as a sequent logo, or a transcription 
# factor binding a motif, and details a probabilistic representation of a short 
# sequence.precision weight matrix is also known as a sequent logo, or a 
# transcription factor binding a motif, and details a probabilistic 
# representation of a short sequence. 

# MatchPVM allows us to search the genome for example for binding site for given
# transcription factor

# 9. pairwiseAlignment() function
# Implements a classic pairwiseAlignment algorithm, that's known in computations 
# biology, either a global pairwiseAlignment or a local pairwiseAlignment like a 
# Smith Waterman or Needleman Bush, a type algorithm
#
# It allows to map millions of reads against a short sequence such as a gene

# 10. TrimLRPatterns() functions
# Allows for trimming off specific patterns on the left and the right of a DNA string set
#
# The use case here is trimming off sequence adapters, but trimLRPatterns has a very
# rich set of functionality, allowing indels and mismatches in the sequence adaptors




