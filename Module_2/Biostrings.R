# Biostrings

# 1. Installing Biostrings Package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::valid()

BiocManager::install(c("BiocParallel", "Rsamtools"), update = TRUE, ask = FALSE, force = TRUE)

BiocManager::install("Biostrings")

# 2. Loading Biostrings Package
library(BiocGenerics)
library(S4Vectors)
library(IRanges)
library(GenomicRanges)
library(Biostrings)

# 3. Representing sequences
dna1 <- DNAString("ACGT-G")
dna1

# 4. DNA string set , a collection of DNA strings
dna2 <- DNAStringSet(c("ACGT", "GTCA", "GCTA"))
dna2

# 5. IUPAC code
IUPAC_CODE_MAP

# 6. Subset DNA gives a subsequence
dna1[2:4]

# 7. Subset a DNA String set -> another string set
dna2[1:2]

# 8. Get a string out of a DNA string set
dna2[[1]]

# 9. Put name on a DNA string set
names(dna2) = paste0("seq", 1:3)
dna2

# 10. Number of bases in the string
width(dna2)

# 11. Sort DNA string set
sort(dna2)

# 12. Reversing a DNA string set
rev(dna2)

# 13. Reversing a DNA string
rev(dna1)

# 14. Reversion of a DNA string set using the reverse() function, reverses each 
#     string in the DNA string set
reverse(dna2)

# 15. Use reverse compliment
reverseComplement(dna2)

# 16. Translate DNA to protein
translate(dna2)

# 17. Alphabet frequency gives a general frequency table of DNA bases
alphabetFrequency(dna2)

# 18. Counting GC content of a string
letterFrequency(dna2, letters = "GC")

# 19. Higher Order Nucleotide frequency e.g. dinucleotides
dinucleotideFrequency(dna2)

# 20. Consensus Matrix
consensusMatrix(dna2)













