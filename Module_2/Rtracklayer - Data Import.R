# Rtracklayer - Data Import

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing rtracklayer and Rsamtools
BiocManager::install(c("rtracklayer", "Rsamtools"), force = TRUE)
BiocManager::install(c("Biostrings", "IRanges"), update = TRUE, ask = FALSE)

# 3. Troubleshooting installing packages
BiocManager::valid()

# 4. Load R libraries
library(rtracklayer)
library(AnnotationHub)
library(Rsamtools)

# 5. Help for Import function
?import

# 6. Set-up the Annotation Hub
ahub <- AnnotationHub()
ahub

# 7. Look at the R data class slot
table(ahub$rdataclass)

# 8. Look at a BIGWIG file on Annotation Hub
ahub.bw = subset(ahub, rdataclass == "BigWigFile" 
                 & species == "Homo sapiens")

# 9. Get the first dataset in this BIGWIG file
bw = ahub.bw[[1]]

# 10. Let's look at BW
bw

# 11. Import the data into memory
import(bw)

# 12. Import (read in) part of the file into memory
gr.chr22 = import(bw, which=GRanges("chr22", ranges = IRanges(1, 10^8)))
gr.chr22

# 13. Import (read in) part of the file into memory and return it as a RLE
# (Run Length Encoded Vector)
rle.chr22 = import(bw, which=GRanges("chr22", ranges = IRanges(1, 10^8)), as = "Rle")
rle.chr22

# 14. Display only RLE for chromosome 22
rle.chr22$chr22

# 15. Looking at the lift over function but to use lift over you need a chain file
# Accessing a chain file on Annotation Hub
ahub.chain = subset(ahub, rdataclass == "ChainFile")
ahub.chain

# 16. Subset further to only contain Human data
ahub.chain = subset(ahub.chain, species == "Homo sapiens")
ahub.chain

# 17. Lift over query for hg18 and hg19
query(ahub.chain, c("hg18","hg19"))

# 18. If we want to convert hg19 into hg18
chain = query(ahub.chain, c("hg18","hg19"))[[1]]

# 19. A set of GRanges
gr.chr22 = import(bw, which=GRanges("chr22", ranges = IRanges(1, 10^8)))
gr.chr22

# 20. Using the liftOver() function
gr.hg18 = liftOver(gr.chr22, chain)

# 21. What comes out of the liftOver function
class(gr.hg18)

# 22. Length of the GRanges list (gr.hg18) = the GRanges list we lifted over (gr.chr22)
length(gr.hg18)

# 23. How often are the GRanges split when converting from hg19 to hg18?
# Warning elementLengths() function out of date, alternative to elementLengths()
elementLengths=elementNROWS

# Then do
table(elementLengths(gr.hg18))

















