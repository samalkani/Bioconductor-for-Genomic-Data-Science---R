# Usecase - Annotation Hub  and Granges, Part 1

# 1. Installing the AnnotationHub package
# BiocManager::install("AnnotationHub")

# 2. Load AnnotationHub
library(AnnotationHub)

# 3. Instantiate an Annotation Hub
ahub = AnnotationHub()

# 4. Subset the Annotation Hub to be specific to Human data
ahub = subset(ahub, species == "Homo sapiens")

# 5. Search for data under two search terms using query function
qhs = query(ahub, c("H3K4me3", "Gm12878"))
qhs

# 6. Get two genome ranges the second element and the fourth element
# 2nd element Rep1, broadPeak = gr1
# 4th element Rep1, narrowPeak = gr2
gr1 = qhs[[2]]
gr2 = qhs[[4]]

# 7. 2nd Element
gr1

# 8. Summary of 2nd element (broad peak)
summary(width(gr1))

# 9. Summary of 4th element (narrow peak)
summary(width(gr2))

# 10. Create a table of all the widths for the narrow peak
table(width(gr2))

# 11. Rename the GRanges (gr2) to peaks
peaks = gr2

# 12. Re-list the genomes
qhs

# 13. Display the Rep1, narrow peak genome
qhs[4]

# 14. Performing a query on the RefSeq database
qhs = query(ahub, "RefSeq")
qhs

# 15. Looked at the genome (clarify why we have 8 data sets)
qhs$genome

# 16. Download the first data set from ahub the hg19 gene
genes = qhs[[1]]
genes

# 17. Use a table to identify how often a gene occurs, and nesting the table function within a  
#     table function, to see how often a single or two transcripts occur
table(table(genes$name))

# 18. Call promoters on your data
prom = promoters(genes)

# 19. Define the width of a promoter
table(width(prom))

# 20. Look at the arguments for promoters
args(promoters)
function (x, upstream = 2000, downstream = 200, ...) 
  NULL











