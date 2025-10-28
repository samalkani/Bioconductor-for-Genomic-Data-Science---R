# Usecase - Annotation Hub  and Granges, Part 2

# 1. Installing the AnnotationHub package
# BiocManager::install("AnnotationHub")

# 2. Load AnnotationHub
library(AnnotationHub)

# 3. Instantiate an Annotation Hub
ahub = AnnotationHub()

# 4. Performing a query on the RefSeq database
qhs = query(ahub, "RefSeq")
qhs

# 5. Download the first data set from ahub the hg19 gene
genes = qhs[[1]]
genes

# 6. Quick Recap - Obtained promoters (Call promoters on your data)
prom = promoters(genes)
prom

# 7. Subset the Annotation Hub to be specific to Human data
ahub = subset(ahub, species == "Homo sapiens")

# 8. Search for data under two search terms using query function
qhs = query(ahub, c("H3K4me3", "Gm12878"))
qhs

# 9. Get two genome ranges the second element and the fourth element
# 2nd element Rep1, broadPeak = gr1
# 4th element Rep1, narrowPeak = gr2
gr1 = qhs[[2]]
gr2 = qhs[[4]]

# 10. Rename the GRanges (gr2) to peaks
peaks = gr2

# 11. Quick - Recap - Histone Modification Peaks
peaks

# 12. Use findOverlaps function to determine whether this particular histone 
#     modification H3K4 trimethylation is enriched in promoters
findOverlaps(prom, peaks)

# 13. Find the number of promoters that have a peak in them avoid double counting 
#     by using the unique function
ov = findOverlaps(prom, peaks)
length(unique(queryHits(ov)))

# 14. Find the number of peaks that have a promoter in them
length(unique(subjectHits(ov)))

# 15. Out of peaks, how many overlap a promoter?
length(subsetByOverlaps(peaks, prom, ignore.strand = TRUE))

# 16. Percentage of peaks that overlap a promoter
length(subsetByOverlaps(peaks, prom, ignore.strand = TRUE)) / length(peaks)

# 17. Percentage of promoters that overlap a peak
length(subsetByOverlaps(prom, peaks, ignore.strand = TRUE)) / length(prom)

# 18. How many bases (in mega-bases) do the peaks really cover?
sum(width(reduce(peaks, ignore.strand = TRUE))) / 10^6

# 19. How many bases (in mega-bases) do the promoters really cover?
sum(width(reduce(prom, ignore.strand = TRUE))) / 10^6

# 20. Intersection of peaks and promoters ignoring strand
sum(width(intersect(peaks, prom, ignore.strand = TRUE))) / 10^6

# 21. Construct a 2 x 2 matrix of zeros
inOut = matrix(0, ncol = 2, nrow = 2)
colnames(inOut) = c("in","out")
rownames(inOut) = c("in","out")
inOut

# 22. Fill the matrix (bases in both peaks and promoters)
inOut[1,1] = sum(width(intersect(peaks, prom, ignore.strand = TRUE)))

# 23. Fill the matrix (bases in peaks and not in promoters)
inOut[1,2] = sum(width(setdiff(peaks, prom, ignore.strand = TRUE)))

# 24. Fill the matrix (bases in promoters and not in peaks)
inOut[2,1] = sum(width(setdiff(prom, peaks, ignore.strand = TRUE)))

# 25. Print Matrix
inOut

# 26. Do Column sums
colSums(inOut)

# 27. Do Row sums
rowSums(inOut)

# 28. Subtract the number of bases in the Sum of inOut from the total number
# of bases in the human genome 3 billion
inOut[2,2] = 3*10^9 - sum(inOut)
inOut

# 29. Use the Fisher Exact test on the completed matrix
# But you get an error, because the integer is > the biggest integer
fisher.test(inOut)$statistic

# 30. Calculate the odds ratio by hand
oddsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2])
oddsRatio

# 31. Assume human genome is only 1.5 billion bases (mappable)
inOut[2,2] = 0
inOut[2,2] = 1.5*10^9 - sum(inOut)
inOut

# 32. Calculate the odds ratio by hand
oddsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2])
oddsRatio














