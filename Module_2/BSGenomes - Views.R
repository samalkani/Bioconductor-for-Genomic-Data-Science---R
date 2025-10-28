# BSGenome - Views

# 1. Load the package
library(Biostrings)

# 2. Loading a specific yeast genome
library(BSgenome.Scerevisiae.UCSC.sacCer2)

# 3. Call the short name of the genome object
Scerevisiae

# 4. Running a match pattern of a DNA string on the yeast genome chromosome 1
dnaseq <- DNAString("ACGTACGT")
vi = matchPattern(dnaseq, Scerevisiae$chrI)
vi

# 5. Get the IRanges out of the views object
ranges(vi)

# 6. Check the coordinates actually gets us the right nucleotides
Scerevisiae$chrI[57932:57939]

# 7. Run functions on the views object as if it were a DNA string set
alphabetFrequency(vi)

# 8. We can shift the view by 10 bases
shift(vi, 10)

# 9. Running vmatchPattern on the entire genome
gr = vmatchPattern(dnaseq, Scerevisiae)
gr

# 10. Instantiate a view object with the views function and the GRanges (gr)
vi2 = Views(Scerevisiae, gr)
vi2

# 11. Compute the GC content of promoters in the yeast genome

# 11. A. Load Annotation Hub
library(AnnotationHub)

# 11. B. Creating a local annotation hub
ahub = AnnotationHub()

# 11. C. Call Annotation Hub
ahub

# 11. D. Query
qh = query(ahub, c("sacCer2", "genes"))
qh

# 11. E. Retrieve the records
genes = ahub[["AH7048"]]

# 11. F. Getting the promoters
prom = promoters(genes)

# 11. G. Let’s look at element number 2 in the promoter
prom

# 11. H. Cutting off anything outside the sequence length of the genome.
prom = trim(prom)
prom

# 11. I. Instantiate a view of these promoters
promViews = Views(Scerevisiae, prom)
promViews

# 11. J. Get GC percentage for the promoters and density plot
gcProm = letterFrequency(promViews, "GC", as.prob = TRUE)
head(gcProm)
tail(gcProm)
plot(density(gcProm))
abline(v = 0.38)
       




