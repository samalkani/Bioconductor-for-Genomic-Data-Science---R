# Short Read

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing "ShortRead"
BiocManager::install(c("ShortRead"), force = TRUE)

# 3. Troubleshooting installing packages
BiocManager::valid()

BiocManager::install(c(
  "BiocParallel", "promises", "Rsamtools"
), update = TRUE, ask = FALSE, force = TRUE)

# 4. Load the "ShortRead" package
library(ShortRead)

# 5. Pointer down to a FASTQ file
fastqDir <- system.file("extdata", "E-MTAB-1147", package = "ShortRead")
fastqPath <- list.files(fastqDir, pattern = ".fastq.gz$", full = TRUE)[1]
reads <- readFastq(fastqPath)
reads

# 6. Alternative instantiate a FASTQ file object
fqFile <- FastqFile(fastqPath)
fqFile

# 7. Reading the FASTQ file
reads <- readFastq(fqFile)

# 8. Accessing individual reads using sread
sread(reads)[1:2]

# 9. Accessing Quality values for each read
quality(reads)[1:2]

# 10. Accessing read names
id(reads)[1:2]

# 11. Converting (coercing) standard quality score into integers from 1 - 40
as(quality(reads), "matrix")[1:2,1:10]

