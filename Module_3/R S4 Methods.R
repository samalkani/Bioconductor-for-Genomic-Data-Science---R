# R S4 Methods

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing Genomic Ranges Packages"
# BiocManager::install("GenomicRanges", force = TRUE)

# 3. Load the Libraries
library(GenomicRanges)

# 4. Example of an S4 method (mimicMethod)
mimicMethod <- function(x) {
  if (is(x, "matrix"))
    method1(x)
  if (is(x, "data.frame"))
    method2(x)
  if (is(x, "IRanges"))
    method3(x)
}

# 5. Look at the printout of as.data.frame
as.data.frame

# 6. Base R printout of as.data.frame
base::as.data.frame
function (x, row.names = NULL, optional = FALSE, ...) 
{
  if (is.null(x)) 
    return(as.data.frame(list()))
  UseMethod("as.data.frame")
}

# 7. Using the showMethods function
showMethods("as.data.frame")

# 8. Displaying the function or method definition
getMethod("as.data.frame", "GenomicRanges")

# 9. Displaying the function or method definition (proper definition)
getMethod("as.data.frame", signature( x = "GenomicRanges"))

# 10. Accessing help for S4 methods or functions
method?"as.data.frame,DataFrame"

# 11. Alternative method for accessing help for S4 methods or functions
?"as.data.frame,DataFrame-method"

# 12. Alternative method for accessing help for S4 methods or functions
?"as.data.frame,GenomicRanges-method"

# 13. Using the showMethods function
showMethods("findOverlaps")

# 14. Displaying the function or method definition (proper definition)
# getMethod("findOverlaps", signature(query = "Ranges", subject = "Ranges"))
getMethod("findOverlaps", signature(query = "GenomicRanges", subject = "GenomicRanges"))

# 15. Getting help
?"findOverlaps,Ranges,Ranges-method"
?"findOverlaps,GenomicRanges,GenomicRanges-method"








