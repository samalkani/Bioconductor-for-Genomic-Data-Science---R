# biomaRt

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing "biomaRt"
BiocManager::install(c("biomaRt"), force = TRUE)

# 3. Troubleshooting installing packages
BiocManager::valid()

# 4. Load Libraries
library(GenomicRanges)
library(biomaRt)
library(Biobase)
library(GenomicFiles)

# 5. Listing a marts
head(listMarts())

# 6. Picking Ensembl database
mart <- useMart("ensembl")
mart

# 7. List datasets in the database
head(listDatasets(mart))

# 8. Selecting a dataset
ensembl <- useDataset("hsapiens_gene_ensembl", mart)
ensembl

# 9. Building a query
values <- c("202763_at","209310_s_at","207500_at")
getBM(attributes = c("ensembl_gene_id", "affy_hg_u133_plus_2"),
                     filters = "affy_hg_u133_plus_2", values = values, mart = ensembl)

# 10. Calling list attributes on Ensembl database
attributes <- listAttributes(ensembl)
nrow(attributes)

# 11. Getting names of the attributes
head(attributes)

# 12. Getting names of the attributes (bottom)
tail(attributes)

# 13. Getting names of the attributes (last 100) output not shown
tail(attributes, n=100)

# 14. Getting names of the attributes (last 500) output not shown
tail(attributes, n=500)

# 15. Listing filters
filters <- listFilters(ensembl)
head(filters)

# 16. Number of filters
nrow(filters)

# 17. Attribute pages
attributePages(ensembl)

# 18. Listing Attributes without homologs
attributes <- listAttributes(ensembl, page = "feature_page")
head(attributes)




