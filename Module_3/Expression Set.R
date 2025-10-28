# Expression Set

# 1. install ALL package

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

BiocManager::install("ALL")

# 2. Load ALL package
library(ALL)

# 3. Print, display ALL package
data(ALL)
ALL

# 4. Details of Experiment given using an accessor function
experimentData(ALL)

# 5. More information can be obtained looking at the help page
?ALL

# 6. Accessing the expression matrix
exprs(ALL)[1:4, 1:4]

# 7. Details of Sample names by the sample names accessor function
head(sampleNames(ALL))

# 8.  Details of Feature names by the feature names accessor function
head(featureNames(ALL))

# 9. Accessing phenotype data from the data frame
head(pData(ALL))

# 10. Accessing specific covariants for modelling purposes
pData(ALL)$sex

# 11. Alternative
ALL$sex

# 12. Selecting the first 5 samples
ALL[, 1:5]

# 13. Selecting the first 10 features
ALL[1:10,]

# 14. Selecting both features and samples
ALL[1:10, 1:5]

# 15. Obtaining Feature data (normally not stored in the object itself)
featureData(ALL)

# 16. Getting the first 5 feature names
ids = featureNames(ALL)[1:5]
ids

# 17. Mapping the Affymetrix ID's into something specific

# 17. A. Install hgu95av2.db package
BiocManager::install("hgu95av2.db")

# 17. B. Loading hgu95av2.db package
library(hgu95av2.db)

# 17. C. Mapping Affymetrix ID's into ENTREZID
as.list(hgu95av2ENTREZID[ids])

# 18. Using PhenoData instead of pData (the recommended one!)
phenoData(ALL)

# 19. Better example
names(pData(ALL))

# 20. VarLabels
varLabels(ALL)

# 21. PhenoData is different to pData
phenoData(ALL)

# vs
head(pData(ALL))

# 22. PhenoData contains a pData slot
head(pData(phenoData(ALL)))









