# The Limma Package

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing "limma, leukemiasEset"
BiocManager::install(c("limma", "leukemiasEset"), force = TRUE)

# 3. Troubleshooting installing packages
BiocManager::valid()

BiocManager::install(c(
  "BiocParallel", "oligo", "promises", "reticulate", "Rsamtools"
), update = TRUE, ask = FALSE, force = TRUE)

# 4. Load Libraries
library(limma)
library(leukemiasEset)

# 5. Displaying the leukemiasEset expression set
data(leukemiasEset)
leukemiasEset

# 6. Profiling the subsets of Leukemia
table(leukemiasEset$LeukemiaType)

# 7. Subset the Leukemia types (select ALL & NoL)
ourData <- leukemiasEset[, leukemiasEset$LeukemiaType %in% c("ALL", "NoL")]

# 8. Leukemia type is a factor
ourData$LeukemiaType

# 9. Getting rid of the three redundant factors in leukemia type not in the subset
ourData$LeukemiaType <- factor(ourData$LeukemiaType)
ourData$LeukemiaType

# 10. Setting up the design matrix
design <- model.matrix(~ ourData$LeukemiaType)
head(design)

# 11. Fitting a linear model using (limma) to all the genes separately
fit <- lmFit(ourData, design)

# 12. Fit the model
fit <- eBayes(fit)

# 13. Look at the top of the table
topTable(fit)

# 14. Getting rid of the three redundant factors in leukemia type not in the subset
ourData$LeukemiaType <- factor(ourData$LeukemiaType)
ourData$LeukemiaType

# 15. Get the first gene from the table
topTable(fit, n = 1)

# 16. Get the gene name
genename <- rownames(topTable(fit, n=1))
genename

# 17. Compute the mean expression inside each group
typeMean <- tapply(exprs(ourData)[genename,], ourData$LeukemiaType, mean)
typeMean

# 18. The log fold change between ALL and NoL
typeMean["NoL"] - typeMean["ALL"]

# 19. Create design matrix 2
design2 <- model.matrix(~ ourData$LeukemiaType - 1)
head(design2)

# 20. Changing the names of the columns to shorter names
colnames(design2) <- c("ALL", "NoL")
head(design2)

# 21. Fit the model (as before)
fit2 <- lmFit(ourData, design2)

# 22. Make contrasts
contrast.matrix <- makeContrasts("ALL-NoL", levels = design2)
contrast.matrix

# 23. Execute contrasts.fit
fit2C <- contrasts.fit(fit2, contrast.matrix)

# 24. Use eBayes command to borrow information across genes
fit2C <- eBayes(fit2C)

# 25. Display results
topTable(fit2C)



















