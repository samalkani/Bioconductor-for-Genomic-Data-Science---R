# Count-based RNA-seq Analysis

# 1. Select an example data set and display its contents
library(airway)
data(airway)
airway

# 2. Counts of features tells us something about the expression level
assay(airway, "counts")[1:3, 1:3]

# 3. Variable of interest - dex
airway$dex

# 4. Change the reference level so the treatment level is NOT the reference level
airway$dex <- relevel(airway$dex, "untrt")
airway$dex

# 5. There is rich information about which gene model was used for each gene:
granges(airway)

# 6. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 7. Installing "DESeq2", "edgeR"
BiocManager::install(c("DESeq2", "edgeR"), force = TRUE)

# 8. Troubleshooting installing packages
BiocManager::valid()

BiocManager::install(c(
  "BiocParallel", "DESeq2", "oligo", "promises", "Rsamtools"
), update = TRUE, ask = FALSE, force = TRUE)

# 9. Load Libraries
library(DESeq2)
library(edgeR)

# 10. Constructing a DEGlist (Limma class which is compatible with edgeR)
dge <- DGEList(counts = assay(airway, "counts"),
                               group = airway$dex)
dge$samples 

# 11. Merge the existing data frame with the DGE samples with the call data from the airway    
# dataset
dge$samples <- merge(dge$samples, as.data.frame(colData(airway)), by = 0)
dge$samples

# 12. No information in the DGE list on genes
head(dge$genes)

# 13. Take the gene names from row ranges in the airway dataset and put into the dataframe and gene slot
dge$genes <- data.frame(name = names(rowRanges(airway)), stringsAsFactors = FALSE)
head(dge$genes)

# 14. First we estimate the normalization factors or effective library sizes
dge <- calcNormFactors(dge)
dge$samples

# 15. Estimating the variability of the data with;

# 15. A. Common dispersion
dge <- estimateGLMCommonDisp(dge)

# 15. B. Tag-wise dispersion
dge <- estimateGLMTagwiseDisp(dge)

# 16. Design Matrix
design <- model.matrix(~dge$samples$group)
design

# 17. Now we do a glmFit(), similar to limma
fit <- glmFit(dge, design)

# 18. Examining which genes are most significant (using GLM likelihood ratio test)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

# 19. Using DESeq package for for analysing differential gene expression

# 19. A. Putting data into a DESeq2 container (like summarized experiment)
dds <- DESeqDataSet(airway, design = ~ dex)

# 19. B. Fitting the model is simple
dds <- DESeq(dds)

# 20. Getting the results of the fitted model
res <- results(dds)
res

# 21. Sort results according to adjusted P value
res <- res[order(res$padj),]
res

















