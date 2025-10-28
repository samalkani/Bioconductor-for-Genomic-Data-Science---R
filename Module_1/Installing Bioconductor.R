# Installing Bioconductor 3.21 for R version 4.5.1

if (!requireNamespace("BiocManager", quietly = TRUE)) 
  
  install.packages("BiocManager")

BiocManager::install(version = "3.21")

# Installing the Limma package

BiocManager::install("limma")

# Checking version of Bioconductor available on R Studio

BiocManager::version()

# Checking whether Bioconductor installation is updated 
# and the installed packages are from the same version of Bioconductor

BiocManager::valid()
