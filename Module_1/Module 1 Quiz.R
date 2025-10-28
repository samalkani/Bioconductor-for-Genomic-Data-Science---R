# Module 1 Quiz

# 1.	Use the Annotation Hub package to obtain data on “CpG Islands” in the 
#     human genome. How many islands exists on the autosomes?

      # A. Installing the AnnotationHub package
      # BiocManager::install("AnnotationHub")

      # B. Load AnnotationHub
        library(AnnotationHub)

      # C. Retrieve the cpg
      ah <- AnnotationHub()
      ah_human <- subset(ah, species == "Homo sapiens")
      ah_human
      ah_human_cpg <- query(ah_human, "CpG Islands")
      ah_human_cpg
      ah_human_cpg_record <- ah_human_cpg[["AH5086"]]
      ah_human_cpg_record

      # D. Extract autosomes
      filter <- c(paste("chr", 1:22, sep=""))
      filter
      split_record <- split(ah_human_cpg_record, seqnames(ah_human_cpg_record))
      split_record
      autosomes <- split_record[filter]
      autosomes
      
      # E. Check the number of autosomes
      unlist(autosomes)
      
# 2.	How many CpG Islands exists on chromosome 4.
      
      # A Print Autosomes
      autosomes[4]
      
      # or B
      autosomes$chr4

# 3.	Obtain the data for the H3K4me3 histone modification for the H1 cell line
#     from Epigenomics roadmap, using AnnotationHub. Subset these regions to 
#     only keep regions mapped to the autosomes (chromosomes 1 to 22). How many
#     bases does these regions cover?
      
      # A. Retrieve the record
      ah_H3K4me3 <- query(ah, c("H3K4me3", "narrowpeak", "E003"))
      ah_H3K4me3
      
      # B. Select the narrowpeak
      ah_H3K4me3_record <- ah_H3K4me3[["AH29884"]]
      ah_H3K4me3_record
      
      # C. Extract autosomes and check the number of regions cover
      split_H3K4me3 <- split(ah_H3K4me3_record, seqnames(ah_H3K4me3_record))
      split_H3K4me3
      H3K4me3_autosomes <- split_H3K4me3[filter]
      H3K4me3_autosomes
      sum(width(unlist(H3K4me3_autosomes)))
      
# 4.	Obtain the data for the H3K27me3 histone modification for the H1 cell line
#     from Epigenomics Roadmap, using the AnnotationHub package. Subset these 
#     regions to only keep regions mapped to the autosomes. In the return data, 
#     each region has an associated "signalValue".
      
#     Question: What is the mean signalValue across all regions on the standard
#     chromosomes?
        
      # A. Retrieve the record
      ah_H3K27me3 <- query(ah, c("H3K27me3", "narrowpeak", "E003"))
      ah_H3K27me3
      ah_H3K27me3_record <- ah_H3K27me3[["AH29892"]]
      ah_H3K27me3_record
      
      # B. Extract autosomes
      split_H3K27me3 <- split(ah_H3K27me3_record, seqnames(ah_H3K27me3_record))
      split_H3K27me3
      H3K27me3_autosomes <- split_H3K27me3[filter]
      H3K27me3_autosomes
      
      # C. Create a subset of extracted autosomes
      ah_H3K27me3_autosomes <- subset(ah_H3K27me3_record, seqnames %in% filter)
      ah_H3K27me3_autosomes
      
      # D. Mean signalValue
      mean_signalValue <- mean(ah_H3K27me3_autosomes$signalValue)
      mean_signalValue
      
# 5.	Bivalent regions are bound by both H3K4me3 and H3K27me3. Using the regions
#     we have obtained above, how many bases on the standard chromosomes are 
#     bivalently marked?
      
      # A. Intersect between two records
      bivalent <- intersect(unlist(H3K4me3_autosomes), unlist(H3K27me3_autosomes))
      bivalent
      sum(width(bivalent))

# 6.	We will examine the extent to which bivalent regions overlap CpG Islands, 
#     how big a fraction (expressed as a number between 0 and 1) of the bivalent
#     regions, overlap one or more CpG Islands
      
      # A. Find bivalent regions overlap CpG Islands
      cpg_autosomes <- autosomes
      cpg_autosomes
      cpg_bivalent <- findOverlaps(bivalent, unlist(cpg_autosomes))
      cpg_bivalent
      
      # B. Calculate the fraction of the bivalent regions overlap CpG Islands
      fraction_bivalent <- length(unique(queryHits(cpg_bivalent)))/length(bivalent)
      fraction_bivalent
      
# 7.	How big a fraction (expressed as a number between 0 and 1) of the bases 
#     which are part of CpG Islands, are also bivalent marked?
      
      # A. Intersect between bivalent regions and CpG Islands
      cpg_bivalent_intersect <- intersect(bivalent, unlist(cpg_autosomes))
      cpg_bivalent_intersect
      
      # B. Calculate the fraction of the bases intersected between CpG Islands and bivalent
      fraction_bivalent_intersect <- 
        sum(width(reduce(cpg_bivalent_intersect)))/sum(width(unlist(cpg_autosomes)))
      fraction_bivalent_intersect
      
# 8.	How many bases are bivalently marked within 10kb of CpG Islands? Tip: 
#     consider using the "resize()"" function.
      
      # A. Extract CpG Islands within 10kb
      cpg_10kb <- 
        resize(unlist(cpg_autosomes), width = 20000 + width(unlist(cpg_autosomes)), 
               fix = "center")
      cpg_10kb
      cpg_10kb_bivalent <- intersect(cpg_10kb, bivalent)
      cpg_10kb_bivalent
      sum(width(cpg_10kb_bivalent))
      
# 9.	How big a fraction (expressed as a number between 0 and 1) of the human 
#     genome is contained in a CpG Island? Tip 1: the object returned by 
#     AnnotationHub contains "seqlengths". Tip 2: you may encounter an integer 
#     overflow. As described in the session on R Basic Types, you can address 
#      this by converting integers to numeric before summing them, "as.numeric()".
      
      # A. Calculate human genome size
      chr_list <- c(paste("chr", 1:22, sep=""))
      chr_list
      genome <- keepSeqlevels(ah_human_cpg_record, chr_list, pruning.mode = "coarse")
      genome
      genome_size <- sum(as.numeric(seqlengths(genome)))
      genome_size
      
      # B. Calculate the fraction of human genome which contained a CpG Island
      cpg_autosomes_size <- sum(as.numeric(width(unlist(cpg_autosomes))))
      cpg_autosomes_size / genome_size
      
# 10.	Compute an odds-ratio for the overlap of bivalent marks with CpG islands.
      
      ## A. Calculate InOut matrix
      overlapMat <- matrix(0,, ncol = 2, nrow = 2)
      overlapMat
      colnames(overlapMat) <- c("in", "out")
      rownames(overlapMat) <- c("in", "out")
      overlapMat[1,1] <- sum(width(cpg_bivalent_intersect))
      overlapMat[1,1]
      overlapMat[1,2] <- sum(width(setdiff(bivalent, unlist(cpg_autosomes))))
      overlapMat[1,2]
      overlapMat[2,1] <- sum(width(setdiff(unlist(cpg_autosomes), bivalent)))
      overlapMat[2,1]
      overlapMat[2,2] <- genome_size - sum(overlapMat)
      overlapMat[2,2]
      overlapMat
      
      ## B. Calculate odds-ratio
      oddsRatio <- overlapMat[1,1] * overlapMat[2,2] / (overlapMat[2,1] * overlapMat[1,2])
      oddsRatio
      
      
          