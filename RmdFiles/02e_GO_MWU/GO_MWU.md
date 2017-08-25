GO\_MWU uses continuous measure of significance (such as fold-change or
-log(p-value) ) to identify GO categories that are significantly
enriches with either up- or down-regulated genes. The advantage - no
need to impose arbitrary significance cutoff.

If the measure is binary (0 or 1) the script will perform a typical "GO
enrichment" analysis based Fisher's exact test: it will show GO
categories over-represented among the genes that have 1 as their
measure.

On the plot, different fonts are used to indicate significance and color
indicates enrichment with either up (red) or down (blue) regulated
genes. No colors are shown for binary measure analysis.

The tree on the plot is hierarchical clustering of GO categories based
on shared genes. Categories with no branch length between them are
subsets of each other.

The fraction next to GO category name indicates the fracton of "good"
genes in it; "good" genes being the ones exceeding the arbitrary
absValue cutoff (option in gomwuPlot). For Fisher's based test, specify
absValue=0.5. This value does not affect statistics and is used for
plotting only.

Stretch the plot manually to match tree to text

Mikhail V. Matz, UT Austin, February 2015; <matz@utexas.edu>

################################################################ 

NOTES: This program drains memory and creates some very large
intermediate files, especially for the biological process catagory.

First, I run the stats from the command line to make sure its working.
Once I've generated the temp files, I comment out then stats portions
and recreate the plots by kniting the rmd file.

    library(ape)

    ## Warning: package 'ape' was built under R version 3.3.2

    source("gomwu.functions.R")

    # set output file for figures 
    knitr::opts_chunk$set(fig.path = '../../figures/02_RNAseq/')

From Experitment 1: Dissociation Test Molecular Function (MF)
-------------------------------------------------------------

    # input files
    input="padjAPAConsistentControlDG.csv" 
    goAnnotations="goAnnotations.tab" 
    goDatabase="go.obo" 
    goDivision="MF" # either MF, or BP, or CC

    # Calculating stats
    #gomwuStats(input, goDatabase, goAnnotations, goDivision,
    #   perlPath="perl", 
    #   largest=0.1,  
    #   smallest=5,   
    #   clusterCutHeight=0.25)  

    # Data viz
    gomwuPlot(input,goAnnotations,goDivision,
        absValue=-log(0.05,10),  
        level1=0.01, 
        level2=0.05, 
        level3=0.001, 
        txtsize=1.4,    
        treeHeight=0.5, 
      colors=c("grey73","lightsalmon","grey73","lightsalmon") 
    )

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

![](../../figures/02_RNAseq/GOConsistentControlDG-1.png)

    ## GO terms dispayed:  20 
    ## "Good genes" accounted for:  162 out of 429 ( 38% )

    # input files
    input="padjAPAConflictControlDG.csv" 
    goAnnotations="goAnnotations.tab" 
    goDatabase="go.obo" 
    goDivision="MF" # either MF, or BP, or CC

    # Calculating stats
    #gomwuStats(input, goDatabase, goAnnotations, goDivision,
    #   perlPath="perl", 
    #   largest=0.1,  
    #   smallest=5,   
    #   clusterCutHeight=0.25)  

    # Data viz
    gomwuPlot(input,goAnnotations,goDivision,
        absValue=-log(0.05,10),  
        level1=0.01, 
        level2=0.05, 
        level3=0.001, 
        txtsize=1.4,    
        treeHeight=0.5, 
      colors=c("grey73","red","grey73","red") 
    )

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

![](../../figures/02_RNAseq/GOConflictControlDG-1.png)

    ## GO terms dispayed:  5 
    ## "Good genes" accounted for:  47 out of 321 ( 15% )

    #library(dplyr)
    #list of all genes. I only care about column 1
    #allgenes <- read.csv("padjAPAConsistentControlDG.csv", header = T)
    #allgenes <- select(allgenes, gene)
    # first, identify genes differentially expressed by region in all experiments
    #intersection <- read.csv/../data/20c_vennDGControlminusConsistent.csv", header=F)
    #names(intersection)[1] <- "gene"      # rename columns
    #intersection$PresAbs <- 1
    #intersection <- full_join(allgenes, intersection)
    #intersection <- intersection %>%
    #   mutate(PresAbs = replace(PresAbs,is.na(PresAbs),0))
    #str(intersection)
    #write.csv(intersection, "./ConsistenMinusConflict.csv", row.names = F)

    # input files
    input="ConsistenMinusConflict.csv" 
    goAnnotations="goAnnotations.tab" 
    goDatabase="go.obo" 
    goDivision="MF" # either MF, or BP, or CC

    # Calculating stats
    #gomwuStats(input, goDatabase, goAnnotations, goDivision,
    #   perlPath="perl", 
    #   largest=0.1,  
    #   smallest=5,   
    #   clusterCutHeight=0.25,
    #   Alternative="g"
    #)

    gomwuPlot(input,goAnnotations,goDivision,
        absValue=0.001,   
        level1=0.05, 
        level2=0.05, 
        level3=0.05,
        txtsize=1.5,    
        treeHeight=0.5 
    )

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab =
    ## "", : the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

![](../../figures/02_RNAseq/GOConsistenMinusConflict-1.png)

    ## GO terms dispayed:  5 
    ## "Good genes" accounted for:  16 out of 51 ( 31% )
