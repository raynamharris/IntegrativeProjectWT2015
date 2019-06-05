RNAseq gene expression analysis with DESeq2
-------------------------------------------

This workflow was modified from the DESeq2 tutorial found at:
<https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf>

First I load a handful of packages for data wrangling, gene expression
analysis, data visualization, and statistics.

    library(dplyr) ## for filtering and selecting rows
    library(plyr) ## for renmaing factors
    library(reshape2) ## for melting dataframe
    library(DESeq2) ## for gene expression analysis
    library(edgeR)  ## for basic read counts status
    library(magrittr) ## to use the weird pipe
    library(genefilter)  ## for PCA fuction
    library(ggplot2)

    ## Functions
    source("functions_RNAseq.R")

    ## set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/02a_makedfs/')

Now, I create data frames from three csv files - count: Contains counts
for all transcripts generated from the program Kallisto. This data can
be reproducibed from the file kallisto.Rmd - geneids: Contains the
ensemble ids and gene names for all the transcripts in the counts data
frame. This file will be used to convert transcipt counts to gene
counts. This file was also created via kallisto.Rmd file - colData: This
file contains all the information I collected for each sample that was
sequenced. Not all columns will be needed, so some are removed later.

    countData <- read.csv("../data/00_CountData.csv", row.names=1, check.names=FALSE )
    colData <- read.csv("../data/IntegrativeWT2015ColData.csv")

In this next section, I tidy the trait data for each sample so that I
can calculate differential gene expression for the colData of interest.
I also remove some samples for reasons described within the code blocks.

    rownames(colData) <- colData$RNAseqID    # set $genoAPAsessionInd as rownames
    colData <- colData[c(1,2,5,7:11)]  #keeping informative volumns
    head(colData)

    ##              RNAseqID   Mouse Region      Group     APA   Conflict
    ## 143A-CA3-1 143A-CA3-1 15-143A    CA3   conflict Trained   Conflict
    ## 143A-DG-1   143A-DG-1 15-143A     DG   conflict Trained   Conflict
    ## 143B-CA1-1 143B-CA1-1 15-143B    CA1    control   Yoked   Conflict
    ## 143B-DG-1   143B-DG-1 15-143B     DG    control   Yoked   Conflict
    ## 143C-CA1-1 143C-CA1-1 15-143C    CA1 consistent Trained NoConflict
    ## 143D-CA1-3 143D-CA1-3 15-143D    CA1    control   Yoked NoConflict
    ##                  APA_Conflict Treatment
    ## 143A-CA3-1   Trained_Conflict  conflict
    ## 143A-DG-1    Trained_Conflict  conflict
    ## 143B-CA1-1     Yoked_Conflict   shocked
    ## 143B-DG-1      Yoked_Conflict   shocked
    ## 143C-CA1-1 Trained_NoConflict   trained
    ## 143D-CA1-3   Yoked_NoConflict     yoked

    levels(colData$APA_Conflict)

    ## [1] "NA_NA"              "Trained_Conflict"   "Trained_NoConflict"
    ## [4] "Yoked_Conflict"     "Yoked_NoConflict"

    # changing the analysis to include homecage
    # colData <- colData %>% dplyr::filter(!grepl("147-|148-", RNAseqID))  # remove 147, and 148:  homecage animals 
    colData$ID <- gsub("[[:punct:]]", "", colData$Mouse) #make a column that thas id without the dash
    colData$APA <- NULL ## delete old APA column
    names(colData)[names(colData)=="APA_Conflict"] <- "APA3" #rename  APA3 to match color scheme
    names(colData)[names(colData)=="Region"] <- "Punch" #rename  region to punch

    # rename factors & group all Control animals into 1 group
    colData$APA2 <- colData$APA3 
    colData$APA2 <- revalue(colData$APA2, c("NA_NA" = "home-cage")) 
    colData$APA2 <- revalue(colData$APA2, c("Trained_Conflict" = "conflict-trained")) 
    colData$APA2 <- revalue(colData$APA2, c("Trained_NoConflict" = "standard-trained")) 
    colData$APA2 <- revalue(colData$APA2, c("Yoked_Conflict" = "conflict-yoked")) 
    colData$APA2 <- revalue(colData$APA2, c("Yoked_NoConflict" = "standard-yoked")) 
    head(colData)

    ##              RNAseqID   Mouse Punch      Group   Conflict
    ## 143A-CA3-1 143A-CA3-1 15-143A   CA3   conflict   Conflict
    ## 143A-DG-1   143A-DG-1 15-143A    DG   conflict   Conflict
    ## 143B-CA1-1 143B-CA1-1 15-143B   CA1    control   Conflict
    ## 143B-DG-1   143B-DG-1 15-143B    DG    control   Conflict
    ## 143C-CA1-1 143C-CA1-1 15-143C   CA1 consistent NoConflict
    ## 143D-CA1-3 143D-CA1-3 15-143D   CA1    control NoConflict
    ##                          APA3 Treatment     ID             APA2
    ## 143A-CA3-1   Trained_Conflict  conflict 15143A conflict-trained
    ## 143A-DG-1    Trained_Conflict  conflict 15143A conflict-trained
    ## 143B-CA1-1     Yoked_Conflict   shocked 15143B   conflict-yoked
    ## 143B-DG-1      Yoked_Conflict   shocked 15143B   conflict-yoked
    ## 143C-CA1-1 Trained_NoConflict   trained 15143C standard-trained
    ## 143D-CA1-3   Yoked_NoConflict     yoked 15143D   standard-yoked

    # reorder 
    colData <- colData[c(1:5,7:9)]

Now, we are ready to calculate differential gene expression using the
DESeq package. For simplicity, I will use the standard nameing of
"countData" and "colData" for the gene counts and gene information,
respectively.

    colData <- colData %>% arrange(RNAseqID) %>% droplevels() #set the coldata to be the countbygene df

    ## colData and countData must contain the exact same sample. I'll use the next three lines to make that happen
    savecols <- as.character(colData$RNAseqID) #select the sample name column that corresponds to row names
    savecols <- as.vector(savecols) # make it a vector
    countData <- countData %>% dplyr::select(one_of(savecols)) # select just the columns that match the samples in colData

    # colData must be factors
    cols = c(1:8)
    colData[,cols] %<>% lapply(function(x) as.factor(as.character(x)))

    # summary data
    colData %>% select(APA2,Punch)  %>%  summary()

    ##                APA2    Punch   
    ##  conflict-trained:14   CA1:17  
    ##  conflict-yoked  :12   CA3:15  
    ##  home-cage       : 6   DG :18  
    ##  standard-trained: 9           
    ##  standard-yoked  : 9

    dim(countData)

    ## [1] 22485    50

Write the two files
-------------------

    write.csv(colData, file = "../data/02a_colData.csv", row.names = F)
    write.csv(countData, file = "../data/02a_countData.csv", row.names = T)

Total Gene Counts Per Sample
----------------------------

this could say something about data before normalization

    ## stats
    counts <- countData
    dim( counts )

    ## [1] 22485    50

    colSums( counts ) / 1e06  # in millions of reads

    ## 143A-CA3-1  143A-DG-1 143B-CA1-1  143B-DG-1 143C-CA1-1 143D-CA1-3 
    ##   3.327867   5.279392   1.719498   2.085031   2.213452   1.091672 
    ##  143D-DG-3 144A-CA1-2 144A-CA3-2  144A-DG-2 144B-CA1-1 144B-CA3-1 
    ##   1.043885   2.980775   0.421165   3.210030   2.555909   1.027388 
    ## 144C-CA1-2 144C-CA3-2  144C-DG-2 144D-CA3-2  144D-DG-2 145A-CA1-2 
    ##   3.298825   1.238998   2.224182   2.323243   4.691568   4.680960 
    ## 145A-CA3-2  145A-DG-2 145B-CA1-1  145B-DG-1 146A-CA1-2 146A-CA3-2 
    ##   0.345619   1.435833   2.020114   1.509310   1.715282   2.756300 
    ##  146A-DG-2 146B-CA1-2 146B-CA3-2  146B-DG-2 146C-CA1-4  146C-DG-4 
    ##   1.201333   1.063417   2.144771   0.116106   1.360004   0.492145 
    ## 146D-CA1-3 146D-CA3-3  146D-DG-3  147-CA1-4  147-CA3-4   147-DG-4 
    ##   0.391369   2.994536   0.090417   0.159069   0.689232   0.139276 
    ## 147C-CA1-3 147C-CA3-3  147C-DG-3 147D-CA3-1  147D-DG-1  148-CA1-2 
    ##   3.072308   5.754581   4.350647   4.624995  11.700703   1.901256 
    ##  148-CA3-2   148-DG-2 148A-CA1-3 148A-CA3-3  148A-DG-3 148B-CA1-4 
    ##   2.343035   2.231849   5.260906   2.676397   4.019062   0.337174 
    ## 148B-CA3-4  148B-DG-4 
    ##   3.486840   0.798668

    table( rowSums( counts ) )[ 1:30 ] # Number of genes with low counts

    ## 
    ##    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
    ## 4169  341  229  217  164  148  120  117  106  117   86   80   81   71   65 
    ##   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29 
    ##   54   61   53   59   56   42   52   47   47   37   43   49   40   36   28

    rowsum <- as.data.frame(colSums( counts ) / 1e06 )
    names(rowsum)[1] <- "millioncounts"
    rowsum$sample <- row.names(rowsum)

    ggplot(rowsum, aes(x=millioncounts)) + 
      geom_histogram(binwidth = 1, colour = "black", fill = "darkgrey") +
      theme_classic() +
      scale_x_continuous(name = "Millions of Gene Counts per Sample") +
      scale_y_continuous(name = "Number of Samples")

![](../figures/02a_makedfs/totalRNAseqcounts-1.png)
