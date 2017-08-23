RNAseq gene expression analysis with DESeq2
-------------------------------------------

This workflow was modified from the DESeq2 tutorial found at: <https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf>

First I load a handful of packages for data wrangling, gene expression analysis, data visualization, and statistics.

``` r
library(dplyr) ## for filtering and selecting rows
library(plyr) ## for renmaing factors
library(reshape2) ## for melting dataframe
library(DESeq2) ## for gene expression analysis
library(edgeR)  ## for basic read counts status
library(magrittr) ## to use the weird pipe
library(genefilter)  ## for PCA fuction

## Functions
source("functions_RNAseq.R")
source("resvalsfunction.R")

## set output file for figures 
knitr::opts_chunk$set(fig.path = '../figures/02_RNAseq/')
```

Now, I create data frames from three csv files - count: Contains counts for all transcripts generated from the program Kallisto. This data can be reproducibed from the file kallisto.Rmd - geneids: Contains the ensemble ids and gene names for all the transcripts in the counts data frame. This file will be used to convert transcipt counts to gene counts. This file was also created via kallisto.Rmd file - Traits: This file contains all the information I collected for each sample that was sequenced. Not all columns will be needed, so some are removed later.

``` r
count <- read.csv("../data/02_count.csv", row.names=1, check.names=FALSE )
geneids <- read.csv("../data/02_geneids.csv")
Traits <- read.csv("../data/02_JA16444samples.csv", sep=",", header = TRUE, stringsAsFactors=FALSE, na.string = "NA")
```

Rather than analyze transcript level counts, I want to examine gene-level counts. In this next section, I do some data wrangling to tidy the data and summarize the counts for all transcripts belonging to the same gene.

``` r
countbygene <- full_join(geneids, count) # merge count and gene id dataframes
```

    ## Joining, by = "id"

``` r
countbygene <- countbygene[-c(1:6,8:12)] # remove unnecesary columns (aka, keep gene name and counts for samples)
countbygene <- melt(countbygene, id=c("gene")) # lenghten 
countbygene$variable <- gsub('.{4}$', '', countbygene$variable) # string split to remove last for characters aka "_S##"
countbygene$variable <- gsub("\\_", "-", countbygene$variable) #replace _ with - in name
countbygene  <- dcast(countbygene, gene ~ variable, value.var= "value", fun.aggregate=mean) # widen by sum
row.names(countbygene) <- countbygene$gene # make gene the row name 
countbygene[1] <- NULL # now remove the gene name from the df
countbygene <- round(countbygene) # round all value to nearest 1s place
```

In this next section, I tidy the trait data for each sample so that I can calculate differential gene expression for the traits of interest. I also remove some samples for reasons described within the code blocks.

``` r
rownames(Traits) <- Traits$RNAseqID    # set $genoAPAsessionInd as rownames
Traits <- Traits[c(1,3,5:6,10:11)]  #keeping informative volumns
Traits <- Traits %>% dplyr::filter(!grepl("100|101|147-|148-|147D-CA1-1|145B-CA3-1|146C-CA3-4", RNAseqID))  
# remove 100, 100, 147, and 148:  homecage animals
# Remove 147D_CA1_1 and 145B_CA3_1: bad samples with no reads.
# Remove 146C-CA3-4: outlier on all pc analyses
Traits$APAconflict <- as.factor(paste(Traits$APA, Traits$Conflict, sep="_")) # adding combinatorial traits columns
Traits$ID <- gsub("[[:punct:]]", "", Traits$Mouse) #make a column that thas id without the dash
row.names(Traits) <- Traits$RNAseqID # make gene the row name 
Traits$APA <- NULL ## delete old APA column
names(Traits)[names(Traits)=="APAconflict"] <- "APA3" #rename  APA3 to match color scheme

# rename factors & group all Control animals into 1 group
Traits$APA <- Traits$APA3 
Traits$APA <- revalue(Traits$APA, c("Trained_Conflict" = "Conflict")) 
Traits$APA <- revalue(Traits$APA, c("Trained_NoConflict" = "Consistent")) 
Traits$APA <- revalue(Traits$APA, c("Yoked_Conflict" = "Control")) 
Traits$APA <- revalue(Traits$APA, c("Yoked_NoConflict" = "Control")) 
head(Traits)
```

    ##              RNAseqID   Mouse   Conflict Punch Slice               APA3
    ## 143A-CA3-1 143A-CA3-1 15-143A   Conflict   CA3     1   Trained_Conflict
    ## 143A-DG-1   143A-DG-1 15-143A   Conflict    DG     1   Trained_Conflict
    ## 143B-CA1-1 143B-CA1-1 15-143B   Conflict   CA1     1     Yoked_Conflict
    ## 143B-DG-1   143B-DG-1 15-143B   Conflict    DG     1     Yoked_Conflict
    ## 143C-CA1-1 143C-CA1-1 15-143C NoConflict   CA1     1 Trained_NoConflict
    ## 143D-CA1-3 143D-CA1-3 15-143D NoConflict   CA1     3   Yoked_NoConflict
    ##                ID        APA
    ## 143A-CA3-1 15143A   Conflict
    ## 143A-DG-1  15143A   Conflict
    ## 143B-CA1-1 15143B    Control
    ## 143B-DG-1  15143B    Control
    ## 143C-CA1-1 15143C Consistent
    ## 143D-CA1-3 15143D    Control

``` r
# reorder 
Traits <- Traits[c(1:5,7,8,6)]
```

Now, we are ready to calculate differential gene expression using the DESeq package. For simplicity, I will use the standard nameing of "countData" and "colData" for the gene counts and gene information, respectively.

``` r
countData <- countbygene #set the countdata to be the countbygene df
colData <- Traits %>% arrange(RNAseqID) #set the coldata to be the countbygene df

## colData and countData must contain the exact same sample. I'll use the next three lines to make that happen
savecols <- as.character(colData$RNAseqID) #select the sample name column that corresponds to row names
savecols <- as.vector(savecols) # make it a vector
countData <- countData %>% dplyr::select(one_of(savecols)) # select just the columns that match the samples in colData


# colData must be factors
cols = c(1:8)
colData[,cols] %<>% lapply(function(x) as.factor(as.character(x)))

# summary data
colData %>% select(APA,Punch)  %>%  summary()
```

    ##          APA     Punch   
    ##  Conflict  :14   CA1:15  
    ##  Consistent: 9   CA3:13  
    ##  Control   :21   DG :16

``` r
colData %>% select(APA3,Punch)  %>%  summary()
```

    ##                  APA3    Punch   
    ##  Trained_Conflict  :14   CA1:15  
    ##  Trained_NoConflict: 9   CA3:13  
    ##  Yoked_Conflict    :12   DG :16  
    ##  Yoked_NoConflict  : 9

``` r
dim(countData)
```

    ## [1] 22485    44

``` r
write.csv(colData, file = "../data/02a_colData.csv", row.names = F)
write.csv(countData, file = "../data/02a_countData.csv", row.names = T)
```
