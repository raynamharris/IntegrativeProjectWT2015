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
knitr::opts_chunk$set(fig.path = '../figures/02_rnaseq/')
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
head(countbygene)
```

    ##               100-CA1-1 100-CA1-2 100-CA1-3 100-CA3-1 100-CA3-4 100-DG-2
    ## 0610007P14Rik        21        78        28        32        40       24
    ## 0610009B22Rik         7        43         6        12         6        8
    ## 0610009L18Rik         3        35         2        11         8        2
    ## 0610009O20Rik        44       225        66        87        88       58
    ## 0610010F05Rik        10        35        11         8        12        7
    ## 0610010K14Rik         2         7         2         1         3        1
    ##               100-DG-3 101-CA1-1 101-CA1-2 101-CA1-3 101-CA3-1 101-CA3-4
    ## 0610007P14Rik       95        58         0         0        38        16
    ## 0610009B22Rik       30        16         2         0        18         5
    ## 0610009L18Rik       40        12         0         0        10         2
    ## 0610009O20Rik      200       165         3         8       100        34
    ## 0610010F05Rik       28        13         1         1         6         6
    ## 0610010K14Rik        9         7         0         0         0         0
    ##               101-DG-3 101-DG-4 143A-CA3-1 143A-DG-1 143B-CA1-1 143B-DG-1
    ## 0610007P14Rik        2        3         42        56         30        24
    ## 0610009B22Rik        1        0         12        17         10         5
    ## 0610009L18Rik        2        8          4         9         10         8
    ## 0610009O20Rik        4       13         85       185         44        72
    ## 0610010F05Rik        0        1         18        19          7        15
    ## 0610010K14Rik        0        2          2         7          1         2
    ##               143C-CA1-1 143D-CA1-3 143D-DG-3 144A-CA1-2 144A-CA3-2
    ## 0610007P14Rik         19         14        22         40         10
    ## 0610009B22Rik         10          0         0         15          4
    ## 0610009L18Rik          2          0         2          9          5
    ## 0610009O20Rik         76         25        38         95         20
    ## 0610010F05Rik          7          5         8         15          2
    ## 0610010K14Rik          2          2         1          3          1
    ##               144A-DG-2 144B-CA1-1 144B-CA3-1 144C-CA1-2 144C-CA3-2
    ## 0610007P14Rik        40         36         17         32         14
    ## 0610009B22Rik         4          7          4          8         12
    ## 0610009L18Rik         0          1          4          2          4
    ## 0610009O20Rik        82         49         24         89         48
    ## 0610010F05Rik        17          9          3         12          7
    ## 0610010K14Rik         2          4          1          3          1
    ##               144C-DG-2 144D-CA3-2 144D-DG-2 145A-CA1-2 145A-CA3-2
    ## 0610007P14Rik        25         22        75         66          2
    ## 0610009B22Rik         8          6        12         18          3
    ## 0610009L18Rik         6          9        13         21          0
    ## 0610009O20Rik        85         46       151         96         12
    ## 0610010F05Rik        12         14        18         22          1
    ## 0610010K14Rik         2          2         6          3          0
    ##               145A-DG-2 145B-CA1-1 145B-CA3-1 145B-DG-1 146A-CA1-2
    ## 0610007P14Rik        20         26          0        10         14
    ## 0610009B22Rik         7          8          0         5          8
    ## 0610009L18Rik         1          3          0         0          8
    ## 0610009O20Rik        52        124          0        46         68
    ## 0610010F05Rik         6          4          0         5          9
    ## 0610010K14Rik         2          2          0         2          3
    ##               146A-CA3-2 146A-DG-2 146B-CA1-2 146B-CA3-2 146B-DG-2
    ## 0610007P14Rik         44        12          6         16         0
    ## 0610009B22Rik          4         3          3         22         0
    ## 0610009L18Rik          9         6          0          2         0
    ## 0610009O20Rik        102        44         16         46         2
    ## 0610010F05Rik         10         6          7         14         1
    ## 0610010K14Rik          1         2          1          1         1
    ##               146C-CA1-4 146C-CA3-4 146C-DG-4 146D-CA1-3 146D-CA3-3
    ## 0610007P14Rik         19          2        11          6         44
    ## 0610009B22Rik          4          4         3          2         12
    ## 0610009L18Rik          9          0         0          0          7
    ## 0610009O20Rik         31          3        10          1         83
    ## 0610010F05Rik          5          1         4          1         18
    ## 0610010K14Rik          1          0         0          1          1
    ##               146D-DG-3 147-CA1-4 147-CA3-4 147-DG-4 147C-CA1-3 147C-CA3-3
    ## 0610007P14Rik         3         1        46        0         34         82
    ## 0610009B22Rik         0         1        10        4          8         15
    ## 0610009L18Rik         0         0         0        0          2         11
    ## 0610009O20Rik         3         0        29        0         58        191
    ## 0610010F05Rik         0         2         0        2         19         41
    ## 0610010K14Rik         0         0         1        0          4          3
    ##               147C-DG-3 147D-CA1-1 147D-CA3-1 147D-DG-1 148-CA1-2
    ## 0610007P14Rik        41          0         40       152        24
    ## 0610009B22Rik        20          0         20        52         3
    ## 0610009L18Rik         3          0          9        67        13
    ## 0610009O20Rik       145          0         88       377        27
    ## 0610010F05Rik        23          0         28        56         9
    ## 0610010K14Rik         2          0          0        16         0
    ##               148-CA3-2 148-DG-2 148A-CA1-3 148A-CA3-3 148A-DG-3
    ## 0610007P14Rik        29       21         68         28        52
    ## 0610009B22Rik        23        8         30         12         8
    ## 0610009L18Rik         0        0         16          7        11
    ## 0610009O20Rik        91      164        162         65       226
    ## 0610010F05Rik        19        7         25         19        22
    ## 0610010K14Rik         2        1          5          4         2
    ##               148B-CA1-4 148B-CA3-4 148B-DG-4
    ## 0610007P14Rik          2         61         8
    ## 0610009B22Rik          0         22         1
    ## 0610009L18Rik          0         11         1
    ## 0610009O20Rik          0         70        18
    ## 0610010F05Rik          0         22         3
    ## 0610010K14Rik          0          4         1

``` r
summary(countbygene)
```

    ##    100-CA1-1          100-CA1-2         100-CA1-3       
    ##  Min.   :    0.00   Min.   :    0.0   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.0   1st Qu.:    0.00  
    ##  Median :    4.00   Median :   12.0   Median :    4.00  
    ##  Mean   :   50.55   Mean   :  147.3   Mean   :   49.58  
    ##  3rd Qu.:   29.00   3rd Qu.:   91.0   3rd Qu.:   30.00  
    ##  Max.   :31720.00   Max.   :95996.0   Max.   :24445.00  
    ##    100-CA3-1          100-CA3-4           100-DG-2       
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    4.00   Median :    5.00   Median :    3.00  
    ##  Mean   :   42.98   Mean   :   53.61   Mean   :   29.28  
    ##  3rd Qu.:   29.00   3rd Qu.:   33.00   3rd Qu.:   19.00  
    ##  Max.   :24878.00   Max.   :42838.00   Max.   :22711.00  
    ##     100-DG-3          101-CA1-1          101-CA1-2       
    ##  Min.   :     0.0   Min.   :     0.0   Min.   :   0.000  
    ##  1st Qu.:     0.0   1st Qu.:     0.0   1st Qu.:   0.000  
    ##  Median :    14.0   Median :     9.0   Median :   0.000  
    ##  Mean   :   135.9   Mean   :   118.7   Mean   :   3.204  
    ##  3rd Qu.:    92.0   3rd Qu.:    58.0   3rd Qu.:   2.000  
    ##  Max.   :100671.0   Max.   :183815.0   Max.   :3478.000  
    ##    101-CA1-3          101-CA3-1          101-CA3-4       
    ##  Min.   :   0.000   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:   0.000   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :   0.000   Median :    5.00   Median :    2.00  
    ##  Mean   :   6.868   Mean   :   60.53   Mean   :   28.46  
    ##  3rd Qu.:   4.000   3rd Qu.:   34.00   3rd Qu.:   15.00  
    ##  Max.   :6174.000   Max.   :86004.00   Max.   :37665.00  
    ##     101-DG-3           101-DG-4         143A-CA3-1         143A-DG-1      
    ##  Min.   :   0.000   Min.   :   0.00   Min.   :    0.00   Min.   :    0.0  
    ##  1st Qu.:   0.000   1st Qu.:   0.00   1st Qu.:    0.00   1st Qu.:    0.0  
    ##  Median :   0.000   Median :   0.00   Median :    5.00   Median :   10.0  
    ##  Mean   :   1.623   Mean   :  13.37   Mean   :   71.18   Mean   :  118.4  
    ##  3rd Qu.:   1.000   3rd Qu.:   8.00   3rd Qu.:   42.00   3rd Qu.:   75.0  
    ##  Max.   :2351.000   Max.   :9988.00   Max.   :23989.00   Max.   :76185.0  
    ##    143B-CA1-1        143B-DG-1          143C-CA1-1      
    ##  Min.   :    0.0   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.0   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    3.0   Median :    4.00   Median :    4.00  
    ##  Mean   :   38.9   Mean   :   45.32   Mean   :   48.31  
    ##  3rd Qu.:   23.0   3rd Qu.:   32.00   3rd Qu.:   30.00  
    ##  Max.   :30026.0   Max.   :21691.00   Max.   :21143.00  
    ##    143D-CA1-3         143D-DG-3          144A-CA1-2      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    1.00   Median :    2.00   Median :    6.00  
    ##  Mean   :   24.05   Mean   :   22.28   Mean   :   66.91  
    ##  3rd Qu.:   15.00   3rd Qu.:   16.00   3rd Qu.:   43.00  
    ##  Max.   :12512.00   Max.   :10111.00   Max.   :44270.00  
    ##    144A-CA3-2         144A-DG-2          144B-CA1-1      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    1.00   Median :    6.00   Median :    4.00  
    ##  Mean   :   10.14   Mean   :   70.08   Mean   :   56.71  
    ##  3rd Qu.:    5.00   3rd Qu.:   45.00   3rd Qu.:   35.00  
    ##  Max.   :12302.00   Max.   :33414.00   Max.   :35177.00  
    ##    144B-CA3-1         144C-CA1-2         144C-CA3-2     
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.0  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.0  
    ##  Median :    1.00   Median :    6.00   Median :    2.0  
    ##  Mean   :   22.54   Mean   :   71.77   Mean   :   28.8  
    ##  3rd Qu.:   13.00   3rd Qu.:   46.00   3rd Qu.:   16.0  
    ##  Max.   :12714.00   Max.   :35028.00   Max.   :26898.0  
    ##    144C-DG-2          144D-CA3-2         144D-DG-2         145A-CA1-2     
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.0   Min.   :    0.0  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.0   1st Qu.:    0.0  
    ##  Median :    4.00   Median :    4.00   Median :    9.0   Median :    8.0  
    ##  Mean   :   48.18   Mean   :   53.78   Mean   :  100.3   Mean   :  105.6  
    ##  3rd Qu.:   32.00   3rd Qu.:   29.00   3rd Qu.:   67.0   3rd Qu.:   62.0  
    ##  Max.   :15607.00   Max.   :46442.00   Max.   :39287.0   Max.   :73533.0  
    ##    145A-CA3-2         145A-DG-2          145B-CA1-1      
    ##  Min.   :   0.000   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:   0.000   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :   0.000   Median :    2.00   Median :    3.00  
    ##  Mean   :   8.004   Mean   :   30.73   Mean   :   45.99  
    ##  3rd Qu.:   4.000   3rd Qu.:   21.00   3rd Qu.:   26.00  
    ##  Max.   :7739.000   Max.   :11567.00   Max.   :39155.00  
    ##    145B-CA3-1         145B-DG-1          146A-CA1-2      
    ##  Min.   :0.000000   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:0.000000   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :0.000000   Median :    3.00   Median :    3.00  
    ##  Mean   :0.003514   Mean   :   32.06   Mean   :   39.06  
    ##  3rd Qu.:0.000000   3rd Qu.:   22.00   3rd Qu.:   24.00  
    ##  Max.   :3.000000   Max.   :15903.00   Max.   :28994.00  
    ##    146A-CA3-2         146A-DG-2          146B-CA1-2       146B-CA3-2      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :   0.0   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:   0.0   1st Qu.:    0.00  
    ##  Median :    5.00   Median :    2.00   Median :   2.0   Median :    4.00  
    ##  Mean   :   67.24   Mean   :   26.33   Mean   :  22.5   Mean   :   46.97  
    ##  3rd Qu.:   36.00   3rd Qu.:   18.00   3rd Qu.:  14.0   3rd Qu.:   29.00  
    ##  Max.   :89329.00   Max.   :13368.00   Max.   :7314.0   Max.   :23297.00  
    ##    146B-DG-2         146C-CA1-4         146C-CA3-4         146C-DG-4      
    ##  Min.   :   0.00   Min.   :    0.00   Min.   :   0.000   Min.   :   0.00  
    ##  1st Qu.:   0.00   1st Qu.:    0.00   1st Qu.:   0.000   1st Qu.:   0.00  
    ##  Median :   0.00   Median :    2.00   Median :   0.000   Median :   0.00  
    ##  Mean   :   2.47   Mean   :   29.48   Mean   :   5.964   Mean   :  10.56  
    ##  3rd Qu.:   1.00   3rd Qu.:   18.00   3rd Qu.:   3.000   3rd Qu.:   7.00  
    ##  Max.   :3901.00   Max.   :15918.00   Max.   :5751.000   Max.   :2745.00  
    ##    146D-CA1-3         146D-CA3-3         146D-DG-3       
    ##  Min.   :   0.000   Min.   :    0.00   Min.   :   0.000  
    ##  1st Qu.:   0.000   1st Qu.:    0.00   1st Qu.:   0.000  
    ##  Median :   0.000   Median :    5.00   Median :   0.000  
    ##  Mean   :   8.644   Mean   :   65.28   Mean   :   1.934  
    ##  3rd Qu.:   4.000   3rd Qu.:   40.00   3rd Qu.:   0.000  
    ##  Max.   :9683.000   Max.   :38276.00   Max.   :3002.000  
    ##    147-CA1-4         147-CA3-4           147-DG-4          147C-CA1-3   
    ##  Min.   :   0.00   Min.   :    0.00   Min.   :   0.000   Min.   :    0  
    ##  1st Qu.:   0.00   1st Qu.:    0.00   1st Qu.:   0.000   1st Qu.:    0  
    ##  Median :   0.00   Median :    0.00   Median :   0.000   Median :    6  
    ##  Mean   :   3.59   Mean   :   15.32   Mean   :   3.098   Mean   :   67  
    ##  3rd Qu.:   1.00   3rd Qu.:    7.00   3rd Qu.:   0.000   3rd Qu.:   43  
    ##  Max.   :5374.00   Max.   :14487.00   Max.   :3772.000   Max.   :37687  
    ##    147C-CA3-3         147C-DG-3          147D-CA1-1        147D-CA3-1     
    ##  Min.   :     0.0   Min.   :    0.00   Min.   :0.00000   Min.   :    0.0  
    ##  1st Qu.:     0.0   1st Qu.:    0.00   1st Qu.:0.00000   1st Qu.:    0.0  
    ##  Median :    13.0   Median :   10.00   Median :0.00000   Median :    8.0  
    ##  Mean   :   134.3   Mean   :   94.22   Mean   :0.00378   Mean   :  105.7  
    ##  3rd Qu.:    79.0   3rd Qu.:   66.00   3rd Qu.:0.00000   3rd Qu.:   59.0  
    ##  Max.   :150301.0   Max.   :46988.00   Max.   :6.00000   Max.   :95754.0  
    ##    147D-DG-1         148-CA1-2          148-CA3-2       
    ##  Min.   :    0.0   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.0   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :   24.0   Median :    2.00   Median :    4.00  
    ##  Mean   :  249.9   Mean   :   41.76   Mean   :   51.06  
    ##  3rd Qu.:  176.0   3rd Qu.:   24.00   3rd Qu.:   31.00  
    ##  Max.   :91299.0   Max.   :24841.00   Max.   :23437.00  
    ##     148-DG-2          148A-CA1-3        148A-CA3-3      
    ##  Min.   :    0.00   Min.   :    0.0   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.0   1st Qu.:    0.00  
    ##  Median :    3.00   Median :    9.0   Median :    4.00  
    ##  Mean   :   47.46   Mean   :  115.2   Mean   :   60.18  
    ##  3rd Qu.:   33.00   3rd Qu.:   72.0   3rd Qu.:   34.00  
    ##  Max.   :15414.00   Max.   :52783.0   Max.   :32891.00  
    ##    148A-DG-3          148B-CA1-4         148B-CA3-4      
    ##  Min.   :    0.00   Min.   :    0.00   Min.   :    0.00  
    ##  1st Qu.:    0.00   1st Qu.:    0.00   1st Qu.:    0.00  
    ##  Median :    9.00   Median :    0.00   Median :    6.00  
    ##  Mean   :   85.29   Mean   :    8.26   Mean   :   76.68  
    ##  3rd Qu.:   60.00   3rd Qu.:    3.00   3rd Qu.:   45.00  
    ##  Max.   :31971.00   Max.   :33665.00   Max.   :37680.00  
    ##    148B-DG-4       
    ##  Min.   :    0.00  
    ##  1st Qu.:    0.00  
    ##  Median :    2.00  
    ##  Mean   :   17.71  
    ##  3rd Qu.:   12.00  
    ##  Max.   :10089.00

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
names(Traits)[names(Traits)=="APAconflict"] <- "APA" #rename APAconflict APA (for simplicity)

# rename factors & group all yoked animals into 1 group
Traits$APA <- revalue(Traits$APA, c("Trained_Conflict" = "Conflict")) 
Traits$APA <- revalue(Traits$APA, c("Trained_NoConflict" = "Same")) 
Traits$APA <- revalue(Traits$APA, c("Yoked_Conflict" = "Yoked")) 
Traits$APA <- revalue(Traits$APA, c("Yoked_NoConflict" = "Yoked")) 
head(Traits)
```

    ##              RNAseqID   Mouse   Conflict Punch Slice      APA     ID
    ## 143A-CA3-1 143A-CA3-1 15-143A   Conflict   CA3     1 Conflict 15143A
    ## 143A-DG-1   143A-DG-1 15-143A   Conflict    DG     1 Conflict 15143A
    ## 143B-CA1-1 143B-CA1-1 15-143B   Conflict   CA1     1    Yoked 15143B
    ## 143B-DG-1   143B-DG-1 15-143B   Conflict    DG     1    Yoked 15143B
    ## 143C-CA1-1 143C-CA1-1 15-143C NoConflict   CA1     1     Same 15143C
    ## 143D-CA1-3 143D-CA1-3 15-143D NoConflict   CA1     3    Yoked 15143D

Now, we are ready to calculate differential gene expression using the DESeq package. For simplicity, I will use the standard nameing of "countData" and "colData" for the gene counts and gene information, respectively.

``` r
countData <- countbygene #set the countdata to be the countbygene df
colData <- Traits %>% arrange(RNAseqID) #set the coldata to be the countbygene df

## colData and countData must contain the exact same sample. I'll use the next three lines to make that happen
savecols <- as.character(colData$RNAseqID) #select the sample name column that corresponds to row names
savecols <- as.vector(savecols) # make it a vector
countData <- countData %>% dplyr::select(one_of(savecols)) # select just the columns that match the samples in colData


# colData must be factors
cols = c(1:4,7)
colData[,cols] %<>% lapply(function(x) as.factor(as.character(x)))
colData$Slice <- as.factor(colData$Slice)
str(colData)
```

    ## 'data.frame':    44 obs. of  7 variables:
    ##  $ RNAseqID: Factor w/ 44 levels "143A-CA3-1","143A-DG-1",..: 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ Mouse   : Factor w/ 18 levels "15-143A","15-143B",..: 1 1 2 2 3 4 4 5 5 5 ...
    ##  $ Conflict: Factor w/ 2 levels "Conflict","NoConflict": 1 1 1 1 2 2 2 1 1 1 ...
    ##  $ Punch   : Factor w/ 3 levels "CA1","CA3","DG": 2 3 1 3 1 1 3 1 2 3 ...
    ##  $ Slice   : Factor w/ 4 levels "1","2","3","4": 1 1 1 1 1 3 3 2 2 2 ...
    ##  $ APA     : Factor w/ 3 levels "Conflict","Same",..: 1 1 3 3 2 3 3 1 1 1 ...
    ##  $ ID      : Factor w/ 18 levels "15143A","15143B",..: 1 1 2 2 3 4 4 5 5 5 ...

Now, we are create differential gene expression object and remove genes with 0 counts. Before filtering, there are 22,485 genes in the object. After filtering genes with 0 counts, we will be left with 17,746 genes that are expressed in a least 1 sample. Then, we can caluate the size factors, estimate gene dispersion estimates, fit a model, test for outliers, and remove outliers.

``` r
## create DESeq object using the factors Punch and APA
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ Punch + APA + Punch*APA)
```

    ## converting counts to integer mode

``` r
## DESeq2 1.3.7 specify the factor levels
dds$Punch <- factor(dds$Punch, levels=c("DG","CA3", "CA1"))
dds$APA <- factor(dds$APA, levels=c("Yoked", "Same", "Conflict"))


dds # view the DESeq object - note numnber of genes
```

    ## class: DESeqDataSet 
    ## dim: 22485 44 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(22485): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(44): 143A-CA3-1 143A-DG-1 ... 148B-CA3-4 148B-DG-4
    ## colData names(7): RNAseqID Mouse ... APA ID

``` r
## DESeq2 1.3.6 Pre-filtering genes with 0 counts
dds <- dds[ rowSums(counts(dds)) > 1, ] 

dds # view the DESeq object - note numnber of genes
```

    ## class: DESeqDataSet 
    ## dim: 17674 44 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(17674): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(44): 143A-CA3-1 143A-DG-1 ... 148B-CA3-4 148B-DG-4
    ## colData names(7): RNAseqID Mouse ... APA ID

``` r
# dim: 17746 45 
# 17,746 genes and 45 samples

## DESeq2 1.4  Differential expression analysi
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 7 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
## for variance stablized gene expression and log transformed data
rld <- rlog(dds, blind=FALSE)
#vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
#vsd.fast <- vst(dds, blind=FALSE)
```

Now, we can calculate the number of differentiall expressed genes for each contrast. See the function in "resvalsfunction.R" for details.

``` r
source("resvalsfunction.R")
contrast1 <- resvals(contrastvector = c("Punch", "CA1", "DG"), mypadj = 0.05) #2497
```

    ## [1] 2497

``` r
contrast2 <- resvals(contrastvector = c("Punch", "CA1", "CA3"), mypadj = 0.05) #1803
```

    ## [1] 1803

``` r
contrast3 <- resvals(contrastvector = c("Punch", "CA3", "DG"), mypadj = 0.05) #3445
```

    ## [1] 3445

``` r
contrast4 <- resvals(contrastvector = c("APA", "Same", "Yoked"), mypadj = 0.05) #95
```

    ## [1] 95

``` r
contrast5 <- resvals(contrastvector = c("APA", "Conflict", "Yoked"), mypadj = 0.05) #42
```

    ## [1] 42

``` r
contrast6 <- resvals(contrastvector = c("APA", "Conflict", "Same"), mypadj = 0.05) # 0 
```

    ## [1] 0

``` r
#create a new DF with the gene counts
## note: contrast1 had 0 differentially expressed genes, so it is not included 
rldpadjs <- assay(rld)
rldpadjs <- cbind(rldpadjs, contrast1, contrast2, contrast3, contrast4, contrast5, contrast6)
rldpadjs <- as.data.frame(rldpadjs)
rldpadjs <- rldpadjs[ , grepl( "padj" , names( rldpadjs ) ) ]
head(rldpadjs)
```

    ##               padjPunchCA1DG padjPunchCA1CA3 padjPunchCA3DG
    ## 0610007P14Rik      0.9959885       1.0000000      0.9953218
    ## 0610009B22Rik      0.8920416       0.4918862      0.1034590
    ## 0610009L18Rik      0.6926063       0.6752603      0.9767398
    ## 0610009O20Rik      0.7325075       0.9860516      0.4348118
    ## 0610010F05Rik      0.7118161       0.1823768      0.3936081
    ## 0610010K14Rik      0.9959885       0.5247152      0.2851455
    ##               padjAPASameYoked padjAPAConflictYoked padjAPAConflictSame
    ## 0610007P14Rik                1                    1                   1
    ## 0610009B22Rik                1                    1                   1
    ## 0610009L18Rik                1                    1                   1
    ## 0610009O20Rik                1                    1                   1
    ## 0610010F05Rik                1                    1                   1
    ## 0610010K14Rik                1                    1                   1

``` r
volcano1 <- respadjfold(contrastvector = c("Punch", "CA1", "DG")) 
volcano2 <- respadjfold(contrastvector = c("Punch", "CA1", "CA3")) 
volcano3 <- respadjfold(contrastvector = c("Punch", "CA3", "DG")) 
volcano4 <- respadjfold(contrastvector = c("APA", "Same", "Yoked")) 
volcano5 <- respadjfold(contrastvector = c("APA", "Conflict", "Yoked")) 
volcano6 <- respadjfold(contrastvector = c("APA", "Conflict", "Same")) 

volcanos <- cbind(volcano1,volcano2, volcano3, volcano4, volcano5, volcano6)
volcanos <- as.data.frame(volcanos)
```

Now let's look at a heat map of the data

``` r
DEGes <- assay(rld)
DEGes <- cbind(DEGes, contrast1, contrast2, contrast3, contrast4, contrast5, contrast6)
DEGes <- as.data.frame(DEGes) # convert matrix to dataframe
DEGes$rownames <- rownames(DEGes)  # add the rownames to the dataframe
DEGes$padjmin <- with(DEGes, pmin(padjPunchCA1DG, padjPunchCA1CA3, padjPunchCA3DG, padjAPASameYoked, padjAPAConflictYoked, padjAPAConflictSame)) # create new col with min padj
DEGes <- DEGes %>% filter(padjmin < 0.05)
rownames(DEGes) <- DEGes$rownames
drop.cols <-colnames(DEGes[,grep("padj|pval|rownames", colnames(DEGes))])
DEGes <- DEGes %>% dplyr::select(-one_of(drop.cols))
DEGes <- as.matrix(DEGes)
DEGes <- DEGes - rowMeans(DEGes)
head(DEGes)
```

    ##                143A-CA3-1   143A-DG-1 143B-CA1-1  143B-DG-1 143C-CA1-1
    ## 1110002E22Rik -1.17062965  1.16379383  0.1336557  1.6190139 -0.6782200
    ## 1110008P14Rik  0.36927931  0.08527737  0.1932882  0.3345603 -0.2305912
    ## 1110012L19Rik -1.01406115  0.15190875 -0.4999942  1.4520809  0.1372077
    ## 1190002N15Rik -1.71094068  0.87439179  0.8886826 -0.8287686  1.3205345
    ## 1700001L19Rik -0.38880277  0.04666205  0.0797948 -0.8249946  0.1009883
    ## 1700001O22Rik  0.03095492 -0.34440929  1.4238938 -1.0345635  1.0126539
    ##                143D-CA1-3  143D-DG-3 144A-CA1-2 144A-CA3-2    144A-DG-2
    ## 1110002E22Rik -0.54696368  1.3063606 -0.5532396 -0.4945070  1.409349738
    ## 1110008P14Rik  0.08968605 -0.2042056  0.1854291  1.1511067 -0.966571550
    ## 1110012L19Rik -0.75569759 -0.1356544  0.8681495  0.8474443 -1.004864720
    ## 1190002N15Rik  0.82641426  0.3743473  0.9836874 -2.0396439  2.043854141
    ## 1700001L19Rik  1.57735013  0.1742612  0.1614121 -0.1651675  0.181071220
    ## 1700001O22Rik  1.18417461  0.3604488  1.3026481 -0.5132948  0.007678619
    ##                144B-CA1-1  144B-CA3-1  144C-CA1-2 144C-CA3-2  144C-DG-2
    ## 1110002E22Rik -0.90331591 -0.79872707 -0.60823782 -0.3420103  1.4161957
    ## 1110008P14Rik -0.91125940  0.32456437 -0.21269379  0.9871236 -0.4374099
    ## 1110012L19Rik  0.05513958  0.06754007  0.62029557 -0.1587617  0.2596652
    ## 1190002N15Rik  1.03636092 -1.22646311  0.33962537 -1.0881270  1.7848510
    ## 1700001L19Rik  0.23195726 -0.54942243  0.08975955 -0.8300424 -0.2683466
    ## 1700001O22Rik  0.80736819 -0.31618449  0.54993501 -0.8464314 -0.8591813
    ##                144D-CA3-2  144D-DG-2  145A-CA1-2   145A-CA3-2  145A-DG-2
    ## 1110002E22Rik -1.04515942  0.9459104 -0.68141493  0.003149581  0.1129408
    ## 1110008P14Rik  0.57648302 -0.0829272 -0.01396450  0.865795418 -0.1037691
    ## 1110012L19Rik -0.14371378  0.3498476  0.01343134  1.181093959 -0.2921593
    ## 1190002N15Rik -0.99788302 -0.3156500  0.98076845 -1.100280568  0.5961078
    ## 1700001L19Rik -0.15054695  0.1614783  0.68929978 -0.073316637 -0.6968575
    ## 1700001O22Rik -0.01174999  0.1163171  0.76595195 -0.037995397 -0.6998790
    ##               145B-CA1-1    145B-DG-1 146A-CA1-2 146A-CA3-2   146A-DG-2
    ## 1110002E22Rik -1.0087776  1.545806825 -0.9723121  0.4703477  1.24339240
    ## 1110008P14Rik -0.9133877  0.104118151 -0.1087641  0.8081349 -0.02833311
    ## 1110012L19Rik -0.8976013  0.250077842  0.6959071  0.2851445  0.47524730
    ## 1190002N15Rik  0.5991391 -0.008861396  0.7765526 -0.8766732  0.82531605
    ## 1700001L19Rik  0.5472203  0.063293967  0.3618480 -0.5915478  0.28772501
    ## 1700001O22Rik  0.8345009 -0.518189707  1.4729829  0.5655674 -0.24786451
    ##               146B-CA1-2  146B-CA3-2  146B-DG-2 146C-CA1-4  146C-DG-4
    ## 1110002E22Rik  0.8727357 -1.03834593  0.3640463  0.6892547 -0.5833415
    ## 1110008P14Rik -1.3376166  0.62794995 -0.1849090 -0.8267128 -0.9006632
    ## 1110012L19Rik -0.1129676 -0.91956311 -0.1624537 -0.8041846  0.2800476
    ## 1190002N15Rik  1.3822158 -0.90607154  1.1001809  1.2744499  1.2563017
    ## 1700001L19Rik  0.9251486  0.08585096 -0.2031513  0.6850137 -0.2845466
    ## 1700001O22Rik  0.3916913 -0.52158600 -0.1824397  1.1653560 -0.2562017
    ##               146D-CA1-3 146D-CA3-3  146D-DG-3 147C-CA1-3 147C-CA3-3
    ## 1110002E22Rik -0.4829639 -1.1312959 -0.1209408 -0.2732693 -0.6465445
    ## 1110008P14Rik  0.2578079 -0.4310780 -0.8579243 -0.5717336  0.2038659
    ## 1110012L19Rik -0.4644786 -0.3874941 -0.1277593 -0.1296321 -0.5829517
    ## 1190002N15Rik -0.1794671 -2.1265122 -1.0563482  0.9527176 -1.0217795
    ## 1700001L19Rik -0.5115538 -0.4761342 -0.1681240  0.3040984 -0.1849206
    ## 1700001O22Rik -0.5019169 -1.1061821 -0.1467142  0.4333845 -1.2396908
    ##                147C-DG-3  147D-CA3-1  147D-DG-1  148A-CA1-3  148A-CA3-3
    ## 1110002E22Rik  1.8955758 -1.24750682  1.2232659 -0.94352026 -0.23573927
    ## 1110008P14Rik -0.4257153  0.83081817 -0.1715680 -0.02561183  1.04656762
    ## 1110012L19Rik  0.3212809 -0.06717396  0.3567257  0.01871763 -0.23589012
    ## 1190002N15Rik  1.1439133 -1.86266575 -0.3923846  1.16426445 -2.18577660
    ## 1700001L19Rik  0.2139276 -0.30761652 -0.1121860  0.69108355  0.08359638
    ## 1700001O22Rik -0.3924831 -1.20579228 -0.8352275  0.74035623 -0.89759923
    ##                 148A-DG-3 148B-CA1-4  148B-CA3-4  148B-DG-4
    ## 1110002E22Rik  1.61172630 -0.4036172 -1.16990160  0.0539803
    ## 1110008P14Rik  0.09149396  0.1377390  0.49551158  0.1808092
    ## 1110012L19Rik -0.18345983 -0.3928404 -0.09674135  0.8831455
    ## 1190002N15Rik  0.27408930 -0.1918781 -1.31171728 -1.3708738
    ## 1700001L19Rik -0.44746094 -0.4353280 -0.31880166  0.2460275
    ## 1700001O22Rik -1.05861891 -0.4230747 -0.02269119  1.0540977

``` r
## the heatmap annotation file
df <- as.data.frame(colData(dds)[,c("Punch","APA")]) ## matrix to df
```

Now lets look at a principle component analysis of the data

``` r
# create the dataframe using my function pcadataframe
pcadata <- pcadataframe(rld, intgroup=c("Punch","APA"), returnData=TRUE)
percentVar <- round(100 * attr(pcadata, "percentVar"))
percentVar
```

    ## [1] 49 21  5  3  2  1  1  1  1

``` r
## statistics
aov1 <- aov(PC1 ~ Punch, data=pcadata)
summary(aov1) 
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)    
    ## Punch        2  16378    8189     254 <2e-16 ***
    ## Residuals   41   1322      32                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TukeyHSD(aov1, which = "Punch") 
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ Punch, data = pcadata)
    ## 
    ## $Punch
    ##               diff        lwr        upr     p adj
    ## CA3-DG  -40.657971 -45.813529 -35.502413 0.0000000
    ## CA1-DG  -39.611585 -44.573891 -34.649278 0.0000000
    ## CA1-CA3   1.046387  -4.185641   6.278415 0.8781712

``` r
aov2 <- aov(PC2 ~ Punch, data=pcadata)
summary(aov2) 
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)    
    ## Punch        2   7310    3655   948.1 <2e-16 ***
    ## Residuals   41    158       4                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TukeyHSD(aov2, which = "Punch") 
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ Punch, data = pcadata)
    ## 
    ## $Punch
    ##              diff       lwr       upr p adj
    ## CA3-DG  -16.60420 -18.38688 -14.82152     0
    ## CA1-DG   15.78048  14.06463  17.49634     0
    ## CA1-CA3  32.38468  30.57556  34.19380     0

``` r
aov3 <- aov(PC3 ~ APA, data=pcadata)
summary(aov3) 
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## APA          2  221.8  110.92   2.634 0.0839 .
    ## Residuals   41 1726.8   42.12                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TukeyHSD(aov3, which = "APA")
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC3 ~ APA, data = pcadata)
    ## 
    ## $APA
    ##                      diff        lwr       upr     p adj
    ## Same-Yoked     -4.5144437 -10.801630 1.7727431 0.2006707
    ## Conflict-Yoked -4.4831349  -9.927998 0.9617286 0.1244768
    ## Conflict-Same   0.0313088  -6.710948 6.7735655 0.9999297

``` r
aov4 <- aov(PC4 ~ APA, data=pcadata)
summary(aov4) 
```

    ##             Df Sum Sq Mean Sq F value  Pr(>F)    
    ## APA          2  382.4  191.18   11.35 0.00012 ***
    ## Residuals   41  690.8   16.85                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TukeyHSD(aov4, which = "APA") 
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC4 ~ APA, data = pcadata)
    ## 
    ## $APA
    ##                     diff       lwr       upr     p adj
    ## Same-Yoked      6.988720  3.012118 10.965321 0.0003232
    ## Conflict-Yoked  4.893361  1.449523  8.337199 0.0036237
    ## Conflict-Same  -2.095358 -6.359788  2.169072 0.4629850

``` r
aov5 <- aov(PC5 ~ APA, data=pcadata)
summary(aov5) 
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## APA          2    6.0   3.018   0.207  0.814
    ## Residuals   41  598.3  14.593

``` r
TukeyHSD(aov5, which = "APA") 
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC5 ~ APA, data = pcadata)
    ## 
    ## $APA
    ##                      diff       lwr      upr     p adj
    ## Same-Yoked      0.7175734 -2.983349 4.418496 0.8850340
    ## Conflict-Yoked -0.3255313 -3.530624 2.879562 0.9669609
    ## Conflict-Same  -1.0431047 -5.011902 2.925692 0.7995031

``` r
lm124 <- lm(PC1+PC2+PC4~APA*Punch, data=pcadata)
summary(lm124)
```

    ## 
    ## Call:
    ## lm(formula = PC1 + PC2 + PC4 ~ APA * Punch, data = pcadata)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -14.7950  -2.6364   0.1244   2.2727  19.8373 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            15.774      2.505   6.297 3.14e-07 ***
    ## APASame                20.708      4.797   4.317 0.000124 ***
    ## APAConflict            16.416      4.039   4.064 0.000259 ***
    ## PunchCA3              -51.557      3.826 -13.474 2.05e-15 ***
    ## PunchCA1              -16.670      3.667  -4.546 6.27e-05 ***
    ## APASame:PunchCA3      -10.584      7.515  -1.408 0.167821    
    ## APAConflict:PunchCA3   -9.754      5.892  -1.655 0.106780    
    ## APASame:PunchCA1      -15.659      6.537  -2.396 0.022076 *  
    ## APAConflict:PunchCA1  -14.350      6.003  -2.390 0.022342 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 7.085 on 35 degrees of freedom
    ## Multiple R-squared:  0.933,  Adjusted R-squared:  0.9177 
    ## F-statistic: 60.97 on 8 and 35 DF,  p-value: < 2.2e-16

``` r
anova(lm124) 
```

    ## Analysis of Variance Table
    ## 
    ## Response: PC1 + PC2 + PC4
    ##           Df  Sum Sq Mean Sq  F value    Pr(>F)    
    ## APA        2  1077.0   538.5  10.7273 0.0002325 ***
    ## Punch      2 22958.2 11479.1 228.6763 < 2.2e-16 ***
    ## APA:Punch  4   448.5   112.1   2.2339 0.0852730 .  
    ## Residuals 35  1756.9    50.2                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
write.csv(colData, file = "../data/02a_colData.csv", row.names = F)
write.csv(countData, file = "../data/02a_countData.csv", row.names = T)
write.csv(rldpadjs, file = "../data/02a_rldpadjs.csv", row.names = T)
write.csv(DEGes, file = "../data/02a_DEGes.csv", row.names = T)
write.csv(df, file = "../data/02a_df.csv", row.names = F)
write.csv(pcadata, file = "../data/02a_pcadata.csv", row.names = F)
write.table(percentVar, file = "../data/02a_percentVar.txt")
```

make volcano plos here.. still perfecting
-----------------------------------------

``` r
res <- results(dds, contrast =c("Punch", "CA1", "DG"), independentFiltering = F)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DG - CA1", xlim=c(-10,10)))
with(subset(res, log2FoldChange>0), points(log2FoldChange, -log10(pvalue), pch=20, col=c("#7570b3")))
with(subset(res, log2FoldChange<0), points(log2FoldChange, -log10(pvalue), pch=20, col=c("#d95f02")))
with(subset(res, padj>.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="grey"))
```

![](../figures/02_rnaseq/volcanoplots-1.png)

``` r
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("Punch"))
```

![](../figures/02_rnaseq/volcanoplots-2.png)

``` r
res <- results(dds, contrast =c("Punch", "CA1", "CA3"), independentFiltering = F)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="CA3-CA1", xlim=c(-10,10)))
with(subset(res, log2FoldChange>0), points(log2FoldChange, -log10(pvalue), pch=20, col=c("#7570b3")))
with(subset(res, log2FoldChange<0), points(log2FoldChange, -log10(pvalue), pch=20, col=c("#1b9e77")))
with(subset(res, padj>.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="grey"))
```

![](../figures/02_rnaseq/volcanoplots-3.png)

``` r
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("Punch"))
```

![](../figures/02_rnaseq/volcanoplots-4.png)

``` r
res <- results(dds, contrast =c("Punch", "CA3", "DG"), independentFiltering = F)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DG - CA3", xlim=c(-10,10)))
with(subset(res, log2FoldChange>0), points(log2FoldChange, -log10(pvalue), pch=20, col=c("#1b9e77")))
with(subset(res, log2FoldChange<0), points(log2FoldChange, -log10(pvalue), pch=20, col=c("#d95f02")))
with(subset(res, padj>.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="grey"))
```

![](../figures/02_rnaseq/volcanoplots-5.png)

``` r
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("Punch"))
```

![](../figures/02_rnaseq/volcanoplots-6.png)

``` r
res <- results(dds, contrast =c("APA", "Conflict", "Same"), independentFiltering = F)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Same - Conflict", xlim=c(-10,10)))
with(subset(res, log2FoldChange>0), points(log2FoldChange, -log10(pvalue), pch=20, col=c("#ca0020")))
with(subset(res, log2FoldChange<0), points(log2FoldChange, -log10(pvalue), pch=20, col=c("#f4a582")))
with(subset(res, padj>.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="grey"))
```

![](../figures/02_rnaseq/volcanoplots-7.png)

``` r
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("APA"))
```

![](../figures/02_rnaseq/volcanoplots-8.png)

``` r
res <- results(dds, contrast =c("APA", "Conflict", "Yoked"), independentFiltering = F)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Yoked - Conflict", xlim=c(-10,10)))
with(subset(res, log2FoldChange>0), points(log2FoldChange, -log10(pvalue), pch=20, col=c("#ca0020")))
with(subset(res, log2FoldChange<0), points(log2FoldChange, -log10(pvalue), pch=20, col=c("#404040")))
with(subset(res, padj>.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="grey"))
```

![](../figures/02_rnaseq/volcanoplots-9.png)

``` r
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("APA"))
```

![](../figures/02_rnaseq/volcanoplots-10.png)

``` r
res <- results(dds, contrast =c("APA", "Same", "Yoked"), independentFiltering = F)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Yoked - Same", xlim=c(-10,10)))
with(subset(res, log2FoldChange>0), points(log2FoldChange, -log10(pvalue), pch=20, col=c("#f4a582")))
with(subset(res, log2FoldChange<0), points(log2FoldChange, -log10(pvalue), pch=20, col=c("#404040")))
with(subset(res, padj>.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="grey"))
```

![](../figures/02_rnaseq/volcanoplots-11.png)

``` r
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("APA"), transform = T)
```

![](../figures/02_rnaseq/volcanoplots-12.png)
