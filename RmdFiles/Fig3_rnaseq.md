RNAseq gene expression analysis with DESeq2
-------------------------------------------

This workflow was modified from the DESeq2 tutorial found at: <https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf>

First I load a handful of packages for data wrangling, gene expression analysis, data visualization, and statistics.

``` r
library(dplyr) ## for filtering and selecting rows
library(plyr) ## for renmaing factors
library(reshape2) ## for melting dataframe
library(magrittr) ## to use the weird pipe
library(gplots) ##for making awesome plots
library(cowplot) ## for some easy to use themes
library(ggplot2) ## for awesome plots!
library(DESeq2) ## for gene expression analysis
library(pheatmap) ## for awesome heatmaps
library(edgeR)  ## for basic read counts status
library(VennDiagram) ## for Venn Diagrams
```

Then, I set all path so that all the figures are saved in the specified subdirectory.

``` r
## set output file for figures 
knitr::opts_chunk$set(fig.path = '../figures/Fig3/')
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
Traits <- Traits %>% dplyr::filter(!grepl("100|101|147-|148-|147D-CA1-1|145B-CA3-1", RNAseqID))  # remove 100, 100, 147, and 148 animals because we aren't interested in these homecage animals that were not trained in the active place avoidance experiement. Remove mice 147D_CA1_1 and 145B_CA3_1 because these were bad samples with no reads
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

    ## 'data.frame':    45 obs. of  7 variables:
    ##  $ RNAseqID: Factor w/ 45 levels "143A-CA3-1","143A-DG-1",..: 1 2 3 4 5 6 7 8 9 10 ...
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
    ## dim: 22485 45 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(22485): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(45): 143A-CA3-1 143A-DG-1 ... 148B-CA3-4 148B-DG-4
    ## colData names(7): RNAseqID Mouse ... APA ID

``` r
## DESeq2 1.3.6 Pre-filtering genes with 0 counts
dds <- dds[ rowSums(counts(dds)) > 1, ] 

dds # view the DESeq object - note numnber of genes
```

    ## class: DESeqDataSet 
    ## dim: 17746 45 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(17746): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(45): 143A-CA3-1 143A-DG-1 ... 148B-CA3-4 148B-DG-4
    ## colData names(7): RNAseqID Mouse ... APA ID

``` r
## DESeq2 1.4  Differential expression analysi
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 18 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
## for variance stablized gene expression and log transformed data
rld <- rlog(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd.fast <- vst(dds, blind=FALSE)
```

Now, we can calculate the number of differentiall expressed genes for each contrast. See the function in "resvalsfunction.R" for details.

``` r
source("resvalsfunction.R")
contrast1 <- resvals(contrastvector = c("Punch", "DG", "CA1"), mypadj = 0.1)
```

    ## [1] 2999

``` r
contrast2 <- resvals(contrastvector = c("Punch", "CA3", "CA1"), mypadj = 0.1)
```

    ## [1] 2204

``` r
contrast3 <- resvals(contrastvector = c("Punch", "DG", "CA3"), mypadj = 0.1)
```

    ## [1] 4110

``` r
contrast4 <- resvals(contrastvector = c("APA", "Yoked", "Same"), mypadj = 0.1)
```

    ## [1] 107

``` r
contrast5 <- resvals(contrastvector = c("APA", "Yoked", "Conflict"), mypadj = 0.1)
```

    ## [1] 59

``` r
contrast6 <- resvals(contrastvector = c("APA", "Same", "Conflict"), mypadj = 0.1)
```

    ## [1] 0

``` r
#create a new DF with the gene counts
## note: contrast1 had 0 differentially expressed genes, so it is not included 
rldpadjs <- assay(rld)
rldpadjs <- cbind(rldpadjs, contrast1, contrast2, contrast3, contrast4, contrast5)
rldpadjs <- as.data.frame(rldpadjs)
rldpadjs <- rldpadjs[ , grepl( "padj" , names( rldpadjs ) ) ]
head(rldpadjs)
```

    ##               padjPunchDGCA1 padjPunchCA3CA1 padjPunchDGCA3
    ## 0610007P14Rik      0.9988301       1.0000000      1.0000000
    ## 0610009B22Rik      0.9058675       0.5192173      0.1194886
    ## 0610009L18Rik      0.6977672       0.6926354      0.9805026
    ## 0610009O20Rik      0.7262958       0.9875592      0.4348464
    ## 0610010F05Rik      0.7090649       0.1991185      0.4357406
    ## 0610010K14Rik      0.9988301       0.5327519      0.2821504
    ##               padjAPAYokedSame padjAPAYokedConflict
    ## 0610007P14Rik                1                    1
    ## 0610009B22Rik                1                    1
    ## 0610009L18Rik                1                    1
    ## 0610009O20Rik                1                    1
    ## 0610010F05Rik                1                    1
    ## 0610010K14Rik                1                    1

Here, the goal is the analyze the distribution of pvalues to see if they are randomly distributed or if that is a tendency towards and increase or decrease of low pvalues. There, I'm showing the pval and adjusted pvale (padj) for all for two-way comparision.

``` r
head(rldpadjs)
```

    ##               padjPunchDGCA1 padjPunchCA3CA1 padjPunchDGCA3
    ## 0610007P14Rik      0.9988301       1.0000000      1.0000000
    ## 0610009B22Rik      0.9058675       0.5192173      0.1194886
    ## 0610009L18Rik      0.6977672       0.6926354      0.9805026
    ## 0610009O20Rik      0.7262958       0.9875592      0.4348464
    ## 0610010F05Rik      0.7090649       0.1991185      0.4357406
    ## 0610010K14Rik      0.9988301       0.5327519      0.2821504
    ##               padjAPAYokedSame padjAPAYokedConflict
    ## 0610007P14Rik                1                    1
    ## 0610009B22Rik                1                    1
    ## 0610009L18Rik                1                    1
    ## 0610009O20Rik                1                    1
    ## 0610010F05Rik                1                    1
    ## 0610010K14Rik                1                    1

``` r
rldpadjslong <- rldpadjs
rldpadjslong$gene <- row.names(rldpadjslong) 
rldpadjslong <- melt(rldpadjslong, id=c("gene"))
head(rldpadjslong)
```

    ##            gene       variable     value
    ## 1 0610007P14Rik padjPunchDGCA1 0.9988301
    ## 2 0610009B22Rik padjPunchDGCA1 0.9058675
    ## 3 0610009L18Rik padjPunchDGCA1 0.6977672
    ## 4 0610009O20Rik padjPunchDGCA1 0.7262958
    ## 5 0610010F05Rik padjPunchDGCA1 0.7090649
    ## 6 0610010K14Rik padjPunchDGCA1 0.9988301

``` r
qplot(value, data=rldpadjslong, geom="histogram") + 
  facet_grid( ~ variable) +
  scale_y_log10()
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 85 rows containing non-finite values (stat_bin).

    ## Warning: Stacking not well defined when ymin != 0

![](../figures/Fig3/pvaluedistribution-1.png)

``` r
ggplot(rldpadjslong, aes(x=value)) + geom_density()
```

    ## Warning: Removed 85 rows containing non-finite values (stat_density).

![](../figures/Fig3/pvaluedistribution-2.png)

``` r
ggplot(rldpadjslong, aes(x=value, fill=variable)) +
  geom_density(alpha=.3) 
```

    ## Warning: Removed 85 rows containing non-finite values (stat_density).

![](../figures/Fig3/pvaluedistribution-3.png)

``` r
ggplot(rldpadjslong, aes(x=value, fill=variable)) +
  geom_density(alpha=.3) +
  scale_y_log10()
```

    ## Warning: Removed 85 rows containing non-finite values (stat_density).

![](../figures/Fig3/pvaluedistribution-4.png)

Now, we count the number of differnetially expressed genes (according to padj) and plot some venn diagrams.

``` r
venn1 <- row.names(rldpadjs[rldpadjs[1] <0.1 & !is.na(rldpadjs[1]),])
venn2 <- row.names(rldpadjs[rldpadjs[2] <0.1 & !is.na(rldpadjs[2]),])
venn3 <- row.names(rldpadjs[rldpadjs[3] <0.1 & !is.na(rldpadjs[3]),])
venn4 <- row.names(rldpadjs[rldpadjs[4] <0.1 & !is.na(rldpadjs[4]),])
venn5 <- row.names(rldpadjs[rldpadjs[5] <0.1 & !is.na(rldpadjs[5]),])

## check order for correctness
candidates <- list("DG vs CA1" = venn1, "Yoked vs Same" = venn4,"DG vs CA3" = venn3,  "CA3 vs CA1" = venn2)

prettyvenn <- venn.diagram(
  scaled=T,
  x = candidates, filename=NULL, 
  col = "black",
  fill = c( "white", "white", "white", "white"),
  alpha = 0.5,
  cex = 1, fontfamily = "sans", #fontface = "bold",
  cat.default.pos = "text",
  #cat.dist = c(0.08, 0.08, 0.08), cat.pos = 1,
  cat.cex = 1, cat.fontfamily = "sans")
grid.draw(prettyvenn)
```

![](../figures/Fig3/venndiagram-1.png)

``` r
plot.new()

candidates <- list("DG vs CA1" = venn1, "DG vs CA3" = venn3,  "CA3 vs CA1" = venn2)

prettyvenn <- venn.diagram(
  scaled=T,
  x = candidates, filename=NULL, 
  col = "black",
  fill = c( "white", "white", "white"),
  alpha = 0.5,
  cex = 1, fontfamily = "sans", #fontface = "bold",
  cat.default.pos = "text",
  #cat.dist = c(0.08, 0.08, 0.08), cat.pos = 1,
  cat.cex = 1, cat.fontfamily = "sans")
grid.draw(prettyvenn)
```

![](../figures/Fig3/venndiagram-2.png)

Now let's look at a heat map of the data

``` r
## set coloring for heat map anntations
df <- as.data.frame(colData(dds)[,c("Punch","APA")]) ## matrix to df
ann_colors = list(
  APA =  c(Yoked = (values=c("#8073ac")), Same = (values=c("#e08214")), Conflict = (values=c("#7f3b08"))),
  Punch =  c(DG = (values=c("#006837")),  CA3 = (values=c("#41ab5d")), 
             CA1 = (values=c("#d9f0a3"))))

DEGes <- assay(rld)
DEGes <- cbind(DEGes, contrast1, contrast2, contrast3, contrast4, contrast5)
DEGes <- as.data.frame(DEGes) # convert matrix to dataframe
DEGes$rownames <- rownames(DEGes)  # add the rownames to the dataframe
DEGes$padjmin <- with(DEGes, pmin(padjPunchDGCA1, padjPunchCA3CA1, padjPunchDGCA3, padjAPAYokedSame, padjAPAYokedConflict)) # create new col with min padj
DEGes <- DEGes %>% filter(padjmin < 0.01)
rownames(DEGes) <- DEGes$rownames
drop.cols <-colnames(DEGes[,grep("padj|pval|rownames", colnames(DEGes))])
DEGes <- DEGes %>% dplyr::select(-one_of(drop.cols))
DEGes <- as.matrix(DEGes)
DEGes <- DEGes - rowMeans(DEGes)
head(DEGes)
```

    ##                143A-CA3-1  143A-DG-1  143B-CA1-1   143B-DG-1  143C-CA1-1
    ## 1110002E22Rik -1.16230774  1.1705090  0.14432205  1.63110398 -0.66763649
    ## 1190002N15Rik -1.72506907  0.8755782  0.88793327 -0.82549312  1.32373075
    ## 1700017B05Rik  1.29813450 -0.2383787 -0.84876206 -1.16200667 -0.65883370
    ## 1700025G04Rik -0.43111905  0.6047283 -0.33863691  0.02097928 -0.09274292
    ## 1810010H24Rik -0.36806542  0.7967341  0.06788807  0.74720221  0.12590565
    ## 1810041L15Rik -0.04800454  1.2368213 -1.69968331  0.78541070 -0.98978547
    ##                143D-CA1-3  143D-DG-3 144A-CA1-2  144A-CA3-2   144A-DG-2
    ## 1110002E22Rik -0.54330515  1.3214462 -0.5395372 -0.50186871  1.42087835
    ## 1190002N15Rik  0.82572960  0.3824383  0.9864699 -2.08732342  2.05439746
    ## 1700017B05Rik  0.02743781 -1.2656349 -0.5637580  0.94786057  0.00303291
    ## 1700025G04Rik -0.80203072  0.3368837 -0.7794041  0.40143012  0.74260373
    ## 1810010H24Rik -0.45798840  0.8737816  0.1875457 -0.29357058 -0.17110780
    ## 1810041L15Rik -1.37081153  0.9329462 -1.1065663  0.05083861  1.43559517
    ##               144B-CA1-1   144B-CA3-1  144C-CA1-2  144C-CA3-2  144C-DG-2
    ## 1110002E22Rik -0.8942817 -0.802574577 -0.59401924 -0.33616232 1.42474894
    ## 1190002N15Rik  1.0372495 -1.240480093  0.34090582 -1.10153393 1.79036094
    ## 1700017B05Rik -0.4518229 -0.307333915 -0.47661869  1.19914121 0.05271724
    ## 1700025G04Rik -0.9835528 -0.307508308 -0.06585606  0.02148736 0.59522838
    ## 1810010H24Rik -0.5654846 -0.000859784  0.07789520  0.24857348 0.75948850
    ## 1810041L15Rik -1.1054536 -0.322546654 -0.62252195  0.56542278 1.26503568
    ##                144D-CA3-2  144D-DG-2 145A-CA1-2  145A-CA3-2  145A-DG-2
    ## 1110002E22Rik -1.04344508  0.9574167 -0.6657968  0.01053546  0.1287196
    ## 1190002N15Rik -1.01325284 -0.3124629  0.9819782 -1.11819161  0.6042550
    ## 1700017B05Rik  0.85948272 -0.2555703 -0.9363856 -0.08544741  0.2783899
    ## 1700025G04Rik  0.12503652  0.7399715 -0.7598524  0.96942707  0.2687021
    ## 1810010H24Rik -0.03471651  0.2408481 -0.1203884 -0.25561606 -0.1391646
    ## 1810041L15Rik -0.02059493  1.2791806 -1.1638068  0.58498224  0.9643158
    ##               145B-CA1-1    145B-DG-1  146A-CA1-2 146A-CA3-2  146A-DG-2
    ## 1110002E22Rik -1.0057109  1.556313069 -0.97101273  0.4790266  1.2586023
    ## 1190002N15Rik  0.6006954 -0.005450833  0.77622333 -0.8847532  0.8351507
    ## 1700017B05Rik -0.3135722 -0.603158148 -0.11444672  0.5592856 -0.5950860
    ## 1700025G04Rik -0.3337989  0.498527522 -0.30895547  0.1440145  0.3725948
    ## 1810010H24Rik -0.5327576  0.328876963  0.06547299 -0.5579117  0.8581657
    ## 1810041L15Rik -1.1155530  1.438908695 -1.53463962 -0.5755530  1.3098625
    ##                146B-CA1-2  146B-CA3-2   146B-DG-2  146C-CA1-4 146C-CA3-4
    ## 1110002E22Rik  0.88401686 -1.03561615  0.38401572  0.70096129 -0.3658572
    ## 1190002N15Rik  1.38487883 -0.91583497  1.12556775  1.27761462  0.1246173
    ## 1700017B05Rik  0.38539037  1.04175448 -0.32052087  0.17758647  0.3427710
    ## 1700025G04Rik -0.35484194  0.07325742  0.19667259 -0.84617479 -0.4568169
    ## 1810010H24Rik -0.01696541 -0.24961470 -0.01369138 -0.09645644 -0.1931296
    ## 1810041L15Rik -0.98434885  0.15102121 -0.58976954 -0.89642435 -0.2225356
    ##                146C-DG-4  146D-CA1-3  146D-CA3-3   146D-DG-3 147C-CA1-3
    ## 1110002E22Rik -0.5885531 -0.48934895 -1.12381112 -0.11663828 -0.2577157
    ## 1190002N15Rik  1.2634222 -0.18075430 -2.14192834 -1.08213907  0.9585774
    ## 1700017B05Rik -0.2862878 -1.46709905  1.72581565 -0.70822971 -0.7567584
    ## 1700025G04Rik  0.6169602 -0.85871033 -0.08087894  0.80634323 -0.4728886
    ## 1810010H24Rik  1.1169687 -0.28670975 -0.57631108  0.01662962 -0.2573655
    ## 1810041L15Rik  1.1958631  0.08952726  0.44251290  0.38603596 -1.2010405
    ##               147C-CA3-3  147C-DG-3 147D-CA3-1  147D-DG-1 148A-CA1-3
    ## 1110002E22Rik -0.6316126  1.9010310 -1.2350377  1.2330545 -0.9265685
    ## 1190002N15Rik -1.0288397  1.1482525 -1.8759059 -0.3891130  1.1655303
    ## 1700017B05Rik  0.6984457 -0.2595952  0.5180043  0.1401251 -0.4115819
    ## 1700025G04Rik  0.3067797  0.7224128  0.5653680  0.3044112 -0.2727652
    ## 1810010H24Rik -0.4623700  0.4051276 -0.4309390  0.6399662 -0.4059048
    ## 1810041L15Rik  0.2208505  1.2363855  0.3938438  1.2096327 -0.9118822
    ##               148A-CA3-3  148A-DG-3 148B-CA1-4  148B-CA3-4   148B-DG-4
    ## 1110002E22Rik -0.2256539  1.6226887 -0.4095695 -1.16175529  0.06600623
    ## 1190002N15Rik -2.2072486  0.2817692 -0.1913985 -1.32352344 -1.38262954
    ## 1700017B05Rik  0.7634976 -0.4472843  1.0062427  1.10535223  0.40370498
    ## 1700025G04Rik -0.2047959  0.6032214 -1.1253419 -0.08870245 -0.07166697
    ## 1810010H24Rik -0.0977259  0.2627586 -0.2325531 -0.58918127 -0.41327959
    ## 1810041L15Rik -0.4364960  0.8592621 -1.0321960 -0.24462694  0.16058534

``` r
pheatmap(DEGes, show_colnames=F, show_rownames = F,
         #annotation_col=df, annotation_colors = ann_colors,
         fontsize = 12, fontsize_row = 10, 
         border_color = "grey60"
)
```

![](../figures/Fig3/heatmap-1.png)

Now lets look at a principle component analysis of the data

``` r
plotPCA(rld, intgroup=c("APA", "Punch"), returnData=TRUE)
```

    ##                   PC1         PC2          group      APA Punch       name
    ## 143A-CA3-1 -17.792212 -15.2604021 Conflict : CA3 Conflict   CA3 143A-CA3-1
    ## 143A-DG-1   30.272394  -1.0154986  Conflict : DG Conflict    DG  143A-DG-1
    ## 143B-CA1-1 -14.570296  16.0530334    Yoked : CA1    Yoked   CA1 143B-CA1-1
    ## 143B-DG-1   26.771138  -1.7319963     Yoked : DG    Yoked    DG  143B-DG-1
    ## 143C-CA1-1 -14.901443  18.0914921     Same : CA1     Same   CA1 143C-CA1-1
    ## 143D-CA1-3 -15.163203  17.5024305    Yoked : CA1    Yoked   CA1 143D-CA1-3
    ## 143D-DG-3   26.985562  -1.1966857     Yoked : DG    Yoked    DG  143D-DG-3
    ## 144A-CA1-2 -10.002146  11.7128574 Conflict : CA1 Conflict   CA1 144A-CA1-2
    ## 144A-CA3-2 -13.386024 -12.9115544 Conflict : CA3 Conflict   CA3 144A-CA3-2
    ## 144A-DG-2   31.076199  -0.2770855  Conflict : DG Conflict    DG  144A-DG-2
    ## 144B-CA1-1 -16.237056  17.7079186    Yoked : CA1    Yoked   CA1 144B-CA1-1
    ## 144B-CA3-1 -10.207559 -12.2683464    Yoked : CA3    Yoked   CA3 144B-CA3-1
    ## 144C-CA1-2 -10.382000  13.6851480     Same : CA1     Same   CA1 144C-CA1-2
    ## 144C-CA3-2 -15.353031 -13.6441111     Same : CA3     Same   CA3 144C-CA3-2
    ## 144C-DG-2   32.994388   0.6484256      Same : DG     Same    DG  144C-DG-2
    ## 144D-CA3-2 -16.476016 -16.0102716    Yoked : CA3    Yoked   CA3 144D-CA3-2
    ## 144D-DG-2   30.847370  -1.0912677     Yoked : DG    Yoked    DG  144D-DG-2
    ## 145A-CA1-2 -15.574394  18.0431051 Conflict : CA1 Conflict   CA1 145A-CA1-2
    ## 145A-CA3-2  -7.963753  -9.8483263 Conflict : CA3 Conflict   CA3 145A-CA3-2
    ## 145A-DG-2   28.237011   0.3650217  Conflict : DG Conflict    DG  145A-DG-2
    ## 145B-CA1-1 -14.750868  17.0780543    Yoked : CA1    Yoked   CA1 145B-CA1-1
    ## 145B-DG-1   28.933376  -1.1425224     Yoked : DG    Yoked    DG  145B-DG-1
    ## 146A-CA1-2 -14.667132  16.1055163 Conflict : CA1 Conflict   CA1 146A-CA1-2
    ## 146A-CA3-2 -14.345862 -15.6450410 Conflict : CA3 Conflict   CA3 146A-CA3-2
    ## 146A-DG-2   28.274114  -0.4573622  Conflict : DG Conflict    DG  146A-DG-2
    ## 146B-CA1-2 -11.670416  15.8397533    Yoked : CA1    Yoked   CA1 146B-CA1-2
    ## 146B-CA3-2 -17.980826 -14.9389457    Yoked : CA3    Yoked   CA3 146B-CA3-2
    ## 146B-DG-2    7.737648  -1.4794107     Yoked : DG    Yoked    DG  146B-DG-2
    ## 146C-CA1-4 -12.287822  17.9126305     Same : CA1     Same   CA1 146C-CA1-4
    ## 146C-CA3-4 -15.458569 -35.8629446     Same : CA3     Same   CA3 146C-CA3-4
    ## 146C-DG-4   26.503359   0.3968615      Same : DG     Same    DG  146C-DG-4
    ## 146D-CA1-3 -10.305140  11.5331103    Yoked : CA1    Yoked   CA1 146D-CA1-3
    ## 146D-CA3-3 -15.013848 -15.1077783    Yoked : CA3    Yoked   CA3 146D-CA3-3
    ## 146D-DG-3    9.177987  -2.7725621     Yoked : DG    Yoked    DG  146D-DG-3
    ## 147C-CA1-3 -14.136545  16.8557902     Same : CA1     Same   CA1 147C-CA1-3
    ## 147C-CA3-3 -15.181895 -14.8710062     Same : CA3     Same   CA3 147C-CA3-3
    ## 147C-DG-3   32.130716   0.1262738      Same : DG     Same    DG  147C-DG-3
    ## 147D-CA3-1 -16.985316 -16.0287160    Yoked : CA3    Yoked   CA3 147D-CA3-1
    ## 147D-DG-1   33.087908  -1.2145155     Yoked : DG    Yoked    DG  147D-DG-1
    ## 148A-CA1-3 -11.943096  15.8531531 Conflict : CA1 Conflict   CA1 148A-CA1-3
    ## 148A-CA3-3 -18.249825 -16.1718546 Conflict : CA3 Conflict   CA3 148A-CA3-3
    ## 148A-DG-3   28.737964  -0.9842093  Conflict : DG Conflict    DG  148A-DG-3
    ## 148B-CA1-4 -11.769117  13.1580821    Yoked : CA1    Yoked   CA1 148B-CA1-4
    ## 148B-CA3-4 -17.594947 -15.3859638    Yoked : CA3    Yoked   CA3 148B-CA3-4
    ## 148B-DG-4    8.583223  -1.3502797     Yoked : DG    Yoked    DG  148B-DG-4

``` r
pcadata <- plotPCA(rld, intgroup=c("APA", "Punch"), returnData=TRUE)
percentVar <- round(100 * attr(pcadata, "percentVar"))


ggplot(pcadata, aes(PC1, PC2, color=Punch, shape=APA)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +  
  #stat_ellipse(level = 0.95, (aes(color=Punch)),size=1.5)   + 
  scale_color_manual(values=c("#006837", "#41ab5d", "#d9f0a3")) + 
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18))
```

![](../figures/Fig3/pca-1.png)

Now for some basic stats about the read and gene counts

``` r
## stats

counts <- countData
dim( counts )
```

    ## [1] 22485    45

``` r
colSums( counts ) / 1e06  # in millions of reads
```

    ## 143A-CA3-1  143A-DG-1 143B-CA1-1  143B-DG-1 143C-CA1-1 143D-CA1-3 
    ##   1.600535   2.662497   0.874614   1.019113   1.086358   0.540738 
    ##  143D-DG-3 144A-CA1-2 144A-CA3-2  144A-DG-2 144B-CA1-1 144B-CA3-1 
    ##   0.500935   1.504457   0.228021   1.575766   1.275137   0.506698 
    ## 144C-CA1-2 144C-CA3-2  144C-DG-2 144D-CA3-2  144D-DG-2 145A-CA1-2 
    ##   1.613785   0.647568   1.083336   1.209306   2.254320   2.375356 
    ## 145A-CA3-2  145A-DG-2 145B-CA1-1  145B-DG-1 146A-CA1-2 146A-CA3-2 
    ##   0.179967   0.690882   1.034066   0.720798   0.878270   1.511881 
    ##  146A-DG-2 146B-CA1-2 146B-CA3-2  146B-DG-2 146C-CA1-4 146C-CA3-4 
    ##   0.591933   0.506014   1.056001   0.055549   0.662938   0.134106 
    ##  146C-DG-4 146D-CA1-3 146D-CA3-3  146D-DG-3 147C-CA1-3 147C-CA3-3 
    ##   0.237419   0.194359   1.467877   0.043490   1.506436   3.020727 
    ##  147C-DG-3 147D-CA3-1  147D-DG-1 148A-CA1-3 148A-CA3-3  148A-DG-3 
    ##   2.118624   2.377445   5.618550   2.590815   1.353197   1.917857 
    ## 148B-CA1-4 148B-CA3-4  148B-DG-4 
    ##   0.185637   1.724144   0.398258

``` r
table( rowSums( counts ) )[ 1:30 ] # Number of genes with low counts
```

    ## 
    ##    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
    ## 4326  413  339  250  217  199  156  146  128  105   98   83   88   96   81 
    ##   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29 
    ##   60   73   77   72   66   59   70   59   54   66   55   53   49   40   39

``` r
rowsum <- as.data.frame(colSums( counts ) / 1e06 )
names(rowsum)[1] <- "millioncounts"
rowsum$sample <- row.names(rowsum)

ggplot(rowsum, aes(x=millioncounts)) + 
  geom_histogram(binwidth = 1, colour = "black", fill = "darkgrey") +
  theme_classic() +
  scale_x_continuous(name = "Millions of Gene Counts per Sample",
                     breaks = seq(0, 8, 1),
                     limits=c(0, 8)) +
  scale_y_continuous(name = "Number of Samples")
```

![](../figures/Fig3/stats-1.png)
