    library("tidyverse") 
    library("forcats")  

    #knitr::opts_chunk$set(cache = T)

Disclaimer
----------

Hello! If you are viewing this page then hopefully you want to access
some really large data files containing raw RNA transcript counts and
estimates of transcripts per million.

To obtain the data analyzed in this markdown file, I [ran
kallisto](https://github.com/raynamharris/IntegrativeProjectWT2015/blob/master/UNIXworkflow/04_kallisto.md)
on the Stampede Cluster at the Texas Advanced Computing Facility. The
data are exported as abunance files in a subdirectory for every sample.

These files and this analysis will take up considerable space and time.
If you want to run the analysis, first download
`GSE100225_IntegrativeWT2015` from this GitHub repository:
<a href="https://github.com/raynamharris/MouseHippocampusRNAseqData" class="uri">https://github.com/raynamharris/MouseHippocampusRNAseqData</a>
and save them in this repo in \`../data/.

Importing data from many Kallisto files into a single data frame
----------------------------------------------------------------

The kallisto output gives you read counts for sample in an abundance
file for every single sample. This portion of the code goes through and
finds each samples’ abundance.tsv file, extracts the data, and combines
it all into a dataframe. The `counts` file is unnormalized, but the
`tpm` is the data after being normalized by transcripts per million.
This script was developed with assistance from Anna Batthenhouse and
Dennis Whylie.

Rather than examine unique transcripts, my analyses will focus on
gene-level exprrssion. I use some string splitting to take the very long
transcript identifying and create a `geneids` file that has all the
database identifiers for each transcript. Then, I’ll save the dount
data.

    ## this will create lists of all the samples
    kallistoDirs = dir("../data/GSE100225_IntegrativeWT2015/")
    kallistoDirs = kallistoDirs[!grepl("\\.(R|py|pl|sh|xlsx?|txt|tsv|csv|org|md|obo|png|jpg|pdf)$",
    kallistoDirs, ignore.case=TRUE)]

    path <- c("../data/GSE100225_IntegrativeWT2015/")

    kallistoFiles = paste0(path, kallistoDirs, "/abundance.tsv")
    names(kallistoFiles) = kallistoDirs
    if(file.exists(kallistoFiles))
      kallistoData = lapply(
      kallistoFiles,
      read.table,
      sep = "\t",
      row.names = 1,
      header = TRUE
    )

    ## this for loop uses the reduce function to make two data frame with counts or tpm from all the samples. note, only counts are used

    ids = Reduce(f=union, x=lapply(kallistoData, rownames))
    if (all(sapply(kallistoData, function(x) {all(rownames(x)==ids)}))) {
      count = data.frame(
        id = ids,
        sapply(kallistoData, function(x) {x$est_counts}),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
      tpm = data.frame(
        id = ids,
        sapply(kallistoData, function(x) {x$tpm}),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    }

Gene ids
--------

This takes one column of information from the transcriptome and splits
it up to make a dataframe with tons of useful gene information

    head(count)

    ##                                                                                                                                                    id
    ## 1      ENSMUST00000070533.4|ENSMUSG00000051951.5|OTTMUSG00000026353.2|OTTMUST00000065166.1|Xkr4-001|Xkr4|3634|UTR5:1-150|CDS:151-2094|UTR3:2095-3634|
    ## 2                        ENSMUST00000208660.1|ENSMUSG00000025900.11|OTTMUSG00000049985.3|OTTMUST00000145515.1|Rp1-003|Rp1|4170|UTR5:1-54|CDS:55-4170|
    ## 3       ENSMUST00000027032.5|ENSMUSG00000025900.11|OTTMUSG00000049985.3|OTTMUST00000127195.2|Rp1-001|Rp1|6869|UTR5:1-127|CDS:128-6415|UTR3:6416-6869|
    ## 4                                                                            ENSMUST00000194992.6|ENSMUSG00000109048.1|-|-|Rp1-201|Rp1|858|CDS:1-858|
    ## 5 ENSMUST00000027035.9|ENSMUSG00000025902.13|OTTMUSG00000050014.7|OTTMUST00000127245.2|Sox17-001|Sox17|3127|UTR5:1-1082|CDS:1083-2342|UTR3:2343-3127|
    ## 6   ENSMUST00000195555.1|ENSMUSG00000025902.13|OTTMUSG00000050014.7|OTTMUST00000127249.1|Sox17-005|Sox17|1977|UTR5:1-635|CDS:636-1511|UTR3:1512-1977|
    ##   142C_CA1 142C_DG 143A-CA3-1 143A-DG-1 143B-CA1-1 143B-DG-1  143C_CA1
    ## 1   95.000     200         98        81         29        35 473.00000
    ## 2    0.000       0          0         0          0         0   0.00000
    ## 3    0.000       0          0         0          0         0   0.00000
    ## 4    0.000       0          0         0          0         0   0.00000
    ## 5    0.000       0          0         0          0         0   0.00000
    ## 6   27.363       0          0         0          0         0   4.89368
    ##   143C_DG 143C-CA1-1 143D-CA1-3 143D-DG-3 144A-CA1-2 144A-CA3-2 144A-DG-2
    ## 1      60         99         39        12         73          1        88
    ## 2       0          0          0         0          0          0         0
    ## 3       0          0          0         0          0          0         0
    ## 4       0          0          0         0          0          0         0
    ## 5       0          0          0         0          0          0         0
    ## 6       2          0          0         0          9          0        15
    ##   144B-CA1-1 144B-CA3-1 144C-CA1-2 144C-CA3-2 144C-DG-2 144D-CA3-2
    ## 1         47         23        173         17        60         60
    ## 2          0          0          0          0         0          0
    ## 3          0          0          0          0         0          0
    ## 4          0          0          0          0         0          0
    ## 5          0          0          0          2         0          0
    ## 6          0          0          1          0         0          4
    ##   144D-DG-2 145A-CA1-2 145A-CA3-2 145A-DG-2 145B-CA1-1 145B-DG-1
    ## 1  115.0000        133          3        46         60        27
    ## 2    0.0000          0          0         0          0         0
    ## 3    3.0000          0          0         0          0         0
    ## 4    0.0000          0          0         0          0         0
    ## 5    5.1248         44          0         0          0         0
    ## 6    0.0000          0          0         0          0         0
    ##   146A-CA1-2 146A-CA3-2 146A-DG-2 146B-CA1-2 146B-CA3-2 146B-DG-2
    ## 1         20          8        10         48         48         1
    ## 2          0          0         0          0          0         0
    ## 3          0          0         0          0          0         0
    ## 4          0          0         0          0          0         0
    ## 5          0          0         0          0          0         0
    ## 6          0         26         0          9          0         0
    ##   146C-CA1-4 146C-CA3-4 146C-DG-4 146D-CA1-3 146D-CA3-3 146D-DG-3
    ## 1   93.00000          1         2         15         53         0
    ## 2    0.00000          0         2          0          0         0
    ## 3    0.00000          0         0          0          0         0
    ## 4    0.00000          0         0          0          0         0
    ## 5    0.00000          0         0          0          0         0
    ## 6    4.37439          0         0          0          1         0
    ##   147-CA1-4 147-CA3-4 147-DG-4 147C-CA1-3 147C-CA3-3 147C-DG-3 147D-CA3-1
    ## 1         1        11        0  106.00000   237.0000 104.00000        160
    ## 2         0         0        0    0.00000     0.0000   0.00000          0
    ## 3         0         0        0    0.00000     0.0000   0.00000          0
    ## 4         0         0        0    0.00000     0.0000   0.00000          0
    ## 5         0         0        0    5.09227    10.8636   0.00000          0
    ## 6         0         0        0    1.90773     0.0000   7.10077          0
    ##   147D-DG-1 148-CA1-2 148-CA3-2 148-DG-2 148A-CA1-3 148A-CA3-3 148A-DG-3
    ## 1 227.00000        54        49       34        253         39       108
    ## 2   0.00000         0         0        0          0          7         0
    ## 3   0.00000         0         0        0          0          0         0
    ## 4   0.00000         0         0        0          0          0         0
    ## 5   1.80308         0         0        0          0          0         0
    ## 6   0.00000         0         0        0          0          1         0
    ##   148B-CA1-4 148B-CA3-4 148B-DG-4
    ## 1          5         73        17
    ## 2          0          0         0
    ## 3          0          0         0
    ## 4          0          0         0
    ## 5          3          1         1
    ## 6          0          0         0

    ## make a dataframe with the parts of the gene id as columns
    geneids <- count[c(1)] 
    geneids$ENSMUST <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 1)
    geneids$ENSMUSG <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 2)
    geneids$OTTMUSG <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 3)
    geneids$OTTMUST <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 4)
    geneids$transcript <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 5)
    geneids$gene <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 6)
    geneids$length <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 7)
    geneids$structure1 <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 8)
    geneids$structure2 <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 9)
    geneids$structure3 <- sapply(strsplit(as.character(geneids$id),'\\|'), "[", 10)
    geneids$transcript_lenght <- as.factor(paste(geneids$transcript, geneids$length, sep="_"))
    names(geneids)

    ##  [1] "id"                "ENSMUST"           "ENSMUSG"          
    ##  [4] "OTTMUSG"           "OTTMUST"           "transcript"       
    ##  [7] "gene"              "length"            "structure1"       
    ## [10] "structure2"        "structure3"        "transcript_lenght"

    # save slim version for joining immediately
    geneidgene <- geneids %>% select(id, gene)
    geneidtranscript <- geneids %>% select(id, transcript_lenght)

TPM
---

My downstream analyses don’t use TPM, but some people find that useful,
so I made one and saved it as a df for you.

    # tpm to tpmbytranscript - note: not used for subsequent analyses
    tpmbytranscript <-  full_join(geneidtranscript, tpm) %>% select(-id)   # merge tpm and genids

    ## Joining, by = "id"

    row.names(tpmbytranscript) <- tpmbytranscript$transcript_lenght ## make gene the row name
    tpmbytranscript[1] <- NULL ## make gene the row name
    tpmbytranscript <- round(tpmbytranscript) #round all value to nearest 1s place
    head(tpmbytranscript,3)

    ##               142C_CA1 142C_DG 143A-CA3-1 143A-DG-1 143B-CA1-1 143B-DG-1
    ## Xkr4-001_3634       13      38         19         9         10        10
    ## Rp1-003_4170         0       0          0         0          0         0
    ## Rp1-001_6869         0       0          0         0          0         0
    ##               143C_CA1 143C_DG 143C-CA1-1 143D-CA1-3 143D-DG-3 144A-CA1-2
    ## Xkr4-001_3634       48      11         28         21         7         14
    ## Rp1-003_4170         0       0          0          0         0          0
    ## Rp1-001_6869         0       0          0          0         0          0
    ##               144A-CA3-2 144A-DG-2 144B-CA1-1 144B-CA3-1 144C-CA1-2
    ## Xkr4-001_3634          1        19         11         13         34
    ## Rp1-003_4170           0         0          0          0          0
    ## Rp1-001_6869           0         0          0          0          0
    ##               144C-CA3-2 144C-DG-2 144D-CA3-2 144D-DG-2 145A-CA1-2
    ## Xkr4-001_3634          7        18         13        16         16
    ## Rp1-003_4170           0         0          0         0          0
    ## Rp1-001_6869           0         0          0         0          0
    ##               145A-CA3-2 145A-DG-2 145B-CA1-1 145B-DG-1 146A-CA1-2
    ## Xkr4-001_3634          4        22         17        12          6
    ## Rp1-003_4170           0         0          0         0          0
    ## Rp1-001_6869           0         0          0         0          0
    ##               146A-CA3-2 146A-DG-2 146B-CA1-2 146B-CA3-2 146B-DG-2
    ## Xkr4-001_3634          1         5         32         14         5
    ## Rp1-003_4170           0         0          0          0         0
    ## Rp1-001_6869           0         0          0          0         0
    ##               146C-CA1-4 146C-CA3-4 146C-DG-4 146D-CA1-3 146D-CA3-3
    ## Xkr4-001_3634         46          2         3         21         11
    ## Rp1-003_4170           0          0         2          0          0
    ## Rp1-001_6869           0          0         0          0          0
    ##               146D-DG-3 147-CA1-4 147-CA3-4 147-DG-4 147C-CA1-3 147C-CA3-3
    ## Xkr4-001_3634         0         3         8        0         22         22
    ## Rp1-003_4170          0         0         0        0          0          0
    ## Rp1-001_6869          0         0         0        0          0          0
    ##               147C-DG-3 147D-CA3-1 147D-DG-1 148-CA1-2 148-CA3-2 148-DG-2
    ## Xkr4-001_3634        16         19        13        18        13       10
    ## Rp1-003_4170          0          0         0         0         0        0
    ## Rp1-001_6869          0          0         0         0         0        0
    ##               148A-CA1-3 148A-CA3-3 148A-DG-3 148B-CA1-4 148B-CA3-4
    ## Xkr4-001_3634         30          8        19          9         13
    ## Rp1-003_4170           0          1         0          0          0
    ## Rp1-001_6869           0          0         0          0          0
    ##               148B-DG-4
    ## Xkr4-001_3634        12
    ## Rp1-003_4170          0
    ## Rp1-001_6869          0

Counts per gene
---------------

This sums the counts for all transcripts then rounds the measure to give
a df that will be used fof all down stream analysis.

    # count to countbygene
    countbygene <- full_join(geneidgene, count) %>% select(-id)  %>% 
      pivot_longer(-gene, names_to = "samples", values_to = "counts") %>%  
      pivot_wider(
        names_from = samples, 
        values_from = counts,
        values_fn = list(counts = sum))  %>% 
      arrange(gene)

    ## Joining, by = "id"

    countbygene <- as.data.frame(countbygene)
    row.names(countbygene) <- countbygene$gene ## make gene the row name
    countbygene[1] <- NULL ## make gene the row name
    countbygene <- round(countbygene) #round all value to nearest 1s place
    head(countbygene,3)

    ##               142C_CA1 142C_DG 143A-CA3-1 143A-DG-1 143B-CA1-1 143B-DG-1
    ## 0610007P14Rik      120      88         85       112         60        48
    ## 0610009B22Rik       62      27         24        34         21        10
    ## 0610009L18Rik        0       0          4         9         10         8
    ##               143C_CA1 143C_DG 143C-CA1-1 143D-CA1-3 143D-DG-3 144A-CA1-2
    ## 0610007P14Rik       83      32         38         28        43         80
    ## 0610009B22Rik       44      19         19          0         1         30
    ## 0610009L18Rik        6       0          2          0         2          9
    ##               144A-CA3-2 144A-DG-2 144B-CA1-1 144B-CA3-1 144C-CA1-2
    ## 0610007P14Rik         21        80         72         34         63
    ## 0610009B22Rik          9         9         14          9         15
    ## 0610009L18Rik          5         0          1          4          2
    ##               144C-CA3-2 144C-DG-2 144D-CA3-2 144D-DG-2 145A-CA1-2
    ## 0610007P14Rik         28        49         43       150        133
    ## 0610009B22Rik         24        16         13        23         36
    ## 0610009L18Rik          4         6          9        13         21
    ##               145A-CA3-2 145A-DG-2 145B-CA1-1 145B-DG-1 146A-CA1-2
    ## 0610007P14Rik          5        41         51        21         29
    ## 0610009B22Rik          6        14         15        10         15
    ## 0610009L18Rik          0         1          3         0          8
    ##               146A-CA3-2 146A-DG-2 146B-CA1-2 146B-CA3-2 146B-DG-2
    ## 0610007P14Rik         87        23         13         33         1
    ## 0610009B22Rik          9         6          6         43         1
    ## 0610009L18Rik          9         6          0          2         0
    ##               146C-CA1-4 146C-CA3-4 146C-DG-4 146D-CA1-3 146D-CA3-3
    ## 0610007P14Rik         38          4        22         11         87
    ## 0610009B22Rik          7          7         6          3         23
    ## 0610009L18Rik          9          0         0          0          7
    ##               146D-DG-3 147-CA1-4 147-CA3-4 147-DG-4 147C-CA1-3 147C-CA3-3
    ## 0610007P14Rik         6         2        92        0         69        164
    ## 0610009B22Rik         0         2        21        9         17         30
    ## 0610009L18Rik         0         0         0        0          2         11
    ##               147C-DG-3 147D-CA3-1 147D-DG-1 148-CA1-2 148-CA3-2 148-DG-2
    ## 0610007P14Rik        82         79       305        49        58       42
    ## 0610009B22Rik        39         41       105         6        46       15
    ## 0610009L18Rik         3          9        67        13         0        0
    ##               148A-CA1-3 148A-CA3-3 148A-DG-3 148B-CA1-4 148B-CA3-4
    ## 0610007P14Rik        135         55       104          5        122
    ## 0610009B22Rik         59         24        15          0         45
    ## 0610009L18Rik         16          7        11          0         11
    ##               148B-DG-4
    ## 0610007P14Rik        16
    ## 0610009B22Rik         2
    ## 0610009L18Rik         1

Prep for DESeq
--------------

DESeq is the tool I’ll be using downstream, so I use “colData” to
describe the meta data that will be used in downstream analyses

    # clean col Data
    colData <- read_csv("../data/IntegrativeWT2015ColData.csv") %>%
      mutate(ID = gsub("[[:punct:]]", "", Mouse)) %>%
      filter(APA_Conflict != "NA_NA") %>%
      mutate(subfield = Region) %>%
      mutate(treatment = fct_recode(APA_Conflict,
                                    "standard.yoked" = "Yoked_NoConflict",
                                    "standard.trained" = "Trained_NoConflict",
                                    "conflict.yoked" = "Yoked_Conflict",
                                    "conflict.trained" = "Trained_Conflict")) %>%
      mutate(training = fct_collapse(treatment,
                                          trained = c("standard.trained", "conflict.trained"),
                                          yoked = c("standard.yoked", "conflict.yoked"))) %>%
      select(RNAseqID,ID,subfield, treatment, training) %>%
      arrange(RNAseqID) %>%
      droplevels() 

    ## Parsed with column specification:
    ## cols(
    ##   RNAseqID = col_character(),
    ##   Mouse = col_character(),
    ##   year = col_double(),
    ##   Genotype = col_character(),
    ##   Region = col_character(),
    ##   jobnumber = col_character(),
    ##   Group = col_character(),
    ##   APA = col_character(),
    ##   Conflict = col_character(),
    ##   APA_Conflict = col_character(),
    ##   Treatment = col_character()
    ## )

    head(colData,3)

    ## # A tibble: 3 x 5
    ##   RNAseqID   ID     subfield treatment        training
    ##   <chr>      <chr>  <chr>    <fct>            <fct>   
    ## 1 143A-CA3-1 15143A CA3      conflict.trained trained 
    ## 2 143A-DG-1  15143A DG       conflict.trained trained 
    ## 3 143B-CA1-1 15143B CA1      conflict.yoked   yoked

Remove count data for samples that we aren’t analyzing.

    ## colData and countData must contain the exact same samples. 
    savecols <- as.character(colData$RNAseqID) #select the rowsname 
    savecols <- as.vector(savecols) # make it a vector
    countData <- countbygene %>% dplyr::select(one_of(savecols)) # select just the columns 
    head(countData,3)  

    ##               143A-CA3-1 143A-DG-1 143B-CA1-1 143B-DG-1 143C-CA1-1
    ## 0610007P14Rik         85       112         60        48         38
    ## 0610009B22Rik         24        34         21        10         19
    ## 0610009L18Rik          4         9         10         8          2
    ##               143D-CA1-3 143D-DG-3 144A-CA1-2 144A-CA3-2 144A-DG-2
    ## 0610007P14Rik         28        43         80         21        80
    ## 0610009B22Rik          0         1         30          9         9
    ## 0610009L18Rik          0         2          9          5         0
    ##               144B-CA1-1 144B-CA3-1 144C-CA1-2 144C-CA3-2 144C-DG-2
    ## 0610007P14Rik         72         34         63         28        49
    ## 0610009B22Rik         14          9         15         24        16
    ## 0610009L18Rik          1          4          2          4         6
    ##               144D-CA3-2 144D-DG-2 145A-CA1-2 145A-CA3-2 145A-DG-2
    ## 0610007P14Rik         43       150        133          5        41
    ## 0610009B22Rik         13        23         36          6        14
    ## 0610009L18Rik          9        13         21          0         1
    ##               145B-CA1-1 145B-DG-1 146A-CA1-2 146A-CA3-2 146A-DG-2
    ## 0610007P14Rik         51        21         29         87        23
    ## 0610009B22Rik         15        10         15          9         6
    ## 0610009L18Rik          3         0          8          9         6
    ##               146B-CA1-2 146B-CA3-2 146B-DG-2 146C-CA1-4 146C-DG-4
    ## 0610007P14Rik         13         33         1         38        22
    ## 0610009B22Rik          6         43         1          7         6
    ## 0610009L18Rik          0          2         0          9         0
    ##               146D-CA1-3 146D-CA3-3 146D-DG-3 147C-CA1-3 147C-CA3-3
    ## 0610007P14Rik         11         87         6         69        164
    ## 0610009B22Rik          3         23         0         17         30
    ## 0610009L18Rik          0          7         0          2         11
    ##               147C-DG-3 147D-CA3-1 147D-DG-1 148A-CA1-3 148A-CA3-3
    ## 0610007P14Rik        82         79       305        135         55
    ## 0610009B22Rik        39         41       105         59         24
    ## 0610009L18Rik         3          9        67         16          7
    ##               148A-DG-3 148B-CA1-4 148B-CA3-4 148B-DG-4
    ## 0610007P14Rik       104          5        122        16
    ## 0610009B22Rik        15          0         45         2
    ## 0610009L18Rik        11          0         11         1

Save files
----------

    write.csv(geneids, "../data/00_geneids.csv", row.names=F)
    write.csv(tpmbytranscript, "../data/00_tpmbytranscript.csv", row.names=T)
    write.csv(colData, file = "../data/00_colData.csv", row.names = F)
    write.csv(countData, file = "../data/00_countData.csv", row.names = T)  

Session Info
------------

    sessionInfo()

    ## R version 3.6.0 (2019-04-26)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Mojave 10.14.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] forcats_0.4.0   stringr_1.4.0   dplyr_0.8.3     purrr_0.3.3    
    ## [5] readr_1.3.1     tidyr_1.0.0     tibble_2.1.3    ggplot2_3.2.1  
    ## [9] tidyverse_1.3.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_0.2.5 xfun_0.9         haven_2.2.0      lattice_0.20-38 
    ##  [5] colorspace_1.4-1 vctrs_0.2.0      generics_0.0.2   htmltools_0.3.6 
    ##  [9] yaml_2.2.0       utf8_1.1.4       rlang_0.4.1      pillar_1.4.2    
    ## [13] glue_1.3.1       withr_2.1.2      DBI_1.0.0        dbplyr_1.4.2    
    ## [17] modelr_0.1.5     readxl_1.3.1     lifecycle_0.1.0  munsell_0.5.0   
    ## [21] gtable_0.3.0     cellranger_1.1.0 rvest_0.3.5      evaluate_0.14   
    ## [25] knitr_1.24       fansi_0.4.0      broom_0.5.2      Rcpp_1.0.2      
    ## [29] scales_1.0.0     backports_1.1.4  jsonlite_1.6     fs_1.3.1        
    ## [33] hms_0.5.2        digest_0.6.20    stringi_1.4.3    grid_3.6.0      
    ## [37] cli_1.1.0        tools_3.6.0      magrittr_1.5     lazyeval_0.2.2  
    ## [41] crayon_1.3.4     pkgconfig_2.0.2  zeallot_0.1.0    xml2_1.2.2      
    ## [45] reprex_0.3.0     lubridate_1.7.4  assertthat_0.2.1 rmarkdown_1.15  
    ## [49] httr_1.4.1       rstudioapi_0.10  R6_2.4.0         nlme_3.1-140    
    ## [53] compiler_3.6.0

    citation("tidyverse") 

    ## 
    ##   Wickham et al., (2019). Welcome to the tidyverse. Journal of
    ##   Open Source Software, 4(43), 1686,
    ##   https://doi.org/10.21105/joss.01686
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Article{,
    ##     title = {Welcome to the {tidyverse}},
    ##     author = {Hadley Wickham and Mara Averick and Jennifer Bryan and Winston Chang and Lucy D'Agostino McGowan and Romain François and Garrett Grolemund and Alex Hayes and Lionel Henry and Jim Hester and Max Kuhn and Thomas Lin Pedersen and Evan Miller and Stephan Milton Bache and Kirill Müller and Jeroen Ooms and David Robinson and Dana Paige Seidel and Vitalie Spinu and Kohske Takahashi and Davis Vaughan and Claus Wilke and Kara Woo and Hiroaki Yutani},
    ##     year = {2019},
    ##     journal = {Journal of Open Source Software},
    ##     volume = {4},
    ##     number = {43},
    ##     pages = {1686},
    ##     doi = {10.21105/joss.01686},
    ##   }
