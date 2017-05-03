    library(dplyr) ## for filtering and selecting rows
    library(plyr) ## for renmaing factors
    library(reshape2) ## for melting dataframe
    library(DESeq2) ## for gene expression analysis
    library(edgeR)  ## for basic read counts status
    library(magrittr) ## to use the weird pipe
    library(genefilter)  ## for PCA fuction
    library(ggplot2) ## for awesome plots!
    library(pheatmap) ## awesome heatmaps


    ## Functions
    source("functions_RNAseq.R")
    source("resvalsfunction.R")
    source("figureoptions.R")

    ## set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/02_RNAseq/')

CA1 only differential gene expression
-------------------------------------

Zero genes have padj &lt; 0.05

    countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)
    colData <- read.csv("../data/02a_colData.csv", header = T)
    colData$APA <- factor(colData$APA, levels=c("Control", "Consistent", "Conflict"))

    colData <- colData %>% 
      filter(Punch == "CA1") 
    savecols <- as.character(colData$RNAseqID) 
    savecols <- as.vector(savecols) 
    countData <- countData %>% dplyr::select(one_of(savecols)) 

    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ APA )
    dds$APA <- factor(dds$APA, levels=c("Control", "Consistent", "Conflict"))
    dds <- dds[ rowSums(counts(dds)) > 1, ] 
    dds # dim: 16467 15  

    ## class: DESeqDataSet 
    ## dim: 16467 15 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(16467): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(15): 143B-CA1-1 143C-CA1-1 ... 148A-CA1-3 148B-CA1-4
    ## colData names(7): RNAseqID Mouse ... APA ID

    dds <- DESeq(dds)

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 334 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    rld <- rlog(dds, blind=FALSE)

    contrast4 <- resvals(contrastvector = c("APA", "Consistent", "Control"), mypadj = 0.05) #0

    ## [1] 0

    contrast5 <- resvals(contrastvector = c("APA", "Conflict", "Control"), mypadj = 0.05) #0

    ## [1] 0

    contrast6 <- resvals(contrastvector = c("APA", "Conflict", "Consistent"), mypadj = 0.05) # 0 

    ## [1] 0

CA3 only differential gene expression
-------------------------------------

    countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)
    colData <- read.csv("../data/02a_colData.csv", header = T)
    colData$APA <- factor(colData$APA, levels=c("Control", "Consistent", "Conflict"))

    colData <- colData %>% 
      filter(Punch == "CA3") 
    savecols <- as.character(colData$RNAseqID) 
    savecols <- as.vector(savecols) 
    countData <- countData %>% dplyr::select(one_of(savecols)) 

    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ APA )
    dds$APA <- factor(dds$APA, levels=c("Control", "Consistent", "Conflict"))
    dds <- dds[ rowSums(counts(dds)) > 1, ] 
    dds # dim: 16467 15  

    ## class: DESeqDataSet 
    ## dim: 16208 13 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(16208): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(13): 143A-CA3-1 144A-CA3-2 ... 148A-CA3-3 148B-CA3-4
    ## colData names(7): RNAseqID Mouse ... APA ID

    dds <- DESeq(dds)

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    rld <- rlog(dds, blind=FALSE)

    contrast4 <- resvals(contrastvector = c("APA", "Consistent", "Control"), mypadj = 0.05) #0

    ## [1] 0

    contrast5 <- resvals(contrastvector = c("APA", "Conflict", "Control"), mypadj = 0.05) #0

    ## [1] 0

    contrast6 <- resvals(contrastvector = c("APA", "Conflict", "Consistent"), mypadj = 0.05) # 3 

    ## [1] 3

    countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)
    colData <- read.csv("../data/02a_colData.csv", header = T)
    colData$APA <- factor(colData$APA, levels=c("Control", "Consistent", "Conflict"))

    colData <- colData %>% 
      filter(Punch == "DG") 
    savecols <- as.character(colData$RNAseqID) 
    savecols <- as.vector(savecols) 
    countData <- countData %>% dplyr::select(one_of(savecols)) 

    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ APA )
    dds$APA <- factor(dds$APA, levels=c("Control", "Consistent", "Conflict"))
    dds <- dds[ rowSums(counts(dds)) > 1, ] 
    dds # dim: 16467 15  

    ## class: DESeqDataSet 
    ## dim: 16658 16 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(16658): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(16): 143A-DG-1 143B-DG-1 ... 148A-DG-3 148B-DG-4
    ## colData names(7): RNAseqID Mouse ... APA ID

    dds <- DESeq(dds)

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 217 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    rld <- rlog(dds, blind=FALSE)

    contrast4 <- resvals(contrastvector = c("APA", "Consistent", "Control"), mypadj = 0.05) #101

    ## [1] 101

    contrast5 <- resvals(contrastvector = c("APA", "Conflict", "Control"), mypadj = 0.05) #39

    ## [1] 39

    contrast6 <- resvals(contrastvector = c("APA", "Conflict", "Consistent"), mypadj = 0.05) # 0 

    ## [1] 0

    res <- results(dds, contrast =c("APA", "Consistent", "Control"), independentFiltering = F)
    with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DG Control - Consistent"))
    with(subset(res, log2FoldChange>0), points(log2FoldChange, -log10(pvalue), pch=20, col=c("#f4a582")))
    with(subset(res, log2FoldChange<0), points(log2FoldChange, -log10(pvalue), pch=20, col=c("#404040")))
    with(subset(res, padj>.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="grey"))

![](../figures/02_RNAseq/DGonly-1.png)

    res <- results(dds, contrast =c("APA", "Conflict", "Control"), independentFiltering = F)
    with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DG Control - Conflict"))
    with(subset(res, log2FoldChange>0), points(log2FoldChange, -log10(pvalue), pch=20, col=c("#ca0020")))
    with(subset(res, log2FoldChange<0), points(log2FoldChange, -log10(pvalue), pch=20, col=c("#404040")))
    with(subset(res, padj>.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="grey"))

![](../figures/02_RNAseq/DGonly-2.png)

    pcadata <- pcadataframe(rld, intgroup=c("Punch","APA"), returnData=TRUE)
    percentVar <- round(100 * attr(pcadata, "percentVar"))
    percentVar

    ## [1] 31 20 11  6  6  5  4  3  3

    pcadata <- pcadataframe(rld, intgroup=c("Punch","APA"), returnData=TRUE)
    percentVar <- round(100 * attr(pcadata, "percentVar"))
    aov1 <- aov(PC1 ~ APA, data=pcadata)
    summary(aov1) 

    ##             Df Sum Sq Mean Sq F value  Pr(>F)   
    ## APA          2  891.8   445.9   7.525 0.00675 **
    ## Residuals   13  770.3    59.3                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov1, which = "APA")

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ APA, data = pcadata)
    ## 
    ## $APA
    ##                           diff       lwr       upr     p adj
    ## Consistent-Control  -16.911835 -30.67185 -3.151817 0.0163925
    ## Conflict-Control    -13.437796 -25.02477 -1.850824 0.0230507
    ## Conflict-Consistent   3.474039 -11.36916 18.317242 0.8129880

    pcadata$Punch <- factor(pcadata$APA, levels=c("Control", "Consistent", "Conflict"))
    plotPCs(pcadata, 1, 2, aescolor = pcadata$APA, colorname = "APA", aesshape = pcadata$APA, shapename = "APA",  colorvalues = colorvalAPA)

    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.

    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.

![](../figures/02_RNAseq/DGonly-3.png)

    DEGes <- assay(rld)
    DEGes <- cbind(DEGes, contrast4, contrast5, contrast6)
    DEGes <- as.data.frame(DEGes) # convert matrix to dataframe
    DEGes$rownames <- rownames(DEGes)  # add the rownames to the dataframe
    DEGes$padjmin <- with(DEGes, pmin(padjAPAConsistentControl, padjAPAConflictControl, padjAPAConflictConsistent)) # create new col with min padj
    DEGes <- DEGes %>% filter(padjmin < 0.05)
    rownames(DEGes) <- DEGes$rownames
    drop.cols <-colnames(DEGes[,grep("padj|pval|rownames", colnames(DEGes))])
    DEGes <- DEGes %>% dplyr::select(-one_of(drop.cols))
    DEGes <- as.matrix(DEGes)
    DEGes <- DEGes - rowMeans(DEGes)
    head(DEGes)

    ##               143A-DG-1  143B-DG-1   143D-DG-3 144A-DG-2 144C-DG-2
    ## 1190002N15Rik 0.4455570 -1.2157030 -0.03475106  1.597700 1.3458024
    ## A830010M20Rik 0.1635088 -0.6078429 -0.33894975  1.139229 0.9555086
    ## Acan          0.8068617 -0.6044171 -0.55024770  1.027837 0.5205708
    ## Ankrd28       0.4531462 -0.3024109 -0.32540927  1.014265 0.4702013
    ## Arc           1.1604161 -0.7219387 -0.32721112  1.377827 1.9855491
    ## Arid5b        0.7316017 -0.4228392 -0.14309363  0.995551 0.8789801
    ##                144D-DG-2   145A-DG-2  145B-DG-1   146A-DG-2   146B-DG-2
    ## 1190002N15Rik -0.7292517  0.16765646 -0.4307213  0.38895475  0.66573864
    ## A830010M20Rik -0.4386699 -0.07480234 -0.4229666  0.27064657 -0.13171398
    ## Acan          -0.8016263  0.49549568 -0.1804890  0.01206221 -0.52452040
    ## Ankrd28       -0.2616318 -0.37720419 -0.6823164 -0.38907739 -0.05811206
    ## Arc           -1.8027171  0.09813946 -0.3341147  0.96853380 -0.84340626
    ## Arid5b        -0.4248934 -0.35016818 -0.2597698  0.04681877 -0.48094854
    ##               146C-DG-4  146D-DG-3 147C-DG-3  147D-DG-1   148A-DG-3
    ## 1190002N15Rik 0.8285891 -1.1320414 0.7129418 -0.8091523 -0.14966254
    ## A830010M20Rik 1.0159530 -0.7108966 0.6500946 -0.7030243 -0.13199863
    ## Acan          0.7003147 -0.4840008 1.0180239 -1.2202602  0.29333251
    ## Ankrd28       1.0828383 -1.1383214 1.0082297 -0.6743511  0.41821674
    ## Arc           0.9282108 -1.1627532 1.0877477 -2.3179783  0.15267970
    ## Arid5b        0.1811024 -0.9836882 0.7750028 -0.2809820  0.03062427
    ##                148B-DG-4
    ## 1190002N15Rik -1.6516569
    ## A830010M20Rik -0.6340759
    ## Acan          -0.5089366
    ## Ankrd28       -0.2380624
    ## Arc           -0.2489839
    ## Arid5b        -0.2932981

    df <- as.data.frame(colData(dds)[,c("APA", "Punch")]) #
    rownames(df) <- names(countData)
    ann_colors = ann_colors1 
    DEGes <- as.matrix(DEGes) 
    paletteLength <- 30
    myBreaks <- c(seq(min(DEGes), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(DEGes)/paletteLength, max(DEGes), length.out=floor(paletteLength/2)))

    pheatmap(DEGes, show_colnames=F, show_rownames = F,
             annotation_col=df, annotation_colors = ann_colors,
             treeheight_row = 0, treeheight_col = 25,
             fontsize = 11, 
             #width=4.5, height=3,
             border_color = "grey60" ,
             color = colorpalette,
             cellwidth = 8, 
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation" 
             )

![](../figures/02_RNAseq/DGonly-4.png)
