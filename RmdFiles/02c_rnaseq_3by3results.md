The figures made from this script were compiled in Adobe.

<img src="../figures/02_RNAseq_3by3/02_RNAseq-01.png" width="1370" />

    library(ggplot2) ## for awesome plots!

    ## Warning: package 'ggplot2' was built under R version 3.3.2

    library(cowplot) ## for some easy to use themes
    library(dplyr) ## for filtering and selecting rows
    library(car) ## stats

    ## Warning: package 'car' was built under R version 3.3.2

    library(VennDiagram) ## venn diagrams
    library(pheatmap) ## awesome heatmaps
    library(viridis) # for awesome color pallette
    library(reshape2) ## for melting dataframe
    library(DESeq2) ## for gene expression analysis

    ## Warning: package 'DESeq2' was built under R version 3.3.2

    ## Warning: package 'S4Vectors' was built under R version 3.3.2

    ## Warning: package 'GenomicRanges' was built under R version 3.3.2

    ## Warning: package 'GenomeInfoDb' was built under R version 3.3.2

    library(edgeR)  ## for basic read counts status

    ## Warning: package 'edgeR' was built under R version 3.3.2

    ## Warning: package 'limma' was built under R version 3.3.2

    library(magrittr) ## to use the weird pipe
    library(genefilter)  ## for PCA fuction
    library(ggrepel) ## for labeling volcano plot

    ## Warning: package 'ggrepel' was built under R version 3.3.2

    library(colorblindr)

    ## load functions 
    source("figureoptions.R")
    source("functions_RNAseq.R")


    ## set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/02_RNAseq_3by3/')

    colData <- read.csv("../data/02a_colData.csv", header = T)
    countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)

Now, we are create differential gene expression object and remove genes
with 0 counts. Before filtering, there are 22,485 genes in the object.
After filtering genes with 0 counts, we will be left with 17,746 genes
that are expressed in a least 1 sample. Then, we can caluate the size
factors, estimate gene dispersion estimates, fit a model, test for
outliers, and remove outliers.

    ## create DESeq object using the factors Punch and APA
    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ Punch + APA + Punch*APA)

    ## it appears that the last variable in the design formula, 'APA',
    ##   has a factor level, 'Control', which is not the reference level. we recommend
    ##   to use factor(...,levels=...) or relevel() to set this as the reference level
    ##   before proceeding. for more information, please see the 'Note on factor levels'
    ##   in vignette('DESeq2').

    ## DESeq2 1.3.7 specify the factor levels
    dds$Punch <- factor(dds$Punch, levels=c("DG","CA3", "CA1"))
    dds$APA <- factor(dds$APA, levels=c("Control", "Consistent", "Conflict"))


    dds # view the DESeq object - note numnber of genes

    ## class: DESeqDataSet 
    ## dim: 22485 44 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(22485): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(44): 143A-CA3-1 143A-DG-1 ... 148B-CA3-4 148B-DG-4
    ## colData names(8): RNAseqID Mouse ... APA APA2

    ## DESeq2 1.3.6 Pre-filtering genes with 0 counts
    dds <- dds[ rowSums(counts(dds)) > 1, ] 

    dds # view the DESeq object - note numnber of genes

    ## class: DESeqDataSet 
    ## dim: 17674 44 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(17674): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(44): 143A-CA3-1 143A-DG-1 ... 148B-CA3-4 148B-DG-4
    ## colData names(8): RNAseqID Mouse ... APA APA2

    # dim: 17674 44 
    # 17,674 genes and 44 samples

    ## DESeq2 1.4  Differential expression analysi
    dds <- DESeq(dds)

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

    ## for variance stablized gene expression and log transformed data
    rld <- rlog(dds, blind=FALSE)
    #vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
    #vsd.fast <- vst(dds, blind=FALSE)

    contrast1 <- resvals(contrastvector = c("Punch", "CA1", "DG"), mypval = 0.05) #2497

    ## [1] 2497

    contrast2 <- resvals(contrastvector = c("Punch", "CA1", "CA3"), mypval = 0.05) #1803

    ## [1] 1803

    contrast3 <- resvals(contrastvector = c("Punch", "CA3", "DG"), mypval = 0.05) #3445

    ## [1] 3445

    contrast4 <- resvals(contrastvector = c("APA", "Consistent", "Control"), mypval = 0.05) #95

    ## [1] 95

    contrast5 <- resvals(contrastvector = c("APA", "Conflict", "Control"), mypval = 0.05) #42

    ## [1] 42

    contrast6 <- resvals(contrastvector = c("APA", "Conflict", "Consistent"), mypval = 0.05) # 0 

    ## [1] 0

    #create a new DF with the gene counts
    ## note: contrast1 had 0 differentially expressed genes, so it is not included 
    rldpadjs <- assay(rld)
    rldpadjs <- cbind(rldpadjs, contrast1, contrast2, contrast3, contrast4, contrast5, contrast6)
    rldpadjs <- as.data.frame(rldpadjs)
    rldpadjs <- rldpadjs[ , grepl( "padj" , names( rldpadjs ) ) ]
    head(rldpadjs)

    ##               padjPunchCA1DG padjPunchCA1CA3 padjPunchCA3DG
    ## 0610007P14Rik      0.9959885       1.0000000      0.9953218
    ## 0610009B22Rik      0.8920416       0.4918862      0.1034590
    ## 0610009L18Rik      0.6926063       0.6752603      0.9767398
    ## 0610009O20Rik      0.7325075       0.9860516      0.4348118
    ## 0610010F05Rik      0.7118161       0.1823768      0.3936081
    ## 0610010K14Rik      0.9959885       0.5247152      0.2851455
    ##               padjAPAConsistentControl padjAPAConflictControl
    ## 0610007P14Rik                        1                      1
    ## 0610009B22Rik                        1                      1
    ## 0610009L18Rik                        1                      1
    ## 0610009O20Rik                        1                      1
    ## 0610010F05Rik                        1                      1
    ## 0610010K14Rik                        1                      1
    ##               padjAPAConflictConsistent
    ## 0610007P14Rik                         1
    ## 0610009B22Rik                         1
    ## 0610009L18Rik                         1
    ## 0610009O20Rik                         1
    ## 0610010F05Rik                         1
    ## 0610010K14Rik                         1

Supplementray P-value distributions
-----------------------------------

Here, the goal is the analyze the distribution of pvalues to see if they
are randomly distributed or if that is a tendency towards and increase
or decrease of low pvalues. There, I'm showing the pval and adjusted
pvale (padj) for all for two-way comparision.

Venn Diagrams of DEgenes
------------------------

Now, we count the number of differnetially expressed genes (according to
padj) and plot some venn diagrams.

    venn1 <- row.names(rldpadjs[rldpadjs[1] <0.05 & !is.na(rldpadjs[1]),])
    venn2 <- row.names(rldpadjs[rldpadjs[2] <0.05 & !is.na(rldpadjs[2]),])
    venn3 <- row.names(rldpadjs[rldpadjs[3] <0.05 & !is.na(rldpadjs[3]),])
    venn4 <- row.names(rldpadjs[rldpadjs[4] <0.05 & !is.na(rldpadjs[4]),])
    venn5 <- row.names(rldpadjs[rldpadjs[5] <0.50 & !is.na(rldpadjs[5]),])

    candidates <- list("Control vs Consistent" = venn4, "DG vs CA3" = venn3  ,"Control vs Conflict" = venn5, "DG vs CA1" = venn1 )

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

![](../figures/02_RNAseq_3by3/venndiagram4-1.png)

    ## other 3 way

    candidates <- list("DG vs CA1" = venn1, "Control vs Consistent" = venn4,"DG vs CA3" = venn3,  "CA3 vs CA1" = venn2)

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

![](../figures/02_RNAseq_3by3/venndiagram3-1.png)

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

![](../figures/02_RNAseq_3by3/venndiagram3brainregions-1.png)

### PCA

    pcadata <- pcadataframe(rld, intgroup=c("Punch","APA"), returnData=TRUE)
    percentVar <- round(100 * attr(pcadata, "percentVar"))
    percentVar

    ## [1] 49 21  5  3  2  1  1  1  1

    aov1 <- aov(PC1 ~ APA, data=pcadata)
    summary(aov1) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## APA          2     51    25.4   0.059  0.943
    ## Residuals   41  17649   430.5

    TukeyHSD(aov1, which = "APA") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ APA, data = pcadata)
    ## 
    ## $APA
    ##                          diff       lwr      upr     p adj
    ## Consistent-Control  1.6866401 -18.41353 21.78681 0.9773221
    ## Conflict-Control    2.3622622 -15.04500 19.76952 0.9418272
    ## Conflict-Consistent 0.6756221 -20.87941 22.23066 0.9968027

    aov2 <- aov(PC2 ~ APA, data=pcadata)
    summary(aov2) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## APA          2    170    85.2   0.479  0.623
    ## Residuals   41   7297   178.0

    TukeyHSD(aov2, which = "APA") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ APA, data = pcadata)
    ## 
    ## $APA
    ##                          diff        lwr       upr     p adj
    ## Consistent-Control   4.141820  -8.782801 17.066442 0.7177431
    ## Conflict-Control    -1.312513 -12.505564  9.880538 0.9562190
    ## Conflict-Consistent -5.454333 -19.314446  8.405779 0.6078845

    pcadata$APA <- factor(pcadata$APA, levels=c("Control", "Consistent", "Conflict"))

    plotPCs(pcadata, 1, 2, aescolor = pcadata$Punch, colorname = "Punch", aesshape = pcadata$APA, shapename = "APA",  colorvalues = colorvalPunch)

    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.
    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.

![](../figures/02_RNAseq_3by3/pca-1.png)

    plotPCs(pcadata, 4, 2, aescolor = pcadata$APA, colorname = "APA", aesshape = pcadata$Punch, shapename = "Punch",  colorvalues = colorvalAPA)

    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.
    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.

![](../figures/02_RNAseq_3by3/pca-2.png)

    # pdf the same pca plots descripbed above of the above
    pdf(file="../figures/02_RNAseq_3by3/PCA12.pdf", width=4.5, height=3)
    PCA12 <- plotPCs(pcadata, 1, 2, aescolor = pcadata$Punch, colorname = " ", aesshape = pcadata$APA, shapename = " ",  colorvalues = colorvalPunch)
    plot(PCA12)

    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.
    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.

    dev.off()

    ## quartz_off_screen 
    ##                 2

    pdf(file="../figures/02_RNAseq_3by3/PCA42.pdf", width=4.5, height=3)
    PCA42 <- plotPCs(pcadata, 4, 2, aescolor = pcadata$APA, colorname = "APA2", aesshape = pcadata$Punch, shapename = "Punch",  colorvalues = colorvalAPA)
    plot(PCA42)

    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.
    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.

    dev.off()

    ## quartz_off_screen 
    ##                 2

Heatmap
-------

    DEGes <- assay(rld)
    DEGes <- cbind(DEGes, contrast1, contrast2, contrast3, contrast4, contrast5, contrast6)
    DEGes <- as.data.frame(DEGes) # convert matrix to dataframe
    DEGes$rownames <- rownames(DEGes)  # add the rownames to the dataframe
    DEGes$padjmin <- with(DEGes, pmin(padjPunchCA1DG, padjPunchCA1CA3, padjPunchCA3DG, padjAPAConsistentControl, padjAPAConflictControl, padjAPAConflictConsistent)) # create new col with min padj
    DEGes <- DEGes %>% filter(padjmin < 0.05)
    rownames(DEGes) <- DEGes$rownames
    drop.cols <-colnames(DEGes[,grep("padj|pval|rownames", colnames(DEGes))])
    DEGes <- DEGes %>% dplyr::select(-one_of(drop.cols))
    DEGes <- as.matrix(DEGes)
    DEGes <- DEGes - rowMeans(DEGes)
    head(DEGes)

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

    ## the heatmap annotation file
    df <- as.data.frame(colData(dds)[,c("Punch","APA")]) ## matrix to df


    rownames(df) <- names(countData)
    ann_colors <- ann_colors2 #use 

    # make sure the data is a matrix
    DEGes <- as.matrix(DEGes) 

    # set color breaks
    paletteLength <- 30
    myBreaks <- c(seq(min(DEGes), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(DEGes)/paletteLength, max(DEGes), length.out=floor(paletteLength/2)))

    pheatmap(DEGes, show_colnames=T, show_rownames = F,
             annotation_col=df, annotation_colors = ann_colors,
             treeheight_row = 0, treeheight_col = 25,
             fontsize = 11, 
             #width=4.5, height=3,
             border_color = "grey60" ,
             color = viridis(30),
             cellwidth = 8, 
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation" 
             )

![](../figures/02_RNAseq_3by3/heatmap-1.png)

    # for adobe
    pheatmap(DEGes, show_colnames=F, show_rownames = F,
             annotation_col=df, annotation_colors = ann_colors,
             treeheight_row = 0, treeheight_col = 50,
             fontsize = 10, 
             #width=4.5, height=3,
             border_color = "grey60" ,
             color = viridis(30),
             cellwidth = 8, 
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation",
             filename = "../figures/02_RNAseq_3by3/HeatmapPadj-1.pdf"
             )



    rownames(df) <- names(countData)
    ann_colors <- ann_colors2 #use 

    # make sure the data is a matrix
    DEGes <- as.matrix(DEGes) 

    # set color breaks
    paletteLength <- 30
    myBreaks <- c(seq(min(DEGes), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(DEGes)/paletteLength, max(DEGes), length.out=floor(paletteLength/2)))

    pheatmap(DEGes, show_colnames=T, show_rownames = F,
             annotation_col=df, annotation_colors = ann_colors,
             treeheight_row = 0, treeheight_col = 25,
             fontsize = 11, 
             #width=4.5, height=3,
             border_color = "grey60" ,
             color = viridis(30),
             cellwidth = 8, 
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation" 
             )

![](../figures/02_RNAseq_3by3/heatmap-2.png)

    # for adobe
    pheatmap(DEGes, show_colnames=F, show_rownames = F,
             annotation_col=df, annotation_colors = ann_colors,
             treeheight_row = 0, treeheight_col = 50,
             fontsize = 10, 
             #width=4.5, height=3,
             border_color = "grey60" ,
             color = viridis(30),
             cellwidth = 8, 
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation",
             filename = "../figures/02_RNAseq_3by3/HeatmapPadj-1.pdf"
             )

Volcanos plots and and gene lists
---------------------------------

    # gene lists
    res <- results(dds, contrast =c("Punch", "CA1", "DG"), independentFiltering = F)
    resOrdered <- res[order(res$padj),]
    head(resOrdered, 10)

    ## log2 fold change (MLE): Punch CA1 vs DG 
    ## Wald test p-value: Punch CA1 vs DG 
    ## DataFrame with 10 rows and 6 columns
    ##          baseMean log2FoldChange     lfcSE      stat       pvalue
    ##         <numeric>      <numeric> <numeric> <numeric>    <numeric>
    ## Fam163b 632.63230      -3.918682 0.2476783 -15.82166 2.205934e-56
    ## Pitpnm2 157.76399      -2.710742 0.1724904 -15.71532 1.187738e-55
    ## Ncald   127.62509      -4.453148 0.2862203 -15.55846 1.394069e-54
    ## Map4     98.71267       2.187816 0.1451653  15.07121 2.504944e-51
    ## Actr3b  147.01919       1.850844 0.1274328  14.52408 8.527759e-48
    ## Prkcg   534.81768       2.627146 0.1839476  14.28203 2.832248e-46
    ## Stum    393.20190       1.646700 0.1166105  14.12137 2.804944e-45
    ## Wfs1    186.06663       5.102024 0.3663791  13.92553 4.432272e-44
    ## Pex5l    46.20933       3.275083 0.2388209  13.71355 8.423359e-43
    ## Pde2a   104.87140       2.216737 0.1623222  13.65640 1.849049e-42
    ##                 padj
    ##            <numeric>
    ## Fam163b 3.897665e-52
    ## Pitpnm2 1.049307e-51
    ## Ncald   8.210599e-51
    ## Map4    1.106497e-47
    ## Actr3b  3.013540e-44
    ## Prkcg   8.340498e-43
    ## Stum    7.080080e-42
    ## Wfs1    9.789226e-41
    ## Pex5l   1.653693e-39
    ## Pde2a   3.267085e-39

    data <- data.frame(gene = row.names(res), pvalue = -log10(res$padj), lfc = res$log2FoldChange)
    data <- na.omit(data)
    data <- data %>%
      mutate(color = ifelse(data$lfc > 0 & data$pvalue > 1.3, 
                            yes = "CA1", 
                            no = ifelse(data$lfc < 0 & data$pvalue > 1.3, 
                                        yes = "DG", 
                                        no = "none")))
    top_labelled <- top_n(data, n = 3, wt = pvalue)

    # Color corresponds to fold change directionality
    colored <- ggplot(data, aes(x = lfc, y = pvalue)) + 
      geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
      theme_bw(base_size = 16) + # clean up theme
      theme(legend.position = "none") + # remove legend 
      xlim(c(-10, 10)) +  
      ylim(c(0, 100)) +   
      scale_color_manual(values = c("CA1" = "#7570b3",
                                    "DG" = "#d95f02", 
                                    "none" = "#bdbdbd")) + theme(panel.grid.minor=element_blank(),
               panel.grid.major=element_blank()) + 
      theme(axis.title.x = element_blank())+ 
      theme(axis.title.y = element_blank()) + 
      geom_text_repel(data = top_labelled, 
                              mapping = aes(label = gene), 
                              size = 3,
                              fontface = 'bold', 
                              color = 'black',
                              box.padding = unit(0.5, "lines"),
                              point.padding = unit(0.5, "lines"))

    colored

![](../figures/02_RNAseq_3by3/volcanos-1.png)

    #cvd_grid(colored) # to view plot for color blind 
    pdf(file="../figures/02_RNAseq_3by3/CA1DG.pdf", width=3, height=3)
    plot(colored)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    res <- results(dds, contrast =c("Punch", "CA1", "CA3"), independentFiltering = F)
    resOrdered <- res[order(res$padj),]
    head(resOrdered, 10)

    ## log2 fold change (MLE): Punch CA1 vs CA3 
    ## Wald test p-value: Punch CA1 vs CA3 
    ## DataFrame with 10 rows and 6 columns
    ##         baseMean log2FoldChange     lfcSE      stat        pvalue
    ##        <numeric>      <numeric> <numeric> <numeric>     <numeric>
    ## Itpka   712.2031       2.879520 0.1316754  21.86832 5.203945e-106
    ## Doc2b   351.7638       6.385453 0.3440064  18.56202  6.520976e-77
    ## Syn2    367.3440      -2.000931 0.1180102 -16.95558  1.750175e-64
    ## Nptxr   943.4033      -2.654819 0.1620732 -16.38037  2.641336e-60
    ## Bcr     343.0037       2.329822 0.1439082  16.18964  5.967757e-59
    ## Fibcd1  361.2475       7.244385 0.4473197  16.19510  5.461311e-59
    ## Celf4   315.2192      -1.593626 0.1098523 -14.50698  1.094343e-47
    ## Ncald   127.6251      -4.202885 0.2906670 -14.45945  2.185272e-47
    ## Bcl11b  134.7735       3.373066 0.2341499  14.40559  4.772594e-47
    ## Wfs1    186.0666       5.343254 0.3725146  14.34375  1.165872e-46
    ##                 padj
    ##            <numeric>
    ## Itpka  9.194851e-102
    ## Doc2b   5.760956e-73
    ## Syn2    1.030794e-60
    ## Nptxr   1.166744e-56
    ## Bcr     1.757405e-55
    ## Fibcd1  1.757405e-55
    ## Celf4   2.762279e-44
    ## Ncald   4.826446e-44
    ## Bcl11b  9.369662e-44
    ## Wfs1    2.059980e-43

    data <- data.frame(gene = row.names(res), pvalue = -log10(res$padj), lfc = res$log2FoldChange)
    data <- na.omit(data)
    head(data)

    ##            gene      pvalue         lfc
    ## 1 0610007P14Rik 0.000000000  0.02802707
    ## 2 0610009B22Rik 0.308135333 -0.75347715
    ## 3 0610009L18Rik 0.170528796 -0.95713263
    ## 4 0610009O20Rik 0.006100336  0.11554405
    ## 5 0610010F05Rik 0.739030388 -0.73711178
    ## 6 0610010K14Rik 0.280076346  1.06103762

    data <- data %>%
      mutate(color = ifelse(data$lfc > 0 & data$pvalue > 1.3, 
                            yes = "CA1", 
                            no = ifelse(data$lfc < 0 & data$pvalue > 1.3, 
                                        yes = "CA3", 
                                        no = "none")))
    top_labelled <- top_n(data, n = 3, wt = pvalue)
    colored <- ggplot(data, aes(x = lfc, y = pvalue)) + 
      geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
      theme_bw(base_size = 16) + # clean up theme
      theme(legend.position = "none") + # remove legend 
      xlim(c(-10, 10)) +  
      ylim(c(0, 100)) +  
      scale_color_manual(values = c("CA1" = "#7570b3",
                                    "CA3" = "#1b9e77", 
                                    "none" = "#bdbdbd")) + theme(panel.grid.minor=element_blank(),
               panel.grid.major=element_blank()) + 
      theme(axis.title.x = element_blank()) + 
      theme(axis.title.y = element_blank()) + 
      geom_text_repel(data = top_labelled, 
                              mapping = aes(label = gene), 
                              size = 3,
                              fontface = 'bold', 
                              color = 'black',
                              box.padding = unit(0.5, "lines"),
                              point.padding = unit(0.5, "lines"))

    colored

    ## Warning: Removed 1 rows containing missing values (geom_text_repel).

![](../figures/02_RNAseq_3by3/volcanos-2.png)

    pdf(file="../figures/02_RNAseq_3by3/CA1CA3.pdf", width=3, height=3)
    plot(colored)

    ## Warning: Removed 1 rows containing missing values (geom_text_repel).

    dev.off()

    ## quartz_off_screen 
    ##                 2

    res <- results(dds, contrast =c("Punch", "CA3", "DG"), independentFiltering = F)
    resOrdered <- res[order(res$padj),]
    head(resOrdered, 10)

    ## log2 fold change (MLE): Punch CA3 vs DG 
    ## Wald test p-value: Punch CA3 vs DG 
    ## DataFrame with 10 rows and 6 columns
    ##           baseMean log2FoldChange      lfcSE      stat       pvalue
    ##          <numeric>      <numeric>  <numeric> <numeric>    <numeric>
    ## Fam163b  632.63230      -5.515648 0.26028028 -21.19119 1.151616e-99
    ## Pitpnm2  157.76399      -2.974671 0.16632932 -17.88423 1.565051e-71
    ## C1ql3    288.38305      -7.669488 0.43971862 -17.44181 3.973016e-68
    ## Doc2b    351.76381      -5.886124 0.33939153 -17.34317 2.221555e-67
    ## Adcy1   2632.63884      -3.554641 0.22457803 -15.82809 1.991737e-56
    ## Gnao1    164.33213       1.531840 0.09677608  15.82870 1.972523e-56
    ## Syn2     367.34401       1.729652 0.11324656  15.27333 1.151397e-52
    ## Tiam1     45.93648      -3.906561 0.26275706 -14.86758 5.351455e-50
    ## Fam84a   296.31225      -1.923974 0.13124976 -14.65888 1.182120e-48
    ## Ahi1     241.16807      -1.626489 0.11177026 -14.55207 5.665820e-48
    ##                 padj
    ##            <numeric>
    ## Fam163b 2.034791e-95
    ## Pitpnm2 1.382644e-67
    ## C1ql3   2.339974e-64
    ## Doc2b   9.813162e-64
    ## Adcy1   5.865333e-53
    ## Gnao1   5.865333e-53
    ## Syn2    2.906291e-49
    ## Tiam1   1.181936e-46
    ## Fam84a  2.320764e-45
    ## Ahi1    1.001094e-44

    data <- data.frame(gene = row.names(res), pvalue = -log10(res$padj), lfc = res$log2FoldChange)
    data <- na.omit(data)
    head(data)

    ##            gene     pvalue         lfc
    ## 1 0610007P14Rik 0.00203650  0.01805968
    ## 2 0610009B22Rik 0.98523153  1.11049704
    ## 3 0610009L18Rik 0.01022111  0.11022248
    ## 4 0610009O20Rik 0.36169871 -0.39180738
    ## 5 0610010F05Rik 0.40493601  0.39920496
    ## 6 0610010K14Rik 0.54493355 -1.15405296

    data <- data %>%
      mutate(color = ifelse(data$lfc > 0 & data$pvalue > 1.3, 
                            yes = "CA3", 
                            no = ifelse(data$lfc < 0 & data$pvalue > 1.3, 
                                        yes = "DG", 
                                        no = "none")))
    top_labelled <- top_n(data, n = 3, wt = pvalue)
    colored <- ggplot(data, aes(x = lfc, y = pvalue)) + 
      geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
      theme_bw(base_size = 16) + # clean up theme
      theme(legend.position = "none") + # remove legend 
      xlim(c(-10, 10)) +  
      ylim(c(0, 100)) +  
      scale_color_manual(values = c("CA3" = "#1b9e77",
                                    "DG" = "#d95f02", 
                                    "none" = "#bdbdbd")) + theme(panel.grid.minor=element_blank(),
               panel.grid.major=element_blank()) + 
      theme(axis.title.x = element_blank()) + 
      theme(axis.title.y = element_blank()) + 
      geom_text_repel(data = top_labelled, 
                              mapping = aes(label = gene), 
                              size = 3,
                              fontface = 'bold', 
                              color = 'black',
                              box.padding = unit(0.5, "lines"),
                              point.padding = unit(0.5, "lines"))

    colored

![](../figures/02_RNAseq_3by3/volcanos-3.png)

    #cvd_grid(colored)
    pdf(file="../figures/02_RNAseq_3by3/DGCA3.pdf", width=3, height=3)
    plot(colored)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    res <- results(dds, contrast =c("APA", "Consistent", "Control"), independentFiltering = F)
    resOrdered <- res[order(res$padj),]
    head(resOrdered, 10)

    ## log2 fold change (MLE): APA Consistent vs Control 
    ## Wald test p-value: APA Consistent vs Control 
    ## DataFrame with 10 rows and 6 columns
    ##         baseMean log2FoldChange     lfcSE      stat       pvalue
    ##        <numeric>      <numeric> <numeric> <numeric>    <numeric>
    ## Plk2   697.03923      2.0089874 0.2656720  7.561910 3.971929e-14
    ## Frmd6  116.84596      3.0501052 0.4208546  7.247408 4.248219e-13
    ## Sgk1    27.18260      2.1618090 0.3096130  6.982294 2.903992e-12
    ## Arc    436.24502      2.4535807 0.3641878  6.737130 1.615458e-11
    ## Fbxo33  74.39061      2.8002536 0.4141262  6.761837 1.362532e-11
    ## Smad7   25.81154      2.9386801 0.4333667  6.781047 1.193080e-11
    ## Dnaja1 106.22420      0.9956356 0.1556782  6.395470 1.600539e-10
    ## Fosl2  298.06084      2.3892668 0.3796235  6.293781 3.098256e-10
    ## Ubc    469.30786      1.0276446 0.1633569  6.290793 3.158479e-10
    ## Homer1  42.54046      2.6587543 0.4283628  6.206781 5.408076e-10
    ##                padj
    ##           <numeric>
    ## Plk2   7.018002e-10
    ## Frmd6  3.753089e-09
    ## Sgk1   1.710354e-08
    ## Arc    4.757254e-08
    ## Fbxo33 4.757254e-08
    ## Smad7  4.757254e-08
    ## Dnaja1 4.039989e-07
    ## Fosl2  6.200797e-07
    ## Ubc    6.200797e-07
    ## Homer1 9.450620e-07

    data <- data.frame(gene = row.names(res),
                       pvalue = -log10(res$padj), 
                       lfc = res$log2FoldChange)
    data <- na.omit(data)
    data <- data %>%
      mutate(color = ifelse(data$lfc > 0 & data$pvalue > 1.3, 
                            yes = "Consistent", 
                            no = ifelse(data$lfc < 0 & data$pvalue > 1.3, 
                                        yes = "Control", 
                                        no = "none")))
    top_labelled <- top_n(data, n = 3, wt = pvalue)
    colored <- ggplot(data, aes(x = lfc, y = pvalue)) + 
      geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
      theme_bw(base_size = 16) + # clean up theme
      theme(legend.position = "none") + # remove legend 
      xlim(c(-10, 10)) +  
      ylim(c(0, 10)) + 
      scale_color_manual(values = c("Consistent" = "#f4a582",
                                    "Control" = "#bababa", 
                                    "none" = "#bdbdbd")) + theme(panel.grid.minor=element_blank(),
               panel.grid.major=element_blank()) + 
      theme(axis.title.x = element_blank())+ 
      theme(axis.title.y = element_blank()) + 
      geom_text_repel(data = top_labelled, 
                              mapping = aes(label = gene), 
                              size = 3,
                              fontface = 'bold', 
                              color = 'black',
                              box.padding = unit(0.5, "lines"),
                              point.padding = unit(0.5, "lines"))

    colored

![](../figures/02_RNAseq_3by3/volcanos-4.png)

    pdf(file="../figures/02_RNAseq_3by3/consistentcontrol.pdf", width=3, height=3)
    plot(colored)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    ####

    res <- results(dds, contrast =c("APA", "Conflict", "Control"), independentFiltering = F)
    resOrdered <- res[order(res$padj),]
    head(resOrdered, 10)

    ## log2 fold change (MLE): APA Conflict vs Control 
    ## Wald test p-value: APA Conflict vs Control 
    ## DataFrame with 10 rows and 6 columns
    ##          baseMean log2FoldChange     lfcSE      stat       pvalue
    ##         <numeric>      <numeric> <numeric> <numeric>    <numeric>
    ## Plk2    697.03923       1.412284 0.2249393  6.278511 3.418308e-10
    ## Arc     436.24502       1.839986 0.3087950  5.958600 2.544079e-09
    ## Smad7    25.81154       2.254145 0.3747463  6.015123 1.797506e-09
    ## Sgk1     27.18260       1.566089 0.2673891  5.856969 4.713923e-09
    ## Slc16a1  51.34226       1.933710 0.3347108  5.777257 7.592832e-09
    ## Errfi1   75.94706       1.463933 0.2637376  5.550718 2.844984e-08
    ## Frmd6   116.84596       1.965433 0.3590953  5.473290 4.417557e-08
    ## Fosl2   298.06084       1.723689 0.3213667  5.363621 8.157011e-08
    ## Ptgs2    66.22778       2.033994 0.3816851  5.328984 9.876359e-08
    ## Fbxo33   74.39061       1.835648 0.3540087  5.185319 2.156453e-07
    ##                 padj
    ##            <numeric>
    ## Plk2    6.039808e-06
    ## Arc     1.498378e-05
    ## Smad7   1.498378e-05
    ## Sgk1    2.082258e-05
    ## Slc16a1 2.683155e-05
    ## Errfi1  8.378004e-05
    ## Frmd6   1.115055e-04
    ## Fosl2   1.801578e-04
    ## Ptgs2   1.938949e-04
    ## Fbxo33  3.810237e-04

    data <- data.frame(gene = row.names(res),
                       pvalue = -log10(res$padj), 
                       lfc = res$log2FoldChange)
    data <- na.omit(data)
    data <- data %>%
      mutate(color = ifelse(data$lfc > 0 & data$pvalue > 1.3, 
                            yes = "Conflict", 
                            no = ifelse(data$lfc < 0 & data$pvalue > 1.3, 
                                        yes = "Conflict", 
                                        no = "none")))
    top_labelled <- top_n(data, n = 3, wt = pvalue)
    # Color corresponds to fold change directionality
    colored <- ggplot(data, aes(x = lfc, y = pvalue)) + 
      geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
      theme_bw(base_size = 16) + # clean up theme
      theme(legend.position = "none") + # remove legend 
      xlim(c(-10, 10)) +  
      ylim(c(0, 10)) +  
      scale_color_manual(values = c("Conflict" = "#ca0020",
                                    "Control" = "#404040", 
                                    "none" = "#bdbdbd")) + theme(panel.grid.minor=element_blank(),
               panel.grid.major=element_blank()) + 
      theme(axis.title.x = element_blank())+ 
      theme(axis.title.y = element_blank()) + 
      geom_text_repel(data = top_labelled, 
                              mapping = aes(label = gene), 
                              size = 3,
                              fontface = 'bold', 
                              color = 'black',
                              box.padding = unit(0.5, "lines"),
                              point.padding = unit(0.5, "lines"))

    colored

![](../figures/02_RNAseq_3by3/volcanos-5.png)

    pdf(file="../figures/02_RNAseq_3by3/conflictcontrol.pdf", width=3, height=3)
    plot(colored)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    #write.csv(rldpadjs, file = "../data/02c_rldpadjs.csv", row.names = T)
    #write.csv(DEGes, file = "../data/02c_DEGes.csv", row.names = T)
    #write.csv(df, file = "../data/02c_df.csv", row.names = F)
    #write.csv(pcadata, file = "../data/02c_pcadata.csv", row.names = F)
    #write.table(percentVar, file = "../data/02c_percentVar.txt")
