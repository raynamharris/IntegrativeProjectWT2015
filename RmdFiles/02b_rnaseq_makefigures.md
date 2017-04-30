    library(ggplot2) ## for awesome plots!

    ## Warning: package 'ggplot2' was built under R version 3.3.2

    library(cowplot) ## for some easy to use themes
    library(dplyr) ## for filtering and selecting rows
    library(car) ## stats

    ## Warning: package 'car' was built under R version 3.3.2

    library(VennDiagram) ## venn diagrams
    library(pheatmap) ## awesome heatmaps

    ## load functions 
    source("figureoptions.R")
    source("functions_RNAseq.R")

    ## set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/02_rnaseq/')

    colData <- read.csv("../data/02a_colData.csv", header = T)
    countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)
    rldpadjs <- read.csv("../data/02a_rldpadjs.csv", header = T, row.names = 1)
    DEGes <- read.csv("../data/02a_DEGes.csv", header = T, check.names = F, row.names = 1)
    df <- read.csv("../data/02a_df.csv", header = T)
    pcadata <- read.csv("../data/02a_pcadata.csv", header = T)

P-value distributions
---------------------

Here, the goal is the analyze the distribution of pvalues to see if they
are randomly distributed or if that is a tendency towards and increase
or decrease of low pvalues. There, I'm showing the pval and adjusted
pvale (padj) for all for two-way comparision.

    ggplot(rldpadjs, aes(x = padjPunchCA1DG)) + geom_histogram(binwidth = 0.05) + scale_y_log10()

    ## Warning: Removed 5 rows containing non-finite values (stat_bin).

![](../figures/02_rnaseq/pvaluedistribution-1.png)

    ggplot(rldpadjs, aes(x = padjPunchCA1CA3)) + geom_histogram(binwidth = 0.05) + scale_y_log10()

    ## Warning: Removed 5 rows containing non-finite values (stat_bin).

![](../figures/02_rnaseq/pvaluedistribution-2.png)

    ggplot(rldpadjs, aes(x = padjPunchCA3DG)) + geom_histogram(binwidth = 0.05) + scale_y_log10()

    ## Warning: Removed 5 rows containing non-finite values (stat_bin).

![](../figures/02_rnaseq/pvaluedistribution-3.png)

    ggplot(rldpadjs, aes(x = padjAPASameYoked)) + geom_histogram(binwidth = 0.05) + scale_y_log10()

    ## Warning: Removed 5 rows containing non-finite values (stat_bin).

![](../figures/02_rnaseq/pvaluedistribution-4.png)

    ggplot(rldpadjs, aes(x = padjAPAConflictYoked)) + geom_histogram(binwidth = 0.05) + scale_y_log10()

    ## Warning: Removed 5 rows containing non-finite values (stat_bin).

![](../figures/02_rnaseq/pvaluedistribution-5.png)

Venn Diagrams of DEgenes
------------------------

Now, we count the number of differnetially expressed genes (according to
padj) and plot some venn diagrams.

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

![](../figures/02_rnaseq/venndiagram-1.png)

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

![](../figures/02_rnaseq/venndiagram-2.png)

Heatmaps
--------

    source("figureoptions.R")

    rownames(df) <- names(countData)

    ann_colors = ann_colors1 #use 

    # make sure the data is a matrix
    DEGes <- as.matrix(DEGes) 

    # set color breaks
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
             #cellwidth = 12, 
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation" 
             )

![](../figures/02_rnaseq/heatmap-1.png)

    # for adobe
    pheatmap(DEGes, show_colnames=F, show_rownames = F,
             annotation_col=df, annotation_colors = ann_colors,
             treeheight_row = 0, treeheight_col = 25,
             fontsize = 11, 
             #width=4.5, height=3,
             border_color = "grey60" ,
             color = colorpalette,
             #cellwidth = 12, 
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation",
             filename = "../figures/02_RNAseq/HeatmapPadj-1.pdf"
             )

    ## pca plots

    ## I haven't been able to get this to work , so I hand coded the vector from a previous analysis of this
    #percentVar <- round(100 * attr(pcadata, "percentVar"))
    percentVar <- c(49,21,5,3,2,1,1,1,1)


    # separates brain regions
    plotPCs(pcadata, 1, 2, aescolor = pcadata$Punch, colorname = " ", aesshape = pcadata$APA, shapename = " ",  colorvalues = colorvalPunch)

    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.
    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.

![](../figures/02_rnaseq/pca-1.png)

    plotPCs(pcadata, 1, 2, aescolor = pcadata$APA, colorname = "APA", aesshape = pcadata$Punch, shapename = "Punch",  colorvalues = colorvalAPA)

    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.
    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.

![](../figures/02_rnaseq/pca-2.png)

    # PC4 significant for training
    plotPCs(pcadata, 2, 4, aescolor = pcadata$APA, colorname = "APA", aesshape = pcadata$Punch, shapename = "Punch",  colorvalues = colorvalAPA)

    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.
    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.

![](../figures/02_rnaseq/pca-3.png)

    # pdf the same pca plots descripbed above of the above
    pdf(file="../figures/02_RNAseq/PCA12.pdf", width=4.5, height=3)
    PCA12 <- plotPCs(pcadata, 1, 2, aescolor = pcadata$Punch, colorname = " ", aesshape = pcadata$APA, shapename = " ",  colorvalues = colorvalPunch)
    plot(PCA12)

    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.
    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.

    dev.off()

    ## quartz_off_screen 
    ##                 2

    pdf(file="../figures/02_RNAseq/PCA25.pdf", width=4.5, height=3)
    PCA25 <- plotPCs(pcadata, 2,4, aescolor = pcadata$APA, colorname = " ", aesshape = pcadata$Punch, shapename = " ",  colorvalues = colorvalAPA)
    plot(PCA25)

    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.
    ## Don't know how to automatically pick scale for object of type data.frame. Defaulting to continuous.

    dev.off()

    ## quartz_off_screen 
    ##                 2

Total Gene Counts Per Sample
----------------------------

    ## stats
    counts <- countData
    dim( counts )

    ## [1] 22485    44

    colSums( counts ) / 1e06  # in millions of reads

    ## 143A-CA3-1  143A-DG-1 143B-CA1-1  143B-DG-1 143C-CA1-1 143D-CA1-3 
    ##   1.600535   2.662497   0.874614   1.019113   1.086358   0.540738 
    ##  143D-DG-3 144A-CA1-2 144A-CA3-2  144A-DG-2 144B-CA1-1 144B-CA3-1 
    ##   0.500935   1.504457   0.228021   1.575766   1.275137   0.506698 
    ## 144C-CA1-2 144C-CA3-2  144C-DG-2 144D-CA3-2  144D-DG-2 145A-CA1-2 
    ##   1.613785   0.647568   1.083336   1.209306   2.254320   2.375356 
    ## 145A-CA3-2  145A-DG-2 145B-CA1-1  145B-DG-1 146A-CA1-2 146A-CA3-2 
    ##   0.179967   0.690882   1.034066   0.720798   0.878270   1.511881 
    ##  146A-DG-2 146B-CA1-2 146B-CA3-2  146B-DG-2 146C-CA1-4  146C-DG-4 
    ##   0.591933   0.506014   1.056001   0.055549   0.662938   0.237419 
    ## 146D-CA1-3 146D-CA3-3  146D-DG-3 147C-CA1-3 147C-CA3-3  147C-DG-3 
    ##   0.194359   1.467877   0.043490   1.506436   3.020727   2.118624 
    ## 147D-CA3-1  147D-DG-1 148A-CA1-3 148A-CA3-3  148A-DG-3 148B-CA1-4 
    ##   2.377445   5.618550   2.590815   1.353197   1.917857   0.185637 
    ## 148B-CA3-4  148B-DG-4 
    ##   1.724144   0.398258

    table( rowSums( counts ) )[ 1:30 ] # Number of genes with low counts

    ## 
    ##    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
    ## 4392  419  329  252  212  196  152  143  124  103   93   78   90   97   79 
    ##   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29 
    ##   60   71   83   67   65   61   70   55   56   62   56   50   44   37   40

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

![](../figures/02_rnaseq/totalRNAseqcounts-1.png)
