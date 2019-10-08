Subfield analysis
-----------------

This script is used to identify treatement differences within each
subfield, generate volcano plots, venn diagrams, and tables for
subsequent GO analyses. The final mutlipanel figures for the manuscript
have been inserted just below the subheadings.

    library(tidyverse)
    library(cowplot) ## for some easy to use themes
    library(DESeq2) ## for gene expression analysis
    library(UpSetR)
    #devtools::install_github("clauswilke/ggtextures")
    library(ggtextures)
    library(magick)
    library(ggtext) # for markdown in plots

    library(BiocParallel)
    register(MulticoreParam(6))

    ## load functions 
    source("figureoptions.R")
    source("functions_RNAseq.R")

    ## set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/02c_rnaseqSubfield/', cache = T)

Get varience stabilized gene expression for each tissue
-------------------------------------------------------

    a.colData <- read.csv("../data/02a_colData.csv", header = T)
    a.countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)

    returndds <- function(mytissue){
      print(mytissue)
      colData <- a.colData %>% 
        filter(Punch %in% c(mytissue))  %>% 
      droplevels()
      
      savecols <- as.character(colData$RNAseqID) 
      savecols <- as.vector(savecols) 
      countData <- a.countData %>% dplyr::select(one_of(savecols)) 

      ## create DESeq object using the factors Punch and APA
      dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ APA2)

      dds <- dds[ rowSums(counts(dds)) > 1, ]  # Pre-filtering genes with 0 counts
      dds <- DESeq(dds, parallel = TRUE)
      return(dds)
    }

    DGdds <- returndds("DG") 

    ## [1] "DG"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    CA1dds <- returndds("CA1") 

    ## [1] "CA1"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    CA3dds <- returndds("CA3") 

    ## [1] "CA3"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    returnvsds <- function(mydds, vsdfilename){
      dds <- mydds
      vsd <- vst(dds, blind=FALSE) ## variance stabilized
      print(head(assay(vsd),3))
      return(write.csv(assay(vsd), file = vsdfilename, row.names = T))
    }

    returnvsds(DGdds, "../data/02c_DGvsd.csv")

    ##               143A-DG-1 143B-DG-1 143D-DG-3 144A-DG-2 144C-DG-2 144D-DG-2
    ## 0610007P14Rik  6.061813  6.093253  6.546473  6.191509  6.081483  6.330372
    ## 0610009B22Rik  5.427332  5.303459  4.916098  5.155434  5.469507  5.309580
    ## 0610009L18Rik  5.028120  5.230032  5.047701  4.596165  5.136024  5.134864
    ##               145A-DG-2 145B-DG-1 146A-DG-2 146B-DG-2 146C-DG-4 146D-DG-3
    ## 0610007P14Rik  6.255339  5.796661  5.981589  5.666827  6.642897  7.198395
    ## 0610009B22Rik  5.599669  5.436754  5.323551  5.666827  5.727726  4.596165
    ## 0610009L18Rik  4.869394  4.596165  5.323551  4.596165  4.596165  4.596165
    ##               147-DG-4 147C-DG-3 147D-DG-1 148-DG-2 148A-DG-3 148B-DG-4
    ## 0610007P14Rik 4.596165  5.976081  6.169871 5.962466  6.186452  6.033484
    ## 0610009B22Rik 7.598426  5.566137  5.548592 5.431813  5.226147  5.122685
    ## 0610009L18Rik 4.596165  4.869872  5.361845 4.596165  5.136778  4.969498

    returnvsds(CA1dds, "../data/02c_CA1vsd.csv")

    ##               143B-CA1-1 143C-CA1-1 143D-CA1-3 144A-CA1-2 144B-CA1-1
    ## 0610007P14Rik   6.880698   6.419536   6.666046   6.684314   6.722455
    ## 0610009B22Rik   6.271651   6.106656   5.319977   6.174401   5.957844
    ## 0610009L18Rik   5.982821   5.578039   5.319977   5.792716   5.491744
    ##               144C-CA1-2 145A-CA1-2 145B-CA1-1 146A-CA1-2 146B-CA1-2
    ## 0610007P14Rik   6.477328   6.733923   6.693372   6.430347   6.267586
    ## 0610009B22Rik   5.896135   6.076698   6.084280   6.127811   5.969877
    ## 0610009L18Rik   5.531571   5.900652   5.664972   5.913479   5.319977
    ##               146C-CA1-4 146D-CA1-3 147-CA1-4 147C-CA1-3 148-CA1-2
    ## 0610007P14Rik   6.733555   6.817016  6.305772   6.569382  6.682223
    ## 0610009B22Rik   5.946313   6.126776  6.305772   5.954570  5.812180
    ## 0610009L18Rik   6.028620   5.319977  5.319977   5.539188  6.040488
    ##               148A-CA1-3 148B-CA1-4
    ## 0610007P14Rik   6.642458   6.481732
    ## 0610009B22Rik   6.210959   5.319977
    ## 0610009L18Rik   5.789298   5.319977

    returnvsds(CA3dds, "../data/02c_CA3vsd.csv")

    ##               143A-CA3-1 144A-CA3-2 144B-CA3-1 144C-CA3-2 144D-CA3-2
    ## 0610007P14Rik   7.229103   7.692534   7.389598   7.189574   7.086379
    ## 0610009B22Rik   6.710399   7.171219   6.779068   7.111954   6.649230
    ## 0610009L18Rik   6.350172   6.906144   6.554700   6.519607   6.557669
    ##               145A-CA3-2 146A-CA3-2 146B-CA3-2 146D-CA3-3 147-CA3-4
    ## 0610007P14Rik   7.013722   7.407009   6.976367   7.296329  8.544319
    ## 0610009B22Rik   7.097826   6.532640   7.096058   6.727478  7.370690
    ## 0610009L18Rik   6.099026   6.532640   6.318145   6.447628  6.099026
    ##               147C-CA3-3 147D-CA3-1 148-CA3-2 148A-CA3-3 148B-CA3-4
    ## 0610007P14Rik   7.325639   7.050746  7.202076   7.113575   7.414902
    ## 0610009B22Rik   6.636480   6.790515  7.086099   6.776859   6.915261
    ## 0610009L18Rik   6.425657   6.425417  6.099026   6.467472   6.506630

Volcano plots
-------------

    DGconsyokcons <-  plot.cons.yokcons(DGdds, "DG", "DG: standard yoked v. trained") 

    ## [1] "DG"
    ## 
    ## out of 17033 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 69, 0.41%
    ## LFC < 0 (down)     : 2, 0.012%
    ## outliers [1]       : 27, 0.16%
    ## low counts [2]     : 4285, 25%
    ## (mean count < 3)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/02c_rnaseqSubfield/volcanos-1.png)

    CA1consyokcons <-  plot.cons.yokcons(CA1dds, "CA1", "CA1: standard yoked v. trained")  

    ## [1] "CA1"
    ## 
    ## out of 16896 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 183, 1.1%
    ## LFC < 0 (down)     : 105, 0.62%
    ## outliers [1]       : 31, 0.18%
    ## low counts [2]     : 5556, 33%
    ## (mean count < 6)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/02c_rnaseqSubfield/volcanos-2.png)

    CA1yoked <-  plot.yokconf.yokcons(CA1dds, "CA1", "CA1: standard v. conflict yoked")  

    ## [1] "CA1"
    ## 
    ## out of 16896 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 203, 1.2%
    ## LFC < 0 (down)     : 136, 0.8%
    ## outliers [1]       : 31, 0.18%
    ## low counts [2]     : 4905, 29%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/02c_rnaseqSubfield/volcanos-3.png)

    volcanos <- plot_grid(DGconsyokcons  ,CA1consyokcons  , CA1yoked , nrow = 1,
                          labels = c("(d)", "(e)", "(f)"), label_size = 8) 
    volcanos

![](../figures/02c_rnaseqSubfield/volcanos-4.png)

    pdf(file="../figures/02c_rnaseqSubfield/volcanos.pdf", width=6.65, height=2)
    plot(volcanos)    
    dev.off()

    ## quartz_off_screen 
    ##                 2

candidate gnees
===============

    betterPlotCounts <- function(mygene, mydds, mysubfield){
      df <- plotCounts(mydds, mygene, intgroup = "APA2",  transform = F, replaced = F, returnData = T)
      names(df) <- c("count", "treatment")
      df$treatment <- factor(df$treatment, levels = c("home.cage","standard.yoked","standard.trained",
                                                      "conflict.yoked", "conflict.trained"))
      
      #df <- df %>% filter(treatment %in% c("home.cage","standard.trained"))
      
      ggplot(df, aes(x = treatment, y = count)) +
        geom_boxplot(aes(fill = treatment)) + 
        geom_point() +
        labs(subtitle = paste(mysubfield, " *", mygene, "*",  sep = "")) +
        theme_classic() +
        theme(plot.subtitle  = element_markdown(),
              legend.position = "none", 
              panel.grid.major  = element_blank(),  # remove major gridlines
              panel.grid.minor  = element_blank()) + # remove minor gridlines) +
        scale_fill_manual(values = fivegroups) 

    }

    betterPlotCounts("Prkcz", DGdds, "DG")

![](../figures/02c_rnaseqSubfield/DG_boxplots-1.png)

    betterPlotCounts("Dusp16", DGdds, "DG")

![](../figures/02c_rnaseqSubfield/DG_boxplots-2.png)

    betterPlotCounts("Thbs1", DGdds, "DG")

![](../figures/02c_rnaseqSubfield/DG_boxplots-3.png)

    betterPlotCounts("Slc16a1", DGdds, "DG")

![](../figures/02c_rnaseqSubfield/DG_boxplots-4.png)

    betterPlotCounts("Hes5", DGdds, "DG")

![](../figures/02c_rnaseqSubfield/DG_boxplots-5.png)

    betterPlotCounts("Rtl1", DGdds, "DG")

![](../figures/02c_rnaseqSubfield/DG_boxplots-6.png)

    betterPlotCounts("Nlrp3", DGdds, "DG")

![](../figures/02c_rnaseqSubfield/DG_boxplots-7.png)

    betterPlotCounts("Kcnc2", DGdds, "DG")

![](../figures/02c_rnaseqSubfield/DG_boxplots-8.png)

    betterPlotCounts("1110008F13Rik", DGdds, "DG")

![](../figures/02c_rnaseqSubfield/DG_boxplots-9.png)

    a <- betterPlotCounts("Irs1", DGdds, "DG") + labs(title = "wald = -6, p = 1.577782e-06")
    b <- betterPlotCounts("Dab2ip", DGdds, "DG")  + labs(title = "wald =  5, p = 3.227358e-05")
    plot_grid(a,b)

![](../figures/02c_rnaseqSubfield/DG_boxplots-10.png)

    betterPlotCounts("Prkcz", CA1dds, "CA1")

![](../figures/02c_rnaseqSubfield/CA1_boxplots-1.png)

    betterPlotCounts("Grm1", CA1dds, "CA1")

![](../figures/02c_rnaseqSubfield/CA1_boxplots-2.png)

    betterPlotCounts("Grin2b", CA1dds, "CA1")

![](../figures/02c_rnaseqSubfield/CA1_boxplots-3.png)

    betterPlotCounts("Foxj3", CA1dds, "CA1")

![](../figures/02c_rnaseqSubfield/CA1_boxplots-4.png)

    betterPlotCounts("Grik3", CA1dds, "CA1")

![](../figures/02c_rnaseqSubfield/CA1_boxplots-5.png)

    betterPlotCounts("0610010K14Rik", CA1dds, "CA1")

![](../figures/02c_rnaseqSubfield/CA1_boxplots-6.png)

    a <- betterPlotCounts("Prkcz", CA1dds, "CA1")
    b <- betterPlotCounts("Prkcz", CA3dds, "CA3")
    c <- betterPlotCounts("Prkcz", DGdds, "DG")

    d <- betterPlotCounts("Prkci", CA1dds, "CA1")
    e <- betterPlotCounts("Prkci", CA3dds, "CA3")
    f <- betterPlotCounts("Prkci", DGdds, "DG")

    plot_grid(a,b,c,d,e,f)

![](../figures/02c_rnaseqSubfield/pkrcs-1.png)

genes that are correlated with number of entrances
--------------------------------------------------

Requires anlaysis of `04_integration.Rmd` first.

    plotCounts(DGdds, "Acan", intgroup = "APA2", normalized = TRUE, main="Acan in DG")

![](../figures/02c_rnaseqSubfield/DGcorrelations-1.png)

    plotCounts(DGdds, "Amigo2", intgroup = "APA2", normalized = TRUE, main="Amigo2 in DG")

![](../figures/02c_rnaseqSubfield/DGcorrelations-2.png)

    plotCounts(DGdds, "Armcx5", intgroup = "APA2", normalized = TRUE, main="Armcx5 in DG")

![](../figures/02c_rnaseqSubfield/DGcorrelations-3.png)

    plotCounts(DGdds, "Ptgs2", intgroup = "APA2", normalized = TRUE, main="Ptgs2 in DG")

![](../figures/02c_rnaseqSubfield/DGcorrelations-4.png)

    plotCounts(DGdds, "Rgs2", intgroup = "APA2", normalized = TRUE, main="Rgs2 in DG")

![](../figures/02c_rnaseqSubfield/DGcorrelations-5.png)

    plotCounts(DGdds, "Syt4", intgroup = "APA2", normalized = TRUE, main="Syt4 in DG")

![](../figures/02c_rnaseqSubfield/DGcorrelations-6.png)

Upset plots
-----------

What genes overlap within cetain comparisons?

    a.colData <- read.csv("../data/02a_colData.csv", header = T)
    a.countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)

    eachsubfield <- levels(a.colData$Punch)

    listofDEGs <- function(group1, group2){
      res <- results(dds, contrast = c("APA2", group1, group2), independentFiltering = T)
      
      print(paste(group1,group2, sep = " vs "))
      print(summary(res))
      
      data <- data.frame(gene = row.names(res),
                         lfc = res$log2FoldChange,
                         padj = res$padj,
                         tissue = i,
                         comparison = paste(group1, group2, sep = "-"))
      data <- data %>% dplyr::filter(padj < 0.1) %>% droplevels()
      return(data)
    }

    for(i in eachsubfield){
      
      colData <- a.colData %>% 
        dplyr::filter(Punch == i)  %>%
        droplevels()
      print(i)
      
      savecols <- as.character(colData$RNAseqID) 
      savecols <- as.vector(savecols) 
      countData <- a.countData %>% dplyr::select(one_of(savecols)) 

    ## create DESeq object using the factors Punch and APA
    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ APA2)

    dds # view the DESeq object - note numnber of genes
    dds <- dds[ rowSums(counts(dds)) > 1, ]  # Pre-filtering genes with 0 counts
    dds <- DESeq(dds, parallel = TRUE) # Differential expression analysis

    A <- listofDEGs("standard.trained","standard.yoked")
    B <- listofDEGs("conflict.trained","standard.trained")
    C <- listofDEGs("conflict.trained","conflict.yoked")
    D <- listofDEGs("conflict.yoked","standard.yoked")

    E <- listofDEGs("home.cage","standard.yoked")
    G <- listofDEGs("home.cage","standard.trained")
    H <- listofDEGs("home.cage","conflict.yoked")
    I <- listofDEGs("home.cage","conflict.trained")

    all <- rbind(A,B,C,D,E,G,H,I)

    write.csv(all, file = paste("../data/02c_",i,"forupset.csv", sep = ""), row.names = F)
    }

    ## [1] "CA1"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    ## [1] "standard.trained vs standard.yoked"
    ## 
    ## out of 16896 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 183, 1.1%
    ## LFC < 0 (down)     : 105, 0.62%
    ## outliers [1]       : 31, 0.18%
    ## low counts [2]     : 5556, 33%
    ## (mean count < 6)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "conflict.trained vs standard.trained"
    ## 
    ## out of 16896 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 31, 0.18%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "conflict.trained vs conflict.yoked"
    ## 
    ## out of 16896 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1, 0.0059%
    ## LFC < 0 (down)     : 1, 0.0059%
    ## outliers [1]       : 31, 0.18%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "conflict.yoked vs standard.yoked"
    ## 
    ## out of 16896 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 203, 1.2%
    ## LFC < 0 (down)     : 136, 0.8%
    ## outliers [1]       : 31, 0.18%
    ## low counts [2]     : 4905, 29%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "home.cage vs standard.yoked"
    ## 
    ## out of 16896 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 441, 2.6%
    ## LFC < 0 (down)     : 658, 3.9%
    ## outliers [1]       : 31, 0.18%
    ## low counts [2]     : 4905, 29%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "home.cage vs standard.trained"
    ## 
    ## out of 16896 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 209, 1.2%
    ## LFC < 0 (down)     : 420, 2.5%
    ## outliers [1]       : 31, 0.18%
    ## low counts [2]     : 4580, 27%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "home.cage vs conflict.yoked"
    ## 
    ## out of 16896 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 257, 1.5%
    ## LFC < 0 (down)     : 561, 3.3%
    ## outliers [1]       : 31, 0.18%
    ## low counts [2]     : 4905, 29%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "home.cage vs conflict.trained"
    ## 
    ## out of 16896 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 236, 1.4%
    ## LFC < 0 (down)     : 492, 2.9%
    ## outliers [1]       : 31, 0.18%
    ## low counts [2]     : 4905, 29%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "CA3"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    ## [1] "standard.trained vs standard.yoked"
    ## 
    ## out of 16583 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1, 0.006%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 9, 0.054%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "conflict.trained vs standard.trained"
    ## 
    ## out of 16583 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 9, 0.054%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "conflict.trained vs conflict.yoked"
    ## 
    ## out of 16583 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 9, 0.054%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "conflict.yoked vs standard.yoked"
    ## 
    ## out of 16583 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1, 0.006%
    ## LFC < 0 (down)     : 1, 0.006%
    ## outliers [1]       : 9, 0.054%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "home.cage vs standard.yoked"
    ## 
    ## out of 16583 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 11, 0.066%
    ## LFC < 0 (down)     : 38, 0.23%
    ## outliers [1]       : 9, 0.054%
    ## low counts [2]     : 4177, 25%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "home.cage vs standard.trained"
    ## 
    ## out of 16583 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 49, 0.3%
    ## LFC < 0 (down)     : 260, 1.6%
    ## outliers [1]       : 9, 0.054%
    ## low counts [2]     : 4497, 27%
    ## (mean count < 5)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "home.cage vs conflict.yoked"
    ## 
    ## out of 16583 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 2, 0.012%
    ## outliers [1]       : 9, 0.054%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "home.cage vs conflict.trained"
    ## 
    ## out of 16583 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 26, 0.16%
    ## LFC < 0 (down)     : 170, 1%
    ## outliers [1]       : 9, 0.054%
    ## low counts [2]     : 4497, 27%
    ## (mean count < 5)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "DG"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    ## [1] "standard.trained vs standard.yoked"
    ## 
    ## out of 17033 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 69, 0.41%
    ## LFC < 0 (down)     : 2, 0.012%
    ## outliers [1]       : 27, 0.16%
    ## low counts [2]     : 4285, 25%
    ## (mean count < 3)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "conflict.trained vs standard.trained"
    ## 
    ## out of 17033 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 27, 0.16%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "conflict.trained vs conflict.yoked"
    ## 
    ## out of 17033 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 1, 0.0059%
    ## outliers [1]       : 27, 0.16%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "conflict.yoked vs standard.yoked"
    ## 
    ## out of 17033 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 3, 0.018%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 27, 0.16%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "home.cage vs standard.yoked"
    ## 
    ## out of 17033 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 908, 5.3%
    ## LFC < 0 (down)     : 352, 2.1%
    ## outliers [1]       : 27, 0.16%
    ## low counts [2]     : 5274, 31%
    ## (mean count < 5)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "home.cage vs standard.trained"
    ## 
    ## out of 17033 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 834, 4.9%
    ## LFC < 0 (down)     : 411, 2.4%
    ## outliers [1]       : 27, 0.16%
    ## low counts [2]     : 4944, 29%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "home.cage vs conflict.yoked"
    ## 
    ## out of 17033 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 875, 5.1%
    ## LFC < 0 (down)     : 480, 2.8%
    ## outliers [1]       : 27, 0.16%
    ## low counts [2]     : 4944, 29%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL
    ## [1] "home.cage vs conflict.trained"
    ## 
    ## out of 17033 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1055, 6.2%
    ## LFC < 0 (down)     : 550, 3.2%
    ## outliers [1]       : 27, 0.16%
    ## low counts [2]     : 4285, 25%
    ## (mean count < 3)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

Standard deviation stuff
------------------------

    res <- results(DGdds, contrast =c("APA2", "home.cage", "standard.trained"))
    res

    ## log2 fold change (MLE): APA2 home.cage vs standard.trained 
    ## Wald test p-value: APA2 home.cage vs standard.trained 
    ## DataFrame with 17033 rows and 6 columns
    ##                       baseMean     log2FoldChange             lfcSE
    ##                      <numeric>          <numeric>         <numeric>
    ## 0610007P14Rik  33.106847648943 -0.824688805172559 0.706201212652936
    ## 0610009B22Rik 15.4308333541882    2.4559198020875   1.0053171058589
    ## 0610009L18Rik 2.01213758864991  -1.16221292431944  2.27283316314352
    ## 0610009O20Rik 42.3702664233309  0.711085181407761 0.592961039759206
    ## 0610010F05Rik 64.7038994019559   1.07957936173565 0.629057203565766
    ## ...                        ...                ...               ...
    ## Zxdc          43.1486819277236 -0.766976885312159 0.594530233960452
    ## Zyg11b        321.833700366766   1.08433885348256 0.510532304484217
    ## Zyx           95.5520024625261   1.61913942562895 0.693930838538127
    ## Zzef1         150.425199831691 -0.500013592342612 0.491422612638784
    ## Zzz3          141.469097182289    1.5276616995255 0.575200802445279
    ##                             stat              pvalue               padj
    ##                        <numeric>           <numeric>          <numeric>
    ## 0610007P14Rik  -1.16778163276513   0.242894858709327  0.610408782238628
    ## 0610009B22Rik   2.44293048210819  0.0145685432626862  0.123429619198798
    ## 0610009L18Rik -0.511349861998668    0.60910609613768                 NA
    ## 0610009O20Rik   1.19921062890831   0.230446055372204  0.595424576296302
    ## 0610010F05Rik   1.71618631122279  0.0861279478114974  0.351162213907659
    ## ...                          ...                 ...                ...
    ## Zxdc           -1.29005531006044   0.197031454883037  0.549395614740788
    ## Zyg11b          2.12393778798787  0.0336753460502013  0.204522462998817
    ## Zyx             2.33328645407937   0.019633115916954  0.148566276154516
    ## Zzef1          -1.01748185672144   0.308924258994894  0.681214700547789
    ## Zzz3            2.65587546649995 0.00791028201436956 0.0840650411077759

    summary(res)

    ## 
    ## out of 17033 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 834, 4.9%
    ## LFC < 0 (down)     : 411, 2.4%
    ## outliers [1]       : 27, 0.16%
    ## low counts [2]     : 4944, 29%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    resSig <- subset(res, padj < 0.1)
    resSig <- as.data.frame(resSig)
    resSig$gene <- row.names(resSig)
    resSig <- resSig %>% arrange(stat)
    head(resSig)

    ##    baseMean log2FoldChange     lfcSE      stat       pvalue         padj
    ## 1  84.21645      -5.272785 0.8205097 -6.426232 1.308060e-10 1.577782e-06
    ## 2  39.01042      -7.231125 1.3104090 -5.518220 3.424502e-08 7.645806e-05
    ## 3 152.69849      -3.779588 0.6943040 -5.443708 5.218263e-08 7.645806e-05
    ## 4 168.14676      -3.709184 0.7158185 -5.181738 2.198273e-07 1.473087e-04
    ## 5  38.29160      -5.778719 1.1493714 -5.027721 4.963422e-07 2.394752e-04
    ## 6  39.27817      -3.905596 0.7979443 -4.894572 9.851984e-07 3.304778e-04
    ##     gene
    ## 1   Irs1
    ## 2   Acan
    ## 3  Frmd6
    ## 4  Ptgs2
    ## 5 Adgrg1
    ## 6   Mios

    head(assays(DGdds)[["mu"]])

    ##               143A-DG-1 143B-DG-1 143D-DG-3  144A-DG-2  144C-DG-2
    ## 0610007P14Rik 125.51220 39.128224 39.202876  74.495857  60.363965
    ## 0610009B22Rik  28.51032 11.155148  5.196671  16.921864  20.521916
    ## 0610009L18Rik  10.30368  3.794863  3.465592   6.115592   2.805252
    ## 0610009O20Rik 198.49054 61.377565 35.248745 117.811039  70.217834
    ## 0610010F05Rik 198.97150 85.100062 41.396385 118.096504 110.742907
    ## 0610010K14Rik  61.80884 41.932264 12.555485  36.685695  19.763325
    ##               144D-DG-2 145A-DG-2 145B-DG-1 146A-DG-2 146B-DG-2  146C-DG-4
    ## 0610007P14Rik 178.41307 35.011671 27.475921 29.107447 1.6642736 13.2096875
    ## 0610009B22Rik  23.65015  7.952962  7.833168  6.611807 0.4744713  4.4908928
    ## 0610009L18Rik  15.77198  2.874215  2.664761  2.389519 0.1614101  0.6138845
    ## 0610009O20Rik 160.41774 55.369003 43.099455 46.031802 2.6106236 15.3660489
    ## 0610010F05Rik 188.39577 55.503167 59.757443 46.143341 3.6196325 24.2343126
    ## 0610010K14Rik  57.14026 17.241596 29.444924 14.334044 1.7835402  4.3248873
    ##               146D-DG-3   147-DG-4 147C-DG-3 147D-DG-1   148-DG-2
    ## 0610007P14Rik 2.7413665 1.18006267 118.43938 449.80613  34.988778
    ## 0610009B22Rik 0.3633912 3.89858326  40.26580  59.62558 115.592730
    ## 0610009L18Rik 0.2423409 0.04340046   5.50415  39.76353   1.286821
    ## 0610009O20Rik 2.4648633 3.98005226 137.77354 404.43720 118.008281
    ## 0610010F05Rik 2.8947535 8.10373203 217.28728 474.97402 240.275108
    ## 0610010K14Rik 0.8779760 0.38435899  38.77737 144.05917  11.396218
    ##               148A-DG-3 148B-DG-4
    ## 0610007P14Rik  97.52312 14.248047
    ## 0610009B22Rik  22.15255  4.062006
    ## 0610009L18Rik   8.00597  1.381851
    ## 0610009O20Rik 154.22737 22.349863
    ## 0610010F05Rik 154.60108 30.988110
    ## 0610010K14Rik  48.02554 15.269103

    head(assays(DGdds)[["cooks"]])

    ##                  143A-DG-1  143B-DG-1    143D-DG-3   144A-DG-2  144C-DG-2
    ## 0610007P14Rik 0.0047563777 0.04134201 0.0054250301 0.002074481 0.03493602
    ## 0610009B22Rik 0.0350000133 0.01279839 0.2592337018 0.146187763 0.09452981
    ## 0610009L18Rik 0.0006873141 0.08492622 0.0096689407 0.039273461 0.12798056
    ## 0610009O20Rik 0.0068402553 0.07395699 0.0093284481 0.122813378 0.14066383
    ## 0610010F05Rik 0.0305325964 0.15504542 0.1911136419 0.016748112 0.02292438
    ## 0610010K14Rik 0.0087830510 0.10107064 0.0005578724 0.043109661 0.02444902
    ##                 144D-DG-2   145A-DG-2  145B-DG-1    146A-DG-2   146B-DG-2
    ## 0610007P14Rik 0.019881864 0.009108214 0.03809356 0.0127633548 0.004315805
    ## 0610009B22Rik 0.001341341 0.194068436 0.06220278 0.0023478043 0.010501794
    ## 0610009L18Rik 0.002299207 0.013285044 0.05950413 0.0657668325 0.006331534
    ## 0610009O20Rik 0.010003784 0.003757888 0.00919852 0.0018009626 0.003139963
    ## 0610010F05Rik 0.078280582 0.002013498 0.09952414 0.0002934111 0.013390475
    ## 0610010K14Rik 0.003282354 0.026030270 0.01280175 0.0575034420 0.220773428
    ##                 146C-DG-4   146D-DG-3    147-DG-4   147C-DG-3   147D-DG-1
    ## 0610007P14Rik 0.195497379 0.069518483 0.047851744 0.107318778 0.085969880
    ## 0610009B22Rik 0.043310414 0.005762166 1.278932710 0.002906652 1.647032927
    ## 0610009L18Rik 0.026093660 0.007348606 0.009552561 0.026837584 0.037065486
    ## 0610009O20Rik 0.118262917 0.002212013 0.372001978 0.011024654 0.015260052
    ## 0610010F05Rik 0.088207959 0.006035332 1.343856544 0.035886249 0.004695359
    ## 0610010K14Rik 0.004265768 0.013422750 0.015760008 0.060694479 0.020048654
    ##                 148-DG-2   148A-DG-3   148B-DG-4
    ## 0610007P14Rik  1.2908053 0.001752103 0.006749729
    ## 0610009B22Rik 11.7200413 0.084232385 0.088854014
    ## 0610009L18Rik  0.3106720 0.005781759 0.003018016
    ## 0610009O20Rik  5.7323780 0.306253729 0.046181651
    ## 0610010F05Rik  4.6360407 0.011684718 0.052072898
    ## 0610010K14Rik  0.3476747 0.091077411 0.035367588

    sizeFactors(DGdds)

    ##  143A-DG-1  143B-DG-1  143D-DG-3  144A-DG-2  144C-DG-2  144D-DG-2 
    ## 4.11991222 1.68610385 0.83727353 2.44531116 1.75103389 3.81044847 
    ##  145A-DG-2  145B-DG-1  146A-DG-2  146B-DG-2  146C-DG-4  146D-DG-3 
    ## 1.14925089 1.18398567 0.95544594 0.07171647 0.38318574 0.05854860 
    ##   147-DG-4  147C-DG-3  147D-DG-1   148-DG-2  148A-DG-3  148B-DG-4 
    ## 0.06062865 3.43568179 9.60671264 1.79763552 3.20117640 0.61397337

    head(coef(DGdds))

    ##               Intercept APA2_conflict.yoked_vs_conflict.trained
    ## 0610007P14Rik  4.929070                             -0.39262591
    ## 0610009B22Rik  2.790799                             -0.06485422
    ## 0610009L18Rik  1.322474                             -0.15211976
    ## 0610009O20Rik  5.590313                             -0.40436674
    ## 0610010F05Rik  5.593804                              0.06359049
    ## 0610010K14Rik  3.907128                              0.72916779
    ##               APA2_home.cage_vs_conflict.trained
    ## 0610007P14Rik                         -0.6463504
    ## 0610009B22Rik                          3.2160078
    ## 0610009L18Rik                         -1.8047637
    ## 0610009O20Rik                          0.4463309
    ## 0610010F05Rik                          1.4686384
    ## 0610010K14Rik                         -1.2427450
    ##               APA2_standard.trained_vs_conflict.trained
    ## 0610007P14Rik                                 0.1783384
    ## 0610009B22Rik                                 0.7600880
    ## 0610009L18Rik                                -0.6425508
    ## 0610009O20Rik                                -0.2647543
    ## 0610010F05Rik                                 0.3890591
    ## 0610010K14Rik                                -0.4105809
    ##               APA2_standard.yoked_vs_conflict.trained
    ## 0610007P14Rik                            0.6200464503
    ## 0610009B22Rik                           -0.1569818073
    ## 0610009L18Rik                            0.7268569386
    ## 0610009O20Rik                           -0.1945837943
    ## 0610010F05Rik                            0.0338575601
    ## 0610010K14Rik                           -0.0006527513
