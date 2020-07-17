Subfield analysis
-----------------

This script is used to identify treatement differences within each
subfield, generate volcano plots, venn diagrams, and tables for
subsequent GO analyses. The final mutlipanel figures for the manuscript
have been inserted just below the subheadings.

    library(tidyverse)
    library(cowplot) ## for some easy to use themes
    library(DESeq2) ## for gene expression analysis
    library(png)
    library(grid)
    library(scales)
    library(apaTables) #  for ANOVA tables

    library(BiocParallel)
    register(MulticoreParam(6))

    ## load functions 
    source("figureoptions.R")
    source("functions_RNAseq.R")

    ## set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/03_rnaseqSubfield/', cache = T)

Wrangle data
------------

    # prep col data,
    outliers <- c("146D-DG-3", "145A-CA3-2", "146B-DG-2", "146D-CA1-3", "148B-CA1-4")

    colData <- read.csv("../data/00_colData.csv", header = T) %>%
      filter(!RNAseqID %in% outliers)
    colData$training <- factor(colData$training, levels = levelstraining)
    colData$treatment <- factor(colData$treatment, levels = levelstreatment)

    # remove outliers
    savecols <- as.character(colData$RNAseqID) #select the rowsname 
    savecols <- as.vector(savecols) # make it a vector

    countData <- read.csv("../data/00_countData.csv", 
                          header = T, check.names = F, row.names = 1) %>%
      dplyr::select(one_of(savecols)) # select just the columns 
    head(countData)

    ##               143A-CA3-1 143A-DG-1 143B-CA1-1 143B-DG-1 143C-CA1-1 143D-CA1-3
    ## 0610007P14Rik         85       112         60        48         38         28
    ## 0610009B22Rik         24        34         21        10         19          0
    ## 0610009L18Rik          4         9         10         8          2          0
    ## 0610009O20Rik         85       185         44        72         76         25
    ## 0610010F05Rik        142       155         54       117         57         39
    ## 0610010K14Rik         24        74         14        21         23         21
    ##               143D-DG-3 144A-CA1-2 144A-CA3-2 144A-DG-2 144B-CA1-1 144B-CA3-1
    ## 0610007P14Rik        43         80         21        80         72         34
    ## 0610009B22Rik         1         30          9         9         14          9
    ## 0610009L18Rik         2          9          5         0          1          4
    ## 0610009O20Rik        38         95         20        82         49         24
    ## 0610010F05Rik        61        119         17       138         75         25
    ## 0610010K14Rik        12         32         10        20         41         11
    ##               144C-CA1-2 144C-CA3-2 144C-DG-2 144D-CA3-2 144D-DG-2 145A-CA1-2
    ## 0610007P14Rik         63         28        49         43       150        133
    ## 0610009B22Rik         15         24        16         13        23         36
    ## 0610009L18Rik          2          4         6          9        13         21
    ## 0610009O20Rik         89         48        85         46       151         96
    ## 0610010F05Rik         95         59        97        111       140        179
    ## 0610010K14Rik         38         14        24         17        62         34
    ##               145A-DG-2 145B-CA1-1 145B-DG-1 146A-CA1-2 146A-CA3-2 146A-DG-2
    ## 0610007P14Rik        41         51        21         29         87        23
    ## 0610009B22Rik        14         15        10         15          9         6
    ## 0610009L18Rik         1          3         0          8          9         6
    ## 0610009O20Rik        52        124        46         68        102        44
    ## 0610010F05Rik        52         30        41         69         83        45
    ## 0610010K14Rik        24         27        24         33         12        23
    ##               146B-CA1-2 146B-CA3-2 146C-CA1-4 146C-DG-4 146D-CA3-3 147C-CA1-3
    ## 0610007P14Rik         13         33         38        22         87         69
    ## 0610009B22Rik          6         43          7         6         23         17
    ## 0610009L18Rik          0          2          9         0          7          2
    ## 0610009O20Rik         16         46         31        10         83         58
    ## 0610010F05Rik         53        110         41        32        148        149
    ## 0610010K14Rik         11         11          7         5         10         46
    ##               147C-CA3-3 147C-DG-3 147D-CA3-1 147D-DG-1 148A-CA1-3 148A-CA3-3
    ## 0610007P14Rik        164        82         79       305        135         55
    ## 0610009B22Rik         30        39         41       105         59         24
    ## 0610009L18Rik         11         3          9        67         16          7
    ## 0610009O20Rik        191       145         88       377        162         65
    ## 0610010F05Rik        328       185        222       446        198        149
    ## 0610010K14Rik         31        27          0       173         60         40
    ##               148A-DG-3 148B-CA3-4 148B-DG-4
    ## 0610007P14Rik       104        122        16
    ## 0610009B22Rik        15         45         2
    ## 0610009L18Rik        11         11         1
    ## 0610009O20Rik       226         70        18
    ## 0610010F05Rik       176        177        23
    ## 0610010K14Rik        17         39        10

Get varience stabilized gene expression for each tissue
-------------------------------------------------------

    # DEGs with looking at all four treatments individually
    DGdds <- returnddstreatment("DG") 

    ## [1] "DG"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    CA3dds <- returnddstreatment("CA3") 

    ## [1] "CA3"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    CA1dds <- returnddstreatment("CA1") 

    ## [1] "CA1"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    # DEGs with looking at all grouped trained and yoked
    DGdds2 <- returnddstraining("DG") 

    ## [1] "DG"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    ## -- replacing outliers and refitting for 58 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    CA3dds2 <- returnddstraining("CA3") 

    ## [1] "CA3"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    CA1dds2 <- returnddstraining("CA1") 

    ## [1] "CA1"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    ## -- replacing outliers and refitting for 82 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    savevsds(DGdds2, "../data/03_DG_vsdtraining.csv")

    ##               143A-DG-1 143B-DG-1 143D-DG-3 144A-DG-2 144C-DG-2 144D-DG-2
    ## 0610007P14Rik  7.153718  7.167228  7.569987  7.271891  7.178889  7.395072
    ## 0610009B22Rik  6.607383  6.495129  6.172178  6.378164  6.648238  6.509850
    ## 0610009L18Rik  6.268598  6.433471  6.282645  5.904182  6.362953  6.360970
    ##               145A-DG-2 145B-DG-1 146A-DG-2 146C-DG-4 147C-DG-3 147D-DG-1
    ## 0610007P14Rik  7.326240  6.920197  7.075780  7.654996  7.065253  7.242216
    ## 0610009B22Rik  6.756751  6.612735  6.514596  6.858188  6.715943  6.707290
    ## 0610009L18Rik  6.135124  5.904182  6.514596  5.904182  6.132069  6.548655
    ##               148A-DG-3 148B-DG-4
    ## 0610007P14Rik  7.247632  7.131677
    ## 0610009B22Rik  6.430112  6.349605
    ## 0610009L18Rik  6.355221  6.219767

    savevsds(CA3dds2, "../data/03_CA3_vsdtraining.csv")

    ##               143A-CA3-1 144A-CA3-2 144B-CA3-1 144C-CA3-2 144D-CA3-2 146A-CA3-2
    ## 0610007P14Rik   7.159328   7.688310   7.347158   7.097481   6.973766   7.346091
    ## 0610009B22Rik   6.518149   7.068934   6.600198   7.003315   6.435071   6.284311
    ## 0610009L18Rik   6.064367   6.746786   6.320031   6.272773   6.320996   6.284311
    ##               146B-CA3-2 146D-CA3-3 147C-CA3-3 147D-CA3-1 148A-CA3-3 148B-CA3-4
    ## 0610007P14Rik   6.842126   7.230928   7.244900   6.929800   7.016913   7.377410
    ## 0610009B22Rik   6.988761   6.533834   6.410424   6.609981   6.600412   6.768893
    ## 0610009L18Rik   6.021975   6.184369   6.150573   6.155501   6.212184   6.259681

    savevsds(CA1dds2, "../data/03_CA1_vsdtraining.csv")

    ##               143B-CA1-1 143C-CA1-1 143D-CA1-3 144A-CA1-2 144B-CA1-1 144C-CA1-2
    ## 0610007P14Rik   7.626135   7.205231   7.426464   7.430562   7.477774   7.258673
    ## 0610009B22Rik   7.063836   6.916634   6.195554   6.965957   6.775834   6.723175
    ## 0610009L18Rik   6.799392   6.431682   6.195554   6.621009   6.351612   6.389144
    ##               145A-CA1-2 145B-CA1-1 146A-CA1-2 146B-CA1-2 146C-CA1-4 147C-CA1-3
    ## 0610007P14Rik   7.492491   7.451376   7.197799   7.077403   7.494980   7.335370
    ## 0610009B22Rik   6.886631   6.891557   6.923209   6.799601   6.768535   6.772299
    ## 0610009L18Rik   6.725450   6.509225   6.729555   6.195554   6.844058   6.394538
    ##               148A-CA1-3
    ## 0610007P14Rik   7.411631
    ## 0610009B22Rik   7.012540
    ## 0610009L18Rik   6.625123

Results to compare with volcano plots
-------------------------------------

    print("DG")

    ## [1] "DG"

    res_summary_subfield(DGdds2, c("training", "trained", "yoked"))

    ## [1] "training" "trained"  "yoked"   
    ## [1] 214
    ## 
    ## out of 17006 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 177, 1%
    ## LFC < 0 (down)     : 37, 0.22%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 6929, 41%
    ## (mean count < 13)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    res_summary_subfield(DGdds, c("treatment", "conflict.trained", "standard.trained"))

    ## [1] "treatment"        "conflict.trained" "standard.trained"
    ## [1] 0
    ## 
    ## out of 17011 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 20, 0.12%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    res_summary_subfield(DGdds, c("treatment", "conflict.yoked", "standard.yoked"))

    ## [1] "treatment"      "conflict.yoked" "standard.yoked"
    ## [1] 3
    ## 
    ## out of 17011 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 3, 0.018%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 20, 0.12%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    print("CA3")

    ## [1] "CA3"

    res_summary_subfield(CA3dds2, c("training", "trained", "yoked"))

    ## [1] "training" "trained"  "yoked"   
    ## [1] 0
    ## 
    ## out of 16497 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 27, 0.16%
    ## low counts [2]     : 5, 0.03%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    res_summary_subfield(CA3dds, c("treatment", "conflict.trained", "standard.trained"))

    ## [1] "treatment"        "conflict.trained" "standard.trained"
    ## [1] 0
    ## 
    ## out of 16502 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 11, 0.067%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    res_summary_subfield(CA3dds, c("treatment", "conflict.yoked", "standard.yoked"))

    ## [1] "treatment"      "conflict.yoked" "standard.yoked"
    ## [1] 2
    ## 
    ## out of 16502 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1, 0.0061%
    ## LFC < 0 (down)     : 1, 0.0061%
    ## outliers [1]       : 11, 0.067%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    print("CA1")

    ## [1] "CA1"

    res_summary_subfield(CA1dds2, c("training", "trained", "yoked"))

    ## [1] "training" "trained"  "yoked"   
    ## [1] 16
    ## 
    ## out of 16846 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1, 0.0059%
    ## LFC < 0 (down)     : 15, 0.089%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 2619, 16%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    res_summary_subfield(CA1dds, c("treatment", "conflict.trained", "standard.trained"))

    ## [1] "treatment"        "conflict.trained" "standard.trained"
    ## [1] 0
    ## 
    ## out of 16852 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 32, 0.19%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    res_summary_subfield(CA1dds, c("treatment", "conflict.yoked", "standard.yoked"))

    ## [1] "treatment"      "conflict.yoked" "standard.yoked"
    ## [1] 917
    ## 
    ## out of 16852 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 545, 3.2%
    ## LFC < 0 (down)     : 372, 2.2%
    ## outliers [1]       : 32, 0.19%
    ## low counts [2]     : 4892, 29%
    ## (mean count < 5)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

Volcano plots
-------------

    # create data frame for making volcanos plots
    DGa <-  calculateDEGs(DGdds, "DG", "treatment", "standard.trained", "standard.yoked") 
    DGb <-  calculateDEGs(DGdds, "DG", "treatment", "conflict.trained", "conflict.yoked") 
    DGc <-  calculateDEGs(DGdds, "DG", "treatment", "conflict.trained", "standard.trained") 
    DGd <-  calculateDEGs(DGdds, "DG", "treatment", "conflict.yoked", "standard.yoked") 
    DGe <- calculateDEGs(DGdds2, "DG", "training", "trained", "yoked") 

    CA3a <-  calculateDEGs(CA3dds, "CA3", "treatment", "standard.trained", "standard.yoked") 
    CA3b <-  calculateDEGs(CA3dds, "CA3", "treatment", "conflict.trained", "conflict.yoked") 
    CA3c <-  calculateDEGs(CA3dds, "CA3", "treatment", "conflict.trained", "standard.trained") 
    CA3d <-  calculateDEGs(CA3dds, "CA3", "treatment", "conflict.yoked", "standard.yoked") 
    CA3e <- calculateDEGs(CA3dds2, "CA3", "training", "trained", "yoked") 

    CA1a <-  calculateDEGs(CA1dds, "CA1", "treatment", "standard.trained", "standard.yoked") 
    CA1b <-  calculateDEGs(CA1dds, "CA1", "treatment", "conflict.trained", "conflict.yoked") 
    CA1c <-  calculateDEGs(CA1dds, "CA1", "treatment", "conflict.trained", "standard.trained") 
    CA1d <-  calculateDEGs(CA1dds, "CA1", "treatment", "conflict.yoked", "standard.yoked") 
    CA1e <- calculateDEGs(CA1dds2, "CA1", "training", "trained", "yoked") 

    # save df with DEGs

    allDEG <- rbind(DGa, DGb, DGc, DGd, DGe,
                   CA3a, CA3b, CA3c, CA3d, CA3e, 
                   CA1a, CA1b, CA1c, CA1d, CA1e) %>% 
      dplyr::filter(direction != "NS") %>%
        dplyr::mutate(lfc = round(lfc, 2),
                      padj = scientific(padj, digits = 3),
                      logpadj = round(logpadj, 2)) %>%
      arrange(tissue, comparison, gene)

pca analysis and bar plots functions
------------------------------------

    a <- plotPCs(DGdds, "DG") 
    b <- plot.volcano(DGa, "\ns. trained vs s. yoked") + labs(y = "-log10(p-value)")
    c <- plot.volcano(DGb, "\nc. trained vs c. yoked")  
    d <- plot.volcano(DGc, "\ns. trained vs c. trained")  
    e <- plot.volcano(DGd, "\ns. yoked vs. c. yoked")  
    f <- plot.volcano(DGe, "\nyoked vs. trained")  

    g <- plotPCs(CA3dds, "CA3")
    h <- plot.volcano(CA3a, " ") + labs(y = "-log10(p-value)")
    i <- plot.volcano(CA3b, " ")  
    j <- plot.volcano(CA3c, " ")  
    k <- plot.volcano(CA3d, " ")  
    l <- plot.volcano(CA3e, " ")  


    m <- plotPCs(CA1dds, "CA1")
    n <- plot.volcano(CA1a, " ") + labs(y = "-log10(p-value)")
    o <- plot.volcano(CA1b, " ")  
    p <- plot.volcano(CA1c, " ")   
    q <- plot.volcano(CA1d, " ")  
    r <- plot.volcano(CA1e, " ")  

    legend <- get_legend(a + theme(legend.position = "bottom", 
                                   legend.title = element_blank()) +
                           guides(color = guide_legend(nrow = 2)))

    mainplot <- plot_grid(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r, 
                          nrow = 3, rel_widths = c(1,1,0.8,0.8,0.8,0.8),
                          labels = c("(a)", "(b)", "", "", "", "(c)", 
                                     "", "", "", "", "", "", 
                                     "", "", "", "", "", ""),
                          label_size = 8)

    ## Warning in MASS::cov.trob(data[, vars]): Probable convergence failure

    fig3 <- plot_grid(mainplot, legend, ncol = 1, rel_heights = c(1, 0.1))
    fig3

![](../figures/03_rnaseqSubfield/volcanos-1.png)

save files
----------

    pdf(file="../figures/03_rnaseqSubfield/volcanos.pdf", width=6.69, height=6)
    plot(fig3)    
    dev.off()

    ## quartz_off_screen 
    ##                 2

    pdf(file="../figures/fig-3.pdf", width=6.69, height=6)
    plot(fig3)    
    dev.off()

    ## quartz_off_screen 
    ##                 2

    suppltable4 <- allDEG %>% filter(tissue == "DG"  & comparison == "yoked vs. trained")
    head(suppltable4)

    ##   tissue          gene   lfc     padj logpadj        comparison direction
    ## 1     DG 1190002N15Rik  1.78 3.60e-04    3.44 yoked vs. trained   trained
    ## 2     DG 2410002F23Rik -0.51 9.26e-02    1.03 yoked vs. trained     yoked
    ## 3     DG A830010M20Rik  1.47 1.68e-07    6.78 yoked vs. trained   trained
    ## 4     DG         Abhd2  0.68 1.33e-02    1.88 yoked vs. trained   trained
    ## 5     DG          Acan  1.96 5.25e-10    9.28 yoked vs. trained   trained
    ## 6     DG       Adamts1  1.75 1.28e-02    1.89 yoked vs. trained   trained

    suppltable5 <- allDEG %>% filter(tissue != "DG"  & comparison != "yoked vs. trained")
    head(suppltable5)

    ##   tissue          gene   lfc     padj logpadj
    ## 1    CA1          Chd8 -0.77 8.19e-02    1.09
    ## 2    CA1 D430019H16Rik -0.64 8.19e-02    1.09
    ## 3    CA1         Fgfr1 -0.66 8.19e-02    1.09
    ## 4    CA1       Gm20390  2.84 1.06e-03    2.97
    ## 5    CA1        Itga10 -2.17 8.19e-02    1.09
    ## 6    CA1       Ppp1r10 -0.86 6.80e-02    1.17
    ##                            comparison        direction
    ## 1 conflict.yoked vs. conflict.trained   conflict.yoked
    ## 2 conflict.yoked vs. conflict.trained   conflict.yoked
    ## 3 conflict.yoked vs. conflict.trained   conflict.yoked
    ## 4 conflict.yoked vs. conflict.trained conflict.trained
    ## 5 conflict.yoked vs. conflict.trained   conflict.yoked
    ## 6 conflict.yoked vs. conflict.trained   conflict.yoked

    write_csv(suppltable4, "../data/suppltable-4.csv")
    write_csv(suppltable5, "../data/suppltable-5.csv")

    citation("DESeq2") 

    ## 
    ##   Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change
    ##   and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550
    ##   (2014)
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Article{,
    ##     title = {Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2},
    ##     author = {Michael I. Love and Wolfgang Huber and Simon Anders},
    ##     year = {2014},
    ##     journal = {Genome Biology},
    ##     doi = {10.1186/s13059-014-0550-8},
    ##     volume = {15},
    ##     issue = {12},
    ##     pages = {550},
    ##   }

    citation("png")

    ## 
    ## To cite package 'png' in publications use:
    ## 
    ##   Simon Urbanek (2013). png: Read and write PNG images. R package
    ##   version 0.1-7. https://CRAN.R-project.org/package=png
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {png: Read and write PNG images},
    ##     author = {Simon Urbanek},
    ##     year = {2013},
    ##     note = {R package version 0.1-7},
    ##     url = {https://CRAN.R-project.org/package=png},
    ##   }
    ## 
    ## ATTENTION: This citation information has been auto-generated from the
    ## package DESCRIPTION file and may need manual editing, see
    ## 'help("citation")'.

    citation("grid")

    ## 
    ## The 'grid' package is part of R.  To cite R in publications use:
    ## 
    ##   R Core Team (2019). R: A language and environment for statistical
    ##   computing. R Foundation for Statistical Computing, Vienna, Austria.
    ##   URL https://www.R-project.org/.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {R: A Language and Environment for Statistical Computing},
    ##     author = {{R Core Team}},
    ##     organization = {R Foundation for Statistical Computing},
    ##     address = {Vienna, Austria},
    ##     year = {2019},
    ##     url = {https://www.R-project.org/},
    ##   }
    ## 
    ## We have invested a lot of time and effort in creating R, please cite it
    ## when using it for data analysis. See also 'citation("pkgname")' for
    ## citing R packages.

    citation("BiocParallel")

    ## 
    ## To cite package 'BiocParallel' in publications use:
    ## 
    ##   Martin Morgan, Valerie Obenchain, Michel Lang, Ryan Thompson and
    ##   Nitesh Turaga (2019). BiocParallel: Bioconductor facilities for
    ##   parallel evaluation. R package version 1.18.0.
    ##   https://github.com/Bioconductor/BiocParallel
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {BiocParallel: Bioconductor facilities for parallel evaluation},
    ##     author = {Martin Morgan and Valerie Obenchain and Michel Lang and Ryan Thompson and Nitesh Turaga},
    ##     year = {2019},
    ##     note = {R package version 1.18.0},
    ##     url = {https://github.com/Bioconductor/BiocParallel},
    ##   }
