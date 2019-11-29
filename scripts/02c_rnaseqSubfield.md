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

    a.colData$training <- factor(a.colData$training, levels = c("yoked", "trained"))
    a.colData$treatment <- factor(a.colData$treatment, levels = c("standard.yoked","standard.trained",
                                                      "conflict.yoked", "conflict.trained"))

    DGdds <- returnddstreatment("DG") 

    ## [1] "DG"

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

    CA3dds <- returnddstreatment("CA3") 

    ## [1] "CA3"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    DGdds2 <- returndds2("DG") 

    ## [1] "DG"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    ## -- replacing outliers and refitting for 54 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    CA1dds2 <- returndds2("CA1") 

    ## [1] "CA1"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    ## -- replacing outliers and refitting for 98 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    CA3dds2 <- returndds2("CA3") 

    ## [1] "CA3"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    ## -- replacing outliers and refitting for 53 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    returnvsds(DGdds2, "../data/02c_DGvsd.csv")

    ##               143A-DG-1 143B-DG-1 143D-DG-3 144A-DG-2 144C-DG-2 144D-DG-2
    ## 0610007P14Rik  6.497231  6.520310  6.927199  6.622006  6.516323  6.739843
    ## 0610009B22Rik  5.933401  5.821822  5.483341  5.697426  5.971474  5.830043
    ## 0610009L18Rik  5.582826  5.757564  5.597896  5.205323  5.677909  5.676574
    ##               145A-DG-2 145B-DG-1 146A-DG-2 146B-DG-2 146C-DG-4 146D-DG-3
    ## 0610007P14Rik  6.669012  6.257587  6.419110  6.134992  7.017216  7.479295
    ## 0610009B22Rik  6.084135  5.939702  5.838619  6.134992  6.195432  5.205323
    ## 0610009L18Rik  5.443572  5.205323  5.838619  5.205323  5.205323  5.205323
    ##               147C-DG-3 147D-DG-1 148A-DG-3 148B-DG-4
    ## 0610007P14Rik  6.418946  6.592025  6.603472  6.467815
    ## 0610009B22Rik  6.054862  6.038992  5.754020  5.664124
    ## 0610009L18Rik  5.444089  5.874552  5.675946  5.530425

    ## class: DESeqTransform 
    ## dim: 17011 16 
    ## metadata(36): version version ... version version
    ## assays(1): ''
    ## rownames(17011): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(23): baseMean baseVar ... replace dispFit
    ## colnames(16): 143A-DG-1 143B-DG-1 ... 148A-DG-3 148B-DG-4
    ## colData names(7): RNAseqID ID ... sizeFactor replaceable

    returnvsds(CA1dds2, "../data/02c_CA1vsd.csv")

    ##               143B-CA1-1 143C-CA1-1 143D-CA1-3 144A-CA1-2 144B-CA1-1
    ## 0610007P14Rik   7.415967   7.026659   7.234953   7.242002   7.278148
    ## 0610009B22Rik   6.904844   6.767590   6.124712   6.819405   6.644173
    ## 0610009L18Rik   6.666449   6.334833   6.124712   6.507772   6.264241
    ##               144C-CA1-2 145A-CA1-2 145B-CA1-1 146A-CA1-2 146B-CA1-2
    ## 0610007P14Rik   7.075661   7.291669   7.252495   7.032181   6.904449
    ## 0610009B22Rik   6.595087   6.743767   6.747191   6.782460   6.657871
    ## 0610009L18Rik   6.297128   6.599020   6.404817   6.606986   6.124712
    ##               146C-CA1-4 146D-CA1-3 147C-CA1-3 148A-CA1-3 148B-CA1-4
    ## 0610007P14Rik   7.290895   7.357361   7.149545   7.213656   7.068189
    ## 0610009B22Rik   6.636289   6.782490   6.641399   6.854020   6.124712
    ## 0610009L18Rik   6.703932   6.124712   6.302770   6.507437   6.124712

    ## class: DESeqTransform 
    ## dim: 16852 15 
    ## metadata(36): version version ... version version
    ## assays(1): ''
    ## rownames(16852): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(23): baseMean baseVar ... replace dispFit
    ## colnames(15): 143B-CA1-1 143C-CA1-1 ... 148A-CA1-3 148B-CA1-4
    ## colData names(7): RNAseqID ID ... sizeFactor replaceable

    returnvsds(CA3dds2, "../data/02c_CA3vsd.csv")

    ##               143A-CA3-1 144A-CA3-2 144B-CA3-1 144C-CA3-2 144D-CA3-2
    ## 0610007P14Rik   7.053774   7.561649   7.233758   6.996898   6.880131
    ## 0610009B22Rik   6.452876   6.973831   6.531009   6.908274   6.375915
    ## 0610009L18Rik   6.030451   6.670471   6.269295   6.224806   6.269568
    ##               145A-CA3-2 146A-CA3-2 146B-CA3-2 146D-CA3-3 147C-CA3-3
    ## 0610007P14Rik   6.784496   7.233784   6.757449   7.122555   7.137734
    ## 0610009B22Rik   6.879922   6.236353   6.895220   6.468228   6.353976
    ## 0610009L18Rik   5.734766   6.236353   5.991499   6.142415   6.111529
    ##               147D-CA3-1 148A-CA3-3 148B-CA3-4
    ## 0610007P14Rik   6.840862   6.923812   7.262408
    ## 0610009B22Rik   6.540744   6.532445   6.689017
    ## 0610009L18Rik   6.116202   6.169424   6.213015

    ## class: DESeqTransform 
    ## dim: 16502 13 
    ## metadata(36): version version ... version version
    ## assays(1): ''
    ## rownames(16502): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(23): baseMean baseVar ... replace dispFit
    ## colnames(13): 143A-CA3-1 144A-CA3-2 ... 148A-CA3-3 148B-CA3-4
    ## colData names(7): RNAseqID ID ... sizeFactor replaceable

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

    plot.volcano <- function(mydds, whichtissue, whichfactor, up, down, mycolors){

      # calculate DEG results
      res <- results(mydds, contrast =c(whichfactor, up, down),
                     independentFiltering = T, alpha = 0.1)
       
      # create dataframe with pvalues and lfc
      data <- data.frame(gene = row.names(res),
                         padj = res$padj, 
                         logpadj = -log10(res$padj),
                         lfc = res$log2FoldChange)
      data <- na.omit(data)
      data <- data %>%
        dplyr::mutate(direction = ifelse(data$lfc > 0 & data$padj < 0.1, 
                                         yes = up, no = ifelse(data$lfc < 0 & data$padj < 0.1, 
                                                     yes = down, no = "NS")))
      data$direction <- factor(data$direction, levels = c(down, "NS", up))
      
      #save df that only has the values for significant genes IF there are DEGs
      toprint <- data %>% filter(direction != "NS")
      myfilename = paste("../data/02c_DEGs", whichtissue, whichfactor, down, up, "csv" ,sep = ".")
      
      # for some reason i can only return the df or return the plot...
      #if (dim(toprint)[1] != 0) { return(write.csv(toprint, myfilename))}
      
      # plot volcanos
      volcano <- data %>%
        ggplot(aes(x = lfc, y = logpadj)) + 
        geom_point(aes(color = direction), size = 1, alpha = 0.75, na.rm = T) + 
         theme_ms() +
        scale_color_manual(values = mycolors,
                           name = " ",
                           drop = FALSE) +
        ylim(c(0,12.5)) +  
        xlim(c(-8,8)) +
        labs(y = NULL, x = "log fold change", subtitle = " ")  +
        theme(legend.position = "top",
              legend.spacing.x = unit(-0.1, 'cm'),
              legend.margin=margin(t=-0.25, r=0, b=0, l=0, unit="cm"),
              panel.grid = element_blank()) 
      return(volcano)
      

    }

    # usage: mydds, whichfactor, up, down, mycolors, mysubfield
    DGa <- plot.volcano(DGdds2, "DG", "training", "trained", "yoked", volcano7) 
    DGb <-  plot.volcano(DGdds, "DG", "treatment", "conflict.trained", "standard.trained",  volcano2) 
    DGc <-  plot.volcano(DGdds, "DG", "treatment", "conflict.yoked", "standard.yoked",  volcano6) 

    CA3a <- plot.volcano(CA3dds2, "CA3", "training", "trained", "yoked", volcano7) 
    CA3b <-  plot.volcano(CA3dds, "CA3","treatment", "conflict.trained", "standard.trained", volcano2) 
    CA3c <-  plot.volcano(CA3dds, "CA3","treatment", "conflict.yoked", "standard.yoked", volcano6) 

    CA1a <- plot.volcano(CA1dds2, "CA1","training", "trained", "yoked", volcano7) 
    CA1b <-  plot.volcano(CA1dds, "CA1","treatment", "conflict.trained", "standard.trained", volcano2) 
    CA1c <-  plot.volcano(CA1dds, "CA1","treatment", "conflict.yoked", "standard.yoked", volcano6) 

    volcanos <- plot_grid(DGa + labs(x= NULL, y = "DG \n -log10(p)") + theme(legend.position = "none"),
                          DGb + labs(x= NULL) + theme(legend.position = "none"), 
                          DGc + labs(x= NULL) + theme(legend.position = "none"),
                          
                          CA3a + theme(legend.position = "none") + labs(x= NULL,  y = "CA3 \n -log10(p)"),
                          CA3b + theme(legend.position = "none") +  labs(x= NULL ), 
                          CA3c + theme(legend.position = "none") +  labs(x= NULL ), 
                          
                          CA1a + theme(legend.position = "bottom") +  labs( y = "CA1 \n -log10(p)"), 
                          CA1b + theme(legend.position = "bottom"), 
                          CA1c + theme(legend.position = "bottom"), 
                          nrow = 3, 
                          rel_heights =  c(1,1,1.25) ,
                          rel_widths = c(1.2,1,1) ,
                          labels =  c("a", "d", "g",
                                      "b", "e", "h",
                                      "c", "f", "i"),
                          label_size = 8) 
    volcanos

![](../figures/02c_rnaseqSubfield/volcanos-1.png)

    comparisons <- png::readPNG("../figures/00_schematics/figure_DEGcomparisons.png")
    comparisons <- ggdraw() +  draw_image(comparisons, scale = 1)

    right <- plot_grid(comparisons,volcanos, nrow = 2, rel_heights = c(0.1,1))

    circuit <- png::readPNG("../figures/00_schematics/figure_hippocircuit.png")
    circuit <- ggdraw() +  draw_image(circuit, scale = 1)

    left <- plot_grid(circuit, right, nrow = 1, rel_widths = c(0.1,1))
    left

![](../figures/02c_rnaseqSubfield/volcanos-2.png)

    pdf(file="../figures/02c_rnaseqSubfield/volcanos.pdf", width=6.69, height=6)
    plot(left)    
    dev.off()

    ## quartz_off_screen 
    ##                 2

    pdf(file="../figures/figure_3.pdf", width=6.69, height=6)
    plot(left)    
    dev.off()

    ## quartz_off_screen 
    ##                 2

    citation("DESeq2") ## for gene expression analysis

    ## 
    ##   Love, M.I., Huber, W., Anders, S. Moderated estimation of fold
    ##   change and dispersion for RNA-seq data with DESeq2 Genome
    ##   Biology 15(12):550 (2014)
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
    ## ATTENTION: This citation information has been auto-generated from
    ## the package DESCRIPTION file and may need manual editing, see
    ## 'help("citation")'.

    citation("grid")

    ## 
    ## The 'grid' package is part of R.  To cite R in publications use:
    ## 
    ##   R Core Team (2019). R: A language and environment for
    ##   statistical computing. R Foundation for Statistical Computing,
    ##   Vienna, Austria. URL https://www.R-project.org/.
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
    ## We have invested a lot of time and effort in creating R, please
    ## cite it when using it for data analysis. See also
    ## 'citation("pkgname")' for citing R packages.

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
