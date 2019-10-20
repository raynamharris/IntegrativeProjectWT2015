    library(tidyverse)
    library(forcats)
    library(cowplot) ## for some easy to use themes
    library(DESeq2) ## for gene expression analysis

    library(BiocParallel)
    register(MulticoreParam(6))

    ## load functions 
    source("figureoptions.R")
    source("functions_RNAseq.R")

    ## set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/02f_trainedvyoked/', cache = F)

    a.colData <- read.csv("../data/02a_colData.csv", header = T)
    a.countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)

    a.colData <- a.colData %>%
      mutate(combinedgroups = fct_collapse(Treatment,
                                           trained = c("conflict", "trained"),
                                           yoked = c("shocked", "yoked")))
    a.colData$combinedgroups <- factor(a.colData$combinedgroups, levels = c("yoked", "trained"))

    returndds2 <- function(mytissue){
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
                                  design = ~ combinedgroups)

      dds <- dds[ rowSums(counts(dds)) > 1, ]  # Pre-filtering genes with 0 counts
      dds <- DESeq(dds, parallel = TRUE)
      return(dds)
    }

    DGdds <- returndds2("DG") 

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

    CA3dds <- returndds2("CA3") 

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

    CA1dds <- returndds2("CA1") 

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

    returnvsds2 <- function(mydds, vsdfilename){
      dds <- mydds
      vsd <- vst(dds, blind=FALSE) ## variance stabilized
      print(head(assay(vsd),3))
      myvsd <- assay(vsd)
      myvsd <- as.data.frame(myvsd)
      return(myvsd)
    }

    DGvsd <- returnvsds2(DGdds)

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

    CA3vsd <- returnvsds2(CA3dds)

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

    CA1vsd <- returnvsds2(CA1dds)

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

    print("DG")

    ## [1] "DG"

    res_summary_subfield(DGdds, c("combinedgroups", "trained", "yoked"))

    ## [1] "combinedgroups" "trained"        "yoked"         
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
    ## log2 fold change (MLE): combinedgroups trained vs yoked 
    ## Wald test p-value: combinedgroups trained vs yoked 
    ## DataFrame with 5 rows and 6 columns
    ##                baseMean   log2FoldChange             lfcSE
    ##               <numeric>        <numeric>         <numeric>
    ## Smad7  171.392871064045 2.52915045290283 0.300055947411442
    ## Sgk1   341.089572273562 1.86437254406297  0.25434038097016
    ## Fzd5   26.8401177227407 3.21028367726372  0.45125074708102
    ## Acan   50.8597490321187 1.97329538674358  0.27993362211107
    ## Errfi1  196.30327794802  1.6626205039142 0.252431976397565
    ##                    stat               pvalue                 padj
    ##               <numeric>            <numeric>            <numeric>
    ## Smad7  8.42892958703734  3.4884140023949e-17 3.51701899721454e-13
    ## Sgk1   7.33022627768143 2.29764457026629e-13 1.15824262787124e-09
    ## Fzd5   7.11419027675832  1.1257175900809e-12 3.78316158106521e-09
    ## Acan    7.0491546240938 1.80008035791982e-12 4.53710254213691e-09
    ## Errfi1 6.58641004060308 4.50588351541571e-11 9.08566352048423e-08

    res_summary_subfield(CA3dds, c("combinedgroups", "trained", "yoked"))

    ## [1] "combinedgroups" "trained"        "yoked"         
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
    ## log2 fold change (MLE): combinedgroups trained vs yoked 
    ## Wald test p-value: combinedgroups trained vs yoked 
    ## DataFrame with 5 rows and 6 columns
    ##                  baseMean   log2FoldChange             lfcSE
    ##                 <numeric>        <numeric>         <numeric>
    ## AW011738 9.17276939908517 2.75526103220329 0.744391234682857
    ## Ccl4     20.1742899726882 2.17694386734611 0.608574059817134
    ## Cldn11   224.152887308636 1.14258100212196 0.299522313314391
    ## Gm9830    5.2365659464215  3.8500083682701  1.07167002084399
    ## Hbb-bt   13.3724700406441 3.05189602672508 0.846081438997779
    ##                      stat               pvalue              padj
    ##                 <numeric>            <numeric>         <numeric>
    ## AW011738 3.70136146669855 0.000214445727305785 0.715204813790057
    ## Ccl4     3.57712234399252 0.000347397602326682 0.715204813790057
    ## Cldn11   3.81467740910062 0.000136360945687631 0.715204813790057
    ## Gm9830   3.59253155671748 0.000327480981566366 0.715204813790057
    ## Hbb-bt   3.60709488006283 0.000309644457522496 0.715204813790057

    res_summary_subfield(CA1dds, c("combinedgroups", "trained", "yoked"))

    ## [1] "combinedgroups" "trained"        "yoked"         
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
    ## log2 fold change (MLE): combinedgroups trained vs yoked 
    ## Wald test p-value: combinedgroups trained vs yoked 
    ## DataFrame with 5 rows and 6 columns
    ##                 baseMean    log2FoldChange             lfcSE
    ##                <numeric>         <numeric>         <numeric>
    ## Glcci1  27.9431361099374  1.98275767136008  0.37236542908503
    ## Igf2bp2  8908.7711010914 -3.85162178226654 0.796911109754504
    ## Pde6a   16883.7614383138 -3.38908219974529 0.719364141096424
    ## Stox2   4433.92096042267 -3.79057173334957 0.812056110081023
    ## Gnaz    4292.66609719128 -3.11102225995498 0.674161232170265
    ##                      stat               pvalue                padj
    ##                 <numeric>            <numeric>           <numeric>
    ## Glcci1   5.32476303246537 1.01084730905324e-07 0.00143873897497547
    ## Igf2bp2 -4.83318871467743 1.34363305890186e-06 0.00956196466367509
    ## Pde6a   -4.71121926452964 2.46239128592003e-06  0.0108291598219097
    ## Stox2   -4.66786923501058  3.0433948772317e-06  0.0108291598219097
    ## Gnaz    -4.61465612601298 3.93746694525317e-06  0.0112083934063577

    listofDEGstrainedvyoked <- function(mydds, myitssue){
      res <- results(mydds, contrast = c("combinedgroups", "trained", "yoked"), independentFiltering = T)
      
      print(paste(myitssue, "trained vs yoked", sep = " "))
      
      data <- data.frame(gene = row.names(res),
                         lfc = res$log2FoldChange,
                         padj = res$padj,
                         tissue = myitssue,
                         comparison = paste("trained", "yoked", sep = "-"))
      data <- data %>% dplyr::filter(padj < 0.1) %>% droplevels()
      print(head(data))
      return(data)
      
    }

    DGDEGs <- listofDEGstrainedvyoked(DGdds, "DG")

    ## [1] "DG trained vs yoked"
    ##            gene       lfc         padj tissue    comparison
    ## 1 1190002N15Rik 1.6390074 2.451662e-04     DG trained-yoked
    ## 2 A830010M20Rik 1.5256191 7.889227e-07     DG trained-yoked
    ## 3         Abhd2 0.8633871 1.533267e-02     DG trained-yoked
    ## 4          Acan 1.9732954 4.537103e-09     DG trained-yoked
    ## 5       Adamts1 1.8771209 1.877775e-03     DG trained-yoked
    ## 6         Adrb1 0.9763798 3.110123e-02     DG trained-yoked

    CA1DEGs <- listofDEGstrainedvyoked(CA1dds, "CA1")

    ## [1] "CA1 trained vs yoked"
    ##     gene       lfc        padj tissue    comparison
    ## 1  Ahdc1 -1.607378 0.065020388    CA1 trained-yoked
    ## 2   Bmt2 -2.478179 0.056625662    CA1 trained-yoked
    ## 3  Ctcfl -2.612072 0.036325436    CA1 trained-yoked
    ## 4 Fn3krp -1.427451 0.048013887    CA1 trained-yoked
    ## 5 Glcci1  1.982758 0.001438739    CA1 trained-yoked
    ## 6   Gnaz -3.111022 0.011208393    CA1 trained-yoked
