New approach. Behavior-centric analysis figures first.

    library(tidyverse) ## for respahing data

    ## ── Attaching packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.2.1     ✔ purrr   0.3.2
    ## ✔ tibble  2.1.3     ✔ dplyr   0.8.1
    ## ✔ tidyr   0.8.3     ✔ stringr 1.4.0
    ## ✔ readr   1.3.1     ✔ forcats 0.4.0

    ## ── Conflicts ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

    library(cowplot) ## for some easy to use themes

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    library(DESeq2)

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter,
    ##     Find, get, grep, grepl, intersect, is.unsorted, lapply, Map,
    ##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
    ##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
    ##     setdiff, sort, table, tapply, union, unique, unsplit, which,
    ##     which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## Loading required package: BiocParallel

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following object is masked from 'package:purrr':
    ## 
    ##     simplify

    ## The following objects are masked from 'package:base':
    ## 
    ##     aperm, apply, rowsum

    library("png")
    library("grid")

    library(BiocParallel)
    register(MulticoreParam(6))

    source("functions_RNAseq.R")
    source("figureoptions.R")

    knitr::opts_chunk$set(echo = TRUE, fig.path = '../figures/fig1/', message = F, results = F)

    # all behavior data
    behav <- read.csv("../data/01a_behavior.csv") 

    # make mouse name
    behav$mouse <- sapply(strsplit(as.character(behav$ID),"15"), "[", 2)

    behav <-  behav %>% 
      select(mouse,APA2, TrainSessionCombo, TrainSessionComboNum, Time1stEntr, Path1stEntr, pTimeTarget, NumEntrances)

    # subset to standard or conflict paradigm only
    standard <- behav %>% filter(APA2 %in% c("standard-yoked", "standard-trained")) 
    standard$APA2 <- factor(standard$APA2, levels = c("standard-yoked", "standard-trained"))

    conflict <- behav %>% filter(APA2 %in% c("conflict-yoked", "conflict-trained")) 
    conflict$APA2 <- factor(conflict$APA2, levels = c("conflict-yoked", "conflict-trained"))


    # gather and summarize

    calculatemeandev <- function(mydf){
      
      mydf <- mydf %>% gather(behavior, measure, Time1stEntr:NumEntrances)
      mydf$behavior <- factor(mydf$behavior, levels = c("Path1stEntr", "Time1stEntr", "pTimeTarget", "NumEntrances"))

      meandev <- mydf %>%
        dplyr::group_by(APA2, TrainSessionComboNum, behavior) %>%
        dplyr::summarise(m = mean(measure), 
                       se = sd(measure)/sqrt(length(measure)))
      return(meandev)
    }

    standard.meandev <- calculatemeandev(standard)
    conflict.meandev <- calculatemeandev(conflict)


    plotmeansd <- function(mymeandev, mybehavior, myylab, mycolors){
      
      mymeandev %>% 
        filter(behavior == mybehavior) %>% 
        ggplot(aes(x=, TrainSessionComboNum, y=m, color=APA2)) + 
          geom_errorbar(aes(ymin=m-se, ymax=m+se, color=APA2), width=.1) +
          geom_point(size = 1.5) +
          geom_line() +
        theme_minimal(base_size = 8) + 
        theme(legend.position = "none") + 
        scale_color_manual(values = mycolors) +
        scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                           labels = c( "P", "T1", "T2", "T3",
                                       "Rt", "T4", "T5", "T6", "Rn")) +
        labs(x = NULL, y = myylab) 
      }
      
    # Time1stEntr, Path1stEntr, pTimeTarget, NumEntrances
    a <- plotmeansd(standard.meandev, "Time1stEntr", "Time 1st entr. (s)" , trainedcolors)
    b <- plotmeansd(standard.meandev, "Path1stEntr", "Path 1st entr. (m)" , trainedcolors)
    c <- plotmeansd(standard.meandev, "pTimeTarget", "Prop. time in zone" , trainedcolors)
    d <- plotmeansd(standard.meandev, "NumEntrances", "Num. of entr." , trainedcolors)

    h <- plotmeansd(conflict.meandev, "Time1stEntr", "Time 1st entr. (s)" , conflictcolors)
    i <- plotmeansd(conflict.meandev, "Path1stEntr", "Path 1st entr. (m)" , conflictcolors)
    j <- plotmeansd(conflict.meandev, "pTimeTarget", "Prop. time in zone" , conflictcolors)
    k <- plotmeansd(conflict.meandev, "NumEntrances", "Num. of entr." , conflictcolors)

    a.colData <- read.csv("../data/02a_colData.csv", header = T)
    a.countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)

    trained <- c("standard.yoked", "standard.trained")
    conflict <- c("conflict.yoked", "conflict.trained")
    trainedconflict <- c("standard.trained", "conflict.trained")
    yokedyoked <- c("standard.yoked", "conflict.yoked")


    # trained
    DGdds <- returndds("DG", trained) 
    CA1dds <- returndds("CA1", trained) 
    CA3dds <- returndds("CA3", trained) 

    e <-  plot.cons.yokcons(DGdds, "DG", "DG") 

![](../figures/fig1/DESeq2-1.png)

    f <-  plot.cons.yokcons(CA3dds, "CA3", "CA3") 

![](../figures/fig1/DESeq2-2.png)

    g <-  plot.cons.yokcons(CA1dds, "CA1", "CA1") 

![](../figures/fig1/DESeq2-3.png)

    efg <- plot_grid(e + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                     f + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                     g, nrow = 3)

    # conflict
    DGdds <- returndds("DG", conflict) 
    CA1dds <- returndds("CA1", conflict) 
    CA3dds <- returndds("CA3", conflict) 

    l <-  plot.conf.yokconf(DGdds, "DG", "DG") 

![](../figures/fig1/DESeq2-4.png)

    m <-  plot.conf.yokconf(CA3dds, "CA3", "CA3") 

![](../figures/fig1/DESeq2-5.png)

    n <-  plot.conf.yokconf(CA1dds, "CA1", "CA1") 

![](../figures/fig1/DESeq2-6.png)

    lmn <- plot_grid(l + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                     m + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                     n, nrow = 3)

    # trained v conflict
    DGdds <- returndds("DG", trainedconflict) 
    CA1dds <- returndds("CA1", trainedconflict) 
    CA3dds <- returndds("CA3", trainedconflict) 

    o <-  plot.conf.cons(DGdds, "DG") 

![](../figures/fig1/DESeq2-7.png)

    p <-  plot.conf.cons(CA3dds, "CA3") 

![](../figures/fig1/DESeq2-8.png)

    q <-  plot.conf.cons(CA1dds, "CA1") 

![](../figures/fig1/DESeq2-9.png)

    opq <- plot_grid(o + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                     p + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                     q, nrow = 3)

    # yoked v yoked

    DGdds <- returndds("DG", yokedyoked) 
    CA1dds <- returndds("CA1", yokedyoked) 
    CA3dds <- returndds("CA3", yokedyoked) 

    r <-  plot.yokconf.yokcons(DGdds, "DG") 

![](../figures/fig1/DESeq2-10.png)

    s <-  plot.yokconf.yokcons(CA3dds, "CA3") 

![](../figures/fig1/DESeq2-11.png)

    t <-  plot.yokconf.yokcons(CA1dds, "CA1") 

![](../figures/fig1/DESeq2-12.png)

    rst <- plot_grid(r + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                     s + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                     t, nrow = 3)

    schematicTrained <- ggdraw() +  draw_image("../figures/figure_fig1a.png")
    schematicConflict <- ggdraw() +  draw_image("../figures/figure_fig2a.png")

    a14 <- plot_grid(a,b,d,d, ncol = 1)
    b14 <- plot_grid(h + theme(axis.title.y = element_blank()),
                i + theme(axis.title.y = element_blank()),
                j + theme(axis.title.y = element_blank()),
                k + theme(axis.title.y = element_blank()), ncol = 1)
    c13 <- plot_grid(e + theme(axis.title.x = element_blank()),
                        f + theme(axis.title.x = element_blank()),
                     g, ncol = 1)
    d13 <- plot_grid(l + theme(axis.title = element_blank()),
                     m + theme( axis.title = element_blank()),
                     n+ theme( axis.title.y = element_blank()),
                     ncol = 1)

    opqrst <- plot_grid(opq,rst, ncol = 2, labels = c( "standard v\nconflict", "yoked v \nyoked"),label_size = 8) 


    fig1top <- plot_grid(a14, b14, c13, d13, opqrst, nrow = 1,
              labels = c("standard training", "conflict training", "standard\nyoked v trained", "conflict\nyoked v trained", NULL),
              label_size = 8,
              rel_widths = c(1,0.9,0.75,0.75,1))
    fig1top

![](../figures/fig1/fig1-1.png)
