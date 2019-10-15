New approach. Behavior-centric analysis figures first.

    library(tidyverse) ## for respahing data

    ## ── Attaching packages ───────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.2.1     ✔ purrr   0.3.2
    ## ✔ tibble  2.1.3     ✔ dplyr   0.8.1
    ## ✔ tidyr   0.8.3     ✔ stringr 1.4.0
    ## ✔ readr   1.3.1     ✔ forcats 0.4.0

    ## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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

    knitr::opts_chunk$set(echo = TRUE, fig.path = '../figures/fig1/')

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


    top <- plot_grid(h,i,j,k, ncol = 1)
    top

![](../figures/fig1/behav-1.png)

    a.colData <- read.csv("../data/02a_colData.csv", header = T)
    a.countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)

    trained <- c("standard.yoked", "standard.trained")
    conflict <- c("conflict.yoked", "conflict.trained")
    trainedconflict <- c("standard.trained", "conflict.trained")
    yokedyoked <- c("standard.yoked", "conflict.yoked")


    # trained
    DGdds <- returndds("DG", trained) 

    ## [1] "DG"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    CA1dds <- returndds("CA1", trained) 

    ## [1] "CA1"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    CA3dds <- returndds("CA3", trained) 

    ## [1] "CA3"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    e <-  plot.cons.yokcons(DGdds, "DG", "DG") 

    ## [1] "DG"
    ## 
    ## out of 16461 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 106, 0.64%
    ## LFC < 0 (down)     : 4, 0.024%
    ## outliers [1]       : 25, 0.15%
    ## low counts [2]     : 2873, 17%
    ## (mean count < 2)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/fig1/DESeq2-1.png)

    f <-  plot.cons.yokcons(CA3dds, "CA3", "CA3") 

    ## [1] "CA3"
    ## 
    ## out of 15699 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 21, 0.13%
    ## LFC < 0 (down)     : 3, 0.019%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 9740, 62%
    ## (mean count < 88)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/fig1/DESeq2-2.png)

    g <-  plot.cons.yokcons(CA1dds, "CA1", "CA1") 

    ## [1] "CA1"
    ## 
    ## out of 15918 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 305, 1.9%
    ## LFC < 0 (down)     : 196, 1.2%
    ## outliers [1]       : 13, 0.082%
    ## low counts [2]     : 4628, 29%
    ## (mean count < 6)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/fig1/DESeq2-3.png)

    efg <- plot_grid(e + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                     f + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                     g, nrow = 3)
    efg

![](../figures/fig1/DESeq2-4.png)

    # conflict
    DGdds <- returndds("DG", conflict) 

    ## [1] "DG"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    CA1dds <- returndds("CA1", conflict) 

    ## [1] "CA1"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    CA3dds <- returndds("CA3", conflict) 

    ## [1] "CA3"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    l <-  plot.conf.yokconf(DGdds, "DG", "DG") 

    ## [1] "DG"
    ## 
    ## out of 16252 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 30, 0.18%
    ## LFC < 0 (down)     : 13, 0.08%
    ## outliers [1]       : 28, 0.17%
    ## low counts [2]     : 7538, 46%
    ## (mean count < 21)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/fig1/DESeq2-5.png)

    m <-  plot.conf.yokconf(CA3dds, "CA3", "CA3") 

    ## [1] "CA3"
    ## 
    ## out of 15884 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 24, 0.15%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/fig1/DESeq2-6.png)

    n <-  plot.conf.yokconf(CA1dds, "CA1", "CA1") 

    ## [1] "CA1"
    ## 
    ## out of 16170 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1, 0.0062%
    ## LFC < 0 (down)     : 6, 0.037%
    ## outliers [1]       : 33, 0.2%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/fig1/DESeq2-7.png)

    lmn <- plot_grid(l + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                     m + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                     n, nrow = 3)
    lmn

![](../figures/fig1/DESeq2-8.png)

    # trained v conflict
    DGdds <- returndds("DG", trainedconflict) 

    ## [1] "DG"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    CA1dds <- returndds("CA1", trainedconflict) 

    ## [1] "CA1"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    CA3dds <- returndds("CA3", trainedconflict) 

    ## [1] "CA3"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    o <-  plot.conf.cons(DGdds, "DG") 

    ## [1] "DG"
    ## 
    ## out of 16556 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 24, 0.14%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/fig1/DESeq2-9.png)

    p <-  plot.conf.cons(CA3dds, "CA3") 

    ## [1] "CA3"
    ## 
    ## out of 15721 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 16, 0.1%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/fig1/DESeq2-10.png)

    q <-  plot.conf.cons(CA1dds, "CA1") 

    ## [1] "CA1"
    ## 
    ## out of 16595 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 223, 1.3%
    ## LFC < 0 (down)     : 252, 1.5%
    ## outliers [1]       : 45, 0.27%
    ## low counts [2]     : 12184, 73%
    ## (mean count < 136)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/fig1/DESeq2-11.png)

    opq <- plot_grid(o + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                     p + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                     q, nrow = 3)
    opq

![](../figures/fig1/DESeq2-12.png)

    # yoked v yoked

    DGdds <- returndds("DG", yokedyoked) 

    ## [1] "DG"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    CA1dds <- returndds("CA1", yokedyoked) 

    ## [1] "CA1"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    CA3dds <- returndds("CA3", yokedyoked) 

    ## [1] "CA3"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    r <-  plot.yokconf.yokcons(DGdds, "DG") 

    ## [1] "DG"
    ## 
    ## out of 16135 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1, 0.0062%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 18, 0.11%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/fig1/DESeq2-13.png)

    s <-  plot.yokconf.yokcons(CA3dds, "CA3") 

    ## [1] "CA3"
    ## 
    ## out of 15856 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1, 0.0063%
    ## LFC < 0 (down)     : 3, 0.019%
    ## outliers [1]       : 5, 0.032%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/fig1/DESeq2-14.png)

    t <-  plot.yokconf.yokcons(CA1dds, "CA1") 

    ## [1] "CA1"
    ## 
    ## out of 15085 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 20, 0.13%
    ## LFC < 0 (down)     : 4, 0.027%
    ## outliers [1]       : 12, 0.08%
    ## low counts [2]     : 6722, 45%
    ## (mean count < 14)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/fig1/DESeq2-15.png)

    rst <- plot_grid(r + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                     s + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                     t, nrow = 3)
    rst

![](../figures/fig1/DESeq2-16.png)

    opqrst <- plot_grid(opq,rst, ncol = 2)

    schematicTrained <- ggdraw() +  draw_image("../figures/figure_fig1a.png")
    schematicConflict <- ggdraw() +  draw_image("../figures/figure_fig2a.png")


    left <- plot_grid(schematicTrained,a,b,c, ncol = 1)
    right <- plot_grid(e,f,g, ncol = 1)

    anotherleft <- plot_grid(schematicConflict,h,i,j, ncol = 1)
    anotherright <- plot_grid(l,m,n, ncol = 1)

    farleft <- plot_grid(left,right, ncol = 2)
    farright <- plot_grid(anotherleft,anotherright, ncol = 2)



    plot_grid(farleft,farright, opqrst, nrow = 1)

![](../figures/fig1/fig1-1.png)
