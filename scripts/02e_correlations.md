    library(tidyverse) 

    ## ── Attaching packages ──────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.2.1     ✔ purrr   0.3.2
    ## ✔ tibble  2.1.3     ✔ dplyr   0.8.1
    ## ✔ tidyr   0.8.3     ✔ stringr 1.4.0
    ## ✔ readr   1.3.1     ✔ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

    library(corrplot)

    ## corrplot 0.84 loaded

    library(cowplot)

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    knitr::opts_chunk$set(fig.path = '../figures/02e_correlations/', cache = F)

For this analysis, I want to explor correlations between a behavioral
measure and gene expression.

    # import behavior data, create mouse id, select relvant samples
    behav <- read.csv("../data/01a_behavior.csv") 
    behav$mouse <- sapply(strsplit(as.character(behav$ID),"15"), "[", 2)
    behav <- behav %>% filter(APA2 %in% c("conflict-trained", "standard-trained"),
                                      TrainSession == "Retention") %>% 
                               select(mouse,Time1stEntr,pTimeTarget) 
    behav$mouse <- as.character(behav$mouse)
    head(behav)

    ##   mouse Time1stEntr pTimeTarget
    ## 1  140C      599.97      0.0021
    ## 2  141C       30.53      0.0909
    ## 3  142C      482.43      0.0445
    ## 4  143C      516.47      0.0039
    ## 5  144C      599.97      0.0000
    ## 6  145C       68.53      0.1336

    # import varance stabilized data and fix sample names
    DGvsd <- read.csv("../data/02c_DGvsd.csv", stringsAsFactors = F, check.names = F, row.names = 1) 
    names(DGvsd) <- sapply(strsplit(names(DGvsd),"-"), "[", 1)

    # import col Data 
    a.colData <- read.csv("../data/02a_colData.csv", header = T)
    a.colData$mouse <- sapply(strsplit(as.character(a.colData$RNAseqID),"-"), "[", 1)

    # subset data by subfield and training and create list of samples for substting
    DGcols <- a.colData %>% 
      filter(Punch == "DG", APA2 %in% c("conflict.trained", "standard.trained")) %>% 
      droplevels()
    DGtrained <- DGcols$mouse

    # keep only trained varance stabilized data
    DGvsd <- DGvsd %>% select(DGtrained)

    # transform and set rownames
    DGvsd <- as.data.frame(t(DGvsd))
    DGvsd$mouse <- row.names(DGvsd)
    DGvsdbehav <- left_join(behav, DGvsd) %>% drop_na()

    ## Joining, by = "mouse"

    row.names(DGvsdbehav) <- DGvsdbehav$mouse
    DGvsdbehav$mouse <- NULL
    #head(DGvsdbehav)

    # subset for genes that are significant
    DGsig <- read.csv("../data/02c_DGforupset.csv", stringsAsFactors = F) %>% select(gene) %>% droplevels()
    DGsig <- DGsig$gene
    DGvsdSig <- DGvsd %>% select(mouse,DGsig) %>% arrange(mouse)

    # join significant genes and time to first entrance
    DGsigbehav <- left_join(behav,DGvsdSig) %>% drop_na()  %>% arrange(mouse)

    ## Joining, by = "mouse"

    DGsigbehav <- as.data.frame(DGsigbehav)

    row.names(DGsigbehav) <- DGsigbehav$mouse 
    DGsigbehav$mouse <- NULL
    #head(DGsigbehav)

    # subset for sanes genes
    candidates <- c("Gria1", "Gria2", "Grin1", "Grin2a", "Grin2d",  "Prkcz" , "Prkci", "Wwc1")
    DGvsdCan <- DGvsd %>% select(mouse,candidates)
    DGcan <- left_join(behav,DGvsdCan) %>% drop_na()

    ## Joining, by = "mouse"

    DGcan <- as.data.frame(DGcan)
    row.names(DGcan) <- DGcan$mouse 
    DGcan$mouse <- NULL
    #head(DGcan)

    M <- cor(DGvsdbehav)

    ## Warning in cor(DGvsdbehav): the standard deviation is zero

    M <- as.data.frame(M)
    M$rownames <- row.names(M)

    Mslim <- M %>% filter(pTimeTarget > 0.95 | pTimeTarget < -0.95)
    row.names(Mslim) <- Mslim$rownames
    colstokeep <- Mslim$rownames
    Mslim <- Mslim %>% select(colstokeep)
    Mslim <- as.matrix(Mslim)
    corrplot.mixed(Mslim, number.cex = .7)

![](../figures/02e_correlations/corrplot-1.png)

    Mslim <- M %>% filter(Time1stEntr > 0.95 | Time1stEntr < -0.95)
    row.names(Mslim) <- Mslim$rownames
    colstokeep <- Mslim$rownames
    Mslim <- Mslim %>% select(colstokeep)
    Mslim <- as.matrix(Mslim)
    corrplot.mixed(Mslim, number.cex = .7)

![](../figures/02e_correlations/corrplot-2.png)

    M <- cor(DGcan)
    M <- as.data.frame(M)
    M$rownames <- row.names(M)
    row.names(M) <- M$rownames
    greatthan05 <- M$rownames
    M <- M %>% select(greatthan05)
    M <- as.matrix(M)
    corrplot.mixed(M, order = "hclust")

![](../figures/02e_correlations/corrplot-3.png)

    zDG <- ggplot(DGcan, aes(x = Time1stEntr, y = Prkcz)) +
      geom_point(aes( color = DGcols$Treatment))  +
      scale_color_manual(values = c("#f4a582", "#ca0020")) +
      theme(legend.position = "none") +
      geom_smooth(method='lm', color = "grey") +
      labs(title = "DG")

    iDG <- ggplot(DGcan, aes(x = Time1stEntr, y = Prkci )) +
      geom_point(aes(color = DGcols$Treatment)) + 
      scale_color_manual(values = c("#f4a582", "#ca0020")) +
      theme(legend.position = "none") +
      geom_smooth(method='lm', color = "grey")  +
      labs(title = "DG")


    plot_grid(zDG, iDG, 
              labels = c("R2 = -0.45" ,"R2 = 0.10"), 
              label_size = 10)

![](../figures/02e_correlations/pkcs-1.png)

    # retention plot
    behav <- read.csv("../data/01a_behavior.csv") 
    behav$mouse <- sapply(strsplit(as.character(behav$ID),"15"), "[", 2)
    retention <- behav %>% filter(APA2 %in% c("conflict-trained", "standard-trained"),
                                     TrainSession == "Retention") %>% 
                               select(mouse,APA2, Time1stEntr, Path1stEntr, pTimeTarget)
    retention$APA2 <- factor(retention$APA2, levels = c("standard-trained", "conflict-trained"))

    a <- ggplot(retention, aes(y = Time1stEntr, x = APA2, color = APA2)) +
      geom_point() + 
      scale_color_manual(values = c("#ca0020", "#f4a582")) +
      theme(legend.position = "none") +
      labs(y = "reten. time to 1st entr. (s)", x = NULL)

    b <- ggplot(retention, aes(y =pTimeTarget, x = APA2, color = APA2)) +
      geom_point() + 
      scale_color_manual(values = c("#ca0020", "#f4a582")) +
      theme(legend.position = "none") +
      labs(y = "reten. prop. time in shock zone", x = NULL)

    plot_grid(a,b)

![](../figures/02e_correlations/behav-1.png)
