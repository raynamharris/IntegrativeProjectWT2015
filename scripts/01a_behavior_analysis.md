This markdown file is used for behavioral data wrangling, statistical
analysis, and data visualization. Figures from this analysis were
assembled into this multi-panel plot using Adobe Illustrator. Files used
to create the individual figures are saved in the data subdirectory with
the prefix 01a.

Setup
-----

    ## load libraries 
    library(tidyverse) ## for respahing data
    library(plyr) ## for renmaing factors
    library(reshape2) ## for melting dataframe
    library(cowplot) ## for some easy to use themes
    library(ggfortify) # pca
    library(factoextra)  ## pca with vectors
    library(FactoMineR) # more pca
    library(car) ## stats
    library(pheatmap)  # for pretty heatmap
    library(viridis) # for awesome color pallette
    library(kableExtra) # for better markdown tables
    library(ggpubr) # for stats on figures

    ## load user-written functions 
    source("functions_behavior.R")
    source("figureoptions.R")

    ## set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/01_behavior/')

Sample sizes
------------

The ‘treatment’ column describes the four behavioral treatment groups.  
The ‘TrainSessionCombo’ column describes the behvioral training
sessions. Here I filter by a single session to calculte the number of
mice.

    ## import output from video tracker program 
    behavior <- read.csv("../data/01_behaviordata.csv", header = T)

    # rename to match RNAseq data, select useful variables
    behavior <- behavior %>%
      mutate(treatment = fct_recode(APA2,
                                    "standard.yoked" = "YokedSame",
                                    "standard.trained" = "Same",
                                    "conflict.yoked" = "YokedConflict",
                                    "conflict.trained" = "Conflict")) %>%
      mutate(training = fct_collapse(treatment,
                                          trained = c("standard.trained", "conflict.trained"),
                                          yoked = c("standard.yoked", "conflict.yoked"))) %>%
      select(ID,Day,treatment, training,TrainSessionCombo,TrainSessionComboNum, ShockOnOff, PairedPartner,
             SdevSpeedArena:Speed2) %>%
      arrange(ID) 

    ## rename columns 
    colnames(behavior)[colnames(behavior)=="pTimeTarget"] <- "pTimeShockZone"
    colnames(behavior)[colnames(behavior)=="TimeTarget"] <- "TimeShockZone"


    # set levels
    behavior$treatment <- factor(behavior$treatment, levels = c("standard.yoked", "standard.trained",
                                                              "conflict.yoked", "conflict.trained"))
    behavior$training <- factor(behavior$training, levels = c("yoked", "trained"))

    # sample sizes
    behavior %>% 
      filter(TrainSessionCombo == "Hab") %>%
      select(treatment)  %>%  summary()

    ##             treatment
    ##  standard.yoked  :8  
    ##  standard.trained:8  
    ##  conflict.yoked  :9  
    ##  conflict.trained:9

    head(behavior)

    ##       ID Day        treatment training TrainSessionCombo
    ## 1 15140A   1 conflict.trained  trained               Hab
    ## 2 15140A   1 conflict.trained  trained                T1
    ## 3 15140A   1 conflict.trained  trained                T2
    ## 4 15140A   1 conflict.trained  trained                T3
    ## 5 15140A   2 conflict.trained  trained            Retest
    ## 6 15140A   2 conflict.trained  trained             T4_C1
    ##   TrainSessionComboNum ShockOnOff PairedPartner SdevSpeedArena
    ## 1                    1        Off        15140B           3.07
    ## 2                    2         On        15140B           2.78
    ## 3                    3         On        15140B           2.68
    ## 4                    4         On        15140B           2.78
    ## 5                    5         On        15140B           3.11
    ## 6                    6         On        15140B           2.52
    ##   Linearity.Arena. NumEntrances Time1stEntr Path1stEntr Speed1stEntr.cm.s.
    ## 1           0.4790           28       24.63        1.09               4.56
    ## 2           0.4016            6        9.83        0.62              16.42
    ## 3           0.3170            2      118.37        3.17               2.31
    ## 4           0.3122            3      256.53        7.48               4.26
    ## 5           0.2895            1      432.07       10.56               9.38
    ## 6           0.3107           10        0.87        0.00              -1.00
    ##   Dist1stEntr.m. NumShock MaxTimeAvoid Time2ndEntr Path2ndEntr
    ## 1           1.12       52           53       59.97        2.59
    ## 2           0.30        7          327       18.30        1.23
    ## 3           0.11        3          312      287.63        8.54
    ## 4           0.17        3          256      447.80       12.74
    ## 5           0.06        1          432      599.97       15.66
    ## 6           0.56       13          447       25.90        0.75
    ##   Speed2ndEntr TimeShockZone pTimeShockZone pTimeCCW pTimeOPP pTimeCW
    ## 1         7.85        94.665         0.2277   0.2583   0.1788  0.3352
    ## 2         6.53         8.433         0.0211   0.6961   0.2049  0.0779
    ## 3         3.73         3.366         0.0092   0.6413   0.3245  0.0250
    ## 4         1.56         2.498         0.0069   0.5790   0.4018  0.0123
    ## 5        -1.00         1.067         0.0026   0.2945   0.6300  0.0729
    ## 6        16.19        17.735         0.0339   0.0195   0.1561  0.7905
    ##   RayleigLength RayleigAngle PolarAvgVal PolarSdVal PolarMinVal
    ## 1          0.11       330.67      163.83     106.39      0.0157
    ## 2          0.65       112.66      198.80      58.56      0.0004
    ## 3          0.78       124.87      214.13      41.85      0.0000
    ## 4          0.80       128.39      218.33      38.58      0.0000
    ## 5          0.72       159.36      246.61      46.73      0.0001
    ## 6          0.67       257.90      172.75      54.93      0.0005
    ##   PolarMinBin Min50.RngLoBin Min50.RngHiBin PolarMaxVal PolarMaxBin
    ## 1         160             60            280      0.0458         280
    ## 2         330            130             90      0.0976         100
    ## 3           0            150            120      0.1233         120
    ## 4           0            150            120      0.1175         140
    ## 5         340            170            120      0.0911         150
    ## 6         120            280            240      0.1003         280
    ##   Max50.RngLoBin Max50.RngHiBin AnnularMinVal AnnularMinBin AnnularMaxVal
    ## 1            230             40        0.0005          19.4        0.2807
    ## 2             70            140        0.0098          19.4        0.2773
    ## 3            100            160        0.0005           3.5        0.4294
    ## 4            100            160        0.0063           8.5        0.4233
    ## 5            110            190        0.0100           3.5        0.3537
    ## 6            230            300        0.0003           3.5        0.4390
    ##   AnnularMaxBin AnnularAvg AnnularSd AnnularSkewnes AnnularKurtosis
    ## 1          13.2      11.80     23.93           0.88            3.13
    ## 2          16.6      15.09     20.26           1.81            6.70
    ## 3          18.0      16.94     10.91           1.87            8.91
    ## 4          18.0      16.58     16.83           2.84           12.51
    ## 5          18.0      16.62     15.14           2.42           11.83
    ## 6          18.0      17.44      9.49           0.98            4.65
    ##       Speed1     Speed2
    ## 1 0.04425497 0.04318826
    ## 2 0.06307223 0.06721311
    ## 3 0.02678043 0.02969092
    ## 4 0.02915838 0.02845020
    ## 5 0.02444048 0.02610131
    ## 6 0.00000000 0.02895753

Number of shocks
----------------

The values in the column “NumShock” are actually measures of the number
of entraces into the shock zone. Because, that’s what the software
records. For standard.trained and conflict.trained animals, the number
of shocks equals equals the number of entraces. However, for yoked
individuals, the number of entrances does not equal the number of
shocks. For them, the number of shocks is equal to their
standard.trained or conflict.trained trained partner.

    # supset beahvior to keep only factors and num shocks
    numshocks <- behavior %>%
      select(ID, TrainSessionCombo, treatment, NumShock) 

    # widen datafram, and sum total
    numshocks <- spread(numshocks, key=TrainSessionCombo, value= NumShock)
    numshocks$sums <- rowSums(numshocks[sapply(numshocks, is.numeric)])
    head(numshocks)

    ##       ID        treatment Hab Retention Retest T1 T2 T3 T4_C1 T5_C2 T6_C3
    ## 1 15140A conflict.trained  52         9      1  7  3  3    13     6     2
    ## 2 15140B   conflict.yoked  55        33     71 96 30 71    87    31    32
    ## 3 15140C standard.trained  62         0     10  6  7  8     3     8     0
    ## 4 15140D   standard.yoked  61        41     34 58 32 22    32    54    48
    ## 5 15141C standard.trained  44        21      7  8 19 11    10     7    12
    ## 6 15141D   standard.yoked  55        52     50 54 48 27    47    36    34
    ##   sums
    ## 1   96
    ## 2  506
    ## 3  104
    ## 4  382
    ## 5  139
    ## 6  403

    # delete values for yoked animals
    numshocks <- numshocks %>%
      filter(treatment %in% c("standard.trained", "conflict.trained")) %>%
      droplevels()

    # create a tempdataframe with dupclicate values for yoked
    numshockstemp <- numshocks
    levels(numshockstemp$treatment) 

    ## [1] "standard.trained" "conflict.trained"

    levels(numshockstemp$treatment) <- c("standard.yoked","conflict.yoked")
    levels(numshockstemp$treatment) 

    ## [1] "standard.yoked" "conflict.yoked"

    # combine the two and plot

    realnumshocks <- rbind(numshocks, numshockstemp)
    levels(realnumshocks$treatment) 

    ## [1] "standard.trained" "conflict.trained" "standard.yoked"  
    ## [4] "conflict.yoked"

    realnumshocks$treatment <- factor(realnumshocks$treatment, levels = c("standard.yoked", "standard.trained", "conflict.yoked", "conflict.trained"))
    levels(realnumshocks$treatment) <- c("standard\nyoked", "standard\ntrained", "conflict\nyoked", "conflict\ntrained")

    realnumshocks %>%
      dplyr::group_by(treatment) %>%
      dplyr::summarise(meanshocks = mean(sums, na.rm = TRUE))

    ## # A tibble: 4 x 2
    ##   treatment           meanshocks
    ##   <fct>                    <dbl>
    ## 1 "standard\nyoked"          96 
    ## 2 "standard\ntrained"        96 
    ## 3 "conflict\nyoked"         125.
    ## 4 "conflict\ntrained"       125.

    # define what levels to compare for stats

    a <- ggplot(realnumshocks, aes(x = treatment, y = sums, fill = treatment)) +
      geom_boxplot(outlier.size = 0.5) +
      theme_ms() +
      scale_fill_manual(values = colorvalAPA00,
                        name = NULL) +
      labs(x = "treatment", subtitle = " ", y = "Total Shocks") +
        theme(axis.text.x=element_text(angle=60, vjust = 1, hjust = 1),
              legend.position = "none") 
    a

![](../figures/01_behavior/shockentrplot-1.png)

Vizualizing Mean and Standard error for num entrace and time 1st entrance
=========================================================================

To make the point and line graphs, I must create and join some data
frames, then I have a function that makes four plots with specific
titles, y labels and limits.

    dfb <- behavior %>%
      dplyr::group_by(treatment, TrainSessionComboNum) %>%
      dplyr::summarise(m = mean(NumEntrances), 
                       se = sd(NumEntrances)/sqrt(length(NumEntrances))) %>%
      dplyr::mutate(measure = "Number of target zone entrances")

    dfc <- behavior %>%
      dplyr::group_by(treatment, TrainSessionComboNum) %>%
      dplyr::mutate(minutes = Time1stEntr/60) %>%
      dplyr::summarise(m = mean(minutes), 
                       se = sd(minutes)/sqrt(length(minutes))) %>%
      dplyr::mutate(measure = "Time to 1st target zone entrance (min)")

    dfd <- behavior %>%
      dplyr::group_by(treatment, TrainSessionComboNum) %>%
      dplyr::summarise(m = mean(pTimeShockZone), 
                       se = sd(pTimeShockZone)/sqrt(length(pTimeShockZone))) %>%
      dplyr::mutate(measure = "Proportion of time in target zone")

    fourmeasures <- rbind(dfb,dfc,dfd)
    head(fourmeasures)

    ## # A tibble: 6 x 5
    ## # Groups:   treatment [1]
    ##   treatment      TrainSessionComboN…     m    se measure                   
    ##   <fct>                        <int> <dbl> <dbl> <chr>                     
    ## 1 standard.yoked                   1  31.4 2.32  Number of target zone ent…
    ## 2 standard.yoked                   2  21.4 2.02  Number of target zone ent…
    ## 3 standard.yoked                   3  15.4 1.40  Number of target zone ent…
    ## 4 standard.yoked                   4  14.5 2.01  Number of target zone ent…
    ## 5 standard.yoked                   5  16.9 0.875 Number of target zone ent…
    ## 6 standard.yoked                   6  15   1.56  Number of target zone ent…

    # see https://cran.r-project.org/web/packages/cowplot/vignettes/shared_legends.html for share legends


    b <- meansdplots(dfb, "NumEntrances" ,  c(0,10,20,30), c(0, 35)) + theme(legend.justification = "center")
    c <- meansdplots(dfc, "Time1stEntr.m (min)",  c(0,2,4,6,8), c(0, 8))
    d <- meansdplots(dfd, "pTimeShockZone", c(0,.12,.25,.37), c(0, .37 ))

    fourplots <- plot_grid(a + theme(legend.position="none"),
               b + theme(legend.position="none"),
               c + theme(legend.position="none"),
               d + theme(legend.position="none"),
               #align = 'vh',
               labels = c("(a)", "(b)", "(c)", "(d)"),
               nrow = 1,
               label_size = 8
               )
    fourplots

![](../figures/01_behavior/fourmeasures-1.png)

### individual variation in behavior

    behavior %>% 
        ggplot(aes(x=, TrainSessionComboNum, y=NumEntrances, color=ID)) + 
        geom_line() +
        geom_point(size = 1.5) +
        labs(subtitle = " ") +
        scale_x_continuous(name= "training session", 
                           breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                           labels = c( "P", "T1", "T2", "T3",
                                       "Rt", "T4", "T5", "T6", "Rn")) +
        theme_ms() +
      facet_wrap(~treatment)

![](../figures/01_behavior/individualbehavior-1.png)

    behavior %>% 
        filter(ID %in% c("15148B", "15143B","15148A", "15143A", "15145A"))  %>%
        ggplot(aes(x=, TrainSessionComboNum, y=NumEntrances, color=treatment, shape = ID)) + 
        geom_line() +
        geom_point(size = 1.5) +
        labs(subtitle = " ") +
        scale_x_continuous(name= "training session", 
                           breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                           labels = c( "P", "T1", "T2", "T3",
                                       "Rt", "T4", "T5", "T6", "Rn")) +
        theme_ms() +
      scale_color_manual(values = treatmentcolors)

![](../figures/01_behavior/individualbehavior-2.png)

### Principle component analysis

Next, I next reduced the dimentionality of the data with a PCA anlaysis.

    retention <- behavior %>% filter(TrainSessionCombo %in% c("Retention") ) %>% droplevels()

    pcadf <- makepcadf(behavior)

    pcadfretention <- pcadf %>% filter(TrainSessionComboNum == 9) %>% 
      group_by(treatment) %>% 
      dplyr::summarize(avePC1 = mean(PC1),
                       avePC2 = mean(PC2),
                       sePC1 = sd(PC1)/sqrt(length(PC1)),
                       sePC2 = sd(PC2)/sqrt(length(PC2)))
    pcadfretention

    ## # A tibble: 4 x 5
    ##   treatment        avePC1  avePC2 sePC1 sePC2
    ##   <fct>             <dbl>   <dbl> <dbl> <dbl>
    ## 1 standard.yoked     3.50 -2.19   0.452 0.594
    ## 2 standard.trained  -4.06 -0.645  0.982 0.292
    ## 3 conflict.yoked     2.50 -1.14   0.508 0.534
    ## 4 conflict.trained  -2.28 -0.0208 1.14  0.604

    e <- ggplot(pcadf, aes(x = PC1, y = PC2, color = treatment, fill = treatment)) +

      geom_point(data = pcadf, aes(alpha = TrainSessionComboNum)) + 
      geom_point(data = pcadfretention, aes(x = avePC1, y = avePC2), size = 4) +
      theme_ms() +
        scale_fill_manual(guide = 'none',values = colorvalAPA00) +
      scale_color_manual(guide = 'none',values = colorvalAPA00) +
      scale_alpha_continuous( breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                           labels = c( "P", "T1", "T2", "T3",
                                       "Rt", "T4", "T5", "T6", "Rn")) +
      theme(legend.position = "none") +
      labs(y = "PC2: 9.7% variance explained", x = "PC1: 35.7% variance explained",
           subtitle = " ") 
    e

![](../figures/01_behavior/PCA-1.png)

    # get contributions
    df <- behavior %>% select(SdevSpeedArena:Speed2)
    res.pca <- PCA(df,  graph = FALSE)
    # Visualize eigenvalues/variances
    fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))

![](../figures/01_behavior/PCA-2.png)

    # Contributions of variables to PC1
    f <- fviz_contrib_rmh(res.pca, choice = "var", axes = 1, top = 8, 
                     ylab = "PC1 % contributions", xlab = "estimates of memory", subtitle = " ") +
      theme_ms() + theme(axis.text.x = element_text(angle=45, hjust = 1))
    # Contributions of variables to PC2
    g <- fviz_contrib_rmh(res.pca, choice = "var", axes = 2, top = 8, 
                     ylab = "PC2 % contributions" , xlab = "estimates of activity", subtitle = " ") +
      theme_ms() + theme(axis.text.x = element_text(angle=45, hjust = 1))

    threeplots <- plot_grid(e,f,g, labels = c("(e)", "(f)", "(g)"),
               nrow = 1,
               label_size = 8,
              rel_widths = c(0.5,0.25,0.25))

    threeplots

![](../figures/01_behavior/PCA-3.png)

    behaviorfig <- plot_grid(fourplots, threeplots, nrow = 2)
    behaviorfig

![](../figures/01_behavior/behaviorfig-1.png)

    pdf(file="../figures/01_behavior/behaviorfig.pdf", width=7, height=3.5)
    plot(behaviorfig)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    write.csv(behavior, file = "../data/01a_behavior.csv", row.names = FALSE)
    write.csv(retention, file = "../data/01a_retention.csv", row.names = FALSE)
    write.csv(fourmeasures, file = "../data/01a_fourmeasures.csv", row.names = FALSE)
    write.csv(pcadf, file = "../data/01a_pcadf.csv", row.names = FALSE)
