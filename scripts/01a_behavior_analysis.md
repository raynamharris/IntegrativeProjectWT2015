This markdown file is used for behavioral data wrangling, statistical
analysis, and data visualization. Figures from this analysis were
assembled into this multi-panel plot using Adobe Illustrator. Files used
to create the individual figures are saved in the data subdirectory with
the prefix 01a.

Setup
-----

    ## load libraries 
    library(tidyverse) ## for respahing data
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

    grouped <- numshocks %>% group_by(treatment)
    summarise(grouped, mean=mean(sums), sd=sd(sums))

    ## # A tibble: 4 x 3
    ##   treatment         mean    sd
    ##   <fct>            <dbl> <dbl>
    ## 1 standard.yoked    473.  83.7
    ## 2 standard.trained   96   29.1
    ## 3 conflict.yoked    480.  87.4
    ## 4 conflict.trained  125.  32.5

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

    trainedonly <- realnumshocks %>% filter(treatment %in% c("conflict\ntrained", "standard\ntrained"))
    summary(trainedonly)

    ##        ID                treatment       Hab          Retention    
    ##  15140A : 1   standard\nyoked  :0   Min.   :44.00   Min.   : 0.00  
    ##  15140C : 1   standard\ntrained:8   1st Qu.:52.00   1st Qu.: 0.00  
    ##  15141C : 1   conflict\nyoked  :0   Median :57.00   Median : 6.00  
    ##  15142A : 1   conflict\ntrained:9   Mean   :56.24   Mean   :13.06  
    ##  15142C : 1                         3rd Qu.:61.00   3rd Qu.:28.00  
    ##  15143A : 1                         Max.   :66.00   Max.   :40.00  
    ##  (Other):11                                                        
    ##      Retest         T1               T2           T3        
    ##  Min.   : 0   Min.   : 4.000   Min.   : 1   Min.   : 1.000  
    ##  1st Qu.: 1   1st Qu.: 6.000   1st Qu.: 2   1st Qu.: 1.000  
    ##  Median : 2   Median : 7.000   Median : 3   Median : 1.000  
    ##  Mean   : 3   Mean   : 8.471   Mean   : 4   Mean   : 3.294  
    ##  3rd Qu.: 5   3rd Qu.:10.000   3rd Qu.: 4   3rd Qu.: 3.000  
    ##  Max.   :10   Max.   :25.000   Max.   :19   Max.   :15.000  
    ##                                                             
    ##      T4_C1           T5_C2            T6_C3             sums      
    ##  Min.   : 0.00   Min.   : 0.000   Min.   : 0.000   Min.   : 65.0  
    ##  1st Qu.: 3.00   1st Qu.: 1.000   1st Qu.: 1.000   1st Qu.: 92.0  
    ##  Median :13.00   Median : 4.000   Median : 2.000   Median :111.0  
    ##  Mean   :13.12   Mean   : 5.118   Mean   : 5.118   Mean   :111.4  
    ##  3rd Qu.:24.00   3rd Qu.: 6.000   3rd Qu.: 9.000   3rd Qu.:123.0  
    ##  Max.   :33.00   Max.   :28.000   Max.   :23.000   Max.   :200.0  
    ## 

    summary(aov(sums ~ treatment, data=trainedonly))

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## treatment    1   3589    3589   3.739 0.0723 .
    ## Residuals   15  14401     960                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

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

    dfa <- behavior %>%
      dplyr::group_by(treatment, TrainSessionComboNum) %>%
      dplyr::summarise(m = mean(Dist1stEntr.m.), 
                       se = sd(Dist1stEntr.m.)/sqrt(length(Dist1stEntr.m.))) %>%
      dplyr::mutate(measure = "Distance 1st entr")


    fourmeasures <- rbind(dfb,dfc,dfd, dfa)


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

    meansdplots(dfa, "Dist1stEntr.m." ,  c(0,0.5,1,1.5), c(0, 1.5)) + theme(legend.position = "none")

![](../figures/01_behavior/fourmeasures-1.png)

    hab <- behavior %>% filter(TrainSessionCombo == "Hab")
    hab

    ##        ID Day        treatment training TrainSessionCombo
    ## 1  15140A   1 conflict.trained  trained               Hab
    ## 2  15140B   1   conflict.yoked    yoked               Hab
    ## 3  15140C   1 standard.trained  trained               Hab
    ## 4  15140D   1   standard.yoked    yoked               Hab
    ## 5  15141C   1 standard.trained  trained               Hab
    ## 6  15141D   1   standard.yoked    yoked               Hab
    ## 7  15142A   1 conflict.trained  trained               Hab
    ## 8  15142B   1   conflict.yoked    yoked               Hab
    ## 9  15142C   1 standard.trained  trained               Hab
    ## 10 15142D   1   standard.yoked    yoked               Hab
    ## 11 15143A   1 conflict.trained  trained               Hab
    ## 12 15143B   1   conflict.yoked    yoked               Hab
    ## 13 15143C   1 standard.trained  trained               Hab
    ## 14 15143D   1   standard.yoked    yoked               Hab
    ## 15 15144A   1 conflict.trained  trained               Hab
    ## 16 15144B   1   conflict.yoked    yoked               Hab
    ## 17 15144C   1 standard.trained  trained               Hab
    ## 18 15144D   1   standard.yoked    yoked               Hab
    ## 19 15145A   1 conflict.trained  trained               Hab
    ## 20 15145B   1   conflict.yoked    yoked               Hab
    ## 21 15145C   1 standard.trained  trained               Hab
    ## 22 15145D   1   standard.yoked    yoked               Hab
    ## 23 15146A   1 conflict.trained  trained               Hab
    ## 24 15146B   1   conflict.yoked    yoked               Hab
    ## 25 15146C   1 standard.trained  trained               Hab
    ## 26 15146D   1   standard.yoked    yoked               Hab
    ## 27 15147A   1 conflict.trained  trained               Hab
    ## 28 15147B   1   conflict.yoked    yoked               Hab
    ## 29 15147C   1 standard.trained  trained               Hab
    ## 30 15147D   1   standard.yoked    yoked               Hab
    ## 31 15148A   1 conflict.trained  trained               Hab
    ## 32 15148B   1   conflict.yoked    yoked               Hab
    ## 33 15148C   1 conflict.trained  trained               Hab
    ## 34 15148D   1   conflict.yoked    yoked               Hab
    ##    TrainSessionComboNum ShockOnOff PairedPartner SdevSpeedArena
    ## 1                     1        Off        15140B           3.07
    ## 2                     1        Off        15140A           3.02
    ## 3                     1        Off        15140D           3.19
    ## 4                     1        Off        15140C           4.18
    ## 5                     1        Off        15141D           3.96
    ## 6                     1        Off        15141C           3.92
    ## 7                     1        Off        15142B           4.07
    ## 8                     1        Off        15142A           3.49
    ## 9                     1        Off        15142D           4.07
    ## 10                    1        Off        15142C           1.98
    ## 11                    1        Off        15143B           3.12
    ## 12                    1        Off        15143A           3.67
    ## 13                    1        Off        15143D           3.44
    ## 14                    1        Off        15143C           3.85
    ## 15                    1        Off        15144B           4.22
    ## 16                    1        Off        15144A           4.01
    ## 17                    1        Off        15144D           4.14
    ## 18                    1        Off        15144C           3.57
    ## 19                    1        Off        15145B           4.66
    ## 20                    1        Off        15145A           3.29
    ## 21                    1        Off        15145D           3.98
    ## 22                    1        Off        15145C           4.15
    ## 23                    1        Off        15146B           4.73
    ## 24                    1        Off        15146A           3.16
    ## 25                    1        Off        15146D           4.36
    ## 26                    1        Off        15146C           3.06
    ## 27                    1        Off        15147B           3.50
    ## 28                    1        Off        15147A           3.43
    ## 29                    1        Off        15147D           3.15
    ## 30                    1        Off        15147C           3.40
    ## 31                    1        Off        15148B           3.61
    ## 32                    1        Off        15148A           3.99
    ## 33                    1        Off        15148D           3.00
    ## 34                    1        Off        15148C           3.10
    ##    Linearity.Arena. NumEntrances Time1stEntr Path1stEntr
    ## 1            0.4790           28       24.63        1.09
    ## 2            0.4984           35       20.10        0.61
    ## 3            0.4678           28        0.53        0.00
    ## 4            0.5405           35        8.80        0.47
    ## 5            0.4678           23       14.30        0.84
    ## 6            0.5142           29        2.13        0.15
    ## 7            0.4975           34       17.43        0.95
    ## 8            0.4905           30       16.17        1.11
    ## 9            0.5227           34        0.50        0.00
    ## 10           0.3420           19       24.33        0.28
    ## 11           0.4609           29       12.43        0.76
    ## 12           0.4991           36        8.43        0.56
    ## 13           0.5004           31        0.53        0.00
    ## 14           0.5556           38       30.47        1.81
    ## 15           0.5421           33       20.23        0.99
    ## 16           0.5125           33       11.63        0.70
    ## 17           0.5015           42        1.33        0.15
    ## 18           0.5157           37        0.57        0.00
    ## 19           0.5416           39        3.10        0.17
    ## 20           0.4684           34        0.53        0.00
    ## 21           0.5037           31        9.27        0.72
    ## 22           0.5281           34       23.33        1.85
    ## 23           0.5499           37        1.97        0.01
    ## 24           0.4618           33        3.40        0.28
    ## 25           0.5104           35       15.80        1.18
    ## 26           0.4857           25        1.83        0.12
    ## 27           0.5008           32        6.10        0.48
    ## 28           0.4931           30       11.20        0.85
    ## 29           0.4407           25        6.87        0.35
    ## 30           0.4944           34        8.33        0.52
    ## 31           0.4780           28        0.50        0.00
    ## 32           0.5293           35       14.60        1.01
    ## 33           0.4104           25        2.73        0.17
    ## 34           0.4637           25       12.50        0.46
    ##    Speed1stEntr.cm.s. Dist1stEntr.m. NumShock MaxTimeAvoid Time2ndEntr
    ## 1                4.56           1.12       52           53       59.97
    ## 2                7.29           1.41       55           57       23.37
    ## 3               -1.00           1.14       62           62       11.30
    ## 4                9.58           1.08       61           45       18.47
    ## 5               17.47           0.81       44           90       25.57
    ## 6                3.76           0.97       55           77       12.07
    ## 7               18.08           1.09       60           40       41.20
    ## 8               18.95           1.12       57           58       29.10
    ## 9               -1.00           1.07       57           55       10.30
    ## 10               0.70           1.28       52           65       61.10
    ## 11              11.74           1.16       56           78       26.40
    ## 12              16.22           1.23       75           36       19.00
    ## 13              -1.00           1.15       52           46        7.30
    ## 14              15.95           1.22       59           30       45.27
    ## 15               8.10           1.02       64           48       33.00
    ## 16              14.26           1.08       63           47       22.63
    ## 17              11.79           1.33       64           33        5.77
    ## 18              -1.00           1.34       71           33       18.70
    ## 19               6.23           1.11       59           44       13.20
    ## 20              -1.00           1.31       55           44        8.90
    ## 21              11.98           0.99       54           42       18.27
    ## 22              12.72           1.10       55           50       29.73
    ## 23               7.27           1.00       66           50       10.83
    ## 24               2.66           1.29       65           35       18.60
    ## 25              15.05           1.04       58           42       26.67
    ## 26               8.98           1.04       51           59       27.10
    ## 27               8.84           1.17       48           56       13.73
    ## 28              12.71           1.15       56           69       29.13
    ## 29              20.14           1.03       49           51       24.90
    ## 30               9.04           1.32       59           50       24.23
    ## 31              -1.00           1.04       61           59        9.20
    ## 32              14.81           1.20       51           63       48.20
    ## 33               5.81           1.07       50           82       14.30
    ## 34              10.08           1.10       61           53       26.37
    ##    Path2ndEntr Speed2ndEntr TimeShockZone pTimeShockZone pTimeCCW pTimeOPP
    ## 1         2.59         7.85        94.665         0.2277   0.2583   0.1788
    ## 2         0.69         3.30        98.704         0.2476   0.2118   0.2462
    ## 3         0.80         9.86       114.934         0.2772   0.2821   0.2146
    ## 4         1.41        12.52       107.802         0.2634   0.2157   0.2261
    ## 5         1.58        17.50        83.731         0.2069   0.3228   0.1702
    ## 6         0.76         7.28       102.096         0.2473   0.2709   0.2433
    ## 7         2.33        13.86        97.965         0.2498   0.2483   0.2645
    ## 8         1.81        10.13       109.031         0.2565   0.2910   0.2232
    ## 9         0.93        14.60        96.867         0.2343   0.1785   0.2963
    ## 10        0.85         1.13       105.467         0.2676   0.4001   0.2192
    ## 11        1.52         5.98       105.065         0.2768   0.1913   0.2475
    ## 12        1.24        18.41       134.705         0.3323   0.3049   0.1820
    ## 13        0.31         3.93        96.040         0.2497   0.3902   0.1998
    ## 14        3.02        13.07        95.336         0.2326   0.2924   0.2684
    ## 15        1.90         9.99       122.296         0.3015   0.2639   0.2049
    ## 16        1.43        13.00       110.933         0.2744   0.2100   0.2393
    ## 17        0.33         6.20       115.464         0.2997   0.3136   0.2083
    ## 18        1.33         3.98       135.666         0.3363   0.2732   0.1923
    ## 19        1.27        16.78       116.634         0.2874   0.3077   0.1949
    ## 20        0.53         4.07       108.667         0.2898   0.2351   0.2174
    ## 21        1.26         9.85       113.835         0.2874   0.2815   0.2088
    ## 22        2.53        16.03       105.238         0.2601   0.3029   0.2160
    ## 23        0.87        21.44       112.933         0.2949   0.2624   0.2254
    ## 24        1.22         5.78       116.636         0.3035   0.2205   0.2324
    ## 25        2.36        14.58       101.073         0.2656   0.3268   0.1812
    ## 26        1.30         9.59        98.338         0.2516   0.2806   0.2741
    ## 27        0.93        13.00        81.968         0.2032   0.2868   0.2464
    ## 28        1.75         7.06       104.203         0.2670   0.3403   0.2173
    ## 29        1.68         9.10        88.002         0.2180   0.3152   0.2532
    ## 30        1.19        10.24       109.605         0.2776   0.3790   0.1582
    ## 31        0.73        11.03       117.834         0.3046   0.2771   0.2026
    ## 32        2.82         1.66        88.033         0.2259   0.3034   0.2352
    ## 33        1.01         6.63        97.966         0.2554   0.2866   0.2395
    ## 34        1.16         8.85       126.207         0.3237   0.2894   0.1942
    ##    pTimeCW RayleigLength RayleigAngle PolarAvgVal PolarSdVal PolarMinVal
    ## 1   0.3352          0.11       330.67      163.83     106.39      0.0157
    ## 2   0.2945          0.08       286.89      175.79     108.65      0.0159
    ## 3   0.2261          0.10        46.70      174.07      99.19      0.0179
    ## 4   0.2948          0.04       278.87      186.30     106.66      0.0152
    ## 5   0.3001          0.05        11.68      175.07     104.81      0.0136
    ## 6   0.2385          0.07       150.85      192.84     100.46      0.0129
    ## 7   0.2374          0.05       176.15      183.49     103.24      0.0157
    ## 8   0.2293          0.10        67.03      177.62      98.68      0.0165
    ## 9   0.2908          0.10       250.81      184.80     109.29      0.0164
    ## 10  0.1131          0.29        90.12      176.53      83.72      0.0057
    ## 11  0.2844          0.11       304.23      171.00     107.92      0.0129
    ## 12  0.1808          0.19        43.44      169.33      94.02      0.0168
    ## 13  0.1603          0.21        72.85      170.69      90.15      0.0098
    ## 14  0.2066          0.09       107.35      181.61      98.31      0.0182
    ## 15  0.2298          0.07        10.31      175.04     102.45      0.0187
    ## 16  0.2764          0.06       313.08      172.21     105.77      0.0161
    ## 17  0.1785          0.18        32.80      157.62      92.62      0.0104
    ## 18  0.1983          0.18        29.57      164.30      95.82      0.0123
    ## 19  0.2100          0.16        44.06      167.00      95.53      0.0130
    ## 20  0.2577          0.06       349.09      173.60     104.05      0.0154
    ## 21  0.2224          0.06        63.19      174.45     100.26      0.0180
    ## 22  0.2210          0.07        73.23      176.63      99.60      0.0149
    ## 23  0.2173          0.05        44.77      171.46     100.05      0.0172
    ## 24  0.2437          0.04       304.28      187.76     104.43      0.0167
    ## 25  0.2264          0.16        76.83      174.89      94.60      0.0156
    ## 26  0.1936          0.10        91.59      177.04      96.46      0.0138
    ## 27  0.2636          0.04       101.53      182.83     102.89      0.0159
    ## 28  0.1754          0.19        58.06      168.12      92.25      0.0136
    ## 29  0.2136          0.12       109.55      186.18      97.08      0.0165
    ## 30  0.1852          0.21        49.64      161.13      90.76      0.0113
    ## 31  0.2157          0.10        33.40      168.93      98.97      0.0142
    ## 32  0.2355          0.09        67.80      176.73      99.20      0.0123
    ## 33  0.2185          0.10        67.22      171.38      97.63      0.0148
    ## 34  0.1927          0.15        30.01      166.92      96.75      0.0149
    ##    PolarMinBin Min50.RngLoBin Min50.RngHiBin PolarMaxVal PolarMaxBin
    ## 1          160             60            280      0.0458         280
    ## 2          110             20            230      0.0443         330
    ## 3          170            150             10      0.0410          10
    ## 4          300            260            100      0.0393         260
    ## 5          180             50            250      0.0456         250
    ## 6          310            230             80      0.0437          80
    ## 7          200            160              0      0.0459         160
    ## 8          350            150             10      0.0455          20
    ## 9           80            330            190      0.0439         260
    ## 10         280            140             30      0.0569         110
    ## 11          90             30            250      0.0457         300
    ## 12         230             90            330      0.0473          80
    ## 13         200            160             40      0.0484          70
    ## 14         230            180             40      0.0365         170
    ## 15         270            140            350      0.0427         350
    ## 16          70              0            210      0.0378         310
    ## 17         250             70            310      0.0456         320
    ## 18         220            100            340      0.0530          30
    ## 19         180            100            330      0.0457          50
    ## 20         170             30            240      0.0462          30
    ## 21         160            120            330      0.0379         120
    ## 22         300            140              0      0.0557         110
    ## 23         260            140            350      0.0420         140
    ## 24         280            180             10      0.0459          10
    ## 25         170            140             20      0.0481         140
    ## 26         240            190             50      0.0440         190
    ## 27         350            220             60      0.0390          70
    ## 28         220            130             10      0.0432          50
    ## 29         340            190             50      0.0407         110
    ## 30         150            120             10      0.0486          60
    ## 31         270             70            290      0.0415         350
    ## 32           0            330            150      0.0411         330
    ## 33         240            150             10      0.0387          40
    ## 34         270             90            320      0.0438          30
    ##    Max50.RngLoBin Max50.RngHiBin AnnularMinVal AnnularMinBin AnnularMaxVal
    ## 1             230             40        0.0005          19.4        0.2807
    ## 2             170            350        0.0189          18.0        0.2089
    ## 3             310            120        0.0069          18.0        0.3405
    ## 4             110            300        0.0224          18.0        0.3004
    ## 5             210             40        0.0024          18.0        0.3270
    ## 6              50            230        0.0086          18.0        0.2525
    ## 7              10            200        0.0204          18.0        0.4025
    ## 8             350            160        0.0034          18.0        0.3046
    ## 9             140            310        0.0082          18.0        0.3731
    ## 10             30            160        0.0047          16.6        0.3551
    ## 11            240             50        0.0131          18.0        0.3490
    ## 12            320            110        0.0014          18.0        0.3649
    ## 13             20            170        0.0058          18.0        0.3355
    ## 14             30            200        0.0011          18.0        0.3134
    ## 15            280            100        0.0131          18.0        0.4132
    ## 16            240             60        0.0005          18.0        0.2699
    ## 17            270             70        0.0071          18.0        0.3614
    ## 18            330            120        0.0020          18.0        0.3367
    ## 19            310            110        0.0183          18.0        0.3314
    ## 20            240             60        0.0002          18.0        0.3134
    ## 21            340            160        0.0159          18.0        0.3335
    ## 22            340            150        0.0029          18.0        0.3221
    ## 23             60            260        0.0003          19.4        0.3499
    ## 24            300            130        0.0026          19.4        0.3379
    ## 25             20            170        0.0200          18.0        0.2856
    ## 26             60            240        0.0178          18.0        0.2585
    ## 27             10            190        0.0120          18.0        0.3075
    ## 28              0            150        0.0008          19.4        0.2849
    ## 29             10            170        0.0136          18.0        0.2717
    ## 30              0            140        0.0143          18.0        0.2737
    ## 31            270             80        0.0021          19.4        0.3481
    ## 32            150              0        0.0093          18.0        0.2734
    ## 33             10            180        0.0023          18.0        0.3001
    ## 34            290             90        0.0102          18.0        0.2719
    ##    AnnularMaxBin AnnularAvg AnnularSd AnnularSkewnes AnnularKurtosis
    ## 1           13.2      11.80     23.93           0.88            3.13
    ## 2           13.2      11.79     24.51           0.82            3.01
    ## 3           13.2      11.95     22.62           1.14            3.89
    ## 4           15.0      13.26     19.14           1.27            4.80
    ## 5           13.2      12.88     19.02           1.45            5.19
    ## 6           15.0      12.87     21.12           1.22            4.20
    ## 7           15.0      13.76     17.94           1.73            6.39
    ## 8           15.0      12.86     19.70           1.50            5.29
    ## 9           15.0      13.30     20.53           1.64            5.43
    ## 10           8.5       9.83     19.41           0.33            2.81
    ## 11          13.2      13.32     18.90           1.73            6.43
    ## 12          13.2      13.12     16.98           1.59            5.85
    ## 13          15.0      12.86     19.00           1.58            5.59
    ## 14          13.2      12.39     20.13           1.20            3.98
    ## 15          15.0      13.78     14.93           1.93            8.16
    ## 16          13.2      12.03     20.23           1.05            3.81
    ## 17          13.2      13.23     16.13           1.79            8.03
    ## 18          13.2      12.52     20.29           1.29            4.47
    ## 19          15.0      13.67     18.54           1.80            6.79
    ## 20          13.2      11.39     22.05           1.12            3.67
    ## 21          15.0      13.36     17.72           1.58            6.01
    ## 22          13.2      12.78     18.32           1.47            5.19
    ## 23          15.0      13.95     18.06           1.77            7.29
    ## 24          15.0      13.44     20.27           1.59            5.99
    ## 25          13.2      13.26     19.24           1.36            5.11
    ## 26          13.2      11.77     24.09           0.81            2.88
    ## 27          15.0      13.25     18.80           1.49            5.54
    ## 28          15.0      13.49     20.26           1.49            5.53
    ## 29          13.2      12.38     23.04           1.09            3.46
    ## 30          15.0      12.61     20.08           1.19            4.46
    ## 31          13.2      13.18     18.46           1.53            6.04
    ## 32          15.0      12.76     20.15           1.25            4.65
    ## 33          13.2      12.42     19.20           1.05            3.58
    ## 34          13.2      12.37     20.16           1.02            3.90
    ##         Speed1     Speed2
    ## 1  0.044254974 0.04318826
    ## 2  0.030348259 0.02952503
    ## 3  0.000000000 0.07079646
    ## 4  0.053409091 0.07634001
    ## 5  0.058741259 0.06179116
    ## 6  0.070422535 0.06296603
    ## 7  0.054503729 0.05655340
    ## 8  0.068645640 0.06219931
    ## 9  0.000000000 0.09029126
    ## 10 0.011508426 0.01391162
    ## 11 0.061142397 0.05757576
    ## 12 0.066429419 0.06526316
    ## 13 0.000000000 0.04246575
    ## 14 0.059402691 0.06671085
    ## 15 0.048937222 0.05757576
    ## 16 0.060189166 0.06319046
    ## 17 0.112781955 0.05719237
    ## 18 0.000000000 0.07112299
    ## 19 0.054838710 0.09621212
    ## 20 0.000000000 0.05955056
    ## 21 0.077669903 0.06896552
    ## 22 0.079297042 0.08509923
    ## 23 0.005076142 0.08033241
    ## 24 0.082352941 0.06559140
    ## 25 0.074683544 0.08848894
    ## 26 0.065573770 0.04797048
    ## 27 0.078688525 0.06773489
    ## 28 0.075892857 0.06007552
    ## 29 0.050946143 0.06746988
    ## 30 0.062424970 0.04911267
    ## 31 0.000000000 0.07934783
    ## 32 0.069178082 0.05850622
    ## 33 0.062271062 0.07062937
    ## 34 0.036800000 0.04398938

    summary(aov(NumEntrances ~ treatment, data=hab))

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment    3    7.0   2.338   0.085  0.967
    ## Residuals   30  820.8  27.358

### individual variation in behavior

    problemsamples1 <- behavior %>%   filter(ID %in% c("15148A", "15143A")) 
    problemsamples2 <- behavior %>%   filter(ID %in% c("15148B", "15143B")) 

    behav1 <- behavior %>% filter(treatment %in% "conflict.trained")
    behav2 <- behavior %>% filter(treatment %in% "conflict.yoked")


    ggplot(data = behav1, aes(x=TrainSessionComboNum, y= Time1stEntr, color=ID)) + 
      geom_point(data = problemsamples1, aes(x=TrainSessionComboNum, y=Time1stEntr, shape=ID, size = 1.5)) +
         geom_line(data = behav1, aes(x=TrainSessionComboNum, y=Time1stEntr)) +
        # geom_line() +
        labs(subtitle = " ") +
        scale_x_continuous(name= "training session", 
                           breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                           labels = c( "P", "T1", "T2", "T3",
                                       "Rt", "T4", "T5", "T6", "Rn")) +
        theme_ms() +
      facet_wrap(~treatment)

![](../figures/01_behavior/individualbehavior-1.png)

    ggplot(data = behav1, aes(x=TrainSessionComboNum, y= NumEntrances, color=ID)) + 
      geom_point(data = problemsamples1, aes(x=TrainSessionComboNum, y=NumEntrances, shape=ID, size = 1.5)) +
         geom_line(data = behav1, aes(x=TrainSessionComboNum, y=NumEntrances)) +
        # geom_line() +
        labs(subtitle = " ") +
        scale_x_continuous(name= "training session", 
                           breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                           labels = c( "P", "T1", "T2", "T3",
                                       "Rt", "T4", "T5", "T6", "Rn")) +
        theme_ms() +
      facet_wrap(~treatment)

![](../figures/01_behavior/individualbehavior-2.png)

    ggplot(data = behav1, aes(x=TrainSessionComboNum, y= pTimeShockZone, color=ID)) + 
      geom_point(data = problemsamples1, aes(x=TrainSessionComboNum, y=pTimeShockZone, shape=ID, size = 1.5)) +
         geom_line(data = behav1, aes(x=TrainSessionComboNum, y=pTimeShockZone)) +
        # geom_line() +
        labs(subtitle = " ") +
        scale_x_continuous(name= "training session", 
                           breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                           labels = c( "P", "T1", "T2", "T3",
                                       "Rt", "T4", "T5", "T6", "Rn")) +
        theme_ms() +
      facet_wrap(~treatment)

![](../figures/01_behavior/individualbehavior-3.png)

    ggplot(data = behav1, aes(x=TrainSessionComboNum, y= Time2ndEntr, color=ID)) + 
      geom_point(data = problemsamples1, aes(x=TrainSessionComboNum, y=Time2ndEntr, shape=ID, size = 1.5)) +
         geom_line(data = behav1, aes(x=TrainSessionComboNum, y=Time2ndEntr)) +
        # geom_line() +
        labs(subtitle = " ") +
        scale_x_continuous(name= "training session", 
                           breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                           labels = c( "P", "T1", "T2", "T3",
                                       "Rt", "T4", "T5", "T6", "Rn")) +
        theme_ms() +
      facet_wrap(~treatment)

![](../figures/01_behavior/individualbehavior-4.png)

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
