    ## load libraries 
    library(tidyverse) ## for respahing data
    library(cowplot) ## for some easy to use themes
    library(factoextra)  ## pca with vectors
    library(FactoMineR) # more pca
    library(apaTables) #  for ANOVA tables

    library(png) # for ading images to plots
    library(grid)  # for ading images to plots

    ## load user-written functions 
    source("functions_behavior.R")
    source("figureoptions.R")

    ## set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/01_behavior/')

Sample sizes
------------

The ‘treatment’ column describes the four behavioral treatment groups.  
The ‘trial’ column describes the behvioral training sessions. Here I
filter by a single session to calculte the number of mice.

    ## import output from video tracker program 
    behavior <- read.csv("../data/00_behaviordata.csv", header = T)

    # set levels
    behavior$treatment <- factor(behavior$treatment, levels = levelstreatment)
    behavior$training <- factor(behavior$training, levels = levelstraining)

    # sample sizes
    behavior %>% 
      filter(trial == "Hab") %>%
      select(treatment)  %>%  summary()

    ##             treatment
    ##  standard.yoked  :8  
    ##  standard.trained:8  
    ##  conflict.yoked  :9  
    ##  conflict.trained:9

    head(behavior)

    ##       ID Day        treatment training  trial trialNum ShockOnOff
    ## 1 15140A   1 conflict.trained  trained    Hab        1        Off
    ## 2 15140A   1 conflict.trained  trained     T1        2         On
    ## 3 15140A   1 conflict.trained  trained     T2        3         On
    ## 4 15140A   1 conflict.trained  trained     T3        4         On
    ## 5 15140A   2 conflict.trained  trained Retest        5         On
    ## 6 15140A   2 conflict.trained  trained  T4_C1        6         On
    ##   PairedPartner TotalPath.Arena. SpeedArena.cm.s sd.Speed.Arena.
    ## 1        15140B            22.68            3.78            3.07
    ## 2        15140B            19.36            3.23            2.78
    ## 3        15140B            15.01            2.50            2.68
    ## 4        15140B            14.39            2.40            2.78
    ## 5        15140B            14.04            2.34            3.11
    ## 6        15140B            12.50            2.08            2.52
    ##   Linearity.Arena. NumEntrances Time1stEntr Path1stEntr Speed1stEntr.cm.s.
    ## 1           0.4790           28       24.63        1.09               4.56
    ## 2           0.4016            6        9.83        0.62              16.42
    ## 3           0.3170            2      118.37        3.17               2.31
    ## 4           0.3122            3      256.53        7.48               4.26
    ## 5           0.2895            1      432.07       10.56               9.38
    ## 6           0.3107           10        0.87        0.00              -1.00
    ##   Entr.Dist.1.m. NumShock MaxTimeAvoid MaxPathAvoid Time2ndEntr
    ## 1           1.12        0           53         2.15       59.97
    ## 2           0.30        7          327        11.70       18.30
    ## 3           0.11        3          312         4.98      287.63
    ## 4           0.17        3          256         7.48      447.80
    ## 5           0.06        1          432        10.56      599.97
    ## 6           0.56       13          447         1.82       25.90
    ##   Path2ndEntr Speed2ndEntr TimeShockZone pTimeShockZone pTimeCCW pTimeOPP
    ## 1        2.59         7.85        94.665         0.2277   0.2583   0.1788
    ## 2        1.23         6.53         8.433         0.0211   0.6961   0.2049
    ## 3        8.54         3.73         3.366         0.0092   0.6413   0.3245
    ## 4       12.74         1.56         2.498         0.0069   0.5790   0.4018
    ## 5       15.66        -1.00         1.067         0.0026   0.2945   0.6300
    ## 6        0.75        16.19        17.735         0.0339   0.0195   0.1561
    ##   pTimeCW RayleigLength RayleigAngle Min50.RngLoBin AnnularSkewnes
    ## 1  0.3352          0.11       330.67             60           0.88
    ## 2  0.0779          0.65       112.66            130           1.81
    ## 3  0.0250          0.78       124.87            150           1.87
    ## 4  0.0123          0.80       128.39            150           2.84
    ## 5  0.0729          0.72       159.36            170           2.42
    ## 6  0.7905          0.67       257.90            280           0.98
    ##   AnnularKurtosis ShockPerEntrance
    ## 1            3.13         0.000000
    ## 2            6.70         1.166667
    ## 3            8.91         1.500000
    ## 4           12.51         1.000000
    ## 5           11.83         1.000000
    ## 6            4.65         1.300000

    names(behavior[9:34])

    ##  [1] "TotalPath.Arena."   "SpeedArena.cm.s"    "sd.Speed.Arena."   
    ##  [4] "Linearity.Arena."   "NumEntrances"       "Time1stEntr"       
    ##  [7] "Path1stEntr"        "Speed1stEntr.cm.s." "Entr.Dist.1.m."    
    ## [10] "NumShock"           "MaxTimeAvoid"       "MaxPathAvoid"      
    ## [13] "Time2ndEntr"        "Path2ndEntr"        "Speed2ndEntr"      
    ## [16] "TimeShockZone"      "pTimeShockZone"     "pTimeCCW"          
    ## [19] "pTimeOPP"           "pTimeCW"            "RayleigLength"     
    ## [22] "RayleigAngle"       "Min50.RngLoBin"     "AnnularSkewnes"    
    ## [25] "AnnularKurtosis"    "ShockPerEntrance"

    dfshocks <- behavior %>%
      dplyr::group_by(treatment, trialNum) %>%
      dplyr::summarise(m = mean(NumShock), 
                       se = sd(NumShock)/sqrt(length(NumShock))) %>%
      dplyr::mutate(m = round(m,0)) %>%
      dplyr::mutate(measure = "NumShock")
    dfshocks

    ## # A tibble: 36 x 5
    ## # Groups:   treatment [4]
    ##    treatment        trialNum     m    se measure 
    ##    <fct>               <int> <dbl> <dbl> <chr>   
    ##  1 standard.yoked          1     0 0     NumShock
    ##  2 standard.yoked          2     7 0.854 NumShock
    ##  3 standard.yoked          3     5 2.12  NumShock
    ##  4 standard.yoked          4     4 1.37  NumShock
    ##  5 standard.yoked          5     4 1.20  NumShock
    ##  6 standard.yoked          6     4 1.66  NumShock
    ##  7 standard.yoked          7     4 1.33  NumShock
    ##  8 standard.yoked          8     4 2.04  NumShock
    ##  9 standard.yoked          9     0 0     NumShock
    ## 10 standard.trained        1     0 0     NumShock
    ## # … with 26 more rows

    numshocks <- ggplot(dfshocks,  aes(x=, trialNum, y=m, color=treatment, label = m)) + 
        geom_errorbar(aes(ymin=m-se, ymax=m+se, color=treatment), width=.1) +
        geom_line() +
        geom_point(size = 1) +
        labs(y = "NumShocks") +
        scale_x_continuous(name= "trial", 
                           breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                           labels = c( "P", "T1", "T2", "T3",
                                       "Rt", "T4", "T5", "T6", "Rn")) +
        scale_alpha_continuous( breaks = c(1, 2, 3)) +
        theme_ms() +
        scale_color_manual(values = treatmentcolors,
                           name  = NULL)  +
        theme(legend.position = "none",
              legend.justification=c(0,0),
              legend.text=element_text(size=5),
              strip.text = element_blank(),
              axis.text.y = element_blank()) +
      facet_wrap(~treatment, nrow = 4) +
      scale_y_continuous(breaks = c(0,10,20), limits = c(0,30)) +
      geom_text(vjust= -1, size=2.5)

    numshocks

![](../figures/01_behavior/experimentaldesign-1.png)

    paradigm <- png::readPNG("../figures/00_schematics/figure_1a.png")
    paradigm <- ggdraw() +  draw_image(paradigm, scale = 1)

    subfields <- png::readPNG("../figures/00_schematics/figure_1c.png")
    subfields <- ggdraw() +  draw_image(subfields, scale = 1)


    experimentaldesign <- plot_grid(paradigm, numshocks, subfields, rel_widths = c(2,1,0.8), nrow = 1,
              label_size = 8, labels = c("(a)", "(b) ", "(c)"))
    experimentaldesign

![](../figures/01_behavior/experimentaldesign-2.png)

    pdf(file="../figures/01_behavior/experimentaldesign.pdf", width=6.69, height=2.4)
    plot(experimentaldesign)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    pdf(file="../figures/figure_1.pdf", width=6.69, height=2.4)
    plot(experimentaldesign)
    dev.off()

    ## quartz_off_screen 
    ##                 2

Vizualizing Mean and Standard error for avoidance behaviors
===========================================================

To make the point and line graphs, I must create and join some data
frames, then I have a function that makes four plots with specific
titles, y labels and limits.

    dfa <- behavior %>%
      dplyr::group_by(treatment, trialNum) %>%
      dplyr::summarise(m = mean(NumEntrances), 
                       se = sd(NumEntrances)/sqrt(length(NumEntrances))) %>%
      dplyr::mutate(measure = "NumEntrances")

    dfb <- behavior %>%
      dplyr::group_by(treatment, trialNum) %>%
      dplyr::mutate(minutes = Time1stEntr/60) %>%
      dplyr::summarise(m = mean(minutes), 
                       se = sd(minutes)/sqrt(length(minutes))) %>%
      dplyr::mutate(measure = "Time1stEntr.min")

    dfc <- behavior %>%
      dplyr::group_by(treatment, trialNum) %>%
      dplyr::summarise(m = mean(pTimeShockZone), 
                       se = sd(pTimeShockZone)/sqrt(length(pTimeShockZone))) %>%
      dplyr::mutate(measure = "pTimeShockZone")

    avoidancedf <- rbind(dfa, dfb,dfc)
    head(avoidancedf)

    ## # A tibble: 6 x 5
    ## # Groups:   treatment [1]
    ##   treatment      trialNum     m    se measure     
    ##   <fct>             <int> <dbl> <dbl> <chr>       
    ## 1 standard.yoked        1  31.4 2.32  NumEntrances
    ## 2 standard.yoked        2  21.4 2.02  NumEntrances
    ## 3 standard.yoked        3  15.4 1.40  NumEntrances
    ## 4 standard.yoked        4  14.5 2.01  NumEntrances
    ## 5 standard.yoked        5  16.9 0.875 NumEntrances
    ## 6 standard.yoked        6  15   1.56  NumEntrances

    a <- meansdplots(dfa, "NumEntrances" ,  c(0,10,20,30), c(0, 35)) 
    b <- meansdplots(dfb, "Time1stEntr.min",  c(0,2,4,6,8), c(0, 8))
    c <- meansdplots(dfc, "pTimeShockZone", c(0,.12,.25,.37), c(0, .37 ))

    a

![](../figures/01_behavior/behavmeanstdev-1.png)

    b

![](../figures/01_behavior/behavmeanstdev-2.png)

    c

![](../figures/01_behavior/behavmeanstdev-3.png)

    avoidancebehaviors <- plot_grid(a + theme(legend.position = "none"),
                           b + theme(legend.position = "none"), 
                           c + theme(legend.position = "none"), nrow = 1,
                           label_size = 8,
                           labels = c("(a)", "(b)", "(c)"))
    avoidancebehaviors

![](../figures/01_behavior/behavmeanstdev-4.png)

### Principle component analysis

Next, I next reduced the dimentionality of the data with a PCA anlaysis.

    pca.all <- makepcadf(behavior)

    retention <- behavior %>% filter(trialNum == 9)
    pca.Rn <- makepcadf(retention)

    pca.Rn.summary <- pca.all %>% filter(trialNum == 9) %>% 
      group_by(treatment) %>% 
      dplyr::summarize(avePC1 = mean(PC1),
                       avePC2 = mean(PC2),
                       sePC1 = sd(PC1)/sqrt(length(PC1)),
                       sePC2 = sd(PC2)/sqrt(length(PC2)))


    d <- ggplot(pca.all, aes(x = PC1, y = PC2, color = treatment, fill = treatment)) +

      geom_point(data = pca.all, aes(alpha = Day)) + 
      geom_point(data = pca.Rn.summary, aes(x = avePC1, y = avePC2), size = 4) +
      theme_ms() +
        scale_fill_manual(guide = 'none',values = treatmentcolors) +
      scale_color_manual(guide = 'none',values = treatmentcolors) +
      scale_alpha_continuous(breaks = c(1, 2, 3)) +
      theme(legend.position = "none") +
      labs( x = "PC1: 38.3% variance explained",
            y = "PC2: 16.7% \n variance explained",
           subtitle = " ") 
    d

![](../figures/01_behavior/PCA-1.png)

    # get contributions
    df <- behavior %>% select(TotalPath.Arena.:AnnularKurtosis)
    res.pca <- PCA(df,  graph = FALSE)
    # Visualize eigenvalues/variances
    fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))

![](../figures/01_behavior/PCA-2.png)

    # Contributions of variables to PC1
    e <- fviz_contrib_rmh(res.pca, choice = "var", axes = 1, top = 8, 
                     ylab = "PC1 % contrib.", xlab = "estimates of memory", subtitle = " ") +
      theme_ms() + theme(axis.text.x = element_text(angle=45, hjust = 1))
    # Contributions of variables to PC2
    f <- fviz_contrib_rmh(res.pca, choice = "var", axes = 2, top = 8, 
                     ylab = "PC2 % contrib." , xlab = "estimates of activity", subtitle = " ") +
      theme_ms() + theme(axis.text.x = element_text(angle=45, hjust = 1))

    e

![](../figures/01_behavior/PCA-3.png)

    f

![](../figures/01_behavior/PCA-4.png)

    pcaplots <- plot_grid(d,e,f, labels = c("(d)", "(e)", "(f)"),
               nrow = 1,
               label_size = 8)
    pcaplots

![](../figures/01_behavior/PCA-5.png)

    avoidance <- plot_grid(avoidancebehaviors, pcaplots, nrow = 2)
    avoidance

![](../figures/01_behavior/avoidance-1.png)

    pdf(file="../figures/01_behavior/avoidance.pdf", width=6.69, height=3.5)
    plot(avoidance)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    pdf(file="../figures/figure_2.pdf", width=6.69, height=3.5)
    plot(avoidance)
    dev.off()

    ## quartz_off_screen 
    ##                 2

now all the stats
-----------------

Maybe it should be done like this: Q1: are the groups different? 1-way
ANOVA of groups on Pre Q2: are the groups different during initial
training T1-T3? 2-way ANOVA of groups X trial Q3: do the groups differ
in initial recall? 1-way ANOVA of groups on Rt Q4: Do the groups differ
in subsequent training? T4-T6 2-way ANOVA of groups X trial Q5: do the
groups differ in subsequent recall? 1-way ANOVA of groups on Rn

    # twoway anova table function

    twowayANOVAfor3measures <- function(mydata, mydescription){
      apa1 <- apa.aov.table(aov(NumEntrances ~ treatment * trial, data=mydata))
      apa1df <- as.data.frame(apa1$table_body)
      totaldf <- apa1df[5, 3]
      apa1df$df <- paste(apa1df$df, ", " , totaldf, sep = "")
      apa1df$ANOVA <- "NumEntrances ~ treatment * trial"
      apa1df

      apa2 <- apa.aov.table(aov(pTimeShockZone ~ treatment * trial, data=mydata))
      apa2df <- as.data.frame(apa2$table_body)
      apa2df$df <- paste(apa2df$df, ", " , totaldf, sep = "")
      apa2df$ANOVA <- "pTimeShockZone ~ treatment * trial"

      apa3 <- apa.aov.table(aov(Time1stEntr ~ treatment * trial, data=mydata))
      apa3df <- as.data.frame(apa3$table_body)
      apa3df$df <- paste(apa3df$df, ", " , totaldf, sep = "")
      apa3df$ANOVA <- "Time1stEntr ~ treatment * trial"
      apa3df

      apa123 <- as.data.frame(rbind(apa1df,apa2df,apa3df))
      apa123$trials <- mydescription
      apa123 <- apa123 %>%
        select(trials, ANOVA, Predictor, df, "F", p) %>%
        filter(!Predictor %in% c("(Intercept)", "Error"))

      return(apa123)
    }


    onewayANOVAfor3measures <- function(mydata, whichtrial, mydescription){
      
      mydata <- mydata %>% filter(trial == whichtrial)
      
      apa1 <- apa.aov.table(aov(NumEntrances ~ treatment , data=mydata))
      apa1df <- as.data.frame(apa1$table_body)
      totaldf <- apa1df[3, 3]
      apa1df$df <- paste(apa1df$df, ", " , totaldf, sep = "")
      apa1df$ANOVA <- "NumEntrances ~ treatment"
      apa1df

      apa2 <- apa.aov.table(aov(pTimeShockZone ~ treatment , data=mydata))
      apa2df <- as.data.frame(apa2$table_body)
      apa2df$df <- paste(apa2df$df, ", " , totaldf, sep = "")
      apa2df$ANOVA <- "pTimeShockZone ~ treatment"

      apa3 <- apa.aov.table(aov(Time1stEntr ~ treatment , data=mydata))
      apa3df <- as.data.frame(apa3$table_body)
      apa3df$df <- paste(apa3df$df, ", " , totaldf, sep = "")
      apa3df$ANOVA <- "Time1stEntr ~ treatment"
      apa3df

      apa123 <- as.data.frame(rbind(apa1df,apa2df,apa3df))
      apa123$trials <- mydescription
      apa123 <- apa123 %>%
        select(trials, ANOVA, Predictor, df, "F", p) %>%
        filter(!Predictor %in% c("(Intercept)", "Error"))

      return(apa123)
    }

    # Q1. Are groups different at pre? No.
    Q1 <- onewayANOVAfor3measures(behavior, "Hab", "Pre-training (Pre)")
      
      
    # Q2. Are the groups different during initial training T1-T3? Yes (sometimes alone, sometime interaction)
    T1T2T3 <-  behavior %>% filter(trial %in% c("T1", "T2", "T3"))
    Q2 <- twowayANOVAfor3measures(T1T2T3, "Initial training (T1 - T3)")


    # Q3 Do the groups differ in initial recall? Yes.
    Q3 <- onewayANOVAfor3measures(behavior, "Retest", "Initial recall (Rt)")

    # Q4 Do the groups differ in subsequent training? Yes
    T4T5T6 <-  behavior %>% filter(trial %in% c("T4_C1", "T5_C2", "T6_C3"))
    Q4 <- twowayANOVAfor3measures(T4T5T6,  "Conflict training  (T4 - T6)")


    # Q5  Do the groups differ in subsequent recall? Yes
    Q5 <- onewayANOVAfor3measures(behavior, "Retention", "Conflict recall (Rn)")


    # more stats. didn't bother with a function

    PC1all <- apa.aov.table(aov(PC1 ~ treatment , data=pca.all))
    PC1all <- as.data.frame(PC1all$table_body)
    totaldf <- PC1all[3, 3]
    PC1all$df <- paste(PC1all$df, ", " , totaldf, sep = "")
    PC1all$ANOVA <- "PC1 ~ treatment"
    PC1all$trials <- "All trials"

    PC2all <- apa.aov.table(aov(PC2 ~ treatment , data=pca.all))
    PC2all <- as.data.frame(PC2all$table_body)
    totaldf <- PC2all[3, 3]
    PC2all$df <- paste(PC2all$df, ", " , totaldf, sep = "")
    PC2all$ANOVA <- "PC2 ~ treatment"
    PC2all$trials <- "All trails"

    PC1rn <- apa.aov.table(aov(PC1 ~ treatment , data=pca.Rn))
    PC1rn <- as.data.frame(PC1rn$table_body)
    totaldf <- PC1rn[3, 3]
    PC1rn$df <- paste(PC1rn$df, ", " , totaldf, sep = "")
    PC1rn$ANOVA <- "PC1 ~ treatment"
    PC1rn$trials <- "Retention (Rn)"

    PC2rn <- apa.aov.table(aov(PC2 ~ treatment , data=pca.Rn))
    PC2rn <- as.data.frame(PC2rn$table_body)
    totaldf <- PC2rn[3, 3]
    PC2rn$df <- paste(PC2rn$df, ", " , totaldf, sep = "")
    PC2rn$ANOVA <- "PC2 ~ treatment"
    PC2rn$trials <- "Retention"

    PC.APA <- as.data.frame(rbind(PC1all, PC2all, PC1rn, PC2rn))
    PC.APA <- PC.APA %>% 
        select(trials, ANOVA, Predictor, df, "F", p) %>%
        filter(!Predictor %in% c("(Intercept)", "Error"))

    Q12345 <- as.data.frame(rbind(Q1,Q2,Q3,Q4,Q5, PC.APA))
    Q12345

    ##                          trials                              ANOVA
    ## 1            Pre-training (Pre)           NumEntrances ~ treatment
    ## 2            Pre-training (Pre)         pTimeShockZone ~ treatment
    ## 3            Pre-training (Pre)            Time1stEntr ~ treatment
    ## 4    Initial training (T1 - T3)   NumEntrances ~ treatment * trial
    ## 5    Initial training (T1 - T3)   NumEntrances ~ treatment * trial
    ## 6    Initial training (T1 - T3)   NumEntrances ~ treatment * trial
    ## 7    Initial training (T1 - T3) pTimeShockZone ~ treatment * trial
    ## 8    Initial training (T1 - T3) pTimeShockZone ~ treatment * trial
    ## 9    Initial training (T1 - T3) pTimeShockZone ~ treatment * trial
    ## 10   Initial training (T1 - T3)    Time1stEntr ~ treatment * trial
    ## 11   Initial training (T1 - T3)    Time1stEntr ~ treatment * trial
    ## 12   Initial training (T1 - T3)    Time1stEntr ~ treatment * trial
    ## 13          Initial recall (Rt)           NumEntrances ~ treatment
    ## 14          Initial recall (Rt)         pTimeShockZone ~ treatment
    ## 15          Initial recall (Rt)            Time1stEntr ~ treatment
    ## 16 Conflict training  (T4 - T6)   NumEntrances ~ treatment * trial
    ## 17 Conflict training  (T4 - T6)   NumEntrances ~ treatment * trial
    ## 18 Conflict training  (T4 - T6)   NumEntrances ~ treatment * trial
    ## 19 Conflict training  (T4 - T6) pTimeShockZone ~ treatment * trial
    ## 20 Conflict training  (T4 - T6) pTimeShockZone ~ treatment * trial
    ## 21 Conflict training  (T4 - T6) pTimeShockZone ~ treatment * trial
    ## 22 Conflict training  (T4 - T6)    Time1stEntr ~ treatment * trial
    ## 23 Conflict training  (T4 - T6)    Time1stEntr ~ treatment * trial
    ## 24 Conflict training  (T4 - T6)    Time1stEntr ~ treatment * trial
    ## 25         Conflict recall (Rn)           NumEntrances ~ treatment
    ## 26         Conflict recall (Rn)         pTimeShockZone ~ treatment
    ## 27         Conflict recall (Rn)            Time1stEntr ~ treatment
    ## 28                   All trials                    PC1 ~ treatment
    ## 29                   All trails                    PC2 ~ treatment
    ## 30               Retention (Rn)                    PC1 ~ treatment
    ## 31                    Retention                    PC2 ~ treatment
    ##            Predictor     df      F    p
    ## 1          treatment  3, 30   0.09 .967
    ## 2          treatment  3, 30   0.78 .512
    ## 3          treatment  3, 30   0.80 .506
    ## 4          treatment  3, 90  26.42 .000
    ## 5              trial  2, 90   5.88 .004
    ## 6  treatment x trial  6, 90   0.72 .631
    ## 7          treatment  3, 90  48.26 .000
    ## 8              trial  2, 90   0.88 .419
    ## 9  treatment x trial  6, 90   0.59 .736
    ## 10         treatment  3, 90   0.02 .997
    ## 11             trial  2, 90   0.11 .895
    ## 12 treatment x trial  6, 90   3.24 .006
    ## 13         treatment  3, 30  45.44 .000
    ## 14         treatment  3, 30 129.49 .000
    ## 15         treatment  3, 30   7.96 .000
    ## 16         treatment  3, 90   9.37 .000
    ## 17             trial  2, 90   0.09 .913
    ## 18 treatment x trial  6, 90   3.43 .004
    ## 19         treatment  3, 90  25.26 .000
    ## 20             trial  2, 90   0.03 .972
    ## 21 treatment x trial  6, 90   1.54 .174
    ## 22         treatment  3, 90   6.03 .001
    ## 23             trial  2, 90   0.05 .954
    ## 24 treatment x trial  6, 90   1.74 .120
    ## 25         treatment  3, 30  18.01 .000
    ## 26         treatment  3, 30  26.90 .000
    ## 27         treatment  3, 30   5.97 .003
    ## 28         treatment 3, 302  91.83 .000
    ## 29         treatment 3, 302  10.76 .000
    ## 30         treatment  3, 30  17.69 .000
    ## 31         treatment  3, 30   0.39 .761

save files
----------

    # supp table 1
    write.csv(behavior, file = "../data/01_behavior.csv", row.names = FALSE)

    # supp table 2
    write.csv(avoidancedf, file = "../data/01_avoidance.csv", row.names = FALSE)

    # table 1
    write.csv(Q12345, "../data/01_APA.csv", row.names = F)

    # unsure wether/how to include
    write.csv(pca.all, file = "../data/01_pca.all.csv", row.names = FALSE)
    write.csv(pca.Rn.summary, file = "../data/01_pca.Rn.summary.csv", row.names = FALSE)
    write.csv(pca.Rn, file = "../data/01_pca.Rn.csv", row.names = FALSE)

    citation("tidyverse")

    ## 
    ##   Wickham et al., (2019). Welcome to the tidyverse. Journal of
    ##   Open Source Software, 4(43), 1686,
    ##   https://doi.org/10.21105/joss.01686
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Article{,
    ##     title = {Welcome to the {tidyverse}},
    ##     author = {Hadley Wickham and Mara Averick and Jennifer Bryan and Winston Chang and Lucy D'Agostino McGowan and Romain François and Garrett Grolemund and Alex Hayes and Lionel Henry and Jim Hester and Max Kuhn and Thomas Lin Pedersen and Evan Miller and Stephan Milton Bache and Kirill Müller and Jeroen Ooms and David Robinson and Dana Paige Seidel and Vitalie Spinu and Kohske Takahashi and Davis Vaughan and Claus Wilke and Kara Woo and Hiroaki Yutani},
    ##     year = {2019},
    ##     journal = {Journal of Open Source Software},
    ##     volume = {4},
    ##     number = {43},
    ##     pages = {1686},
    ##     doi = {10.21105/joss.01686},
    ##   }

    citation("cowplot")  

    ## 
    ## To cite package 'cowplot' in publications use:
    ## 
    ##   Claus O. Wilke (2019). cowplot: Streamlined Plot Theme and Plot
    ##   Annotations for 'ggplot2'. R package version 0.9.4.
    ##   https://CRAN.R-project.org/package=cowplot
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {cowplot: Streamlined Plot Theme and Plot Annotations for 'ggplot2'},
    ##     author = {Claus O. Wilke},
    ##     year = {2019},
    ##     note = {R package version 0.9.4},
    ##     url = {https://CRAN.R-project.org/package=cowplot},
    ##   }

    citation("factoextra")   

    ## 
    ## To cite package 'factoextra' in publications use:
    ## 
    ##   Alboukadel Kassambara and Fabian Mundt (2017). factoextra:
    ##   Extract and Visualize the Results of Multivariate Data Analyses.
    ##   R package version 1.0.5.
    ##   https://CRAN.R-project.org/package=factoextra
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {factoextra: Extract and Visualize the Results of Multivariate Data Analyses},
    ##     author = {Alboukadel Kassambara and Fabian Mundt},
    ##     year = {2017},
    ##     note = {R package version 1.0.5},
    ##     url = {https://CRAN.R-project.org/package=factoextra},
    ##   }

    citation("FactoMineR")  

    ## 
    ## To cite FactoMineR in publications use:
    ## 
    ##   Sebastien Le, Julie Josse, Francois Husson (2008). FactoMineR:
    ##   An R Package for Multivariate Analysis. Journal of Statistical
    ##   Software, 25(1), 1-18. 10.18637/jss.v025.i01
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Article{,
    ##     title = {{FactoMineR}: A Package for Multivariate Analysis},
    ##     author = {S\'ebastien L\^e and Julie Josse and Fran\c{c}ois Husson},
    ##     journal = {Journal of Statistical Software},
    ##     year = {2008},
    ##     volume = {25},
    ##     number = {1},
    ##     pages = {1--18},
    ##     doi = {10.18637/jss.v025.i01},
    ##   }

    citation("apaTables")  

    ## 
    ## To cite package 'apaTables' in publications use:
    ## 
    ##   David Stanley (2018). apaTables: Create American Psychological
    ##   Association (APA) Style Tables. R package version 2.0.5.
    ##   https://CRAN.R-project.org/package=apaTables
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {apaTables: Create American Psychological Association (APA) Style Tables},
    ##     author = {David Stanley},
    ##     year = {2018},
    ##     note = {R package version 2.0.5},
    ##     url = {https://CRAN.R-project.org/package=apaTables},
    ##   }
