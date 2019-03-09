This markdown file is used for behavioral data wrangling, statistical
analysis, and data visualization. Figures from this analysis were
assembled into this multi-panel plot using Adobe Illustrator. Files used
to create the individual figures are saved in the data subdirectory with
the prefix 01a.

<img src="../figures/figure_2.png" width="1370" />

Setup
-----

    ## load libraries 
    library(tidyr) ## for respahing data
    library(plyr) ## for renmaing factors
    library(dplyr) ## for filtering and selecting rows
    library(reshape2) ## for melting dataframe
    library(ggplot2) ## for awesome plots!
    library(cowplot) ## for some easy to use themes
    library(factoextra)  ##pca with vectors
    library(car) ## stats
    library(pheatmap)  # for pretty heatmap
    library(viridis) # for awesome color pallette
    library(kableExtra) # for better markdown tables

    ## load user-written functions 
    source("functions_behavior.R")
    source("figureoptions.R")

    ## set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/01_behavior/')

Sample sizes
------------

The ‘APA2’ column describes the four behavioral treatment groups.  
The ‘TrainSessionCombo’ column describes the behvioral training
sessions. Here I filter by a single session to calculte the number of
mice.

    ## import output from video tracker program 
    behavior <- read.csv("../data/01_behaviordata.csv", header = T)

    ## set level for APA2 then renmae
    behavior$APA2 <- factor(behavior$APA2, levels = c("YokedSame", "Same", "YokedConflict","Conflict"))
    levels(behavior$APA2) <-  c("yoked-consistent" ,"consistent", "yoked-conflict", "conflict")

    # sample sizes
    behavior %>% 
      filter(TrainSessionCombo == "Retention") %>%
      select(APA2)  %>%  summary()

    ##                APA2  
    ##  yoked-consistent:8  
    ##  consistent      :8  
    ##  yoked-conflict  :9  
    ##  conflict        :9

    # keep only relevant columns
    behavior_slim <- behavior[,c(15,16,14,20:58)] 

Vizualizing Mean and Standard error for num entrace and time 1st entrance
=========================================================================

To make the point and line graphs, I must create and merge some data
frames

    ## number of entrances
    meannumentr <- dplyr::summarise(group_by(behavior, APA2, TrainSessionComboNum), m = mean(NumEntrances), se = sd(NumEntrances)/sqrt(length(NumEntrances)))
    meannumentr$APA2 <- factor(meannumentr$APA2, levels = c("yoked-consistent" ,"consistent", "yoked-conflict", "conflict"))

    numentr <- ggplot(meannumentr, 
                      aes(x=, TrainSessionComboNum, y=m, color=APA2)) + 
        geom_errorbar(aes(ymin=m-se, ymax=m+se, color=APA2), width=.1) +
        geom_point(size = 2) +
       geom_line() +
      labs(subtitle = "A. Number of target zone entrances") +
       scale_y_continuous(name= "Counts") +
        scale_x_continuous(name= NULL, 
                           breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                           labels = NULL) +
      theme_cowplot(font_size = 7, line_size = 0.25) +
      background_grid(major = "y", minor = "y") +
      scale_color_manual(values = colorvalAPA00)  +
      theme(legend.title=element_blank(),
            legend.position="none",
            legend.text=element_text(size=7)) 
    numentr

![](../figures/01_behavior/twobehaviors-1.png)

    pdf(file="../figures/01_behavior/numentr.pdf", width=2.55, height=2)
    plot(numentr)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    # time spent in target
    target <- dplyr::summarise(group_by(behavior, APA2, TrainSessionComboNum), m = mean(pTimeTarget), se = sd(pTimeTarget)/sqrt(length(pTimeTarget)))
    target$APA2 <- factor(target$APA2, levels = c("yoked-consistent" ,"consistent", "yoked-conflict", "conflict"))

    tagetplot <- ggplot(target, 
                      aes(x=, TrainSessionComboNum, y=m, color=APA2)) + 
        geom_errorbar(aes(ymin=m-se, ymax=m+se, color=APA2), width=.1) +
        geom_point(size = 2) +
       geom_line() +
      labs(subtitle = "B. Time spent in the target zone") +
       scale_y_continuous(name= "Proportion",
                          breaks = c(0, .12, .25, .37)) +
        scale_x_continuous(name= NULL, 
                           breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                           labels = NULL) +
      theme_cowplot(font_size = 7, line_size = 0.25) +
      background_grid(major = "y", minor = "y") +
      scale_color_manual(values = colorvalAPA00)  +
      theme(legend.title=element_blank(),
            legend.position="none",
            legend.text=element_text(size=7)) 
    tagetplot

![](../figures/01_behavior/twobehaviors-2.png)

    pdf(file="../figures/01_behavior/tagetplot.pdf", width=2.55, height=2)
    plot(tagetplot)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    ## time first entrance
    Time1Entr <- dplyr::summarise(group_by(behavior, APA2, TrainSessionComboNum), m = mean(Time1stEntr), se = sd(Time1stEntr)/sqrt(length(Time1stEntr)))
    Time1Entr$APA2 <- factor(Time1Entr$APA2, levels = c("yoked-consistent" ,"consistent", "yoked-conflict", "conflict"))


    timeentr <- ggplot(Time1Entr, 
                      aes(x=, TrainSessionComboNum, y=m, color=APA2)) + 
        geom_errorbar(aes(ymin=m-se, ymax=m+se, color=APA2), width=.1) +
        geom_point(size = 2) +
       geom_line() +
      labs(subtitle = "C. Time to first entrance into the target zone") +
       scale_y_continuous(name= "Time (s)") +
        scale_x_continuous(name= "Training session", 
                           breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                           labels = c( "P", "T1", "T2", "T3",
                                       "Rt", "T4", "T5", "T6", "Rn")) +
      theme_cowplot(font_size = 7, line_size = 0.25) +
      background_grid(major = "y", minor = "y") +
      scale_color_manual(values = colorvalAPA00)  +
      theme(legend.title=element_blank(),
            legend.position="none", 
            legend.text=element_text(size=5)) 

    timeentr

![](../figures/01_behavior/twobehaviors-3.png)

    pdf(file="../figures/01_behavior/time1entr.pdf", width=2.55, height=2.25)
    plot(timeentr)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    ### Propostion of time in opposite zone
    timeopp <- dplyr::summarise(group_by(behavior, APA2, TrainSessionComboNum), m = mean(pTimeOPP), se = sd(pTimeOPP)/sqrt(length(pTimeOPP)))
    timeopp$APA2 <- factor(timeopp$APA2, levels = c("yoked-consistent" ,"consistent", "yoked-conflict", "conflict"))


    timeoppplot <- ggplot(timeopp, 
                      aes(x=, TrainSessionComboNum, y=m, color=APA2)) + 
        geom_errorbar(aes(ymin=m-se, ymax=m+se, color=APA2), width=.1) +
        geom_point(size = 2) +
       geom_line() +
      labs(subtitle = "D. Time spent opposite the target zone") +
       scale_y_continuous(name= "Proportion",
                         breaks = c(0, .25, .5, .75),
                         limits = c(0,.75)) +
        scale_x_continuous(name= "Training session", 
                           breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                           labels = c( "P", "T1", "T2", "T3",
                                       "Rt", "T4", "T5", "T6", "Rn")) +
      theme_cowplot(font_size = 7, line_size = 0.25) +
      background_grid(major = "y", minor = "y") +
      scale_color_manual(values = colorvalAPA00)  +
      theme(legend.title=element_blank(),
            legend.position="none", 
            legend.text=element_text(size=5)) 
    timeoppplot

![](../figures/01_behavior/twobehaviors-4.png)

    pdf(file="../figures/01_behavior/timeoppplot.pdf", width=2.55, height=2.25)
    plot(timeoppplot)
    dev.off()

    ## quartz_off_screen 
    ##                 2

Proportion of time spent
------------------------

Now let’s look at how much time they spend in each proportion

    proptime <- behavior_slim %>%
      dplyr::select(APA2, TrainSessionCombo, pTimeOPP, pTimeTarget, pTimeCCW, pTimeCW) %>%
      tidyr::gather(quadrant, proportion, -APA2, -TrainSessionCombo)

    proptime$quadrant <- as.factor(proptime$quadrant)
    proptime$quadrant <- factor(proptime$quadrant, levels = c("pTimeTarget", "pTimeCCW", "pTimeOPP", "pTimeCW"))
    proptime$TrainSessionCombo <- factor(proptime$TrainSessionCombo, 
                                         levels = c("Hab" ,"T1", "T2", "T3", "Retest", 
                                                    "T4_C1", "T5_C2", "T6_C3","Retention"))

    timespent1 <- proptime %>%
      ggplot(aes(x = TrainSessionCombo, y = proportion, fill = quadrant)) + 
      geom_bar(position = "fill",stat = "identity")  +
      facet_wrap(~APA2, nrow=2) +
      theme_cowplot(font_size = 7, line_size = 0.25) +  
      theme(legend.title=element_blank(),
            legend.position="bottom",
            legend.text=element_text(size=7),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold")) +
      scale_y_continuous(name= "Proportion",
                         breaks = c(0.25,0.50, 0.75)) +
      scale_x_discrete(name="Training Session", 
                         labels = c( "P", "T1", "T2", "T3",
                                       "Rt", "C1", "C2","C3", 
                                       "Rn")) +
      scale_fill_manual(values = c("#f03b20", "#e5f5e0" ,"#a1d99b", "#31a354"),
                         labels=c("Target zone", "Counter-clockwise zone", "Opposite zone", "Clockwise zone")) + 
      geom_hline(yintercept=c(0.25,0.50, 0.75), color="black" , 
                 linetype="dashed",line_size = 0.25) +
      labs(subtitle = "Propotion of time spent in each quadrant") 

    ## Warning: Ignoring unknown parameters: line_size

    timespent1

![](../figures/01_behavior/proptime-1.png)

    pdf(file="../figures/01_behavior/timespent1.pdf",  width=4.25, height=4)
    plot(timespent1)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    meanproplong <- proptime %>%
      filter(TrainSessionCombo != "Hab") %>%
      group_by(APA2, quadrant) %>%
      summarize(mean_prop = mean(proportion)) 
    meanproplong$mean_prop <- round(meanproplong$mean_prop,2)
    head(meanproplong)

    ## # A tibble: 6 x 3
    ## # Groups:   APA2 [2]
    ##   APA2             quadrant    mean_prop
    ##   <fct>            <fct>           <dbl>
    ## 1 yoked-consistent pTimeTarget      0.27
    ## 2 yoked-consistent pTimeCCW         0.27
    ## 3 yoked-consistent pTimeOPP         0.24
    ## 4 yoked-consistent pTimeCW          0.22
    ## 5 consistent       pTimeTarget      0.01
    ## 6 consistent       pTimeCCW         0.36

    summary(aov(data = meanproplong, mean_prop ~ APA2*quadrant))

    ##               Df Sum Sq Mean Sq
    ## APA2           3 0.0000 0.00000
    ## quadrant       3 0.1384 0.04612
    ## APA2:quadrant  9 0.1447 0.01607

    meanpropwide <- spread(meanproplong, quadrant, mean_prop)
    kable(meanpropwide)

<table>
<thead>
<tr>
<th style="text-align:left;">
APA2
</th>
<th style="text-align:right;">
pTimeTarget
</th>
<th style="text-align:right;">
pTimeCCW
</th>
<th style="text-align:right;">
pTimeOPP
</th>
<th style="text-align:right;">
pTimeCW
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
yoked-consistent
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.22
</td>
</tr>
<tr>
<td style="text-align:left;">
consistent
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
yoked-conflict
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.27
</td>
</tr>
<tr>
<td style="text-align:left;">
conflict
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.16
</td>
</tr>
</tbody>
</table>

Hierarchical clusering of time series behavioral data
-----------------------------------------------------

Here I use heirarhical cluster to identify patterns in the behavioral
data. On the y axis see three distinct clusters of behaviors that are 1)
higher in trained animals, 2) higher in yoked animals, and 3) measures
of speed.

    ## create scaled data frame
    behavior_slim_heat <- behavior_slim %>%
      filter(TrainSessionCombo != "Hab")

    behavior_slim_heat$RayleigAngle <- NULL
    behavior_slim_heat$PolarMinBin <- NULL
    scaledaveragedata <- as.data.frame(makescaledaveragedata(behavior_slim_heat))

    ## make annotation df and ann_colors for pheatmap
    ann_cols <- as.data.frame(makecolumnannotations(scaledaveragedata))
    ann_colors = ann_colors_APA2

    # set color breaks
    paletteLength <- 30
    myBreaks <- c(seq(min(scaledaveragedata), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(scaledaveragedata)/paletteLength, max(scaledaveragedata), length.out=floor(paletteLength/2)))

    ## pheatmap for markdown
    pheatmap(scaledaveragedata, show_colnames=T, show_rownames = T,
             annotation_col=ann_cols, 
             annotation_colors = ann_colors,
             treeheight_row = 0, treeheight_col = 50,
             border_color = "grey60" ,
             color = viridis(30),
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation" ,
             clustering_distance_rows = "correlation"
             )

![](../figures/01_behavior/pheatmap2-1.png)

    # pheatmapfor adobe
    pheatmap(scaledaveragedata, show_colnames=F, show_rownames = F,
             annotation_col=ann_cols, annotation_colors = ann_colors,
             annotation_names_col = F,
             treeheight_row = 0, treeheight_col = 15,
             fontsize = 6, 
             border_color = "grey60" ,
             color = viridis(30),
              width = 3, height = 2,
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation",
             filename = "../figures/01_behavior/pheatmap2.pdf",
             legend = TRUE,
             annotation_legend = FALSE
             )

### Principle component analysis

Next, I next reduced the dimentionality of the data with a PCA anlaysis.

    longdata <- makelongdata(behavior)
    Z <- longdata[,3:362]
    Z <- Z[,apply(Z, 2, var, na.rm=TRUE) != 0]
    pc = prcomp(Z, scale.=TRUE)
    loadings <- pc$rotation
    scores <- pc$x

    scoresdf <- makepcadf(behavior) #create the df of pcas
    rotationdf <- mkrotationdf(behavior) #loadings for specific factors
    behaviormatrix <- behavior[c(20:58)]  # for 2nd pca analysis
    scoresdf$PC1 <- scoresdf$PC1 * -1
    scoresdf$APA2 <- factor(scoresdf$APA2, levels = c("yoked-consistent" ,"consistent", "yoked-conflict", "conflict"))


    ## data wraningly for pca anlysis
    behaviormatrix %>% 
      scale() %>%                 # scale to 0 mean and unit variance
      prcomp() ->                 # do PCA
      pca                         # store result as `pca`
    percent <- round(100*pca$sdev^2/sum(pca$sdev^2),2)
    perc_data <- data.frame(percent=percent, PC=1:length(percent))
    res.pca <- prcomp(behaviormatrix,  scale = TRUE)

    # plot of percent contribution
    ggplot(perc_data, aes(x=PC, y=percent)) + 
      geom_bar(stat="identity") + 
      geom_text(aes(label=round(percent, 2)), size=4, vjust=-.5) + 
      xlim(0, 10)

    ## Warning: Removed 29 rows containing missing values (position_stack).

    ## Warning: Removed 29 rows containing missing values (geom_text).

![](../figures/01_behavior/PCA-1.png)

    ## print anova and TukeyHSD stats for first 6 PCs
    j <- 0
    for (i in (scoresdf[,c(1:6)])){
      j <- j+1
      print(paste("PC", j, sep = " "))
      myaov <- aov(i ~ APA2, data=scoresdf)
      print(summary(myaov))
      print(TukeyHSD(myaov, which = "APA2"))
    }

    ## [1] "PC 1"
    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## APA2         3   3195  1065.0   73.05 6.84e-14 ***
    ## Residuals   30    437    14.6                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = i ~ APA2, data = scoresdf)
    ## 
    ## $APA2
    ##                                        diff        lwr        upr
    ## consistent-yoked-consistent      20.9070831  15.716020  26.098146
    ## yoked-conflict-yoked-consistent   2.8124017  -2.232405   7.857208
    ## conflict-yoked-consistent        20.6555262  15.610720  25.700333
    ## yoked-conflict-consistent       -18.0946814 -23.139488 -13.049875
    ## conflict-consistent              -0.2515569  -5.296364   4.793250
    ## conflict-yoked-conflict          17.8431245  12.948943  22.737306
    ##                                     p adj
    ## consistent-yoked-consistent     0.0000000
    ## yoked-conflict-yoked-consistent 0.4409343
    ## conflict-yoked-consistent       0.0000000
    ## yoked-conflict-consistent       0.0000000
    ## conflict-consistent             0.9990895
    ## conflict-yoked-conflict         0.0000000
    ## 
    ## [1] "PC 2"
    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## APA2         3  171.6   57.20   3.091 0.0419 *
    ## Residuals   30  555.2   18.51                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = i ~ APA2, data = scoresdf)
    ## 
    ## $APA2
    ##                                       diff        lwr       upr     p adj
    ## consistent-yoked-consistent     -0.8999943 -6.7486749  4.948686 0.9748977
    ## yoked-conflict-yoked-consistent  4.4517329 -1.2321630 10.135629 0.1669142
    ## conflict-yoked-consistent        3.4421998 -2.2416961  9.126096 0.3688648
    ## yoked-conflict-consistent        5.3517272 -0.3321688 11.035623 0.0707214
    ## conflict-consistent              4.3421941 -1.3417018 10.026090 0.1836316
    ## conflict-yoked-conflict         -1.0095331 -6.5237221  4.504656 0.9589440
    ## 
    ## [1] "PC 3"
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## APA2         3   94.9   31.64   2.001  0.135
    ## Residuals   30  474.3   15.81               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = i ~ APA2, data = scoresdf)
    ## 
    ## $APA2
    ##                                      diff        lwr      upr     p adj
    ## consistent-yoked-consistent     -2.477541 -7.8832385 2.928157 0.6031763
    ## yoked-conflict-yoked-consistent -1.147686 -6.4010795 4.105708 0.9330967
    ## conflict-yoked-consistent        2.053796 -3.1995983 7.307190 0.7141711
    ## yoked-conflict-consistent        1.329855 -3.9235387 6.583249 0.9007412
    ## conflict-consistent              4.531336 -0.7220575 9.784730 0.1102806
    ## conflict-yoked-conflict          3.201481 -1.8950595 8.298022 0.3373455
    ## 
    ## [1] "PC 4"
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## APA2         3   60.0   20.00   1.305  0.291
    ## Residuals   30  459.9   15.33               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = i ~ APA2, data = scoresdf)
    ## 
    ## $APA2
    ##                                       diff       lwr      upr     p adj
    ## consistent-yoked-consistent      1.3080022 -4.015060 6.631065 0.9082391
    ## yoked-conflict-yoked-consistent -2.3451756 -7.518263 2.827911 0.6115099
    ## conflict-yoked-consistent       -0.8791771 -6.052264 4.293910 0.9667025
    ## yoked-conflict-consistent       -3.6531778 -8.826265 1.519909 0.2412559
    ## conflict-consistent             -2.1871792 -7.360266 2.985908 0.6623510
    ## conflict-yoked-conflict          1.4659986 -3.552633 6.484630 0.8564632
    ## 
    ## [1] "PC 5"
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## APA2         3   10.5   3.496   0.233  0.873
    ## Residuals   30  450.8  15.027               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = i ~ APA2, data = scoresdf)
    ## 
    ## $APA2
    ##                                       diff       lwr      upr     p adj
    ## consistent-yoked-consistent     -0.6565012 -5.926822 4.613819 0.9863545
    ## yoked-conflict-yoked-consistent -1.4904797 -6.612311 3.631351 0.8578281
    ## conflict-yoked-consistent       -0.3613765 -5.483207 4.760454 0.9974400
    ## yoked-conflict-consistent       -0.8339785 -5.955809 4.287852 0.9705093
    ## conflict-consistent              0.2951246 -4.826706 5.416956 0.9985986
    ## conflict-yoked-conflict          1.1291031 -3.839803 6.098009 0.9255676
    ## 
    ## [1] "PC 6"
    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## APA2         3  123.5   41.18    4.13 0.0145 *
    ## Residuals   30  299.1    9.97                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = i ~ APA2, data = scoresdf)
    ## 
    ## $APA2
    ##                                       diff        lwr       upr     p adj
    ## consistent-yoked-consistent     -3.6169924 -7.9100582 0.6760735 0.1228341
    ## yoked-conflict-yoked-consistent -3.2507693 -7.4228793 0.9213408 0.1703278
    ## conflict-yoked-consistent        0.6641045 -3.5080056 4.8362146 0.9723522
    ## yoked-conflict-consistent        0.3662231 -3.8058870 4.5383331 0.9951112
    ## conflict-consistent              4.2810969  0.1089868 8.4532069 0.0426272
    ## conflict-yoked-conflict          3.9148738 -0.1326675 7.9624151 0.0608419

    pca12elipse <- ggplot(scoresdf, aes(PC1,PC2, color=APA2)) +
        geom_point(size=2.5, alpha = 0.7) +
        xlab(paste0("PC 1: ", percent[1],"% variance")) +
        ylab(paste0("PC 2: ", percent[2],"% variance")) +
        stat_ellipse(level = 0.95, (aes(color=APA2)),size=0.25) + 
        scale_colour_manual(values=c(colorvalAPA00)) + 
        theme_cowplot(font_size = 7, line_size = 0.25) +
        theme(legend.position="none") +
        ylim(-22,20)
    pca12elipse

![](../figures/01_behavior/PCA-2.png)

    pdf(file="../figures/01_behavior/pca12elipse.pdf",  width=2, height=2)
    plot(pca12elipse)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    # PCA with contributions
    res.pca <- prcomp(behaviormatrix, scale = TRUE)
    fviz_pca_var(res.pca,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE,     # Avoid text overlapping
                 select.var = list(contrib = 10))

![](../figures/01_behavior/PCA-3.png)

### Comparing Consistent and Conflict behaviors during the T4/C1 training session

    filtered <- behavior_slim %>% filter(TrainSessionCombo == "T4_C1", APA != "control") 
    exp_factors <- as.data.frame(filtered[,1])
    exp_nums <- filtered[,c(4:42)]
    exp_factors$APA2 <- factor(filtered$APA2, levels = c("consistent", "conflict"))

    # Levene's test for normality
    #for(y in names(exp_nums)){
    #  ymod <- leveneTest(exp_nums[[y]] ~ exp_factors$APA2)
    #  cat(paste('\nDependent var:', y, '\n'))
    #  print(ymod)
    #}

    # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # *** Speed1, Path2ndEntr, Time2ndEntr, Path2ndEntr, Time2ndEntr
    # ** Path1stEntr, Time1stEntr
    # *  Max50.RngHiBin ,  PolarMaxBin  , PolarMinVal, RayleigAngle, pTimeCW, pTimeOPP,
    # .  Speed2, Min50.RngLoBin , TimeTarget, NumShock, NumEntrances
    #    AnnularKurtosis, AnnularSkewnes, AnnularSd, AnnularMaxBin, AnnularMaxVal, AnnularMinBin, AnnularMinVal, Max50.RngLoBin, RayleigLength PolarMaxVal, Min50.RngHiBin , PolarMinBin, PolarSdVal, PolarAvgVal, RayleigLength, pTimeTarget, Speed2ndEntr, MaxTimeAvoid, Dist1stEntr.m, Speed1stEntr.cm.s, Linearity.Arena, SdevSpeedArena

    for(y in names(exp_nums)){
      ymod <- wilcox.test(exp_nums[[y]] ~ exp_factors$APA2 )
      cat(paste('\nDependent var:', y, '\n'))
      print(ymod)
    }

    ## Warning in wilcox.test.default(x = c(3.21, 2.82, 2.2, 2.64, 2.88, 2.5, 3, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: SdevSpeedArena 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 33, p-value = 0.8098
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(0.3341, 0.2786, 0.3786, 0.3168,
    ## 0.3158, : cannot compute exact p-value with ties

    ## 
    ## Dependent var: Linearity.Arena. 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 52.5, p-value = 0.1234
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(3L, 10L, 0L, 0L, 2L, 4L, 13L, 2L), :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: NumEntrances 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 4.5, p-value = 0.002804
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(257.83, 39.7, 599.97, 599.97, 249.4, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: Time1stEntr 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 72, p-value = 0.0006258
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(8.59, 0.92, 14.26, 16.08, 7.75,
    ## 3.19, : cannot compute exact p-value with ties

    ## 
    ## Dependent var: Path1stEntr 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 72, p-value = 0.0005879
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(9.54, 2.04, -1, -1, 13.37, 2.38,
    ## 1.96, : cannot compute exact p-value with ties

    ## 
    ## Dependent var: Speed1stEntr.cm.s. 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 48, p-value = 0.258
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(0.15, 0.57, 0, 0, 0.11, 0.24, 0.7, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: Dist1stEntr.m. 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 6, p-value = 0.004506
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(3L, 10L, 0L, 0L, 3L, 5L, 13L, 2L), :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: NumShock 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 2.5, p-value = 0.001456
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(257L, 195L, 599L, 599L, 293L, 161L, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: MaxTimeAvoid 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 49.5, p-value = 0.2104
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(321.23, 55, 599.97, 599.97, 543.07, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: Time2ndEntr 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 72, p-value = 0.0006306
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: Path2ndEntr 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 72, p-value = 8.227e-05
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(1.61, 1.97, -1, -1, 1.82, 2.47,
    ## 2.17, : cannot compute exact p-value with ties

    ## 
    ## Dependent var: Speed2ndEntr 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 23.5, p-value = 0.2473
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(4.5, 8.064, 0, 0, 4.601, 4.599,
    ## 10.632, : cannot compute exact p-value with ties

    ## 
    ## Dependent var: TimeTarget 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 1, p-value = 0.0008944
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(0.0107, 0.0203, 0, 0, 0.0116,
    ## 0.0124, : cannot compute exact p-value with ties

    ## 
    ## Dependent var: pTimeTarget 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 3, p-value = 0.001753
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: pTimeCCW 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 37, p-value = 0.9626
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: pTimeOPP 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 64, p-value = 0.005512
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: pTimeCW 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 20, p-value = 0.1388
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(0.52, 0.61, 0.86, 0.7, 0.68, 0.67, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: RayleigLength 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 40.5, p-value = 0.6998
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: RayleigAngle 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 38, p-value = 0.8884
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: PolarAvgVal 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 57, p-value = 0.0464
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: PolarSdVal 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 11, p-value = 0.01522
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(1e-04, 1e-04, 0, 0, 0, 0, 0, 0), y =
    ## c(5e-04, : cannot compute exact p-value with ties

    ## 
    ## Dependent var: PolarMinVal 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 10, p-value = 0.00948
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(10L, 320L, 0L, 0L, 0L, 310L, 0L, 0L:
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: PolarMinBin 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 22, p-value = 0.1828
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(250L, 160L, 190L, 200L, 150L, 160L, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: Min50.RngLoBin 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 36.5, p-value = 1
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(180L, 100L, 170L, 150L, 110L, 110L, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: Min50.RngHiBin 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 44, p-value = 0.4688
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: PolarMaxVal 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 30, p-value = 0.6058
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(200L, 120L, 170L, 180L, 120L, 130L, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: PolarMaxBin 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 34.5, p-value = 0.9231
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(160L, 90L, 150L, 140L, 100L, 100L, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: Max50.RngLoBin 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 40.5, p-value = 0.6998
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(260L, 180L, 200L, 220L, 170L, 180L, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: Max50.RngHiBin 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 35, p-value = 0.9615
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(6e-04, 0.0096, 0.0022, 0.005, 9e-04, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: AnnularMinVal 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 40, p-value = 0.736
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(3.5, 3.5, 8.5, 8.5, 19.4, 19.4, 3.5, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: AnnularMinBin 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 34, p-value = 0.8797
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: AnnularMaxVal 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 24, p-value = 0.2766
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(18, 16.6, 16.6, 16.6, 16.6, 15, 18, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: AnnularMaxBin 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 40.5, p-value = 0.6614
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(16.18, 15.81, 16.7, 16.39, 16.18,
    ## 14.94, : cannot compute exact p-value with ties

    ## 
    ## Dependent var: AnnularAvg 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 41, p-value = 0.6648
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: AnnularSd 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 34, p-value = 0.8884
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: AnnularSkewnes 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 34, p-value = 0.8884
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: AnnularKurtosis 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 33, p-value = 0.8148
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(0.03331652639336,
    ## 0.0231738035264484, : cannot compute exact p-value with ties

    ## 
    ## Dependent var: Speed1 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 48, p-value = 0.2655
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: Speed2 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 6, p-value = 0.002468
    ## alternative hypothesis: true location shift is not equal to 0

    # *** Path2ndEntr 
    # **  Speed2, PolarMinVal, pTimeOPP, pTimeTarget, TimeTarget, Time2ndEntr, NumShock
    # **  Dist1stEntr.m., Path1stEntr , Time1stEntr, NumEntrances
    # *   PolarSdVal, PolarAvgVal 
    # .    
    #     Speed1, AnnularKurtosis, AnnularSkewnes, AnnularSd, AnnularAvg, AnnularMaxBin,
    #     AnnularMaxVal, AnnularMinBin, AnnularMinVal, Max50.RngHiBin , PolarMaxBin, 
    #     PolarMaxVal, Min50.RngHiBin, Min50.RngLoBin, PolarMinBin, RayleigAngle
    #     RayleigLength, pTimeCW, pTimeCCW, Speed2ndEntr, MaxTimeAvoid, Speed1stEntr.cm.s., #     Linearity.Arena., SdevSpeedArena 
     
    par(mfrow=c(3,3))
    for(y in names(exp_nums)){
      ymod <- boxplot(exp_nums[[y]] ~ exp_factors$APA2,
                   main = y,
                   xlab = "T4/C1")
    }

![](../figures/01_behavior/T4consistentconflict-1.png)![](../figures/01_behavior/T4consistentconflict-2.png)![](../figures/01_behavior/T4consistentconflict-3.png)![](../figures/01_behavior/T4consistentconflict-4.png)![](../figures/01_behavior/T4consistentconflict-5.png)

### Comparing Consistent and Conflict behaviors during the T6/C3 training session

    filtered <- behavior_slim %>% filter(TrainSessionCombo == "T6_C3", APA != "control") 
    exp_factors <- as.data.frame(filtered[,1])
    exp_nums <- filtered[,c(4:42)]
    exp_factors$APA2 <- factor(filtered$APA2, levels = c("consistent", "conflict"))

    for(y in names(exp_nums)){
      ymod<- wilcox.test(exp_nums[[y]] ~ exp_factors$APA2 )
      cat(paste('\nDependent var:', y, '\n'))
      print(ymod)
    }

    ## Warning in wilcox.test.default(x = c(2.37, 3.08, 1.96, 2.32, 2.33, 2.55, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: SdevSpeedArena 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 21.5, p-value = 0.1774
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: Linearity.Arena. 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 36, p-value = 1
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(0L, 12L, 0L, 0L, 1L, 13L, 8L, 0L), :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: NumEntrances 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 23.5, p-value = 0.2421
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(599.97, 3.13, 599.97, 599.97, 46.53, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: Time1stEntr 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 42, p-value = 0.5944
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: Path1stEntr 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 41, p-value = 0.673
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(-1, 2.1, -1, -1, 1.47, 1.96, 2.17, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: Speed1stEntr.cm.s. 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 11, p-value = 0.01769
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(0, 0.68, 0, 0, 0.07, 0.84, 0.43, 0), :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: Dist1stEntr.m. 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 25, p-value = 0.3081
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(0L, 12L, 0L, 0L, 1L, 13L, 8L, 0L), :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: NumShock 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 23.5, p-value = 0.2421
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(599L, 113L, 599L, 599L, 553L, 86L, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: MaxTimeAvoid 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 47, p-value = 0.3093
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(599.97, 28.37, 599.97, 599.97,
    ## 599.97, : cannot compute exact p-value with ties

    ## 
    ## Dependent var: Time2ndEntr 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 44.5, p-value = 0.4163
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: Path2ndEntr 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 31, p-value = 0.673
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(-1, 2.17, -1, -1, -1, 1.56, 1.95, -1:
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: Speed2ndEntr 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 26.5, p-value = 0.3605
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(0, 6, 0, 0, 0.667, 9.401, 7.833, 0), :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: TimeTarget 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 21.5, p-value = 0.175
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(0, 0.016, 0, 0, 0.002, 0.0258,
    ## 0.0201, : cannot compute exact p-value with ties

    ## 
    ## Dependent var: pTimeTarget 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 27, p-value = 0.4102
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: pTimeCCW 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 41, p-value = 0.673
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: pTimeOPP 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 31, p-value = 0.673
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(0.3233, 0.0443, 0.0213, 0.3027, 0, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: pTimeCW 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 26, p-value = 0.3603
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(0.78, 0.67, 0.87, 0.68, 0.92, 0.72, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: RayleigLength 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 42.5, p-value = 0.5635
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: RayleigAngle 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 35, p-value = 0.9626
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: PolarAvgVal 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 68, p-value = 0.0009872
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: PolarSdVal 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 18, p-value = 0.09272
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(0, 0, 0, 0, 0, 0, 1e-04, 0), y =
    ## c(1e-04, : cannot compute exact p-value with ties

    ## 
    ## Dependent var: PolarMinVal 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 24, p-value = 0.1657
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(0L, 310L, 0L, 0L, 0L, 0L, 340L, 0L), :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: PolarMinBin 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 31, p-value = 0.6121
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(230L, 140L, 190L, 220L, 140L, 120L, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: Min50.RngLoBin 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 34.5, p-value = 0.9229
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(200L, 90L, 160L, 170L, 120L, 70L,
    ## 70L, : cannot compute exact p-value with ties

    ## 
    ## Dependent var: Min50.RngHiBin 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 37.5, p-value = 0.9226
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: PolarMaxVal 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 37, p-value = 0.9626
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(210L, 140L, 190L, 200L, 120L, 70L, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: PolarMaxBin 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 35, p-value = 0.9615
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(180L, 70L, 140L, 160L, 100L, 60L,
    ## 60L, : cannot compute exact p-value with ties

    ## 
    ## Dependent var: Max50.RngLoBin 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 37, p-value = 0.9612
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(240L, 150L, 200L, 240L, 150L, 140L, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: Max50.RngHiBin 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 35.5, p-value = 1
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(0.009, 0.0034, 0.0059, 0.0218,
    ## 8e-04, : cannot compute exact p-value with ties

    ## 
    ## Dependent var: AnnularMinVal 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 50.5, p-value = 0.1777
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(3.5, 8.5, 8.5, 8.5, 8.5, 3.5, 3.5, :
    ## cannot compute exact p-value with ties

    ## 
    ## Dependent var: AnnularMinBin 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 23, p-value = 0.185
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: AnnularMaxVal 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 34, p-value = 0.8884
    ## alternative hypothesis: true location shift is not equal to 0

    ## Warning in wilcox.test.default(x = c(16.6, 16.6, 16.6, 16.6, 16.6, 18,
    ## 18, : cannot compute exact p-value with ties

    ## 
    ## Dependent var: AnnularMaxBin 
    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 33, p-value = 0.7609
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: AnnularAvg 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 26, p-value = 0.3704
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: AnnularSd 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 41, p-value = 0.673
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: AnnularSkewnes 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 32, p-value = 0.743
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: AnnularKurtosis 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 31, p-value = 0.673
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: Speed1 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 20, p-value = 0.1388
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ## Dependent var: Speed2 
    ## 
    ##  Wilcoxon rank sum test
    ## 
    ## data:  exp_nums[[y]] by exp_factors$APA2
    ## W = 22, p-value = 0.1996
    ## alternative hypothesis: true location shift is not equal to 0

    # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # ***  
    # **  PolarAvgVal,
    # *   Speed1stEntr.cm.s. 
    # .   PolarSdVal,   
    #     Speed2, Speed1, AnnularKurtosis, AnnularSkewnes, AnnularSd,
    #     AnnularAvg, AnnularMaxBin, AnnularMaxVal, AnnularMinBin, AnnularMinVal
    #     Max50.RngHiBin, Max50.RngLoBin, PolarMaxBin, PolarMaxVal, Min50.RngHiBin   
    #     Min50.RngLoBin, PolarMinBin, PolarMinVal, RayleigAngle, RayleigLength,
    #     pTimeCW, pTimeOPP, pTimeCCW, pTimeTarget, TimeTarget, Speed2ndEntr, 
    #     Path2ndEntr, Time2ndEntr, MaxTimeAvoid, NumShock, Dist1stEntr.m. 
    #     Path1stEntr, Time1stEntr, NumEntrances, Linearity.Arena., SdevSpeedArena

    par(mfrow=c(3,3))
    for(y in names(exp_nums)){
      ymod <- boxplot(exp_nums[[y]] ~ exp_factors$APA2,
                   main = y,
                   xlab = "T6/C3")
    }

![](../figures/01_behavior/T6consistentconflict-1.png)![](../figures/01_behavior/T6consistentconflict-2.png)![](../figures/01_behavior/T6consistentconflict-3.png)![](../figures/01_behavior/T6consistentconflict-4.png)![](../figures/01_behavior/T6consistentconflict-5.png)

    write.csv(behavior, file = "../data/01a_behavior.csv", row.names = FALSE)
    write.csv(meannumentr, file = "../data/01a_meannumentr.csv", row.names = FALSE)
    write.csv(Time1Entr, file = "../data/01a_Time1Entr.csv", row.names = FALSE)
    write.csv(scoresdf, file = "../data/01a_scoresdf.csv", row.names = FALSE)
    write.csv(rotationdf, file = "../data/01a_rotationdf.csv", row.names = TRUE)
    write.csv(behaviormatrix, file = "../data/01a_behaviormatrix.csv", row.names = TRUE)
