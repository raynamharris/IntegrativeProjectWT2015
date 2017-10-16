After wrangling the behaivoral data in the previous script ([01a\_beahvior\_create\_dfs.Rmd](./01a_beahvior_create_dfs.Rmd)), conducted the analyses descirbed in this scipt. All plots are generate as **.png** files for markdown viewing and as **.pdf** files for incorporation to Adobe Illustrator.

<img src="../figures/figures-01.png" width="1370" />

``` r
library(ggplot2) ## for awesome plots!
library(cowplot) ## for some easy to use themes
library(plyr)
library(dplyr) ## for filtering and selecting rows
library(factoextra)  ##pca with vectors
library(car) ## stats
library(superheat) # for kmeans clustered heatmap
library(pheatmap)  # for pretty heatmap
library(viridis) # for awesome color pallette
library(reshape2) ## for melting dataframe
library(tidyr) ## for respahing data


## load user-written functions 
source("functions_behavior.R")
source("figureoptions.R")

## set output file for figures 
knitr::opts_chunk$set(fig.path = '../figures/01_behavior/')
```

To help make production of these figures more reproducible, I first import some intermediate data files that I cleaned and manipulated for data vizualization. I also relevel factors here to overide the defalut alphabetical plotting.

``` r
behavior <- read.csv("../data/01a_behavior.csv", header = T)
threeplots <- read.csv("../data/01a_threeplots.csv", header = T)
scoresdf <- read.csv("../data/01a_scoresdf.csv", header = T)
rotationdf <- read.csv("../data/01a_rotationdf.csv", header = T, row.names = 1)
behaviormatrix <- read.csv("../data/01a_behaviormatrix.csv", header = T, row.names = 1)

#set factor levels
scoresdf$APA <- NULL
scoresdf$APA2 <- factor(scoresdf$APA2, levels = c("yoked-consistent" ,"consistent", "yoked-conflict", "conflict"))
threeplots$APA2 <- factor(threeplots$APA2, levels = c("yoked-consistent" ,"consistent", "yoked-conflict", "conflict"))
threeplots$measure <- factor(threeplots$measure, levels = c("Number of Entrances" , "Max Avoidance Time", "Speed"))
```

Figure 1B: Standard vizualization of mean avoidance beavior
-----------------------------------------------------------

First, I visualze the group mean and standard error for the time it takes before an individual mouse enters the spatial region marked "schock zone" or equivilent (Fig. 1A).

``` r
behaviorwrap <- ggplot(threeplots, aes(x=, TrainSessionComboNum, y=m, color=APA2)) + 
    geom_errorbar(aes(ymin=m-se, ymax=m+se, color=APA2), width=.1) +
    geom_point(size = 2) +
   geom_line() +
   scale_y_continuous(name= NULL) +
    scale_x_continuous(name="Training Session", 
                       breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                       labels = c( "Hab.", "T1", "T2", "T3",
                                   "Retest", "T4", "T5", "T6", "Reten.")) +
  theme_cowplot(font_size = 8, line_size = 0.25) +
  #background_grid(major = "y", minor = "y") +
  scale_color_manual(values = colorvalAPA00)  +
  theme(legend.position=c(0.7, 0.8))  +
  theme(legend.title=element_blank()) +
  theme(legend.position="none") +
  facet_wrap(~measure, ncol=1, scales = "free_y")
behaviorwrap
```

![](../figures/01_behavior/behaviorplots-1.png)

``` r
pdf(file="../figures/01_behavior/threebehaviors.pdf", width=3.25, height=3.75)
plot(behaviorwrap)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

Hierarchical clusering of time series behavioral data
-----------------------------------------------------

Here I use heirarhical cluster to identify patterns in the behavioral data. On the y axis see three distinct clusters of behaviors that are 1) higher in trained animals, 2) higher in yoked animals, and 3) measures of speed (Fig. 1C).

``` r
## make annotation df and ann_colors for pheatmap
behavior$RayleigAngle <- NULL
behavior$PolarMinBin <- NULL
scaledaveragedata2 <- as.data.frame(makescaledaveragedata2(behavior))

# create a rubric for the color coding and load the colors from figureoptions.R
df2 <- as.data.frame(makecolumnannotations2(scaledaveragedata2))
ann_colors = ann_colors9  # includes color for session & APA

# set color breaks
paletteLength <- 30
myBreaks <- c(seq(min(scaledaveragedata2), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(scaledaveragedata2)/paletteLength, max(scaledaveragedata2), length.out=floor(paletteLength/2)))

## pheatmap for markdown
pheatmap(scaledaveragedata2, show_colnames=T, show_rownames = T,
         annotation_col=df2, 
         annotation_colors = ann_colors,
         treeheight_row = 0, treeheight_col = 50,
         border_color = "grey60" ,
         color = viridis(30),
         clustering_method="average",
         breaks=myBreaks,
         clustering_distance_cols="correlation" ,
         clustering_distance_rows = "correlation"
         )
```

![](../figures/01_behavior/pheatmap2-1.png)

``` r
# pheatmapfor adobe
pheatmap(scaledaveragedata2, show_colnames=F, show_rownames = F,
         annotation_col=df2, annotation_colors = ann_colors,
         annotation_names_col = F,
         treeheight_row = 0, treeheight_col = 25,
         fontsize = 6, 
         border_color = "grey60" ,
         color = viridis(30),
          width = 3.25, height = 3.5,
         clustering_method="average",
         breaks=myBreaks,
         clustering_distance_cols="correlation",
         filename = "../figures/01_behavior/pheatmap2.pdf",
         legend = TRUE,
         annotation_legend = FALSE
         )
```

### Figure 1D1 and 1D2: Principle component analysis

Next, I next reduced the dimentionality of the data with a PCA anlaysis. PC1 explains 35% of the variation in the data. All other PCs explain less than 10% of the variation.

``` r
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
```

    ## Warning: Removed 29 rows containing missing values (position_stack).

    ## Warning: Removed 29 rows containing missing values (geom_text).

![](../figures/01_behavior/PCA-1.png)

### Figure 1D1: 35% of behaivor variance (PC1) separates yoked from trained

PC1 encompases differences between yoked trained indivdual but does not significantly differ between consistent and conflict trained aniamls. To confirm statistical significance of this visual pattern, we conducted a two-way treatment x region ANOVA and confirmed a significant effect of region (F2,31= 101.39; p = 2.5e-14). Post hoc Tukey tests confirmed conflict = consistent &lt; control). The major contibutors to this variation are number of shocks and distance to first entrance.

``` r
## statistics
aov1 <- aov(PC1 ~ APA2, data=scoresdf)
summary(aov1) # p = 1.01e-13
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## APA2         3   3247  1082.5   70.92 1.01e-13 ***
    ## Residuals   30    458    15.3                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TukeyHSD(aov1, which = "APA2") # p<< 0.001 for both control comparisions
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ APA2, data = scoresdf)
    ## 
    ## $APA2
    ##                                        diff        lwr        upr
    ## consistent-yoked-consistent      21.0488034  15.737228  26.360379
    ## yoked-conflict-yoked-consistent   2.8209039  -2.341020   7.982828
    ## conflict-yoked-consistent        20.8378476  15.675924  25.999771
    ## yoked-conflict-consistent       -18.2278994 -23.389823 -13.065976
    ## conflict-consistent              -0.2109558  -5.372879   4.950968
    ## conflict-yoked-conflict          18.0169437  13.009142  23.024745
    ##                                     p adj
    ## consistent-yoked-consistent     0.0000000
    ## yoked-conflict-yoked-consistent 0.4582117
    ## conflict-yoked-consistent       0.0000000
    ## yoked-conflict-consistent       0.0000000
    ## conflict-consistent             0.9994975
    ## conflict-yoked-conflict         0.0000000

``` r
aov2 <- aov(PC2 ~ APA2, data=scoresdf)
summary(aov2) # p = 0.0295 *
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## APA2         3  188.7   62.89   3.427 0.0295 *
    ## Residuals   30  550.6   18.35                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TukeyHSD(aov2, which = "APA2") # p > 0.05
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ APA2, data = scoresdf)
    ## 
    ## $APA2
    ##                                       diff        lwr       upr     p adj
    ## consistent-yoked-consistent     -1.0965242 -6.9208797  4.727831 0.9556089
    ## yoked-conflict-yoked-consistent  4.4155199 -1.2447361 10.075776 0.1695436
    ## conflict-yoked-consistent        3.7560908 -1.9041652  9.416347 0.2911558
    ## yoked-conflict-consistent        5.5120442 -0.1482119 11.172300 0.0585081
    ## conflict-consistent              4.8526151 -0.8076410 10.512871 0.1134571
    ## conflict-yoked-conflict         -0.6594291 -6.1506841  4.831826 0.9877336

``` r
summary(aov(PC3 ~ APA2, data=scoresdf)) # p = 0.117
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## APA2         3  102.5   34.17   2.135  0.117
    ## Residuals   30  480.0   16.00

``` r
summary(aov(PC3 ~ APA2, data=scoresdf))
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## APA2         3  102.5   34.17   2.135  0.117
    ## Residuals   30  480.0   16.00

``` r
summary(aov(PC4 ~ APA2, data=scoresdf))
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## APA2         3   61.9   20.63   1.337  0.281
    ## Residuals   30  463.0   15.44

``` r
summary(aov(PC5 ~ APA2, data=scoresdf))
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## APA2         3   10.5   3.484   0.227  0.877
    ## Residuals   30  461.4  15.381

``` r
summary(aov(PC6 ~ APA2, data=scoresdf))
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## APA2         3  136.3   45.44   4.447 0.0106 *
    ## Residuals   30  306.6   10.22                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TukeyHSD((aov(PC6 ~ APA2, data=scoresdf)), which = "APA2") # p > 0.05
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC6 ~ APA2, data = scoresdf)
    ## 
    ## $APA2
    ##                                       diff         lwr       upr     p adj
    ## consistent-yoked-consistent     -3.7852251 -8.13141662 0.5609664 0.1053573
    ## yoked-conflict-yoked-consistent -3.0836300 -7.30736900 1.1401089 0.2160376
    ## conflict-yoked-consistent        0.9506748 -3.27306411 5.1744138 0.9274591
    ## yoked-conflict-consistent        0.7015950 -3.52214391 4.9253340 0.9687918
    ## conflict-consistent              4.7358999  0.51216098 8.9596389 0.0233390
    ## conflict-yoked-conflict          4.0343049 -0.06332378 8.1319336 0.0548778

``` r
scoresdf$wrap <- "Principal Component Analyses of Behavior"

pca12 <- ggplot(scoresdf, aes(PC1,PC2, color=APA2)) +
    geom_point(size=3, alpha = 0.7) +
    xlab(paste0("PC 1: ", percent[1],"% variance")) +
    ylab(paste0("PC 2: ", percent[2],"% variance")) +
    scale_colour_manual(values=c(colorvalAPA00)) + 
    theme_cowplot(font_size = 8, line_size = 0.25) +
      theme(legend.position="none") +
  facet_wrap(~wrap)
pca12
```

![](../figures/01_behavior/PCAplots-1.png)

``` r
pdf(file="../figures/01_behavior/pca12.pdf",  width=1.5, height=1.5)
plot(pca12)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
pca16 <- ggplot(scoresdf, aes(PC1,PC6, color=APA2)) +
    geom_point(size=3, alpha = 0.7) +
    xlab(paste0("PC 1: ", percent[1],"% variance")) +
    ylab(paste0("PC 6: ", percent[6],"% variance")) +
    stat_ellipse(level = 0.95, (aes(color=APA2)),size=0.25) + 
    scale_colour_manual(values=c(colorvalAPA00)) + 
    theme_cowplot(font_size = 8, line_size = 0.25) +
      theme(legend.position="none") +
  facet_wrap(~wrap)
pca16
```

![](../figures/01_behavior/PCAplots-2.png)

``` r
pdf(file="../figures/01_behavior/pca16.pdf",  width=3.25, height=2.25)
plot(pca16)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

Here are some stats modeling combinatorail PCs. I'm really not sure if this makes sense. I should probabaly model some behavior time interaction....

``` r
lm1 <- lm(PC1~APA2, data=scoresdf)
summary(lm1)
```

    ## 
    ## Call:
    ## lm(formula = PC1 ~ APA2, data = scoresdf)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -7.026 -2.518  0.098  2.252 10.180 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         -11.215      1.381  -8.119 4.60e-09 ***
    ## APA2consistent       21.049      1.953  10.775 7.82e-12 ***
    ## APA2yoked-conflict    2.821      1.898   1.486    0.148    
    ## APA2conflict         20.838      1.898  10.977 5.00e-12 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 3.907 on 30 degrees of freedom
    ## Multiple R-squared:  0.8764, Adjusted R-squared:  0.8641 
    ## F-statistic: 70.92 on 3 and 30 DF,  p-value: 1.009e-13

``` r
lm16 <- lm(PC1+PC6~APA2, data=scoresdf)
summary(lm16)
```

    ## 
    ## Call:
    ## lm(formula = PC1 + PC6 ~ APA2, data = scoresdf)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -11.7560  -3.6630   0.7667   2.6693  10.6940 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         -9.7600     1.7434  -5.598 4.30e-06 ***
    ## APA2consistent      17.2636     2.4655   7.002 8.82e-08 ***
    ## APA2yoked-conflict  -0.2627     2.3960  -0.110    0.913    
    ## APA2conflict        21.7885     2.3960   9.094 3.99e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 4.931 on 30 degrees of freedom
    ## Multiple R-squared:  0.8242, Adjusted R-squared:  0.8066 
    ## F-statistic: 46.87 on 3 and 30 DF,  p-value: 1.945e-11

``` r
lm136 <- lm(PC1+PC3+PC6~APA2, data=scoresdf)
summary(lm136)
```

    ## 
    ## Call:
    ## lm(formula = PC1 + PC3 + PC6 ~ APA2, data = scoresdf)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -10.637  -3.346  -1.297   3.727  14.601 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          -9.379      1.838  -5.103 1.74e-05 ***
    ## APA2consistent       14.746      2.599   5.674 3.48e-06 ***
    ## APA2yoked-conflict   -1.594      2.526  -0.631    0.533    
    ## APA2conflict         23.917      2.526   9.469 1.60e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 5.198 on 30 degrees of freedom
    ## Multiple R-squared:  0.8287, Adjusted R-squared:  0.8115 
    ## F-statistic: 48.36 on 3 and 30 DF,  p-value: 1.322e-11

getting file name for a time spent heatmap
------------------------------------------

``` r
y <- behavior %>%
  distinct(ID, APA2)
x <- read.csv("~/Github/BehavEphyRNAseq/data/behavior/APA_2013-2016.csv", header = T)
x <- x %>%
  distinct(ID, filename)

z <- inner_join(x,y)
```

    ## Joining, by = "ID"

    ## Warning in inner_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
    ## factors with different levels, coercing to character vector

``` r
write.csv(z, "filnames.csv", row.names = F)
```
