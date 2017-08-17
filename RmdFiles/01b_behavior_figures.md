After wrangling the behaivoral data in the previous script ([01a\_beahvior\_create\_dfs.Rmd](./01a_beahvior_create_dfs.Rmd)), conducted the analyses descirbed in this scipt. All plots are generate as **.png** files for markdown viewing and as **.pdf** files for incorporation to Adobe Illustrator.

<img src="../figures/01_behavior/01_avoidancebehvaior-01.png" width="1370" />

``` r
library(ggplot2) ## for awesome plots!
library(cowplot) ## for some easy to use themes
library(dplyr) ## for filtering and selecting rows
library(factoextra)  ##pca with vectors
library(car) ## stats
library(superheat) # for kmeans clustered heatmap
library(pheatmap)  # for pretty heatmap
library(viridis) # for awesome color pallette

## load user-written functions 
source("functions_behavior.R")
source("figureoptions.R")

## set output file for figures 
knitr::opts_chunk$set(fig.path = '../figures/01_behavior/')
```

To help make production of these figures more reproducible, I first import some intermediate data files that I cleaned and manipulated for data vizualization. I also relevel factors here to overide the defalut alphabetical plotting.

``` r
behavior <- read.csv("../data/01a_behavior.csv", header = T)
retention <- read.csv("../data/01a_retention.csv", header = T) 
behaviorsummaryNum <- read.csv("../data/01a_behaviorsummaryNum.csv", header = T)
behaviorsummaryNumAPA2 <- read.csv("../data/01a_behaviorsummaryNumAPA2.csv", header = T)
scaledaveragedata <- read.csv("../data/01a_scaledaveragedata.csv", header = T, row.names = 1)
columnannotations <- read.csv("../data/01a_columnannotations.csv", header = T, row.names = 1)
scoresdf <- read.csv("../data/01a_scoresdf.csv", header = T)
rotationdf <- read.csv("../data/01a_rotationdf.csv", header = T, row.names = 1)
behaviormatrix <- read.csv("../data/01a_behaviormatrix.csv", header = T, row.names = 1)

#set factor levels
behavior$APA <- factor(behavior$APA, levels = c("control", "consistent", "conflict"))

behaviorsummaryNum$APA <- factor(behaviorsummaryNum$APA, levels = c("control", "consistent", "conflict"))
behaviorsummaryNumAPA2$APA2 <- factor(behaviorsummaryNumAPA2$APA2, levels = c("yoked-consistent", "yoked-conflict" ,"consistent", "conflict"))
scoresdf$APA <- factor(scoresdf$APA, levels = c("control", "consistent", "conflict"))
```

Figure 1B: Standard vizualization of mean avoidance beavior
-----------------------------------------------------------

First, I visualze the group mean and standard error for the time it takes before an individual mouse enters the spatial region marked "schock zone" or equivilent (Fig. 1A).

``` r
# plotting mean and se for time to total number of entrances
numentrance1 <- ggplot(behaviorsummaryNum, aes(x=, TrainSessionComboNum, y=m, color=APA)) + 
    geom_errorbar(aes(ymin=m-se, ymax=m+se, color=APA), width=.1) +
    geom_point(size = 2) +
   geom_line() +
    scale_y_continuous(name="Number of Entrances") +
    scale_x_continuous(name = NULL, 
                       breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                       labels=c("1" = "Hab.", "2" = "T1", "3" = "T2", 
                                "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                                "7" = "T5/C2", "8" = "T6/C3", "9"= "Reten.")) +
  theme_cowplot(font_size = 14, line_size = 0.5) +
  background_grid(major = "none", minor = "none") +
  scale_color_manual(values = colorvalAPA) + 
  theme(legend.position=c(0.8, 0.8))  + 
  theme(legend.title=element_blank())
numentrance1
```

![](../figures/01_behavior/numentrance-1.png)

``` r
pdf(file="../figures/01_behavior/numentrance1.pdf", width=6, height=3)
plot(numentrance1)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
numentrance2 <- ggplot(behaviorsummaryNumAPA2, aes(x=, TrainSessionComboNum, y=m, color=APA2)) + 
    geom_errorbar(aes(ymin=m-se, ymax=m+se, color=APA2), width=.1) +
    geom_point(size = 2) +
   geom_line() +
    scale_y_continuous(name="Number of Entrances") +
    scale_x_continuous(name = NULL, 
                       breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                       labels=c("1" = "Hab.", "2" = "T1", "3" = "T2", 
                                "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                                "7" = "T5/C2", "8" = "T6/C3", "9"= "Reten.")) +
  theme_cowplot(font_size = 14, line_size = 0.5) +
  background_grid(major = "none", minor = "none") +
  scale_color_manual(values = colorvalAPA2) + 
  theme(legend.position=c(0.8, 0.8))  + 
  theme(legend.title=element_blank())
numentrance2
```

![](../figures/01_behavior/numentrance-2.png)

``` r
pdf(file="../figures/01_behavior/numentrance2.pdf", width=6, height=3)
plot(numentrance2)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

### Supplemnentary time series plot

Here I visualze the individual data points for each annimal then use a linar model to fit a 95% confidence internal for the time it take a mouse to enter the "schock zone" or equivilent .

``` r
# plotting all data points and linear model smoothing for number of entrances
numentrance3 <- onebehavior(data=behavior, 
                            xcol="TrainSessionComboNum", ycol="NumEntrances",
                  yaxislabel="Number of Entrances",
                  colorcode="APA")

numentrance2
```

![](../figures/01_behavior/numentrance3-1.png)

``` r
pdf(file="../figures/01_behavior/numentrance3.pdf", width=6, height=3)
plot(numentrance3)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

### Supplmentary heatmap

The next image shows how all the behaviors measured change over time. Here, the data are normalized to a z-score with more positive values shown the viridis color scheme with yellow being positive and purple being negative. Each row contains value for each behavioral measurement. Each column is the average value for a group of animals as specific by APA group and training session.

This particular graph clusters the data according to a kmean of 3 rather than by heriarchical clustering. I use the R package superheat.

``` r
scaledaveragedatatranposed <- t(scaledaveragedata)

superheat(scaledaveragedata,
          # change the size of the labels
          left.label.size = 0.3, 
          bottom.label.size = 0.4,
          bottom.label.text.angle = 90, 
          # cluster rows and add dendrogram
          pretty.order.cols = TRUE,
          n.clusters.rows = 3,
          left.label = 'variable',
          heat.lim = c(-1.5, 1.5), 
          extreme.values.na = FALSE)
```

![](../figures/01_behavior/superheat-1.png)

Figure 1C: Hierarchical clusering of time series behavioral data
----------------------------------------------------------------

Here I use heirarhical cluster to identify patterns in the behavioral data. On the y axis see three distinct clusters of behaviors that are 1) higher in trained animals, 2) higher in yoked animals, and 3) measures of speed (Fig. 1C). The row dendrogram helps visuzlise the interaction between treatment group and (grey: control, pink: consistent, red: conflict) and training session (white: habitutation, shades of brown: training sessino 1-3, dark down: retest, shades of purple: training session 4-6 or conflict session 1-3, dark purple: retention).

``` r
## make annotation df and ann_colors for pheatmap
df <- columnannotations
ann_colors = session_colors 

# set color breaks
paletteLength <- 30
myBreaks <- c(seq(min(scaledaveragedata), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(scaledaveragedata)/paletteLength, max(scaledaveragedata), length.out=floor(paletteLength/2)))

## pheatmap for markdown
pheatmap(scaledaveragedata, show_colnames=F, show_rownames = T,
         annotation_col=df, annotation_colors = ann_colors,
         treeheight_row = 0, treeheight_col = 50,
         fontsize = 8, 
         #width=4.5, height=3,
         border_color = "grey60" ,
         color = viridis(30),
         cellwidth = 10, 
         clustering_method="average",
         breaks=myBreaks,
         clustering_distance_cols="correlation" 
         )
```

![](../figures/01_behavior/pheatmap-1.png)

``` r
# pheatmapfor adobe
pheatmap(scaledaveragedata, show_colnames=F, show_rownames = T,
         annotation_col=df, annotation_colors = ann_colors,
         treeheight_row = 0, treeheight_col = 50,
         fontsize = 8, 
         border_color = "grey60" ,
         color = viridis(30),
         cellwidth = 6, 
         clustering_method="average",
         breaks=myBreaks,
         clustering_distance_cols="correlation",
         filename = "../figures/01_behavior/pheatmap.pdf"
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
percent <- 100*pca$sdev^2/sum(pca$sdev^2)
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
pca12 <- makepcaplotwithpercent(data=scoresdf,xcol="PC1",ycol="PC2",colorcode="APA", newxlab = "PC1 (35.7%)", newylab = "PC2 (9.7%)")
pca12
```

![](../figures/01_behavior/PC12-1.png)

``` r
pdf(file="../figures/01_behavior/pca12.pdf", width=3, height=3)
plot(pca12)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
## statistics
aov1 <- aov(PC1 ~ APA, data=scoresdf)
summary(aov1) # p = 2.53e-14
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## APA          2   3214  1606.8   101.3 2.53e-14 ***
    ## Residuals   31    492    15.9                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TukeyHSD(aov1, which = "APA") # p<< 0.001 for both control comparisions
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ APA, data = scoresdf)
    ## 
    ## $APA
    ##                            diff        lwr        upr     p adj
    ## consistent-control  -19.5553836 -23.757552 -15.353215 0.0000000
    ## conflict-control    -19.3444279 -23.384730 -15.304125 0.0000000
    ## conflict-consistent   0.2109558  -4.551502   4.973413 0.9934703

``` r
fviz12 <- fviz_pca_var(res.pca, select.var = list(contrib = 2), axes = c(1, 2))
fviz12
```

![](../figures/01_behavior/PC12-2.png)

``` r
pdf(file="../figures/01_behavior/fviz12.pdf", width=3, height=3)
plot(fviz12)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
pca36 <-makepcaplotwithpercent(data=scoresdf,xcol="PC3",ycol="PC6",colorcode="APA", newxlab = "PC3 (8.4%)", newylab = "PC6 (3.6%)")
pca36
```

![](../figures/01_behavior/PC36-1.png)

``` r
pdf(file="../figures/01_behavior/pca36.pdf", width=3, height=3)
plot(pca36)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
fviz36 <- fviz_pca_var(res.pca, select.var = list(contrib = 2), axes = c(3, 6))
fviz36
```

![](../figures/01_behavior/PC36-2.png)

``` r
pdf(file="../figures/01_behavior/fviz36.pdf", width=3, height=3)
plot(fviz36)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
aov3 <- aov(PC3 ~ APA, data=scoresdf)
summary(aov3) # p = 0.0633
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## APA          2   95.0   47.50    3.02 0.0633 .
    ## Residuals   31  487.5   15.73                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TukeyHSD(aov3, which = "APA") # p = 0.0557503 for conflict-consistent 
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC3 ~ APA, data = scoresdf)
    ## 
    ## $APA
    ##                          diff         lwr      upr     p adj
    ## consistent-control  -1.812712 -5.99726370 2.371839 0.5417745
    ## conflict-control     2.833634 -1.18973029 6.856997 0.2090687
    ## conflict-consistent  4.646346 -0.09614557 9.388837 0.0557503

``` r
aov6 <- aov(PC6 ~ APA, data=scoresdf)
summary(aov6) # p = 0.0226
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## APA          2   96.1   48.03   4.293 0.0226 *
    ## Residuals   31  346.9   11.19                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TukeyHSD(aov6, which = "APA") # p = 0.0175577 for conflict-consistent 
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC6 ~ APA, data = scoresdf)
    ## 
    ## $APA
    ##                          diff        lwr      upr     p adj
    ## consistent-control  -2.152715 -5.6824153 1.376985 0.3043094
    ## conflict-control     2.583185 -0.8105525 5.976922 0.1633660
    ## conflict-consistent  4.735900  0.7355730 8.736227 0.0175577

PC2 and 9 also difference by p &lt; 0.01.

``` r
aov2 <- aov(PC2 ~ APA, data=scoresdf)
summary(aov2) # p = 0.0906
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## APA          2  106.1   53.05   2.597 0.0906 .
    ## Residuals   31  633.2   20.42                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TukeyHSD(aov2, which = "APA") # p = 0.0852578 for conflict-consistent 
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ APA, data = scoresdf)
    ## 
    ## $APA
    ##                          diff        lwr       upr     p adj
    ## consistent-control  -3.434152 -8.2030827  1.334778 0.1955003
    ## conflict-control     1.418463 -3.1667701  6.003695 0.7290903
    ## conflict-consistent  4.852615 -0.5521726 10.257403 0.0852578

``` r
aov9 <- aov(PC9 ~ APA, data=scoresdf)
summary(aov9) # p =  0.0584
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## APA          2  62.55   31.28   3.117 0.0584 .
    ## Residuals   31 311.07   10.04                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TukeyHSD(aov9, which = "APA") # p = 0.0722837 for conflict-consistent 
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC9 ~ APA, data = scoresdf)
    ## 
    ## $APA
    ##                           diff        lwr      upr     p adj
    ## consistent-control  -0.7833752 -4.1260684 2.559318 0.8334263
    ## conflict-control     2.7406860 -0.4732478 5.954620 0.1065599
    ## conflict-consistent  3.5240612 -0.2643244 7.312447 0.0722837

Here are some stats modeling combinatorail PCs. I'm really not sure if this makes sense. I should probabaly model some behavior time interaction....

``` r
lm1 <- lm(PC1~APA, data=scoresdf)
summary(lm1)
```

    ## 
    ## Call:
    ## lm(formula = PC1 ~ APA, data = scoresdf)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -11.5071  -2.5761   0.0012   2.0463   7.0257 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     9.7219     0.9658   10.07 2.74e-11 ***
    ## APAconsistent -19.5554     1.7074  -11.45 1.14e-12 ***
    ## APAconflict   -19.3444     1.6416  -11.78 5.54e-13 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 3.982 on 31 degrees of freedom
    ## Multiple R-squared:  0.8673, Adjusted R-squared:  0.8588 
    ## F-statistic: 101.3 on 2 and 31 DF,  p-value: 2.531e-14

``` r
lm16 <- lm(PC1+PC6~APA, data=scoresdf)
summary(lm16)
```

    ## 
    ## Call:
    ## lm(formula = PC1 + PC6 ~ APA, data = scoresdf)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -12.4439  -3.4814  -0.5956   3.6301  13.1110 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      9.545      1.341   7.119 5.33e-08 ***
    ## APAconsistent  -21.708      2.370  -9.160 2.50e-10 ***
    ## APAconflict    -16.761      2.279  -7.356 2.79e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 5.528 on 31 degrees of freedom
    ## Multiple R-squared:  0.7717, Adjusted R-squared:  0.7569 
    ## F-statistic: 52.38 on 2 and 31 DF,  p-value: 1.143e-10

``` r
lm136 <- lm(PC1+PC3+PC6~APA, data=scoresdf)
summary(lm136)
```

    ## 
    ## Call:
    ## lm(formula = PC1 + PC3 + PC6 ~ APA, data = scoresdf)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -10.755  -4.705  -1.010   3.253  17.750 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      9.221      1.659   5.559 4.32e-06 ***
    ## APAconsistent  -23.521      2.932  -8.021 4.67e-09 ***
    ## APAconflict    -13.928      2.819  -4.940 2.55e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 6.839 on 31 degrees of freedom
    ## Multiple R-squared:  0.6935, Adjusted R-squared:  0.6737 
    ## F-statistic: 35.07 on 2 and 31 DF,  p-value: 1.096e-08
