``` r
library(ggplot2) ## for awesome plots!
```

    ## Warning: package 'ggplot2' was built under R version 3.3.2

``` r
library(cowplot) ## for some easy to use themes
library(dplyr) ## for filtering and selecting rows
library(factoextra)
```

    ## Warning: package 'factoextra' was built under R version 3.3.2

``` r
## load functions 
source("functions_behavior.R")
source("figureoptions.R")

## set output file for figures 
knitr::opts_chunk$set(fig.path = '../figures/01_behavior/')
```

``` r
behavior <- read.csv("../data/01a_behavior.csv", header = T)
retention <- read.csv("../data/01a_retention.csv", header = T) 
behaviorsummaryTime <- read.csv("../data/01a_behaviorsummaryTime.csv", header = T)
behaviorsummaryNum <- read.csv("../data/01a_behaviorsummaryNum.csv", header = T)
scaledaveragedata <- read.csv("../data/01a_scaledaveragedata.csv", header = T, row.names = 1)
columnannotations <- read.csv("../data/01a_columnannotations.csv", header = T, row.names = 1)
scoresdf <- read.csv("../data/01a_scoresdf.csv", header = T)
rotationdf <- read.csv("../data/01a_rotationdf.csv", header = T, row.names = 1)
behaviormatrix <- read.csv("../data/01a_behaviormatrix.csv", header = T, row.names = 1)
```

``` r
behavior$APA <- factor(behavior$APA, levels = c("control", "consistent", "conflict"))
behaviorsummaryTime$APA <- factor(behaviorsummaryTime$APA, levels = c("control", "consistent", "conflict"))
behaviorsummaryNum$APA <- factor(behaviorsummaryNum$APA, levels = c("control", "consistent", "conflict"))
scoresdf$APA <- factor(scoresdf$APA, levels = c("control", "consistent", "conflict"))
```

``` r
A <- ggplot(behaviorsummaryTime, aes(x=APA, y=m, color=APA)) + 
    geom_errorbar(aes(ymin=m-se, ymax=m+se, color=APA), width=.1) +
    geom_point(size = 2) +
    scale_y_continuous(name="Time to First Entrance (s)") +
    scale_x_discrete(name=NULL) +
  theme_cowplot(font_size = 14, line_size = 1) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(values = colorvalAPA) +
theme(legend.justification=c(1,0), legend.position=c(1,0), legend.title = element_blank()) + 
theme(axis.ticks = element_blank(), axis.text.x = element_blank())

B <- ggplot(behaviorsummaryNum, aes(x=APA, y=m, color=APA)) + 
    geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.1) +
    geom_point(size = 2) + 
    scale_y_continuous(name="Number of Entrances") +
    scale_x_discrete(name=NULL) +
  theme_cowplot(font_size = 14, line_size = 1) +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(values = colorvalAPA) +
theme(legend.justification=c(1,1), legend.position=c(1,1), legend.title = element_blank()) + theme(axis.ticks = element_blank(), axis.text.x = element_blank())


plot_grid(A,B, nrow=1, labels=c("A", "B"), rel_widths = c(1, 1))
```

![](../figures/01_behavior/avoidancebehavior-1.png)

``` r
meansem <- plot_grid(A,B, nrow=1, labels=c("A", "B"), rel_widths = c(1, 1))


pdf(file="../figures/01_behavior/avoidancebehavior-1.pdf", width=6, height=3)
plot(meansem)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
A <- myboxplotlegendtop(data = behavior,xcol = "TrainSessionCombo", 
                ycol = "Time1stEntr", colorcode = "APA", session ="Retention",
                yaxislabel="\n Time to 1st Entrance (s)")

B <- myboxplotlegendbottom(data = behavior, xcol = "TrainSessionCombo", 
                ycol = "NumEntrances", colorcode = "APA", 
                session ="Retention",
                yaxislabel="\n Number of Entrances") 

plot_grid(A,B, nrow=1,  labels=c("A", "B"))
```

![](../figures/01_behavior/avoidancebehavior-2.png)

``` r
boxplot <- plot_grid(A,B, nrow=1, labels=c("A", "B"))

pdf(file="../figures/01_behavior/avoidancebehavior-2.pdf", width=8.5, height=3)
plot(boxplot)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
onebehavior(data=behavior, xcol="TrainSessionComboNum", ycol="pTimeOPP",
                  yaxislabel=" Proportion of time spent\n opposite the shock zone",
                  colorcode="APA")
```

![](../figures/01_behavior/avoidancebehavior-3.png)

``` r
timeopp <- onebehavior(data=behavior, xcol="TrainSessionComboNum", ycol="pTimeOPP",
                  yaxislabel=" Proportion of time spent\n opposite the shock zone",
                  colorcode="APA")

pdf(file="../figures/01_behavior/avoidancebehavior-3.pdf", width=6, height=4)
plot(timeopp)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
## hab to retest
habtoretest <- behavior %>% 
  filter(TrainSessionCombo %in% c("Hab", "Retest", "T1", "T3", "T2")) %>% droplevels()
onebehaviorhabtoretest(data=habtoretest, xcol="TrainSessionComboNum", ycol="pTimeOPP",
                  yaxislabel=" Proportion of time spent\n opposite the shock zone",
                  colorcode="APA")
```

    ## `geom_smooth()` using method = 'loess'

![](../figures/01_behavior/avoidancebehavior-4.png)

``` r
## c4 to rentetion
c4toretention <- behavior %>% 
  filter(TrainSessionCombo %in% c("Retention", "T4_C1", "T5_C2", "T6_C3")) %>% droplevels()
onebehaviorc4toRentention(data=c4toretention, xcol="TrainSessionComboNum", ycol="pTimeOPP",
                  yaxislabel=" Proportion of time spent\n opposite the shock zone",
                  colorcode="APA")
```

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric =
    ## parametric, : pseudoinverse used at 5.985

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric =
    ## parametric, : neighborhood radius 2.015

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric =
    ## parametric, : reciprocal condition number 4.896e-17

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric =
    ## parametric, : There are other near singularities as well. 4.0602

    ## Warning in predLoess(object$y, object$x, newx = if
    ## (is.null(newdata)) object$x else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : pseudoinverse used
    ## at 5.985

    ## Warning in predLoess(object$y, object$x, newx = if
    ## (is.null(newdata)) object$x else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : neighborhood radius
    ## 2.015

    ## Warning in predLoess(object$y, object$x, newx = if
    ## (is.null(newdata)) object$x else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : reciprocal
    ## condition number 4.896e-17

    ## Warning in predLoess(object$y, object$x, newx = if
    ## (is.null(newdata)) object$x else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : There are other
    ## near singularities as well. 4.0602

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric =
    ## parametric, : pseudoinverse used at 5.985

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric =
    ## parametric, : neighborhood radius 2.015

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric =
    ## parametric, : reciprocal condition number 4.5798e-17

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric =
    ## parametric, : There are other near singularities as well. 4.0602

    ## Warning in predLoess(object$y, object$x, newx = if
    ## (is.null(newdata)) object$x else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : pseudoinverse used
    ## at 5.985

    ## Warning in predLoess(object$y, object$x, newx = if
    ## (is.null(newdata)) object$x else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : neighborhood radius
    ## 2.015

    ## Warning in predLoess(object$y, object$x, newx = if
    ## (is.null(newdata)) object$x else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : reciprocal
    ## condition number 4.5798e-17

    ## Warning in predLoess(object$y, object$x, newx = if
    ## (is.null(newdata)) object$x else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : There are other
    ## near singularities as well. 4.0602

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric =
    ## parametric, : pseudoinverse used at 5.985

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric =
    ## parametric, : neighborhood radius 2.015

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric =
    ## parametric, : reciprocal condition number 4.896e-17

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric =
    ## parametric, : There are other near singularities as well. 4.0602

    ## Warning in predLoess(object$y, object$x, newx = if
    ## (is.null(newdata)) object$x else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : pseudoinverse used
    ## at 5.985

    ## Warning in predLoess(object$y, object$x, newx = if
    ## (is.null(newdata)) object$x else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : neighborhood radius
    ## 2.015

    ## Warning in predLoess(object$y, object$x, newx = if
    ## (is.null(newdata)) object$x else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : reciprocal
    ## condition number 4.896e-17

    ## Warning in predLoess(object$y, object$x, newx = if
    ## (is.null(newdata)) object$x else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : There are other
    ## near singularities as well. 4.0602

![](../figures/01_behavior/avoidancebehavior-5.png)

``` r
levels(behavior$TrainSessionCombo)
```

    ## [1] "Hab"       "Retention" "Retest"    "T1"        "T2"        "T3"       
    ## [7] "T4_C1"     "T5_C2"     "T6_C3"

### Heatmap

The next image shows how all the behaviors measured change over time. Here, the data are normalized to a z-score with more positive values shown in red and negative values show in blue. Each row contains value for each behavioral measurement. Each column is the average value for a group of animals as specific by APA group (purple, orange, brown) and training session (from white to black according to increasing time spend in the active place avoidance group).

``` r
## see the makesessionheatmap documentataion for data tidying and plot specifications
#makesessionheatmap(behavior)
```

### Principle component analysis (PCA)

Given the correlational structure of the data, I next reduced the dimentionality with a PCA anlaysis. You can see that PC1 speparates trained and untraned animals (D,E) but neither PC2 (D) nor PC3 (E) separate same and conflict aniamls. Elipses show 95% confidence interval.

``` r
makepcaplot(data=scoresdf,xcol="PC1",ycol="PC2",colorcode="APA")
```

![](../figures/01_behavior/PCA-1.png)

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
makepcaplot(data=scoresdf,xcol="PC1",ycol="PC3",colorcode="APA")
```

![](../figures/01_behavior/PCA-2.png)

``` r
makepcaplot(data=scoresdf,xcol="PC1",ycol="PC4",colorcode="APA")
```

![](../figures/01_behavior/PCA-3.png)

``` r
makepcaplot(data=scoresdf,xcol="PC1",ycol="PC5",colorcode="APA")
```

![](../figures/01_behavior/PCA-4.png)

``` r
makepcaplot(data=scoresdf,xcol="PC1",ycol="PC6",colorcode="APA")
```

![](../figures/01_behavior/PCA-5.png)

``` r
## statistics
aov6 <- aov(PC6 ~ APA, data=scoresdf)
summary(aov6) # p = 0.0584 
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## APA          2   96.1   48.03   4.293 0.0226 *
    ## Residuals   31  346.9   11.19                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TukeyHSD(aov6, which = "APA") # consistent-conflict p= 0.0175577
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

``` r
makepcaplot(data=scoresdf,xcol="PC1",ycol="PC7",colorcode="APA")
```

![](../figures/01_behavior/PCA-6.png)

``` r
makepcaplot(data=scoresdf,xcol="PC1",ycol="PC8",colorcode="APA")
```

![](../figures/01_behavior/PCA-7.png)

``` r
makepcaplot(data=scoresdf,xcol="PC1",ycol="PC9",colorcode="APA")
```

![](../figures/01_behavior/PCA-8.png)

``` r
## statistics
aov9 <- aov(PC9 ~ APA, data=scoresdf)
summary(aov9) # p = 0.0584
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## APA          2  62.55   31.28   3.117 0.0584 .
    ## Residuals   31 311.07   10.04                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
TukeyHSD(aov9, which = "APA") # consistent-conflict p= 0.0722
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

``` r
makepcaplot(data=scoresdf,xcol="PC1",ycol="PC6",colorcode="APA")
```

![](../figures/01_behavior/PCA-9.png)

``` r
makepcaplot(data=scoresdf,xcol="PC1",ycol="PC9",colorcode="APA")
```

![](../figures/01_behavior/PCA-10.png)

``` r
res.pca <- prcomp(behaviormatrix,  scale = TRUE)
fviz_pca_var(res.pca, select.var = list(contrib = 8), axes = c(1, 2))
```

![](../figures/01_behavior/fviz_pca-1.png)

``` r
fviz_pca_var(res.pca, select.var = list(contrib = 8), axes = c(1, 6))
```

![](../figures/01_behavior/fviz_pca-2.png)

``` r
fviz_pca_var(res.pca, select.var = list(contrib = 8), axes = c(1, 9))
```

![](../figures/01_behavior/fviz_pca-3.png)

``` r
library(ggfortify)
```

    ## Warning: package 'ggfortify' was built under R version 3.3.2

``` r
df <- behaviormatrix
autoplot(prcomp(df))
```

![](../figures/01_behavior/otherpca-1.png)

``` r
autoplot(prcomp(df), data = behavior, colour = 'APA', loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5, scale = 0)
```

![](../figures/01_behavior/otherpca-2.png)

``` r
library(cluster)
autoplot(clara(behaviormatrix, 3))
```

![](../figures/01_behavior/otherpca-3.png)

``` r
autoplot(fanny(behaviormatrix, 2), frame = TRUE)
```

    ## Warning in fanny(behaviormatrix, 2): the memberships are all very close to
    ## 1/k. Maybe decrease 'memb.exp' ?

![](../figures/01_behavior/otherpca-4.png)

``` r
# capture the rotation matrix in a data frame
rotation_data <- rotationdf
# define a pleasing arrow style
arrow_style <- arrow(length = unit(0.05, "inches"),
                     type = "closed")
# now plot, using geom_segment() for arrows and geom_text for labels
ggplot(rotation_data) + 
  geom_segment(aes(xend=PC1, yend=PC2), x=0, y=0, arrow=arrow_style) + 
  geom_text(aes(x=PC1, y=PC2, label=variable), hjust=0, size=3, color='red')
```

![](../figures/01_behavior/wilkepca-1.png)