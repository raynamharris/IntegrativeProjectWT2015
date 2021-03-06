    library(tidyverse)

    ## ── Attaching packages ──────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    library(corrr)
    library(cowplot)

    ## 
    ## ********************************************************

    ## Note: As of version 1.0.0, cowplot does not change the

    ##   default ggplot2 theme anymore. To recover the previous

    ##   behavior, execute:
    ##   theme_set(theme_cowplot())

    ## ********************************************************

    library(Hmisc) # for correlations with pvalue

    ## Loading required package: lattice

    ## Loading required package: survival

    ## Loading required package: Formula

    ## 
    ## Attaching package: 'Hmisc'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     src, summarize

    ## The following objects are masked from 'package:base':
    ## 
    ##     format.pval, units

    library(ggrepel)


    source("./figureoptions.R")
    source("./functions_RNAseq.R")

Sample information and PC1
--------------------------

    # read the sample data, set levels, join iwth behvior PCA data
    colData <- read.csv("../data/00_colData.csv", row.names = 1, stringsAsFactors = T) 
    pca.Rn <- read_csv("../data/suppltable-3.csv") %>% dplyr::filter(trialNum == 9)

    ## Parsed with column specification:
    ## cols(
    ##   ID = col_character(),
    ##   treatment = col_character(),
    ##   trial = col_character(),
    ##   trialNum = col_double(),
    ##   PC1 = col_double(),
    ##   PC2 = col_double()
    ## )

    colData <- left_join(colData, pca.Rn)

    ## Joining, by = c("ID", "treatment")

    ## Warning: Column `ID` joining factor and character vector, coercing into
    ## character vector

    ## Warning: Column `treatment` joining factor and character vector, coercing into
    ## character vector

    head(colData)

    ##       ID subfield        treatment training     trial trialNum        PC1
    ## 1 15143A      CA3 conflict.trained  trained Retention        9 -0.2275039
    ## 2 15143A       DG conflict.trained  trained Retention        9 -0.2275039
    ## 3 15143B      CA1   conflict.yoked    yoked Retention        9 -3.1436627
    ## 4 15143B       DG   conflict.yoked    yoked Retention        9 -3.1436627
    ## 5 15143C      CA1 standard.trained  trained Retention        9  5.8860390
    ## 6 15143D      CA1   standard.yoked    yoked Retention        9 -2.7532719
    ##           PC2
    ## 1  3.03543738
    ## 2  3.03543738
    ## 3 -0.48834291
    ## 4 -0.48834291
    ## 5  0.33235063
    ## 6 -0.07584809

    # read all count data prep to join with sample colData

    combinePCvsd <- function(filename, whichsubfield, whichgenes){
      
      
      vsd <- read.csv(filename, row.names = 1, check.names = F) 
      vsd$gene <- row.names(vsd)
      vsd$gene <- str_to_title(vsd$gene)
      vsd <- as.data.frame(vsd)
      row.names(vsd) <- vsd$gene
      vsd$gene <- NULL
      vsd <- as.data.frame(t(vsd))
      vsd$sample <- row.names(vsd)
      vsd$mouse <- sapply(strsplit(as.character(vsd$sample),"\\-"), "[", 1)
      vsd$ID <- paste(15, vsd$mouse, sep = "")
      
      vsd <- left_join(colData, vsd) %>%
        select(ID:training, 
               PC1, PC2,
               whichgenes) %>%
        filter(subfield == whichsubfield) 
      print(head(vsd))
      return(vsd)
    }

    candidategenes <- c("Prkcz", "Prkci", "Wwc1", "Grin1", "Gria1",
             "Pick1", "Nsf", "Fmr1", "Camk2a",
             "Fos", "Fosl2", "Npas4", "Arc")

    vsdDG <- combinePCvsd("../data/03_DG_vsdtraining.csv", "DG", candidategenes)

    ## Joining, by = "ID"

    ##       ID subfield        treatment training        PC1         PC2    Prkcz
    ## 1 15143A       DG conflict.trained  trained -0.2275039  3.03543738 9.052238
    ## 2 15143B       DG   conflict.yoked    yoked -3.1436627 -0.48834291 8.922782
    ## 3 15143D       DG   standard.yoked    yoked -2.7532719 -0.07584809 9.060078
    ## 4 15144A       DG conflict.trained  trained  6.7041815 -0.07853719 9.004395
    ## 5 15144C       DG standard.trained  trained  7.0499369 -1.78499206 8.903003
    ## 6 15144D       DG   standard.yoked    yoked -3.3026284  1.17314374 9.094801
    ##      Prkci     Wwc1    Grin1    Gria1    Pick1      Nsf     Fmr1   Camk2a
    ## 1 7.719338 8.301099 10.68066 11.52905 7.530504 11.49573 8.027732 13.21496
    ## 2 7.548660 8.215442 10.80562 11.42916 7.327508 11.42863 8.205092 12.63324
    ## 3 8.006286 8.350489 10.83083 11.63562 7.717724 11.41348 8.200184 12.89590
    ## 4 7.852252 8.451785 10.69996 11.75007 7.474781 11.61321 8.401654 13.45862
    ## 5 7.888498 8.410608 10.73685 11.67325 7.272327 11.62802 8.219559 13.01407
    ## 6 7.736955 8.438542 10.69857 11.75503 7.611018 11.72119 8.012141 13.18540
    ##         Fos     Fosl2     Npas4       Arc
    ## 1 11.273104 10.047064 11.707078 11.248244
    ## 2 10.859684  9.277070  9.895664  9.360658
    ## 3 10.440561  9.137341 10.057463  9.740559
    ## 4 11.515130 11.374861 12.733872 11.465318
    ## 5 11.765583 11.054132 12.631030 12.124175
    ## 6  8.213153  8.489472  8.594899  8.468157

    vsdCA1 <- combinePCvsd("../data/03_CA3_vsdtraining.csv", "CA1", candidategenes)

    ## Joining, by = "ID"

    ##       ID subfield        treatment training       PC1         PC2    Prkcz
    ## 1 15143B      CA1   conflict.yoked    yoked -3.143663 -0.48834291       NA
    ## 2 15143C      CA1 standard.trained  trained  5.886039  0.33235063       NA
    ## 3 15143D      CA1   standard.yoked    yoked -2.753272 -0.07584809       NA
    ## 4 15144A      CA1 conflict.trained  trained  6.704182 -0.07853719 8.480993
    ## 5 15144B      CA1   conflict.yoked    yoked -2.347673 -0.02836475 9.020478
    ## 6 15144C      CA1 standard.trained  trained  7.049937 -1.78499206 9.330134
    ##      Prkci     Wwc1    Grin1    Gria1    Pick1      Nsf     Fmr1   Camk2a
    ## 1       NA       NA       NA       NA       NA       NA       NA       NA
    ## 2       NA       NA       NA       NA       NA       NA       NA       NA
    ## 3       NA       NA       NA       NA       NA       NA       NA       NA
    ## 4 7.257581 7.604937 10.42402 10.87217 7.515542 11.42011 7.561053 12.29037
    ## 5 7.789876 8.568090 10.77667 11.30591 7.449523 11.67896 7.789876 12.94293
    ## 6 7.491501 7.761925 10.40600 11.09883 7.283242 11.49317 7.952061 12.78681
    ##         Fos    Fosl2    Npas4      Arc
    ## 1        NA       NA       NA       NA
    ## 2        NA       NA       NA       NA
    ## 3        NA       NA       NA       NA
    ## 4  9.855929 7.976064 7.515542 9.208575
    ## 5 10.205795 7.525519 9.020478 9.821415
    ## 6 10.384651 7.824398 7.264073 9.413381

    vsdCA3 <- combinePCvsd("../data/03_DG_vsdtraining.csv", "CA3", candidategenes)

    ## Joining, by = "ID"

    ##       ID subfield        treatment training        PC1         PC2    Prkcz
    ## 1 15143A      CA3 conflict.trained  trained -0.2275039  3.03543738 9.052238
    ## 2 15144A      CA3 conflict.trained  trained  6.7041815 -0.07853719 9.004395
    ## 3 15144B      CA3   conflict.yoked    yoked -2.3476730 -0.02836475       NA
    ## 4 15144C      CA3 standard.trained  trained  7.0499369 -1.78499206 8.903003
    ## 5 15144D      CA3   standard.yoked    yoked -3.3026284  1.17314374 9.094801
    ## 6 15145A      CA3 conflict.trained  trained  6.9392690  0.25866181 8.988535
    ##      Prkci     Wwc1    Grin1    Gria1    Pick1      Nsf     Fmr1   Camk2a
    ## 1 7.719338 8.301099 10.68066 11.52905 7.530504 11.49573 8.027732 13.21496
    ## 2 7.852252 8.451785 10.69996 11.75007 7.474781 11.61321 8.401654 13.45862
    ## 3       NA       NA       NA       NA       NA       NA       NA       NA
    ## 4 7.888498 8.410608 10.73685 11.67325 7.272327 11.62802 8.219559 13.01407
    ## 5 7.736955 8.438542 10.69857 11.75503 7.611018 11.72119 8.012141 13.18540
    ## 6 7.965716 8.296659 10.85650 11.69118 7.260001 11.63924 8.056599 12.98003
    ##         Fos     Fosl2     Npas4       Arc
    ## 1 11.273104 10.047064 11.707078 11.248244
    ## 2 11.515130 11.374861 12.733872 11.465318
    ## 3        NA        NA        NA        NA
    ## 4 11.765583 11.054132 12.631030 12.124175
    ## 5  8.213153  8.489472  8.594899  8.468157
    ## 6 10.314514  9.471430 11.342022 10.184857

    # save only genes in all
    DGgenes <- names(vsdDG)
    CA3genes <- names(vsdCA3)
    CA1genes <- names(vsdCA1)
    savecols <- intersect(DGgenes, CA3genes)
    savecols <- intersect(savecols, CA1genes)


    vsdDG <- vsdDG %>% dplyr::select(one_of(savecols)) 
    vsdCA1 <- vsdCA1 %>% dplyr::select(one_of(savecols)) 
    vsdCA3 <- vsdCA3 %>% dplyr::select(one_of(savecols)) 

    allvsd <- rbind(vsdDG, vsdCA1)
    allvsd <- rbind(allvsd, vsdCA3)
    head(allvsd)

    ##       ID subfield        treatment training        PC1         PC2    Prkcz
    ## 1 15143A       DG conflict.trained  trained -0.2275039  3.03543738 9.052238
    ## 2 15143B       DG   conflict.yoked    yoked -3.1436627 -0.48834291 8.922782
    ## 3 15143D       DG   standard.yoked    yoked -2.7532719 -0.07584809 9.060078
    ## 4 15144A       DG conflict.trained  trained  6.7041815 -0.07853719 9.004395
    ## 5 15144C       DG standard.trained  trained  7.0499369 -1.78499206 8.903003
    ## 6 15144D       DG   standard.yoked    yoked -3.3026284  1.17314374 9.094801
    ##      Prkci     Wwc1    Grin1    Gria1    Pick1      Nsf     Fmr1   Camk2a
    ## 1 7.719338 8.301099 10.68066 11.52905 7.530504 11.49573 8.027732 13.21496
    ## 2 7.548660 8.215442 10.80562 11.42916 7.327508 11.42863 8.205092 12.63324
    ## 3 8.006286 8.350489 10.83083 11.63562 7.717724 11.41348 8.200184 12.89590
    ## 4 7.852252 8.451785 10.69996 11.75007 7.474781 11.61321 8.401654 13.45862
    ## 5 7.888498 8.410608 10.73685 11.67325 7.272327 11.62802 8.219559 13.01407
    ## 6 7.736955 8.438542 10.69857 11.75503 7.611018 11.72119 8.012141 13.18540
    ##         Fos     Fosl2     Npas4       Arc
    ## 1 11.273104 10.047064 11.707078 11.248244
    ## 2 10.859684  9.277070  9.895664  9.360658
    ## 3 10.440561  9.137341 10.057463  9.740559
    ## 4 11.515130 11.374861 12.733872 11.465318
    ## 5 11.765583 11.054132 12.631030 12.124175
    ## 6  8.213153  8.489472  8.594899  8.468157

    summary(lm( PC1 ~ Arc, data = vsdDG))

    ## 
    ## Call:
    ## lm(formula = PC1 ~ Arc, data = vsdDG)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -4.7195 -1.8538  0.1447  1.4708  5.6309 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  -29.184      6.810  -4.286 0.001059 ** 
    ## Arc            2.994      0.660   4.536 0.000682 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.748 on 12 degrees of freedom
    ##   (2 observations deleted due to missingness)
    ## Multiple R-squared:  0.6316, Adjusted R-squared:  0.6009 
    ## F-statistic: 20.58 on 1 and 12 DF,  p-value: 0.0006823

    summary(lm( PC2 ~ Arc, data = vsdDG))

    ## 
    ## Call:
    ## lm(formula = PC2 ~ Arc, data = vsdDG)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.9226 -1.0067 -0.5923  0.5110  3.2901 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)   9.8100     3.7627   2.607   0.0229 *
    ## Arc          -0.8948     0.3647  -2.454   0.0304 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.519 on 12 degrees of freedom
    ##   (2 observations deleted due to missingness)
    ## Multiple R-squared:  0.3341, Adjusted R-squared:  0.2786 
    ## F-statistic:  6.02 on 1 and 12 DF,  p-value: 0.03039

    cor.test(vsdDG$PC1, vsdDG$Arc, method = c("pearson"))

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  vsdDG$PC1 and vsdDG$Arc
    ## t = 4.5362, df = 12, p-value = 0.0006823
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.4568067 0.9322321
    ## sample estimates:
    ##       cor 
    ## 0.7947587

    cor.test(vsdDG$PC2, vsdDG$Arc, method = c("pearson"))

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  vsdDG$PC2 and vsdDG$Arc
    ## t = -2.4536, df = 12, p-value = 0.03039
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.84839573 -0.06839018
    ## sample estimates:
    ##        cor 
    ## -0.5779963

    iconDG <- png::readPNG("../figures/00_schematics/DG.png")
    iconDG <-  grid::rasterGrob(iconDG, interpolate=TRUE)

    vsdDG$training <- factor(vsdDG$training, levels = levelstraining)

    a <- ggplot(vsdDG, aes(x = Arc, y = PC1)) +
       geom_point(aes( color = training)) + 
       geom_smooth(method = "lm", color = "grey") +
       scale_color_manual(values = allcolors) +
      theme_ms() +
       theme(legend.position = "bottom",
             axis.title.x = element_text(face = "italic"),
             legend.title = element_blank(), 
             legend.key.size = unit(0.25, "cm"))  +
      labs(subtitle = "r = 0.81, p = 0.0002") +
       annotation_custom(iconDG, ymin = 6, ymax = 11, xmin = 7.5, xmax = 9)
    a

    ## `geom_smooth()` using formula 'y ~ x'

    ## Warning: Removed 2 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 2 rows containing missing values (geom_point).

![](../figures/04_correlations/Arc-1.png)
