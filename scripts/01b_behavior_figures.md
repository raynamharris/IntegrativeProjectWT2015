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
# list of all behaviors measured
names(behavior[c(19:57)])
```

    ##  [1] "TrainSessionComboNum" "SdevSpeedArena"       "Linearity.Arena."    
    ##  [4] "NumEntrances"         "Time1stEntr"          "Path1stEntr"         
    ##  [7] "Speed1stEntr.cm.s."   "Dist1stEntr.m."       "NumShock"            
    ## [10] "MaxTimeAvoid"         "Time2ndEntr"          "Path2ndEntr"         
    ## [13] "Speed2ndEntr"         "TimeTarget"           "pTimeTarget"         
    ## [16] "pTimeCCW"             "pTimeOPP"             "pTimeCW"             
    ## [19] "RayleigLength"        "RayleigAngle"         "PolarAvgVal"         
    ## [22] "PolarSdVal"           "PolarMinVal"          "PolarMinBin"         
    ## [25] "Min50.RngLoBin"       "Min50.RngHiBin"       "PolarMaxVal"         
    ## [28] "PolarMaxBin"          "Max50.RngLoBin"       "Max50.RngHiBin"      
    ## [31] "AnnularMinVal"        "AnnularMinBin"        "AnnularMaxVal"       
    ## [34] "AnnularMaxBin"        "AnnularAvg"           "AnnularSd"           
    ## [37] "AnnularSkewnes"       "AnnularKurtosis"      "Speed1"

``` r
threeplots <- read.csv("../data/01a_threeplots.csv", header = T)
scoresdf <- read.csv("../data/01a_scoresdf.csv", header = T)
rotationdf <- read.csv("../data/01a_rotationdf.csv", header = T, row.names = 1)
behaviormatrix <- read.csv("../data/01a_behaviormatrix.csv", header = T, row.names = 1)

trained <- behavior %>% filter(TrainGroup == "trained") %>% droplevels()

#set factor levels
behavior$APA2 <- factor(behavior$APA2, levels = c("yoked-consistent" ,"consistent", "yoked-conflict", "conflict"))
scoresdf$APA <- NULL
scoresdf$APA2 <- factor(scoresdf$APA2, levels = c("yoked-consistent" ,"consistent", "yoked-conflict", "conflict"))
threeplots$APA2 <- factor(threeplots$APA2, levels = c("yoked-consistent" ,"consistent", "yoked-conflict", "conflict"))
threeplots$measure <- factor(threeplots$measure, levels = c("Number of Entrances" , "Max Avoidance Time", "Speed"))
```

Summary Statistics
==================

``` r
# sample sizes
behavior %>% 
  filter(TrainSession == "Retention") %>%
  select(APA2)  %>%  summary()
```

    ##                APA2  
    ##  yoked-consistent:8  
    ##  consistent      :8  
    ##  yoked-conflict  :9  
    ##  conflict        :9

``` r
names(behavior[c(19:57)])
```

    ##  [1] "TrainSessionComboNum" "SdevSpeedArena"       "Linearity.Arena."    
    ##  [4] "NumEntrances"         "Time1stEntr"          "Path1stEntr"         
    ##  [7] "Speed1stEntr.cm.s."   "Dist1stEntr.m."       "NumShock"            
    ## [10] "MaxTimeAvoid"         "Time2ndEntr"          "Path2ndEntr"         
    ## [13] "Speed2ndEntr"         "TimeTarget"           "pTimeTarget"         
    ## [16] "pTimeCCW"             "pTimeOPP"             "pTimeCW"             
    ## [19] "RayleigLength"        "RayleigAngle"         "PolarAvgVal"         
    ## [22] "PolarSdVal"           "PolarMinVal"          "PolarMinBin"         
    ## [25] "Min50.RngLoBin"       "Min50.RngHiBin"       "PolarMaxVal"         
    ## [28] "PolarMaxBin"          "Max50.RngLoBin"       "Max50.RngHiBin"      
    ## [31] "AnnularMinVal"        "AnnularMinBin"        "AnnularMaxVal"       
    ## [34] "AnnularMaxBin"        "AnnularAvg"           "AnnularSd"           
    ## [37] "AnnularSkewnes"       "AnnularKurtosis"      "Speed1"

``` r
names(behavior)
```

    ##  [1] "ID"                   "Year"                 "Genotype"            
    ##  [4] "TrainProtocol"        "TrainSequence"        "TrainGroup"          
    ##  [7] "Day"                  "TrainSession"         "ShockOnOff"          
    ## [10] "PairedPartner"        "Experimenter"         "Housing"             
    ## [13] "TestLocation"         "APA"                  "APA2"                
    ## [16] "TrainSessionCombo"    "pair1"                "pair2"               
    ## [19] "TrainSessionComboNum" "SdevSpeedArena"       "Linearity.Arena."    
    ## [22] "NumEntrances"         "Time1stEntr"          "Path1stEntr"         
    ## [25] "Speed1stEntr.cm.s."   "Dist1stEntr.m."       "NumShock"            
    ## [28] "MaxTimeAvoid"         "Time2ndEntr"          "Path2ndEntr"         
    ## [31] "Speed2ndEntr"         "TimeTarget"           "pTimeTarget"         
    ## [34] "pTimeCCW"             "pTimeOPP"             "pTimeCW"             
    ## [37] "RayleigLength"        "RayleigAngle"         "PolarAvgVal"         
    ## [40] "PolarSdVal"           "PolarMinVal"          "PolarMinBin"         
    ## [43] "Min50.RngLoBin"       "Min50.RngHiBin"       "PolarMaxVal"         
    ## [46] "PolarMaxBin"          "Max50.RngLoBin"       "Max50.RngHiBin"      
    ## [49] "AnnularMinVal"        "AnnularMinBin"        "AnnularMaxVal"       
    ## [52] "AnnularMaxBin"        "AnnularAvg"           "AnnularSd"           
    ## [55] "AnnularSkewnes"       "AnnularKurtosis"      "Speed1"              
    ## [58] "Speed2"               "Time1stEntrLog"

``` r
slim1 <- behavior[,c(15,16,14,20:59)]
slim2 <- slim1 %>% filter(TrainSessionCombo == "T4_C1", APA != "control") 
slim3 <- as.data.frame(slim2[,1])
slim4 <- slim2[,c(4:43)]

slim3$APA2 <- factor(slim2$APA2, levels = c("consistent", "conflict"))

 
for(y in names(slim4)){
  ymod<- summary(aov(slim4[[y]] ~ slim3$APA2 ))
  cat(paste('\nDependent var:', y, '\n'))
  print(ymod)
}
```

    ## 
    ## Dependent var: SdevSpeedArena 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1 0.0107 0.01071    0.11  0.745
    ## Residuals   15 1.4629 0.09753               
    ## 
    ## Dependent var: Linearity.Arena. 
    ##             Df  Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1 0.00250 0.00250   2.428   0.14
    ## Residuals   15 0.01544 0.00103               
    ## 
    ## Dependent var: NumEntrances 
    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## slim3$APA2   1  853.3   853.3   17.49 0.000801 ***
    ## Residuals   15  731.7    48.8                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Time1stEntr 
    ##             Df Sum Sq Mean Sq F value  Pr(>F)   
    ## slim3$APA2   1 310090  310090   13.71 0.00213 **
    ## Residuals   15 339379   22625                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Path1stEntr 
    ##             Df Sum Sq Mean Sq F value  Pr(>F)   
    ## slim3$APA2   1  212.8  212.81   14.16 0.00188 **
    ## Residuals   15  225.4   15.03                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Speed1stEntr.cm.s. 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1    0.5    0.54    0.01  0.922
    ## Residuals   15  807.3   53.82               
    ## 
    ## Dependent var: Dist1stEntr.m. 
    ##             Df Sum Sq Mean Sq F value  Pr(>F)   
    ## slim3$APA2   1  1.897  1.8965   16.04 0.00115 **
    ## Residuals   15  1.774  0.1183                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: NumShock 
    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## slim3$APA2   1 1122.2    1122   21.99 0.000291 ***
    ## Residuals   15  765.6      51                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: MaxTimeAvoid 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1  56936   56936   2.224  0.157
    ## Residuals   15 384022   25601               
    ## 
    ## Dependent var: Time2ndEntr 
    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## slim3$APA2   1 472992  472992   20.49 0.000402 ***
    ## Residuals   15 346269   23085                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Path2ndEntr 
    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## slim3$APA2   1  324.9   324.9   19.26 0.000529 ***
    ## Residuals   15  253.1    16.9                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Speed2ndEntr 
    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## slim3$APA2   1  70.98   70.98     4.8 0.0447 *
    ## Residuals   15 221.80   14.79                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: TimeTarget 
    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## slim3$APA2   1 1283.1  1283.1   34.88 2.88e-05 ***
    ## Residuals   15  551.8    36.8                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: pTimeTarget 
    ##             Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## slim3$APA2   1 0.005264 0.005264   22.61 0.000256 ***
    ## Residuals   15 0.003493 0.000233                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: pTimeCCW 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1  0.009 0.00901   0.113  0.741
    ## Residuals   15  1.194 0.07960               
    ## 
    ## Dependent var: pTimeOPP 
    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## slim3$APA2   1 0.3032 0.30322   7.371  0.016 *
    ## Residuals   15 0.6171 0.04114                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: pTimeCW 
    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## slim3$APA2   1 0.3284  0.3284   4.288 0.0561 .
    ## Residuals   15 1.1489  0.0766                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: RayleigLength 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1 0.0235 0.02347   0.855   0.37
    ## Residuals   15 0.4116 0.02744               
    ## 
    ## Dependent var: RayleigAngle 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1    506     506   0.108  0.747
    ## Residuals   15  70164    4678               
    ## 
    ## Dependent var: PolarAvgVal 
    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## slim3$APA2   1   5637    5637   6.647  0.021 *
    ## Residuals   15  12721     848                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: PolarSdVal 
    ##             Df Sum Sq Mean Sq F value  Pr(>F)   
    ## slim3$APA2   1   7426    7426   11.46 0.00407 **
    ## Residuals   15   9717     648                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: PolarMinVal 
    ##             Df    Sum Sq   Mean Sq F value Pr(>F)  
    ## slim3$APA2   1 1.930e-06 1.930e-06   6.062 0.0264 *
    ## Residuals   15 4.775e-06 3.183e-07                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: PolarMinBin 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1  11550   11550    0.71  0.413
    ## Residuals   15 244156   16277               
    ## 
    ## Dependent var: Min50.RngLoBin 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1    356     356    0.07  0.795
    ## Residuals   15  76550    5103               
    ## 
    ## Dependent var: Min50.RngHiBin 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1     21      21   0.004  0.953
    ## Residuals   15  89156    5944               
    ## 
    ## Dependent var: PolarMaxVal 
    ##             Df   Sum Sq   Mean Sq F value Pr(>F)
    ## slim3$APA2   1 0.000006 0.0000055   0.008   0.93
    ## Residuals   15 0.010249 0.0006833               
    ## 
    ## Dependent var: PolarMaxBin 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1   1765    1765   0.268  0.612
    ## Residuals   15  98788    6586               
    ## 
    ## Dependent var: Max50.RngLoBin 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1    284     284   0.051  0.825
    ## Residuals   15  84410    5627               
    ## 
    ## Dependent var: Max50.RngHiBin 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1    971     971   0.188   0.67
    ## Residuals   15  77276    5152               
    ## 
    ## Dependent var: AnnularMinVal 
    ##             Df    Sum Sq   Mean Sq F value Pr(>F)
    ## slim3$APA2   1 1.820e-06 1.820e-06   0.177   0.68
    ## Residuals   15 1.546e-04 1.031e-05               
    ## 
    ## Dependent var: AnnularMinBin 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1    1.4    1.35   0.026  0.874
    ## Residuals   15  783.1   52.21               
    ## 
    ## Dependent var: AnnularMaxVal 
    ##             Df  Sum Sq  Mean Sq F value Pr(>F)
    ## slim3$APA2   1 0.00563 0.005628   0.905  0.357
    ## Residuals   15 0.09329 0.006219               
    ## 
    ## Dependent var: AnnularMaxBin 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1  0.156  0.1556   0.168  0.688
    ## Residuals   15 13.915  0.9277               
    ## 
    ## Dependent var: AnnularAvg 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1  0.007  0.0068   0.015  0.903
    ## Residuals   15  6.646  0.4431               
    ## 
    ## Dependent var: AnnularSd 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1   0.64   0.641   0.059  0.811
    ## Residuals   15 161.71  10.781               
    ## 
    ## Dependent var: AnnularSkewnes 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1   0.00   0.000       0  0.998
    ## Residuals   15  20.59   1.373               
    ## 
    ## Dependent var: AnnularKurtosis 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2   1   13.4   13.38   0.275  0.607
    ## Residuals   15  728.7   48.58               
    ## 
    ## Dependent var: Speed1 
    ##             Df    Sum Sq   Mean Sq F value Pr(>F)  
    ## slim3$APA2   1 0.0005401 0.0005401   4.525 0.0504 .
    ## Residuals   15 0.0017905 0.0001194                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Speed2 
    ##             Df    Sum Sq   Mean Sq F value  Pr(>F)   
    ## slim3$APA2   1 0.0004692 0.0004692    14.4 0.00176 **
    ## Residuals   15 0.0004886 0.0000326                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Time1stEntrLog 
    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## slim3$APA2   1  74.93   74.93   39.13 1.54e-05 ***
    ## Residuals   15  28.72    1.91                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# *** Time1stEntrLog, pTimeTarget, TimeTarget, Path2ndEntr, Time2ndEntr, NumEntrances
# **  Speed2, PolarSdVal, Dist1stEntr.m, Time1stEntr
# *   PolarMinVal, PolarAvgVal, pTimeOPP, Speed2ndEntr
# .   Speed1, pTimeCW
# not AnnularKurtosis, AnnularSkewnes, AnnularSd, AnnularAvg, AnnularMinBin, AnnularMaxBin, Max50.RngHiBin, Max50.RngLoBin, PolarMaxBin, PolarMaxVal, Min50.RngHiBin, Min50.RngLoBin, PolarMinBin, RayleigAngle, pTimeCCW, MaxTimeAvoid, Linearity.Arena, SdevSpeedArena

for(y in names(slim4)){
  ymod <- t.test(slim4[[y]] ~ slim3$APA2 )
  cat(paste('\nDependent var:', y, '\n'))
  print(ymod)
}
```

    ## 
    ## Dependent var: SdevSpeedArena 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -0.32727, df = 13.612, p-value = 0.7484
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.3806588  0.2801033
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 2.697500                 2.747778 
    ## 
    ## 
    ## Dependent var: Linearity.Arena. 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 1.5173, df = 11.627, p-value = 0.1559
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.01071673  0.05930840
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                0.3237625                0.2994667 
    ## 
    ## 
    ## Dependent var: NumEntrances 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -4.3217, df = 12.816, p-value = 0.0008563
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -21.300512  -7.088377
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                  4.25000                 18.44444 
    ## 
    ## 
    ## Dependent var: Time1stEntr 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 3.4769, df = 7.0325, p-value = 0.01023
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##   86.7354 454.4327
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##               278.576250                 7.992222 
    ## 
    ## 
    ## Dependent var: Path1stEntr 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 3.5347, df = 7.0353, p-value = 0.00946
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##   2.351236 11.825709
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                7.2962500                0.2077778 
    ## 
    ## 
    ## Dependent var: Speed1stEntr.cm.s. 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 0.10298, df = 13.056, p-value = 0.9195
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -7.105767  7.817433
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 3.622500                 3.266667 
    ## 
    ## 
    ## Dependent var: Dist1stEntr.m. 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -4.1096, df = 13.783, p-value = 0.001095
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -1.0189237 -0.3194096
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                0.2375000                0.9066667 
    ## 
    ## 
    ## Dependent var: NumShock 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -4.8544, df = 12.519, p-value = 0.0003503
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -23.550278  -9.005278
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                  4.50000                 20.77778 
    ## 
    ## 
    ## Dependent var: MaxTimeAvoid 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 1.4636, df = 12.762, p-value = 0.1675
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -55.52603 287.41492
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 316.5000                 200.5556 
    ## 
    ## 
    ## Dependent var: Time2ndEntr 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 4.2515, df = 7.0392, p-value = 0.003738
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  148.5268 519.8407
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 352.6838                  18.5000 
    ## 
    ## 
    ## Dependent var: Path2ndEntr 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 4.1214, df = 7.0325, p-value = 0.004408
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##   3.738084 13.779138
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                9.3675000                0.6088889 
    ## 
    ## 
    ## Dependent var: Speed2ndEntr 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -2.3113, df = 9.3935, p-value = 0.04499
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -8.0750631 -0.1124369
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                  1.27625                  5.37000 
    ## 
    ## 
    ## Dependent var: TimeTarget 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -6.131, df = 12.076, p-value = 4.953e-05
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -23.58683 -11.22439
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                  4.22450                 21.63011 
    ## 
    ## 
    ## Dependent var: pTimeTarget 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -4.9141, df = 12.776, p-value = 0.0002976
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.05078248 -0.01972863
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##               0.01100000               0.04625556 
    ## 
    ## 
    ## Dependent var: pTimeCCW 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 0.33985, df = 14.976, p-value = 0.7387
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.2431431  0.3353736
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                0.4144375                0.3683222 
    ## 
    ## 
    ## Dependent var: pTimeOPP 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 2.6446, df = 11.683, p-value = 0.02181
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  0.04646268 0.48867621
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                0.4683250                0.2007556 
    ## 
    ## 
    ## Dependent var: pTimeCW 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -2.1871, df = 9.1991, p-value = 0.05589
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.565539722  0.008611944
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                0.1062250                0.3846889 
    ## 
    ## 
    ## Dependent var: RayleigLength 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 0.96335, df = 11.51, p-value = 0.3552
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.09472544  0.24361433
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                0.6700000                0.5955556 
    ## 
    ## 
    ## Dependent var: RayleigAngle 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -0.34369, df = 10.894, p-value = 0.7376
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -80.97483  59.12372
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 145.9400                 156.8656 
    ## 
    ## 
    ## Dependent var: PolarAvgVal 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 2.5868, df = 14.922, p-value = 0.0207
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##   6.408734 66.557932
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 228.1000                 191.6167 
    ## 
    ## 
    ## Dependent var: PolarSdVal 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -3.5564, df = 10.106, p-value = 0.005127
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -68.07182 -15.67680
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 53.39125                 95.26556 
    ## 
    ## 
    ## Dependent var: PolarMinVal 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -2.6199, df = 8.0648, p-value = 0.03045
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.0012682916 -0.0000817084
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                  2.5e-05                  7.0e-04 
    ## 
    ## 
    ## Dependent var: PolarMinBin 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -0.8282, df = 13.006, p-value = 0.4225
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -188.43861   83.99416
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                  80.0000                 132.2222 
    ## 
    ## 
    ## Dependent var: Min50.RngLoBin 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -0.2737, df = 12.336, p-value = 0.7888
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -81.92026  63.58693
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 167.5000                 176.6667 
    ## 
    ## 
    ## Dependent var: Min50.RngHiBin 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -0.061894, df = 11.224, p-value = 0.9517
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -81.05321  76.60876
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 120.0000                 122.2222 
    ## 
    ## 
    ## Dependent var: PolarMaxVal 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 0.089917, df = 14.849, p-value = 0.9296
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.02591349  0.02819405
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                0.1011625                0.1000222 
    ## 
    ## 
    ## Dependent var: PolarMaxBin 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -0.54112, df = 10.943, p-value = 0.5993
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -103.5134   62.6801
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 136.2500                 156.6667 
    ## 
    ## 
    ## Dependent var: Max50.RngLoBin 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -0.23477, df = 11.078, p-value = 0.8187
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -84.95305  68.56416
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 106.2500                 114.4444 
    ## 
    ## 
    ## Dependent var: Max50.RngHiBin 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -0.45091, df = 11.962, p-value = 0.6601
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -88.31571  58.03793
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 183.7500                 198.8889 
    ## 
    ## 
    ## Dependent var: AnnularMinVal 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 0.41672, df = 14.069, p-value = 0.6832
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.002716933  0.004028044
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##              0.003100000              0.002444444 
    ## 
    ## 
    ## Dependent var: AnnularMinBin 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -0.16031, df = 14.456, p-value = 0.8749
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -8.105574  6.975018
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 10.71250                 11.27778 
    ## 
    ## 
    ## Dependent var: AnnularMaxVal 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -0.94398, df = 14.133, p-value = 0.361
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.11920758  0.04629924
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                0.4272125                0.4636667 
    ## 
    ## 
    ## Dependent var: AnnularMaxBin 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 0.40551, df = 13.899, p-value = 0.6913
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.8227699  1.2061033
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 16.92500                 16.73333 
    ## 
    ## 
    ## Dependent var: AnnularAvg 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 0.12252, df = 13.961, p-value = 0.9042
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.6604011  0.7404011
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                    16.29                    16.25 
    ## 
    ## 
    ## Dependent var: AnnularSd 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 0.23896, df = 12.607, p-value = 0.815
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -3.138132  3.915910
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 11.07000                 10.68111 
    ## 
    ## 
    ## Dependent var: AnnularSkewnes 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -0.0027527, df = 13.842, p-value = 0.9978
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -1.193181  1.190125
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 2.046250                 2.047778 
    ## 
    ## 
    ## Dependent var: AnnularKurtosis 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -0.53377, df = 14.684, p-value = 0.6015
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -8.889479  5.334201
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 11.07125                 12.84889 
    ## 
    ## 
    ## Dependent var: Speed1 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 2.2403, df = 9.6588, p-value = 0.04988
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  7.094273e-06 2.257782e-02
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##               0.02562134               0.01432888 
    ## 
    ## 
    ## Dependent var: Speed2 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = -3.9218, df = 12.807, p-value = 0.001801
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.016332075 -0.004718381
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##               0.02610536               0.03663059 
    ## 
    ## 
    ## Dependent var: Time1stEntrLog 
    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  slim4[[y]] by slim3$APA2
    ## t = 6.4273, df = 13.627, p-value = 1.797e-05
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  2.798961 5.613415
    ## sample estimates:
    ## mean in group consistent   mean in group conflict 
    ##                 5.255859                 1.049671

``` r
for(y in names(slim4)){
  ymod <- leveneTest(slim4[[y]] ~ slim3$APA2 )
  cat(paste('\nDependent var:', y, '\n'))
  print(ymod)
}
```

    ## 
    ## Dependent var: SdevSpeedArena 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.6373 0.4372
    ##       15               
    ## 
    ## Dependent var: Linearity.Arena. 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  1.0011 0.3329
    ##       15               
    ## 
    ## Dependent var: NumEntrances 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  3.5807 0.07792 .
    ##       15                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Time1stEntr 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value   Pr(>F)   
    ## group  1  10.719 0.005125 **
    ##       15                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Path1stEntr 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value   Pr(>F)   
    ## group  1  15.271 0.001399 **
    ##       15                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Speed1stEntr.cm.s. 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0635 0.8044
    ##       15               
    ## 
    ## Dependent var: Dist1stEntr.m. 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  2.1913 0.1595
    ##       15               
    ## 
    ## Dependent var: NumShock 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  3.1029 0.09852 .
    ##       15                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: MaxTimeAvoid 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.1734  0.683
    ##       15               
    ## 
    ## Dependent var: Time2ndEntr 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value    Pr(>F)    
    ## group  1  23.588 0.0002093 ***
    ##       15                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Path2ndEntr 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value    Pr(>F)    
    ## group  1  22.408 0.0002664 ***
    ##       15                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Speed2ndEntr 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  2.5895 0.1284
    ##       15               
    ## 
    ## Dependent var: TimeTarget 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  3.8701 0.06793 .
    ##       15                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: pTimeTarget 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  2.4332 0.1396
    ##       15               
    ## 
    ## Dependent var: pTimeCCW 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.6844  0.421
    ##       15               
    ## 
    ## Dependent var: pTimeOPP 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.5979 0.4514
    ##       15               
    ## 
    ## Dependent var: pTimeCW 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  5.4947 0.03326 *
    ##       15                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: RayleigLength 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  1.2091 0.2888
    ##       15               
    ## 
    ## Dependent var: RayleigAngle 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  4.7896 0.04487 *
    ##       15                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: PolarAvgVal 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1   0.033 0.8583
    ##       15               
    ## 
    ## Dependent var: PolarSdVal 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  2.7169 0.1201
    ##       15               
    ## 
    ## Dependent var: PolarMinVal 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  6.9164 0.01894 *
    ##       15                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: PolarMinBin 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1   0.011  0.918
    ##       15               
    ## 
    ## Dependent var: Min50.RngLoBin 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  4.1584 0.05945 .
    ##       15                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Min50.RngHiBin 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1    1.82 0.1973
    ##       15               
    ## 
    ## Dependent var: PolarMaxVal 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.2287 0.6394
    ##       15               
    ## 
    ## Dependent var: PolarMaxBin 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  6.0328 0.02671 *
    ##       15                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Max50.RngLoBin 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  2.8794 0.1104
    ##       15               
    ## 
    ## Dependent var: Max50.RngHiBin 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  6.5891 0.02147 *
    ##       15                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: AnnularMinVal 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.2773 0.6062
    ##       15               
    ## 
    ## Dependent var: AnnularMinBin 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0054 0.9426
    ##       15               
    ## 
    ## Dependent var: AnnularMaxVal 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0502 0.8257
    ##       15               
    ## 
    ## Dependent var: AnnularMaxBin 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.4136 0.5299
    ##       15               
    ## 
    ## Dependent var: AnnularAvg 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.1464 0.7074
    ##       15               
    ## 
    ## Dependent var: AnnularSd 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0508 0.8248
    ##       15               
    ## 
    ## Dependent var: AnnularSkewnes 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1   0.685 0.4208
    ##       15               
    ## 
    ## Dependent var: AnnularKurtosis 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  1.3212 0.2684
    ##       15               
    ## 
    ## Dependent var: Speed1 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value    Pr(>F)    
    ## group  1  19.839 0.0004639 ***
    ##       15                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Speed2 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  3.0982 0.09875 .
    ##       15                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Time1stEntrLog 
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  3.1531 0.09607 .
    ##       15                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# *** Speed1, Path2ndEntr, Time2ndEntr, Path2ndEntr, Time2ndEntr
# ** Path1stEntr, Time1stEntr
# * Max50.RngHiBin ,  PolarMaxBin  , PolarMinVal, RayleigAngle, pTimeCW, pTimeOPP,
# . Time1stEntrLog, Speed2, Min50.RngLoBin , TimeTarget, NumShock, NumEntrances
# AnnularKurtosis, AnnularSkewnes, AnnularSd, AnnularMaxBin, AnnularMaxVal, AnnularMinBin, AnnularMinVal, Max50.RngLoBin, RayleigLength PolarMaxVal, Min50.RngHiBin , PolarMinBin, PolarSdVal, PolarAvgVal, RayleigLength, pTimeTarget, Speed2ndEntr, MaxTimeAvoid, Dist1stEntr.m, Speed1stEntr.cm.s, Linearity.Arena, SdevSpeedArena

par(mfrow=c(3,3))
hist(slim4$Time1stEntrLog)
boxplot(slim4$Time1stEntrLog ~ slim3$APA2)

for(y in names(slim4)){
  ymod <- hist(slim4[[y]],
               main = y, xlab = NULL,
               col = c())
}
```

![](../figures/01_behavior/unnamed-chunk-1-1.png)![](../figures/01_behavior/unnamed-chunk-1-2.png)![](../figures/01_behavior/unnamed-chunk-1-3.png)![](../figures/01_behavior/unnamed-chunk-1-4.png)

``` r
for(y in names(slim4)){
  ymod <- boxplot(slim4[[y]] ~ slim3$APA2,
               main = y,
               xlab = NULL)
}
```

![](../figures/01_behavior/unnamed-chunk-1-5.png)![](../figures/01_behavior/unnamed-chunk-1-6.png)![](../figures/01_behavior/unnamed-chunk-1-7.png)![](../figures/01_behavior/unnamed-chunk-1-8.png)![](../figures/01_behavior/unnamed-chunk-1-9.png)![](../figures/01_behavior/unnamed-chunk-1-10.png)

Repeated measures anova
-----------------------

ges: Generalized Eta-Squared measure of effect size (see in references below: Bakeman, 2005).
Sphericity: assumption of the test

``` r
# sample sizes
slim1 <- behavior[,c(15,16,6,20:59)] # drop frivolous columns
slim2 <- slim1 %>% 
  filter(TrainGroup == "trained") 
slim3 <- slim2[,c(1:3)] # extract experimental design variables
slim4 <- slim2[,c(4:43)]

for(y in names(slim4)){
  ymod<- summary(aov(slim4[[y]] ~ slim3$APA2 * slim3$TrainSessionCombo ))
  cat(paste('\nDependent var:', y, '\n'))
  print(ymod)
}
```

    ## 
    ## Dependent var: SdevSpeedArena 
    ##                                     Df Sum Sq Mean Sq F value Pr(>F)    
    ## slim3$APA2                           1  0.359  0.3588   2.119  0.148    
    ## slim3$TrainSessionCombo              8 22.373  2.7966  16.516 <2e-16 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8  1.298  0.1623   0.958  0.471    
    ## Residuals                          135 22.860  0.1693                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Linearity.Arena. 
    ##                                     Df Sum Sq Mean Sq F value Pr(>F)    
    ## slim3$APA2                           1 0.0132 0.01319   6.963 0.0093 ** 
    ## slim3$TrainSessionCombo              8 0.4562 0.05702  30.097 <2e-16 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8 0.0155 0.00193   1.020 0.4240    
    ## Residuals                          135 0.2558 0.00189                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: NumEntrances 
    ##                                     Df Sum Sq Mean Sq F value   Pr(>F)    
    ## slim3$APA2                           1    148   147.6   5.838 0.017018 *  
    ## slim3$TrainSessionCombo              8  11253  1406.7  55.653  < 2e-16 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8    813   101.6   4.019 0.000262 ***
    ## Residuals                          135   3412    25.3                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Time1stEntr 
    ##                                     Df  Sum Sq Mean Sq F value   Pr(>F)
    ## slim3$APA2                           1   58341   58341   1.677    0.198
    ## slim3$TrainSessionCombo              8 1736981  217123   6.242 7.38e-07
    ## slim3$APA2:slim3$TrainSessionCombo   8  406492   50812   1.461    0.177
    ## Residuals                          135 4695862   34784                 
    ##                                       
    ## slim3$APA2                            
    ## slim3$TrainSessionCombo            ***
    ## slim3$APA2:slim3$TrainSessionCombo    
    ## Residuals                             
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Path1stEntr 
    ##                                     Df Sum Sq Mean Sq F value  Pr(>F)    
    ## slim3$APA2                           1   21.3   21.29   0.929   0.337    
    ## slim3$TrainSessionCombo              8 1124.5  140.57   6.132 9.8e-07 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8  293.5   36.69   1.601   0.130    
    ## Residuals                          135 3094.5   22.92                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Speed1stEntr.cm.s. 
    ##                                     Df Sum Sq Mean Sq F value   Pr(>F)    
    ## slim3$APA2                           1    1.3    1.33   0.060    0.807    
    ## slim3$TrainSessionCombo              8 1141.2  142.64   6.465 4.15e-07 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8   72.7    9.09   0.412    0.912    
    ## Residuals                          135 2978.8   22.07                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Dist1stEntr.m. 
    ##                                     Df Sum Sq Mean Sq F value  Pr(>F)    
    ## slim3$APA2                           1  0.394  0.3935   5.855 0.01686 *  
    ## slim3$TrainSessionCombo              8 11.633  1.4542  21.637 < 2e-16 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8  1.841  0.2301   3.424 0.00129 ** 
    ## Residuals                          135  9.073  0.0672                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: NumShock 
    ##                                     Df Sum Sq Mean Sq F value  Pr(>F)    
    ## slim3$APA2                           1    399     399   7.862 0.00579 ** 
    ## slim3$TrainSessionCombo              8  38859    4857  95.757 < 2e-16 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8   1188     149   2.928 0.00478 ** 
    ## Residuals                          135   6848      51                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: MaxTimeAvoid 
    ##                                     Df  Sum Sq Mean Sq F value   Pr(>F)
    ## slim3$APA2                           1   71883   71883   3.340   0.0698
    ## slim3$TrainSessionCombo              8 1641566  205196   9.534 2.12e-10
    ## slim3$APA2:slim3$TrainSessionCombo   8  178024   22253   1.034   0.4138
    ## Residuals                          135 2905642   21523                 
    ##                                       
    ## slim3$APA2                         .  
    ## slim3$TrainSessionCombo            ***
    ## slim3$APA2:slim3$TrainSessionCombo    
    ## Residuals                             
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Time2ndEntr 
    ##                                     Df  Sum Sq Mean Sq F value   Pr(>F)
    ## slim3$APA2                           1  128822  128822   3.032   0.0839
    ## slim3$TrainSessionCombo              8 3124476  390559   9.193 4.75e-10
    ## slim3$APA2:slim3$TrainSessionCombo   8  682401   85300   2.008   0.0500
    ## Residuals                          135 5735214   42483                 
    ##                                       
    ## slim3$APA2                         .  
    ## slim3$TrainSessionCombo            ***
    ## slim3$APA2:slim3$TrainSessionCombo *  
    ## Residuals                             
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Path2ndEntr 
    ##                                     Df Sum Sq Mean Sq F value   Pr(>F)    
    ## slim3$APA2                           1     46   45.62   1.573   0.2119    
    ## slim3$TrainSessionCombo              8   1967  245.92   8.480 2.65e-09 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8    511   63.91   2.204   0.0308 *  
    ## Residuals                          135   3915   29.00                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Speed2ndEntr 
    ##                                     Df Sum Sq Mean Sq F value   Pr(>F)    
    ## slim3$APA2                           1   94.4   94.40   4.410   0.0376 *  
    ## slim3$TrainSessionCombo              8 1427.9  178.48   8.338 3.75e-09 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8  239.6   29.95   1.399   0.2023    
    ## Residuals                          135 2889.9   21.41                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: TimeTarget 
    ##                                     Df Sum Sq Mean Sq F value Pr(>F)    
    ## slim3$APA2                           1    736     736   5.968 0.0159 *  
    ## slim3$TrainSessionCombo              8 143803   17975 145.845 <2e-16 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8   1823     228   1.849 0.0732 .  
    ## Residuals                          135  16639     123                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: pTimeTarget 
    ##                                     Df Sum Sq Mean Sq F value Pr(>F)    
    ## slim3$APA2                           1 0.0022 0.00215   3.430 0.0662 .  
    ## slim3$TrainSessionCombo              8 0.9194 0.11493 183.304 <2e-16 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8 0.0060 0.00075   1.193 0.3075    
    ## Residuals                          135 0.0846 0.00063                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: pTimeCCW 
    ##                                     Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2                           1  0.000 0.00019   0.004  0.951
    ## slim3$TrainSessionCombo              8  0.423 0.05293   1.064  0.392
    ## slim3$APA2:slim3$TrainSessionCombo   8  0.340 0.04250   0.854  0.557
    ## Residuals                          135  6.718 0.04976               
    ## 
    ## Dependent var: pTimeOPP 
    ##                                     Df Sum Sq Mean Sq F value  Pr(>F)    
    ## slim3$APA2                           1  0.114 0.11426   3.266   0.073 .  
    ## slim3$TrainSessionCombo              8  1.686 0.21078   6.025 1.3e-06 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8  0.332 0.04145   1.185   0.313    
    ## Residuals                          135  4.723 0.03499                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: pTimeCW 
    ##                                     Df Sum Sq Mean Sq F value   Pr(>F)    
    ## slim3$APA2                           1 0.0932 0.09321   4.590 0.033950 *  
    ## slim3$TrainSessionCombo              8 0.6046 0.07557   3.722 0.000582 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8 0.4004 0.05005   2.465 0.015932 *  
    ## Residuals                          135 2.7414 0.02031                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: RayleigLength 
    ##                                     Df Sum Sq Mean Sq F value Pr(>F)    
    ## slim3$APA2                           1  0.002  0.0022   0.133  0.716    
    ## slim3$TrainSessionCombo              8  5.418  0.6772  41.100 <2e-16 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8  0.152  0.0190   1.153  0.333    
    ## Residuals                          135  2.224  0.0165                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: RayleigAngle 
    ##                                     Df Sum Sq Mean Sq F value Pr(>F)  
    ## slim3$APA2                           1    344     344   0.135 0.7144  
    ## slim3$TrainSessionCombo              8  38237    4780   1.869 0.0698 .
    ## slim3$APA2:slim3$TrainSessionCombo   8  11154    1394   0.545 0.8205  
    ## Residuals                          135 345189    2557                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: PolarAvgVal 
    ##                                     Df Sum Sq Mean Sq F value   Pr(>F)    
    ## slim3$APA2                           1  64074   64074   92.17  < 2e-16 ***
    ## slim3$TrainSessionCombo              8 104175   13022   18.73  < 2e-16 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8  84282   10535   15.15 1.02e-15 ***
    ## Residuals                          135  93846     695                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: PolarSdVal 
    ##                                     Df Sum Sq Mean Sq F value   Pr(>F)    
    ## slim3$APA2                           1   5254    5254  14.087 0.000258 ***
    ## slim3$TrainSessionCombo              8  38612    4827  12.942 9.96e-14 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8  14326    1791   4.802 3.25e-05 ***
    ## Residuals                          135  50347     373                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: PolarMinVal 
    ##                                     Df   Sum Sq   Mean Sq F value Pr(>F)
    ## slim3$APA2                           1 0.000009 0.0000091   2.887 0.0916
    ## slim3$TrainSessionCombo              8 0.003313 0.0004142 130.932 <2e-16
    ## slim3$APA2:slim3$TrainSessionCombo   8 0.000040 0.0000050   1.571 0.1390
    ## Residuals                          135 0.000427 0.0000032               
    ##                                       
    ## slim3$APA2                         .  
    ## slim3$TrainSessionCombo            ***
    ## slim3$APA2:slim3$TrainSessionCombo    
    ## Residuals                             
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: PolarMinBin 
    ##                                     Df  Sum Sq Mean Sq F value   Pr(>F)
    ## slim3$APA2                           1   33884   33884   2.112 0.148463
    ## slim3$TrainSessionCombo              8  465482   58185   3.627 0.000749
    ## slim3$APA2:slim3$TrainSessionCombo   8  302140   37767   2.354 0.021112
    ## Residuals                          135 2165800   16043                 
    ##                                       
    ## slim3$APA2                            
    ## slim3$TrainSessionCombo            ***
    ## slim3$APA2:slim3$TrainSessionCombo *  
    ## Residuals                             
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Min50.RngLoBin 
    ##                                     Df Sum Sq Mean Sq F value Pr(>F)  
    ## slim3$APA2                           1   1617    1617   0.636 0.4266  
    ## slim3$TrainSessionCombo              8  39605    4951   1.947 0.0579 .
    ## slim3$APA2:slim3$TrainSessionCombo   8  11219    1402   0.552 0.8157  
    ## Residuals                          135 343247    2543                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Min50.RngHiBin 
    ##                                     Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2                           1   2343    2343   0.579  0.448
    ## slim3$TrainSessionCombo              8  51499    6437   1.591  0.133
    ## slim3$APA2:slim3$TrainSessionCombo   8  21186    2648   0.655  0.730
    ## Residuals                          135 546094    4045               
    ## 
    ## Dependent var: PolarMaxVal 
    ##                                     Df  Sum Sq  Mean Sq F value   Pr(>F)
    ## slim3$APA2                           1 0.00102 0.001023   2.044    0.155
    ## slim3$TrainSessionCombo              8 0.06305 0.007881  15.746 3.14e-16
    ## slim3$APA2:slim3$TrainSessionCombo   8 0.00257 0.000321   0.641    0.742
    ## Residuals                          135 0.06757 0.000500                 
    ##                                       
    ## slim3$APA2                            
    ## slim3$TrainSessionCombo            ***
    ## slim3$APA2:slim3$TrainSessionCombo    
    ## Residuals                             
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: PolarMaxBin 
    ##                                     Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2                           1   1431    1431   0.368  0.545
    ## slim3$TrainSessionCombo              8  19899    2487   0.639  0.744
    ## slim3$APA2:slim3$TrainSessionCombo   8  10347    1293   0.332  0.952
    ## Residuals                          135 525540    3893               
    ## 
    ## Dependent var: Max50.RngLoBin 
    ##                                     Df Sum Sq Mean Sq F value Pr(>F)
    ## slim3$APA2                           1     41      41   0.011  0.917
    ## slim3$TrainSessionCombo              8  39835    4979   1.333  0.232
    ## slim3$APA2:slim3$TrainSessionCombo   8   5956     745   0.199  0.991
    ## Residuals                          135 504321    3736               
    ## 
    ## Dependent var: Max50.RngHiBin 
    ##                                     Df Sum Sq Mean Sq F value Pr(>F)  
    ## slim3$APA2                           1   1297    1297   0.489 0.4856  
    ## slim3$TrainSessionCombo              8  48847    6106   2.302 0.0241 *
    ## slim3$APA2:slim3$TrainSessionCombo   8   9156    1144   0.431 0.9005  
    ## Residuals                          135 358147    2653                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: AnnularMinVal 
    ##                                     Df   Sum Sq   Mean Sq F value   Pr(>F)
    ## slim3$APA2                           1 0.000047 4.668e-05   0.879 0.350206
    ## slim3$TrainSessionCombo              8 0.001676 2.094e-04   3.943 0.000322
    ## slim3$APA2:slim3$TrainSessionCombo   8 0.000157 1.961e-05   0.369 0.935264
    ## Residuals                          135 0.007171 5.312e-05                 
    ##                                       
    ## slim3$APA2                            
    ## slim3$TrainSessionCombo            ***
    ## slim3$APA2:slim3$TrainSessionCombo    
    ## Residuals                             
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: AnnularMinBin 
    ##                                     Df Sum Sq Mean Sq F value   Pr(>F)    
    ## slim3$APA2                           1      7    7.47   0.225    0.636    
    ## slim3$TrainSessionCombo              8   1618  202.21   6.077 1.13e-06 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8    204   25.53   0.767    0.632    
    ## Residuals                          135   4492   33.28                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: AnnularMaxVal 
    ##                                     Df Sum Sq  Mean Sq F value  Pr(>F)    
    ## slim3$APA2                           1 0.0108 0.010841   1.778   0.185    
    ## slim3$TrainSessionCombo              8 0.2473 0.030908   5.069 1.6e-05 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8 0.0265 0.003314   0.544   0.822    
    ## Residuals                          135 0.8231 0.006097                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: AnnularMaxBin 
    ##                                     Df Sum Sq Mean Sq F value Pr(>F)    
    ## slim3$APA2                           1   0.72   0.722   0.946  0.333    
    ## slim3$TrainSessionCombo              8 135.40  16.925  22.161 <2e-16 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8   3.29   0.411   0.538  0.826    
    ## Residuals                          135 103.10   0.764                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: AnnularAvg 
    ##                                     Df Sum Sq Mean Sq F value Pr(>F)    
    ## slim3$APA2                           1   1.31   1.307   3.308 0.0712 .  
    ## slim3$TrainSessionCombo              8 159.57  19.946  50.474 <2e-16 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8   1.17   0.146   0.370 0.9346    
    ## Residuals                          135  53.35   0.395                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: AnnularSd 
    ##                                     Df Sum Sq Mean Sq F value   Pr(>F)    
    ## slim3$APA2                           1   14.9   14.92   1.364    0.245    
    ## slim3$TrainSessionCombo              8 1219.1  152.39  13.928 1.25e-14 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8   84.1   10.52   0.961    0.469    
    ## Residuals                          135 1477.0   10.94                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: AnnularSkewnes 
    ##                                     Df Sum Sq Mean Sq F value  Pr(>F)   
    ## slim3$APA2                           1   3.07  3.0663   4.535 0.03503 * 
    ## slim3$TrainSessionCombo              8  14.32  1.7901   2.647 0.00996 **
    ## slim3$APA2:slim3$TrainSessionCombo   8   5.00  0.6256   0.925 0.49791   
    ## Residuals                          135  91.28  0.6762                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: AnnularKurtosis 
    ##                                     Df Sum Sq Mean Sq F value  Pr(>F)   
    ## slim3$APA2                           1    498   498.1   9.342 0.00270 **
    ## slim3$TrainSessionCombo              8   1392   174.0   3.264 0.00197 **
    ## slim3$APA2:slim3$TrainSessionCombo   8    373    46.7   0.875 0.53929   
    ## Residuals                          135   7198    53.3                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Speed1 
    ##                                     Df  Sum Sq   Mean Sq F value   Pr(>F)
    ## slim3$APA2                           1 0.00005 0.0000501   0.191    0.663
    ## slim3$TrainSessionCombo              8 0.01253 0.0015663   5.973 1.48e-06
    ## slim3$APA2:slim3$TrainSessionCombo   8 0.00125 0.0001559   0.594    0.781
    ## Residuals                          135 0.03540 0.0002622                 
    ##                                       
    ## slim3$APA2                            
    ## slim3$TrainSessionCombo            ***
    ## slim3$APA2:slim3$TrainSessionCombo    
    ## Residuals                             
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Speed2 
    ##                                     Df   Sum Sq  Mean Sq F value Pr(>F)
    ## slim3$APA2                           1 0.000502 0.000502   5.549 0.0199
    ## slim3$TrainSessionCombo              8 0.025590 0.003199  35.391 <2e-16
    ## slim3$APA2:slim3$TrainSessionCombo   8 0.000773 0.000097   1.069 0.3886
    ## Residuals                          135 0.012202 0.000090               
    ##                                       
    ## slim3$APA2                         *  
    ## slim3$TrainSessionCombo            ***
    ## slim3$APA2:slim3$TrainSessionCombo    
    ## Residuals                             
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dependent var: Time1stEntrLog 
    ##                                     Df Sum Sq Mean Sq F value   Pr(>F)    
    ## slim3$APA2                           1   1.43    1.43   0.659    0.418    
    ## slim3$TrainSessionCombo              8 290.39   36.30  16.756  < 2e-16 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8  87.18   10.90   5.031 1.77e-05 ***
    ## Residuals                          135 292.45    2.17                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(slim4$NumEntrances ~ slim3$APA2 * slim3$TrainSessionCombo))
```

    ##                                     Df Sum Sq Mean Sq F value   Pr(>F)    
    ## slim3$APA2                           1    148   147.6   5.838 0.017018 *  
    ## slim3$TrainSessionCombo              8  11253  1406.7  55.653  < 2e-16 ***
    ## slim3$APA2:slim3$TrainSessionCombo   8    813   101.6   4.019 0.000262 ***
    ## Residuals                          135   3412    25.3                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# sample sizes
behavior %>% filter(TrainSession == "Retention") %>% select(APA2)  %>%  summary()
```

    ##                APA2  
    ##  yoked-consistent:8  
    ##  consistent      :8  
    ##  yoked-conflict  :9  
    ##  conflict        :9

``` r
library(ez)
ezANOVA(trained, dv = NumEntrances, wid = ID, detailed = F,
        within = TrainSessionCombo, between = APA2)
```

    ## Warning: Data is unbalanced (unequal N per group). Make sure you specified
    ## a well-considered value for the type argument to ezANOVA().

    ## $ANOVA
    ##                   Effect DFn DFd         F            p p<.05        ges
    ## 2                   APA2   1  15  2.312766 1.491132e-01       0.04145403
    ## 3      TrainSessionCombo   8 120 68.753449 3.402258e-41     * 0.76732992
    ## 4 APA2:TrainSessionCombo   8 120  4.965579 2.514143e-05     * 0.19236717
    ## 
    ## $`Mauchly's Test for Sphericity`
    ##                   Effect            W            p p<.05
    ## 3      TrainSessionCombo 8.316293e-05 8.646959e-10     *
    ## 4 APA2:TrainSessionCombo 8.316293e-05 8.646959e-10     *
    ## 
    ## $`Sphericity Corrections`
    ##                   Effect       GGe        p[GG] p[GG]<.05      HFe
    ## 3      TrainSessionCombo 0.3996539 8.396272e-18         * 0.520595
    ## 4 APA2:TrainSessionCombo 0.3996539 3.740327e-03         * 0.520595
    ##          p[HF] p[HF]<.05
    ## 3 1.574820e-22         *
    ## 4 1.336868e-03         *

Standard vizualization of mean avoidance beavior
------------------------------------------------

First, I visualze the group mean and standard error for the time it takes before an individual mouse enters the spatial region marked "schock zone" or equivilent.

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
TukeyHSD((aov(PC3 ~ APA2, data=scoresdf)), which = "APA2") 
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC3 ~ APA2, data = scoresdf)
    ## 
    ## $APA2
    ##                                      diff        lwr      upr     p adj
    ## consistent-yoked-consistent     -2.517549 -7.9556882 2.920591 0.5954357
    ## yoked-conflict-yoked-consistent -1.331357 -6.6162791 3.953564 0.9020015
    ## conflict-yoked-consistent        2.128797 -3.1561244 7.413719 0.6950837
    ## yoked-conflict-consistent        1.186191 -4.0987304 6.471113 0.9280078
    ## conflict-consistent              4.646346 -0.6385758 9.931268 0.1006931
    ## conflict-yoked-conflict          3.460155 -1.6669725 8.587282 0.2773302

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
pca16 <- ggplot(scoresdf, aes(PC1,PC2, color=APA2)) +
    geom_point(size=3, alpha = 0.7) +
    xlab(paste0("PC 1: ", percent[1],"% variance")) +
    ylab(paste0("PC 2: ", percent[2],"% variance")) +
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

``` r
pca12 <- ggplot(scoresdf, aes(PC1,PC2, color=APA2)) +
    geom_point(size=3, alpha = 0.7) +
    xlab(paste0("PC 1: ", percent[1],"% variance")) +
    ylab(paste0("PC 2: ", percent[2],"% variance")) +
    stat_ellipse(level = 0.95, (aes(color=APA2)),size=0.25) + 
    scale_colour_manual(values=c(colorvalAPA00)) + 
    theme_cowplot(font_size = 8, line_size = 0.25) +
      theme(legend.position="none") +
  facet_wrap(~wrap)
pca12
```

![](../figures/01_behavior/PCAplots-3.png)

``` r
library(factoextra)
res.pca <- prcomp(behaviormatrix, scale = TRUE)
fviz_eig(res.pca)
```

![](../figures/01_behavior/PCAplots-4.png)

``` r
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             select.var = list(contrib = 10))
```

![](../figures/01_behavior/PCAplots-5.png)

``` r
fviz_pca_biplot(res.pca, label ="var")
```

![](../figures/01_behavior/PCAplots-6.png)

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
