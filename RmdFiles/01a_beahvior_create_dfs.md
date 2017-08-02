``` r
## load libraries 
#library(magrittr) ## to use the weird pipe
library(tidyr) ## for respahing data
```

    ## Warning: package 'tidyr' was built under R version 3.3.2

``` r
library(plyr) ## for renmaing factors
library(dplyr) ## for filtering and selecting rows
library(reshape2) ## for melting dataframe

## load functions 
source("functions_behavior.R")
```

``` r
## read intermediate data (raw data from video tracker program analyzed in matlab)
behavior <- read.csv("../data/01_behaviordata.csv", header = T)
## relevel factors
behavior$APA <- factor(behavior$APA, levels = c("Yoked", "Same", "Conflict"))
behavior$APA2 <- factor(behavior$APA2, levels = c("YokedSame", "YokedConflict","Same", "Conflict"))
str(behavior$APA)
```

    ##  Factor w/ 3 levels "Yoked","Same",..: 2 2 2 2 2 2 2 2 2 2 ...

``` r
## revalue factors
levels(behavior$APA) <- c("control","consistent","conflict")
levels(behavior$APA2) <- c("yoked-consistent","yoked-conflict","consistent","conflict")

## log transformation 
behavior$Time1stEntrLog <- log(behavior$Time1stEntr) 
summary(behavior)
```

    ##        ID           Year      Genotype TrainProtocol        TrainSequence
    ##  15140A :  9   Min.   :2015   WT:306   3-day:306     train-conflict:162  
    ##  15140B :  9   1st Qu.:2015                          train-train   :144  
    ##  15140C :  9   Median :2015                                              
    ##  15140D :  9   Mean   :2015                                              
    ##  15141C :  9   3rd Qu.:2015                                              
    ##  15141D :  9   Max.   :2015                                              
    ##  (Other):252                                                             
    ##    TrainGroup       Day           TrainSession ShockOnOff PairedPartner
    ##  trained:153   Min.   :1.000   Hab      : 34   Off: 68    15140A :  9  
    ##  yoked  :153   1st Qu.:1.000   Retention: 34   On :238    15140B :  9  
    ##                Median :2.000   Retest   : 34              15140C :  9  
    ##                Mean   :1.676   T1       : 34              15140D :  9  
    ##                3rd Qu.:2.000   T2       : 34              15141C :  9  
    ##                Max.   :4.000   T3       : 34              15141D :  9  
    ##                                (Other)  :102              (Other):252  
    ##  Experimenter       Housing    TestLocation         APA     
    ##  Maddy:306    AnimalRoom:306   MBL:306      control   :153  
    ##                                             consistent: 72  
    ##                                             conflict  : 81  
    ##                                                             
    ##                                                             
    ##                                                             
    ##                                                             
    ##                APA2    TrainSessionCombo              pair1    
    ##  yoked-consistent:72   Hab      : 34     15140A_Hab      :  1  
    ##  yoked-conflict  :81   Retention: 34     15140A_Retention:  1  
    ##  consistent      :72   Retest   : 34     15140A_Retest   :  1  
    ##  conflict        :81   T1       : 34     15140A_T1       :  1  
    ##                        T2       : 34     15140A_T2       :  1  
    ##                        T3       : 34     15140A_T3       :  1  
    ##                        (Other)  :102     (Other)         :300  
    ##               pair2     TrainSessionComboNum SdevSpeedArena 
    ##  15140A_Hab      :  1   Min.   :1            Min.   :0.570  
    ##  15140A_Retention:  1   1st Qu.:3            1st Qu.:2.170  
    ##  15140A_Retest   :  1   Median :5            Median :2.525  
    ##  15140A_T1       :  1   Mean   :5            Mean   :2.550  
    ##  15140A_T2       :  1   3rd Qu.:7            3rd Qu.:3.000  
    ##  15140A_T3       :  1   Max.   :9            Max.   :4.730  
    ##  (Other)         :300                                       
    ##  Linearity.Arena.  NumEntrances    Time1stEntr       Path1stEntr    
    ##  Min.   :0.1360   Min.   : 0.00   Min.   :  0.330   Min.   : 0.000  
    ##  1st Qu.:0.2764   1st Qu.: 4.00   1st Qu.:  9.178   1st Qu.: 0.220  
    ##  Median :0.3191   Median :12.00   Median : 24.050   Median : 0.615  
    ##  Mean   :0.3261   Mean   :13.01   Mean   : 96.018   Mean   : 2.475  
    ##  3rd Qu.:0.3701   3rd Qu.:18.00   3rd Qu.: 67.025   3rd Qu.: 1.785  
    ##  Max.   :0.5556   Max.   :42.00   Max.   :600.000   Max.   :16.510  
    ##                                                                     
    ##  Speed1stEntr.cm.s. Dist1stEntr.m.     NumShock       MaxTimeAvoid   
    ##  Min.   :-1.000     Min.   :0.000   Min.   :  0.00   Min.   : 30.00  
    ##  1st Qu.: 1.610     1st Qu.:0.240   1st Qu.:  5.00   1st Qu.: 73.75  
    ##  Median : 2.040     Median :0.765   Median : 33.50   Median :122.00  
    ##  Mean   : 5.084     Mean   :0.738   Mean   : 32.67   Mean   :204.45  
    ##  3rd Qu.: 7.655     3rd Qu.:1.140   3rd Qu.: 55.00   3rd Qu.:312.75  
    ##  Max.   :31.730     Max.   :2.420   Max.   :113.00   Max.   :599.00  
    ##                                                                      
    ##   Time2ndEntr      Path2ndEntr      Speed2ndEntr      TimeTarget     
    ##  Min.   :  4.83   Min.   : 0.080   Min.   :-1.000   Min.   :  0.000  
    ##  1st Qu.: 35.04   1st Qu.: 1.040   1st Qu.: 1.472   1st Qu.:  4.807  
    ##  Median : 69.72   Median : 1.700   Median : 2.115   Median : 68.215  
    ##  Mean   :163.92   Mean   : 4.323   Mean   : 4.450   Mean   : 62.075  
    ##  3rd Qu.:149.00   3rd Qu.: 4.037   3rd Qu.: 6.465   3rd Qu.:106.965  
    ##  Max.   :600.00   Max.   :17.960   Max.   :27.180   Max.   :224.426  
    ##                                                                      
    ##   pTimeTarget         pTimeCCW         pTimeOPP         pTimeCW       
    ##  Min.   :0.00000   Min.   :0.0129   Min.   :0.0105   Min.   :0.00000  
    ##  1st Qu.:0.01137   1st Qu.:0.1958   1st Qu.:0.2173   1st Qu.:0.09052  
    ##  Median :0.14200   Median :0.2635   Median :0.2833   Median :0.19810  
    ##  Mean   :0.14915   Mean   :0.2983   Mean   :0.3557   Mean   :0.19685  
    ##  3rd Qu.:0.26360   3rd Qu.:0.3498   3rd Qu.:0.4826   3rd Qu.:0.27005  
    ##  Max.   :0.54510   Max.   :0.9787   Max.   :0.9578   Max.   :0.90630  
    ##                                                                       
    ##  RayleigLength     RayleigAngle     PolarAvgVal       PolarSdVal    
    ##  Min.   :0.0200   Min.   :  0.71   Min.   : 75.32   Min.   : 24.05  
    ##  1st Qu.:0.1300   1st Qu.: 98.82   1st Qu.:169.01   1st Qu.: 54.48  
    ##  Median :0.3100   Median :143.76   Median :182.34   Median : 91.83  
    ##  Mean   :0.3959   Mean   :154.91   Mean   :190.08   Mean   : 81.55  
    ##  3rd Qu.:0.6800   3rd Qu.:189.31   3rd Qu.:222.35   3rd Qu.:101.66  
    ##  Max.   :0.9200   Max.   :358.16   Max.   :280.00   Max.   :156.58  
    ##                                                                     
    ##   PolarMinVal        PolarMinBin    Min50.RngLoBin  Min50.RngHiBin
    ##  Min.   :0.000000   Min.   :  0.0   Min.   :  0.0   Min.   :  0   
    ##  1st Qu.:0.000100   1st Qu.:  0.0   1st Qu.:120.0   1st Qu.: 80   
    ##  Median :0.006950   Median :160.0   Median :160.0   Median :130   
    ##  Mean   :0.007258   Mean   :150.2   Mean   :162.4   Mean   :151   
    ##  3rd Qu.:0.014150   3rd Qu.:270.0   3rd Qu.:200.0   3rd Qu.:210   
    ##  Max.   :0.021500   Max.   :350.0   Max.   :350.0   Max.   :350   
    ##                                                                   
    ##   PolarMaxVal       PolarMaxBin    Max50.RngLoBin  Max50.RngHiBin 
    ##  Min.   :0.03650   Min.   :  0.0   Min.   :  0.0   Min.   :  0.0  
    ##  1st Qu.:0.04582   1st Qu.:100.0   1st Qu.: 70.0   1st Qu.:130.0  
    ##  Median :0.05795   Median :150.0   Median :120.0   Median :180.0  
    ##  Mean   :0.07174   Mean   :156.1   Mean   :143.3   Mean   :176.6  
    ##  3rd Qu.:0.09450   3rd Qu.:200.0   3rd Qu.:210.0   3rd Qu.:220.0  
    ##  Max.   :0.17710   Max.   :350.0   Max.   :350.0   Max.   :350.0  
    ##                                                                   
    ##  AnnularMinVal      AnnularMinBin   AnnularMaxVal    AnnularMaxBin  
    ##  Min.   :0.000100   Min.   : 3.50   Min.   :0.1817   Min.   : 3.50  
    ##  1st Qu.:0.001100   1st Qu.: 3.50   1st Qu.:0.3004   1st Qu.:15.00  
    ##  Median :0.004200   Median :11.10   Median :0.3602   Median :16.60  
    ##  Mean   :0.009289   Mean   :12.11   Mean   :0.3765   Mean   :15.28  
    ##  3rd Qu.:0.010575   3rd Qu.:18.00   3rd Qu.:0.4281   3rd Qu.:16.60  
    ##  Max.   :0.094900   Max.   :19.40   Max.   :0.7857   Max.   :18.00  
    ##                                                                     
    ##    AnnularAvg      AnnularSd     AnnularSkewnes   AnnularKurtosis 
    ##  Min.   : 4.91   Min.   : 4.84   Min.   :-1.380   Min.   : 1.470  
    ##  1st Qu.:13.27   1st Qu.:11.54   1st Qu.: 1.120   1st Qu.: 4.065  
    ##  Median :15.31   Median :15.57   Median : 1.590   Median : 6.880  
    ##  Mean   :14.62   Mean   :16.25   Mean   : 1.662   Mean   : 8.466  
    ##  3rd Qu.:16.29   3rd Qu.:20.26   3rd Qu.: 2.158   3rd Qu.:10.185  
    ##  Max.   :17.59   Max.   :41.95   Max.   : 4.150   Max.   :65.450  
    ##                                                                   
    ##      Speed1            Speed2         Time1stEntrLog  
    ##  Min.   :0.00000   Min.   :0.002154   Min.   :-1.109  
    ##  1st Qu.:0.01780   1st Qu.:0.022085   1st Qu.: 2.217  
    ##  Median :0.02400   Median :0.026601   Median : 3.180  
    ##  Mean   :0.02690   Mean   :0.030936   Mean   : 3.157  
    ##  3rd Qu.:0.03106   3rd Qu.:0.034506   3rd Qu.: 4.205  
    ##  Max.   :0.11278   Max.   :0.100000   Max.   : 6.397  
    ## 

``` r
## subset to view just the data from the retention session
retention <- behavior %>% filter(TrainSessionCombo == "Retention") %>%  droplevels() 
```

``` r
behaviorsummaryTime <- summarise(group_by(behavior, APA, TrainSessionComboNum), m = mean(Time1stEntr), se = sd(Time1stEntr)/sqrt(length(Time1stEntr)))
behaviorsummaryTime
```

    ## Source: local data frame [27 x 4]
    ## Groups: APA [?]
    ## 
    ##           APA TrainSessionComboNum        m       se
    ##        <fctr>                <int>    <dbl>    <dbl>
    ## 1     control                    1 11.66765 2.172617
    ## 2     control                    2 12.99706 4.081765
    ## 3     control                    3 41.59588 7.201730
    ## 4     control                    4 28.79882 5.581521
    ## 5     control                    5 24.91706 6.749053
    ## 6     control                    6 30.78412 4.045806
    ## 7     control                    7 15.08588 2.535963
    ## 8     control                    8 34.96000 8.202296
    ## 9     control                    9 41.13412 9.478989
    ## 10 consistent                    1  6.14125 2.268263
    ## # ... with 17 more rows

``` r
behaviorsummaryNum <- summarise(group_by(behavior, APA, TrainSessionComboNum), m = mean(NumEntrances), se = sd(NumEntrances)/sqrt(length(NumEntrances)))
behaviorsummaryNum
```

    ## Source: local data frame [27 x 4]
    ## Groups: APA [?]
    ## 
    ##           APA TrainSessionComboNum        m        se
    ##        <fctr>                <int>    <dbl>     <dbl>
    ## 1     control                    1 31.88235 1.2153502
    ## 2     control                    2 21.05882 1.1129934
    ## 3     control                    3 14.41176 1.2006702
    ## 4     control                    4 15.11765 1.2715274
    ## 5     control                    5 16.47059 0.9930555
    ## 6     control                    6 14.94118 1.0017286
    ## 7     control                    7 14.17647 0.8838223
    ## 8     control                    8 15.76471 1.3787500
    ## 9     control                    9 15.41176 1.0916084
    ## 10 consistent                    1 31.12500 2.1333868
    ## # ... with 17 more rows

``` r
## see the makesessionheatmap documentataion for data tidying and plot specifications
scaledaveragedata <- as.data.frame(makescaledaveragedata(behavior))
summary(scaledaveragedata)
columnannotations <- as.data.frame(makecolumnannotations(scaledaveragedata))
summary(columnannotations)
scaledaveragedata <- scaledaveragedata[-1,]
```

``` r
longdata <- makelongdata(behavior)
Z <- longdata[,3:371]
Z <- Z[,apply(Z, 2, var, na.rm=TRUE) != 0]
pc = prcomp(Z, scale.=TRUE)
loadings <- pc$rotation
scores <- pc$x
scoresdf <- makepcadf(behavior) #create the df of pcas
rotationdf <- mkrotationdf(behavior) #loadings for specific factors
```

    ## 'data.frame':    360 obs. of  35 variables:
    ##  $ PC1     : num  -0.0441 -0.049 -0.0233 -0.0414 -0.017 ...
    ##  $ PC2     : num  0.0231 -0.0202 0.0262 -0.0455 0.0466 ...
    ##  $ PC3     : num  0.0757 0.0738 0.0351 0.0183 0.012 ...
    ##  $ PC4     : num  -0.0987 -0.1047 -0.0887 -0.071 -0.0571 ...
    ##  $ PC5     : num  -0.02711 -0.02321 -0.02452 -0.00852 -0.07279 ...
    ##  $ PC6     : num  -0.0215 -0.0179 -0.0305 -0.0114 0.0113 ...
    ##  $ PC7     : num  -0.01731 -0.06523 0.06213 -0.10985 0.00896 ...
    ##  $ PC8     : num  -0.1057 -0.0675 -0.1347 0.0086 0.0485 ...
    ##  $ PC9     : num  0.0408 0.0561 0.0339 0.0312 0.1213 ...
    ##  $ PC10    : num  0.0115 0.0532 -0.0177 0.0375 -0.0291 ...
    ##  $ PC11    : num  -0.13857 -0.0563 -0.13212 0.00426 -0.13492 ...
    ##  $ PC12    : num  -0.0358 -0.0116 -0.0637 -0.0132 -0.0184 ...
    ##  $ PC13    : num  -0.00952 -0.03178 0.0117 0.07186 -0.06115 ...
    ##  $ PC14    : num  0.02667 0.05893 -0.00494 0.10288 0.04849 ...
    ##  $ PC15    : num  -0.0486 -0.0975 0.0156 -0.1 0.017 ...
    ##  $ PC16    : num  0.0126 0.0338 0.0188 0.0495 -0.0434 ...
    ##  $ PC17    : num  0.02118 0.02026 0.07371 -0.06109 0.00136 ...
    ##  $ PC18    : num  0.0474 0.0369 0.0461 0.0291 0.0213 ...
    ##  $ PC19    : num  0.00154 -0.03977 -0.02195 -0.09976 -0.03083 ...
    ##  $ PC20    : num  -0.006 0.0534 -0.1055 0.0741 -0.0121 ...
    ##  $ PC21    : num  0.0838 0.0272 0.0806 -0.0904 0.0954 ...
    ##  $ PC22    : num  0.01047 0.09824 0.06281 -0.02868 -0.00431 ...
    ##  $ PC23    : num  0.03905 -0.01752 -0.00796 -0.07779 0.08059 ...
    ##  $ PC24    : num  -0.0424 -0.0148 -0.0056 0.0215 0.1458 ...
    ##  $ PC25    : num  0.0334 0.0125 -0.0428 -0.0477 0.0806 ...
    ##  $ PC26    : num  -0.0887 -0.0712 -0.0675 -0.1032 0.0442 ...
    ##  $ PC27    : num  -0.013009 0.008914 0.000234 0.081901 -0.047246 ...
    ##  $ PC28    : num  -0.0156 0.0372 -0.07 0.0401 -0.0193 ...
    ##  $ PC29    : num  0.04799 0.00863 0.000895 0.045915 -0.044004 ...
    ##  $ PC30    : num  -0.04214 -0.03311 0.04827 0.01359 0.00935 ...
    ##  $ PC31    : num  -0.00415 -0.0528 0.05691 -0.04775 0.10618 ...
    ##  $ PC32    : num  -0.0168 0.0344 -0.1252 0.0473 -0.1981 ...
    ##  $ PC33    : num  0.0338 0.0536 0.066 0.0211 -0.0472 ...
    ##  $ PC34    : num  -0.22674 0.00302 0.33786 -0.06677 0.08603 ...
    ##  $ variable: Factor w/ 360 levels "Hab_AnnularAvg",..: 1 2 3 4 5 6 7 8 9 10 ...

``` r
behaviormatrix <- behavior[c(20:58)]  # for 2nd pca analysis
```

``` r
write.csv(behavior, file = "../data/01a_behavior.csv", row.names = FALSE)
write.csv(retention, file = "../data/01a_retention.csv", row.names = FALSE)
write.csv(behaviorsummaryTime, file = "../data/01a_behaviorsummaryTime.csv", row.names = FALSE)
write.csv(behaviorsummaryNum, file = "../data/01a_behaviorsummaryNum.csv", row.names = FALSE)
write.csv(scaledaveragedata, file = "../data/01a_scaledaveragedata.csv", row.names = TRUE)
write.csv(columnannotations, file = "../data/01a_columnannotations.csv", row.names = TRUE)
write.csv(scoresdf, file = "../data/01a_scoresdf.csv", row.names = FALSE)
write.csv(rotationdf, file = "../data/01a_rotationdf.csv", row.names = TRUE)
write.csv(behaviormatrix, file = "../data/01a_behaviormatrix.csv", row.names = TRUE)
```
