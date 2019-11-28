    ## load libraries 
    library(tidyverse) ## for respahing data

    ## set output file for figures 
    knitr::opts_chunk$set(echo = TRUE)

    # because this csv file has some rows with text, I have to specific which columns should be double
    mycols <- cols(
      "Year"  = col_integer(),  "ID"  = col_character(), "Genotype" = col_character(), "TrainProtocol"  = col_character(), 
    "TrainSequence"  = col_factor(),     "TrainGroup"  = col_factor(),        "Day"   = col_integer(),  "TrainSession"   = col_factor(),    
    "ShockOnOff"  = col_character(),       "PairedPartner" = col_character(),     "Experimenter"  = col_character(),     "Housing"  = col_character(),         
     "TestLocation" = col_character(),      "filename" = col_character(),          "p-miss" = col_double(),            "TotalTime(s)"  = col_double(),     
    "TotalPath(Arena)" = col_double(),   "Speed (Arena)" = col_double(),      "sd Speed(Arena)" = col_double(),    "Linearity(Arena)"  = col_double(), 
     "TotalPath(m)"= col_double(),        "#Entrances"  = col_double(),        "Time1stEntr" = col_double(),        "Path1stEntr"  = col_double(),      
     "Speed1stEntr(cm/s)" = col_double(), "Entr/Dist(1/m)" = col_double(),     "#Shock"  = col_double(),            "Time1stShock"  = col_double(),     
     "Path1stShock"  = col_double(),      "Speed (cm/s)"   = col_double(),     "SdSpeed"  = col_double(),           "Linearity" = col_double(),         
     "MaxTimeAvoid" = col_double(),       "MaxPathAvoid"  = col_double(),      "Time2ndEntr"   = col_double(),      "Path2ndEntr" = col_double(),       
    "Speed2ndEntr"   = col_double(),     "TimeTarget" = col_double(),         "pTimeTarget"  = col_double(),       "pTimeCCW" = col_double(),          
     "pTimeOPP"    = col_double(),        "pTimeCW" = col_double(),            "RayleigLength" = col_double(),      "RayleigAngle"   = col_double(),    
     "PolarAvgVal"  = col_double(),       "PolarSdVal"= col_double(),          "PolarMinVal"  = col_double(),       "PolarMinBin"   = col_double(),     
     "Min50%RngLoBin" = col_double(),     "Min50%RngHiBin" = col_double(),     "PolarMaxVal"  = col_double(),       "PolarMaxBin" = col_double(),       
     "Max50%RngLoBin" = col_double(),     "Max50%RngHiBin" = col_double(),     "AnnularMinVal"    = col_double(),   "AnnularMinBin"  = col_double(),    
     "AnnularMaxVal"  = col_double(),     "AnnularMaxBin" = col_double(),      "AnnularAvg"  = col_double(),        "AnnularSd"  = col_double(),        
     "AnnularSkewnes"  = col_double(),    "AnnularKurtosis"= col_double() 
    )

    # import raw data
    behav <- read_csv("../data/00_Data2013_2016_forAnalysis.csv", col_types = mycols) 

    ## Warning: 650 parsing failures.
    ## row              col   expected           actual                                       file
    ##  97 Year             an integer Year             '../data/00_Data2013_2016_forAnalysis.csv'
    ##  97 Day              an integer Day              '../data/00_Data2013_2016_forAnalysis.csv'
    ##  97 p-miss           a double   p-miss           '../data/00_Data2013_2016_forAnalysis.csv'
    ##  97 TotalTime(s)     a double   TotalTime(s)     '../data/00_Data2013_2016_forAnalysis.csv'
    ##  97 TotalPath(Arena) a double   TotalPath(Arena) '../data/00_Data2013_2016_forAnalysis.csv'
    ## ... ................ .......... ................ ..........................................
    ## See problems(...) for more details.

    behav <- as.data.frame(behav)
    names(behav)

    ##  [1] "Year"               "ID"                 "Genotype"          
    ##  [4] "TrainProtocol"      "TrainSequence"      "TrainGroup"        
    ##  [7] "Day"                "TrainSession"       "ShockOnOff"        
    ## [10] "PairedPartner"      "Experimenter"       "Housing"           
    ## [13] "TestLocation"       "filename"           "p-miss"            
    ## [16] "TotalTime(s)"       "TotalPath(Arena)"   "Speed (Arena)"     
    ## [19] "sd Speed(Arena)"    "Linearity(Arena)"   "TotalPath(m)"      
    ## [22] "#Entrances"         "Time1stEntr"        "Path1stEntr"       
    ## [25] "Speed1stEntr(cm/s)" "Entr/Dist(1/m)"     "#Shock"            
    ## [28] "Time1stShock"       "Path1stShock"       "Speed (cm/s)"      
    ## [31] "SdSpeed"            "Linearity"          "MaxTimeAvoid"      
    ## [34] "MaxPathAvoid"       "Time2ndEntr"        "Path2ndEntr"       
    ## [37] "Speed2ndEntr"       "TimeTarget"         "pTimeTarget"       
    ## [40] "pTimeCCW"           "pTimeOPP"           "pTimeCW"           
    ## [43] "RayleigLength"      "RayleigAngle"       "PolarAvgVal"       
    ## [46] "PolarSdVal"         "PolarMinVal"        "PolarMinBin"       
    ## [49] "Min50%RngLoBin"     "Min50%RngHiBin"     "PolarMaxVal"       
    ## [52] "PolarMaxBin"        "Max50%RngLoBin"     "Max50%RngHiBin"    
    ## [55] "AnnularMinVal"      "AnnularMinBin"      "AnnularMaxVal"     
    ## [58] "AnnularMaxBin"      "AnnularAvg"         "AnnularSd"         
    ## [61] "AnnularSkewnes"     "AnnularKurtosis"

    # remove some columns that are not behavioral measures
    columnstoremove <- c("PolarAvgVal" , "PolarSdVal", "PolarMinVal" ,
                         "PolarMinBin", "Max50%RngLoBin", "Min50%RngHiBin" , 
                         "PolarMaxVal", "PolarMaxBin", "Max50%RngLoBin", "Max50%RngHiBin",
                         "p-miss", "TotalTime(s)", "TotalPath(m)",  "Time1stShock",  
                         "Path1stShock", "Speed (cm/s)", "SdSpeed", "Linearity",
                         "AnnularMinVal" , "AnnularMinBin" , "AnnularMaxVal" , 
                         "AnnularMaxBin", "AnnularAvg"  , "AnnularSd" )

    # rename some columns
    names(behav)[names(behav)=="#Shock"] <- "NumShock"
    names(behav)[names(behav)=="#Entrances"] <- "NumEntrances"
    names(behav)[names(behav)=="pTimeTarget"] <- "pTimeShockZone"
    names(behav)[names(behav)=="TimeTarget"] <- "TimeShockZone"
    names(behav)[names(behav)=="Speed (Arena)"] <- "SpeedArena.cm.s"


    # keep only 2015 animals and remove unnecessary columns
    behav <- behav %>% dplyr::filter(Year == "2015") %>% 
                       dplyr::select(-columnstoremove) %>% 
                       
                      droplevels()

    # check which measures are included
    measures <- behav %>% select(`TotalPath(Arena)`:AnnularKurtosis)
    names(measures)

    ##  [1] "TotalPath(Arena)"   "SpeedArena.cm.s"    "sd Speed(Arena)"   
    ##  [4] "Linearity(Arena)"   "NumEntrances"       "Time1stEntr"       
    ##  [7] "Path1stEntr"        "Speed1stEntr(cm/s)" "Entr/Dist(1/m)"    
    ## [10] "NumShock"           "MaxTimeAvoid"       "MaxPathAvoid"      
    ## [13] "Time2ndEntr"        "Path2ndEntr"        "Speed2ndEntr"      
    ## [16] "TimeShockZone"      "pTimeShockZone"     "pTimeCCW"          
    ## [19] "pTimeOPP"           "pTimeCW"            "RayleigLength"     
    ## [22] "RayleigAngle"       "Min50%RngLoBin"     "AnnularSkewnes"    
    ## [25] "AnnularKurtosis"

    # make new factors
    behav <- behav %>%
      mutate(APA = paste(TrainSequence, TrainGroup, sep = "-")) %>%
      mutate(treatment = fct_recode(APA,
                                    "standard.yoked" = "train-train-yoked",
                                    "standard.trained" = "train-train-trained",
                                    "conflict.yoked" = "train-conflict-yoked",
                                    "conflict.trained" = "train-conflict-trained")) %>%
      mutate(training = fct_collapse(treatment,
                                          trained = c("standard.trained", "conflict.trained"),
                                          yoked = c("standard.yoked", "conflict.yoked")))  %>%
      mutate(trial = fct_collapse(TrainSession,
                                          "T4_C1" = c("T4", "C1"),
                                          "T5_C2" = c("T5", "C2"),
                                          "T6_C3" = c("T6", "C3"),
                                          "Hab"  = "Hab",
                                          "T1" = "T1",
                                          "T2"= "T2",
                                          "T3"= "T3",
                                          "Retest" = "Retest",
                                          "Retention" = "Retention"))

    behav$trial <- factor(behav$trial, 
                                      levels = c("Hab", "T1","T2","T3","Retest",
                                       "T4_C1","T5_C2","T6_C3","Retention"))
    behav$trialNum <- as.numeric(behav$trial)

    behav <- behav %>%  dplyr::select(ID,Day,treatment, training,trial,
                                     trialNum, ShockOnOff, PairedPartner,
                                     `TotalPath(Arena)`:AnnularKurtosis) %>%
                        dplyr::arrange(ID)

    head(behav)

    ##       ID Day        treatment training  trial trialNum ShockOnOff
    ## 1 15140A   1 conflict.trained  trained    Hab        1        Off
    ## 2 15140A   1 conflict.trained  trained     T1        2         On
    ## 3 15140A   1 conflict.trained  trained     T2        3         On
    ## 4 15140A   1 conflict.trained  trained     T3        4         On
    ## 5 15140A   2 conflict.trained  trained Retest        5         On
    ## 6 15140A   2 conflict.trained  trained  T4_C1        6         On
    ##   PairedPartner TotalPath(Arena) SpeedArena.cm.s sd Speed(Arena)
    ## 1        15140B            22.68            3.78            3.07
    ## 2        15140B            19.36            3.23            2.78
    ## 3        15140B            15.01            2.50            2.68
    ## 4        15140B            14.39            2.40            2.78
    ## 5        15140B            14.04            2.34            3.11
    ## 6        15140B            12.50            2.08            2.52
    ##   Linearity(Arena) NumEntrances Time1stEntr Path1stEntr Speed1stEntr(cm/s)
    ## 1           0.4790           28       24.63        1.09               4.56
    ## 2           0.4016            6        9.83        0.62              16.42
    ## 3           0.3170            2      118.37        3.17               2.31
    ## 4           0.3122            3      256.53        7.48               4.26
    ## 5           0.2895            1      432.07       10.56               9.38
    ## 6           0.3107           10        0.87        0.00              -1.00
    ##   Entr/Dist(1/m) NumShock MaxTimeAvoid MaxPathAvoid Time2ndEntr
    ## 1           1.12       52           53         2.15       59.97
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
    ##   pTimeCW RayleigLength RayleigAngle Min50%RngLoBin AnnularSkewnes
    ## 1  0.3352          0.11       330.67             60           0.88
    ## 2  0.0779          0.65       112.66            130           1.81
    ## 3  0.0250          0.78       124.87            150           1.87
    ## 4  0.0123          0.80       128.39            150           2.84
    ## 5  0.0729          0.72       159.36            170           2.42
    ## 6  0.7905          0.67       257.90            280           0.98
    ##   AnnularKurtosis
    ## 1            3.13
    ## 2            6.70
    ## 3            8.91
    ## 4           12.51
    ## 5           11.83
    ## 6            4.65

    ## Fix the NumShocks problem but it doesn't know the the yoked got shocked when the trained got shocked rather than they enterer the shock zone.
    # Steps:
    # 1. Create a "shocks1" df that have only ID, PairedPartner, training, trialNum, NumShock, ShockOnOff
    # 2. Separate the tarined and yoked mice into two dataframes, relabel columsn so the two dfs can be rejoined
    # 3. Creat a "shocks2" df by join the two dataframes and replace all NumShock value for all the yoked animals, and set it to 0 when shock is off
    # 4. Sparate the trained and yoked again and rename so they can be join into a final "shocks3" datafram by 
    # 5. Join the behavior data with the new shocks data. Replace old NumShock with new NumShock

    shocks1 <- behav %>% select(ID, PairedPartner, training, trialNum, NumShock, ShockOnOff) 

    trained <- shocks1 %>% filter(training == "trained") %>% 
      rename(trained = ID, yoked = PairedPartner, ttraining = training,  tNumShock = NumShock)
    yoked <- shocks1 %>% filter(training == "yoked") %>% 
      rename(yoked = ID, trained = PairedPartner, ytraining = training, yNumShock = NumShock)

    shocks2 <- full_join(trained,yoked) %>% select(trained, yoked, ShockOnOff, trialNum, tNumShock) 

    ## Joining, by = c("trained", "yoked", "trialNum", "ShockOnOff")

    shocks2$tNumShock <- ifelse(shocks2$ShockOnOff == "Off", 0, shocks2$tNumShock )

    trained <- shocks2 %>% select(trained, trialNum, tNumShock) %>% rename(ID = trained )
    yoked <- shocks2 %>% select(yoked, trialNum, tNumShock) %>% rename(ID = yoked )
    shocks3 <- rbind(trained, yoked)

    behav <- full_join(behav, shocks3) %>% mutate(NumShock = tNumShock) %>% select(-tNumShock)

    ## Joining, by = c("ID", "trialNum")

    # add new measure 
    behav <- behav %>%  mutate(ShockPerEntrance = NumShock / NumEntrances) %>%  
      mutate(ShockPerEntrance = ifelse(is.na(ShockPerEntrance), 0, ShockPerEntrance))

    write.csv(behav, "../data/00_behaviordata.csv", row.names = F)
