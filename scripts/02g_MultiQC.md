    # read meta data for plotting
    colData <- read.csv("../data/02a_colData.csv")
    colData <- colData %>% select(RNAseqID, treatment)  

    # stats from fastqc 
    fastqc <- read.csv(file = "../data/multiqc/multiqc_fastqc.csv")
    fastqc$read <- ifelse(grepl("R1", fastqc$Sample.Name), "R1", "R2")  # make read columns
    fastqc$mouse <- sapply(strsplit(as.character(fastqc$Sample.Name),"\\_"), "[", 1)
    fastqc$subfield <- sapply(strsplit(as.character(fastqc$Sample.Name),"\\_"), "[", 2)
    fastqc$section <- sapply(strsplit(as.character(fastqc$Sample.Name),"\\_"), "[", 3)
    fastqc$RNAseqID <- paste(fastqc$mouse, fastqc$subfield, fastqc$section, sep = "-")
    fastqc <- left_join(colData,fastqc) %>% 
      select(RNAseqID, treatment, subfield, read, QualityFiltered, Dups, GC, Length, MillionReads)  %>% 
      filter(QualityFiltered == "No")

    ## Joining, by = "RNAseqID"

    ## Warning: Column `RNAseqID` joining factor and character vector, coercing
    ## into character vector

    fastqc$treatment <- factor(fastqc$treatment, levels = levelstreatment)
    fastqc$subfield <- factor(fastqc$subfield, levels = levelssubfield)
    head(fastqc)

    ##     RNAseqID        treatment subfield read QualityFiltered Dups   GC
    ## 1 143A-CA3-1 conflict.trained      CA3   R1              No 0.81 0.49
    ## 2 143A-CA3-1 conflict.trained      CA3   R2              No 0.71 0.48
    ## 3  143A-DG-1 conflict.trained       DG   R1              No 0.72 0.50
    ## 4  143A-DG-1 conflict.trained       DG   R2              No 0.63 0.49
    ## 5 143B-CA1-1   conflict.yoked      CA1   R1              No 0.71 0.48
    ## 6 143B-CA1-1   conflict.yoked      CA1   R2              No 0.62 0.48
    ##   Length MillionReads
    ## 1    150          5.8
    ## 2    150          5.8
    ## 3    150          7.9
    ## 4    150          7.9
    ## 5    150          2.9
    ## 6    150          2.9

    # stats from kallisto
    kallisto <- read.csv(file = "../data/multiqc/multiqc_kallisto.csv")
    kallisto$mouse <- sapply(strsplit(as.character(kallisto$sample),"\\_"), "[", 1)
    kallisto$subfield <- sapply(strsplit(as.character(kallisto$sample),"\\_"), "[", 2)
    kallisto$section <- sapply(strsplit(as.character(kallisto$sample),"\\_"), "[", 3)
    kallisto$RNAseqID <- paste(kallisto$mouse, kallisto$subfield, kallisto$section, sep = "-")
    kallisto <- left_join(colData,kallisto) %>% 
      select(RNAseqID, treatment, subfield, QC, bp, fracalign, millalign) %>% 
      filter(QC == "raw")

    ## Joining, by = "RNAseqID"

    ## Warning: Column `RNAseqID` joining factor and character vector, coercing
    ## into character vector

    kallisto$treatment <- factor(kallisto$treatment, levels = levelstreatment)
    kallisto$subfield <- factor(kallisto$subfield, levels = levelssubfield)
    head(kallisto)

    ##     RNAseqID        treatment subfield  QC    bp fracalign millalign
    ## 1 143A-CA3-1 conflict.trained      CA3 raw 201.1      0.60       3.5
    ## 2  143A-DG-1 conflict.trained       DG raw 198.9      0.69       5.4
    ## 3 143B-CA1-1   conflict.yoked      CA1 raw 200.6      0.62       1.8
    ## 4  143B-DG-1   conflict.yoked       DG raw 200.0      0.57       2.2
    ## 5 143C-CA1-1 standard.trained      CA1 raw 198.4      0.66       2.3
    ## 6 143D-CA1-3   standard.yoked      CA1 raw 197.0      0.37       1.2

    multiqc <- left_join(fastqc, kallisto) 

    ## Joining, by = c("RNAseqID", "treatment", "subfield")

    multiqc$treatment <-   fct_recode(multiqc$treatment, "standard\nyoked" = "standard.yoked", 
                                      "standard\ntrained" = "standard.trained",
                                      "conflict\nyoked" = "conflict.yoked", 
                                      "conflict\ntrained" = "conflict.trained")
    head(multiqc)

    ##     RNAseqID         treatment subfield read QualityFiltered Dups   GC
    ## 1 143A-CA3-1 conflict\ntrained      CA3   R1              No 0.81 0.49
    ## 2 143A-CA3-1 conflict\ntrained      CA3   R2              No 0.71 0.48
    ## 3  143A-DG-1 conflict\ntrained       DG   R1              No 0.72 0.50
    ## 4  143A-DG-1 conflict\ntrained       DG   R2              No 0.63 0.49
    ## 5 143B-CA1-1   conflict\nyoked      CA1   R1              No 0.71 0.48
    ## 6 143B-CA1-1   conflict\nyoked      CA1   R2              No 0.62 0.48
    ##   Length MillionReads  QC    bp fracalign millalign
    ## 1    150          5.8 raw 201.1      0.60       3.5
    ## 2    150          5.8 raw 201.1      0.60       3.5
    ## 3    150          7.9 raw 198.9      0.69       5.4
    ## 4    150          7.9 raw 198.9      0.69       5.4
    ## 5    150          2.9 raw 200.6      0.62       1.8
    ## 6    150          2.9 raw 200.6      0.62       1.8

    a <- ggplot(multiqc, aes(x = treatment, y = MillionReads, color = treatment)) +
      geom_boxplot() + geom_point() +
      scale_color_manual(values = treatmentcolors2) +  
      #facet_wrap(~subfield, nrow = 3) +
      labs(y = "Total Reads (millions)", x = NULL, subtitle = " ") +
      theme_ms() + theme(legend.position = "none")


    b <- ggplot(multiqc, aes(x = treatment, y = millalign, color = treatment)) +
      geom_boxplot() + geom_point() +
      scale_color_manual(values = treatmentcolors2) +  
      labs(y = "Pseudo-aligned Reads (millions)", x = NULL, subtitle = " ") +
      theme_ms() + theme(legend.position = "none")

    c <- ggplot(multiqc, aes(x = treatment, y = fracalign, color = treatment)) +
             geom_boxplot() + geom_point() +
      scale_color_manual(values = treatmentcolors2) +  
      labs(y = "Fraction aligned (millions)", x = NULL, subtitle = " ") +
      theme_ms() + theme(legend.position = "none")


    plot_grid( a, b , c , 
              labels = c( "(a)", "(b)", "(c)"),
              label_size = 8, nrow = 1)

![](../figures/02g_MultiQC/multiqc-1.png)

    summary(multiqc)

    ##    RNAseqID                    treatment   subfield     read          
    ##  Length:88          standard\nyoked  :18   DG :32   Length:88         
    ##  Class :character   standard\ntrained:18   CA3:26   Class :character  
    ##  Mode  :character   conflict\nyoked  :24   CA1:30   Mode  :character  
    ##                     conflict\ntrained:28                              
    ##                                                                       
    ##                                                                       
    ##  QualityFiltered      Dups              GC             Length   
    ##  No :88          Min.   :0.4700   Min.   :0.4400   Min.   :150  
    ##  Yes: 0          1st Qu.:0.6400   1st Qu.:0.4700   1st Qu.:150  
    ##                  Median :0.7300   Median :0.4800   Median :150  
    ##                  Mean   :0.7326   Mean   :0.4782   Mean   :150  
    ##                  3rd Qu.:0.8400   3rd Qu.:0.4900   3rd Qu.:150  
    ##                  Max.   :0.9600   Max.   :0.5400   Max.   :150  
    ##   MillionReads             QC           bp          fracalign    
    ##  Min.   : 1.800   filtertrim: 0   Min.   :161.4   Min.   :0.050  
    ##  1st Qu.: 3.750   raw       :88   1st Qu.:196.7   1st Qu.:0.250  
    ##  Median : 5.250                   Median :198.9   Median :0.480  
    ##  Mean   : 6.900                   Mean   :198.6   Mean   :0.428  
    ##  3rd Qu.: 7.425                   3rd Qu.:201.3   3rd Qu.:0.620  
    ##  Max.   :37.900                   Max.   :214.8   Max.   :0.690  
    ##    millalign     
    ##  Min.   : 0.100  
    ##  1st Qu.: 1.275  
    ##  Median : 2.200  
    ##  Mean   : 2.584  
    ##  3rd Qu.: 3.400  
    ##  Max.   :12.100

    mean(multiqc$MillionReads)

    ## [1] 6.9

    sd(multiqc$MillionReads)

    ## [1] 6.345059

    mean(multiqc$millalign)

    ## [1] 2.584091

    sd(multiqc$millalign)

    ## [1] 2.102647

    mean(multiqc$fracalign)

    ## [1] 0.4279545

    sd(multiqc$fracalign)

    ## [1] 0.2093842

    summary(aov(MillionReads ~ treatment, multiqc))

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment    3    179   59.55   1.505  0.219
    ## Residuals   84   3324   39.57

    summary(aov(millalign ~ treatment, multiqc))

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## treatment    3   31.8    10.6   2.524 0.0631 .
    ## Residuals   84  352.8     4.2                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    summary(aov(fracalign ~ treatment, multiqc))

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## treatment    3  0.300 0.10003   2.391 0.0744 .
    ## Residuals   84  3.514 0.04183                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    a <- ggplot(multiqc, aes(x = subfield, y = MillionReads, color = subfield)) +
      geom_boxplot() + geom_point() +
      scale_color_manual(values = colorvalsubfield) +  
      #facet_wrap(~subfield, nrow = 3) +
      labs(y = "Total Reads (millions)", x = NULL, subtitle = " ") +
      theme_ms() + theme(legend.position = "none")


    b <- ggplot(multiqc, aes(x = subfield, y = millalign, color = subfield)) +
      geom_boxplot() + geom_point() +
      scale_color_manual(values = colorvalsubfield) +  
      labs(y = "Pseudo-aligned Reads (millions)", x = NULL, subtitle = " ") +
      theme_ms() + theme(legend.position = "none")

    c <- ggplot(multiqc, aes(x = subfield, y = fracalign, color = subfield)) +
             geom_boxplot() + geom_point() +
      scale_color_manual(values = colorvalsubfield) +  
      labs(y = "Fraction aligned (millions)", x = NULL, subtitle = " ") +
      theme_ms() + theme(legend.position = "none")


    plot_grid( a, b , c , 
              labels = c( "(a)", "(b)", "(c)"),
              label_size = 8, nrow = 1)

![](../figures/02g_MultiQC/multiqc-2.png)

    summary(aov(MillionReads ~ subfield, multiqc))

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## subfield     2    228  113.78   2.953 0.0576 .
    ## Residuals   85   3275   38.53                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    summary(aov(millalign ~ subfield, multiqc))

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## subfield     2    8.1   4.034   0.911  0.406
    ## Residuals   85  376.6   4.430

    summary(aov(fracalign ~ subfield, multiqc))

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## subfield     2  0.039 0.01943   0.438  0.647
    ## Residuals   85  3.775 0.04442
