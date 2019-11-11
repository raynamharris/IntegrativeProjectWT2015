    library(tidyverse) 

    ## ── Attaching packages ───────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.2.1     ✔ purrr   0.3.2
    ## ✔ tibble  2.1.3     ✔ dplyr   0.8.1
    ## ✔ tidyr   0.8.3     ✔ stringr 1.4.0
    ## ✔ readr   1.3.1     ✔ forcats 0.4.0

    ## ── Conflicts ──────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

    library(corrplot)

    ## corrplot 0.84 loaded

    library(cowplot)

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    library(corrr)
    #devtools::install_github("clauswilke/ggtext")
    library(ggtext)

    source("figureoptions.R")

    theme_ms <- function () { 
      theme_classic(base_size = 8) +
        theme(
          panel.grid.major  = element_blank(),  # remove major gridlines
          panel.grid.minor  = element_blank(),  # remove minor gridlines
          plot.title = element_text(hjust = 0.5, face = "bold") # center & bold 
        )
    }

    theme2 <- function () { 
      theme_minimal(base_size = 7) +
        theme(
          panel.grid.major  = element_blank(),  # remove major gridlines
          panel.grid.minor  = element_blank() # remove minor gridlines
        )
    }

    knitr::opts_chunk$set(fig.path = '../figures/02e_correlations/', cache = F)

For this analysis, I want to explor correlations between a behavioral
measure and gene expression.

    # import behavior data, create mouse id, select relvant samples
    behav <- read.csv("../data/01a_behavior.csv") 
    behav$mouse <- sapply(strsplit(as.character(behav$ID),"15"), "[", 2)
    behav <- behav %>% filter( #treatment %in% c("conflict.trained", "standard.trained"),
                                      TrainSessionCombo == "Retention") %>% 
                               select(mouse,Entr.Dist.1.m.,MaxTimeAvoid) 
    head(behav)

    ##   mouse Entr.Dist.1.m. MaxTimeAvoid
    ## 1  140A           0.22          337
    ## 2  140B           0.61          111
    ## 3  140C           0.00          599
    ## 4  140D           0.82          115
    ## 5  141C           0.51          172
    ## 6  141D           1.29           66

    pcadata <- read_csv("../data/01a_pcadf.csv") %>%
      filter(#treatment %in% c("conflict.trained", "standard.trained"),
             TrainSessionComboNum == 9) %>%
      select(ID,PC1,PC2) 

    ## Parsed with column specification:
    ## cols(
    ##   ID = col_character(),
    ##   treatment = col_character(),
    ##   TrainSessionComboNum = col_double(),
    ##   PC1 = col_double(),
    ##   PC2 = col_double(),
    ##   PC3 = col_double(),
    ##   PC4 = col_double(),
    ##   PC5 = col_double(),
    ##   PC6 = col_double(),
    ##   PC7 = col_double(),
    ##   PC8 = col_double(),
    ##   PC9 = col_double(),
    ##   PC10 = col_double()
    ## )

    pcadata$mouse <- sapply(strsplit(as.character(pcadata$ID),"15"), "[", 2)
    pcadata$ID <- NULL
    head(pcadata)

    ## # A tibble: 6 x 3
    ##       PC1    PC2 mouse
    ##     <dbl>  <dbl> <chr>
    ## 1  2.24   -1.22  140A 
    ## 2 -0.791  -1.60  140B 
    ## 3  6.43    0.226 140C 
    ## 4 -1.21    1.71  140D 
    ## 5 -0.0494  1.12  141C 
    ## 6 -2.93    2.21  141D

    DG_DEGs <- read.csv("../data/02f_DG_DEGs_vsd.csv", row.names = 1, check.names = F)
    DG_DEGs <- as.data.frame(t(DG_DEGs))
    DG_DEGs$sample <- row.names(DG_DEGs)
    DG_DEGs$mouse <- sapply(strsplit(as.character(DG_DEGs$sample),"\\-"), "[", 1)
    DG_DEGs <- DG_DEGs %>% select(mouse,`1190002N15RIK`:ZFP869)
    head(DG_DEGs)[1:5]

    ##           mouse 1190002N15RIK A830010M20RIK      ABHD2       ACAN
    ## 143A-DG-1  143A    0.37666940     0.3321828  0.4370642  0.6682061
    ## 143B-DG-1  143B   -0.79007413    -0.4803026 -0.1164208 -0.4768060
    ## 143D-DG-3  143D   -0.01317465    -0.2575089 -0.1110985 -0.4930643
    ## 144A-DG-2  144A    1.40407413     1.4069082  1.1339066  0.8794021
    ## 144C-DG-2  144C    1.16933394     1.2067769  0.8232341  0.4271263
    ## 144D-DG-2  144D   -0.47594844    -0.2925787  0.1279683 -0.5497050

    DEGsPCA <- left_join(DG_DEGs, pcadata)

    ## Joining, by = "mouse"

    DEGsPCAbeahv <- left_join(DEGsPCA, behav)

    ## Joining, by = "mouse"

    DEGsPCAbeahv <- as.data.frame(DEGsPCAbeahv)
    row.names(DEGsPCAbeahv) <- DEGsPCAbeahv$mouse
    DEGsPCAbeahv$mouse <- NULL
    DEGsPCAbeahv <- as.matrix(DEGsPCAbeahv)

    DGcor <- DEGsPCAbeahv %>% correlate() %>% rearrange() %>%  shave()

    ## 
    ## Correlation method: 'pearson'
    ## Missing treated using: 'pairwise.complete.obs'

    ## Registered S3 method overwritten by 'seriation':
    ##   method         from 
    ##   reorder.hclust gclus

    DGcor

    ## # A tibble: 218 x 219
    ##    rowname  SNX18  PCDH8  NPAS4  KCNJ2 ERRFI1   SGK1  NFIL3 TIPARP   PLK2
    ##    <chr>    <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
    ##  1 SNX18   NA     NA     NA     NA     NA     NA     NA     NA     NA    
    ##  2 PCDH8    0.946 NA     NA     NA     NA     NA     NA     NA     NA    
    ##  3 NPAS4    0.943  0.920 NA     NA     NA     NA     NA     NA     NA    
    ##  4 KCNJ2    0.944  0.927  0.952 NA     NA     NA     NA     NA     NA    
    ##  5 ERRFI1   0.955  0.906  0.974  0.964 NA     NA     NA     NA     NA    
    ##  6 SGK1     0.961  0.920  0.953  0.925  0.940 NA     NA     NA     NA    
    ##  7 NFIL3    0.932  0.902  0.983  0.962  0.961  0.919 NA     NA     NA    
    ##  8 TIPARP   0.935  0.898  0.967  0.924  0.956  0.936  0.959 NA     NA    
    ##  9 PLK2     0.954  0.902  0.973  0.932  0.970  0.977  0.943  0.951 NA    
    ## 10 GADD45G  0.936  0.932  0.920  0.955  0.947  0.870  0.933  0.904  0.885
    ## # … with 208 more rows, and 209 more variables: GADD45G <dbl>,
    ## #   SLC16A1 <dbl>, FAM107B <dbl>, FBXO33 <dbl>, KCNF1 <dbl>, KITL <dbl>,
    ## #   GPR19 <dbl>, MEST <dbl>, PTGS2 <dbl>, JUN <dbl>, APAF1 <dbl>,
    ## #   FZD5 <dbl>, CWC25 <dbl>, PAK6 <dbl>, SMAD7 <dbl>, ZDBF2 <dbl>,
    ## #   CXADR <dbl>, TSC22D2 <dbl>, EGR4 <dbl>, A830010M20RIK <dbl>,
    ## #   ARID5B <dbl>, PLK3 <dbl>, PIGA <dbl>, KCNA4 <dbl>, IRF2BP2 <dbl>,
    ## #   HMGCR <dbl>, RFX2 <dbl>, SIAH2 <dbl>, BTG2 <dbl>, AHR <dbl>,
    ## #   ARL13B <dbl>, ARL4A <dbl>, FZD4 <dbl>, MARCH11 <dbl>, C2CD4B <dbl>,
    ## #   FOSB <dbl>, ING2 <dbl>, RASD1 <dbl>, JUNB <dbl>, LONRF1 <dbl>,
    ## #   NR4A3 <dbl>, KLF2 <dbl>, LEMD3 <dbl>, PPP1R15A <dbl>, FOS <dbl>,
    ## #   DUSP16 <dbl>, CTNND1 <dbl>, MYC <dbl>, EGR3 <dbl>, PEG10 <dbl>,
    ## #   SLC2A3 <dbl>, SH2D3C <dbl>, FRMD6 <dbl>, SCG2 <dbl>, FOXG1 <dbl>,
    ## #   FOSL2 <dbl>, GM13889 <dbl>, ERF <dbl>, LCMT2 <dbl>, MED7 <dbl>,
    ## #   CNNM1 <dbl>, ABHD2 <dbl>, TRIB1 <dbl>, STMN4 <dbl>, GMEB2 <dbl>,
    ## #   LMNA <dbl>, ANKRD28 <dbl>, DYRK2 <dbl>, ODC1 <dbl>, RASL11A <dbl>,
    ## #   PER1 <dbl>, KLF6 <dbl>, RGS2 <dbl>, B3GNT2 <dbl>, ACAN <dbl>,
    ## #   LRRTM2 <dbl>, KDM6B <dbl>, OTUD1 <dbl>, HECA <dbl>, ADAMTS1 <dbl>,
    ## #   ARC <dbl>, EIF5 <dbl>, FERMT2 <dbl>, IRS2 <dbl>, LBH <dbl>,
    ## #   IL16 <dbl>, RASL11B <dbl>, DUSP14 <dbl>, MN1 <dbl>, SLC25A25 <dbl>,
    ## #   NAF1 <dbl>, ZFP654 <dbl>, ZFP275 <dbl>, HS6ST1 <dbl>, DNAJB1 <dbl>,
    ## #   NR4A2 <dbl>, POU3F3 <dbl>, NPTX2 <dbl>, `1190002N15RIK` <dbl>,
    ## #   ZFP869 <dbl>, …

    DGcorSlim <- correlate(DEGsPCAbeahv) %>%  
      focus(PC1)   %>%  
      arrange(PC1)

    ## 
    ## Correlation method: 'pearson'
    ## Missing treated using: 'pairwise.complete.obs'

    DGcorSlim

    ## # A tibble: 217 x 2
    ##    rowname           PC1
    ##    <chr>           <dbl>
    ##  1 Entr.Dist.1.m. -0.913
    ##  2 COQ2           -0.608
    ##  3 PC2            -0.589
    ##  4 ZFP207         -0.568
    ##  5 TUBB4A         -0.559
    ##  6 NXF1           -0.547
    ##  7 GPI1           -0.545
    ##  8 PXN            -0.468
    ##  9 LRRC45         -0.443
    ## 10 EEF1E1         -0.439
    ## # … with 207 more rows

    DGcorarranged <- fashion(DGcorSlim) %>% arrange(desc(PC1))
    head(DGcorarranged)

    ##          rowname  PC1
    ## 1 Entr.Dist.1.m. -.91
    ## 2           COQ2 -.61
    ## 3            PC2 -.59
    ## 4         ZFP207 -.57
    ## 5         TUBB4A -.56
    ## 6           NXF1 -.55

    mousetreatment <- read_csv("../data/01a_behavior.csv")  %>%
      filter(TrainSessionCombo == "Retention") %>%
      mutate(mouse = sapply(strsplit(as.character(ID),"15"), "[", 2)) %>%
      select(mouse, treatment)

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   ID = col_character(),
    ##   treatment = col_character(),
    ##   training = col_character(),
    ##   TrainSessionCombo = col_character(),
    ##   ShockOnOff = col_character(),
    ##   PairedPartner = col_character()
    ## )

    ## See spec(...) for full column specifications.

    DEGsPCAbeahvDf <- as.data.frame(DEGsPCAbeahv)
    DEGsPCAbeahvDf$mouse <- row.names(DEGsPCAbeahvDf)
    DEGsPCAbeahvDf <- left_join(mousetreatment, DEGsPCAbeahvDf)

    ## Joining, by = "mouse"

    summary(lm(ARC ~ PC1, data = DEGsPCAbeahvDf))

    ## 
    ## Call:
    ## lm(formula = ARC ~ PC1, data = DEGsPCAbeahvDf)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.2473 -0.4199  0.1534  0.3754  1.4703 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -0.23810    0.18931  -1.258 0.229064    
    ## PC1          0.22541    0.04435   5.083 0.000167 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7337 on 14 degrees of freedom
    ##   (18 observations deleted due to missingness)
    ## Multiple R-squared:  0.6486, Adjusted R-squared:  0.6235 
    ## F-statistic: 25.84 on 1 and 14 DF,  p-value: 0.0001669

    summary(lm(FOSL2 ~ PC1, data = DEGsPCAbeahvDf))

    ## 
    ## Call:
    ## lm(formula = FOSL2 ~ PC1, data = DEGsPCAbeahvDf)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.7580 -0.3245  0.2459  0.5141  0.8199 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -0.21603    0.20112  -1.074 0.300926    
    ## PC1          0.20452    0.04711   4.341 0.000678 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7795 on 14 degrees of freedom
    ##   (18 observations deleted due to missingness)
    ## Multiple R-squared:  0.5737, Adjusted R-squared:  0.5433 
    ## F-statistic: 18.84 on 1 and 14 DF,  p-value: 0.0006778

    a <- ggplot(DEGsPCAbeahvDf, aes(x = PC1, y = ARC)) +
      geom_point(aes(color = treatment)) + geom_smooth(method='lm', colour = "darkgrey")  +
      labs(y = "ARC expression in DG \n (variance stabilized)", 
           x = "behavior PC1 \n increasing avoidance behavior -->",
           subtitle = "R2 = 0.6235, p = 0.00016693") +
      theme_ms() +
      scale_color_manual(values = treatmentcolors) + theme(legend.position = "none")
      
    b <- ggplot(DEGsPCAbeahvDf, aes(x = PC1, y = FOSL2)) +
      geom_point(aes(color = treatment)) + geom_smooth(method='lm', colour = "darkgrey") +
      labs(y = "FOSL2 expression in DG \n (variance stabilized)", 
           x = "behavior PC1 \n increasing avoidance behavior -->",
           subtitle = "R2 = 0.5433, p = 0.0006778") +
        theme_ms() +
      scale_color_manual(values = treatmentcolors) + theme(legend.position = "none")

    plot_grid(a,b)

    ## Warning: Removed 18 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 18 rows containing missing values (geom_point).

    ## Warning: Removed 18 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 18 rows containing missing values (geom_point).

![](../figures/02e_correlations/DEGvPCA-1.png)

    allDGvsd <- read.csv("../data/02c_DGvsd.csv", check.names = F, stringsAsFactors = F, row.names = 1)
    allDGvsd <- as.data.frame(t(allDGvsd))

    allDGvsd$ID <- row.names(allDGvsd)
    allDGvsd$mouse <- sapply(strsplit(as.character(allDGvsd$ID),"-"), "[", 1)
    allDGvsd <- left_join(mousetreatment, allDGvsd)

    ## Joining, by = "mouse"

    allDGvsd <- left_join(pcadata, allDGvsd)

    ## Joining, by = "mouse"

    summary(lm(Prkcz ~ PC1, data = allDGvsd))

    ## 
    ## Call:
    ## lm(formula = Prkcz ~ PC1, data = allDGvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.18054 -0.05397 -0.02670  0.01022  0.28031 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  8.423056   0.032702 257.568   <2e-16 ***
    ## PC1         -0.010868   0.007661  -1.419    0.178    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1267 on 14 degrees of freedom
    ##   (18 observations deleted due to missingness)
    ## Multiple R-squared:  0.1257, Adjusted R-squared:  0.06323 
    ## F-statistic: 2.012 on 1 and 14 DF,  p-value: 0.1779

    summary(lm(Camk2a ~ PC1, data = allDGvsd))

    ## 
    ## Call:
    ## lm(formula = Camk2a ~ PC1, data = allDGvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.59854 -0.20560  0.03381  0.26081  0.46591 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 12.33416    0.08606 143.321   <2e-16 ***
    ## PC1          0.02282    0.02016   1.132    0.277    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3335 on 14 degrees of freedom
    ##   (18 observations deleted due to missingness)
    ## Multiple R-squared:  0.08385,    Adjusted R-squared:  0.01841 
    ## F-statistic: 1.281 on 1 and 14 DF,  p-value: 0.2767

    c <- ggplot(allDGvsd, aes(x = PC1, y = Prkcz)) +
      geom_point(aes(color = treatment)) + geom_smooth(method='lm', colour = "darkgrey") +
      labs(y = "PRKCZ expression in DG \n (variance stabilized)", 
           x = "behavior PC1 \n increasing avoidance behavior -->",
           subtitle = "R2 = 0.06323 , p = 0.1779") +
        theme_ms() +
      scale_color_manual(values = treatmentcolors) + theme(legend.position = "none")

    summary(lm(Camk2a ~ PC1, data = allDGvsd))

    ## 
    ## Call:
    ## lm(formula = Camk2a ~ PC1, data = allDGvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.59854 -0.20560  0.03381  0.26081  0.46591 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 12.33416    0.08606 143.321   <2e-16 ***
    ## PC1          0.02282    0.02016   1.132    0.277    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3335 on 14 degrees of freedom
    ##   (18 observations deleted due to missingness)
    ## Multiple R-squared:  0.08385,    Adjusted R-squared:  0.01841 
    ## F-statistic: 1.281 on 1 and 14 DF,  p-value: 0.2767

    d <- ggplot(allDGvsd, aes(x = PC1, y = Camk2a)) +
      geom_point(aes(color = treatment)) + geom_smooth(method='lm', colour = "darkgrey") +
      labs(y = "CAMK2A expression in DG \n (variance stabilized)", 
           x = "behavior PC1 \n increasing avoidance behavior -->",
           subtitle = "R2 = 0.01841, p = 0.2767") +
        theme_ms() +
      scale_color_manual(values = treatmentcolors) + theme(legend.position = "none")

    plot_grid(a + theme(axis.title.x = element_blank()),
              b + theme(axis.title.x = element_blank()),
              c,d, rel_heights = c(0.45, 0.55))

    ## Warning: Removed 18 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 18 rows containing missing values (geom_point).

    ## Warning: Removed 18 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 18 rows containing missing values (geom_point).

    ## Warning: Removed 18 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 18 rows containing missing values (geom_point).

    ## Warning: Removed 18 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 18 rows containing missing values (geom_point).

![](../figures/02e_correlations/DEGnCandidatesvPCA-1.png)

    plotPCAvcandidates <- function(filename){
      allvsd <- read.csv(filename, check.names = F, stringsAsFactors = F, row.names = 1)
      allvsd <- as.data.frame(t(allvsd))

      allvsd$ID <- row.names(allvsd)
      allvsd$mouse <- sapply(strsplit(as.character(allvsd$ID),"-"), "[", 1)
      allvsd <- left_join(mousetreatment, allvsd)
      allvsd <- left_join(pcadata, allvsd)
      
      print(summary(lm(Arc ~ PC1, data = allvsd)))
      print(summary(lm(Camk2a ~ PC1, data = allvsd)))
      print(summary(lm(Gria1 ~ PC1, data = allvsd)))
      print(summary(lm(Grin1 ~ PC1, data = allvsd)))
      print(summary(lm(Prkci ~ PC1, data = allvsd)))
      print(summary(lm(Prkcz ~ PC1, data = allvsd)))

      candidatevsdPCA <- allvsd %>% select(mouse, treatment, PC1, PC2, mouse,
                                           Arc, Camk2a, Gria1, Grin1, Prkci, Prkcz)
      candidatevsdPCA <- gather(candidatevsdPCA, gene, expression, Arc:Prkcz, factor_key=TRUE) %>% drop_na()

      p <- ggplot(candidatevsdPCA, aes(x = PC1, y = expression)) +
        geom_point(aes(color = treatment)) + geom_smooth(method='lm', colour = "darkgrey") +
        labs(y = "gene expression", 
            x = "PC1 (estimate of increasing avoidance behavior)") +
          theme_ms() +
        scale_color_manual(values = treatmentcolors) + theme(legend.position = "none") +
        facet_wrap(~gene, scales = "free_y", nrow = 1) +
        theme(axis.text.x = element_blank()) +
        scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
      return(p)
    }

    p1 <- plotPCAvcandidates("../data/02c_DGvsd.csv")

    ## Joining, by = "mouse"
    ## Joining, by = "mouse"

    ## 
    ## Call:
    ## lm(formula = Arc ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.2360 -0.4151  0.1499  0.3651  1.4514 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  9.24279    0.18632  49.607  < 2e-16 ***
    ## PC1          0.22243    0.04365   5.096 0.000163 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7221 on 14 degrees of freedom
    ##   (18 observations deleted due to missingness)
    ## Multiple R-squared:  0.6497, Adjusted R-squared:  0.6247 
    ## F-statistic: 25.97 on 1 and 14 DF,  p-value: 0.0001629
    ## 
    ## 
    ## Call:
    ## lm(formula = Camk2a ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.59854 -0.20560  0.03381  0.26081  0.46591 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 12.33416    0.08606 143.321   <2e-16 ***
    ## PC1          0.02282    0.02016   1.132    0.277    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3335 on 14 degrees of freedom
    ##   (18 observations deleted due to missingness)
    ## Multiple R-squared:  0.08385,    Adjusted R-squared:  0.01841 
    ## F-statistic: 1.281 on 1 and 14 DF,  p-value: 0.2767
    ## 
    ## 
    ## Call:
    ## lm(formula = Gria1 ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.85004 -0.03374  0.04554  0.14516  0.33743 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 10.87872    0.07531  144.46   <2e-16 ***
    ## PC1          0.02240    0.01764    1.27    0.225    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2919 on 14 degrees of freedom
    ##   (18 observations deleted due to missingness)
    ## Multiple R-squared:  0.1033, Adjusted R-squared:  0.03922 
    ## F-statistic: 1.612 on 1 and 14 DF,  p-value: 0.2249
    ## 
    ## 
    ## Call:
    ## lm(formula = Grin1 ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.31391 -0.08040  0.05118  0.08605  0.18433 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 10.145897   0.034253  296.20   <2e-16 ***
    ## PC1          0.006020   0.008024    0.75    0.466    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1328 on 14 degrees of freedom
    ##   (18 observations deleted due to missingness)
    ## Multiple R-squared:  0.03865,    Adjusted R-squared:  -0.03002 
    ## F-statistic: 0.5628 on 1 and 14 DF,  p-value: 0.4656
    ## 
    ## 
    ## Call:
    ## lm(formula = Prkci ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -2.02573 -0.06494  0.08851  0.18068  1.26874 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  7.04318    0.17069  41.262 5.05e-16 ***
    ## PC1          0.00982    0.03999   0.246     0.81    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.6615 on 14 degrees of freedom
    ##   (18 observations deleted due to missingness)
    ## Multiple R-squared:  0.004289,   Adjusted R-squared:  -0.06683 
    ## F-statistic: 0.06031 on 1 and 14 DF,  p-value: 0.8096
    ## 
    ## 
    ## Call:
    ## lm(formula = Prkcz ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.18054 -0.05397 -0.02670  0.01022  0.28031 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  8.423056   0.032702 257.568   <2e-16 ***
    ## PC1         -0.010868   0.007661  -1.419    0.178    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1267 on 14 degrees of freedom
    ##   (18 observations deleted due to missingness)
    ## Multiple R-squared:  0.1257, Adjusted R-squared:  0.06323 
    ## F-statistic: 2.012 on 1 and 14 DF,  p-value: 0.1779

    p1

![](../figures/02e_correlations/CandidatesvPCA-1.png)

    p2 <- plotPCAvcandidates("../data/02c_CA3vsd.csv")

    ## Joining, by = "mouse"
    ## Joining, by = "mouse"

    ## 
    ## Call:
    ## lm(formula = Arc ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.62570 -0.21396  0.01013  0.25956  0.42689 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  9.15432    0.09113 100.451   <2e-16 ***
    ## PC1         -0.02140    0.02053  -1.043     0.32    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3119 on 11 degrees of freedom
    ##   (21 observations deleted due to missingness)
    ## Multiple R-squared:  0.08993,    Adjusted R-squared:  0.007192 
    ## F-statistic: 1.087 on 1 and 11 DF,  p-value: 0.3195
    ## 
    ## 
    ## Call:
    ## lm(formula = Camk2a ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.4935 -0.2405 -0.0896  0.2586  0.6890 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 12.58608    0.10033 125.446   <2e-16 ***
    ## PC1         -0.03908    0.02260  -1.729    0.112    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3434 on 11 degrees of freedom
    ##   (21 observations deleted due to missingness)
    ## Multiple R-squared:  0.2138, Adjusted R-squared:  0.1423 
    ## F-statistic: 2.991 on 1 and 11 DF,  p-value: 0.1117
    ## 
    ## 
    ## Call:
    ## lm(formula = Gria1 ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.39159 -0.18338  0.09158  0.17690  0.31951 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 10.82061    0.07170 150.906   <2e-16 ***
    ## PC1         -0.03533    0.01615  -2.188   0.0512 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2454 on 11 degrees of freedom
    ##   (21 observations deleted due to missingness)
    ## Multiple R-squared:  0.3032, Adjusted R-squared:  0.2399 
    ## F-statistic: 4.787 on 1 and 11 DF,  p-value: 0.05116
    ## 
    ## 
    ## Call:
    ## lm(formula = Grin1 ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.45126 -0.10158 -0.01727  0.09343  0.53051 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 10.49912    0.08413 124.795   <2e-16 ***
    ## PC1         -0.05300    0.01895  -2.797   0.0174 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.288 on 11 degrees of freedom
    ##   (21 observations deleted due to missingness)
    ## Multiple R-squared:  0.4156, Adjusted R-squared:  0.3625 
    ## F-statistic: 7.824 on 1 and 11 DF,  p-value: 0.01737
    ## 
    ## 
    ## Call:
    ## lm(formula = Prkci ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.46060 -0.13626 -0.05721  0.10958  0.82920 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  7.54286    0.09433  79.959   <2e-16 ***
    ## PC1         -0.01149    0.02125  -0.541    0.599    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3229 on 11 degrees of freedom
    ##   (21 observations deleted due to missingness)
    ## Multiple R-squared:  0.02592,    Adjusted R-squared:  -0.06264 
    ## F-statistic: 0.2927 on 1 and 11 DF,  p-value: 0.5993
    ## 
    ## 
    ## Call:
    ## lm(formula = Prkcz ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.40544 -0.09400 -0.04928  0.13286  0.46953 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  8.85333    0.07637 115.920   <2e-16 ***
    ## PC1         -0.02403    0.01720  -1.397     0.19    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2614 on 11 degrees of freedom
    ##   (21 observations deleted due to missingness)
    ## Multiple R-squared:  0.1506, Adjusted R-squared:  0.07342 
    ## F-statistic: 1.951 on 1 and 11 DF,  p-value: 0.19

    p2

![](../figures/02e_correlations/CandidatesvPCA-2.png)

    p3 <- plotPCAvcandidates("../data/02c_CA1vsd.csv")

    ## Joining, by = "mouse"
    ## Joining, by = "mouse"

    ## 
    ## Call:
    ## lm(formula = Arc ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.74215 -0.31878 -0.03302  0.31946  0.74213 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 10.25803    0.12156  84.389   <2e-16 ***
    ## PC1          0.02039    0.02674   0.763    0.459    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4315 on 13 degrees of freedom
    ##   (19 observations deleted due to missingness)
    ## Multiple R-squared:  0.04281,    Adjusted R-squared:  -0.03082 
    ## F-statistic: 0.5814 on 1 and 13 DF,  p-value: 0.4594
    ## 
    ## 
    ## Call:
    ## lm(formula = Camk2a ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.30203 -0.15118 -0.09120  0.09606  0.44466 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 12.67926    0.06426 197.300   <2e-16 ***
    ## PC1         -0.01784    0.01413  -1.262    0.229    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2281 on 13 degrees of freedom
    ##   (19 observations deleted due to missingness)
    ## Multiple R-squared:  0.1091, Adjusted R-squared:  0.04061 
    ## F-statistic: 1.593 on 1 and 13 DF,  p-value: 0.2291
    ## 
    ## 
    ## Call:
    ## lm(formula = Gria1 ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.31026 -0.06263  0.01453  0.05067  0.25455 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 11.566033   0.037230 310.667   <2e-16 ***
    ## PC1          0.007768   0.008188   0.949     0.36    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1321 on 13 degrees of freedom
    ##   (19 observations deleted due to missingness)
    ## Multiple R-squared:  0.06474,    Adjusted R-squared:  -0.007202 
    ## F-statistic: 0.8999 on 1 and 13 DF,  p-value: 0.3601
    ## 
    ## 
    ## Call:
    ## lm(formula = Grin1 ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.24085 -0.06279 -0.01132  0.07823  0.20505 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 11.043113   0.035748 308.915   <2e-16 ***
    ## PC1          0.005040   0.007863   0.641    0.533    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1269 on 13 degrees of freedom
    ##   (19 observations deleted due to missingness)
    ## Multiple R-squared:  0.03064,    Adjusted R-squared:  -0.04393 
    ## F-statistic: 0.4109 on 1 and 13 DF,  p-value: 0.5326
    ## 
    ## 
    ## Call:
    ## lm(formula = Prkci ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.58737 -0.10151  0.04483  0.19160  0.30661 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  7.62171    0.07434 102.531   <2e-16 ***
    ## PC1          0.01169    0.01635   0.715    0.487    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2638 on 13 degrees of freedom
    ##   (19 observations deleted due to missingness)
    ## Multiple R-squared:  0.03782,    Adjusted R-squared:  -0.0362 
    ## F-statistic: 0.511 on 1 and 13 DF,  p-value: 0.4874
    ## 
    ## 
    ## Call:
    ## lm(formula = Prkcz ~ PC1, data = allvsd)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.42541 -0.04775  0.00652  0.04935  0.29667 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  8.778114   0.046397 189.196   <2e-16 ***
    ## PC1         -0.004126   0.010205  -0.404    0.693    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1647 on 13 degrees of freedom
    ##   (19 observations deleted due to missingness)
    ## Multiple R-squared:  0.01242,    Adjusted R-squared:  -0.06355 
    ## F-statistic: 0.1635 on 1 and 13 DF,  p-value: 0.6925

    p3

![](../figures/02e_correlations/CandidatesvPCA-3.png)

    plot_grid(p1 + theme(axis.title.x = element_blank()),
              p2 + theme(axis.title.x = element_blank(),
                         strip.background = element_blank(), strip.text = element_blank()),
              p3 + theme(strip.background = element_blank(), strip.text = element_blank()), 
              nrow = 3, labels = c("DG", "CA3", "CA1"), label_size =  6,
              rel_heights = c(0.35,0.3,0.35))

![](../figures/02e_correlations/CandidatesvPCA-4.png)
