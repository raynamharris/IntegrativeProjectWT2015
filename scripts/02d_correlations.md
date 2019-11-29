    library(tidyverse)
    library(corrr)
    library(cowplot)

    source("./figureoptions.R")
    source("./functions_RNAseq.R")

    # read the sample data, set levels, join iwth behvior PCA data
    colData <- read.csv("../data/02a_colData.csv", row.names = 1, stringsAsFactors = T)
    colData <- colData %>% filter(subfield == "DG")
    pca.Rn <- read_csv("../data/01a_pca.all.csv") %>% filter(trialNum == 9)

    ## Parsed with column specification:
    ## cols(
    ##   ID = col_character(),
    ##   treatment = col_character(),
    ##   trialNum = col_double(),
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

    pca.Rn <- pca.Rn %>% select(ID:PC2)
    colData <- left_join(colData, pca.Rn)

    ## Joining, by = c("ID", "treatment")

    ## Warning: Column `ID` joining factor and character vector, coercing into
    ## character vector

    ## Warning: Column `treatment` joining factor and character vector, coercing
    ## into character vector

    head(colData)

    ##       ID subfield        treatment training trialNum        PC1
    ## 1 15143A       DG conflict.trained  trained        9 -0.2511636
    ## 2 15143B       DG   conflict.yoked    yoked        9 -3.5595365
    ## 3 15143D       DG   standard.yoked    yoked        9 -3.0605115
    ## 4 15144A       DG conflict.trained  trained        9  6.6921089
    ## 5 15144C       DG standard.trained  trained        9  7.0436756
    ## 6 15144D       DG   standard.yoked    yoked        9 -3.4999867
    ##           PC2
    ## 1  3.02981466
    ## 2 -0.47451513
    ## 3 -0.07227782
    ## 4 -0.02270055
    ## 5 -1.72946383
    ## 6  1.16766175

    # read all count data
    vsd <- read.csv("../data/02c_DGvsd.csv", row.names = 1, check.names = F) 
    vsd$gene <- row.names(vsd)
    vsd$gene <- toupper(vsd$gene)
    vsd <- as.data.frame(vsd)
    row.names(vsd) <- vsd$gene
    vsd$gene <- NULL
    head(vsd)

    ##               143A-DG-1 143B-DG-1 143D-DG-3 144A-DG-2 144C-DG-2 144D-DG-2
    ## 0610007P14RIK  6.381145  6.405587  6.834813  6.513173  6.401366  6.637581
    ## 0610009B22RIK  5.781124  5.661825  5.299160  5.528657  5.821796  5.670620
    ## 0610009L18RIK  5.405848  5.593055  5.422004  5.000789  5.507750  5.506319
    ## 0610009O20RIK  6.736671  6.691641  6.736719  6.530386  6.799892  6.642508
    ## 0610010F05RIK  6.603589  7.093384  7.135233  6.931921  6.907778  6.587132
    ## 0610010K14RIK  6.136536  5.949913  6.015194  5.782488  5.999853  6.084773
    ##               145A-DG-2 145B-DG-1 146A-DG-2 146B-DG-2 146C-DG-4 146D-DG-3
    ## 0610007P14RIK  6.562833  6.126760  6.298333  5.996243  6.929322  7.411893
    ## 0610009B22RIK  5.942028  5.787857  5.679794  5.996243  6.060620  5.000789
    ## 0610009L18RIK  5.256496  5.000789  5.679794  5.000789  5.000789  5.000789
    ## 0610009O20RIK  6.739906  6.623046  6.747840  6.383172  6.350262  6.793559
    ## 0610010F05RIK  6.739906  6.540236  6.765447  7.085062  7.264347  6.493162
    ## 0610010K14RIK  6.218650  6.200426  6.298333  7.252267  5.971751  5.000789
    ##               147C-DG-3 147D-DG-1 148A-DG-3 148B-DG-4
    ## 0610007P14RIK  6.298159  6.481476  6.493581  6.349976
    ## 0610009B22RIK  5.910805  5.893873  5.589261  5.492981
    ## 0610009L18RIK  5.257050  5.718223  5.505647  5.349661
    ## 0610009O20RIK  6.687072  6.631927  7.108382  6.425888
    ## 0610010F05RIK  6.880090  6.759923  6.892155  6.595557
    ## 0610010K14RIK  5.761732  6.136024  5.626697  6.081280

    # prep to join with sample colData
    vsd <- as.data.frame(t(vsd))
    vsd$sample <- row.names(vsd)
    vsd$mouse <- sapply(strsplit(as.character(vsd$sample),"\\-"), "[", 1)
    vsd$ID <- paste(15, vsd$mouse, sep = "")

    vsd <- left_join(colData, vsd)

    ## Joining, by = "ID"

    head(vsd)[17015:17020]

    ##     ZYG11B      ZYX    ZZEF1     ZZZ3    sample mouse
    ## 1 8.319205 7.203704 7.829641 7.742970 143A-DG-1  143A
    ## 2 8.118246 7.258360 7.711706 7.569642 143B-DG-1  143B
    ## 3 8.715286 6.853523 7.987434 7.135233 143D-DG-3  143D
    ## 4 9.110958 7.532625 8.240355 7.760152 144A-DG-2  144A
    ## 5 8.724875 7.128516 7.955995 7.634300 144C-DG-2  144C
    ## 6 8.814969 6.671645 7.908857 7.526267 144D-DG-2  144D

    forcorall <-  vsd %>% select(PC1:ZZZ3)

    corrrmat <- correlate(forcorall, diagonal = 1) 

    ## 
    ## Correlation method: 'pearson'
    ## Missing treated using: 'pairwise.complete.obs'

    head(corrrmat)

    ## # A tibble: 6 x 17,014
    ##   rowname     PC1     PC2 `0610007P14RIK` `0610009B22RIK` `0610009L18RIK`
    ##   <chr>     <dbl>   <dbl>           <dbl>           <dbl>           <dbl>
    ## 1 PC1      1      -0.589          -0.0955          0.331          -0.146 
    ## 2 PC2     -0.589   1               0.0157         -0.199           0.0177
    ## 3 061000… -0.0955  0.0157          1              -0.595          -0.173 
    ## 4 061000…  0.331  -0.199          -0.595           1               0.0266
    ## 5 061000… -0.146   0.0177         -0.173           0.0266          1     
    ## 6 061000…  0.0325 -0.163           0.164          -0.342           0.465 
    ## # … with 17,008 more variables: `0610009O20RIK` <dbl>,
    ## #   `0610010F05RIK` <dbl>, `0610010K14RIK` <dbl>, `0610012G03RIK` <dbl>,
    ## #   `0610030E20RIK` <dbl>, `0610037L13RIK` <dbl>, `0610040J01RIK` <dbl>,
    ## #   `1110002E22RIK` <dbl>, `1110004E09RIK` <dbl>, `1110004F10RIK` <dbl>,
    ## #   `1110008F13RIK` <dbl>, `1110008L16RIK` <dbl>, `1110008P14RIK` <dbl>,
    ## #   `1110012L19RIK` <dbl>, `1110017D15RIK` <dbl>, `1110032A03RIK` <dbl>,
    ## #   `1110032F04RIK` <dbl>, `1110034G24RIK` <dbl>, `1110037F02RIK` <dbl>,
    ## #   `1110038F14RIK` <dbl>, `1110051M20RIK` <dbl>, `1110059E24RIK` <dbl>,
    ## #   `1110059G10RIK` <dbl>, `1110065P20RIK` <dbl>, `1190002N15RIK` <dbl>,
    ## #   `1190005I06RIK` <dbl>, `1190007I07RIK` <dbl>, `1300017J02RIK` <dbl>,
    ## #   `1500009C09RIK` <dbl>, `1500009L16RIK` <dbl>, `1500011B03RIK` <dbl>,
    ## #   `1500011K16RIK` <dbl>, `1500015O10RIK` <dbl>, `1520401A03RIK` <dbl>,
    ## #   `1600002H07RIK` <dbl>, `1600002K03RIK` <dbl>, `1600012H06RIK` <dbl>,
    ## #   `1600014C10RIK` <dbl>, `1700001C19RIK` <dbl>, `1700001K19RIK` <dbl>,
    ## #   `1700001L19RIK` <dbl>, `1700001O22RIK` <dbl>, `1700001P01RIK` <dbl>,
    ## #   `1700003E16RIK` <dbl>, `1700003F12RIK` <dbl>, `1700006E09RIK` <dbl>,
    ## #   `1700007G11RIK` <dbl>, `1700007K13RIK` <dbl>, `1700008O03RIK` <dbl>,
    ## #   `1700010I14RIK` <dbl>, `1700011E24RIK` <dbl>, `1700011H14RIK` <dbl>,
    ## #   `1700011M02RIK` <dbl>, `1700013D24RIK` <dbl>, `1700013F07RIK` <dbl>,
    ## #   `1700014D04RIK` <dbl>, `1700015E13RIK` <dbl>, `1700016D06RIK` <dbl>,
    ## #   `1700016K19RIK` <dbl>, `1700017B05RIK` <dbl>, `1700019A02RIK` <dbl>,
    ## #   `1700019B03RIK` <dbl>, `1700019D03RIK` <dbl>, `1700019O17RIK` <dbl>,
    ## #   `1700020A23RIK` <dbl>, `1700020D05RIK` <dbl>, `1700020L24RIK` <dbl>,
    ## #   `1700020N01RIK` <dbl>, `1700021F05RIK` <dbl>, `1700022I11RIK` <dbl>,
    ## #   `1700023F06RIK` <dbl>, `1700024P16RIK` <dbl>, `1700025G04RIK` <dbl>,
    ## #   `1700027J19RIK` <dbl>, `1700028J19RIK` <dbl>, `1700028K03RIK` <dbl>,
    ## #   `1700028P14RIK` <dbl>, `1700029F12RIK` <dbl>, `1700029H14RIK` <dbl>,
    ## #   `1700029I15RIK` <dbl>, `1700029J07RIK` <dbl>, `1700030J22RIK` <dbl>,
    ## #   `1700030K09RIK` <dbl>, `1700034I23RIK` <dbl>, `1700034J05RIK` <dbl>,
    ## #   `1700037C18RIK` <dbl>, `1700037H04RIK` <dbl>, `1700040L02RIK` <dbl>,
    ## #   `1700047I17RIK2` <dbl>, `1700048O20RIK` <dbl>, `1700057G04RIK` <dbl>,
    ## #   `1700061G19RIK` <dbl>, `1700064H15RIK` <dbl>, `1700066B19RIK` <dbl>,
    ## #   `1700066M21RIK` <dbl>, `1700067K01RIK` <dbl>, `1700088E04RIK` <dbl>,
    ## #   `1700092M07RIK` <dbl>, `1700093K21RIK` <dbl>, `1700102P08RIK` <dbl>, …

    summary(lm( PC1 ~ ARC, data = vsd))

    ## 
    ## Call:
    ## lm(formula = PC1 ~ ARC, data = vsd)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -4.6977 -1.7429  0.2442  1.2672  5.6700 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -26.6287     5.4718  -4.867 0.000249 ***
    ## ARC           2.9210     0.5732   5.096 0.000163 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.617 on 14 degrees of freedom
    ## Multiple R-squared:  0.6497, Adjusted R-squared:  0.6247 
    ## F-statistic: 25.97 on 1 and 14 DF,  p-value: 0.0001629

    summary(lm( PC2 ~ ARC, data = vsd))

    ## 
    ## Call:
    ## lm(formula = PC2 ~ ARC, data = vsd)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.0549 -0.9078 -0.1875  0.5348  3.2513 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)   9.8519     2.9665   3.321  0.00505 **
    ## ARC          -0.9469     0.3107  -3.047  0.00870 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.419 on 14 degrees of freedom
    ## Multiple R-squared:  0.3988, Adjusted R-squared:  0.3558 
    ## F-statistic: 9.285 on 1 and 14 DF,  p-value: 0.008699

    cor.test(vsd$PC1, vsd$ARC, method = c("pearson"))

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  vsd$PC1 and vsd$ARC
    ## t = 5.0961, df = 14, p-value = 0.0001629
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.5168971 0.9301213
    ## sample estimates:
    ##       cor 
    ## 0.8060654

    cor.test(vsd$PC2, vsd$ARC, method = c("pearson"))

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  vsd$PC2 and vsd$ARC
    ## t = -3.0471, df = 14, p-value = 0.008699
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.8584581 -0.1976259
    ## sample estimates:
    ##        cor 
    ## -0.6314698

    a <- ggplot(vsd, aes(x = ARC, y = PC1)) +
       geom_point(aes( color = treatment)) + 
       geom_smooth(method = "lm", color = "grey") +
       scale_color_manual(values = treatmentcolors) +
      theme_ms() +
       theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank())  +
      labs(subtitle = "R2 = 0.81, p = 0.0002")
    b <- ggplot(vsd, aes(x = ARC, y = PC2)) +
       geom_point(aes( color = treatment)) + 
       geom_smooth(method = "lm", color = "grey") +
       scale_color_manual(values = treatmentcolors) +
      theme_ms() +
       theme(legend.position = "none")  +
      labs(subtitle = "R2 = -0.63, p = 0.009 ")

    left <- plot_grid(a,b, nrow = 2, labels = c("a","b"), label_size = 8)
    left

![](../figures/02e_correlations/corrrplots-1.png)

    slim <- corrrmat %>% 
      focus(PC1, PC2)  %>% 
      arrange(desc(PC1)) %>% 
      filter(PC1 > 0.78) 
    slim

    ## # A tibble: 11 x 3
    ##    rowname    PC1    PC2
    ##    <chr>    <dbl>  <dbl>
    ##  1 NAF1     0.859 -0.693
    ##  2 PTGS2    0.839 -0.640
    ##  3 RGS2     0.833 -0.605
    ##  4 HIST1H1D 0.818 -0.539
    ##  5 COL10A1  0.813 -0.574
    ##  6 ARC      0.806 -0.631
    ##  7 HSPB3    0.801 -0.277
    ##  8 NPAS4    0.798 -0.637
    ##  9 FZD5     0.797 -0.666
    ## 10 ACAN     0.789 -0.542
    ## 11 AREG     0.784 -0.432

    topcorrrs <- slim$rowname

    p1 <- corrrmat %>% 
      focus(PC1, PC2, topcorrrs,  mirror = TRUE) %>% 
       focus(PC1, PC2) %>%
      mutate(rowname = reorder(rowname, PC1)) %>%
      ggplot(aes(rowname, PC1, fill = PC1)) +
        geom_col() + coord_flip() +
      scale_fill_gradient2(low = "#67a9cf",  high = "#ef8a62", midpoint = 0) +
      theme_ms() +
      theme(legend.position = "none") +
      #ylim(-1,1) +
      labs(x = NULL, y = "Correlation to PC1") +
      labs(subtitle = " ")
    p1 

![](../figures/02e_correlations/corrrplots-2.png)

    p2 <- corrrmat %>% 
      focus(PC1, PC2, topcorrrs,  mirror = TRUE) %>% 
      rearrange() %>% 
      network_plot(colors = c("#67a9cf", "white", "#ef8a62"),
                   min_cor = .6, curved = T, legend = T,
                   repel = TRUE) + 
      theme(legend.position = "bottom") +
      labs(subtitle = " ")

    ## Registered S3 method overwritten by 'seriation':
    ##   method         from 
    ##   reorder.hclust gclus

    p2 

![](../figures/02e_correlations/corrrplots-3.png)

    right <- plot_grid(p1,p2,  nrow = 1, labels = c("c","d"), label_size = 8)
    right

![](../figures/02e_correlations/corrrplots-4.png)

    fig4 <- plot_grid(left, right, nrow = 1, rel_widths = c(1,2))
    fig4

![](../figures/02e_correlations/corrrplots-5.png)

    pdf(file="../figures/02e_correlations/corrrplots.pdf", width=6.69, height=3.5)
    plot(fig4)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    pdf(file="../figures/figure_4.pdf", width=6.69, height=3.5)
    plot(fig4)
    dev.off()

    ## quartz_off_screen 
    ##                 2

what genes genes correlated with PC1? what gene are correlated with PC1 and are differentially epxressed
--------------------------------------------------------------------------------------------------------

    # 57 genes are correlated with PC1 > 0.7
    PC1corrs <- corrrmat %>% 
      focus(PC1)  %>% 
      arrange(desc(PC1)) %>% 
      filter(PC1 > .7) %>%
       mutate_each(funs=toupper) %>%
      mutate(gene = rowname) %>% 
      arrange(gene) %>% 
      select(gene,PC1)   
    PC1corrs$gene

    ##  [1] "ABRA"     "ACAN"     "ADAMTS1"  "AIM1L"    "AMIGO2"   "ARC"     
    ##  [7] "AREG"     "ARL4D"    "ARL5B"    "ARMCX5"   "ATF3"     "BDNF"    
    ## [13] "BTG2"     "CHST14"   "COL10A1"  "CPEB4"    "DUSP8"    "EGR1"    
    ## [19] "EGR4"     "ERRFI1"   "FBXO33"   "FOSL2"    "FRMD6"    "FZD5"    
    ## [25] "GM10269"  "HIST1H1D" "HIST1H3I" "HOMER1"   "HOXC4"    "HSPB3"   
    ## [31] "IGHD"     "KCNK10"   "MEST"     "NAF1"     "NEXN"     "NFIL3"   
    ## [37] "NPAS4"    "NR4A3"    "PCDH8"    "PELI1"    "PER1"     "PLK2"    
    ## [43] "PPIC"     "PTGS2"    "RASD1"    "RGS2"     "SGK1"     "SLC16A1" 
    ## [49] "SLC25A25" "SMAD7"    "SPTY2D1"  "SYT4"     "TIPARP"   "TNFRSF23"
    ## [55] "TRIB1"    "UBC"      "ZFP804B"

    # 217 DEGs in DG
    dim(DG_DEGs)

    ## [1] 214  17

    head(DG_DEGs$gene)

    ## [1] "1190002N15RIK" "A830010M20RIK" "ABHD2"         "ACAN"         
    ## [5] "ADAMTS1"       "ADRB1"

    # 41 found  in both correlate and deg datasets
    PC1corrsDEGS <- inner_join(PC1corrs, DG_DEGs) %>% 
      arrange(gene) %>% select(gene,PC1)   

    ## Joining, by = "gene"

    length(PC1corrsDEGS$gene)

    ## [1] 41

    PC1corrsDEGS$gene

    ##  [1] "ACAN"     "ADAMTS1"  "AMIGO2"   "ARC"      "ARL4D"    "ARL5B"   
    ##  [7] "ARMCX5"   "ATF3"     "BDNF"     "BTG2"     "CPEB4"    "DUSP8"   
    ## [13] "EGR1"     "EGR4"     "ERRFI1"   "FBXO33"   "FOSL2"    "FRMD6"   
    ## [19] "FZD5"     "HOMER1"   "MEST"     "NAF1"     "NFIL3"    "NPAS4"   
    ## [25] "NR4A3"    "PCDH8"    "PELI1"    "PER1"     "PLK2"     "PTGS2"   
    ## [31] "RASD1"    "RGS2"     "SGK1"     "SLC16A1"  "SLC25A25" "SMAD7"   
    ## [37] "SPTY2D1"  "SYT4"     "TIPARP"   "TRIB1"    "UBC"

    ## these 41 genes were used for a go analysis using shinygo http://bioinformatics.sdstate.edu/go/
    ## results are stored as 
      #data/02e_DEGPC1enrichmentBP.csv
        #data/02e_DEGPC1enrichmentCC.csv
        #data/02e_DEGPC1enrichmentMF.csv

    BP <- read_csv("../data/02e_DEGPC1enrichmentBP.csv", n_max = 5) %>% mutate(Domain = "BP") 

    ## Parsed with column specification:
    ## cols(
    ##   `Enrichment FDR` = col_double(),
    ##   `Genes in list` = col_double(),
    ##   `Total genes` = col_double(),
    ##   `Functional Category` = col_character(),
    ##   Genes = col_character()
    ## )

    CC <- read_csv("../data/02e_DEGPC1enrichmentCC.csv", n_max = 5) %>% mutate(Domain = "CC") 

    ## Parsed with column specification:
    ## cols(
    ##   `Enrichment FDR` = col_double(),
    ##   `Genes in list` = col_double(),
    ##   `Total genes` = col_double(),
    ##   `Functional Category` = col_character(),
    ##   Genes = col_character()
    ## )

    MF <- read_csv("../data/02e_DEGPC1enrichmentMF.csv", n_max = 5) %>% mutate(Domain = "MF") 

    ## Parsed with column specification:
    ## cols(
    ##   `Enrichment FDR` = col_double(),
    ##   `Genes in list` = col_double(),
    ##   `Total genes` = col_double(),
    ##   `Functional Category` = col_character(),
    ##   Genes = col_character()
    ## )

    GO <- rbind(BP, CC, MF) %>% select(Domain, `Functional Category`, `Total genes`, `Genes in list`, Genes)
    GO

    ## # A tibble: 15 x 5
    ##    Domain `Functional Category`  `Total genes` `Genes in list` Genes       
    ##    <chr>  <chr>                          <dbl>           <dbl> <chr>       
    ##  1 BP     Memory                           145               9 BDNF ARC PT…
    ##  2 BP     Tissue development              1855              21 ERRFI1 BDNF…
    ##  3 BP     Learning or memory               286              10 BDNF ARC PT…
    ##  4 BP     Cognition                        317              10 BDNF ARC PT…
    ##  5 BP     Behavior                         693              13 BDNF ARC PT…
    ##  6 CC     Neuron projection               1486              13 HOMER1 CPEB…
    ##  7 CC     Cell junction                   1095              10 FZD5 FRMD6 …
    ##  8 CC     Dendrite                         736               8 HOMER1 BDNF…
    ##  9 CC     Dendritic tree                   738               8 HOMER1 BDNF…
    ## 10 CC     Neuron part                     1933              13 HOMER1 CPEB…
    ## 11 MF     Regulatory region nuc…           949              11 PER1 ATF3 N…
    ## 12 MF     Transcription regulat…           946              11 PER1 ATF3 N…
    ## 13 MF     RNA polymerase II reg…           782               9 ATF3 NR4A3 …
    ## 14 MF     DNA-binding transcrip…           738               9 EGR1 NPAS4 …
    ## 15 MF     RNA polymerase II reg…           788               9 ATF3 NR4A3 …

    # alphabetize
    GO$Genes <- GO$Genes %>% str_split(., ' ') %>% lapply(., 'sort') %>%  lapply(., 'paste', collapse=' ') %>% unlist(.)
    GO

    ## # A tibble: 15 x 5
    ##    Domain `Functional Category`  `Total genes` `Genes in list` Genes       
    ##    <chr>  <chr>                          <dbl>           <dbl> <chr>       
    ##  1 BP     Memory                           145               9 ARC BDNF EG…
    ##  2 BP     Tissue development              1855              21 ACAN ARC AR…
    ##  3 BP     Learning or memory               286              10 ARC BDNF BT…
    ##  4 BP     Cognition                        317              10 ARC BDNF BT…
    ##  5 BP     Behavior                         693              13 ARC BDNF BT…
    ##  6 CC     Neuron projection               1486              13 ACAN ARC BD…
    ##  7 CC     Cell junction                   1095              10 ARC CPEB4 F…
    ##  8 CC     Dendrite                         736               8 ARC BDNF CP…
    ##  9 CC     Dendritic tree                   738               8 ARC BDNF CP…
    ## 10 CC     Neuron part                     1933              13 ACAN ARC BD…
    ## 11 MF     Regulatory region nuc…           949              11 ATF3 EGR1 E…
    ## 12 MF     Transcription regulat…           946              11 ATF3 EGR1 E…
    ## 13 MF     RNA polymerase II reg…           782               9 ATF3 EGR1 E…
    ## 14 MF     DNA-binding transcrip…           738               9 ATF3 BTG2 E…
    ## 15 MF     RNA polymerase II reg…           788               9 ATF3 EGR1 E…

    citation("corrr")

    ## 
    ## To cite package 'corrr' in publications use:
    ## 
    ##   Edgar Ruiz, Simon Jackson and Jorge Cimentada (2019). corrr:
    ##   Correlations in R. R package version 0.4.0.
    ##   https://CRAN.R-project.org/package=corrr
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {corrr: Correlations in R},
    ##     author = {Edgar Ruiz and Simon Jackson and Jorge Cimentada},
    ##     year = {2019},
    ##     note = {R package version 0.4.0},
    ##     url = {https://CRAN.R-project.org/package=corrr},
    ##   }

    write.csv(slim, "../data/02e_top10PC1correlations.csv", row.names = F)
    write.csv(PC1corrs, "../data/02e_PC1correlations.csv", row.names = F)
    write.csv(GO, "../data/02e_PC1corrTopGOterms.csv", row.names = F)
